# Architecture

## Module overview

```
ops_post/
├── __main__.py          CLI entry point (same-name convention)
├── tcl_parser.py        TCL model parser + .postdata generator
├── model.py             Data model (dataclasses)
├── mpco_reader.py       File I/O (HDF5 + .postdata)
├── ops_elements.py      Element type registry (61 types)
├── mesh_builder.py      PyVista mesh construction
├── result_processor.py  Result extraction + mapping
├── gui.py               PyQt5 tabbed UI (Results + View), in-place
│                        mesh updates, persistent HDF5, export
└── utils.py             Arcball interactor, von Mises, utilities
```

## Workflow

```
                         ops-post my_model
                               |
                    ┌──────────┴──────────┐
                    v                     v
              my_model.tcl          my_model.mpco
                    |                     |
                    v                     v
             tcl_parser.py         mpco_reader.py
             (parse nodes,         (read HDF5:
              elements,             nodes, elements,
              sections,             sections, results)
              geomTransf)                |
                    |                    |
                    v                    |
         my_model.mpco.postdata         |
         (local axes, beam              |
          profiles, elem info)          |
                    |                    |
                    └────────┬───────────┘
                             v
                         ModelState
                             |
              ┌──────────────┼──────────────┐
              v              v              v
        MeshBuilder   ResultProcessor   MainWindow
        (PyVista        (scalar         (Qt + PyVista
         meshes)         arrays)         viewport)
```

## Data flow

1. **`__main__.py`** resolves the base name, checks if `.postdata` needs
   regeneration, calls `tcl_parser.parse_tcl()` + `write_postdata()` if
   needed.

2. **`tcl_parser.parse_tcl(path)`** reads the TCL model, following all
   `source` commands. Returns a `TclModel` with nodes, elements, sections,
   materials, and geomTransf definitions. Element parsing is dispatched
   by the `ops_elements` registry -- the element's `ElemFamily` determines
   how many nodes to read and what flags to extract.

3. **`tcl_parser.write_postdata(model, path)`** computes local axes
   (rotation matrices -> quaternions), derives beam cross-section profiles
   from section properties, and writes the `.mpco.postdata` file.

4. **`mpco_reader.read_mpco(path)`** reads the HDF5 file and the
   `.mpco.postdata` companion, returning a `ModelState`.

5. **`MeshBuilder(model)`** converts the model into PyVista meshes.

6. **`ResultProcessor(model, builder)`** extracts per-step scalar data from
   HDF5 and maps it onto the meshes.

7. **`MainWindow`** ties it together with a tabbed control panel (Results
   and View tabs), in-place mesh updates during animation, and export
   controls. The viewport has a fixed size (1200x800 default, adjustable
   via View tab spinners). All figure settings apply immediately -- no
   Apply button.

## Key dataclasses (model.py)

| Class              | Purpose                                              |
|--------------------|------------------------------------------------------|
| `NodeData`         | Node IDs, coordinates, ID-to-index map               |
| `ElementGroup`     | One element type: IDs, connectivity, HDF5 key        |
| `SectionAssignment`| Fiber data, thickness, element-to-section mapping     |
| `BeamProfile`      | 2D cross-section vertices for beam visualization      |
| `ResultDescriptor` | Describes one result (category, name, sub-groups)     |
| `ResultSubGroup`   | One HDF5 data group: elem IDs, GP/fiber/comp counts  |
| `StageInfo`        | Stage index and its step keys                         |
| `ModelState`       | Top-level container for the entire model              |

## Element registry (ops_elements.py)

Maps 61 OpenSees element type strings to `ElemInfo` objects. This single
registry drives **three** systems:

1. **TCL parsing**: determines how many nodes to read and which flags
   to extract (`-local`, `-orient`, geomTransf).
2. **Mesh building**: provides shape functions, GP coordinates, and node
   counts for surface/extrusion/GP cloud construction.
3. **Result processing**: provides extrapolation matrices for GP-to-node
   mapping.

```python
info = lookup("ASDShellQ4")
info.family      # ElemFamily.SHELL_QUAD
info.num_nodes   # 4
info.num_gp      # 4
info.gp_coords   # (4, 2) array at +/-1/sqrt(3)
info.shape_fn    # shape_q4(xi, eta) -> (4,)
info.extrap_fn   # extrapolation_matrix_q4() -> (4, 4)
info.is_surface  # True
```

Adding a new element type to `ELEMENT_REGISTRY` automatically handles
TCL parsing, mesh building, and result display.

## TCL parser (tcl_parser.py)

Parses OpenSees TCL commands:

| Command | What is extracted |
|---------|-------------------|
| `model basic -ndm -ndf` | Dimensionality |
| `node` | Node ID, coordinates |
| `element` | ID, type, node IDs, section tag, -local, -orient, geomTransf |
| `geomTransf` | ID, type, vecxz vector |
| `section` | ID, type, LayeredShell layers |
| `nDMaterial` / `uniaxialMaterial` | ID, type |
| `source` | Follows included files recursively |

Element parsing is **registry-driven**: `lookup(elem_type)` returns the
`ElemFamily`, which determines the parser:

| Family | Parser | Extracts |
|--------|--------|----------|
| `SHELL_QUAD`, `SHELL_TRI`, `QUAD_2D`, `TRI_2D` | `_parse_surface_element` | N nodes + section tag + `-local` |
| `BEAM` | `_parse_beam` | 2 nodes + geomTransf + inline A/Iy/Iz |
| `TRUSS`, `LINK` | `_parse_two_node` | 2 nodes |
| `ZERO_LENGTH` | `_parse_zero_length` | 2 nodes + `-orient` |
| Unknown | `_parse_generic_element` | Consecutive integer node IDs |

## Render update strategy

The GUI uses two update paths:

### Full rebuild (`_full_rebuild`)
Triggered by: visibility toggles, colormap change, result/component change,
grid spacing change.

- Clears the scene.
- Rebuilds all meshes from base geometry + current displacement.
- Stores mesh and actor references in `_meshes` / `_actors` dicts.

### Fast update (`_fast_update`)
Triggered by: time step change, displacement scale change.

- Updates `.points[:]` arrays in-place on existing meshes (no clear/add).
- Updates scalar arrays for contour results.
- Only GP spheres are removed and re-added (glyph geometry depends on values).

## HDF5 access

A single `h5py.File` handle is opened in `MainWindow.__init__` and kept
for the session (persistent HDF5 handle). All `read_step_data()` calls
pass this handle to avoid repeated file open/close. The handle is closed
in `closeEvent`.

## GP sphere zero mode

When displaying GP spheres, the zero mode controls sphere radius scaling:

- **"Zero at 0"**: radius is zero when the value is zero (for signed
  quantities like stress).
- **"Zero at min"**: radius is zero at the minimum value (for damage-like
  quantities that range from 0 upward).

This setting is in the View tab under Figure controls and applies
immediately.

## Displacement interpolation

GP cloud displacement uses vectorized numpy:

```python
# Precomputed in build_gauss_point_cloud():
_gp_node_indices  # (N_gp, max_nodes) int array
_gp_weights       # (N_gp, max_nodes) float array

# At render time:
d = displacements[node_indices]          # (N_gp, max_nodes, 3)
interp = np.einsum("ij,ijk->ik", wts, d) # (N_gp, 3)
```
