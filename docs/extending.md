# Extending ops-post

## Adding a new element type

Adding an element to `ELEMENT_REGISTRY` in `ops_elements.py` is the
**only change needed**. The registry drives TCL parsing, mesh building,
and result processing automatically.

### Step 1: If the element fits an existing topology

Just add it to `ELEMENT_REGISTRY`:

```python
"MyNewShellQ4": _Q4,
```

This gives it Q4 shape functions, 2x2 Gauss points, 4-node surface
parsing, and Q4 extrapolation.

### Step 2: If it needs a new topology

Define the GP coordinates, shape functions, and extrapolation matrix,
then create a new `ElemInfo`:

```python
GP_TRI_6 = np.array([...])  # natural coordinates

def shape_t6(xi, eta):
    """6-node triangle shape functions."""
    ...
    return N  # (6,) array

_T6 = ElemInfo(
    family=ElemFamily.SHELL_TRI,
    num_nodes=6, num_gp=3,
    gp_coords=GP_TRI_3,
    shape_fn=shape_t6,
    extrap_fn=None,
    face_nodes=[0, 1, 2, 3, 4, 5],
)

ELEMENT_REGISTRY["MyTriElement"] = _T6
```

### What happens automatically

Once the element is in the registry:

| System | Behavior |
|--------|----------|
| **TCL parser** | Uses `ElemFamily` to pick the right parser (surface, beam, truss, zero-length) and reads the correct number of nodes |
| **Postdata generator** | Computes local axes from `-local`, `geomTransf`, or `-orient` flags |
| **Mesh builder** | Builds surface faces, extrusion, fiber edges, and GP cloud using the registry's shape functions and GP coordinates |
| **Result processor** | Uses the registry's extrapolation matrix for GP-to-node mapping |

### Adding beam inline property extraction

If a new beam type has inline section properties in a different argument
layout, add it to `_BEAM_INLINE_PROPS` in `tcl_parser.py`:

```python
_BEAM_INLINE_PROPS = {
    "Timoshenko":        (7, 9, 10),   # A at token[7], Iy at [9], Iz at [10]
    "elasticBeamColumn": (5, 9, 10),
    "myNewBeam":         (6, 8, 9),    # new
}
```

## Adding a new result type

Nodal results are auto-discovered from `RESULTS/ON_NODES/` in the HDF5 file.
Element results are auto-discovered from `RESULTS/ON_ELEMENTS/`.

To customize component labels for a new nodal result, add it to
`NODAL_COMPONENTS` in `result_processor.py`:

```python
NODAL_COMPONENTS = {
    "DISPLACEMENT": ["Ux", "Uy", "Uz", "|U|"],
    "REACTION_FORCE": ["Fx", "Fy", "Fz", "|F|"],
    "VELOCITY": ["Vx", "Vy", "Vz", "|V|"],  # new
}
```

Element result component names are read from the HDF5 `META/COMPONENTS`
field automatically.

## Adding a new derived quantity

To add a derived scalar (like von Mises), modify `result_processor.py`:

1. Add the name to the component list in `get_available_results()`:
   ```python
   if sg0.num_components >= 6:
       comps.append("Tresca")
   ```

2. Handle it in `_extract_gp_fiber_component()` and
   `extract_element_result_for_gp_cloud()`, using the component index
   to detect when the derived quantity is selected (index >= `n_comp`).

## Adding a new visualization layer

1. Build the mesh in `MeshBuilder` (e.g. `build_my_layer()`).

2. In `MainWindow.__init__`, call the builder and store the base mesh.

3. In `_full_rebuild`, add it with `add_mesh(...)` and store in
   `self._meshes["my_layer"]` / `self._actors["my_layer"]`.

4. In `_fast_update`, update its `.points[:]` from the base + displacement.

5. Add a visibility checkbox in the **Results tab** display options if
   needed.

## Customizing the View tab

The View tab in `gui.py` contains camera presets and figure controls.
Key points for extension:

- **GP zero mode**: controlled by a dropdown in the Figure section of
  the View tab. To add a new zero mode, update the dropdown items and
  the sphere radius calculation in `result_processor.py`.
- **Scale range clamping**: min/max spinboxes in the Figure section.
  To add more clamping options, extend the spinbox layout and the
  scalar array clamping logic in `_apply_scale_range()`.
- **All settings apply immediately**: there is no Apply button. Connect
  new controls to their update methods directly (e.g.
  `spinbox.valueChanged.connect(self._full_rebuild)`).

## Customizing interaction

The arcball interactor is in `utils.py` (`ArcballInteractorStyle`).
It uses the Shoemake arcball algorithm:

- Mouse positions are projected onto a virtual unit sphere.
- A single rotation is computed from drag-start to current position.
- The rotation is applied to the saved initial camera state (no drift).

To change the rotation pivot from world origin to another point, modify
`_arcball_rotate()` to rotate around a different center.

To add left-click pick/inspect, implement `_on_left_press` and use
VTK's cell picker to identify the clicked element.
