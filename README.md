# ops-post

Interactive 3D post-processor for plotting OpenSees masonry wall results.

Reads OpenSees TCL model files and HDF5 `.mpco` results to produce
publication-ready plots of masonry wall models — stress contours,
displacement fields, fiber results, and animated pushover sequences.

## Installation

```bash
cd ops-post
pip install .
```

## Usage

ops-post uses a **same-name convention**: the TCL model and MPCO results
share a base name.

```
my_model.tcl     <-- OpenSees model (single file or sources others)
my_model.mpco    <-- analysis results
```

Then just run:

```bash
ops-post my_model
```

This will:

1. Parse `my_model.tcl` (follows all `source` commands automatically)
2. Generate `my_model.mpco.postdata` (local axes, beam profiles, element info)
3. Open the interactive 3D viewer with results from `my_model.mpco`

On subsequent runs, the `.postdata` is reused unless the `.tcl` is newer.

You can also pass the `.mpco` or `.tcl` extension explicitly, or launch
a file dialog with just `ops-post`.

## Controls

| Input              | Action                          |
|--------------------|---------------------------------|
| Right-click drag   | Arcball rotation around origin  |
| Middle-click drag  | Pan                             |
| Scroll             | Zoom (perspective + orthographic) |

## Features

### Tabbed interface
- **Results tab**: result/component selection, element display mode,
  displacement scale, visibility toggles, time step animation.
- **View tab**: camera presets, figure controls, colormap, scale bar
  settings, export.

### Result display
- **Nodal results** (displacement, reaction force) as smooth contours on the
  mid-surface mesh.
- **Element material results** as contours (GP-averaged or extrapolated) or
  as colored Gauss-point spheres with selectable GP and fiber layer.
- **Element fiber results** (`section.fiber.stress/strain`) shown as scaled
  spheres at every GP x fiber location. Sphere radius is proportional to
  the absolute value. No averaging -- every integration point shows its
  real value.
- **GP sphere zero mode**: "Zero at 0" for stress (+/- values) or
  "Zero at min" for damage-like quantities (0-to-max).
- **Von Mises** stress offered when 3+ components are available, using the
  5-component formula when applicable.

### Visualization layers
- Transparent shell extrusion showing section thickness.
- Mid-surface quads/triangles for contour plots.
- Fiber layer edge lines (auto-hidden when GP spheres are active).
- Beam elements rendered as 3D hexahedra from cross-section profiles
  (derived from TCL section properties).
- Floor grid at the model base, adaptive spacing, toggleable.

### Camera
- Standard views: Front, Back, Left, Right, Top, Bottom, Iso.
- Perspective / Orthographic projection toggle.
- True Shoemake arcball rotation -- no drift on circular mouse paths.

### Publication-ready figures
- Fixed viewport size (1200x800 default), adjustable via spinners.
- Aspect ratio presets (16:9, 4:3, 3:2, 1:1, 2:1).
- Scale bar: vertical/horizontal orientation, position presets
  (Top-Left, Top-Right, Bottom-Left, Bottom-Right), adjustable width,
  height, font size, and custom title override.
- Scale range clamping: optional min/max with value clamping.
- Colormaps: jet, coolwarm, viridis, plasma, turbo, RdYlBu_r, hot_r,
  YlOrRd, inferno_r.
- All settings apply immediately -- no Apply button needed.

### Export
- **Screenshot**: PNG (1x, 2x, 3x, 4x resolution), SVG vector, PDF vector.
- **Animation recording**: choose file (GIF or MP4) first, then
  auto-plays from current step, auto-stops and saves at the last step.
  MP4 uses bundled ffmpeg via imageio (falls back to GIF if unavailable).

### Supported elements
The element registry (`ops_elements.py`) covers 61 OpenSees element types.
The same registry drives both **TCL parsing** and **mesh visualization** --
adding a new element type to the registry is all that's needed.

| Family       | Examples                                    | Nodes | GPs   |
|--------------|---------------------------------------------|-------|-------|
| Shell quad   | ASDShellQ4, ShellMITC4, ShellDKGQ           | 4     | 2x2   |
| Shell tri    | ASDShellT3, ShellDKGT                       | 3     | 3     |
| Shell Q9     | ShellMITC9                                  | 9     | 3x3   |
| 2D quad      | quad, SSPquad, enhancedQuad                 | 4     | 2x2   |
| 2D tri       | Tri31                                       | 3     | 1     |
| Beam/column  | ElasticTimoshenkoBeam, forceBeamColumn, ...  | 2     | --    |
| Truss        | truss, corotTruss                           | 2     | --    |
| Link         | twoNodeLink                                 | 2     | --    |
| Zero-length  | zeroLength, zeroLengthSection, ...           | 2     | --    |

Unknown elements are handled by guessing from the connectivity column count.

### Performance
- **In-place mesh updates**: stepping through time only updates point
  positions and scalar arrays -- no scene rebuild.
- **Vectorized GP displacement**: all Gauss-point positions are interpolated
  in one `numpy.einsum` call.
- **Persistent HDF5 handle**: the `.mpco` file is opened once and kept open
  for the session, avoiding repeated open/close overhead during animation.

### Automatic detection
- **Up-axis**: detected from the node bounding box (largest extent). Floor
  grid and fallback normals orient accordingly.
- **Component names**: read from the HDF5 `META/COMPONENTS` field.
- **Grid spacing**: scales with model dimensions (order of magnitude).

## Dependencies

| Package    | Purpose                         | License                |
|------------|---------------------------------|------------------------|
| h5py       | Read `.mpco` HDF5 files         | BSD-3-Clause           |
| numpy      | Array operations                | BSD-3-Clause           |
| scipy      | Quaternion conversion           | BSD-3-Clause           |
| pyvista    | 3D mesh rendering               | MIT                    |
| pyvistaqt  | PyVista-Qt integration          | MIT                    |
| PyQt5      | GUI framework                   | GPL-3.0                |
| imageio    | MP4 video export (bundled ffmpeg) | BSD-2-Clause         |

Transitive dependencies:

| Package         | Purpose                   | License                     |
|-----------------|---------------------------|-----------------------------|
| VTK             | 3D rendering engine       | BSD-3-Clause                |
| Pillow          | GIF export, image I/O     | HPND (permissive)           |
| imageio-ffmpeg  | Bundled ffmpeg binary     | BSD-2-Clause (ffmpeg: LGPL) |

All licenses are compatible with GPL-3.0. PyQt5 is itself GPL-3.0,
which is why ops-post uses the same license. All dependencies are
installed automatically by `pip install .`.

**OpenSees** (BSD-3-Clause) is the analysis engine that produces `.mpco`
result files. ops-post does not link to or include any OpenSees code --
it reads the `.mpco` HDF5 output and parses the TCL model input
independently. The `.mpco` format is an HDF5 data layout, not a
separately licensed specification.

## Project structure

```
ops-post/
├── pyproject.toml
├── LICENSE              GPL-3.0-or-later
└── src/
    └── ops_post/
        ├── __main__.py          Entry point (same-name convention)
        ├── tcl_parser.py        TCL model parser + .postdata generator
        ├── gui.py               PyQt5 main window (Results + View tabs)
        ├── mesh_builder.py      PyVista mesh construction
        ├── model.py             Dataclasses (nodes, elements, sections, results)
        ├── mpco_reader.py       HDF5 + .postdata file reader
        ├── ops_elements.py      Element type registry (61 types)
        ├── result_processor.py  Result extraction + GP/node mapping
        └── utils.py             Arcball interactor, von Mises, utilities
```

## Acknowledgements

This project was developed with the help of [Claude Code](https://claude.ai/claude-code).

## License

ops-post is licensed under the [GNU General Public License v3.0](LICENSE).

You are free to use, modify, and distribute this software. If you distribute
modified versions, you must also make your source code available under the
same license. See the [LICENSE](LICENSE) file for full terms.
