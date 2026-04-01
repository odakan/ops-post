# ops-post

Interactive 3D post-processor for plotting OpenSees masonry wall results.

<div style="display: flex; gap: 1em; margin: 1.5em 0;">
<div style="flex:1; padding: 1em; border: 1px solid #ddd; border-radius: 8px;">
<strong>Parses</strong><br>
OpenSees TCL model files for wall geometry and local axes
</div>
<div style="flex:1; padding: 1em; border: 1px solid #ddd; border-radius: 8px;">
<strong>Plots</strong><br>
Stress contours, displacement fields, fiber results, pushover animations
</div>
<div style="flex:1; padding: 1em; border: 1px solid #ddd; border-radius: 8px;">
<strong>Exports</strong><br>
Publication-ready PNG, SVG, PDF figures and GIF/MP4 animations
</div>
</div>

## Quick start

```bash
cd ops-post
pip install .
ops-post my_model
```

ops-post finds `my_model.tcl` and `my_model.mpco` by the same-name
convention. See [Getting Started](getting-started.md) for details.

## Why ops-post?

OpenSees lacks a built-in interactive post-processor for visualizing
masonry wall results. **ops-post** is a standalone tool that works with
any Python installation and gives you:

- Automatic TCL model parsing -- no manual companion files
- Transparent shell extrusion showing section thickness
- Fiber-level Gauss point visualization with scaled spheres
- True arcball rotation (no drift)
- Animated time stepping with in-place mesh updates
- Auto-detected up-axis and floor grid
- **Tabbed interface**: Results tab for result/component selection and
  display options; View tab for camera presets, figure controls,
  colormap, scale bar, and export
- **Publication-ready figures**: fixed viewport size, aspect ratio
  presets, scale bar customization, scale range clamping, multiple
  colormaps, PNG/SVG/PDF export, GIF/MP4 animation recording
- **GP sphere zero mode**: choose "Zero at 0" for stress-like quantities
  or "Zero at min" for damage-like quantities

## Built with

| Package | Purpose |
|---------|---------|
| [h5py](https://www.h5py.org/) | Read `.mpco` HDF5 files |
| [PyVista](https://docs.pyvista.org/) | 3D mesh rendering |
| [PyQt5](https://www.riverbankcomputing.com/software/pyqt/) | GUI framework |
| [NumPy](https://numpy.org/) / [SciPy](https://scipy.org/) | Numerics |
| [imageio](https://imageio.readthedocs.io/) | MP4 video export (bundled ffmpeg) |

## License

ops-post is licensed under the
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).
PyQt5 is itself GPL-3.0, which is why ops-post uses the same license.
All other dependencies (h5py, NumPy, SciPy, PyVista, imageio) use
BSD/MIT-compatible licenses. See the README for the full dependency
license table.
