# Getting Started

## Requirements

- Python 3.8 or later
- An OpenSees model (`.tcl` file) and its results (`.mpco` file)

## Installation

### From source

```bash
git clone https://github.com/user/ops-post.git
cd ops-post
pip install .
```

### Editable install (for development)

```bash
pip install -e .
```

This installs all dependencies automatically: h5py, numpy, scipy, pyvista,
pyvistaqt, PyQt5.

## Same-name convention

ops-post expects the TCL model and MPCO results to share a base name:

```
my_model.tcl      <-- OpenSees model (entry point)
my_model.mpco     <-- analysis results (HDF5)
```

The `.tcl` file can `source` other files (nodes.tcl, elements.tcl, etc.)
-- ops-post follows all `source` commands automatically. If your model
uses a multi-file layout with `main.tcl` as the entry point, rename or
create a single-file version with the same base name as the `.mpco`.

!!! tip
    Use the included `merge_tcl.py` utility to merge a multi-file model
    into one: `python merge_tcl.py main.tcl my_model.tcl`

## Running

```bash
# Just the base name -- finds .tcl and .mpco automatically
ops-post my_model

# Or with any extension -- the base name is extracted
ops-post my_model.mpco
ops-post my_model.tcl

# File dialog (no arguments)
ops-post
```

## What happens

```
ops-post my_model
        |
        v
  Find my_model.tcl ──► Parse TCL ──► Generate my_model.mpco.postdata
        |                                  (local axes, beam profiles,
        |                                   element info)
        v
  Find my_model.mpco ──► Read HDF5 ──► Open 3D viewer
                              |
                              v
                    Read .postdata ──► Local axes, beam profiles
```

1. **First run**: parses the TCL, generates `.postdata`, opens the viewer.
2. **Subsequent runs**: skips parsing (`.postdata` is up to date), opens
   directly.
3. **After editing TCL**: `.postdata` is regenerated automatically because
   the `.tcl` file is newer.

## The .postdata file

Generated automatically from your TCL model. Contains:

| Section | Content |
|---------|---------|
| `*LOCAL_AXES` | Quaternion orientation per element |
| `*BEAM_PROFILE` | Rectangular cross-section vertices per beam section |
| `*BEAM_PROFILE_ASSIGNMENT` | Element-to-profile mapping |
| `*ELEMENT_INFO` | Element type and section name per element |

You never need to create or edit this file manually.

!!! note
    The `.mpco` file can be large (50-500 MB for typical models). ops-post
    uses lazy step loading and a persistent HDF5 handle, so only one time
    step is in memory at a time.

## First steps

1. Run `ops-post my_model` -- the viewer opens with the default view.
2. Use **right-click drag** to rotate, **middle-click drag** to pan,
   **scroll** to zoom.
3. Select a result category (Nodal / Element) and component from the
   right panel.
4. Use the time step slider or Play button to animate.

See [Controls & Display](user-guide/controls.md) for the full reference.
