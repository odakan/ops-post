# Postdata File Format

The `.mpco.postdata` file is a plain-text companion to the `.mpco` result
file. It is **generated automatically** by parsing the OpenSees TCL model
files and contains geometry metadata that the HDF5 file does not store.

## Generation

```bash
# Automatic (generated on first run if .tcl exists)
ops-post my_model

# The parser follows all source commands in the TCL:
# my_model.tcl -> source nodes.tcl -> source elements.tcl -> ...
```

The `.postdata` file is regenerated when the `.tcl` file is newer.

## File structure

```
#elem_id q.w q.x q.y q.z
*LOCAL_AXES
1 0.0 0.7071 0.0 -0.7071
2 0.0 0.7071 0.0 -0.7071
...
#profile_id n_vertices\n  y z (per vertex)
*BEAM_PROFILE
-1 4
-0.0285 -0.0285
0.0285 -0.0285
0.0285 0.0285
-0.0285 0.0285
...
#elem_id profile_id
*BEAM_PROFILE_ASSIGNMENT
1 -1
2 -2
...
#elem_id elem_type section_name
*ELEMENT_INFO
1 ElasticTimoshenkoBeam
239 ASDShellQ4 LayeredShell_10
...
```

## Sections

### *LOCAL_AXES

One line per element. Quaternion (scalar-first: qw, qx, qy, qz) defining
the element's local coordinate frame.

```
elem_id  qw  qx  qy  qz
```

**How local axes are computed from TCL:**

| Element family | Source |
|----------------|--------|
| Shell/quad/tri with `-local vx vy vz` | Local-x direction + element normal |
| Beam with `geomTransf ... vecxzX vecxzY vecxzZ` | Element axis + vecxz plane vector |
| Zero-length with `-orient x1 x2 x3 y1 y2 y3` | Explicit local-x and local-y |
| Other | Computed from node positions (fallback) |

### *BEAM_PROFILE

One block per beam cross-section profile. Each block starts with the
profile ID and vertex count, followed by 2D vertices (y, z in local
coordinates).

```
profile_id  n_vertices
y0  z0
y1  z1
...
```

**How profiles are derived from TCL:**

For beams with inline section properties (`ElasticTimoshenkoBeam`,
`elasticBeamColumn`), the cross-section dimensions are computed from
the area and moment of inertia assuming a rectangular section:

```
h = sqrt(12 * Iz / A)
b = A / h
```

For beams referencing a `section Elastic`, the same formula is applied
using the section's `A` and `Iz` values.

Profile IDs are negative (e.g. `-1`, `-2`) for inline beam properties
to avoid collision with TCL section IDs.

### *BEAM_PROFILE_ASSIGNMENT

Maps each beam element to its cross-section profile.

```
elem_id  profile_id
```

### *ELEMENT_INFO

One line per element with its type and section name.

```
elem_id  elem_type  section_name
```

The section name is derived from the TCL section definition
(e.g. `LayeredShell_10` for `section LayeredShell 10 ...`).

## Relationship to .mpco

The `.mpco` file stores:

- Node coordinates and element connectivity
- Section assignments and fiber data
- All result data (displacements, stresses, strains)

The `.postdata` file adds:

- Local coordinate frame orientations (needed for shell normal direction
  and beam cross-section orientation)
- Beam cross-section profiles (needed for 3D beam visualization)
- Element metadata linking elements to their TCL definitions

Together they provide everything needed for visualization.
