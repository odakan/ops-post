# Element Registry

The element registry (`ops_elements.py`) maps OpenSees element type names to
their finite element topology. This is used throughout ops-post to
determine how to build meshes, compute Gauss point positions, and
extrapolate results to nodes.

## API

### Lookup

```python
from ops_post.ops_elements import lookup, lookup_or_guess

# Exact match
info = lookup("ASDShellQ4")  # returns ElemInfo or None

# Match with fallback to connectivity guess
info = lookup_or_guess("ASDShellQ4", num_conn_cols=4)
```

### ElemInfo

Each entry in the registry is an `ElemInfo` object:

| Attribute | Type | Description |
|-----------|------|-------------|
| `family` | `ElemFamily` | Element family enum |
| `num_nodes` | `int` | Number of element nodes |
| `num_gp` | `int` | Number of Gauss points |
| `gp_coords` | `ndarray (n_gp, 2)` | GP natural coordinates |
| `shape_fn` | `callable` | Shape functions: `(xi, eta) -> (n_nodes,)` |
| `extrap_fn` | `callable` | Returns extrapolation matrix: `() -> (n_nodes, n_gp)` |
| `face_nodes` | `list[int]` | Node indices for the display face |

Convenience properties:

| Property | True for |
|----------|----------|
| `is_surface` | Shells, 2D quads/tris |
| `is_beam` | Beams, trusses |
| `is_solid` | Bricks, tetrahedra |
| `is_point` | Zero-length, links |

### ElemFamily enum

```python
from ops_post.ops_elements import ElemFamily

ElemFamily.SHELL_QUAD   # 4- or 9-node shell
ElemFamily.SHELL_TRI    # 3-node shell triangle
ElemFamily.QUAD_2D      # 2D continuum quad
ElemFamily.TRI_2D       # 2D continuum triangle
ElemFamily.BRICK        # 8- or 20-node brick
ElemFamily.TET          # 4- or 10-node tetrahedron
ElemFamily.BEAM         # 2-node beam/column
ElemFamily.TRUSS        # 2-node truss
ElemFamily.LINK         # 2-node link
ElemFamily.ZERO_LENGTH  # 2-node zero-length
```

## Shape functions

```python
from ops_post.ops_elements import shape_q4, shape_t3, shape_q9

N = shape_q4(xi, eta)   # (4,) at (xi, eta) in [-1,1]^2
N = shape_t3(xi, eta)   # (3,) at (xi, eta) in area coordinates
N = shape_q9(xi, eta)   # (9,) at (xi, eta) in [-1,1]^2
```

## Gauss point coordinates

Pre-defined arrays of natural coordinates:

| Constant | Points | Domain |
|----------|--------|--------|
| `GP_QUAD_1x1` | 1 | [-1,1]^2 |
| `GP_QUAD_2x2` | 4 | [-1,1]^2 at +/-1/sqrt(3) |
| `GP_QUAD_3x3` | 9 | [-1,1]^2 at +/-sqrt(3/5), 0 |
| `GP_TRI_1` | 1 | Triangle centroid |
| `GP_TRI_3` | 3 | Triangle 3-point rule |
| `GP_BRICK_2x2x2` | 8 | [-1,1]^3 |
| `GP_BRICK_3x3x3` | 27 | [-1,1]^3 |
| `GP_TET_1` | 1 | Tetrahedron centroid |
| `GP_TET_4` | 4 | Tetrahedron 4-point rule |

## Extrapolation matrices

Map Gauss point values to corner node values:

```python
from ops_post.ops_elements import (
    extrapolation_matrix_q4,      # () -> (4, 4)
    extrapolation_matrix_t3_3gp,  # () -> (3, 3)
    extrapolation_matrix_t3_1gp,  # () -> (3, 1)
)
```

Usage:

```python
E = extrapolation_matrix_q4()
nodal_values = gp_values @ E.T  # (n_elems, 4)
```
