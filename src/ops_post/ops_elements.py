"""OpenSees element type registry.

Maps OpenSees element names to their topology: shape family,
number of nodes, number of Gauss points, and natural GP coordinates.
"""

import numpy as np
from enum import Enum


class ElemFamily(Enum):
    SHELL_QUAD = "shell_quad"      # 4-node shell (Q4) or 9-node (Q9)
    SHELL_TRI = "shell_tri"        # 3-node shell triangle
    QUAD_2D = "quad_2d"            # 2D continuum quad
    TRI_2D = "tri_2d"              # 2D continuum triangle
    BRICK = "brick"                # 8- or 20-node brick
    TET = "tet"                    # 4- or 10-node tetrahedron
    BEAM = "beam"                  # 2-node beam/column
    TRUSS = "truss"                # 2-node truss
    LINK = "link"                  # 2-node link
    ZERO_LENGTH = "zero_length"    # 2-node zero-length
    UNKNOWN = "unknown"


# ── Standard Gauss point coordinates in natural space ───────────────

_g1 = 1.0 / np.sqrt(3.0)

# Quad/hex natural coordinates
GP_QUAD_1x1 = np.array([[0.0, 0.0]])

GP_QUAD_2x2 = np.array([
    [-_g1, -_g1],
    [+_g1, -_g1],
    [+_g1, +_g1],
    [-_g1, +_g1],
])

_g3 = np.sqrt(3.0 / 5.0)
GP_QUAD_3x3 = np.array([
    [-_g3, -_g3], [0.0, -_g3], [+_g3, -_g3],
    [-_g3,  0.0], [0.0,  0.0], [+_g3,  0.0],
    [-_g3, +_g3], [0.0, +_g3], [+_g3, +_g3],
])

# Triangle natural coordinates (area coordinates mapped to xi, eta)
GP_TRI_1 = np.array([[1.0 / 3.0, 1.0 / 3.0]])

GP_TRI_3 = np.array([
    [1.0 / 6.0, 1.0 / 6.0],
    [2.0 / 3.0, 1.0 / 6.0],
    [1.0 / 6.0, 2.0 / 3.0],
])

# Brick natural coordinates
GP_BRICK_2x2x2 = np.array([
    [s0 * _g1, s1 * _g1, s2 * _g1]
    for s2 in [-1, 1] for s1 in [-1, 1] for s0 in [-1, 1]
])

GP_BRICK_3x3x3 = np.array([
    [x, y, z]
    for z in [-_g3, 0, _g3] for y in [-_g3, 0, _g3] for x in [-_g3, 0, _g3]
])

# Tet natural coordinates
GP_TET_1 = np.array([[0.25, 0.25, 0.25]])

GP_TET_4 = np.array([
    [0.1381966, 0.1381966, 0.1381966],
    [0.5854102, 0.1381966, 0.1381966],
    [0.1381966, 0.5854102, 0.1381966],
    [0.1381966, 0.1381966, 0.5854102],
])


# ── Shape functions ─────────────────────────────────────────────────

def shape_q4(xi, eta):
    """4-node quad shape functions at (xi, eta) in [-1,1]²."""
    return 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ])


def shape_t3(xi, eta):
    """3-node triangle shape functions at (xi, eta) in area coords."""
    return np.array([1.0 - xi - eta, xi, eta])


def shape_q9(xi, eta):
    """9-node quad (serendipity) shape functions at (xi, eta) in [-1,1]²."""
    N = np.zeros(9)
    # Corner nodes
    N[0] = 0.25 * xi * (xi - 1) * eta * (eta - 1)
    N[1] = 0.25 * xi * (xi + 1) * eta * (eta - 1)
    N[2] = 0.25 * xi * (xi + 1) * eta * (eta + 1)
    N[3] = 0.25 * xi * (xi - 1) * eta * (eta + 1)
    # Mid-side nodes
    N[4] = 0.5 * (1 - xi**2) * eta * (eta - 1)
    N[5] = 0.5 * xi * (xi + 1) * (1 - eta**2)
    N[6] = 0.5 * (1 - xi**2) * eta * (eta + 1)
    N[7] = 0.5 * xi * (xi - 1) * (1 - eta**2)
    # Center node
    N[8] = (1 - xi**2) * (1 - eta**2)
    return N


# ── Gauss-to-node extrapolation matrices ────────────────────────────

def extrapolation_matrix_q4():
    """4 GPs (2×2) → 4 corner nodes for Q4."""
    s = np.sqrt(3.0) / 2.0
    a, b, c = 1.0 + s, -0.5, 1.0 - s
    return np.array([
        [a, b, c, b],
        [b, a, b, c],
        [c, b, a, b],
        [b, c, b, a],
    ])


def extrapolation_matrix_t3_3gp():
    """3 GPs → 3 corner nodes for T3 (linear extrapolation)."""
    return np.array([
        [ 5.0/3, -1.0/3, -1.0/3],
        [-1.0/3,  5.0/3, -1.0/3],
        [-1.0/3, -1.0/3,  5.0/3],
    ])


def extrapolation_matrix_t3_1gp():
    """1 GP → 3 corner nodes for T3 (constant)."""
    return np.ones((3, 1))


# ── Element registry ────────────────────────────────────────────────

class ElemInfo:
    """Topology descriptor for one element type."""
    __slots__ = ("family", "num_nodes", "num_gp", "gp_coords",
                 "shape_fn", "extrap_fn", "face_nodes")

    def __init__(self, family, num_nodes, num_gp, gp_coords,
                 shape_fn, extrap_fn, face_nodes):
        self.family = family
        self.num_nodes = num_nodes
        self.num_gp = num_gp
        self.gp_coords = gp_coords      # (num_gp, 2 or 3) natural coords
        self.shape_fn = shape_fn         # callable(xi, eta) → (num_nodes,)
        self.extrap_fn = extrap_fn       # callable() → (num_nodes, num_gp) matrix
        self.face_nodes = face_nodes     # node indices forming the display face

    @property
    def is_surface(self):
        """True for shells, 2D quads/tris — elements with a displayable surface."""
        return self.family in (
            ElemFamily.SHELL_QUAD, ElemFamily.SHELL_TRI,
            ElemFamily.QUAD_2D, ElemFamily.TRI_2D,
        )

    @property
    def is_beam(self):
        return self.family in (ElemFamily.BEAM, ElemFamily.TRUSS)

    @property
    def is_solid(self):
        return self.family in (ElemFamily.BRICK, ElemFamily.TET)

    @property
    def is_point(self):
        return self.family in (ElemFamily.ZERO_LENGTH, ElemFamily.LINK)


# ── The registry ────────────────────────────────────────────────────

_Q4 = ElemInfo(
    family=ElemFamily.SHELL_QUAD, num_nodes=4, num_gp=4,
    gp_coords=GP_QUAD_2x2, shape_fn=shape_q4,
    extrap_fn=extrapolation_matrix_q4,
    face_nodes=[0, 1, 2, 3],
)

_T3 = ElemInfo(
    family=ElemFamily.SHELL_TRI, num_nodes=3, num_gp=3,
    gp_coords=GP_TRI_3, shape_fn=shape_t3,
    extrap_fn=extrapolation_matrix_t3_3gp,
    face_nodes=[0, 1, 2],
)

_T3_1GP = ElemInfo(
    family=ElemFamily.SHELL_TRI, num_nodes=3, num_gp=1,
    gp_coords=GP_TRI_1, shape_fn=shape_t3,
    extrap_fn=extrapolation_matrix_t3_1gp,
    face_nodes=[0, 1, 2],
)

_Q9 = ElemInfo(
    family=ElemFamily.SHELL_QUAD, num_nodes=9, num_gp=9,
    gp_coords=GP_QUAD_3x3, shape_fn=shape_q9,
    extrap_fn=None,  # 9→9 identity-like, not commonly extrapolated
    face_nodes=[0, 1, 2, 3, 4, 5, 6, 7, 8],
)

_Q4_2D = ElemInfo(
    family=ElemFamily.QUAD_2D, num_nodes=4, num_gp=4,
    gp_coords=GP_QUAD_2x2, shape_fn=shape_q4,
    extrap_fn=extrapolation_matrix_q4,
    face_nodes=[0, 1, 2, 3],
)

_Q9_2D = ElemInfo(
    family=ElemFamily.QUAD_2D, num_nodes=9, num_gp=9,
    gp_coords=GP_QUAD_3x3, shape_fn=shape_q9,
    extrap_fn=None,
    face_nodes=[0, 1, 2, 3, 4, 5, 6, 7, 8],
)

_T3_2D = ElemInfo(
    family=ElemFamily.TRI_2D, num_nodes=3, num_gp=1,
    gp_coords=GP_TRI_1, shape_fn=shape_t3,
    extrap_fn=extrapolation_matrix_t3_1gp,
    face_nodes=[0, 1, 2],
)

_BEAM = ElemInfo(
    family=ElemFamily.BEAM, num_nodes=2, num_gp=0,
    gp_coords=None, shape_fn=None,
    extrap_fn=None, face_nodes=[],
)

_TRUSS = ElemInfo(
    family=ElemFamily.TRUSS, num_nodes=2, num_gp=0,
    gp_coords=None, shape_fn=None,
    extrap_fn=None, face_nodes=[],
)

_LINK = ElemInfo(
    family=ElemFamily.LINK, num_nodes=2, num_gp=0,
    gp_coords=None, shape_fn=None,
    extrap_fn=None, face_nodes=[],
)

_ZERO = ElemInfo(
    family=ElemFamily.ZERO_LENGTH, num_nodes=2, num_gp=0,
    gp_coords=None, shape_fn=None,
    extrap_fn=None, face_nodes=[],
)


# Name → ElemInfo.  Keys are the element type strings extracted from
# the HDF5 key (e.g. "ASDShellQ4" from "203-ASDShellQ4[201:0]").
ELEMENT_REGISTRY: dict[str, ElemInfo] = {
    # ── Shells ──
    "ASDShellQ4":       _Q4,
    "ShellMITC4":       _Q4,
    "ShellDKGQ":        _Q4,
    "ShellNLDKGQ":      _Q4,
    "ShellANDeS":       _Q4,
    "shellMITC4Thermal": _Q4,
    "shellNLDKGQThermal": _Q4,
    "ASDShellT3":       _T3,
    "ShellDKGT":        _T3,
    "ShellNLDKGT":      _T3,
    "ShellMITC9":       _Q9,
    "ShellNL":          _Q9,
    # ── 2D Continuum ──
    "quad":             _Q4_2D,
    "quad3d":           _Q4_2D,
    "SSPquad":          _Q4_2D,
    "bbarQuad":         _Q4_2D,
    "enhancedQuad":     _Q4_2D,
    "quadWithSensitivity": _Q4_2D,
    "nineNodeQuad":     _Q9_2D,
    "Tri31":            _T3_2D,
    # ── Beams / Columns ──
    "ElasticTimoshenkoBeam":     _BEAM,
    "ElasticTimoshenkoBeam3d":   _BEAM,
    "ElasticTimoshenkoBeam2d":   _BEAM,
    "elasticBeamColumn":         _BEAM,
    "forceBeamColumn":           _BEAM,
    "forceBeamColumnCBDI":       _BEAM,
    "forceBeamColumnCSBDI":      _BEAM,
    "forceBeamColumnThermal":    _BEAM,
    "forceBeamColumnWarping":    _BEAM,
    "elasticForceBeamColumn":    _BEAM,
    "elasticForceBeamColumnWarping": _BEAM,
    "dispBeamColumn":            _BEAM,
    "dispBeamColumnInt":         _BEAM,
    "dispBeamColumnThermal":     _BEAM,
    "dispBeamColumnWithSensitivity": _BEAM,
    "gradientInelasticBeamColumn": _BEAM,
    "internalBeamColumnElement": _BEAM,
    "AxEqDispBeamColumn2d":      _BEAM,
    "element2dGNL":              _BEAM,
    "MVLEM":                     _BEAM,
    "MVLEM_3D":                  _BEAM,
    "SFI_MVLEM":                 _BEAM,
    "SFI_MVLEM_3D":              _BEAM,
    "E_SFI_MVLEM_3D":            _BEAM,
    # ── Truss ──
    "truss":            _TRUSS,
    "Truss":            _TRUSS,
    "corotTruss":       _TRUSS,
    "corotTrussSection": _TRUSS,
    # ── Link ──
    "twoNodeLink":      _LINK,
    "TwoNodeLink":      _LINK,
    # ── Zero-length ──
    "zeroLength":       _ZERO,
    "ZeroLength":       _ZERO,
    "zeroLengthSection": _ZERO,
    "zeroLengthND":     _ZERO,
    "zeroLengthContact2D": _ZERO,
    "zeroLengthContact3D": _ZERO,
    "zeroLengthContactASDimplex": _ZERO,
    "zeroLengthImpact3D": _ZERO,
    "zeroLengthInterface2D": _ZERO,
    "zeroLengthRocking": _ZERO,
    "CoupledZeroLength": _ZERO,
}


def lookup(elem_type: str) -> ElemInfo:
    """Look up element info by OpenSees element type name.

    Falls back to guessing from connectivity column count if not in registry.
    """
    return ELEMENT_REGISTRY.get(elem_type)


def lookup_or_guess(elem_type: str, num_conn_cols: int) -> ElemInfo:
    """Look up element info, or guess from connectivity width if unknown.

    Args:
        elem_type: element type string from HDF5 key
        num_conn_cols: number of node-id columns in connectivity array
    """
    info = ELEMENT_REGISTRY.get(elem_type)
    if info is not None:
        return info

    # Fallback: guess by node count
    if num_conn_cols == 4:
        return _Q4
    if num_conn_cols == 3:
        return _T3
    if num_conn_cols == 9:
        return _Q9
    if num_conn_cols == 2:
        return _BEAM
    return None
