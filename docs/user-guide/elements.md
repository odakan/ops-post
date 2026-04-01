# Supported Elements

ops-post recognizes 61 OpenSees element types through its element
registry (`ops_elements.py`).

## Shell elements

Displayed as transparent extrusion + mid-surface + fiber layers + GP spheres.

| Element | Nodes | GPs | Notes |
|---------|-------|-----|-------|
| ASDShellQ4 | 4 | 2x2 | General-purpose shell |
| ShellMITC4 | 4 | 2x2 | Mixed interpolation |
| ShellDKGQ | 4 | 2x2 | Discrete Kirchhoff |
| ShellNLDKGQ | 4 | 2x2 | Nonlinear DKG |
| ShellANDeS | 4 | 2x2 | ANDeS formulation |
| shellMITC4Thermal | 4 | 2x2 | Thermal variant |
| shellNLDKGQThermal | 4 | 2x2 | Thermal variant |
| ASDShellT3 | 3 | 3 | Triangular shell |
| ShellDKGT | 3 | 3 | Triangular DKG |
| ShellNLDKGT | 3 | 3 | Triangular nonlinear |
| ShellMITC9 / ShellNL | 9 | 3x3 | 9-node shell |

## 2D continuum elements

Displayed the same as shell elements (surface + extrusion).

| Element | Nodes | GPs |
|---------|-------|-----|
| quad | 4 | 2x2 |
| quad3d | 4 | 2x2 |
| SSPquad | 4 | 2x2 |
| bbarQuad | 4 | 2x2 |
| enhancedQuad | 4 | 2x2 |
| quadWithSensitivity | 4 | 2x2 |
| nineNodeQuad | 9 | 3x3 |
| Tri31 | 3 | 1 |

## Beam / column elements

Displayed as 3D hexahedra using the cross-section profile derived from
TCL section data (A, Iz -> rectangular width and height).

| Element | Nodes |
|---------|-------|
| ElasticTimoshenkoBeam2d / 3d | 2 |
| elasticBeamColumn | 2 |
| forceBeamColumn | 2 |
| forceBeamColumnCBDI / CSBDI | 2 |
| forceBeamColumnThermal | 2 |
| forceBeamColumnWarping | 2 |
| elasticForceBeamColumn | 2 |
| elasticForceBeamColumnWarping | 2 |
| dispBeamColumn | 2 |
| dispBeamColumnInt | 2 |
| dispBeamColumnThermal | 2 |
| dispBeamColumnWithSensitivity | 2 |
| gradientInelasticBeamColumn | 2 |
| internalBeamColumnElement | 2 |
| AxEqDispBeamColumn2d | 2 |
| element2dGNL | 2 |
| MVLEM / MVLEM_3D | 2 |
| SFI_MVLEM / SFI_MVLEM_3D | 2 |
| E_SFI_MVLEM_3D | 2 |

## Truss elements

Displayed as thin beams with a default square cross-section.

| Element | Nodes |
|---------|-------|
| truss / Truss | 2 |
| corotTruss / corotTrussSection | 2 |

## Link and zero-length elements

Recognized but not visualized (no geometry to display).

| Element | Type |
|---------|------|
| twoNodeLink / TwoNodeLink | Link |
| zeroLength / ZeroLength | Zero-length |
| zeroLengthSection | Zero-length |
| zeroLengthND | Zero-length |
| zeroLengthContact2D / 3D | Zero-length |
| zeroLengthContactASDimplex | Zero-length |
| zeroLengthImpact3D | Zero-length |
| zeroLengthInterface2D | Zero-length |
| zeroLengthRocking | Zero-length |
| CoupledZeroLength | Zero-length |

## Unknown elements

If an element type is not in the registry, ops-post falls back to
guessing the topology from the connectivity column count:

| Columns | Assumed type |
|---------|-------------|
| 4 | Shell quad (Q4) |
| 3 | Shell triangle (T3) |
| 9 | 9-node quad (Q9) |
| 2 | Beam |

To add a missing element, see [Extending](../extending.md).
