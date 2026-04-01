# Results

## How results are displayed

ops-post supports two categories of results, each with different
visualization strategies.

### Nodal results

Results stored per node (e.g. displacement, reaction forces).

- Displayed as **smooth contours** on the shell mid-surface mesh.
- Values at shared nodes are mapped directly (no averaging needed).
- Available components depend on the result type:

| Result | Components |
|--------|-----------|
| DISPLACEMENT | Ux, Uy, Uz, \|U\| |
| REACTION_FORCE | Fx, Fy, Fz, \|F\| |

### Element results

Results stored per element at Gauss point / fiber locations.

**Contour mode** (non-fiber results):

- GP values are either averaged across all GPs or extrapolated to
  corner nodes using the element's extrapolation matrix.
- The extrapolated values are averaged at shared nodes and displayed
  as a smooth contour on the mid-surface.
- You select which fiber layer to display.

**GP Spheres mode** (non-fiber results):

- Colored spheres placed at the physical Gauss point location for
  a selected fiber layer.
- Uniform sphere size controlled by the Point size slider.

**Fiber results** (e.g. `section.fiber.stress`):

- All Gauss points and all fiber layers are shown simultaneously.
- Each sphere's **radius** is proportional to its absolute value.
- Each sphere's **color** maps the signed value to the colormap.
- No averaging, no selection -- every integration point shows its
  real value.

## GP sphere zero mode

When displaying GP spheres, the **zero mode** controls how the sphere
radius scale is anchored:

| Mode | Behavior | Use for |
|------|----------|---------|
| **Zero at 0** | Radius is zero when the value is zero; positive and negative values grow outward equally | Stress, strain (signed quantities that can be positive or negative) |
| **Zero at min** | Radius is zero at the minimum value; maximum value gets the largest radius | Damage, energy, or any quantity that ranges from 0 upward |

This setting is in the **View tab** under Figure controls.

## Scale range clamping

Optional **min** and **max** fields in the View tab let you clamp the
color scale:

- Values below the **min** are clamped to the minimum color.
- Values above the **max** are clamped to the maximum color.
- Leave a field empty to use the data's natural range for that end.

This is useful for comparing results across time steps with a fixed
color range, or for filtering out outliers.

## Von Mises stress

Von Mises equivalent stress is available when an element result has
3 or more stress components.

**3 components** (plane stress: S11, S22, S12):

$$\sigma_{vm} = \sqrt{\sigma_{11}^2 - \sigma_{11}\sigma_{22} + \sigma_{22}^2 + 3\tau_{12}^2}$$

**5 components** (full shell: S11, S22, S12, S33, S13):

$$\sigma_{vm} = \sqrt{\sigma_{11}^2 + \sigma_{22}^2 + \sigma_{33}^2 - \sigma_{11}\sigma_{22} - \sigma_{22}\sigma_{33} - \sigma_{11}\sigma_{33} + 3(\tau_{12}^2 + \tau_{13}^2)}$$

## Component naming

Component names are read from the HDF5 `META/COMPONENTS` field. When
`"Unknown"` appears as a component name, ops-post uses generic labels
(`C0`, `C1`, `C2`, ...).

## Displacement scaling

The displacement scale factor applies to all visualization layers:

- Shell mid-surface, extrusion, and fiber edges are displaced using
  nodal displacement values.
- Gauss point spheres are displaced using shape-function interpolation
  at each GP's natural coordinates.
- Beam elements are displaced using their end-node displacements.

A scale factor of 1 shows true displacements. Increase it to
exaggerate deformations for visualization.
