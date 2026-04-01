# MPCO File Format

The `.mpco` file is an HDF5 container that stores the model definition and
all result data from an OpenSees analysis.

## Top-level structure

```
/
├── MODEL_STAGE[1]/
│   ├── MODEL/
│   │   ├── NODES/
│   │   │   ├── ID                    (N,) int
│   │   │   └── COORDINATES           (N, 3) float
│   │   ├── ELEMENTS/
│   │   │   ├── 203-ASDShellQ4[201:0] (M, 5) int  [eid, n1, n2, n3, n4]
│   │   │   ├── 204-ElasticTimoshenkoBeam3d[202:0]  (K, 3) int  [eid, n1, n2]
│   │   │   └── ...
│   │   └── SECTION_ASSIGNMENTS/
│   │       ├── SECTION_0/
│   │       │   ├── ASSIGNMENT         (P, 2) int  [elem_id, fiber_index]
│   │       │   └── FIBER_DATA         (F, 3) float [x, y, weight]
│   │       └── ...
│   └── RESULTS/
│       ├── ON_NODES/
│       │   ├── DISPLACEMENT/
│       │   │   ├── ID                 (N,) int
│       │   │   └── DATA/
│       │   │       ├── STEP_0         (N, 3) float  [Ux, Uy, Uz]
│       │   │       ├── STEP_1
│       │   │       └── ...
│       │   └── REACTION_FORCE/
│       │       └── (same structure)
│       └── ON_ELEMENTS/
│           ├── section.fiber.stress/
│           │   ├── 203-ASDShellQ4[201:0:0]/
│           │   │   ├── ID             (M,) int
│           │   │   ├── META/
│           │   │   │   ├── GAUSS_IDS      (n_gp,) int
│           │   │   │   ├── MULTIPLICITY   (n_gp,) int  [= num_fibers]
│           │   │   │   ├── NUM_COMPONENTS (n_gp,) int
│           │   │   │   └── COMPONENTS     (1,) string
│           │   │   └── DATA/
│           │   │       ├── STEP_0     (M, n_gp * n_fib * n_comp) float
│           │   │       └── ...
│           │   └── 203-ASDShellQ4[201:0:1]/
│           │       └── (same structure, different section)
│           └── section.fiber.strain/
│               └── (same structure)
├── MODEL_STAGE[2]/
│   └── (same structure, different step keys)
└── ...
```

## Element keys

Format: `{geom_id}-{ElementType}[{section_id}:{property_id}]`

Examples:
- `203-ASDShellQ4[201:0]` -- ASDShellQ4 elements with geometry 203
- `204-ElasticTimoshenkoBeam3d[202:0]` -- beam elements with geometry 204

The element type name is extracted with regex: `\d+-(\w+)\[`.

## Element result data layout

Each row in a `DATA/STEP_N` dataset has all Gauss points, fibers, and
components flattened into columns:

```
col = gp * (n_fibers * n_components) + fiber * n_components + component
```

For example, with 4 GPs, 8 fibers, and 5 components (= 160 columns):
- Column 0: GP0, Fiber0, Component0
- Column 1: GP0, Fiber0, Component1
- ...
- Column 5: GP0, Fiber1, Component0
- ...
- Column 40: GP1, Fiber0, Component0
- ...

## META/COMPONENTS string

Format: `0.1.2.3.4.Name0,Name1,Name2,Name3,Name4;...`

One semicolon-separated block per Gauss point (all identical in practice).
Each block has dot-separated component indices followed by comma-separated
names. The parser extracts names from the first block.

When `"Unknown"` appears as a component name, ops-post falls back to
generic labels (`C0`, `C1`, ...).

## What the .mpco does NOT store

The HDF5 file contains node coordinates, element connectivity, section
fiber data, and all result arrays. However, it does **not** store:

- **Element local axes** (shell normal direction, beam orientation)
- **Beam cross-section profiles** (for 3D visualization)

This information is derived from the TCL model files and stored in the
`.mpco.postdata` companion file. See
[Postdata Format](postdata-format.md) for details on how this file is
generated and structured.

## Section thickness detection

Fiber layer y-coordinates encode the through-thickness position:

```python
thickness = (y_max - y_min) / (num_fibers - 1) * num_fibers
```

For a single-fiber section, the thickness is `2 * abs(y)`.
