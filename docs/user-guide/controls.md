# Controls & Display

## Mouse controls

| Input | Action |
|-------|--------|
| Right-click drag | Arcball rotation around world origin |
| Middle-click drag | Pan (translate view) |
| Scroll up/down | Zoom in/out (perspective and orthographic) |

The rotation uses a true Shoemake arcball -- drawing circles with the mouse
returns to exactly the starting orientation (no drift).

## Tabbed control panel

The right dock panel is organized into two tabs: **Results** and **View**.

---

## Results tab

### Result Selection

- **Category**: `Nodal` or `Element`
- **Result**: e.g. `DISPLACEMENT`, `section.fiber.stress`
- **Component**: e.g. `Ux`, `C0`, `von Mises`

When component names in the HDF5 file are `"Unknown"`, ops-post falls
back to generic labels (`C0`, `C1`, `C2`, ...).

### Element Display Mode

Only visible when an Element result is selected.

For **non-fiber** results (e.g. `section.stress`):

- **Contour (averaged)**: extrapolates GP values to nodes, displays as
  a smooth contour on the mid-surface.
- **GP Spheres (fiber-level)**: shows colored spheres at Gauss point
  locations. You can select a specific GP and fiber layer.

For **fiber results** (e.g. `section.fiber.stress`):

- The mode selector, GP dropdown, and fiber dropdown are hidden.
- All GP x fiber points are shown simultaneously as scaled spheres.
- Sphere radius is proportional to the absolute value of the result.
- The **Point size** slider controls the maximum sphere radius.

### Display Options

| Control | Effect |
|---------|--------|
| Disp. Scale | Displacement magnification factor (spinbox + slider) |
| Show Shells | Toggle shell element visibility |
| Show Beams | Toggle beam element visibility |
| Show Edges | Toggle element edge lines |
| Show Extrusion | Toggle transparent thickness extrusion |
| Show Fiber Layers | Toggle fiber layer edge lines (auto-hidden when GP spheres are active) |
| Grid | Toggle floor grid visibility |
| Grid spacing | Floor grid line spacing (in model units) |

!!! note
    Fiber layer edge lines are automatically hidden when GP spheres are
    active to avoid visual clutter.

### Time Step

| Control | Effect |
|---------|--------|
| Step slider | Jump to any analysis step |
| \|< | First step |
| < | Previous step |
| Play/Pause | Animate through all steps |
| > | Next step |
| >\| | Last step |
| FPS | Animation speed (frames per second) |

The **Min/Max** label at the bottom shows the range of the currently
displayed scalar field.

---

## View tab

### Camera

| Control | Effect |
|---------|--------|
| Front / Back / Left / Right / Top / Bottom / Iso | Standard camera presets |
| Orthographic / Perspective | Projection mode toggle |

### Figure

Controls for publication-ready output. All settings apply immediately --
no Apply button needed.

| Control | Effect |
|---------|--------|
| Colormap | Color scheme: jet, coolwarm, viridis, plasma, turbo, RdYlBu_r, hot_r, YlOrRd, inferno_r |
| Viewport W / H | Viewport size in pixels (default 1200 x 800) |
| Aspect presets | 16:9, 4:3, 3:2, 1:1, 2:1 -- adjusts H to match W |
| Scale bar orientation | Vertical or Horizontal |
| Scale bar position | Top-Left, Top-Right, Bottom-Left, Bottom-Right |
| Scale bar W / H | Width and height of the scale bar (fraction of viewport) |
| Scale bar font | Font size for scale bar labels |
| Title override | Custom title text for the scale bar (overrides auto-generated) |
| Scale range min | Optional minimum value -- values below are clamped |
| Scale range max | Optional maximum value -- values above are clamped |
| GP zero mode | "Zero at 0" (stress: +/- around zero) or "Zero at min" (damage: 0-to-max) |

!!! tip
    The viewport size is fixed (not tied to window size). Resize it only
    through the W/H spinners in the View tab to get consistent figure
    dimensions for publications.

### Export

| Control | Effect |
|---------|--------|
| Save Figure | Save current view as PNG (1x/2x/3x/4x resolution), SVG vector, or PDF vector |
| Record Animation | Choose output file (GIF or MP4) first, then auto-plays from current step, auto-stops and saves at the last step |

When recording an animation:

1. Click **Record Animation** and choose a filename (`.gif` or `.mp4`).
2. Playback starts automatically from the current step.
3. Recording stops and saves automatically at the last step.

MP4 uses bundled ffmpeg via imageio. If ffmpeg is unavailable, the
recording falls back to GIF format.

The animation uses the current FPS setting for playback speed.
