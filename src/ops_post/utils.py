"""Shared utility functions."""

import math
import numpy as np
from scipy.spatial.transform import Rotation
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
from vtkmodules.vtkRenderingCore import vtkRenderWindowInteractor


def quaternion_to_rotation_matrix(qw, qx, qy, qz):
    """Convert quaternion (scalar-first) to 3x3 rotation matrix.

    The .cdata file uses scalar-first convention (qw, qx, qy, qz),
    while scipy uses scalar-last (qx, qy, qz, qw).
    """
    r = Rotation.from_quat([qx, qy, qz, qw])
    return r.as_matrix()


def compute_von_mises_shell(s11, s22, s12, s33=None, s13=None):
    """Compute von Mises equivalent stress.

    3 components (plane stress): sqrt(s11² - s11·s22 + s22² + 3·s12²)
    5 components (full shell):   includes s33 and s13 terms
    """
    if s33 is not None and s13 is not None:
        vm = np.sqrt(s11**2 + s22**2 + s33**2
                     - s11*s22 - s22*s33 - s11*s33
                     + 3.0 * (s12**2 + s13**2))
    else:
        vm = np.sqrt(s11**2 - s11 * s22 + s22**2 + 3.0 * s12**2)
    return vm


def compute_magnitude(vectors):
    """Compute Euclidean norm along last axis."""
    return np.linalg.norm(vectors, axis=-1)


def detect_up_axis(coords: np.ndarray) -> int:
    """Return axis index (0=X, 1=Y, 2=Z) with the largest bounding-box extent.

    For a tall wall or frame, the tallest direction is typically 'up'.
    Defaults to Z (2) if extents are within 10% of each other.
    """
    extents = coords.max(axis=0) - coords.min(axis=0)
    max_ext = extents.max()
    if max_ext < 1e-12:
        return 2
    # If all extents are similar (within 10%), default to Z
    if extents.min() / max_ext > 0.9:
        return 2
    return int(np.argmax(extents))


class ArcballInteractorStyle(vtkInteractorStyleTrackballCamera):
    """Custom interaction style:
    - Scroll: zoom (default)
    - Middle mouse: pan/translate
    - Right mouse: arcball rotation around origin (0,0,0)

    The arcball projects the mouse position onto a virtual sphere
    enclosing the model. Dragging rotates that sphere.
    """

    def __init__(self):
        super().__init__()
        self._rotating = False
        self._panning = False
        self._drag_start = None
        self._start_cam_pos = None
        self._start_cam_fp = None
        self._start_cam_up = None

        self.AddObserver("RightButtonPressEvent", self._on_right_press)
        self.AddObserver("RightButtonReleaseEvent", self._on_right_release)
        self.AddObserver("MiddleButtonPressEvent", self._on_middle_press)
        self.AddObserver("MiddleButtonReleaseEvent", self._on_middle_release)
        self.AddObserver("MouseMoveEvent", self._on_mouse_move)
        self.AddObserver("MouseWheelForwardEvent", self._on_scroll_forward)
        self.AddObserver("MouseWheelBackwardEvent", self._on_scroll_backward)
        self.AddObserver("LeftButtonPressEvent", self._on_left_press)
        self.AddObserver("LeftButtonReleaseEvent", self._on_left_release)

    def _on_left_press(self, obj, event):
        # Disable default left-click rotation
        pass

    def _on_left_release(self, obj, event):
        pass

    def _on_right_press(self, obj, event):
        self._rotating = True
        pos = self.GetInteractor().GetEventPosition()
        self._drag_start = pos
        # Snapshot camera state at drag start
        self.FindPokedRenderer(*pos)
        renderer = self.GetCurrentRenderer()
        if renderer:
            cam = renderer.GetActiveCamera()
            self._start_cam_pos = np.array(cam.GetPosition())
            self._start_cam_fp = np.array(cam.GetFocalPoint())
            self._start_cam_up = np.array(cam.GetViewUp())

    def _on_right_release(self, obj, event):
        self._rotating = False
        self._drag_start = None

    def _on_middle_press(self, obj, event):
        self._panning = True
        self.StartPan()

    def _on_middle_release(self, obj, event):
        self._panning = False
        self.EndPan()

    def _zoom(self, factor):
        """Zoom by factor. Handles both perspective (Dolly) and orthographic (ParallelScale)."""
        renderer = self.GetCurrentRenderer()
        if not renderer:
            return
        camera = renderer.GetActiveCamera()
        if camera.GetParallelProjection():
            camera.SetParallelScale(camera.GetParallelScale() / factor)
        else:
            camera.Dolly(factor)
        renderer.ResetCameraClippingRange()
        renderer.GetRenderWindow().Render()

    def _on_scroll_forward(self, obj, event):
        self.FindPokedRenderer(*self.GetInteractor().GetEventPosition())
        self._zoom(1.1)

    def _on_scroll_backward(self, obj, event):
        self.FindPokedRenderer(*self.GetInteractor().GetEventPosition())
        self._zoom(1.0 / 1.1)

    def _on_mouse_move(self, obj, event):
        interactor = self.GetInteractor()
        pos = interactor.GetEventPosition()

        if self._panning:
            self.FindPokedRenderer(*pos)
            self.Pan()
            return

        if self._rotating and self._drag_start is not None:
            self.FindPokedRenderer(*pos)
            self._arcball_rotate(self._drag_start, pos)
            return

    def _arcball_rotate(self, start_pos, cur_pos):
        """True arcball rotation around the world origin (0,0,0).

        Computes a single rotation from drag-start to current mouse position
        and applies it to the saved initial camera state.  No incremental
        composition, so closed mouse paths produce zero net rotation.
        """
        renderer = self.GetCurrentRenderer()
        if renderer is None or self._start_cam_pos is None:
            return

        size = renderer.GetRenderWindow().GetSize()
        w, h = size[0], size[1]
        if w == 0 or h == 0:
            return

        # Project screen coordinates onto a unit sphere
        def to_sphere(sx, sy):
            r = min(w, h)
            x = (2.0 * sx - w) / r
            y = (2.0 * sy - h) / r
            r2 = x * x + y * y
            if r2 <= 1.0:
                z = math.sqrt(1.0 - r2)
            else:
                s = 1.0 / math.sqrt(r2)
                x *= s
                y *= s
                z = 0.0
            return np.array([x, y, z])

        p1 = to_sphere(*start_pos)
        p2 = to_sphere(*cur_pos)

        # Rotation axis and angle in screen space
        axis_scr = np.cross(p1, p2)
        sin_a = np.linalg.norm(axis_scr)
        cos_a = np.dot(p1, p2)
        if sin_a < 1e-10:
            return
        axis_scr /= sin_a
        angle = math.atan2(sin_a, cos_a)

        # Use the INITIAL camera frame for screen-to-world conversion
        pos0 = self._start_cam_pos
        fp0 = self._start_cam_fp
        vup0 = self._start_cam_up

        fwd = fp0 - pos0
        fwd = fwd / np.linalg.norm(fwd)
        right = np.cross(fwd, vup0)
        rlen = np.linalg.norm(right)
        if rlen < 1e-10:
            return
        right /= rlen
        up = np.cross(right, fwd)

        # Screen-to-world: x→right, y→up, z→−fwd
        # Negate angle so dragging right rotates model right (grab-and-pull feel)
        world_axis = axis_scr[0] * right + axis_scr[1] * up - axis_scr[2] * fwd
        rot = Rotation.from_rotvec(-angle * world_axis)

        # Apply rotation to the INITIAL camera state (not current)
        camera = renderer.GetActiveCamera()
        camera.SetPosition(*rot.apply(pos0))
        camera.SetFocalPoint(*rot.apply(fp0))
        camera.SetViewUp(*rot.apply(vup0))
        camera.OrthogonalizeViewUp()
        renderer.ResetCameraClippingRange()
        renderer.GetRenderWindow().Render()
