"""PyQt5 main window with PyVista 3D viewport and control panel."""

import os
import numpy as np
import pyvista as pv
import h5py

from PyQt5.QtWidgets import (
    QMainWindow, QDockWidget, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QComboBox, QSlider, QPushButton, QCheckBox, QGroupBox,
    QSpinBox, QDoubleSpinBox, QSizePolicy, QFileDialog, QTabWidget,
    QLineEdit,
)
from PyQt5.QtCore import Qt, QTimer
from pyvistaqt import QtInteractor

from .mpco_reader import read_mpco
from .mesh_builder import MeshBuilder
from .result_processor import ResultProcessor
from .utils import ArcballInteractorStyle


class MainWindow(QMainWindow):

    def __init__(self, mpco_path: str, parent=None):
        super().__init__(parent)
        self.mpco_path = os.path.abspath(mpco_path)
        self.setWindowTitle(f"MPCO Post-Processor — {os.path.basename(mpco_path)}")

        # Load model
        self.model = read_mpco(self.mpco_path)
        self.builder = MeshBuilder(self.model)
        self.processor = ResultProcessor(self.model, self.builder)

        # Keep HDF5 file open for fast step reads during animation
        self._h5file = h5py.File(self.mpco_path, "r")

        # Build base meshes
        self.shell_surface_base = self.builder.build_shell_surface_mesh()
        self.shell_extrusion_base = self.builder.build_shell_extrusion_mesh()
        self.fiber_layer_edges_base = self.builder.build_fiber_layer_edges()
        self.beam_mesh_base = self.builder.build_beam_mesh()

        # Precompute GP cloud for all fibers (base positions)
        gp_data = self.builder.build_gauss_point_cloud(fiber_idx=None)
        self._gp_points_base = gp_data[0]
        self._gp_elem_ids = gp_data[1]
        self._gp_indices = gp_data[2]
        self._gp_fiber_indices = gp_data[3]

        # Build step list
        self._build_step_list()

        # State
        self._current_step_idx = 0
        self._scale_factor = 1.0
        self._playing = False
        self._show_shells = True
        self._show_beams = True
        self._show_edges = True
        self._show_extrusion = True
        self._show_fibers = True
        self._sphere_radius = 0.01

        # Timer
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._on_timer_tick)

        # Model geometric center — used for initial camera framing
        coords = self.model.nodes.coords
        self._model_center = (coords.min(axis=0) + coords.max(axis=0)) / 2.0

        # Floor grid settings — adaptive default spacing
        self._show_grid = True
        max_span = max(self.model.nodes.coords.max(axis=0) - self.model.nodes.coords.min(axis=0))
        self._grid_spacing = self._auto_grid_spacing(max_span)
        self._rebuild_floor_grid()

        # Actor handles for in-place updates (set during full rebuild)
        self._actors = {}   # name -> actor
        self._meshes = {}   # name -> live mesh reference
        self._needs_full_rebuild = True

        # Setup UI
        self._setup_plotter()
        self._setup_control_panel()
        self._populate_dropdowns()
        self._update_display()

    def closeEvent(self, event):
        """Close HDF5 handle when window is closed."""
        if self._h5file:
            self._h5file.close()
            self._h5file = None
        super().closeEvent(event)

    def _build_step_list(self):
        self._steps = []
        for stage in self.model.stages:
            for sk in stage.step_keys:
                self._steps.append((stage.stage_index, sk))

    def _setup_plotter(self):
        self.plotter_widget = QtInteractor(self)
        self.setCentralWidget(self.plotter_widget)
        self.plotter_widget.setFixedSize(1200, 800)
        self.plotter_widget.set_background("white")
        self.plotter_widget.add_axes()

        # Set custom arcball interaction: right-drag=rotate, middle=pan, scroll=zoom
        style = ArcballInteractorStyle()
        self.plotter_widget.interactor.SetInteractorStyle(style)

    def _setup_control_panel(self):
        dock = QDockWidget("Controls", self)
        dock.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetFloatable)
        dock.setMinimumWidth(300)

        tabs = QTabWidget()

        # =====================================================================
        # TAB 1: Results
        # =====================================================================
        results_tab = QWidget()
        rl_top = QVBoxLayout(results_tab)
        rl_top.setSpacing(6)

        # --- Result Selection ---
        result_group = QGroupBox("Result Selection")
        rl = QVBoxLayout()
        rl.addWidget(QLabel("Category:"))
        self.combo_category = QComboBox()
        rl.addWidget(self.combo_category)
        rl.addWidget(QLabel("Result:"))
        self.combo_result = QComboBox()
        rl.addWidget(self.combo_result)
        rl.addWidget(QLabel("Component:"))
        self.combo_component = QComboBox()
        rl.addWidget(self.combo_component)
        result_group.setLayout(rl)
        rl_top.addWidget(result_group)

        # --- Element Display Mode ---
        self.elem_mode_group = QGroupBox("Element Display Mode")
        eml = QVBoxLayout()
        self.combo_elem_mode = QComboBox()
        self.combo_elem_mode.addItems(["Contour (averaged)", "GP Spheres (fiber-level)"])
        eml.addWidget(self.combo_elem_mode)
        self.lbl_gp = QLabel("Gauss Point:")
        eml.addWidget(self.lbl_gp)
        self.combo_gp = QComboBox()
        eml.addWidget(self.combo_gp)
        self.lbl_fiber = QLabel("Fiber:")
        eml.addWidget(self.lbl_fiber)
        self.combo_fiber = QComboBox()
        eml.addWidget(self.combo_fiber)
        sph_row = QHBoxLayout()
        sph_row.addWidget(QLabel("Point size:"))
        self.spin_sphere = QDoubleSpinBox()
        self.spin_sphere.setRange(1, 50)
        self.spin_sphere.setValue(10)
        self.spin_sphere.setSingleStep(2)
        self.spin_sphere.setDecimals(0)
        sph_row.addWidget(self.spin_sphere)
        eml.addLayout(sph_row)
        self.elem_mode_group.setLayout(eml)
        rl_top.addWidget(self.elem_mode_group)

        # --- Display Options ---
        display_group = QGroupBox("Display")
        dl = QVBoxLayout()
        scale_row = QHBoxLayout()
        scale_row.addWidget(QLabel("Disp. Scale:"))
        self.spin_scale = QSpinBox()
        self.spin_scale.setRange(0, 100000)
        self.spin_scale.setValue(1)
        scale_row.addWidget(self.spin_scale)
        dl.addLayout(scale_row)
        self.slider_scale = QSlider(Qt.Horizontal)
        self.slider_scale.setRange(0, 1000)
        self.slider_scale.setValue(10)
        dl.addWidget(self.slider_scale)
        self.chk_shells = QCheckBox("Show Shells")
        self.chk_shells.setChecked(True)
        dl.addWidget(self.chk_shells)
        self.chk_beams = QCheckBox("Show Beams")
        self.chk_beams.setChecked(True)
        dl.addWidget(self.chk_beams)
        self.chk_edges = QCheckBox("Show Edges")
        self.chk_edges.setChecked(True)
        dl.addWidget(self.chk_edges)
        self.chk_extrusion = QCheckBox("Show Extrusion")
        self.chk_extrusion.setChecked(True)
        dl.addWidget(self.chk_extrusion)
        self.chk_fiber_lines = QCheckBox("Show Fiber Layers")
        self.chk_fiber_lines.setChecked(True)
        dl.addWidget(self.chk_fiber_lines)
        self.chk_grid = QCheckBox("Show Floor Grid")
        self.chk_grid.setChecked(True)
        dl.addWidget(self.chk_grid)
        grid_row = QHBoxLayout()
        grid_row.addWidget(QLabel("Grid spacing:"))
        self.spin_grid = QDoubleSpinBox()
        self.spin_grid.setRange(0.01, 100.0)
        self.spin_grid.setValue(self._grid_spacing)
        self.spin_grid.setSingleStep(0.1)
        self.spin_grid.setDecimals(2)
        grid_row.addWidget(self.spin_grid)
        dl.addLayout(grid_row)
        display_group.setLayout(dl)
        rl_top.addWidget(display_group)

        # --- Time Step ---
        time_group = QGroupBox("Time Step")
        tl = QVBoxLayout()
        self.slider_step = QSlider(Qt.Horizontal)
        self.slider_step.setRange(0, max(len(self._steps) - 1, 0))
        self.slider_step.setValue(0)
        tl.addWidget(self.slider_step)
        self.lbl_step_info = QLabel("Frame 1 / 0")
        tl.addWidget(self.lbl_step_info)
        self.lbl_stage_info = QLabel("")
        tl.addWidget(self.lbl_stage_info)
        btn_row = QHBoxLayout()
        self.btn_first = QPushButton("|<")
        self.btn_back = QPushButton("<")
        self.btn_play = QPushButton("Play")
        self.btn_next = QPushButton(">")
        self.btn_last = QPushButton(">|")
        for btn in [self.btn_first, self.btn_back, self.btn_play, self.btn_next, self.btn_last]:
            btn.setMaximumWidth(50)
            btn_row.addWidget(btn)
        tl.addLayout(btn_row)
        speed_row = QHBoxLayout()
        speed_row.addWidget(QLabel("FPS:"))
        self.spin_fps = QSpinBox()
        self.spin_fps.setRange(1, 30)
        self.spin_fps.setValue(5)
        speed_row.addWidget(self.spin_fps)
        tl.addLayout(speed_row)
        time_group.setLayout(tl)
        rl_top.addWidget(time_group)

        # Min/Max label
        self.lbl_min_max = QLabel("")
        rl_top.addWidget(self.lbl_min_max)

        rl_top.addStretch()
        tabs.addTab(results_tab, "Results")

        # =====================================================================
        # TAB 2: View
        # =====================================================================
        view_tab = QWidget()
        vl_top = QVBoxLayout(view_tab)
        vl_top.setSpacing(6)

        # --- Camera ---
        camera_group = QGroupBox("Camera")
        cl = QVBoxLayout()
        view_row1 = QHBoxLayout()
        self.btn_front = QPushButton("Front")
        self.btn_back_view = QPushButton("Back")
        self.btn_left = QPushButton("Left")
        view_row1.addWidget(self.btn_front)
        view_row1.addWidget(self.btn_back_view)
        view_row1.addWidget(self.btn_left)
        cl.addLayout(view_row1)
        view_row2 = QHBoxLayout()
        self.btn_right = QPushButton("Right")
        self.btn_top = QPushButton("Top")
        self.btn_bottom = QPushButton("Bottom")
        view_row2.addWidget(self.btn_right)
        view_row2.addWidget(self.btn_top)
        view_row2.addWidget(self.btn_bottom)
        cl.addLayout(view_row2)
        view_row3 = QHBoxLayout()
        self.btn_iso = QPushButton("Iso")
        self.btn_persp = QPushButton("Orthographic")
        self.btn_persp.setCheckable(True)
        self.btn_persp.setChecked(False)
        view_row3.addWidget(self.btn_iso)
        view_row3.addWidget(self.btn_persp)
        cl.addLayout(view_row3)
        camera_group.setLayout(cl)
        vl_top.addWidget(camera_group)

        # --- Figure ---
        figure_group = QGroupBox("Figure")
        fl = QVBoxLayout()

        # Colormap
        cmap_row = QHBoxLayout()
        cmap_row.addWidget(QLabel("Colormap:"))
        self.combo_cmap = QComboBox()
        self.combo_cmap.addItems([
            "jet", "coolwarm", "viridis", "plasma", "turbo",
            "RdYlBu_r", "hot_r", "YlOrRd", "inferno_r",
        ])
        cmap_row.addWidget(self.combo_cmap)
        fl.addLayout(cmap_row)

        # Viewport size
        size_row = QHBoxLayout()
        size_row.addWidget(QLabel("Size:"))
        self.spin_fig_w = QSpinBox()
        self.spin_fig_w.setRange(200, 4000)
        self.spin_fig_w.setValue(1200)
        self.spin_fig_w.setSuffix(" px")
        size_row.addWidget(self.spin_fig_w)
        size_row.addWidget(QLabel("x"))
        self.spin_fig_h = QSpinBox()
        self.spin_fig_h.setRange(200, 4000)
        self.spin_fig_h.setValue(800)
        self.spin_fig_h.setSuffix(" px")
        size_row.addWidget(self.spin_fig_h)
        fl.addLayout(size_row)
        ar_row = QHBoxLayout()
        ar_row.addWidget(QLabel("Preset:"))
        self.combo_aspect = QComboBox()
        self.combo_aspect.addItems(["Custom", "16:9", "4:3", "3:2", "1:1", "2:1"])
        ar_row.addWidget(self.combo_aspect)
        fl.addLayout(ar_row)

        # Scale bar controls
        sbar_row = QHBoxLayout()
        sbar_row.addWidget(QLabel("Scale bar:"))
        self.combo_sbar_orient = QComboBox()
        self.combo_sbar_orient.addItems(["Vertical", "Horizontal"])
        sbar_row.addWidget(self.combo_sbar_orient)
        fl.addLayout(sbar_row)
        font_row = QHBoxLayout()
        font_row.addWidget(QLabel("Label size:"))
        self.spin_sbar_fontsize = QSpinBox()
        self.spin_sbar_fontsize.setRange(6, 36)
        self.spin_sbar_fontsize.setValue(12)
        font_row.addWidget(self.spin_sbar_fontsize)
        fl.addLayout(font_row)

        # Scale bar title override
        title_row = QHBoxLayout()
        title_row.addWidget(QLabel("Title:"))
        self.txt_sbar_title = QLineEdit()
        self.txt_sbar_title.setPlaceholderText("(auto)")
        title_row.addWidget(self.txt_sbar_title)
        fl.addLayout(title_row)

        # Scale range clamping
        fl.addWidget(QLabel("Scale range:"))
        min_row = QHBoxLayout()
        self.chk_clamp_min = QCheckBox("Min:")
        self.chk_clamp_min.setChecked(False)
        min_row.addWidget(self.chk_clamp_min)
        self.spin_clamp_min = QDoubleSpinBox()
        self.spin_clamp_min.setRange(-1e20, 1e20)
        self.spin_clamp_min.setDecimals(4)
        self.spin_clamp_min.setValue(0)
        self.spin_clamp_min.setEnabled(False)
        min_row.addWidget(self.spin_clamp_min)
        fl.addLayout(min_row)
        max_row = QHBoxLayout()
        self.chk_clamp_max = QCheckBox("Max:")
        self.chk_clamp_max.setChecked(False)
        max_row.addWidget(self.chk_clamp_max)
        self.spin_clamp_max = QDoubleSpinBox()
        self.spin_clamp_max.setRange(-1e20, 1e20)
        self.spin_clamp_max.setDecimals(4)
        self.spin_clamp_max.setValue(1)
        self.spin_clamp_max.setEnabled(False)
        max_row.addWidget(self.spin_clamp_max)
        fl.addLayout(max_row)

        # GP sphere zero-size mode
        gp_zero_row = QHBoxLayout()
        gp_zero_row.addWidget(QLabel("GP size zero at:"))
        self.combo_gp_zero = QComboBox()
        self.combo_gp_zero.addItems(["Zero value", "Min value"])
        self.combo_gp_zero.setToolTip(
            "Zero value: size=0 at value 0 (for stress, +/- values)\n"
            "Min value: size=0 at min (for damage, 0-to-max values)"
        )
        gp_zero_row.addWidget(self.combo_gp_zero)
        fl.addLayout(gp_zero_row)

        # Scale bar dimensions (fraction of viewport, 0-1)
        bw_row = QHBoxLayout()
        bw_row.addWidget(QLabel("Bar W:"))
        self.spin_sbar_width = QDoubleSpinBox()
        self.spin_sbar_width.setRange(0.01, 1.0)
        self.spin_sbar_width.setSingleStep(0.01)
        self.spin_sbar_width.setDecimals(2)
        self.spin_sbar_width.setValue(0.08)
        bw_row.addWidget(self.spin_sbar_width)
        bw_row.addWidget(QLabel("H:"))
        self.spin_sbar_height = QDoubleSpinBox()
        self.spin_sbar_height.setRange(0.01, 1.0)
        self.spin_sbar_height.setSingleStep(0.01)
        self.spin_sbar_height.setDecimals(2)
        self.spin_sbar_height.setValue(0.45)
        bw_row.addWidget(self.spin_sbar_height)
        fl.addLayout(bw_row)

        # Scale bar position
        pos_row = QHBoxLayout()
        pos_row.addWidget(QLabel("Position:"))
        self.combo_sbar_pos = QComboBox()
        self.combo_sbar_pos.addItems([
            "Top-Left", "Top-Right", "Bottom-Left", "Bottom-Right",
        ])
        pos_row.addWidget(self.combo_sbar_pos)
        fl.addLayout(pos_row)

        figure_group.setLayout(fl)
        vl_top.addWidget(figure_group)

        # --- Export ---
        export_group = QGroupBox("Export")
        el = QVBoxLayout()
        self.btn_screenshot = QPushButton("Save Figure")
        el.addWidget(self.btn_screenshot)
        self.btn_record = QPushButton("Record Animation")
        el.addWidget(self.btn_record)
        self.lbl_record_status = QLabel("")
        el.addWidget(self.lbl_record_status)
        export_group.setLayout(el)
        vl_top.addWidget(export_group)

        vl_top.addStretch()
        tabs.addTab(view_tab, "View")

        # =====================================================================
        dock.setWidget(tabs)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

        # --- Connect signals ---
        self.combo_category.currentIndexChanged.connect(self._on_category_changed)
        self.combo_result.currentIndexChanged.connect(self._on_result_changed)
        self.combo_component.currentIndexChanged.connect(self._on_display_param_changed)
        self.combo_elem_mode.currentIndexChanged.connect(self._on_display_param_changed)
        self.combo_gp.currentIndexChanged.connect(self._on_display_param_changed)
        self.combo_fiber.currentIndexChanged.connect(self._on_display_param_changed)
        self.combo_cmap.currentIndexChanged.connect(self._on_display_param_changed)
        self.spin_sphere.valueChanged.connect(self._on_display_param_changed)
        self.combo_gp_zero.currentIndexChanged.connect(self._on_display_param_changed)
        self.spin_grid.valueChanged.connect(self._on_grid_spacing_changed)

        self.slider_step.valueChanged.connect(self._on_step_slider_changed)
        self.slider_scale.valueChanged.connect(self._on_scale_slider_changed)
        self.spin_scale.valueChanged.connect(self._on_scale_spin_changed)

        self.chk_shells.stateChanged.connect(self._on_visibility_changed)
        self.chk_beams.stateChanged.connect(self._on_visibility_changed)
        self.chk_edges.stateChanged.connect(self._on_visibility_changed)
        self.chk_extrusion.stateChanged.connect(self._on_visibility_changed)
        self.chk_fiber_lines.stateChanged.connect(self._on_visibility_changed)
        self.chk_grid.stateChanged.connect(self._on_visibility_changed)

        self.btn_first.clicked.connect(lambda: self._goto_step(0))
        self.btn_back.clicked.connect(lambda: self._goto_step(self._current_step_idx - 1))
        self.btn_next.clicked.connect(lambda: self._goto_step(self._current_step_idx + 1))
        self.btn_last.clicked.connect(lambda: self._goto_step(len(self._steps) - 1))
        self.btn_play.clicked.connect(self._on_play_pause)
        self.btn_screenshot.clicked.connect(self._on_screenshot)
        self.btn_record.clicked.connect(self._on_record_toggle)
        self.combo_aspect.currentIndexChanged.connect(self._on_aspect_changed)

        # Figure controls — all trigger immediate rebuild
        self.spin_fig_w.valueChanged.connect(self._on_figure_changed)
        self.spin_fig_h.valueChanged.connect(self._on_figure_changed)
        self.combo_sbar_orient.currentIndexChanged.connect(self._on_figure_changed)
        self.spin_sbar_fontsize.valueChanged.connect(self._on_figure_changed)
        self.txt_sbar_title.textChanged.connect(self._on_figure_changed)
        self.spin_sbar_width.valueChanged.connect(self._on_figure_changed)
        self.spin_sbar_height.valueChanged.connect(self._on_figure_changed)
        self.combo_sbar_pos.currentIndexChanged.connect(self._on_figure_changed)
        self.combo_cmap.currentIndexChanged.connect(self._on_figure_changed)
        self.chk_clamp_min.toggled.connect(self._on_figure_changed)
        self.chk_clamp_max.toggled.connect(self._on_figure_changed)
        self.spin_clamp_min.valueChanged.connect(self._on_figure_changed)
        self.spin_clamp_max.valueChanged.connect(self._on_figure_changed)

        self.chk_clamp_min.toggled.connect(self.spin_clamp_min.setEnabled)
        self.chk_clamp_max.toggled.connect(self.spin_clamp_max.setEnabled)

        self.btn_front.clicked.connect(lambda: self._set_view("front"))
        self.btn_back_view.clicked.connect(lambda: self._set_view("back"))
        self.btn_left.clicked.connect(lambda: self._set_view("left"))
        self.btn_right.clicked.connect(lambda: self._set_view("right"))
        self.btn_top.clicked.connect(lambda: self._set_view("top"))
        self.btn_bottom.clicked.connect(lambda: self._set_view("bottom"))
        self.btn_iso.clicked.connect(lambda: self._set_view("iso"))
        self.btn_persp.clicked.connect(self._toggle_perspective)
        self.combo_sbar_orient.currentIndexChanged.connect(self._on_sbar_orient_changed)

        # Recording state
        self._recording = False
        self._record_frames = []

    def _populate_dropdowns(self):
        self._available = self.processor.get_available_results()
        self.combo_category.blockSignals(True)
        self.combo_category.clear()
        categories = []
        if self._available.get("Nodal"):
            categories.append("Nodal")
        if self._available.get("Element"):
            categories.append("Element")
        self.combo_category.addItems(categories)
        self.combo_category.blockSignals(False)
        if categories:
            self._on_category_changed()

    def _on_category_changed(self):
        cat = self.combo_category.currentText()
        results = self._available.get(cat, {})

        self.combo_result.blockSignals(True)
        self.combo_result.clear()
        self.combo_result.addItems(list(results.keys()))
        self.combo_result.blockSignals(False)

        self.elem_mode_group.setVisible(cat == "Element")
        if results:
            self._on_result_changed()

    def _on_result_changed(self):
        cat = self.combo_category.currentText()
        result_name = self.combo_result.currentText()
        info = self._available.get(cat, {}).get(result_name, {})

        self.combo_component.blockSignals(True)
        self.combo_component.clear()
        self.combo_component.addItems(info.get("components", []))
        self.combo_component.blockSignals(False)

        if cat == "Element":
            is_fiber = self._is_fiber_result()

            # For fiber results: hide GP/Fiber selectors and mode combo
            # (all GPs × fibers are shown as scaled spheres)
            self.combo_elem_mode.setVisible(not is_fiber)
            self.lbl_gp.setVisible(not is_fiber)
            self.combo_gp.setVisible(not is_fiber)
            self.lbl_fiber.setVisible(not is_fiber)
            self.combo_fiber.setVisible(not is_fiber)

            if not is_fiber:
                self.combo_gp.blockSignals(True)
                self.combo_gp.clear()
                self.combo_gp.addItems(info.get("gauss_points", []))
                self.combo_gp.blockSignals(False)

                max_fibers = info.get("max_fibers", 1)
                self.combo_fiber.blockSignals(True)
                self.combo_fiber.clear()
                self.combo_fiber.addItems([f"Fiber {i}" for i in range(max_fibers)])
                self.combo_fiber.blockSignals(False)

        self._on_display_param_changed()

    def _on_display_param_changed(self):
        self._sphere_radius = self.spin_sphere.value()
        self._needs_full_rebuild = True
        self._update_display()

    def _on_step_slider_changed(self, value):
        self._current_step_idx = value
        self._update_step_label()
        self._update_display()

    def _on_scale_slider_changed(self, value):
        self.spin_scale.blockSignals(True)
        self.spin_scale.setValue(value)
        self.spin_scale.blockSignals(False)
        self._scale_factor = float(value)
        self._update_display()

    def _on_scale_spin_changed(self, value):
        self.slider_scale.blockSignals(True)
        self.slider_scale.setValue(min(value, self.slider_scale.maximum()))
        self.slider_scale.blockSignals(False)
        self._scale_factor = float(value)
        self._update_display()

    def _on_grid_spacing_changed(self, value):
        self._grid_spacing = value
        self._rebuild_floor_grid()
        self._needs_full_rebuild = True
        self._update_display()

    def _on_aspect_changed(self):
        """Update height from width when an aspect ratio preset is selected."""
        preset = self.combo_aspect.currentText()
        if preset == "Custom":
            return
        ratios = {"16:9": 16/9, "4:3": 4/3, "3:2": 3/2, "1:1": 1, "2:1": 2}
        ratio = ratios.get(preset, 1)
        w = self.spin_fig_w.value()
        self.spin_fig_h.blockSignals(True)
        self.spin_fig_h.setValue(int(w / ratio))
        self.spin_fig_h.blockSignals(False)

    def _on_sbar_orient_changed(self):
        """Swap W and H spinner values when orientation changes."""
        w = self.spin_sbar_width.value()
        h = self.spin_sbar_height.value()
        self.spin_sbar_width.blockSignals(True)
        self.spin_sbar_height.blockSignals(True)
        self.spin_sbar_width.setValue(h)
        self.spin_sbar_height.setValue(w)
        self.spin_sbar_width.blockSignals(False)
        self.spin_sbar_height.blockSignals(False)

    def _on_figure_changed(self):
        """Apply figure settings immediately."""
        w = self.spin_fig_w.value()
        h = self.spin_fig_h.value()
        self.plotter_widget.setFixedSize(w, h)
        self.adjustSize()
        self._needs_full_rebuild = True
        self._update_display()

    def _setup_initial_camera(self):
        """Position camera so the model is centered and fills the view.

        Camera looks from a 3/4 view (front-right, slightly above) with
        the model's geometric center as focal point.
        """
        coords = self.model.nodes.coords
        center = self._model_center
        diag = np.linalg.norm(coords.max(axis=0) - coords.min(axis=0))

        # Camera distance: far enough to see the whole model
        dist = diag * 1.5

        # 3/4 view direction: slightly from the right and above
        up = self.builder._up_vector
        cam_pos = center + np.array([dist * 0.6, -dist * 0.7, dist * 0.4])

        self.plotter_widget.camera.focal_point = center
        self.plotter_widget.camera.position = cam_pos
        self.plotter_widget.camera.up = up
        self.plotter_widget.renderer.ResetCameraClippingRange()

    def _set_view(self, direction):
        """Set camera to a standard orthographic view direction.

        All views look at the model center from a distance that fills the frame.
        """
        center = self._model_center
        coords = self.model.nodes.coords
        diag = np.linalg.norm(coords.max(axis=0) - coords.min(axis=0))
        dist = diag * 1.5
        up = self.builder._up_vector

        # Direction vectors: camera position = center + offset
        views = {
            "front":  np.array([0, -dist, 0]),
            "back":   np.array([0, +dist, 0]),
            "left":   np.array([-dist, 0, 0]),
            "right":  np.array([+dist, 0, 0]),
            "top":    np.array([0, 0, +dist]),
            "bottom": np.array([0, 0, -dist]),
            "iso":    np.array([dist * 0.6, -dist * 0.7, dist * 0.4]),
        }
        # Up vector for top/bottom views needs to change
        up_overrides = {
            "top":    np.array([0, -1, 0]),
            "bottom": np.array([0, 1, 0]),
        }

        offset = views.get(direction, views["front"])
        view_up = up_overrides.get(direction, up)

        self.plotter_widget.camera.focal_point = center
        self.plotter_widget.camera.position = center + offset
        self.plotter_widget.camera.up = view_up
        self.plotter_widget.renderer.ResetCameraClippingRange()
        self.plotter_widget.render()

    def _toggle_perspective(self):
        """Switch between perspective and orthographic projection."""
        cam = self.plotter_widget.renderer.GetActiveCamera()
        if self.btn_persp.isChecked():
            cam.ParallelProjectionOn()
            # Set parallel scale to frame the model
            coords = self.model.nodes.coords
            diag = np.linalg.norm(coords.max(axis=0) - coords.min(axis=0))
            cam.SetParallelScale(diag * 0.6)
            self.btn_persp.setText("Perspective")
        else:
            cam.ParallelProjectionOff()
            self.btn_persp.setText("Orthographic")
        self.plotter_widget.renderer.ResetCameraClippingRange()
        self.plotter_widget.render()

    def _reset_figure_size(self):
        """Remove fixed size constraint (restore resizable viewport)."""
        self.plotter_widget.setMinimumSize(0, 0)
        self.plotter_widget.setMaximumSize(16777215, 16777215)

    @staticmethod
    def _auto_grid_spacing(max_span):
        """Pick grid spacing based on model scale.

        Finds the order of magnitude of the model, then uses that as spacing.
        Models < 10 -> spacing 0.1, < 100 -> spacing 1, < 1000 -> spacing 10, etc.
        """
        if max_span <= 0:
            return 0.1
        order = 10 ** np.floor(np.log10(max_span))
        return order / 10.0

    def _rebuild_floor_grid(self):
        """Rebuild floor grid perpendicular to detected up-axis.

        Always square, sized to cover the model with padding.
        Grid extends in multiples of (spacing * 10) blocks beyond the model.
        """
        coords = self.model.nodes.coords
        up = self.builder._up_axis
        direction = np.zeros(3)
        direction[up] = 1.0

        # Grid center at the floor level
        center = np.zeros(3)
        center[up] = coords[:, up].min()

        # Square grid: use the largest model dimension
        max_span = max(coords[:, i].max() - coords[:, i].min() for i in range(3))
        s = self._grid_spacing
        block = s * 10  # one 10x10 block

        # Grid size: round up to the next full block, with at least 1 block padding
        half = max_span / 2.0 + block
        half = np.ceil(half / block) * block
        size = half * 2
        n = max(1, int(round(size / s)))

        self._floor_grid = pv.Plane(
            center=tuple(center), direction=tuple(direction),
            i_size=size, j_size=size,
            i_resolution=n, j_resolution=n,
        )

    def _on_visibility_changed(self):
        self._show_shells = self.chk_shells.isChecked()
        self._show_beams = self.chk_beams.isChecked()
        self._show_edges = self.chk_edges.isChecked()
        self._show_extrusion = self.chk_extrusion.isChecked()
        self._show_fibers = self.chk_fiber_lines.isChecked()
        self._show_grid = self.chk_grid.isChecked()
        self._needs_full_rebuild = True
        self._update_display()

    def _goto_step(self, idx):
        idx = max(0, min(idx, len(self._steps) - 1))
        self.slider_step.setValue(idx)

    def _on_play_pause(self):
        if self._playing:
            self._playing = False
            self._timer.stop()
            self.btn_play.setText("Play")
        else:
            self._playing = True
            self.btn_play.setText("Pause")
            interval = max(30, int(1000 / self.spin_fps.value()))
            self._timer.start(interval)

    def _on_timer_tick(self):
        next_idx = self._current_step_idx + 1
        if next_idx >= len(self._steps):
            if self._recording:
                # Recording finished — stop play and save
                self._stop_recording()
                return
            next_idx = 0
        self._goto_step(next_idx)

    def _update_step_label(self):
        total = len(self._steps)
        self.lbl_step_info.setText(f"Frame {self._current_step_idx + 1} / {total}")
        if self._steps:
            stage_idx, step_key = self._steps[self._current_step_idx]
            step_num = int(step_key.split("_")[1])
            self.lbl_stage_info.setText(f"Stage {stage_idx}, {step_key} (t={step_num})")

    def _update_display(self):
        """Main render update.

        On full rebuild (visibility/result/colormap change): clears scene and
        re-adds all actors, saving references for fast updates.

        On fast update (step/scale change during animation): only updates
        point positions and scalar arrays in-place — no scene rebuild.
        """
        if not self._steps:
            return

        stage_idx, step_key = self._steps[self._current_step_idx]
        self._update_step_label()

        h5 = self._h5file
        disp = self.processor.get_displacement_data(self.mpco_path, stage_idx, step_key, h5file=h5)
        cat = self.combo_category.currentText()
        cmap = self.combo_cmap.currentText() or "jet"

        if self._needs_full_rebuild:
            self._full_rebuild(stage_idx, step_key, disp, cat, cmap, h5)
            self._needs_full_rebuild = False
        else:
            self._fast_update(stage_idx, step_key, disp, cat, h5)

    def _full_rebuild(self, stage_idx, step_key, disp, cat, cmap, h5):
        """Clear scene and rebuild all actors from scratch."""
        cam_pos = self.plotter_widget.camera_position
        has_camera = cam_pos is not None and self.plotter_widget.renderer.GetActors().GetNumberOfItems() > 0

        self.plotter_widget.renderer.SetDraw(False)
        self.plotter_widget.clear()
        self._actors.clear()
        self._meshes.clear()

        # --- Floor grid ---
        if self._show_grid:
            self.plotter_widget.add_mesh(
                self._floor_grid, color="lightgray", style="wireframe",
                line_width=0.5, opacity=0.3, reset_camera=False,
            )

        # --- Shell extrusion ---
        if self._show_shells and self._show_extrusion and self.shell_extrusion_base.n_cells > 0:
            ext_mesh = self.shell_extrusion_base.copy(deep=True)
            if disp is not None:
                self._apply_disp_inplace(ext_mesh, disp, self.builder.shell_node_to_hex_points)
            self._meshes["extrusion"] = ext_mesh
            self._actors["extrusion"] = self.plotter_widget.add_mesh(
                ext_mesh, color="lightgray", opacity=0.15,
                show_edges=True, edge_color="gray", line_width=0.5,
                style="surface", reset_camera=False,
            )

        # --- Fiber layer lines (hidden when GP spheres are active) ---
        show_gp = cat == "Element" and (self._is_fiber_result() or not self._is_contour_mode())
        if self._show_shells and self._show_fibers and not show_gp and self.fiber_layer_edges_base.n_points > 0:
            fiber_mesh = self.fiber_layer_edges_base.copy(deep=True)
            if disp is not None:
                self._apply_disp_inplace(fiber_mesh, disp, self.builder.shell_node_to_layer_points)
            self._meshes["fibers"] = fiber_mesh
            self._actors["fibers"] = self.plotter_widget.add_mesh(
                fiber_mesh, color="darkgray", line_width=0.5, opacity=0.3,
                reset_camera=False,
            )

        # --- Shell mid-surface ---
        if self._show_shells and self.shell_surface_base.n_cells > 0:
            surf_mesh = self.shell_surface_base.copy(deep=True)
            if disp is not None:
                self._apply_disp_inplace(surf_mesh, disp, self.builder.shell_node_to_layer_points)

            scalar_name = self._get_scalar_label()
            use_contour = (cat == "Nodal") or (cat == "Element" and not self._is_fiber_result())
            surface_scalars = self._extract_surface_scalars(stage_idx, step_key, h5) if use_contour else None

            if surface_scalars is not None and len(surface_scalars) == surf_mesh.n_points:
                surf_mesh.point_data[scalar_name] = surface_scalars
                self._meshes["surface"] = surf_mesh
                clim = self._get_clim()
                mesh_kwargs = dict(
                    scalars=scalar_name, cmap=cmap,
                    show_edges=self._show_edges, nan_color="lightgray",
                    scalar_bar_args=self._scalar_bar_args(scalar_name),
                    opacity=0.9, reset_camera=False,
                )
                if clim is not None:
                    mesh_kwargs["clim"] = clim
                self._actors["surface"] = self.plotter_widget.add_mesh(
                    surf_mesh, **mesh_kwargs,
                )
                self._update_min_max(surface_scalars)
            else:
                self._meshes["surface"] = surf_mesh
                self._actors["surface"] = self.plotter_widget.add_mesh(
                    surf_mesh, color="lightblue",
                    show_edges=self._show_edges, opacity=0.7,
                    reset_camera=False,
                )

        # --- GP Spheres ---
        if self._show_shells and cat == "Element" and (self._is_fiber_result() or not self._is_contour_mode()):
            self._add_gp_spheres(stage_idx, step_key, disp, cmap, h5)

        # --- Beams ---
        if self._show_beams and self.beam_mesh_base.n_cells > 0:
            beam_mesh = self.beam_mesh_base.copy(deep=True)
            if disp is not None:
                self._apply_disp_inplace(beam_mesh, disp, self.builder.beam_node_to_hex_points)
            self._meshes["beams"] = beam_mesh
            self._actors["beams"] = self.plotter_widget.add_mesh(
                beam_mesh, color="sienna",
                show_edges=self._show_edges, opacity=0.7,
                reset_camera=False,
            )

        # Restore camera
        self.plotter_widget.renderer.SetDraw(True)
        if has_camera:
            self.plotter_widget.camera_position = cam_pos
        else:
            self._setup_initial_camera()
        self.plotter_widget.renderer.ResetCameraClippingRange()
        self.plotter_widget.render()
        self._capture_frame()

    def _fast_update(self, stage_idx, step_key, disp, cat, h5):
        """Update point positions and scalars in-place without rebuilding actors."""
        self.plotter_widget.renderer.SetDraw(False)

        # Update mesh point positions
        if "extrusion" in self._meshes:
            pts = self.shell_extrusion_base.points.copy()
            if disp is not None:
                self._apply_disp_to_array(pts, disp, self.builder.shell_node_to_hex_points)
            self._meshes["extrusion"].points[:] = pts

        if "fibers" in self._meshes:
            pts = self.fiber_layer_edges_base.points.copy()
            if disp is not None:
                self._apply_disp_to_array(pts, disp, self.builder.shell_node_to_layer_points)
            self._meshes["fibers"].points[:] = pts

        if "surface" in self._meshes:
            pts = self.shell_surface_base.points.copy()
            if disp is not None:
                self._apply_disp_to_array(pts, disp, self.builder.shell_node_to_layer_points)
            surf = self._meshes["surface"]
            surf.points[:] = pts

            # Update scalars
            use_contour = (cat == "Nodal") or (cat == "Element" and not self._is_fiber_result())
            if use_contour:
                scalars = self._extract_surface_scalars(stage_idx, step_key, h5)
                if scalars is not None and len(scalars) == surf.n_points:
                    scalar_name = self._get_scalar_label()
                    surf.point_data[scalar_name] = scalars
                    self._update_min_max(scalars)

        if "beams" in self._meshes:
            pts = self.beam_mesh_base.points.copy()
            if disp is not None:
                self._apply_disp_to_array(pts, disp, self.builder.beam_node_to_hex_points)
            self._meshes["beams"].points[:] = pts

        # GP spheres need full rebuild each step (different scalars + glyph geometry)
        if "gp_spheres" in self._actors:
            self.plotter_widget.remove_actor(self._actors["gp_spheres"])
            del self._actors["gp_spheres"]
            self._meshes.pop("gp_spheres", None)
        if self._show_shells and cat == "Element" and (self._is_fiber_result() or not self._is_contour_mode()):
            cmap = self.combo_cmap.currentText() or "jet"
            self._add_gp_spheres(stage_idx, step_key, disp, cmap, h5)

        self.plotter_widget.renderer.SetDraw(True)
        self.plotter_widget.renderer.ResetCameraClippingRange()
        self.plotter_widget.render()
        self._capture_frame()

    def _apply_disp_inplace(self, mesh, disp, node_to_points):
        """Apply displacement to a mesh's points in-place."""
        nodes = self.model.nodes
        for nid, pt_indices in node_to_points.items():
            node_idx = nodes.id_to_index.get(nid)
            if node_idx is None:
                continue
            d = disp[node_idx] * self._scale_factor
            for pi in pt_indices:
                if pi < mesh.n_points:
                    mesh.points[pi] += d

    def _apply_disp_to_array(self, pts, disp, node_to_points):
        """Apply displacement to a points array (copy of base)."""
        nodes = self.model.nodes
        for nid, pt_indices in node_to_points.items():
            node_idx = nodes.id_to_index.get(nid)
            if node_idx is None:
                continue
            d = disp[node_idx] * self._scale_factor
            for pi in pt_indices:
                if pi < len(pts):
                    pts[pi] += d

    def _scalar_bar_args(self, title):
        """Build scalar bar args from current figure settings."""
        # Title override
        custom_title = self.txt_sbar_title.text().strip()
        if custom_title:
            title = custom_title

        vertical = self.combo_sbar_orient.currentText() == "Vertical"
        fontsize = self.spin_sbar_fontsize.value()
        bar_w = self.spin_sbar_width.value()
        bar_h = self.spin_sbar_height.value()
        pos_name = self.combo_sbar_pos.currentText()

        n_labels = 8 if vertical else 5
        margin = 0.05

        # Compute position from preset
        if vertical:
            if "Right" in pos_name:
                pos_x = 1.0 - bar_w - margin
            else:
                pos_x = margin
            if "Bottom" in pos_name:
                pos_y = margin
            else:
                pos_y = 1.0 - bar_h - margin
        else:
            if "Right" in pos_name:
                pos_x = 1.0 - bar_w - margin
            else:
                pos_x = margin
            if "Top" in pos_name:
                pos_y = 1.0 - bar_h - margin
            else:
                pos_y = margin

        return {
            "title": title,
            "n_labels": n_labels,
            "vertical": vertical,
            "title_font_size": fontsize + 2,
            "label_font_size": fontsize,
            "width": bar_w,
            "height": bar_h,
            "position_x": pos_x,
            "position_y": pos_y,
            "fmt": "%.3e",
        }

    def _get_clim(self):
        """Get user-defined color limits, or None for auto."""
        has_min = self.chk_clamp_min.isChecked()
        has_max = self.chk_clamp_max.isChecked()
        if not has_min and not has_max:
            return None
        lo = self.spin_clamp_min.value() if has_min else None
        hi = self.spin_clamp_max.value() if has_max else None
        return [lo, hi]

    def _is_contour_mode(self):
        return self.combo_elem_mode.currentIndex() == 0

    def _is_fiber_result(self):
        """Check if current element result is fiber-level (e.g. section.fiber.stress)."""
        result_name = self.combo_result.currentText()
        return "fiber" in result_name.lower()

    def _extract_surface_scalars(self, stage_idx, step_key, h5=None):
        """Extract scalars for the mid-surface mesh."""
        cat = self.combo_category.currentText()
        result_name = self.combo_result.currentText()
        comp_text = self.combo_component.currentText()
        if not result_name or not comp_text:
            return None

        if cat == "Nodal":
            comps = self._available["Nodal"][result_name]["components"]
            comp_idx = comps.index(comp_text) if comp_text in comps else 0
            return self.processor.extract_nodal_result(
                self.mpco_path, stage_idx, step_key, result_name, comp_idx,
                h5file=h5
            )
        else:
            comps = self._available["Element"][result_name]["components"]
            comp_idx = comps.index(comp_text) if comp_text in comps else 0

            gp_text = self.combo_gp.currentText()
            gp = -1 if gp_text == "Average" else int(gp_text.replace("GP ", ""))

            fiber_text = self.combo_fiber.currentText()
            fiber_idx = int(fiber_text.replace("Fiber ", "")) if fiber_text else 0

            return self.processor.extract_element_result_contour(
                self.mpco_path, stage_idx, step_key, result_name,
                gp, fiber_idx, comp_idx, h5file=h5
            )

    def _add_gp_spheres(self, stage_idx, step_key, disp, cmap, h5=None):
        """Add colored spheres at Gauss point positions.

        For fiber results: all GPs × all fibers, sphere radius ~ |value|.
        For non-fiber results: filtered by selected GP/fiber, uniform size.
        Uses vectorized displacement interpolation.
        """
        result_name = self.combo_result.currentText()
        comp_text = self.combo_component.currentText()
        if not result_name or not comp_text:
            return

        comps = self._available["Element"][result_name]["components"]
        comp_idx = comps.index(comp_text) if comp_text in comps else 0

        is_fiber = self._is_fiber_result()

        if is_fiber:
            mask = None  # all points
            pts = self._gp_points_base.copy()
            eids = self._gp_elem_ids
            gp_idx = self._gp_indices
            fib_idx = self._gp_fiber_indices
        else:
            fiber_text = self.combo_fiber.currentText()
            fiber_idx = int(fiber_text.replace("Fiber ", "")) if fiber_text else 0
            mask = self._gp_fiber_indices == fiber_idx

            gp_text = self.combo_gp.currentText()
            if gp_text != "Average":
                gp_sel = int(gp_text.replace("GP ", ""))
                mask = mask & (self._gp_indices == gp_sel)

            if not np.any(mask):
                return

            pts = self._gp_points_base[mask].copy()
            eids = self._gp_elem_ids[mask]
            gp_idx = self._gp_indices[mask]
            fib_idx = self._gp_fiber_indices[mask]

        # Vectorized displacement interpolation
        if disp is not None:
            pts = self.builder.displace_gp_points(pts, disp, self._scale_factor, mask)

        # Extract result values
        scalars = self.processor.extract_element_result_for_gp_cloud(
            self.mpco_path, stage_idx, step_key, result_name, comp_idx,
            eids, gp_idx, fib_idx, h5file=h5
        )

        if scalars is None:
            return

        scalar_name = self._get_scalar_label()

        clim = self._get_clim()

        if is_fiber:
            # Compute sphere radii based on zero-mode setting
            vals = np.nan_to_num(scalars, nan=0.0)
            zero_at_zero = self.combo_gp_zero.currentText() == "Zero value"

            bbox_diag = np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)) if len(pts) > 1 else 1.0
            max_radius = bbox_diag * self._sphere_radius / 1000.0

            if zero_at_zero:
                # Size proportional to |value| — zero at value=0
                abs_vals = np.abs(vals)
                max_val = abs_vals.max()
                if max_val > 0:
                    radii = max_radius * abs_vals / max_val
                else:
                    radii = np.full_like(abs_vals, max_radius * 0.1)
            else:
                # Size proportional to (value - min) — zero at min value
                vmin = vals.min()
                vmax = vals.max()
                span = vmax - vmin
                if span > 0:
                    radii = max_radius * (vals - vmin) / span
                else:
                    radii = np.full_like(vals, max_radius * 0.1)

            radii = np.maximum(radii, max_radius * 0.05)

            cloud = pv.PolyData(pts)
            cloud.point_data[scalar_name] = scalars
            cloud.point_data["_radius"] = radii

            sphere_src = pv.Sphere(radius=1.0, theta_resolution=6, phi_resolution=6)
            glyphs = cloud.glyph(geom=sphere_src, scale="_radius", orient=False)

            mesh_kwargs = dict(
                scalars=scalar_name, cmap=cmap,
                nan_color="lightgray",
                scalar_bar_args=self._scalar_bar_args(scalar_name),
                reset_camera=False,
            )
            if clim is not None:
                mesh_kwargs["clim"] = clim
            self._actors["gp_spheres"] = self.plotter_widget.add_mesh(
                glyphs, **mesh_kwargs,
            )
        else:
            cloud = pv.PolyData(pts)
            cloud.point_data[scalar_name] = scalars

            mesh_kwargs = dict(
                scalars=scalar_name, cmap=cmap,
                render_points_as_spheres=True,
                point_size=int(self._sphere_radius),
                nan_color="lightgray",
                scalar_bar_args=self._scalar_bar_args(scalar_name),
                reset_camera=False,
            )
            if clim is not None:
                mesh_kwargs["clim"] = clim
            self._actors["gp_spheres"] = self.plotter_widget.add_mesh(
                cloud, **mesh_kwargs,
            )

        self._update_min_max(scalars)

    def _get_scalar_label(self):
        cat = self.combo_category.currentText()
        result = self.combo_result.currentText()
        comp = self.combo_component.currentText()
        if cat == "Element":
            if self._is_fiber_result():
                return f"{result} | {comp} [All GPs \u00d7 Fibers]"
            mode = "Contour" if self._is_contour_mode() else "GP Spheres"
            gp = self.combo_gp.currentText()
            fiber = self.combo_fiber.currentText()
            return f"{result} | {comp} | {gp} | {fiber} [{mode}]"
        return f"{result} | {comp}"

    def _update_min_max(self, scalars):
        valid = scalars[~np.isnan(scalars)] if scalars is not None else np.array([])
        if len(valid) > 0:
            self.lbl_min_max.setText(f"Min: {valid.min():.4e}  Max: {valid.max():.4e}")
        else:
            self.lbl_min_max.setText("")

    # --- Export ---

    def _on_screenshot(self):
        """Save current view as raster (PNG at 1x-4x) or vector (SVG/PDF)."""
        from PyQt5.QtWidgets import QFileDialog, QInputDialog
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Figure",
            os.path.join(os.path.dirname(self.mpco_path), "figure.png"),
            "PNG Image (*.png);;SVG Vector (*.svg);;PDF Document (*.pdf)",
        )
        if not path:
            return

        ext = os.path.splitext(path)[1].lower()

        if ext in (".svg", ".pdf"):
            # Vector export via VTK GL2PS — sharp at any zoom, includes scalar bar
            self.plotter_widget.save_graphic(path)
        else:
            # High-res raster — ask for scale factor
            scale, ok = QInputDialog.getItem(
                self, "Resolution", "Scale factor:",
                ["1x (screen)", "2x", "3x", "4x"], 1, False,
            )
            if not ok:
                return
            factor = int(scale[0])

            if factor == 1:
                # Grab Qt widget directly (includes scalar bar, overlays)
                pixmap = self.plotter_widget.grab()
                pixmap.save(path)
            else:
                # Render at higher resolution via VTK
                w, h = self.plotter_widget.window_size
                self.plotter_widget.screenshot(
                    path, window_size=(w * factor, h * factor),
                )

        self.lbl_record_status.setText(f"Saved: {os.path.basename(path)}")

    def _on_record_toggle(self):
        """Start or cancel recording.

        Flow: ask file path first, then start recording + play.
        Recording stops automatically when play reaches the last step.
        """
        if self._recording:
            # Manual cancel
            self._stop_recording()
            return

        # Ask for output file FIRST
        filters = "GIF Animation (*.gif)"
        if self._has_ffmpeg():
            filters += ";;MP4 Video (*.mp4)"

        from PyQt5.QtWidgets import QFileDialog
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Animation",
            os.path.join(os.path.dirname(self.mpco_path), "animation.gif"),
            filters,
        )
        if not path:
            return

        # Start recording
        self._recording = True
        self._record_path = path
        self._record_frames = []
        self.btn_record.setText("Cancel Recording")
        self.lbl_record_status.setText("Recording...")

        # Capture the current frame, then start play
        self._capture_frame()
        if not self._playing:
            self._playing = True
            self.btn_play.setText("Pause")
            interval = max(30, int(1000 / self.spin_fps.value()))
            self._timer.start(interval)

    def _stop_recording(self):
        """Stop recording, stop play, save the file."""
        self._recording = False
        self.btn_record.setText("Record Animation")

        # Stop play
        if self._playing:
            self._playing = False
            self._timer.stop()
            self.btn_play.setText("Play")

        # Save
        self._save_recording()

    def _capture_frame(self):
        """Capture current frame to the recording buffer."""
        if self._recording:
            img = self.plotter_widget.screenshot(return_img=True, window_size=None)
            if img is not None:
                self._record_frames.append(img)
                self.lbl_record_status.setText(
                    f"Recording... {len(self._record_frames)} frames"
                )

    def _save_recording(self):
        """Save recorded frames to the pre-selected file."""
        if not self._record_frames:
            self.lbl_record_status.setText("No frames recorded")
            return

        path = self._record_path
        ext = os.path.splitext(path)[1].lower()
        fps = self.spin_fps.value()

        try:
            if ext == ".mp4":
                self._save_mp4(path, fps)
            else:
                self._save_gif(path, fps)
            self.lbl_record_status.setText(
                f"Saved: {os.path.basename(path)} ({len(self._record_frames)} frames)"
            )
        except Exception as e:
            self.lbl_record_status.setText(f"Error: {e}")
        finally:
            self._record_frames.clear()

    def _save_gif(self, path, fps):
        """Save frames as animated GIF using Pillow."""
        from PIL import Image
        images = [Image.fromarray(f) for f in self._record_frames]
        duration = int(1000 / fps)
        images[0].save(
            path, save_all=True, append_images=images[1:],
            duration=duration, loop=0,
        )

    def _save_mp4(self, path, fps):
        """Save frames as MP4 using imageio + ffmpeg."""
        import imageio
        writer = imageio.get_writer(path, fps=fps)
        for frame in self._record_frames:
            writer.append_data(frame)
        writer.close()

    @staticmethod
    def _has_ffmpeg():
        """Check if ffmpeg is available for MP4 export."""
        try:
            import imageio_ffmpeg
            imageio_ffmpeg.get_ffmpeg_exe()
            return True
        except Exception:
            return False
