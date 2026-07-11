import os
import logging
import json
import traceback
import pyvista as pv
import numpy as np
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QLabel,
    QSlider,
    QHBoxLayout,
    QPushButton,
    QDoubleSpinBox,
    QComboBox,
)
from PyQt6.QtGui import QColor
from PyQt6.QtCore import Qt

# Try to import VDW radii from constants, fallback if needed
try:
    from moleditpy.utils.constants import pt, CPK_COLORS_PV
except ImportError:
    try:
        from rdkit import Chem

        pt = Chem.GetPeriodicTable()
    except Exception:
        pt = None
    CPK_COLORS_PV = {}  # Last resort fallback


# Plugin Metadata
PLUGIN_NAME = "VDW Radii Overlay"
PLUGIN_VERSION = "2026.07.11"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Visualizes VDW radii as a translucent surface overlay using PyVista. Refactored for V3 API."
PLUGIN_CONTEXT = None

SETTINGS_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "vdw_radii_overlay.json"
)

# Base model style shown underneath the translucent VDW overlay.
# "default" -> ball_and_stick (the plugin's original behavior)
# "stick"   -> thin sticks, so the overlay reads more clearly.
BASE_STYLE_OPTIONS = ["Default", "Stick"]
_BASE_STYLE_TO_OVERRIDE = {
    "default": "ball_and_stick",
    "stick": "stick",
}

# Global Settings (State)
_vdw_settings = {
    "occupancy": 0.3,  # Opacity (0.0 - 1.0)
    "resolution": 0.125,  # Voxel spacing in Angstroms
    "base_style": "default",  # "default" or "stick"
}


def load_settings():
    try:
        if os.path.exists(SETTINGS_FILE):
            with open(SETTINGS_FILE, "r") as f:
                saved = json.load(f)
                if "occupancy" in saved:
                    _vdw_settings["occupancy"] = float(saved["occupancy"])
                if "resolution" in saved:
                    _vdw_settings["resolution"] = float(saved["resolution"])
                if "base_style" in saved and saved["base_style"] in _BASE_STYLE_TO_OVERRIDE:
                    _vdw_settings["base_style"] = saved["base_style"]
    except Exception as e:
        logging.warning("Error loading VDW settings: %s", e)


def save_settings():
    try:
        with open(SETTINGS_FILE, "w") as f:
            json.dump(_vdw_settings, f, indent=4)
    except Exception as e:
        logging.warning("Error saving VDW settings: %s", e)


class VDWConfigWindow(QDialog):
    """
    Namespaced configuration window for VDW Overlay.
    Refactored for MoleditPy V3.0 API.
    """

    def __init__(self, context):
        super().__init__(parent=context.get_main_window())
        self.context = context
        self.setWindowTitle("VDW Overlay Settings")
        self.setModal(False)
        self.resize(350, 150)
        self.init_ui()

        # Register window for V3 lifecycle management
        self.context.register_window("main_panel", self)

    def init_ui(self):
        layout = QVBoxLayout()

        # Occupancy Slider
        occ_layout = QHBoxLayout()
        occ_layout.addWidget(QLabel("Occupancy:"))
        self.slider_occ = QSlider(Qt.Orientation.Horizontal)
        self.slider_occ.setRange(0, 100)
        current_occ = _vdw_settings.get("occupancy", 0.3)
        self.slider_occ.setValue(int(current_occ * 100))
        self.slider_occ.valueChanged.connect(self.on_occupancy_slider_changed)
        occ_layout.addWidget(self.slider_occ)

        self.spin_occ = QDoubleSpinBox()
        self.spin_occ.setRange(0.0, 1.0)
        self.spin_occ.setSingleStep(0.05)
        self.spin_occ.setValue(current_occ)
        self.spin_occ.valueChanged.connect(self.on_occupancy_spin_changed)
        occ_layout.addWidget(self.spin_occ)

        layout.addLayout(occ_layout)

        # Resolution Slider
        res_layout = QHBoxLayout()
        res_layout.addWidget(QLabel("Resolution (Å):"))
        self.slider_res = QSlider(Qt.Orientation.Horizontal)
        self.slider_res.setRange(5, 50)  # 0.05 to 0.50
        current_res = _vdw_settings.get("resolution", 0.125)
        self.slider_res.setValue(int(current_res * 100))
        self.slider_res.valueChanged.connect(self.on_resolution_slider_changed)
        res_layout.addWidget(self.slider_res)

        self.spin_res = QDoubleSpinBox()
        self.spin_res.setRange(0.05, 0.50)
        self.spin_res.setSingleStep(0.005)
        self.spin_res.setDecimals(3)
        self.spin_res.setValue(current_res)
        self.spin_res.valueChanged.connect(self.on_resolution_spin_changed)
        res_layout.addWidget(self.spin_res)

        layout.addLayout(res_layout)

        # Base Model Style
        base_layout = QHBoxLayout()
        base_layout.addWidget(QLabel("Base model:"))
        self.combo_base_style = QComboBox()
        self.combo_base_style.addItems(BASE_STYLE_OPTIONS)
        current_base = _vdw_settings.get("base_style", "default")
        idx = 1 if current_base == "stick" else 0
        self.combo_base_style.setCurrentIndex(idx)
        self.combo_base_style.currentIndexChanged.connect(
            self.on_base_style_changed
        )
        base_layout.addWidget(self.combo_base_style)
        layout.addLayout(base_layout)

        # Reset Button
        btn_reset = QPushButton("Reset to Defaults")
        btn_reset.clicked.connect(self.reset_defaults)
        layout.addWidget(btn_reset)

        # Close Button
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        layout.addWidget(btn_close)

        self.setLayout(layout)

    def on_occupancy_slider_changed(self, value):
        val_float = value / 100.0
        self.spin_occ.blockSignals(True)
        self.spin_occ.setValue(val_float)
        self.spin_occ.blockSignals(False)
        self._update_occupancy(val_float)

    def on_occupancy_spin_changed(self, value):
        val_int = int(value * 100)
        self.slider_occ.blockSignals(True)
        self.slider_occ.setValue(val_int)
        self.slider_occ.blockSignals(False)
        self._update_occupancy(value)

    def _update_occupancy(self, value):
        _vdw_settings["occupancy"] = value
        save_settings()
        self.update_view()

    def on_resolution_slider_changed(self, value):
        val_float = value / 100.0
        self.spin_res.blockSignals(True)
        self.spin_res.setValue(val_float)
        self.spin_res.blockSignals(False)
        self._update_resolution(val_float)

    def on_resolution_spin_changed(self, value):
        val_int = int(value * 100)
        self.slider_res.blockSignals(True)
        self.slider_res.setValue(val_int)
        self.slider_res.blockSignals(False)
        self._update_resolution(value)

    def _update_resolution(self, value):
        _vdw_settings["resolution"] = value
        save_settings()
        self.update_view()

    def on_base_style_changed(self, index):
        base_style = "stick" if index == 1 else "default"
        _vdw_settings["base_style"] = base_style
        save_settings()
        self.update_view()

    def reset_defaults(self):
        # Default values
        def_occ = 0.3
        def_res = 0.125
        def_base = "default"

        # Block signals to prevent redundant updates/saves during setting
        self.slider_occ.blockSignals(True)
        self.spin_occ.blockSignals(True)
        self.slider_res.blockSignals(True)
        self.spin_res.blockSignals(True)
        self.combo_base_style.blockSignals(True)

        # Set values
        self.slider_occ.setValue(int(def_occ * 100))
        self.spin_occ.setValue(def_occ)
        self.slider_res.setValue(int(def_res * 100))
        self.spin_res.setValue(def_res)
        self.combo_base_style.setCurrentIndex(0)

        # Unblock
        self.slider_occ.blockSignals(False)
        self.spin_occ.blockSignals(False)
        self.slider_res.blockSignals(False)
        self.spin_res.blockSignals(False)
        self.combo_base_style.blockSignals(False)

        # Update settings and view once
        _vdw_settings["occupancy"] = def_occ
        _vdw_settings["resolution"] = def_res
        _vdw_settings["base_style"] = def_base
        save_settings()
        self.update_view()

    def refresh_ui_values(self):
        """Update UI elements from global settings."""
        occ = _vdw_settings.get("occupancy", 0.3)
        res = _vdw_settings.get("resolution", 0.125)
        base_style = _vdw_settings.get("base_style", "default")

        self.slider_occ.blockSignals(True)
        self.spin_occ.blockSignals(True)
        self.slider_res.blockSignals(True)
        self.spin_res.blockSignals(True)
        self.combo_base_style.blockSignals(True)

        self.slider_occ.setValue(int(occ * 100))
        self.spin_occ.setValue(occ)
        self.slider_res.setValue(int(res * 100))
        self.spin_res.setValue(res)
        self.combo_base_style.setCurrentIndex(1 if base_style == "stick" else 0)

        self.slider_occ.blockSignals(False)
        self.spin_occ.blockSignals(False)
        self.slider_res.blockSignals(False)
        self.spin_res.blockSignals(False)
        self.combo_base_style.blockSignals(False)

    def update_view(self):
        # Trigger redraw via context (modern V3 pattern)
        self.context.refresh_3d_view()


def draw_vdw_overlay(mw, mol):
    """
    3D Rendering Callback for VDW Overlay Style.
    Receives (mw, mol) from the engine.
    """
    # 1. Draw the base model first (Ball & Stick by default, or Stick when
    # the user picked "Stick" in the settings dialog -- see BASE_STYLE_OPTIONS).
    base_style = _vdw_settings.get("base_style", "default")
    base_override = _BASE_STYLE_TO_OVERRIDE.get(base_style, "ball_and_stick")
    # # [DIRECT ACCESS] to manager
    if hasattr(mw, "view_3d_manager"):
        mw.view_3d_manager.draw_standard_3d_style(mol, style_override=base_override)
    else:
        # Emergency fallback for legacy systems
        print("VDW Warning: view_3d_manager not found, using legacy draw path.")
        if hasattr(mw, "view3d"):
            mw.view3d.draw_standard_3d_style(mol, style_override=base_override)

    # 2. Draw VDW Surface Overlay
    if mol and mol.GetNumAtoms() > 0:
        try:
            positions = []
            radii = []
            atom_colors = []

            # Use color overrides from the new namespaced dictionary on the manager
            custom_map = {}
            # # [DIRECT ACCESS] to manager state
            if hasattr(mw, "view_3d_manager"):
                custom_map = getattr(mw.view_3d_manager, "_plugin_color_overrides", {})

            if not custom_map:
                # # [DIRECT ACCESS] legacy
                custom_map = getattr(mw, "custom_atom_colors", {})

            if mol.GetNumConformers() > 0:
                conf = mol.GetConformer()
                for i in range(mol.GetNumAtoms()):
                    atom = mol.GetAtomWithIdx(i)
                    pos = conf.GetAtomPosition(i)
                    positions.append([pos.x, pos.y, pos.z])

                    sym = atom.GetSymbol()
                    # GetRvdw returns standard VDW radius.
                    r = 1.5  # Default
                    if pt:
                        r = pt.GetRvdw(atom.GetAtomicNum())
                    radii.append(r)

                    # Color handling
                    if i in custom_map:
                        val = custom_map[i]
                        # Handling new API (Hex string) vs Legacy (List/Tuple)
                        if isinstance(val, str) and val.startswith("#"):
                            # Convert Hex to RGB [0-1]
                            qc = QColor(val)
                            c = [qc.redF(), qc.greenF(), qc.blueF()]
                        else:
                            # Assume legacy list/tuple
                            c = val
                            # Normalize 0-255 to 0-1 if needed
                            if any(x > 1.0 for x in c):
                                c = [x / 255.0 for x in c]
                    else:
                        c = CPK_COLORS_PV.get(
                            sym, [0.8, 0.8, 0.8]
                        )  # Default grey if missing
                    atom_colors.append(c)

            if positions:
                positions = np.array(positions)
                radii = np.array(radii)
                atom_colors = np.array(atom_colors)

                # --- Generate Merged Surface (SDF) ---

                # 1. Define Grid Bounds
                padding = radii.max() + 1.0
                min_bounds = positions.min(axis=0) - padding
                max_bounds = positions.max(axis=0) + padding

                # Resolution (voxel size in Angstroms)
                res_val = _vdw_settings.get("resolution", 0.125)
                # Clamp to safe limits just in case
                if res_val < 0.01:
                    res_val = 0.01
                spacing = (res_val, res_val, res_val)

                dims = np.ceil((max_bounds - min_bounds) / spacing).astype(int)
                nx, ny, nz = dims

                # Generate 1D coordinates for axes
                xs = np.arange(nx) * spacing[0] + min_bounds[0]
                ys = np.arange(ny) * spacing[1] + min_bounds[1]
                zs = np.arange(nz) * spacing[2] + min_bounds[2]

                values = np.full((nx, ny, nz), np.inf)

                # Calculate SDF values by iterating over atoms and using 3D broadcasting
                for i in range(len(positions)):
                    dx2 = (xs - positions[i, 0]) ** 2
                    dy2 = (ys - positions[i, 1]) ** 2
                    dz2 = (zs - positions[i, 2]) ** 2
                    dist = np.sqrt(
                        dx2[:, np.newaxis, np.newaxis]
                        + dy2[np.newaxis, :, np.newaxis]
                        + dz2[np.newaxis, np.newaxis, :]
                    ) - radii[i]
                    np.minimum(values, dist, out=values)

                grid = pv.ImageData()
                grid.dimensions = dims
                grid.origin = min_bounds
                grid.spacing = spacing

                # Assign computed distance field values using Fortran ordering to match VTK dimensions
                grid.point_data["values"] = values.ravel(order="F")

                # 3. Contour at iso-value 0 to get the surface
                mesh = grid.contour([0], scalars="values")

                # 4. Map Colors to Surface Vertices
                mesh_points = mesh.points
                n_mesh_pts = mesh_points.shape[0]
                if n_mesh_pts > 0:
                    min_dists = np.full(n_mesh_pts, np.inf)
                    nearest_atom_indices = np.zeros(n_mesh_pts, dtype=int)

                    # Determine closest atom for each mesh vertex in a vectorized loop over atoms
                    for i in range(len(positions)):
                        dx = mesh_points[:, 0] - positions[i, 0]
                        dy = mesh_points[:, 1] - positions[i, 1]
                        dz = mesh_points[:, 2] - positions[i, 2]
                        dist = np.sqrt(dx*dx + dy*dy + dz*dz) - radii[i]

                        mask = dist < min_dists
                        min_dists[mask] = dist[mask]
                        nearest_atom_indices[mask] = i

                    mesh.point_data["AtomColors"] = atom_colors[nearest_atom_indices]

                    opacity = _vdw_settings.get("occupancy", 0.3)

                    # # [DIRECT ACCESS] to manager's plotter
                    v3d = getattr(mw, "view_3d_manager", None)
                    if v3d and v3d.plotter:
                        v3d.plotter.add_mesh(
                            mesh,
                            scalars="AtomColors",
                            rgb=True,
                            opacity=opacity,
                            smooth_shading=True,
                            specular=0.2,
                            name="vdw_overlay_mesh",
                        )

        except Exception as e:
            logging.exception("VDW Overlay Error: %s", e)


def open_settings(context):
    """Namespaced singleton launcher for the settings dialog."""
    win = context.get_window("main_panel")
    if win:
        win.show()
        win.raise_()
        win.activateWindow()
        return

    win = VDWConfigWindow(context)
    win.show()


def initialize(context):
    """
    MoleditPy Plugin Entry Point (V3.0)
    """
    global PLUGIN_CONTEXT
    PLUGIN_CONTEXT = context
    load_settings()

    # 1. Register the 3D Drawing Style
    context.register_3d_style("vdw_overlay", draw_vdw_overlay)


def run(mw):
    if hasattr(mw, "host"):
        mw = mw.host
    context = PLUGIN_CONTEXT
    if not context:
        return
    open_settings(context)
