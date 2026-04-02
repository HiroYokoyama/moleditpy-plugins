
import os
import json
import traceback
import pyvista as pv
import numpy as np
import functools
import types
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QLabel, 
                             QSlider, QHBoxLayout, QPushButton, QDoubleSpinBox)
from PyQt6.QtGui import QAction, QColor
from PyQt6.QtCore import Qt, QTimer

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
PLUGIN_VERSION = "2026.04.01"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Visualizes VDW radii as a translucent surface overlay using PyVista. Refactored for V3 API."

SETTINGS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vdw_radii_overlay.json")

# Global Settings (State)
_vdw_settings = {
    "occupancy": 0.3,   # Opacity (0.0 - 1.0)
    "resolution": 0.125 # Voxel spacing in Angstroms
}

def load_settings():
    global _vdw_settings
    try:
        if os.path.exists(SETTINGS_FILE):
            with open(SETTINGS_FILE, 'r') as f:
                saved = json.load(f)
                if "occupancy" in saved:
                    _vdw_settings["occupancy"] = float(saved["occupancy"])
                if "resolution" in saved:
                    _vdw_settings["resolution"] = float(saved["resolution"])
    except Exception as e:
        print(f"Error loading VDW settings: {e}")

def save_settings():
    try:
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(_vdw_settings, f, indent=4)
    except Exception as e:
        print(f"Error saving VDW settings: {e}")

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
        self.slider_res.setRange(5, 50) # 0.05 to 0.50
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

    def reset_defaults(self):
        # Default values
        def_occ = 0.3
        def_res = 0.125
        
        # Block signals to prevent redundant updates/saves during setting
        self.slider_occ.blockSignals(True)
        self.spin_occ.blockSignals(True)
        self.slider_res.blockSignals(True)
        self.spin_res.blockSignals(True)
        
        # Set values
        self.slider_occ.setValue(int(def_occ * 100))
        self.spin_occ.setValue(def_occ)
        self.slider_res.setValue(int(def_res * 100))
        self.spin_res.setValue(def_res)
        
        # Unblock
        self.slider_occ.blockSignals(False)
        self.spin_occ.blockSignals(False)
        self.slider_res.blockSignals(False)
        self.spin_res.blockSignals(False)
        
        # Update settings and view once
        _vdw_settings["occupancy"] = def_occ
        _vdw_settings["resolution"] = def_res
        save_settings()
        self.update_view()

    def refresh_ui_values(self):
        """Update UI elements from global settings."""
        occ = _vdw_settings.get("occupancy", 0.3)
        res = _vdw_settings.get("resolution", 0.125)
        
        self.slider_occ.blockSignals(True)
        self.spin_occ.blockSignals(True)
        self.slider_res.blockSignals(True)
        self.spin_res.blockSignals(True)

        self.slider_occ.setValue(int(occ * 100))
        self.spin_occ.setValue(occ)
        self.slider_res.setValue(int(res * 100))
        self.spin_res.setValue(res)

        self.slider_occ.blockSignals(False)
        self.spin_occ.blockSignals(False)
        self.slider_res.blockSignals(False)
        self.spin_res.blockSignals(False)

    def update_view(self):
        # Trigger redraw if we are in the correct mode
        # # [DIRECT ACCESS] via proxy
        mw = self.context.get_main_window()
        current_style = getattr(mw, 'current_3d_style', None)
        
        if current_style == "vdw_overlay":
            self.context.refresh_3d_view()

def draw_vdw_overlay(mw, mol):
    """
    3D Rendering Callback for VDW Overlay Style.
    Receives (mw, mol) from the engine.
    """
    # 1. Draw standard Ball & Stick first
    # # [DIRECT ACCESS] to manager
    if hasattr(mw, 'view_3d_manager'):
        mw.view_3d_manager.draw_standard_3d_style(mol, style_override='ball_and_stick')
    else:
        # Emergency fallback for legacy systems
        print("VDW Warning: view_3d_manager not found, using legacy draw path.")
        if hasattr(mw, 'view3d'):
             mw.view3d.draw_standard_3d_style(mol, style_override='ball_and_stick')

    # 2. Draw VDW Surface Overlay
    if mol and mol.GetNumAtoms() > 0:
        try:
            positions = []
            radii = []
            atom_colors = []
            
            # Use color overrides from the new namespaced dictionary on the manager
            custom_map = {}
            # # [DIRECT ACCESS] to manager state
            if hasattr(mw, 'view_3d_manager'):
                custom_map = getattr(mw.view_3d_manager, '_plugin_color_overrides', {})
            
            if not custom_map:
                # # [DIRECT ACCESS] legacy
                custom_map = getattr(mw, 'custom_atom_colors', {})
            
            if mol.GetNumConformers() > 0:
                conf = mol.GetConformer()
                for i in range(mol.GetNumAtoms()):
                    atom = mol.GetAtomWithIdx(i)
                    pos = conf.GetAtomPosition(i)
                    positions.append([pos.x, pos.y, pos.z])
                    
                    sym = atom.GetSymbol()
                    # GetRvdw returns standard VDW radius. 
                    r = 1.5 # Default
                    if pt:
                        r = pt.GetRvdw(atom.GetAtomicNum())
                    radii.append(r)
                    
                    # Color handling
                    if i in custom_map:
                        val = custom_map[i]
                        # Handling new API (Hex string) vs Legacy (List/Tuple)
                        if isinstance(val, str) and val.startswith('#'):
                            # Convert Hex to RGB [0-1]
                            qc = QColor(val)
                            c = [qc.redF(), qc.greenF(), qc.blueF()]
                        else:
                            # Assume legacy list/tuple
                            c = val
                            # Normalize 0-255 to 0-1 if needed
                            if any(x > 1.0 for x in c):
                                c = [x/255.0 for x in c]
                    else:
                        c = CPK_COLORS_PV.get(sym, [0.8, 0.8, 0.8]) # Default grey if missing
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
                if res_val < 0.01: res_val = 0.01
                spacing = (res_val, res_val, res_val) 
                
                dims = np.ceil((max_bounds - min_bounds) / spacing).astype(int)
                
                grid = pv.ImageData()
                grid.dimensions = dims
                grid.origin = min_bounds
                grid.spacing = spacing
                
                grid_points = grid.points
                n_points = grid_points.shape[0]
                values = np.empty(n_points)
                
                # Process in chunks of 100k points
                chunk_size = 100000 
                for start_idx in range(0, n_points, chunk_size):
                    end_idx = min(start_idx + chunk_size, n_points)
                    chunk_pts = grid_points[start_idx:end_idx]
                    
                    d = np.linalg.norm(chunk_pts[:, np.newaxis, :] - positions[np.newaxis, :, :], axis=2)
                    d_surface = d - radii[np.newaxis, :]
                    values[start_idx:end_idx] = d_surface.min(axis=1)

                grid.point_data["values"] = values
                
                # 3. Contour at iso-value 0 to get the surface
                mesh = grid.contour([0], scalars="values")
                
                # 4. Map Colors to Surface Vertices
                mesh_points = mesh.points
                if mesh_points.shape[0] > 0:
                    n_mesh_pts = mesh_points.shape[0]
                    mesh_colors = np.zeros((n_mesh_pts, 3))
                    
                    chunk_size = 50000
                    for start_idx in range(0, n_mesh_pts, chunk_size):
                        end_idx = min(start_idx + chunk_size, n_mesh_pts)
                        chunk_pts = mesh_points[start_idx:end_idx]
                        
                        d_center = np.linalg.norm(chunk_pts[:, np.newaxis, :] - positions[np.newaxis, :, :], axis=2)
                        d_surface = d_center - radii[np.newaxis, :]
                        nearest_atom_indices = np.argmin(d_surface, axis=1)
                        mesh_colors[start_idx:end_idx] = atom_colors[nearest_atom_indices]
                    
                    mesh.point_data["AtomColors"] = mesh_colors
                
                    opacity = _vdw_settings.get("occupancy", 0.3)
                    
                    # Assume mw.plotter is available
                    # # [DIRECT ACCESS] to plotter
                    if hasattr(mw, 'plotter'):
                        mw.plotter.add_mesh(
                            mesh,
                            scalars="AtomColors", 
                            rgb=True,
                            opacity=opacity,
                            smooth_shading=True,
                            specular=0.2, 
                            name="vdw_overlay_mesh"
                        )

        except Exception as e:
            print(f"VDW Overlay Error: {e}")
            traceback.print_exc()

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
    load_settings()
    
    # 1. Register the 3D Drawing Style
    context.register_3d_style("vdw_overlay", draw_vdw_overlay)
    
def run(mw):
    if hasattr(mw, 'host'):
        mw = mw.host
    from moleditpy.plugins.plugin_interface import PluginContext
    context = PluginContext(mw.plugin_manager, PLUGIN_NAME)
    open_settings(context)

