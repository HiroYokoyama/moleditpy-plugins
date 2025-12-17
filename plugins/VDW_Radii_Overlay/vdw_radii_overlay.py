import os
import json
import traceback
import pyvista as pv
import numpy as np
import functools
import types
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QLabel, 
                             QSlider, QHBoxLayout, QPushButton)
from PyQt6.QtGui import QAction
from PyQt6.QtCore import Qt, QTimer

# Try to import VDW radii from constants, fallback if needed
try:
    from moleditpy.modules.constants import pt, CPK_COLORS_PV
except ImportError:
    try:
        from modules.constants import pt, CPK_COLORS_PV
    except ImportError:
        from rdkit import Chem
        pt = Chem.GetPeriodicTable()
        CPK_COLORS_PV = {} # Last resort fallback


# Plugin Metadata
__version__="2025.12.17"
__author__="HiroYokoyama"
PLUGIN_NAME = "VDW Radii Overlay"
SETTINGS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vdw_radii_overlay.json")

# Global State
_patches_installed = False
_config_window = None
_vdw_settings = {
    "occupancy": 0.3  # Opacity (0.0 - 1.0)
}

def load_settings():
    global _vdw_settings
    try:
        if os.path.exists(SETTINGS_FILE):
            with open(SETTINGS_FILE, 'r') as f:
                saved = json.load(f)
                # Filter to only keep occupancy
                if "occupancy" in saved:
                    _vdw_settings["occupancy"] = float(saved["occupancy"])
    except Exception as e:
        print(f"Error loading VDW settings: {e}")

def save_settings():
    try:
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(_vdw_settings, f, indent=4)
    except Exception as e:
        print(f"Error saving VDW settings: {e}")

class VDWConfigWindow(QDialog):
    def __init__(self, main_window):
        super().__init__(parent=main_window)
        self.mw = main_window
        self.setWindowTitle("VDW Overlay Settings")
        self.setModal(False)
        self.resize(300, 100)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        
        # Occupancy Slider
        occ_layout = QHBoxLayout()
        occ_layout.addWidget(QLabel("Occupancy:"))
        self.slider_occ = QSlider(Qt.Orientation.Horizontal)
        self.slider_occ.setRange(0, 100)
        current_occ = int(_vdw_settings.get("occupancy", 0.3) * 100)
        self.slider_occ.setValue(current_occ)
        self.slider_occ.valueChanged.connect(self.on_occupancy_changed)
        occ_layout.addWidget(self.slider_occ)
        self.lbl_occ_val = QLabel(f"{current_occ}%")
        occ_layout.addWidget(self.lbl_occ_val)
        layout.addLayout(occ_layout)
        
        # Close Button
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        layout.addWidget(btn_close)
        
        self.setLayout(layout)

    def on_occupancy_changed(self, value):
        opacity = value / 100.0
        _vdw_settings["occupancy"] = opacity
        self.lbl_occ_val.setText(f"{value}%")
        save_settings()
        self.update_view()

    def update_view(self):
        # Trigger redraw if we are in the correct mode
        if hasattr(self.mw, 'current_3d_style') and self.mw.current_3d_style == "vdw_overlay":
            if self.mw.current_mol:
                self.mw.draw_molecule_3d(self.mw.current_mol)

def install_patches(mw):
    global _patches_installed
    if _patches_installed:
        return

    # Patch draw_molecule_3d
    original_draw_3d = mw.draw_molecule_3d
    
    # Save original reference for other plugins or debugging
    if not hasattr(mw, '_vdw_overlay_original_draw_3d'):
        mw._vdw_overlay_original_draw_3d = original_draw_3d

    @functools.wraps(original_draw_3d)
    def patched_draw_3d(mol):
        # Check if we are in "vdw_overlay" style
        is_vdw_style = False
        if hasattr(mw, 'current_3d_style') and mw.current_3d_style == "vdw_overlay":
            is_vdw_style = True
            
            # Temporarily switch to 'ball_and_stick' for the base draw
            # We set the attribute directly to avoid triggering a recursive redraw if there's a setter
            mw.current_3d_style = 'ball_and_stick' 
        
        # Call the original draw (which will now draw ball_and_stick)
        res = original_draw_3d(mol)
        
        # If we were in vdw_overlay, draw the overlay and restore the style name
        if is_vdw_style:
            # Restore the style name so the UI and state remain consistent
            mw.current_3d_style = "vdw_overlay"
            
            if mol and mol.GetNumAtoms() > 0:
                try:
                    positions = []
                    radii = []
                    atom_colors = []
                    
                    # Use custom colors if available, otherwise CPK
                    custom_map = getattr(mw, 'custom_atom_colors', {})
                    
                    if mol.GetNumConformers() > 0:
                        conf = mol.GetConformer()
                        for i in range(mol.GetNumAtoms()):
                            atom = mol.GetAtomWithIdx(i)
                            pos = conf.GetAtomPosition(i)
                            positions.append([pos.x, pos.y, pos.z])
                            
                            sym = atom.GetSymbol()
                            # GetRvdw returns standard VDW radius. 
                            r = pt.GetRvdw(atom.GetAtomicNum())
                            radii.append(r)
                            
                            # Color handling
                            if i in custom_map:
                                c = custom_map[i]
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
                        # Padding = max radius + some buffer
                        padding = radii.max() + 1.0
                        min_bounds = positions.min(axis=0) - padding
                        max_bounds = positions.max(axis=0) + padding
                        
                        # Resolution (voxel size in Angstroms)
                        spacing = (0.125, 0.125, 0.125) 
                        
                        dims = np.ceil((max_bounds - min_bounds) / spacing).astype(int)
                        
                        grid = pv.ImageData()
                        grid.dimensions = dims
                        grid.origin = min_bounds
                        grid.spacing = spacing
                        
                        # Get grid points coordinates
                        # grid.points is (N, 3)
                        grid_points = grid.points
                        
                        # 2. Compute Signed Distance Field (SDF)
                        # Value = min( distance_to_atom_center - atom_radius )
                        # We want the surface where Value = 0
                        
                        n_points = grid_points.shape[0]
                        values = np.empty(n_points)
                        
                        # Process in chunks of 100k points
                        chunk_size = 100000 
                        for start_idx in range(0, n_points, chunk_size):
                            end_idx = min(start_idx + chunk_size, n_points)
                            chunk_pts = grid_points[start_idx:end_idx]
                            
                            # Distance matrix (Chunk x Atoms)
                            # dists[k, i] = distance from point k to atom i
                            d = np.linalg.norm(chunk_pts[:, np.newaxis, :] - positions[np.newaxis, :, :], axis=2)
                            
                            # Subtract radii
                            d_surface = d - radii[np.newaxis, :]
                            
                            # Min value for each point across all atoms
                            values[start_idx:end_idx] = d_surface.min(axis=1)

                        grid.point_data["values"] = values
                        
                        # 3. Contour at iso-value 0 to get the surface
                        mesh = grid.contour([0], scalars="values")
                        
                        # 4. Map Colors to Surface Vertices
                        # Using Point Data enables smooth interpolation (gradients) between colors
                        mesh_points = mesh.points
                        if mesh_points.shape[0] > 0:
                            # Use nearest neighbor search on VERTICES
                            
                            n_mesh_pts = mesh_points.shape[0]
                            mesh_colors = np.zeros((n_mesh_pts, 3))
                            
                            # Chunking for color mapping
                            chunk_size = 50000
                            for start_idx in range(0, n_mesh_pts, chunk_size):
                                end_idx = min(start_idx + chunk_size, n_mesh_pts)
                                chunk_pts = mesh_points[start_idx:end_idx]
                                
                                # Distances: (Chunk, Atoms)
                                # dist[k, i] = distance from point k center to atom i center
                                d_center = np.linalg.norm(chunk_pts[:, np.newaxis, :] - positions[np.newaxis, :, :], axis=2)
                                
                                # VDW-radii aware distance metric
                                d_surface = d_center - radii[np.newaxis, :]
                                
                                # Find index of atom with min surface distance
                                nearest_atom_indices = np.argmin(d_surface, axis=1)
                                
                                # Assign colors
                                mesh_colors[start_idx:end_idx] = atom_colors[nearest_atom_indices]
                            
                            mesh.point_data["AtomColors"] = mesh_colors
                        
                            opacity = _vdw_settings.get("occupancy", 0.3)
                            
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

        return res

    # Bind the patched method to the instance
    mw.draw_molecule_3d = patched_draw_3d 
    _patches_installed = True


def autorun(mw):
    """
    Called automatically by plugin system on startup.
    Injects 'VDW Overlay' into the 3D Style menu on the toolbar.
    """
    load_settings()
    install_patches(mw)
    
    # Use a timer to ensure the main window UI (toolbar) is fully initialized
    # before trying to modify the menu.
    
    def add_menu_item():
        # Locate the "3D Style" tool button on the main toolbar
        style_button = getattr(mw, 'style_button', None)
        
        if style_button and style_button.menu():
            style_menu = style_button.menu()
            
            # Check if already added
            exists = False
            for a in style_menu.actions():
                if a.text() == "VDW Overlay":
                    exists = True
                    break
            
            if not exists:
                # We want to add it to the same action group as the others so it behaves like a radio button
                existing_actions = style_menu.actions()
                group = None
                if existing_actions:
                    group = existing_actions[0].actionGroup()
                
                action = QAction("VDW Overlay", mw)
                action.setCheckable(True)
                
                if group:
                    group.addAction(action)
                
                # Connect to set_3d_style
                action.triggered.connect(lambda: mw.set_3d_style("vdw_overlay"))
                
                style_menu.addAction(action)
    
    # Run after 2000ms to be safe
    QTimer.singleShot(2000, add_menu_item)

def run(mw):
    """
    Called from Extensions menu. Shows configuration.
    """
    load_settings()
    install_patches(mw)
    
    global _config_window
    if _config_window is None:
        _config_window = VDWConfigWindow(mw)
        _config_window.finished.connect(lambda: _cleanup_config())
    
    _config_window.show()
    _config_window.raise_()
    _config_window.activateWindow()

def _cleanup_config():
    global _config_window
    _config_window = None
