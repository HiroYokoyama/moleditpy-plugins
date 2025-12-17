import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
                             QLabel, QColorDialog, QDockWidget, QMessageBox, 
                             QLineEdit, QListWidget, QAbstractItemView, QGroupBox, QDialog)
from PyQt6.QtGui import QColor, QCloseEvent
from PyQt6.QtCore import Qt
import traceback
import sys
import os
import json
import functools # Added for polite patching

# Try importing from the installed package first (pip package structure)
try:
    from moleditpy.modules.constants import CPK_COLORS_PV
except ImportError:
    # Fallback to local 'modules' if running from source or sys.path is set that way
    try:
        from modules.constants import CPK_COLORS_PV
    except ImportError:
        # Final fallback map
        CPK_COLORS_PV = {}

__version__="2025.12.17"
__author__="HiroYokoyama"

PLUGIN_NAME = "Atom Colorizer"

class AtomColorizerWindow(QDialog):
    def __init__(self, main_window):
        super().__init__(parent=main_window)
        self.mw = main_window
        # self.dock = dock_widget # Removed as per instruction
        self.plotter = self.mw.plotter
        
        # Set window properties for modeless behavior
        self.setModal(False)
        self.setWindowTitle(PLUGIN_NAME)
        self.setWindowFlags(Qt.WindowType.Window) # Ensures it has min/max/close buttons
        
        # Initialize current_color as QColor object
        self.current_color = QColor(255, 0, 0) # Default red
        
        self.init_ui()
        
        # Auto-enable 3D Selection (Measurement Mode) if not already active
        try:
            if hasattr(self.mw, 'measurement_mode') and not self.mw.measurement_mode:
                if hasattr(self.mw, 'toggle_measurement_mode'):
                    self.mw.toggle_measurement_mode(True)
                    # Sync UI button state if possible
                    if hasattr(self.mw, 'measurement_action'):
                        self.mw.measurement_action.setChecked(True)
        except Exception as e:
            print(f"Failed to auto-enable 3D selection: {e}")

    def init_ui(self):
        layout = QVBoxLayout()
        
        # Information Label
        info_label = QLabel("Select atoms in the 3D viewer and apply color.")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)

        # Selection Group
        sel_group = QGroupBox("Selection")
        sel_layout = QVBoxLayout()
        
        self.le_indices = QLineEdit()
        self.le_indices.setPlaceholderText("e.g. 0, 1, 5")
        sel_layout.addWidget(self.le_indices)
        
        # 'Get Selection' button removed as per user request (auto-update is active)
        
        # Auto-update timer
        from PyQt6.QtCore import QTimer
        self.sel_timer = QTimer(self)
        self.sel_timer.timeout.connect(self._auto_update_selection)
        self.sel_timer.start(200) # Check every 200ms
        
        sel_group.setLayout(sel_layout)
        layout.addWidget(sel_group)

        # Color Group
        col_group = QGroupBox("Color")
        col_layout = QVBoxLayout()
        
        self.btn_color = QPushButton("Choose Color")
        # Initial style based on self.current_color (QColor object)
        self.btn_color.setStyleSheet(f"background-color: {self.current_color.name()}; color: {'black' if self.current_color.lightness() > 128 else 'white'};")
        self.btn_color.clicked.connect(self.choose_color)
        col_layout.addWidget(self.btn_color)
        
        btn_apply = QPushButton("Apply Color")
        btn_apply.clicked.connect(self.apply_color)
        col_layout.addWidget(btn_apply)
        
        col_group.setLayout(col_layout)
        layout.addWidget(col_group)

        # Reset Group
        reset_group = QGroupBox("Reset")
        reset_layout = QVBoxLayout()
        
        btn_reset = QPushButton("Reset to Element Colors")
        btn_reset.clicked.connect(self.reset_colors)
        reset_layout.addWidget(btn_reset)
        
        reset_group.setLayout(reset_layout)
        layout.addWidget(reset_group)
        
        layout.addStretch()
        
        # Close Button
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)

        self.setLayout(layout)
        
        # Resize window to a reasonable default
        self.resize(300, 400)

    def get_selection_from_viewer(self):
        """
        Get selected atom indices from the main window.
        Only checks 3D selection and Measurement selection. 2D selection is ignored per request.
        """
        indices = set()
        
        # 1. Check direct 3D selection (e.g. from 3D Drag or specific 3D select tools)
        if hasattr(self.mw, 'selected_atoms_3d') and self.mw.selected_atoms_3d:
            indices.update(self.mw.selected_atoms_3d)

        # 2. Check measurement selection (commonly used for picking atoms in 3D)
        if hasattr(self.mw, 'selected_atoms_for_measurement') and self.mw.selected_atoms_for_measurement:
             # selected_atoms_for_measurement might be list of int or objects, typically ints in this internal API
             for item in self.mw.selected_atoms_for_measurement:
                 if isinstance(item, int):
                     indices.add(item)

        # 2D Selection logic removed as per request ("2Dはいらない")

        if not indices:
            # Silent return if auto-updating, or maybe clear? 
            # If we invoke manually, we might want info, but generic message is okay if list is empty.
            pass

        # Update the line edit
        sorted_indices = sorted(list(indices))
        new_text = ",".join(map(str, sorted_indices))
        if self.le_indices.text() != new_text:
             self.le_indices.setText(new_text)

    def _auto_update_selection(self):
        """Timer slot to auto-update selection."""
        # Only update if the user is not actively typing? 
        # For now, just call get_selection_from_viewer which now checks for changes before setting text.
        # However, checking if le_indices has focus might be good.
        if self.le_indices.hasFocus():
            return
        self.get_selection_from_viewer()

    def choose_color(self):
        c = QColorDialog.getColor(initial=self.current_color, title="Select Color")
        if c.isValid():
            self.current_color = c
            # Update button style
            self.btn_color.setStyleSheet(f"background-color: {c.name()}; color: {'black' if c.lightness() > 128 else 'white'};")

    def _update_3d_actor(self):
        """Re-generate glyphs and update the actor to reflect color changes."""
        try:
            # 1. Re-run glyph filter to propagate color changes from glyph_source to mesh
            if hasattr(self.mw, 'glyph_source') and self.mw.glyph_source:
                # Read resolution from settings or default
                try:
                    style = self.mw.current_3d_style
                    if style == 'cpk':
                        resolution = self.mw.settings.get('cpk_resolution', 32)
                    elif style == 'stick':
                        resolution = self.mw.settings.get('stick_resolution', 16)
                    else: # ball_stick
                        resolution = self.mw.settings.get('ball_stick_resolution', 16)
                except Exception:
                    resolution = 16

                glyphs = self.mw.glyph_source.glyph(
                    scale='radii',
                    geom=pv.Sphere(radius=1.0, theta_resolution=resolution, phi_resolution=resolution),
                    orient=False
                )
                
                # 2. Update the actor
                if hasattr(self.mw, 'atom_actor') and self.mw.atom_actor:
                    self.mw.plotter.remove_actor(self.mw.atom_actor)
                
                # Re-add mesh (copying properties logic from main_window_view_3d roughly)
                is_lighting_enabled = self.mw.settings.get('lighting_enabled', True)
                mesh_props = dict(
                    smooth_shading=True,
                    specular=self.mw.settings.get('specular', 0.2),
                    specular_power=self.mw.settings.get('specular_power', 20),
                    lighting=is_lighting_enabled,
                )
                
                self.mw.atom_actor = self.mw.plotter.add_mesh(
                    glyphs, scalars='colors', rgb=True, **mesh_props
                )
                
                self.mw.plotter.render()
        except Exception as e:
            print(f"Error updating 3D actor: {e}")
            traceback.print_exc()

    def apply_color(self):
        txt = self.le_indices.text().strip()
        if not txt:
            QMessageBox.warning(self, "Warning", "No atoms selected. Please select atoms in the 3D viewer first.")
            return

        try:
            str_indices = [x.strip() for x in txt.split(',') if x.strip()]
            target_indices = [int(x) for x in str_indices]
        except ValueError:
            QMessageBox.warning(self, "Error", "Invalid indices format.")
            return

        if not hasattr(self.mw, 'glyph_source') or self.mw.glyph_source is None:
            QMessageBox.warning(self, "Error", "No 3D molecule found (glyph_source is None).")
            return

        
        # 1. Update glyph_source colors
        try:
            if not hasattr(self.mw, 'custom_atom_colors'):
                self.mw.custom_atom_colors = {}
                
            colors = self.mw.glyph_source.point_data['colors']
            
            # Helper to normalize color to whatever format 'colors' is using
            r, g, b = self.current_color.red(), self.current_color.green(), self.current_color.blue()
            
            # Store simple 0-255 list for persistence to avoid numpy type issues in JSON
            stored_color = [r, g, b]
            
            # Check if colors are float (0-1) or uint8 (0-255)
            is_float = (colors.dtype.kind == 'f')
            
            new_color_val = [r/255.0, g/255.0, b/255.0] if is_float else [r, g, b]

            for idx in target_indices:
                if 0 <= idx < len(colors):
                    colors[idx] = new_color_val
                    # Update persistent storage
                    self.mw.custom_atom_colors[idx] = stored_color
            
            # 2. Force update of the actor
            self._update_3d_actor()
            
        except Exception as e:
            # QMessageBox.critical(self, "Error", f"Failed to apply color: {e}")
            traceback.print_exc()

    def reset_colors(self):
        if not hasattr(self.mw, 'glyph_source') or self.mw.glyph_source is None:
            return
        if not self.mw.current_mol:
            return

        try:
            # Clear persistent storage
            if hasattr(self.mw, 'custom_atom_colors'):
                self.mw.custom_atom_colors = {}
                
            colors = self.mw.glyph_source.point_data['colors']
            is_float = (colors.dtype.kind == 'f')

            # Iterate atoms and reset to CPK
            for i in range(self.mw.current_mol.GetNumAtoms()):
                atom = self.mw.current_mol.GetAtomWithIdx(i)
                sym = atom.GetSymbol()
                # Get default color (float 0-1)
                base_col = CPK_COLORS_PV.get(sym, [0.5, 0.5, 0.5])
                
                if is_float:
                    colors[i] = base_col
                else:
                    colors[i] = [int(c*255) for c in base_col]
            
            self._update_3d_actor()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reset colors: {e}")


    def _restore_colors_from_file(self):
        """
        Check if the main window has a valid .pmeprj file open.
        If so, read it manually to find 'custom_atom_colors' and apply them.
        This handles the case where the file was loaded *before* this plugin started.
        """
        # If no file path or not a .pmeprj, ignore
        if not hasattr(self.mw, 'current_file_path') or not self.mw.current_file_path:
            return
        if not self.mw.current_file_path.lower().endswith('.pmeprj'):
            return
        
        # If we already have colors (unlikely if plugin just started, unless double-patched), skip
        if hasattr(self.mw, 'custom_atom_colors') and self.mw.custom_atom_colors:
            return

        try:
            with open(self.mw.current_file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            if "3d_structure" in data and data["3d_structure"]:
                raw_colors = data["3d_structure"].get("custom_atom_colors")
                if raw_colors:
                    custom_colors = {int(k): v for k, v in raw_colors.items()}
                    
                    # Apply to MainWindow
                    self.mw.custom_atom_colors = custom_colors
                    
                    # Force update of 3D actor
                    # We might need to ensure glyph_source is ready; assuming file load populated it.
                    if hasattr(self.mw, 'glyph_source') and self.mw.glyph_source:
                         # We need to manually inject these colors into the polydata
                         colors = self.mw.glyph_source.point_data['colors']
                         is_float = (colors.dtype.kind == 'f')

                         for idx, col_val in custom_colors.items():
                             if 0 <= idx < len(colors):
                                 if is_float:
                                     # If stored as 0-255 but buffer is float 0-1
                                     if any(c > 1.0 for c in col_val):
                                         colors[idx] = [c/255.0 for c in col_val]
                                     else:
                                         colors[idx] = col_val
                                 else:
                                      # If stored as float 0-1 but buffer is uint8
                                      if all(c <= 1.0 for c in col_val):
                                         colors[idx] = [int(c*255) for c in col_val]
                                      else:
                                         colors[idx] = col_val
                         
                         self._update_3d_actor()
                         print(f"Atom Colorizer: Restored {len(custom_colors)} custom colors from file.")
        except Exception as e:
            print(f"Atom Colorizer: Failed to lazy-load colors from file: {e}")
            traceback.print_exc()

# Global reference to keep window alive
_atom_colorizer_window = None
_patches_installed = False

def run(main_window):
    global _atom_colorizer_window, _patches_installed
    
    # Check if this is the first run (patches not installed)
    first_run = not _patches_installed

    # Install patches for persistence
    install_patches(main_window)

    # Check if window already exists
    if _atom_colorizer_window is None:
        _atom_colorizer_window = AtomColorizerWindow(main_window)
        # Handle cleanup when window is closed
        _atom_colorizer_window.finished.connect(lambda: _cleanup_window())
    
    # Only restore from file if this is the first execution
    if first_run:
        _atom_colorizer_window._restore_colors_from_file()

    _atom_colorizer_window.show()
    _atom_colorizer_window.raise_()
    _atom_colorizer_window.activateWindow()

def _cleanup_window():
    global _atom_colorizer_window
    _atom_colorizer_window = None

def install_patches(mw):
    """
    Install monkey patches to core modules to ensure color persistence.
    Checks `_patches_installed` to avoid double patching.
    """
    global _patches_installed
    if _patches_installed:
        return

    # Initialize persistent storage on MainWindow if not present
    if not hasattr(mw, 'custom_atom_colors'):
        mw.custom_atom_colors = {} 

    # --- Patch 1: MainWindowView3d.draw_molecule_3d ---
    # Purpose: Re-apply colors after any redraw (e.g. style change, molecular edit)
    # We patch the instance method on 'mw' which is the entry point for other modules.
    
    original_draw_3d = mw.draw_molecule_3d

    def patched_draw_3d(mol):
        # Call original
        res = original_draw_3d(mol)
        
        # Apply custom colors if they exist
        if hasattr(mw, 'custom_atom_colors') and mw.custom_atom_colors and hasattr(mw, 'glyph_source') and mw.glyph_source:
             try:
                import pyvista as pv # Ensure pyvista is available inside closure if needed
                
                colors = mw.glyph_source.point_data['colors']
                is_float = (colors.dtype.kind == 'f')
                
                # 1. Update the source colors
                for idx, col_val in mw.custom_atom_colors.items():
                    if isinstance(idx, str): idx = int(idx)
                    # Check index bounds
                    if 0 <= idx < len(colors):
                        if is_float:
                            if any(c > 1.0 for c in col_val):
                                colors[idx] = [c/255.0 for c in col_val]
                            else:
                                colors[idx] = col_val
                        else:
                             if all(c <= 1.0 for c in col_val):
                                colors[idx] = [int(c*255) for c in col_val]
                             else:
                                colors[idx] = col_val
                
                # 2. Re-generate the actor (Glyph filter)
                # Mimic common 3D view logic to respect resolution settings
                try:
                    style = getattr(mw, 'current_3d_style', 'cpk')
                    if style == 'cpk':
                        resolution = mw.settings.get('cpk_resolution', 32)
                    elif style == 'stick':
                        resolution = mw.settings.get('stick_resolution', 16)
                    else: # ball_stick
                        resolution = mw.settings.get('ball_stick_resolution', 16)
                except Exception:
                    resolution = 16

                glyphs = mw.glyph_source.glyph(
                    scale='radii',
                    geom=pv.Sphere(radius=1.0, theta_resolution=resolution, phi_resolution=resolution),
                    orient=False
                )
                
                # Remove old actor
                if hasattr(mw, 'atom_actor') and mw.atom_actor:
                    mw.plotter.remove_actor(mw.atom_actor)
                
                # Add new actor
                is_lighting_enabled = mw.settings.get('lighting_enabled', True)
                mesh_props = dict(
                    smooth_shading=True,
                    specular=mw.settings.get('specular', 0.2),
                    specular_power=mw.settings.get('specular_power', 20),
                    lighting=is_lighting_enabled,
                )
                
                mw.atom_actor = mw.plotter.add_mesh(
                    glyphs, scalars='colors', rgb=True, **mesh_props
                )
                
                # Force render
                if hasattr(mw, 'plotter'):
                    mw.plotter.render()

             except Exception as e:
                print(f"Patched draw_3d error: {e}")
                traceback.print_exc()
        return res

    mw.draw_molecule_3d = patched_draw_3d

    # --- Patch 2: MainWindowAppState.create_json_data ---
    # Purpose: Save colors to .pmeprj
    
    original_create_json = mw.create_json_data
    
    def patched_create_json():
        data = original_create_json()
        if hasattr(mw, 'custom_atom_colors') and mw.custom_atom_colors:
            if "3d_structure" in data and data["3d_structure"]:
                data["3d_structure"]["custom_atom_colors"] = mw.custom_atom_colors
        return data

    mw.create_json_data = patched_create_json
    
    # --- Patch 3: MainWindowAppState.load_from_json_data ---
    # Purpose: Load colors from .pmeprj
    
    original_load_json = mw.load_from_json_data
    
    def patched_load_json(json_data):
        # Extract custom colors
        custom_colors = {}
        if "3d_structure" in json_data and json_data["3d_structure"]:
             raw_colors = json_data["3d_structure"].get("custom_atom_colors")
             if raw_colors:
                 # Ensure keys are ints (JSON keys are strings)
                 custom_colors = {int(k): v for k, v in raw_colors.items()}
        
        # Set colors to mw BEFORE calling original logic 
        # (because original logic calls draw_molecule_3d, which uses our patch)
        mw.custom_atom_colors = custom_colors
        
        return original_load_json(json_data)
        
    mw.load_from_json_data = patched_load_json

    # --- Patch 4: MainWindow.clear_all (New / Clear) ---
    # Purpose: Reset colors when starting fresh
    
    original_clear_all = mw.clear_all
    
    @functools.wraps(original_clear_all)
    def patched_clear_all():
        # Reset colors BEFORE clearing logic
        if hasattr(mw, 'custom_atom_colors'):
            mw.custom_atom_colors = {}
        return original_clear_all()
        
    mw.clear_all = patched_clear_all
    
    # --- Patch 5: MainWindow.trigger_conversion (2D -> 3D) ---
    # Purpose: Reset colors when generating new 3D structure
    
    original_trigger_conversion = mw.trigger_conversion
    
    @functools.wraps(original_trigger_conversion)
    def patched_trigger_conversion():
        # Reset colors because structure is being regenerated
        if hasattr(mw, 'custom_atom_colors'):
            mw.custom_atom_colors = {}
        return original_trigger_conversion()
        
    mw.trigger_conversion = patched_trigger_conversion
    
    _patches_installed = True
    print("Atom Colorizer patches installed.")