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

__version__="2026.02.13"
__author__="HiroYokoyama"

PLUGIN_NAME = "Atom Colorizer"

class AtomColorizerWindow(QDialog):
    def __init__(self, main_window):
        super().__init__(parent=main_window)
        self.mw = main_window
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

        # Update the line edit
        sorted_indices = sorted(list(indices))
        new_text = ",".join(map(str, sorted_indices))
        if self.le_indices.text() != new_text:
             self.le_indices.setText(new_text)

    def _auto_update_selection(self):
        """Timer slot to auto-update selection."""
        if self.le_indices.hasFocus():
            return
        self.get_selection_from_viewer()

    def choose_color(self):
        c = QColorDialog.getColor(initial=self.current_color, title="Select Color")
        if c.isValid():
            self.current_color = c
            # Update button style
            self.btn_color.setStyleSheet(f"background-color: {c.name()}; color: {'black' if c.lightness() > 128 else 'white'};")

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

        if not self.mw.current_mol:
            QMessageBox.warning(self, "Error", "No molecule loaded.")
            return

        try:
            # Use the API to set atom colors
            hex_color = self.current_color.name()
            
            for idx in target_indices:
                if 0 <= idx < self.mw.current_mol.GetNumAtoms():
                    # Access via main_window_view_3d proxy
                    if hasattr(self.mw, 'main_window_view_3d'):
                        self.mw.main_window_view_3d.update_atom_color_override(idx, hex_color)
                    else:
                         # Fallback if unproxied (unlikely in this architecture)
                         pass
                    
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply color: {e}")
            traceback.print_exc()

    def reset_colors(self):
        if not self.mw.current_mol:
            return

        try:
            # Clear all color overrides using the API
            for i in range(self.mw.current_mol.GetNumAtoms()):
                if hasattr(self.mw, 'main_window_view_3d'):
                    self.mw.main_window_view_3d.update_atom_color_override(i, None)
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reset colors: {e}")


# Global reference to keep window alive
_atom_colorizer_window = None

def run(mw):
    global _atom_colorizer_window
    
    # Check if window already exists
    if _atom_colorizer_window is None:
        _atom_colorizer_window = AtomColorizerWindow(mw)
        # Handle cleanup when window is closed
        _atom_colorizer_window.finished.connect(lambda: _cleanup_window())
    
    _atom_colorizer_window.show()
    _atom_colorizer_window.raise_()
    _atom_colorizer_window.activateWindow()

def _cleanup_window():
    global _atom_colorizer_window
    _atom_colorizer_window = None


def initialize(context):
    """
    Register plugin save/load handlers for persistence.
    """
    mw = context.get_main_window()
    
    def save_handler():
        """Save color overrides to project file."""
        # _plugin_color_overrides is stored on the MainWindow instance by the API
        if not hasattr(mw, '_plugin_color_overrides'):
            return {}
        
        # Convert color overrides to JSON-serializable format
        return {
            "atom_colors": {str(k): v for k, v in mw._plugin_color_overrides.items()}
        }
    
    def load_handler(data):
        """Load color overrides from project file."""
        if not data:
            return
        
        atom_colors = data.get("atom_colors", {})
        
        # Restore color overrides using the API
        if hasattr(mw, 'main_window_view_3d'):
            for atom_idx_str, hex_color in atom_colors.items():
                try:
                    atom_idx = int(atom_idx_str)
                    mw.main_window_view_3d.update_atom_color_override(atom_idx, hex_color)
                except Exception as e:
                    print(f"Failed to restore color for atom {atom_idx_str}: {e}")
    
    def on_document_reset():
        """Reset colors when a new document is created."""
        if hasattr(mw, '_plugin_color_overrides'):
            mw._plugin_color_overrides.clear()
        global _atom_colorizer_window
        if _atom_colorizer_window:
            _atom_colorizer_window.le_indices.clear()

    # Register handlers
    context.register_save_handler(save_handler)
    context.register_document_reset_handler(on_document_reset)
    context.register_load_handler(load_handler)