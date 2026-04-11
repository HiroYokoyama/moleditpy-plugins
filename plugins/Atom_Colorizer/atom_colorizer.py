from PyQt6.QtWidgets import (QVBoxLayout, QPushButton, QLabel, QColorDialog, 
                             QMessageBox, QLineEdit, QGroupBox, QDialog)
from PyQt6.QtGui import QColor, QCloseEvent
from PyQt6.QtCore import Qt, QTimer
import traceback
import logging

# Try importing from the installed package first (pip package structure)
try:
    from moleditpy.utils.constants import CPK_COLORS_PV
except ImportError:
    # Fallback to local 'modules' if running from source or sys.path is set that way
    try:
        from modules.constants import CPK_COLORS_PV
    except ImportError:
        # Final fallback map
        CPK_COLORS_PV = {}

# Plugin Metadata
PLUGIN_NAME = "Atom Colorizer"
PLUGIN_VERSION = "2026.04.11"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Applies custom colors to atoms in the 3D viewer. Refactored for V3 API."

class AtomColorizerWindow(QDialog):
    """
    Plugin window for applying custom colors to atoms in the 3D scene.
    Refactored for MoleditPy V3.0 API.
    """
    def __init__(self, context):
        super().__init__(parent=context.get_main_window())
        self.context = context
        self._restore_measurement_mode = False
        self._restore_edit_mode = False
        self._forced_measurement_mode = False
        
        # Set window properties
        self.setModal(False)
        self.setWindowTitle(PLUGIN_NAME)
        self.setWindowFlags(Qt.WindowType.Window) 
        
        self.current_color = QColor(255, 0, 0) # Default red
        self.init_ui()
        
        # Register window for V3 lifecycle management
        self.context.register_window("main_panel", self)
        
        # Ensure 3D picking is active while this window is open.
        mw = self.context.get_main_window()
        if mw:
            try:
                if hasattr(mw, "edit_3d_manager"):
                    self._restore_measurement_mode = bool(getattr(mw.edit_3d_manager, "measurement_mode", False))
                    self._restore_edit_mode = bool(getattr(mw.edit_3d_manager, "is_3d_edit_mode", False))
                if hasattr(mw, "edit_3d_manager") and not getattr(mw.edit_3d_manager, "measurement_mode", False):
                    if hasattr(mw, "init_manager") and hasattr(mw.init_manager, "measurement_action"):
                        mw.init_manager.measurement_action.setChecked(True)
                    mw.edit_3d_manager.toggle_measurement_mode(True)
                    self._forced_measurement_mode = True
            except Exception as _e:
                logging.warning("[atom_colorizer.py:65] silenced: %s", _e)
            
        self.context.show_status_message("Atom Colorizer: 3D picking enabled.")

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
        self.resize(300, 400)

    def get_selection_from_viewer(self):
        """Get selected atom indices from the context (2D) and Edit3DManager (3D)."""
        # Fetch 2D indices from core API
        indices = set(self.context.get_selected_atom_indices())
        
        # Fetch 3D indices directly from Edit3DManager due to "JUST PLUGIN" constraint
        mw = self.context.get_main_window()
        if mw and hasattr(mw, "edit_3d_manager"):
            indices_3d = getattr(mw.edit_3d_manager, "selected_atoms_3d", set())
            if isinstance(indices_3d, (set, list)):
                indices.update(indices_3d)
            picked_for_measurement = getattr(mw.edit_3d_manager, "selected_atoms_for_measurement", [])
            if isinstance(picked_for_measurement, (set, list, tuple)):
                indices.update(int(i) for i in picked_for_measurement)
        
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

        mol = self.context.current_molecule
        if not mol:
            QMessageBox.warning(self, "Error", "No molecule loaded.")
            return

        try:
            controller = self.context.get_3d_controller()
            hex_color = self.current_color.name()
            
            for idx in target_indices:
                if 0 <= idx < mol.GetNumAtoms():
                    controller.set_atom_color(idx, hex_color)
            
            self.context.refresh_3d_view()
                    
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply color: {e}")
            traceback.print_exc()

    def reset_colors(self):
        mol = self.context.current_molecule
        if not mol:
            return

        try:
            controller = self.context.get_3d_controller()
            for i in range(mol.GetNumAtoms()):
                controller.set_atom_color(i, None)
            
            self.context.refresh_3d_view()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reset colors: {e}")

    def closeEvent(self, event: QCloseEvent):
        """Restore previous 3D interaction mode when the panel closes."""
        try:
            if self.sel_timer.isActive():
                self.sel_timer.stop()
        except Exception as _e:
            logging.warning("[atom_colorizer.py:216] silenced: %s", _e)

        mw = self.context.get_main_window()
        if mw and hasattr(mw, "edit_3d_manager"):
            try:
                if self._forced_measurement_mode:
                    mw.edit_3d_manager.clear_measurement_selection()
                    if hasattr(mw, "init_manager") and hasattr(mw.init_manager, "measurement_action"):
                        mw.init_manager.measurement_action.setChecked(self._restore_measurement_mode)
                    mw.edit_3d_manager.toggle_measurement_mode(self._restore_measurement_mode)
                if hasattr(mw, "ui_manager"):
                    mw.ui_manager.toggle_3d_edit_mode(self._restore_edit_mode)
            except Exception as _e:
                logging.warning("[atom_colorizer.py:229] silenced: %s", _e)
        super().closeEvent(event)

def launch(context):
    win = context.get_window("main_panel")
    if win:
        win.show()
        win.raise_()
        win.activateWindow()
        return
        
    win = AtomColorizerWindow(context)
    win.show()

def initialize(context):
    """MoleditPy Plugin Entry Point (V3.0)"""
    
    def save_handler():
        """Save color overrides to project file."""
        mw = context.get_main_window()
        v3d = getattr(mw, 'view_3d_manager', None)
        if not v3d or not hasattr(v3d, '_plugin_color_overrides'):
            return {}
        
        return {
            "atom_colors": {str(k): v for k, v in v3d._plugin_color_overrides.items()}
        }
    
    def load_handler(data):
        """Restore color overrides from project file."""
        if not data:
            return
        
        atom_colors = data.get("atom_colors", {})
        controller = context.get_3d_controller()
        for atom_idx_str, hex_color in atom_colors.items():
            try:
                controller.set_atom_color(int(atom_idx_str), hex_color)
            except Exception as _e:
                logging.warning("[atom_colorizer.py:268] silenced: %s", _e)
        context.refresh_3d_view()
    
    def on_document_reset():
        """Reset colors when a new document is created."""
        mw = context.get_main_window()
        v3d = getattr(mw, 'view_3d_manager', None)
        if v3d and hasattr(v3d, '_plugin_color_overrides'):
            v3d._plugin_color_overrides.clear()
        
        win = context.get_window("main_panel")
        if win:
            win.le_indices.clear()

    context.register_save_handler(save_handler)
    context.register_load_handler(load_handler)
    context.register_document_reset_handler(on_document_reset)
    context.show_status_message(f"{PLUGIN_NAME} Loaded.")


def run(mw):
    """Legacy entry point for older plugin loaders."""
    if hasattr(mw, "host"):
        mw = mw.host

    if not hasattr(mw, "plugin_manager"):
        return

    from moleditpy.plugins.plugin_interface import PluginContext

    context = PluginContext(mw.plugin_manager, PLUGIN_NAME)
    launch(context)
