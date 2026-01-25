import os
import json
import traceback
from PyQt6.QtWidgets import (
    QMessageBox, QInputDialog, QFileDialog, QApplication, QColorDialog,
    QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QPushButton, QLabel,
    QWidget, QAbstractItemView, QMenu
)
from PyQt6.QtGui import QColor, QAction
from PyQt6.QtCore import Qt

PLUGIN_NAME = "Settings Saver"
PLUGIN_VERSION = "2026.01.25"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Save, load, and manage settings presets in a unified dialog."

SETTINGS_FILENAME = "settings_saver.json"

def initialize(context):
    """Initialize the Settings Saver plugin."""
    
    # Menu Action under "Settings"
    context.add_menu_action("Settings/Settings Saver...", lambda: open_manager(context))

    # Toolbar Button (Quick Access to Manager)
    context.add_toolbar_action(lambda: open_manager(context), "Settings Manager", icon=None, tooltip="Manage Settings Presets")

def open_manager(context):
    """Open the SettingsManager dialog."""
    mw = context.get_main_window()
    dialog = SettingsSaverDialog(context, mw)
    dialog.exec()

# --- Helper Functions ---

def get_plugin_data_path():
    """Get path to the storage JSON file in the plugin directory."""
    plugin_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(plugin_dir, SETTINGS_FILENAME)

def load_library():
    """Load the full library of presets from disk."""
    path = get_plugin_data_path()
    if not os.path.exists(path):
        return {}
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading settings library: {e}")
        return {}

def save_library(data, parent_window=None):
    """Save the full library of presets to disk."""
    path = get_plugin_data_path()
    try:
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=4)
    except Exception as e:
        if parent_window:
            QMessageBox.critical(parent_window, "Error", f"Failed to save settings library:\n{e}")

# --- Logic for Hot Loading (Reused) ---

def apply_settings_hot(mw):
    """Apply settings immediately to the active window (Hot Loading)."""
    try:
        # Mark settings dirty
        try:
            mw.settings_dirty = True
        except Exception:
            pass

        # 1. Apply 3D Settings
        try:
            if hasattr(mw, 'apply_3d_settings'):
                mw.apply_3d_settings()
        except Exception as e:
            print(f"Error applying 3D settings: {e}")

        # 2. Update CPK Colors
        try:
            if hasattr(mw, 'update_cpk_colors_from_settings'):
                mw.update_cpk_colors_from_settings()
        except Exception as e:
            print(f"Error updating CPK colors: {e}")

        # 3. Refresh Color Dialogs
        try:
            for w in QApplication.topLevelWidgets():
                if hasattr(w, 'refresh_ui'):
                    try: w.refresh_ui()
                    except: pass
        except: pass

        # 4. Redraw 3D Molecule
        try:
            if hasattr(mw, 'current_mol') and mw.current_mol:
                mw.draw_molecule_3d(mw.current_mol)
        except Exception as e:
            print(f"Error redrawing 3D molecule: {e}")

        # 5. Update 2D SCENE
        try:
            if hasattr(mw, 'scene') and mw.scene:
                bg_col_2d = mw.settings.get('background_color_2d', '#FFFFFF')
                mw.scene.setBackgroundBrush(QColor(bg_col_2d))
                for item in mw.scene.items():
                    if hasattr(item, 'update_style'):
                        item.update_style()
                    elif hasattr(item, 'update'):
                        item.update()
                if hasattr(mw, 'view_2d') and mw.view_2d and hasattr(mw.view_2d, 'viewport'):
                    mw.view_2d.viewport().update()
        except Exception as e:
            print(f"Error updating 2D settings: {e}")

    except Exception as e:
        print(f"Hot loading failed: {e}")
        traceback.print_exc()

# --- Main Dialog Class ---

class SettingsSaverDialog(QDialog):
    def __init__(self, context, parent=None):
        super().__init__(parent)
        self.context = context
        self.main_window = parent
        self.setWindowTitle("Settings Saver Manager")
        self.resize(500, 400)
        
        self.library = load_library()
        self.init_ui()

    def init_ui(self):
        layout = QHBoxLayout(self)

        # Left: List of Presets
        left_layout = QVBoxLayout()
        left_layout.addWidget(QLabel("Saved Presets:"))
        
        self.preset_list = QListWidget()
        self.preset_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.preset_list.itemDoubleClicked.connect(self.on_load) # Double click to load
        left_layout.addWidget(self.preset_list)
        
        layout.addLayout(left_layout, stretch=2)

        # Right: Buttons
        right_layout = QVBoxLayout()
        
        self.btn_load = QPushButton("Load Preset")
        self.btn_load.clicked.connect(self.on_load)
        self.btn_load.setToolTip("Load selected preset and apply immediately.")
        right_layout.addWidget(self.btn_load)

        self.btn_save = QPushButton("Save New...")
        self.btn_save.clicked.connect(self.on_save)
        self.btn_save.setToolTip("Save current settings as a new preset.")
        right_layout.addWidget(self.btn_save)
        
        self.btn_delete = QPushButton("Delete")
        self.btn_delete.clicked.connect(self.on_delete)
        right_layout.addWidget(self.btn_delete)
        
        right_layout.addStretch()
        
        # Export Menu Button
        self.btn_export = QPushButton("Export...")
        self.export_menu = QMenu(self)
        self.export_menu.addAction("Export Selected...", self.on_export_single)
        self.export_menu.addAction("Export All...", self.on_export_all)
        self.btn_export.setMenu(self.export_menu)
        right_layout.addWidget(self.btn_export)

        self.btn_import = QPushButton("Import...")
        self.btn_import.clicked.connect(self.on_import)
        right_layout.addWidget(self.btn_import)
        
        right_layout.addSpacing(20)
        
        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.accept)
        right_layout.addWidget(self.btn_close)

        layout.addLayout(right_layout, stretch=1)
        
        self.refresh_list()

    def refresh_list(self):
        self.preset_list.clear()
        names = sorted(list(self.library.keys()))
        self.preset_list.addItems(names)

    def get_selected_name(self):
        items = self.preset_list.selectedItems()
        if not items:
            return None
        return items[0].text()

    # --- Actions ---

    def on_load(self):
        name = self.get_selected_name()
        if not name:
            QMessageBox.warning(self, "Warning", "Please select a preset to load.")
            return

        settings_data = self.library.get(name)
        if hasattr(self.main_window, 'settings'):
            self.main_window.settings.update(settings_data)
            apply_settings_hot(self.main_window)
            self.main_window.statusBar().showMessage(f"Preset '{name}' loaded.")
        else:
            QMessageBox.warning(self, "Error", "Main window missing 'settings'.")

    def on_save(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Enter name for new preset:")
        if not ok or not name.strip():
            return
        name = name.strip()
        
        if name in self.library:
            reply = QMessageBox.question(self, "Overwrite?", f"'{name}' already exists. Overwrite?",
                                         QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if reply != QMessageBox.StandardButton.Yes:
                return

        if hasattr(self.main_window, 'settings'):
            self.library[name] = self.main_window.settings.copy()
            save_library(self.library, self)
            self.refresh_list()
            # Select the new item
            items = self.preset_list.findItems(name, Qt.MatchFlag.MatchExactly)
            if items:
                self.preset_list.setCurrentItem(items[0])
            QMessageBox.information(self, "Saved", f"Preset '{name}' saved.")
        else:
            QMessageBox.warning(self, "Error", "Main window missing 'settings'.")

    def on_delete(self):
        name = self.get_selected_name()
        if not name:
            return
        
        reply = QMessageBox.question(self, "Confirm Delete", f"Are you sure you want to delete '{name}'?",
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if reply == QMessageBox.StandardButton.Yes:
            del self.library[name]
            save_library(self.library, self)
            self.refresh_list()

    def on_export_single(self):
        name = self.get_selected_name()
        if not name:
            QMessageBox.warning(self, "Warning", "Please select a preset to export.")
            return
            
        data = self.library[name]
        path, _ = QFileDialog.getSaveFileName(self, "Export Preset", f"{name}.json", "JSON Files (*.json)")
        if path:
            try:
                with open(path, 'w', encoding='utf-8') as f:
                    json.dump(data, f, indent=4)
                QMessageBox.information(self, "Exported", f"Exported '{name}' to {os.path.basename(path)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    def on_export_all(self):
        if not self.library:
            QMessageBox.information(self, "Export", "No presets to export.")
            return
            
        path, _ = QFileDialog.getSaveFileName(self, "Export All", "settings_library.json", "JSON Files (*.json)")
        if path:
            try:
                with open(path, 'w', encoding='utf-8') as f:
                    json.dump(self.library, f, indent=4)
                QMessageBox.information(self, "Exported", f"All presets exported to {os.path.basename(path)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    def on_import(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import Presets", "", "JSON Files (*.json)")
        if not path:
            return
            
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            if not isinstance(data, dict):
                raise ValueError("Invalid JSON format.")

            # Heuristic: Library vs Single
            first_val = next(iter(data.values())) if data else None
            is_library = isinstance(first_val, dict) # If value is dict, it's a library {name: {settings}}

            if is_library:
                count = 0
                for k, v in data.items():
                    if isinstance(v, dict):
                        self.library[k] = v
                        count += 1
                save_library(self.library, self)
                self.refresh_list()
                QMessageBox.information(self, "Imported", f"Imported {count} presets.")
            else:
                # Single preset
                default_name = os.path.splitext(os.path.basename(path))[0]
                name, ok = QInputDialog.getText(self, "Import Preset", "Name for imported preset:", text=default_name)
                if ok and name.strip():
                    name = name.strip()
                    self.library[name] = data
                    save_library(self.library, self)
                    self.refresh_list()
                    QMessageBox.information(self, "Imported", f"Imported '{name}'.")

        except Exception as e:
             QMessageBox.critical(self, "Error", f"Import failed:\n{e}")
