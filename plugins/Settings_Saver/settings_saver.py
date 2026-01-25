import os
import json
import traceback
from PyQt6.QtWidgets import (
    QMessageBox, QInputDialog, QFileDialog, QApplication, QColorDialog,
    QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QPushButton, QLabel,
    QWidget, QAbstractItemView, QMenu, QGroupBox, QCheckBox, QLineEdit,
    QListWidgetItem
)
from PyQt6.QtGui import QColor, QAction
from PyQt6.QtCore import Qt

PLUGIN_NAME = "Settings Saver"
PLUGIN_VERSION = "2026.01.25"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Save, load, and manage settings presets in a unified dialog."

SETTINGS_FILENAME = "settings_saver.json"

# Global state
PLUGIN_CONTEXT = None
PROJECT_PRESETS = {}
EMBED_SETTINGS = {"enabled": False} 
ORIGINAL_SETTINGS = None
PLUGIN_CONFIG = {"always_save_to_project": False}

def initialize(context):
    """Initialize the Settings Saver plugin."""
    global PLUGIN_CONTEXT, PLUGIN_CONFIG
    PLUGIN_CONTEXT = context
    
    # Load configuration on launch
    library = load_library()
    if "_PLUGIN_CONFIG" in library:
        PLUGIN_CONFIG.update(library["_PLUGIN_CONFIG"])
    
    
    # Apply "Always Save to Project" preference
    if PLUGIN_CONFIG.get("always_save_to_project", False):
        if PLUGIN_CONTEXT:
            mw = PLUGIN_CONTEXT.get_main_window()
            if not mw: mw = QApplication.activeWindow()
            if mw:
                enable_project_mode(mw)
    
    # Menu Action under "Settings"
    context.add_menu_action("Settings/Presets...", lambda: open_manager(context))

    # Register Save Handler (for .pmeprj)
    context.register_save_handler(on_save_project)

    # Register Load Handler (for .pmeprj)
    context.register_load_handler(on_load_project)

    # Register Reset Handler (File -> New)
    context.register_document_reset_handler(on_document_reset)

# --- Project Mode Logic ---

def enable_project_mode(mw):
    """Enable strict Project Mode protection."""
    global EMBED_SETTINGS
    EMBED_SETTINGS["enabled"] = True
    
    # 1. Reset Global Dirty Flag
    if hasattr(mw, 'settings_dirty'):
        mw.settings_dirty = False
    
    # 2. Alias initial_settings (Clean Exit Check)
    if hasattr(mw, 'initial_settings') and hasattr(mw, 'settings'):
        mw.initial_settings = mw.settings
        
    # 3. Monkey-patch save_settings (Strict Protection)
    # Even if User changes settings manually -> dirty=True.
    # We must BLOCK the actual save to file.
    if hasattr(mw, 'save_settings') and not hasattr(mw, '_original_save_settings'):
        mw._original_save_settings = mw.save_settings
        
        # Define proxy
        def proxy_save_settings():
            print("[Settings Saver] Global save blocked by Project Settings mode.")
            return # Block it!
            
        mw.save_settings = proxy_save_settings
        print("[Settings Saver] Enforcing strict Project Mode: Global saving blocked.")

def disable_project_mode(mw, restore_content=True):
    """Disable Project Mode and restore global save functionality."""
    global EMBED_SETTINGS, ORIGINAL_SETTINGS
    EMBED_SETTINGS["enabled"] = False
    
    # 1. Restore Original Save Function
    if hasattr(mw, '_original_save_settings'):
        mw.save_settings = mw._original_save_settings
        del mw._original_save_settings
        print("[Settings Saver] Project Mode disabled: Global saving restored.")
        
    # 2. Restore Original Settings Content & Alias
    if restore_content and ORIGINAL_SETTINGS is not None:
        if hasattr(mw, 'initial_settings'):
            mw.initial_settings = ORIGINAL_SETTINGS.copy()
            
        if hasattr(mw, 'settings'):
            try:
                mw.settings.clear()
                mw.settings.update(ORIGINAL_SETTINGS)
                apply_settings_hot(mw)
            except Exception as e:
                print(f"[Settings Saver] Error restoring original settings: {e}")
                
        ORIGINAL_SETTINGS = None

def on_save_project():
    """Return data to be saved into .pmeprj."""
    # Only save if enabled
    if not EMBED_SETTINGS["enabled"]:
        return None
    
    # Get Main Window to access settings
    mw = None
    if PLUGIN_CONTEXT:
        mw = PLUGIN_CONTEXT.get_main_window()
    if not mw:
        mw = QApplication.activeWindow()

    if mw and hasattr(mw, 'settings'):
        # Re-enforce protections just in case
        enable_project_mode(mw)
            
        return {
            "preset_name": "Project Settings", # Hardcoded name
            "settings": mw.settings
        }
    return None

def on_load_project(data):
    """Restore data from .pmeprj."""
    global PROJECT_PRESETS, EMBED_SETTINGS, ORIGINAL_SETTINGS
    
    if not isinstance(data, dict) or "settings" not in data:
        return

    # Helper: Always use "Project Settings" as the key/display name
    preset_name = "Project Settings" 
    project_settings = data.get("settings", {})

    # Store in global project presets
    PROJECT_PRESETS[preset_name] = project_settings
    
    if PLUGIN_CONTEXT:
        mw = PLUGIN_CONTEXT.get_main_window()
        if mw and hasattr(mw, 'settings'):
            try:
                # BACKUP ORIGINAL SETTINGS if not already backed up
                if ORIGINAL_SETTINGS is None:
                    ORIGINAL_SETTINGS = mw.settings.copy()

                mw.settings.update(project_settings)
                apply_settings_hot(mw)
                
                # Auto-enable strict mode AND protections
                enable_project_mode(mw)
                    
                mw.statusBar().showMessage(f"Applied project settings.")
            except Exception as e:
                print(f"Error auto-applying project settings: {e}")

def on_document_reset():
    """Clear project presets when creating a new file."""
    global PROJECT_PRESETS, EMBED_SETTINGS, ORIGINAL_SETTINGS

    if PLUGIN_CONTEXT:
        mw = PLUGIN_CONTEXT.get_main_window()
        if mw:
            # Disable project mode (restores original function and settings)
            disable_project_mode(mw, restore_content=True)

    PROJECT_PRESETS.clear()
    
    # If Default Preference is ON, re-enable for the new document
    if PLUGIN_CONFIG.get("always_save_to_project", False):
        if PLUGIN_CONTEXT:
             mw = PLUGIN_CONTEXT.get_main_window()
             if mw:
                 enable_project_mode(mw)
    else:
        EMBED_SETTINGS["enabled"] = False

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
        # Mark settings dirty, UNLESS we are in Project Mode
        try:
            if EMBED_SETTINGS.get("enabled", False):
                # If enabled, logic is handled by enable_project_mode elsewhere, 
                # but ensure dirty flag is cleared here too just in case
                mw.settings_dirty = False
            else:
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

        # 6. Push Undo State (Requirement: "push undo when applied", Manual: "Call after modifying")
        try:
            if hasattr(mw, 'push_undo_state'):
                mw.push_undo_state()
        except Exception as e:
            print(f"Error pushing undo state: {e}")

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
        self.resize(600, 500) # Slightly taller
        
        self.library = load_library()
        self.init_ui()

    def init_ui(self):
        layout = QHBoxLayout(self)

        # --- Left: List of Presets ---
        left_layout = QVBoxLayout()
        left_layout.addWidget(QLabel("Saved Presets:"))
        
        self.preset_list = QListWidget()
        self.preset_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.preset_list.itemDoubleClicked.connect(self.on_load) # Double click to load
        left_layout.addWidget(self.preset_list)
        
        layout.addLayout(left_layout, stretch=2)

        # --- Right: Buttons & Options ---
        right_layout = QVBoxLayout()
        
        # Action Buttons
        self.btn_load = QPushButton("Load Preset")
        self.btn_load.clicked.connect(self.on_load)
        self.btn_load.setToolTip("Load selected preset and apply immediately.")
        right_layout.addWidget(self.btn_load)

        self.btn_save = QPushButton("Save New...")
        self.btn_save.clicked.connect(self.on_save)
        self.btn_save.setToolTip("Save current settings as a new global preset.")
        right_layout.addWidget(self.btn_save)
        
        # New: Overwrite Global Default
        self.btn_set_global = QPushButton("Set as Global Default")
        self.btn_set_global.clicked.connect(self.on_save_as_global_default)
        self.btn_set_global.setToolTip("Overwrite your main application 'settings.json' with current settings.")
        right_layout.addWidget(self.btn_set_global)
        
        self.btn_delete = QPushButton("Delete")
        self.btn_delete.clicked.connect(self.on_delete)
        right_layout.addWidget(self.btn_delete)
        
        right_layout.addSpacing(10)
        
        # Project Embedding Section
        project_group = QGroupBox("Project Settings")
        project_layout = QVBoxLayout()
        
        self.chk_embed = QCheckBox("Save Settings to Project (This Project Only)")
        self.chk_embed.setChecked(EMBED_SETTINGS["enabled"])
        self.chk_embed.toggled.connect(self.on_embed_toggled)
        project_layout.addWidget(self.chk_embed)
        
        # Enforce dirty flag logic on open
        if EMBED_SETTINGS["enabled"] and hasattr(self.main_window, 'settings_dirty'):
            self.main_window.settings_dirty = False
            # Ensure alias
            if hasattr(self.main_window, 'initial_settings'):
                self.main_window.initial_settings = self.main_window.settings

        # --- Plugin Preference: Always Save ---
        self.chk_always_save = QCheckBox("Default: Always enable for new projects")
        self.chk_always_save.setToolTip("If checked, 'Save Settings to Project' will be enabled by default upon launch.")
        self.chk_always_save.setChecked(PLUGIN_CONFIG.get("always_save_to_project", False))
        self.chk_always_save.toggled.connect(self.on_always_save_toggled)
        project_layout.addWidget(self.chk_always_save)
        
        project_group.setLayout(project_layout)
        right_layout.addWidget(project_group)
        
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
        
        # 1. User Library Presets (Standard)
        # Filter out _PLUGIN_CONFIG
        names = sorted([k for k in self.library.keys() if not k.startswith("_")])
        for name in names:
            item = QListWidgetItem(name)
            self.preset_list.addItem(item)
            
        # 2. Project Presets (Blue)
        project_names = sorted(list(PROJECT_PRESETS.keys()))
        for name in project_names:
            display_name = name
            item = QListWidgetItem(display_name)
            item.setForeground(QColor("blue"))
            item.setToolTip("This preset is embedded in the current project file.")
            item.setData(Qt.ItemDataRole.UserRole, "project") 
            self.preset_list.addItem(item)
    
    def on_always_save_toggled(self, checked):
        PLUGIN_CONFIG["always_save_to_project"] = checked
        self.library["_PLUGIN_CONFIG"] = PLUGIN_CONFIG
        save_library(self.library, self)

    def get_selected_name(self):
        items = self.preset_list.selectedItems()
        if not items:
            return None
        return items[0].text()

    def is_project_preset(self, item):
        return item.data(Qt.ItemDataRole.UserRole) == "project"

    # --- Event Handlers ---

    def on_embed_toggled(self, checked):
        if hasattr(self.main_window, 'settings_dirty'):
            if checked:
                # ENABLE PROJECT MODE (Strict Protection)
                enable_project_mode(self.main_window)
            else:
                # DISABLE PROJECT MODE (Restore)
                disable_project_mode(self.main_window, restore_content=False)
                # Note: We don't restore content here because user might just want to 
                # stop saving to project but keep current look.
                # However, we DO want to restore the save_settings function.
                self.main_window.settings_dirty = True

    def on_load(self):
        items = self.preset_list.selectedItems()
        if not items:
            QMessageBox.warning(self, "Warning", "Please select a preset to load.")
            return

        item = items[0]
        name = item.text()
        
        settings_data = None
        if self.is_project_preset(item):
            settings_data = PROJECT_PRESETS.get(name)
        else:
            settings_data = self.library.get(name)

        if settings_data and hasattr(self.main_window, 'settings'):
            self.main_window.settings.update(settings_data)
            apply_settings_hot(self.main_window)
            self.main_window.statusBar().showMessage(f"Preset '{name}' loaded.")
        else:
            QMessageBox.warning(self, "Error", "Could not load settings.")

    def on_save_as_global_default(self):
        """Overwrite the main application's global settings.json."""
        reply = QMessageBox.question(
            self, "Confirm Overwrite", 
            "Are you sure you want to update the Global Default settings?\n"
            "This will overwrite 'settings.json' with the current configuration.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
            
        try:
            mw = self.main_window
            if not mw:
                 QMessageBox.warning(self, "Error", "Main window not found.")
                 return

            # Check if strict mode blocked the normal function
            save_func = None
            if hasattr(mw, '_original_save_settings'):
                save_func = mw._original_save_settings
            elif hasattr(mw, 'save_settings'):
                save_func = mw.save_settings
            
            if save_func:
                save_func() # Execute actual save
                
                # Also update our "Initial Settings" baseline so the app doesn't think we have diffs anymore
                if hasattr(mw, 'settings'):
                    global ORIGINAL_SETTINGS
                    if hasattr(mw, 'initial_settings'):
                         mw.initial_settings = mw.settings.copy()
                    # Update our plugin's backup too so "Reset" goes back to this new default
                    ORIGINAL_SETTINGS = mw.settings.copy()
                
                QMessageBox.information(self, "Success", "Global settings has been updated.")
            else:
                QMessageBox.warning(self, "Error", "Could not find save_settings function.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save global settings:\n{e}")

    def on_save(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Enter name for new global preset:")
        if not ok or not name.strip():
            return
        name = name.strip()
        
        if name.startswith("_"):
             QMessageBox.warning(self, "Invalid Name", "Preset names cannot start with underscore '_'.")
             return
        
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
            for it in items:
                if not self.is_project_preset(it):
                    self.preset_list.setCurrentItem(it)
                    break 
            QMessageBox.information(self, "Saved", f"Global preset '{name}' saved.")
        else:
            QMessageBox.warning(self, "Error", "Main window missing 'settings'.")

    def on_delete(self):
        items = self.preset_list.selectedItems()
        if not items:
            return
        item = items[0]
        name = item.text()
        
        if self.is_project_preset(item):
            QMessageBox.information(self, "Info", "Cannot delete project presets here. They are loaded from the .pmeprj file.")
            return

        reply = QMessageBox.question(self, "Confirm Delete", f"Are you sure you want to delete '{name}'?",
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if reply == QMessageBox.StandardButton.Yes:
            del self.library[name]
            save_library(self.library, self)
            self.refresh_list()

    def on_export_single(self):
        items = self.preset_list.selectedItems()
        if not items:
            QMessageBox.warning(self, "Warning", "Please select a preset to export.")
            return
        item = items[0]
        name = item.text()
        
        if self.is_project_preset(item):
            data = PROJECT_PRESETS.get(name)
        else:
            data = self.library.get(name)

        path, _ = QFileDialog.getSaveFileName(self, "Export Preset", f"{name}.json", "JSON Files (*.json)")
        if path:
            try:
                with open(path, 'w', encoding='utf-8') as f:
                    json.dump(data, f, indent=4)
                QMessageBox.information(self, "Exported", f"Exported '{name}' to {os.path.basename(path)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    def on_export_all(self):
        # Filter out _PLUGIN_CONFIG
        export_data = {k: v for k, v in self.library.items() if not k.startswith("_")}
        
        if not export_data:
            QMessageBox.information(self, "Export", "No presets to export.")
            return
            
        path, _ = QFileDialog.getSaveFileName(self, "Export All", "settings_library.json", "JSON Files (*.json)")
        if path:
            try:
                with open(path, 'w', encoding='utf-8') as f:
                    json.dump(export_data, f, indent=4)
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
            is_library = isinstance(first_val, dict)

            if is_library:
                count = 0
                for k, v in data.items():
                    if k.startswith("_"): continue # Skip config
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
                    if name.startswith("_"):
                         QMessageBox.warning(self, "Invalid Name", "Preset names cannot start with underscore '_'.")
                         return
                    self.library[name] = data
                    save_library(self.library, self)
                    self.refresh_list()
                    QMessageBox.information(self, "Imported", f"Imported '{name}'.")

        except Exception as e:
             QMessageBox.critical(self, "Error", f"Import failed:\n{e}")
