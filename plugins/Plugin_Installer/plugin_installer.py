#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plugin installer
Checks for updates of moleditpy plugins.
"""

import sys
import os
import json
import urllib.request
import urllib.error
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QTableWidget, QTableWidgetItem, QPushButton, 
    QHeaderView, QMessageBox, QAbstractItemView, QApplication, QCheckBox,
    QLineEdit, QMenu
)
from PyQt6.QtCore import Qt, QUrl
from PyQt6.QtGui import QDesktopServices, QColor, QIcon
import importlib.metadata

import shutil
import tempfile

# --- Metadata ---
PLUGIN_NAME = "Plugin Installer"
PLUGIN_VERSION = "2026.01.11"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Checks for updates, installs new plugins, and allows manual reinstallation."

REMOTE_JSON_URL = "https://raw.githubusercontent.com/HiroYokoyama/moleditpy-plugins/refs/heads/main/explorer/plugins.json"
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "plugin_installer.json")

# Global flag to ensure startup check runs only once per session
_startup_check_performed = False

def load_settings():
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, 'r') as f:
                return json.load(f)
        except:
            pass
    # Default to empty so we can detect "first run" (missing config)
    return {}

def save_settings(settings):
    try:
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(settings, f)
    except Exception as e:
        print(f"Error saving plugin installer settings: {e}")

def initialize(context):
    """
    Called by PluginManager at startup.
    Checks for updates if enabled in settings.
    """
    global _startup_check_performed
    if _startup_check_performed:
        return
    _startup_check_performed = True

    settings = load_settings()
    mw = context.get_main_window()
    if not mw:
        return

    # Check if "check_at_startup" is explicitly configured
    # if "check_at_startup" not in settings:
    #     # First Run: Ask the user ONCE
    #     from PyQt6.QtCore import QTimer
    #     QTimer.singleShot(2000, lambda: ask_user_permission(mw))
    elif settings.get("check_at_startup"):
        # Configured True: Run check
        from PyQt6.QtCore import QTimer
        QTimer.singleShot(2000, lambda: perform_startup_check(mw))

def ask_user_permission(mw):
    """Ask user if they want to enable automatic updates."""
    reply = QMessageBox.question(
        mw, "Plugin Installer",
        "Do you want to enable automatic plugin update checks at startup?\n\n(You can change this later in the Plugin Installer)",
        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        QMessageBox.StandardButton.No
    )
    
    settings = load_settings()
    if reply == QMessageBox.StandardButton.Yes:
        settings["check_at_startup"] = True
        save_settings(settings)
        # Run the check now since they said yes
        perform_startup_check(mw)
    else:
        settings["check_at_startup"] = False
        save_settings(settings)

def perform_startup_check(mw):
    try:
        # Run silent check
        checker = PluginInstallerWindow(mw, auto_check=True)
        if checker.updates_found:
            # Ask user if they want to see updates or skip
            reply = QMessageBox.question(
                mw, "Plugin Installer", 
                "Updates are available for MoleditPy or installed plugins.\nDo you want to check them now?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.Yes
            )
            
            if reply == QMessageBox.StandardButton.Yes:
                checker.show()
                checker.raise_()
            else:
                checker.deleteLater()
        else:
            checker.deleteLater()
    except Exception as e:
        print(f"Plugin Installer Auto-run failed: {e}")

def run(main_window):
    """
    Entry point for the plugin.
    Display the Plugin Installer window.
    """
    dialog = PluginInstallerWindow(main_window)
    dialog.exec()

class PluginDetailsDialog(QDialog):
    def __init__(self, parent, name, author, version, description, dependencies, local_info, target_file):
        super().__init__(parent)
        self.setWindowTitle(f"{name} - Details")
        self.resize(400, 350)
        self.parent_installer = parent # Reference to PluginInstallerWindow
        self.plugin_name = name
        self.target_file = target_file
        
        layout = QVBoxLayout(self)
        
        # Details
        lbl_name = QLabel(f"<h2>{name}</h2>")
        lbl_author = QLabel(f"<b>Author:</b> {author}")
        lbl_ver = QLabel(f"<b>Version:</b> {version}")
        
        lbl_desc = QLabel(f"<b>Description:</b><br>{description}")
        lbl_desc.setWordWrap(True)
        
        layout.addWidget(lbl_name)
        layout.addWidget(lbl_author)
        layout.addWidget(lbl_ver)
        layout.addSpacing(10)
        layout.addWidget(lbl_desc)
        
        if dependencies:
            layout.addSpacing(10)
            lbl_dep_header = QLabel("<b>Dependencies:</b>")
            layout.addWidget(lbl_dep_header)
            
            installed_deps = []
            missing_deps = []
            
            for dep in dependencies:
                try:
                    importlib.metadata.distribution(dep)
                    installed_deps.append(dep)
                except:
                    missing_deps.append(dep)
            
            # Construct display text
            dep_text_parts = []
            if installed_deps:
                text = ", ".join(installed_deps)
                dep_text_parts.append(f"<span style='color:green'>{text} (Installed)</span>")
            
            if missing_deps:
                text = ", ".join(missing_deps)
                dep_text_parts.append(f"<span style='color:red'>{text} (Missing)</span>")
            
            full_text = "<br>".join(dep_text_parts)
            lbl_deps = QLabel(full_text)
            lbl_deps.setWordWrap(True)
            layout.addWidget(lbl_deps)

            # Single Install Button for ALL missing
            if missing_deps:
                self.missing_deps = missing_deps # Store for button
                btn_install_all = QPushButton(f"Copy install command ({len(missing_deps)})")
                btn_install_all.setStyleSheet("color: blue;")
                btn_install_all.clicked.connect(self.copy_install_all_command)
                layout.addWidget(btn_install_all)

        layout.addStretch()
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        # Delete Button (Only if installed)
        if local_info:
            btn_delete = QPushButton("Delete Plugin")
            btn_delete.setStyleSheet("color: red;")
            btn_delete.clicked.connect(self.on_delete)
            btn_layout.addWidget(btn_delete)
            
        btn_layout.addStretch()
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btn_close.setDefault(True)
        btn_close.setFocus()
        btn_layout.addWidget(btn_close)
        
        layout.addLayout(btn_layout)

    def copy_install_all_command(self):
        if not hasattr(self, 'missing_deps') or not self.missing_deps:
            return
            
        deps_str = " ".join(self.missing_deps)
        cmd = f"pip install {deps_str}"
        QApplication.clipboard().setText(cmd)
        QMessageBox.information(self, "Install Command", 
                                f"Command copied to clipboard:\n\n{cmd}\n\nPlease close MoleditPy and run this terminal command.")

    def on_delete(self):
        # We can reuse the logic in the parent window, but we need to trigger it carefully
        # Or we can just Implement simple delete here using parent's method logic
        
        # Confirmation
        ret = QMessageBox.question(self, "Delete Plugin", 
                                   f"Are you sure you want to delete '{self.plugin_name}'?\nThis cannot be undone.",
                                   QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if ret != QMessageBox.StandardButton.Yes:
            return

        try:
            # Determine what to delete
            if self.target_file and os.path.exists(self.target_file):
                if os.path.basename(self.target_file) == "__init__.py":
                    dir_to_delete = os.path.dirname(self.target_file)
                    shutil.rmtree(dir_to_delete)
                else:
                    os.remove(self.target_file)
            else:
                 QMessageBox.warning(self, "Error", "Plugin file not found.")
                 return

            QMessageBox.information(self, "Success", f"Deleted '{self.plugin_name}'.")
            
            # Close details dialog
            self.accept()
            
            # Refresh parent table
            self.parent_installer.main_window.plugin_manager.discover_plugins(self.parent_installer.main_window)
            self.parent_installer.check_updates()
            
        except Exception as e:
             QMessageBox.warning(self, "Delete Error", f"Failed to delete plugin:\n{e}")

class PluginInstallerWindow(QDialog):
    def __init__(self, main_window, auto_check=False):
        super().__init__(main_window)
        self.setWindowTitle("Plugin Installer")
        self.resize(800, 500)
        self.main_window = main_window
        self.remote_data = []
        self.settings = load_settings()
        self.auto_check_mode = auto_check
        self.updates_found = False
        
        self.init_ui()
        
        if self.auto_check_mode:
            self.check_updates_silent()
        else:
            self.check_updates()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # Header Info
        try:
            from moleditpy.modules.constants import VERSION as APP_VERSION
        except ImportError:
            try:
                from moleditpy_linux.modules.constants import VERSION as APP_VERSION
            except ImportError:
                # Fallback for Linux/other environments: Try adding main script dir to sys.path
                try:
                    sys.path.append(os.path.dirname(os.path.abspath(sys.argv[0])))
                    try:
                        from moleditpy.modules.constants import VERSION as APP_VERSION
                    except ImportError:
                        from moleditpy_linux.modules.constants import VERSION as APP_VERSION
                except ImportError:
                    APP_VERSION = "Unknown"

        # Check PyPI for latest version
        self.latest_app_version = "Checking..."
        
        info_layout = QHBoxLayout()
        
        self.lbl_current_version = QLabel(f"<b>MoleditPy Version:</b> {APP_VERSION}")
        self.lbl_current_version.setStyleSheet("font-size: 14px;")
        info_layout.addWidget(self.lbl_current_version)
        
        info_layout.addSpacing(20)
        
        self.lbl_latest_version = QLabel(f"<b>Latest (PyPI):</b> {self.latest_app_version}")
        self.lbl_latest_version.setStyleSheet("font-size: 14px; color: gray;")
        info_layout.addWidget(self.lbl_latest_version)
        
        # Upgrade Button (Hidden by default)
        self.btn_upgrade_app = QPushButton("Copy upgrade command")
        self.btn_upgrade_app.setVisible(False)
        self.btn_upgrade_app.setStyleSheet("color: blue; font-weight: bold;")
        self.btn_upgrade_app.clicked.connect(self.copy_upgrade_command)
        info_layout.addWidget(self.btn_upgrade_app)
        
        info_layout.addStretch()
        layout.addLayout(info_layout)
        
        # Add a small line
        line = QLabel()
        line.setFrameStyle(QLabel.Shape.HLine | QLabel.Shadow.Sunken)
        layout.addWidget(line)

        # Search Bar
        search_layout = QHBoxLayout()
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search plugins by name or author...")
        self.search_input.textChanged.connect(self.filter_plugins)
        search_layout.addWidget(QLabel("Search:"))
        search_layout.addWidget(self.search_input)
        layout.addLayout(search_layout)

        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(6) # Name, Author, Local Ver, Latest Ver, Status, Action
        self.table.setHorizontalHeaderLabels(["Plugin Name", "Author", "Local Ver", "Latest Ver", "Status", "Action"])
        
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)          # Name (Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)          # Author (Stretch)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents) # Local Ver
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents) # Latest Ver
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents) # Status
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents) # Action
        
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.cellDoubleClicked.connect(self.show_plugin_details)
        layout.addWidget(self.table)
        
        # Settings Checkbox
        self.chk_startup = QCheckBox("Check for updates at startup")
        self.chk_startup.setChecked(self.settings.get("check_at_startup", False))
        self.chk_startup.toggled.connect(self.save_settings_ui)
        layout.addWidget(self.chk_startup)
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        refresh_btn = QPushButton("Refresh")
        refresh_btn.clicked.connect(self.check_updates)
        btn_layout.addWidget(refresh_btn)
        
        btn_layout.addStretch()
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.on_close_clicked) # Changed connection
        btn_layout.addWidget(close_btn)
        
        layout.addLayout(btn_layout)

    def save_settings_ui(self):
        self.settings["check_at_startup"] = self.chk_startup.isChecked()
        save_settings(self.settings)

    def fetch_pypi_version(self):
        try:
            url = "https://pypi.org/pypi/moleditpy/json"
            with urllib.request.urlopen(url, timeout=5) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode('utf-8'))
                    return data['info']['version']
        except Exception as e:
            print(f"Error fetching PyPI version: {e}")
        return "Unknown"

    def fetch_remote_data(self):
        try:
            with urllib.request.urlopen(REMOTE_JSON_URL, timeout=5) as response:
                if response.status == 200:
                    # Use utf-8-sig to handle potential BOM
                    data = json.loads(response.read().decode('utf-8-sig'))
                    return data
        except Exception as e:
            if not self.auto_check_mode:
                QMessageBox.warning(self, "Network Error", f"Failed to fetch plugin data:\n{e}")
            return []
        return []

    def check_updates_silent(self):
        """Perform checking without updating UI labels aggressively, return True if updates found."""
        # Check App Version
        latest = self.fetch_pypi_version()
        try:
            from moleditpy.modules.constants import VERSION as APP_VERSION
        except ImportError:
            try:
                from moleditpy_linux.modules.constants import VERSION as APP_VERSION
            except ImportError:
                APP_VERSION = "0.0.0"
            
        if latest != "Unknown" and latest > APP_VERSION:
             self.updates_found = True

        self.remote_data = self.fetch_remote_data()
        self.populate_table(silent=True)
        
        if self.updates_found:
             self.setWindowTitle("Plugin Installer (Updates Available!)")

    def check_updates(self):
        # Update App Version
        self.lbl_latest_version.setText("<b>Latest (PyPI):</b> Checking...")
        QApplication.processEvents() 
        
        latest = self.fetch_pypi_version()
        self.latest_app_version = latest
        
        try:
            from moleditpy.modules.constants import VERSION as APP_VERSION
        except ImportError:
            try:
                from moleditpy_linux.modules.constants import VERSION as APP_VERSION
            except ImportError:
                APP_VERSION = "0.0.0"

        color = "gray"
        if latest != "Unknown":
            if latest > APP_VERSION:
                 color = "#dc3545" # Red
                 status_text = f"{latest} (Update Available)"
            elif latest == APP_VERSION:
                 color = "green"
                 status_text = f"{latest} (Up to date)"
            else:
                 color = "blue"
                 status_text = f"{latest} (Dev/Newer)"
        else:
             status_text = "Unknown"

        self.lbl_latest_version.setText(f"<b>Latest (PyPI):</b> <span style='color:{color}'>{status_text}</span>")

        self.remote_data = self.fetch_remote_data()
        self.populate_table()

    def populate_table(self, silent=False):
        self.table.setRowCount(0)
        self.updates_found = False
        
        # Create map of remote data
        remote_map = {}
        for entry in self.remote_data:
            name = entry.get('name', '')
            remote_map[name] = entry

        installed_plugins = self.main_window.plugin_manager.plugins
        installed_map = {}
        for p in installed_plugins:
            installed_map[p.get('name')] = p

        # Union of names to show everything (installed + available)
        all_names = sorted(list(set(remote_map.keys()) | set(installed_map.keys())))
        
        for name in all_names:
            remote_info = remote_map.get(name)
            local_info = installed_map.get(name)
            
            local_ver = "-"
            remote_ver = "Unknown"
            status = "Unknown"
            author = "Unknown"
            color = None
            
            is_installed = (local_info is not None)
            
            # Visibility Check: 
            # If NOT installed and remote says strictly invisible, skip it.
            # If installed, show it regardless of remote visibility.
            if remote_info and not is_installed:
                 if not remote_info.get('visible', True):
                      continue

            can_update = False
            can_download = False
            target_file = None
            
            if is_installed:
                local_ver = local_info.get('version', 'Unknown')
                target_file = local_info.get('filepath')

            if remote_info:
                remote_ver = remote_info.get('version', 'Unknown')
                # Prioritize remote author info
                author = remote_info.get('author', 'Unknown')
            
            # If author still unknown, try local
            if author == "Unknown" and is_installed:
                author = local_info.get('author', 'Unknown')

            if is_installed:
                if remote_info:
                    if self.latest_app_version != "Checking..." and self.latest_app_version != "Unknown":
                        try:
                            if self.latest_app_version > APP_VERSION:
                                self.updates_found = True
                                self.btn_upgrade_app.setVisible(True)
                                self.btn_upgrade_app.setText(f"Copy upgrade command (v{self.latest_app_version})")
                            else:
                                self.btn_upgrade_app.setVisible(False)
                        except:
                            pass
                    if local_ver != 'Unknown' and remote_ver != 'Unknown':
                        if remote_ver > local_ver:
                            status = "Update Available"
                            color = QColor("#f8d7da")
                            can_update = True
                            self.updates_found = True
                        elif local_ver == remote_ver:
                            status = "Up to date"
                            color = QColor("#d4edda")
                            can_update = True # Allow manual reinstall
                        else:
                            status = "Newer"
                            color = QColor("#cce5ff")
                            can_update = True # Allow manual reinstall
                    else:
                        status = "Version Mismatch"
                        can_update = True
                else:
                    status = "Not in Registry"
                    color = QColor("#e2e3e5")
            else:
                status = "Not Installed"
                can_download = True
                # Determine target file for new download
                if remote_info and 'downloadUrl' in remote_info:
                    d_url = remote_info['downloadUrl']
                    # Download to plugin root directory (flatten structure)
                    filename = os.path.basename(d_url)
                    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                    target_file = os.path.join(base_dir, filename)

            if silent and not (is_installed and can_update and status == "Update Available"):
                # If silent check, mostly we care about notifying updates.
                pass

            row = self.table.rowCount()
            self.table.insertRow(row)
            
            # Name, Author, Local Ver, Latest Ver, Status, Action
            self.table.setItem(row, 0, QTableWidgetItem(name))
            self.table.setItem(row, 1, QTableWidgetItem(author))
            self.table.setItem(row, 2, QTableWidgetItem(str(local_ver)))
            self.table.setItem(row, 3, QTableWidgetItem(str(remote_ver)))
            
            status_item = QTableWidgetItem(status)
            if color:
                status_item.setBackground(color)
            self.table.setItem(row, 4, status_item)
            
            # Action Button
            if remote_info and 'downloadUrl' in remote_info:
                label = ""
                if can_update:
                    if status == "Update Available":
                        label = "Update"
                    else:
                        label = "Reinstall"
                elif can_download:
                     label = "Install"
                
                if label:
                    btn_action = QPushButton(label)
                    btn_action.setProperty("plugin_name", name)
                    btn_action.setProperty("download_url", remote_info['downloadUrl'])
                    btn_action.setProperty("target_file", target_file)
                    
                    if remote_info:
                        dependencies = remote_info.get('dependencies', [])
                        btn_action.setProperty("dependencies", dependencies)
                    
                    # Store data for double click lookups if needed, or we just look up in map
                    btn_action.clicked.connect(self.on_update_clicked)
                    self.table.setCellWidget(row, 5, btn_action)
            else:
                self.table.setItem(row, 5, QTableWidgetItem(""))

        # Re-apply filter if needed
        if hasattr(self, 'search_input'):
            self.filter_plugins()

    def filter_plugins(self):
        text = self.search_input.text().lower()
        rows = self.table.rowCount()
        for r in range(rows):
            name_item = self.table.item(r, 0)
            author_item = self.table.item(r, 1)
            
            name = name_item.text().lower() if name_item else ""
            author = author_item.text().lower() if author_item else ""
            
            should_show = (text in name) or (text in author)
            self.table.setRowHidden(r, not should_show)



    def show_plugin_details(self, row, col):
        name_item = self.table.item(row, 0)
        if not name_item: return
        name = name_item.text()
        
        # Find info
        remote_info = None
        for entry in self.remote_data:
            if entry.get('name') == name:
                remote_info = entry
                break
        
        # Also check local
        local_info = None
        for p in self.main_window.plugin_manager.plugins:
            if p.get('name') == name:
                local_info = p
                break
        
        description = "No description available."
        author = "Unknown"
        version = "Unknown"
        dependencies = []
        target_file = None
        
        if remote_info:
            description = remote_info.get('description', description)
            author = remote_info.get('author', author)
            version = remote_info.get('version', version)
            dependencies = remote_info.get('dependencies', [])
        
        if local_info:
            if description == "No description available.":
                description = local_info.get('description', description)
            if author == "Unknown":
                author = local_info.get('author', author)
            local_ver = local_info.get('version', 'Unknown')
            target_file = local_info.get('filepath')
            
            if remote_info:
                version = f"Installed: {local_ver}\nLatest: {version}"
            else:
                version = f"Installed: {local_ver}"

        dialog = PluginDetailsDialog(self, name, author, version, description, dependencies, local_info, target_file)
        dialog.exec()

    def on_close_clicked(self):
        self.accept()

    def closeEvent(self, event):
        # Reload plugins on close (X button or Close button)
        if self.main_window and hasattr(self.main_window, 'plugin_manager'):
            self.main_window.plugin_manager.discover_plugins(self.main_window)

            # Hot Replace: Update Main Window Menu to reflect changes
            # if hasattr(self.main_window, 'update_plugin_menu'):
            #     plugin_menu = None
            #     try:
            #         if hasattr(self.main_window, 'menuBar'):
            #             menus = [a.menu() for a in self.main_window.menuBar().actions() if a.menu()]
            #             
            #             # Strategy 1: Find by Name "Plugin" (startswith to handle "Plugins")
            #             for menu in menus:
            #                 title = menu.title().replace('&', '')
            #                 if title.lower().startswith('plugin'):
            #                     plugin_menu = menu
            #                     break
            #             
            #             # Strategy 2: Find by Content (contains "Plugin Manager...")
            #             if not plugin_menu:
            #                 for menu in menus:
            #                     for action in menu.actions():
            #                         if action.text().replace('&', '') == "Plugin Manager...":
            #                             plugin_menu = menu
            #                             break
            #                     if plugin_menu: break
            #     except Exception:
            #         pass
            #     
            #     if not plugin_menu:
            #         # Fallback: Use dummy menu to ensure toolbars/global menus update
            #         # even if we can't update the specific Plugin menu list.
            #         plugin_menu = QMenu()
            #     
            #     try:
            #         self.main_window.update_plugin_menu(plugin_menu)
            #     except Exception as e:
            #         print(f"Failed to update plugin menu: {e}")

        # Refresh PluginManagerWindow if it is open (Sync with Manager Window logic)
        for widget in QApplication.topLevelWidgets():
            if type(widget).__name__ == 'PluginManagerWindow':
                if hasattr(widget, 'refresh_plugin_list'):
                    try:
                        widget.refresh_plugin_list()
                    except Exception as e:
                        print(f"Failed to refresh PluginManagerWindow: {e}")

        super().closeEvent(event)

    def copy_upgrade_command(self):
        package_name = "moleditpy"
        try:
            # Attempt to find the distribution name (PyPI package name) 
            # that provides the 'moleditpy' module.
            # 'packages_distributions' is new in Python 3.10
            try:
                from importlib.metadata import packages_distributions
                dists = packages_distributions()
                if "moleditpy" in dists:
                    # returns a list of dist names, pick the first one
                    package_name = dists["moleditpy"][0]
            except ImportError:
                # Fallback for Python < 3.10
                for dist in importlib.metadata.distributions():
                    try:
                        toplevels = (dist.read_text('top_level.txt') or '').split()
                        if "moleditpy" in toplevels:
                            package_name = dist.metadata['Name']
                            break
                    except Exception:
                        pass
        except Exception as e:
            print(f"Error detecting package name: {e}")

        cmd = f"pip install --upgrade {package_name}"
        QApplication.clipboard().setText(cmd)
        QMessageBox.information(self, "Upgrade Command", 
                                f"Command copied to clipboard:\n\n{cmd}\n\nPlease close MoleditPy and run this command in your terminal.")


    def on_update_clicked(self):
        btn = self.sender()
        if not btn: return
        
        plugin_name = btn.property("plugin_name")
        download_url = btn.property("download_url")
        target_file = btn.property("target_file")
        dependencies = btn.property("dependencies") or []
        missing_deps = []
        user_confirmed_intent = False
        
        # Check dependencies first
        if dependencies:
            
            missing_deps = []
            for dep in dependencies:
                try:
                    importlib.metadata.distribution(dep)
                except:
                    missing_deps.append(dep)
            
            if missing_deps:
                # Show Details Dialog instead of installing
                # We need to reconstruct the arguments for show_plugin_details logic or just call it if we have row/col
                # But we don't know the row here easily without searching.
                # Easier to just instantiate the dialog directly since we have the data
                
                # Fetch full info to pass to dialog
                remote_info = None
                for entry in self.remote_data:
                    if entry.get('name') == plugin_name:
                        remote_info = entry
                        break
                
                # Local info
                local_info = None
                for p in self.main_window.plugin_manager.plugins:
                    if p.get('name') == plugin_name:
                        local_info = p
                        break
                
                description = remote_info.get('description', "No description available.") if remote_info else ""
                author = remote_info.get('author', "Unknown") if remote_info else "Unknown"
                version = remote_info.get('version', "Unknown") if remote_info else "Unknown"
                
                if local_info:
                    local_ver = local_info.get('version', 'Unknown')
                    version = f"Installed: {local_ver}\nLatest: {version}"
                
                # Ask user: Install anyway or View Details?
                ret = QMessageBox.question(self, "Missing Dependencies", 
                                    f"The following dependencies are missing: {', '.join(missing_deps)}.\n\nThe plugin may not work without them.\nDo you want to install it anyway?",
                                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
                
                if ret == QMessageBox.StandardButton.No:
                    # User chose NO, so show details to let them fix it
                    dialog = PluginDetailsDialog(self, plugin_name, author, version, description, dependencies, local_info, target_file)
                    dialog.exec()
                    return
                # If Yes, proceed to standard install confirmation below
                user_confirmed_intent = True

        if not download_url:
            QMessageBox.warning(self, "Error", "Cannot update: Missing download URL.")
            return

        # Confirm
        # Check if installed to decide text
        is_installed = False
        if target_file and os.path.exists(target_file):
             is_installed = True
        
        # If user already confirmed "Install anyway" despite missing deps, skip this generic confirmation
        if not user_confirmed_intent:
            action_verb = "Update" if is_installed else "Install"
            ret = QMessageBox.question(self, f"{action_verb} Plugin", 
                                       f"{action_verb} '{plugin_name}'?\nThis will download and install the plugin.",
                                       QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if ret != QMessageBox.StandardButton.Yes:
                return
        else:
            action_verb = "Update" if is_installed else "Install" # Need this for logging below

        # Resolve URL (handle relative)
        final_url = download_url
        if not download_url.startswith("http"):
            # Construct absolute URL relative to plugins.json location
            base_url = REMOTE_JSON_URL.rsplit('/', 1)[0]
            from urllib.parse import urljoin
            final_url = urljoin(base_url + "/", download_url)

        print(f"{action_verb}ing from: {final_url}")
        
        # Download to Temp File
        try:
            # Create a temporary directory to hold the download with its original name
            temp_dir = tempfile.mkdtemp()
            
            try:
                # 1. Determine filename
                filename = None
                
                # If we are updating an existing single-file plugin, FORCE the usage of the existing filename
                # so that install_plugin overwrites it exactly.
                if target_file and os.path.exists(target_file):
                    file_name_existing = os.path.basename(target_file)
                    # Check if it's a standard .py script (not a package __init__.py)
                    if file_name_existing != "__init__.py" and file_name_existing.endswith(".py"):
                        # If the update is also a .py file
                        if not final_url.lower().endswith(".zip"):
                            filename = file_name_existing
                
                # If valid filename not set by existing target, derive from URL
                if not filename:
                    path = urllib.parse.urlparse(final_url).path
                    filename = os.path.basename(path)
                
                if not filename:
                    # Fallback
                    if final_url.lower().endswith('.zip'):
                        filename = "plugin_update.zip"
                    else:
                        filename = "plugin_update.py"
                
                # Full path for the downloaded file
                download_path = os.path.join(temp_dir, filename)
                
                # Download
                with urllib.request.urlopen(final_url, timeout=15) as response:
                    if response.status == 200:
                        content = response.read()
                        with open(download_path, 'wb') as f:
                            f.write(content)
                        
                        # Install / Update Logic
                        try:
                            # Strict Update Mode:
                            # If we are updating an existing .py file (and the update is .py),
                            # we manually overwrite the TARGET file to preserve its location (e.g. inside a subfolder).
                            # We avoid using plugin_manager.install_plugin() here because it might force-install to the root plugins dir.
                            did_manual_overwrite = False
                            
                            if target_file and os.path.exists(target_file):
                                file_name_existing = os.path.basename(target_file)
                                # Check: Existing is .py (not __init__), Update is .py
                                if (file_name_existing != "__init__.py" and 
                                    file_name_existing.endswith(".py") and 
                                    not final_url.lower().endswith(".zip")):
                                    
                                    try:
                                        print(f"Overwriting existing plugin file: {target_file}")
                                        shutil.copy2(download_path, target_file)
                                        did_manual_overwrite = True
                                    except Exception as e:
                                        print(f"Manual overwrite failed: {e}. Falling back to manager.")
                            
                            if did_manual_overwrite:
                                # Manual overwrite success
                                pass
                            else:
                                # Standard Install (New, ZIP, or Transition)
                            
                                # Preservation Logic for settings.json (Folder directory plugins only)
                                # When installing a ZIP, the old directory is often wiped. We save settings.json to memory.
                                saved_settings_content = None
                                target_dir_for_restore = None

                                if target_file and os.path.exists(target_file):
                                    # Check if it is a folder plugin (target is __init__.py)
                                    if os.path.basename(target_file) == "__init__.py":
                                        folder_dir = os.path.dirname(target_file)
                                        settings_json = os.path.join(folder_dir, "settings.json")
                                        if os.path.exists(settings_json):
                                            try:
                                                with open(settings_json, 'r', encoding='utf-8') as f:
                                                    saved_settings_content = f.read()
                                                target_dir_for_restore = folder_dir
                                            except Exception as e:
                                                print(f"Plugin Installer: Failed to backup settings.json: {e}")

                                self.main_window.plugin_manager.install_plugin(download_path)
                                
                                # Restore settings.json
                                if saved_settings_content and target_dir_for_restore:
                                    # We assume the directory name is preserved or at least we try to write back to the same location
                                    # install_plugin replaces the directory.
                                    if os.path.exists(target_dir_for_restore):
                                        settings_json = os.path.join(target_dir_for_restore, "settings.json")
                                        try:
                                            with open(settings_json, 'w', encoding='utf-8') as f:
                                                f.write(saved_settings_content)
                                            print("Plugin Installer: Successfully restored settings.json")
                                        except Exception as e:
                                            print(f"Plugin Installer: Failed to restore settings.json: {e}")
                                
                                # Post-Install Cleanup (Transitions)
                                if target_file and os.path.exists(target_file):
                                    _, new_ext = os.path.splitext(filename)
                                    new_is_zip = new_ext.lower() == ".zip"
                                    old_is_py = target_file.lower().endswith(".py")
                                    old_is_init = os.path.basename(target_file) == "__init__.py"
                                    
                                    # Case 1: Transition from Single File (.py) to Package (.zip)
                                    if new_is_zip and old_is_py and not old_is_init:
                                        try:
                                            os.remove(target_file)
                                        except: pass
                                            
                                    # Case 2: Transition from Package (Folder) to Single File (.py)
                                    elif not new_is_zip and old_is_init:
                                        try:
                                            extracted_py = os.path.join(os.path.dirname(os.path.dirname(target_file)), filename)
                                            # If install_plugin extracted a single file, it's likely at root/plugins/filename
                                            # We need to be careful not to delete if things are messy, but standard logic applies:
                                            parent_dir = os.path.dirname(target_file)
                                            shutil.rmtree(parent_dir)
                                        except: pass
                            
                            # Success Message
                            # Success Message & Dependency Check
                            verb_past = "updated" if is_installed else "installed"
                            
                            if missing_deps:
                                # Reduce dialogs: Combine Success + Warning
                                QMessageBox.warning(self, "Installation Complete (Dependencies Missing)", 
                                                    f"Successfully {verb_past} '{plugin_name}'.\n\nHowever, you must install the missing dependencies for it to work.\n\nOpening details window...")
                                dialog = PluginDetailsDialog(self, plugin_name, author, version, description, dependencies, local_info, target_file)
                                dialog.exec()
                            else:
                                QMessageBox.information(self, "Success", f"Successfully {verb_past} '{plugin_name}'.")
                            
                            # Reload plugins
                            self.main_window.plugin_manager.discover_plugins(self.main_window)
                            # Refresh our list
                            self.check_updates()
                            
                        except Exception as e:
                            QMessageBox.warning(self, "Installation Error", f"Failed to install plugin:\n{e}")

                    else:
                        QMessageBox.warning(self, "Error", f"Download failed with status: {response.status}")
            
            finally:
                # Clean up temp directory
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
                    
        except Exception as e:
             QMessageBox.warning(self, "Update Error", f"Failed to update plugin:\n{e}")


