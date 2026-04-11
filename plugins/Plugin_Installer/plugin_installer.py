#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plugin installer
Checks for updates of moleditpy plugins.
"""

import sys
import os
import json
import zipfile
import urllib.request
import urllib.error
import logging
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QTableWidget, QTableWidgetItem, QPushButton, 
    QHeaderView, QMessageBox, QAbstractItemView, QApplication, QCheckBox,
    QLineEdit, QMenu
)
from PyQt6.QtCore import Qt, QUrl
from PyQt6.QtGui import QDesktopServices, QColor, QIcon
import importlib.metadata
import hashlib
import re

import shutil
import tempfile

# --- Metadata ---
PLUGIN_NAME = "Plugin Installer"
PLUGIN_VERSION = "2026.04.12"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Checks for updates, installs new plugins, and allows manual reinstallation."

REMOTE_JSON_URL = "https://hiroyokoyama.github.io/moleditpy-plugins/REGISTRY/plugins.json"
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "plugin_installer.json")

# Global flag to ensure startup check runs only once per session
_startup_check_performed = False


def _refresh_plugin_menus(mw):
    """
    After discover_plugins(), strip stale plugin-managed menu/toolbar entries and
    rebuild them from the newly populated plugin_manager registries.

    This allows installs, updates, and removals to take effect immediately without
    requiring an application restart.
    """
    PLUGIN_ACTION_TAG = "plugin_managed"

    def _clean_menu(menu):
        """Recursively remove plugin-managed leaf actions and now-empty submenus."""
        for action in list(menu.actions()):
            submenu = action.menu()
            if submenu is not None:
                _clean_menu(submenu)
                # Remove the submenu entry itself if it has no visible items left
                has_content = any(not a.isSeparator() for a in submenu.actions())
                if not has_content:
                    menu.removeAction(action)
            elif action.data() == PLUGIN_ACTION_TAG:
                menu.removeAction(action)

    try:
        for top_action in list(mw.menuBar().actions()):
            top_menu = top_action.menu()
            if top_menu is not None:
                _clean_menu(top_menu)
    except Exception as e:
        logging.warning("Plugin Installer: menu cleanup error: %s", e)

    if not hasattr(mw, 'init_manager'):
        return

    im = mw.init_manager
    # Allow _add_registered_plugin_actions to re-add a separator before any new
    # top-level plugin menus (the flag is only used for first-time top-level menus,
    # not for the existing Plugin menu, so it is safe to reset).
    try:
        im._plugin_menubar_separator_added = False
    except Exception as e:
        logging.warning("Plugin Installer: could not reset _plugin_menubar_separator_added: %s", e)

    try:
        if hasattr(im, '_add_registered_plugin_actions'):
            im._add_registered_plugin_actions()
    except Exception as e:
        logging.warning("Plugin Installer: menu rebuild error: %s", e)

    try:
        if hasattr(im, '_add_plugin_toolbar_actions'):
            im._add_plugin_toolbar_actions()
    except Exception as e:
        logging.warning("Plugin Installer: toolbar rebuild error: %s", e)

def load_settings():
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, 'r') as f:
                return json.load(f)
        except Exception as e:
            logging.warning("Plugin Installer: failed to load settings: %s", e)
    # Default to empty so we can detect "first run" (missing config)
    return {}

def save_settings(settings):
    try:
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(settings, f)
    except Exception as e:
        logging.warning("Plugin Installer: failed to save settings: %s", e)

def initialize(context):
    """
    Called by PluginManager at startup (and on plugin reload).
    Registers the Plugin menu entry every time, and schedules the startup
    update check once per session.
    """
    global _startup_check_performed

    mw = context.get_main_window()
    if not mw:
        return


    # Schedule the startup update check only once per session.
    if not _startup_check_performed:
        _startup_check_performed = True
        settings = load_settings()
        if settings.get("check_at_startup", None):
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
        logging.warning("Plugin Installer: auto-run failed: %s", e)

def run(main_window):
    if hasattr(main_window, 'host'):
        main_window = main_window.host
    if hasattr(main_window, 'host'):
        main_window = main_window.host
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
                except Exception:
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
        if not getattr(self, 'missing_deps', None):
            return
            
        # Security: sanitize dependency names (allow alphanumeric, underscore, hyphen, period)
        safe_deps = [d for d in self.missing_deps if re.match(r'^[a-zA-Z0-9_\-\.]+$', d)]
        
        if len(safe_deps) != len(self.missing_deps):
            QMessageBox.warning(self, "Security Warning", 
                                "Some dependencies contain invalid characters and were skipped.")
            
        if not safe_deps:
            return

        deps_str = " ".join(safe_deps)
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
            
            # Reload plugins and immediately update Plugin menu + toolbar
            mw = self.parent_installer.main_window
            mw.plugin_manager.discover_plugins(mw)
            _refresh_plugin_menus(mw)
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

        self.check_updates_silent()

    def calculate_sha256(self, filepath):
        sha256_hash = hashlib.sha256()
        with open(filepath, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()

    def compare_versions(self, v1, v2):
        """
        Compare two version strings. 
        Returns 1 if v1 > v2, -1 if v1 < v2, 0 if v1 == v2.
        Handles nested numbers like 2.10.0 > 2.5.0
        """
        try:
            def parse(v):
                # Extract digits, ignore non-numeric suffixes for simple comparison
                return [int(x) for x in re.findall(r'\d+', v)]
            
            p1 = parse(v1)
            p2 = parse(v2)
            
            # Pad with zeros to handle different lengths (e.g. 2.5 vs 2.5.1)
            length = max(len(p1), len(p2))
            p1 += [0] * (length - len(p1))
            p2 += [0] * (length - len(p2))
            
            if p1 > p2: return 1
            if p1 < p2: return -1
            return 0
        except Exception as e:
            logging.warning("Plugin Installer: version comparison failed (%r vs %r): %s", v1, v2, e)
            # Fallback to string comparison
            if v1 > v2: return 1
            if v1 < v2: return -1
            return 0

    def init_ui(self):
        layout = QVBoxLayout(self)

    def _get_package_name(self):
        """Return the running PyPI package name: 'moleditpy' or 'moleditpy-linux'."""
        try:
            from moleditpy.utils.constants import VERSION  # noqa: F401
            return "moleditpy"
        except ImportError:
            pass
        try:
            from moleditpy_linux.modules.constants import VERSION  # noqa: F401
            return "moleditpy-linux"
        except ImportError:
            pass
        return "moleditpy"

    def get_app_version(self):
        """Robustly detect MoleditPy version."""
        try:
            from moleditpy.utils.constants import VERSION as APP_VERSION
            return APP_VERSION
        except ImportError:
            try:
                from moleditpy_linux.modules.constants import VERSION as APP_VERSION
                return APP_VERSION
            except ImportError:
                try:
                    main_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
                    if main_dir not in sys.path:
                        sys.path.append(main_dir)
                    try:
                        from moleditpy.utils.constants import VERSION as APP_VERSION
                    except ImportError:
                        from moleditpy_linux.modules.constants import VERSION as APP_VERSION
                    return APP_VERSION
                except:
                    return "0.0.0"

    def init_ui(self):
        self.latest_app_version = "Checking..."
        layout = QVBoxLayout(self)

        # Header Info
        app_ver = self.get_app_version()
        pkg_name = self._get_package_name()

        info_layout = QHBoxLayout()

        self.lbl_current_version = QLabel(f"<b>{pkg_name} Version:</b> {app_ver}")
        self.lbl_current_version.setStyleSheet("font-size: 14px;")
        info_layout.addWidget(self.lbl_current_version)
        
        info_layout.addSpacing(20)
        
        self.lbl_latest_version = QLabel(f"<b>Latest ({pkg_name} on PyPI):</b> {self.latest_app_version}")
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

        self.btn_update_all = QPushButton("Update All")
        self.btn_update_all.clicked.connect(self.update_all_plugins)
        self.btn_update_all.setEnabled(False)
        btn_layout.addWidget(self.btn_update_all)

        btn_layout.addStretch()

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.on_close_clicked) # Changed connection
        btn_layout.addWidget(close_btn)
        
        layout.addLayout(btn_layout)

    def save_settings_ui(self):
        self.settings["check_at_startup"] = self.chk_startup.isChecked()
        save_settings(self.settings)

    def fetch_pypi_version(self):
        package_name = self._get_package_name()
        try:
            url = f"https://pypi.org/pypi/{package_name}/json"
            with urllib.request.urlopen(url, timeout=5) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode('utf-8'))
                    return data['info']['version']
        except Exception as e:
            logging.warning("Plugin Installer: failed to fetch PyPI version for %s: %s", package_name, e)
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
        self.latest_app_version = latest
        app_ver = self.get_app_version()
            
        if latest != "Unknown" and self.compare_versions(latest, app_ver) > 0:
             self.updates_found = True

        self.remote_data = self.fetch_remote_data()
        self.populate_table(silent=True)
        
        if self.updates_found:
             self.setWindowTitle("Plugin Installer (Updates Available!)")

    def check_updates(self):
        # Update App Version
        pkg_name = self._get_package_name()
        self.lbl_latest_version.setText(f"<b>Latest ({pkg_name} on PyPI):</b> Checking...")
        QApplication.processEvents() 
        
        latest = self.fetch_pypi_version()
        self.latest_app_version = latest
        
        app_ver = self.get_app_version()

        color = "gray"
        if latest != "Unknown":
            comp = self.compare_versions(latest, app_ver)
            if comp > 0:
                 color = "#dc3545" # Red
                 status_text = f"{latest} (Update Available)"
            elif comp == 0:
                 color = "green"
                 status_text = f"{latest} (Up to date)"
            else:
                 color = "blue"
                 status_text = f"{latest} (Dev/Newer)"
        else:
             status_text = "Unknown"

        self.lbl_latest_version.setText(f"<b>Latest ({pkg_name} on PyPI):</b> <span style='color:{color}'>{status_text}</span>")

        self.remote_data = self.fetch_remote_data()
        self.populate_table()

    def populate_table(self, silent=False):
        self.table.setRowCount(0)
        self.updates_found = False
        plugin_updates_found = False

        # Check app version upgrade button visibility once per populate
        app_ver = self.get_app_version()
        if self.latest_app_version != "Checking..." and self.latest_app_version != "Unknown":
            try:
                if self.compare_versions(self.latest_app_version, app_ver) > 0:
                    self.updates_found = True
                    self.btn_upgrade_app.setVisible(True)
                    self.btn_upgrade_app.setText(f"Copy upgrade command (v{self.latest_app_version})")
                else:
                    self.btn_upgrade_app.setVisible(False)
            except Exception as e:
                logging.warning("Plugin Installer: upgrade button version check failed: %s", e)
                self.btn_upgrade_app.setVisible(False)
        
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
            remote_info = remote_map.get(name, None)
            local_info = installed_map.get(name, None)
            
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
                target_file = local_info.get('filepath', None)

            if remote_info:
                remote_ver = remote_info.get('version', 'Unknown')
                # Prioritize remote author info
                author = remote_info.get('author', 'Unknown')
            
            # If author still unknown, try local
            if author == "Unknown" and is_installed:
                author = local_info.get('author', 'Unknown')

            if is_installed:
                if remote_info:
                    
                    if local_ver != 'Unknown' and remote_ver != 'Unknown':
                        comp = self.compare_versions(remote_ver, local_ver)
                        if comp > 0:
                            status = "Update Available"
                            color = QColor("#f8d7da")
                            can_update = True
                            self.updates_found = True
                            plugin_updates_found = True
                        elif comp == 0:
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

        # Enable/disable Update All button — only for plugin updates, not app updates
        if getattr(self, 'btn_update_all', None) is not None:
            self.btn_update_all.setEnabled(plugin_updates_found)

        # Re-apply filter if needed
        if getattr(self, 'search_input', None) is not None:
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



    def update_all_plugins(self):
        """Show a single confirmation dialog then update all plugins with 'Update Available'."""
        rows_to_update = []
        for row in range(self.table.rowCount()):
            status_item = self.table.item(row, 4)
            if status_item and status_item.text() == "Update Available":
                btn = self.table.cellWidget(row, 5)
                if btn and isinstance(btn, QPushButton):
                    rows_to_update.append(btn)

        if not rows_to_update:
            QMessageBox.information(self, "Update All", "No plugins with available updates found.")
            return

        names = [btn.property("plugin_name") for btn in rows_to_update]
        name_list = "\n".join(f"  - {n}" for n in names)
        ret = QMessageBox.question(
            self, "Update All",
            f"Update the following {len(names)} plugin(s)?\n\n{name_list}",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes,
        )
        if ret != QMessageBox.StandardButton.Yes:
            return

        self._batch_updating = True
        succeeded = []
        failed = []
        try:
            for btn in rows_to_update:
                name = btn.property("plugin_name")
                try:
                    btn.click()
                    succeeded.append(name)
                except Exception as e:
                    failed.append(f"{name}: {e}")
        finally:
            self._batch_updating = False

        lines = [f"Updated {len(succeeded)} of {len(names)} plugin(s)."]
        if succeeded:
            lines.append("\nSucceeded:\n" + "\n".join(f"  - {n}" for n in succeeded))
        if failed:
            lines.append("\nFailed:\n" + "\n".join(f"  - {n}" for n in failed))
        QMessageBox.information(self, "Update All Complete", "\n".join(lines))

    def show_plugin_details(self, row, col):
        name_item = self.table.item(row, 0)
        if not name_item: return
        name = name_item.text()
        
        # Find info
        remote_info = None
        for entry in self.remote_data:
            if entry.get('name', None) == name:
                remote_info = entry
                break
        
        # Also check local
        local_info = None
        for p in self.main_window.plugin_manager.plugins:
            if p.get('name', None) == name:
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
            target_file = local_info.get('filepath', None)
            
            if remote_info:
                version = f"Installed: {local_ver}\nLatest: {version}"
            else:
                version = f"Installed: {local_ver}"

        dialog = PluginDetailsDialog(self, name, author, version, description, dependencies, local_info, target_file)
        dialog.exec()

    def on_close_clicked(self):
        self.accept()

    def closeEvent(self, event):
        # Re-discover plugins so any installs/removals that happened this session
        # are reflected, then immediately rebuild the Plugin menu and toolbar.
        if self.main_window and hasattr(self.main_window, 'plugin_manager'):
            try:
                self.main_window.plugin_manager.discover_plugins(self.main_window)
                _refresh_plugin_menus(self.main_window)
            except Exception as e:
                logging.warning("Plugin Installer: failed to refresh menus on close: %s", e)

        # Sync any open PluginManagerWindow
        for widget in QApplication.topLevelWidgets():
            if type(widget).__name__ == 'PluginManagerWindow':
                if hasattr(widget, 'refresh_plugin_list'):
                    try:
                        widget.refresh_plugin_list()
                    except Exception as e:
                        logging.warning("Plugin Installer: failed to sync PluginManagerWindow: %s", e)

        super().closeEvent(event)

    def copy_upgrade_command(self):
        package_name = self._get_package_name()
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
                except Exception:
                    missing_deps.append(dep)
            
            if missing_deps:
                # Show Details Dialog instead of installing
                # We need to reconstruct the arguments for show_plugin_details logic or just call it if we have row/col
                # But we don't know the row here easily without searching.
                # Easier to just instantiate the dialog directly since we have the data
                
                # Fetch full info to pass to dialog
                remote_info = None
                for entry in self.remote_data:
                    if entry.get('name', None) == plugin_name:
                        remote_info = entry
                        break
                
                # Local info
                local_info = None
                for p in self.main_window.plugin_manager.plugins:
                    if p.get('name', None) == plugin_name:
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
        # Also skip if we're in batch-update mode (single dialog already shown in update_all_plugins)
        action_verb = "Update" if is_installed else "Install"
        if not user_confirmed_intent and not getattr(self, '_batch_updating', False):
            ret = QMessageBox.question(self, f"{action_verb} Plugin",
                                       f"{action_verb} '{plugin_name}'?\nThis will download and install the plugin.",
                                       QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if ret != QMessageBox.StandardButton.Yes:
                return

        # Resolve URL (handle relative)
        final_url = download_url
        if not download_url.startswith("http"):
            # Construct absolute URL relative to plugins.json location
            base_url = REMOTE_JSON_URL.rsplit('/', 1)[0]
            from urllib.parse import urljoin
            final_url = urljoin(base_url + "/", download_url)

        logging.info("Plugin Installer: %sing %s from: %s", action_verb, plugin_name, final_url)
        
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
                        
                        # SHA256 Verification
                        # Find remote_info for this plugin to get SHA256
                        remote_info = None
                        for entry in self.remote_data:
                            if entry.get('name', None) == plugin_name:
                                remote_info = entry
                                break
                        
                        remote_sha256 = None
                        if remote_info and 'sha256' in remote_info:
                            remote_sha256 = remote_info['sha256']
                        
                        # Warn if SHA256 not available
                        if not remote_sha256:
                            ret = QMessageBox.warning(self, "Security Warning", 
                                                      f"SHA256 checksum is not available for '{plugin_name}'.\n\nThis means the plugin file's integrity cannot be verified.\n\nDo you want to proceed anyway?",
                                                      QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                                                      QMessageBox.StandardButton.No)
                            if ret != QMessageBox.StandardButton.Yes:
                                return
                        else:
                            # Verify SHA256 checksum
                            calculated_sha256 = self.calculate_sha256(download_path)
                            if calculated_sha256 != remote_sha256:
                                # Show details and ask for support
                                detail_msg = (
                                    f"SHA256 Checksum Mismatch for '{plugin_name}'!\n\n"
                                    f"Expected: {remote_sha256}\n"
                                    f"Calculated: {calculated_sha256}\n\n"
                                    f"This could indicate:\n"
                                    f"• The file was corrupted during download\n"
                                    f"• The file has been tampered with\n"
                                    f"• The plugin registry data is outdated\n\n"
                                    f"Installation has been stopped for your safety.\n\n"
                                    f"Please report this issue to the plugin author or visit:\n"
                                    f"https://github.com/HiroYokoyama/moleditpy-plugins/issues"
                                )
                                QMessageBox.critical(self, "Security Error", detail_msg)
                                return

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
                                        logging.info("Plugin Installer: overwriting existing plugin file: %s", target_file)
                                        shutil.copy2(download_path, target_file)
                                        did_manual_overwrite = True
                                    except Exception as e:
                                        logging.warning("Plugin Installer: manual overwrite failed: %s — falling back to manager", e)
                                
                                # Check: Folder Plugin Overwrite (Target is __init__.py, Source is ZIP)
                                elif (file_name_existing == "__init__.py" and 
                                      final_url.lower().endswith(".zip")):
                                      
                                      try:
                                          target_dir = os.path.dirname(target_file)
                                          logging.info("Plugin Installer: overwriting existing folder plugin: %s", target_dir)
                                          
                                          # Extract and validate (Security Fix)
                                          extract_temp = os.path.join(temp_dir, "extracted")
                                          os.makedirs(extract_temp, exist_ok=True)
                                          with zipfile.ZipFile(download_path, 'r') as z:
                                              for member in z.infolist():
                                                  # Security check for directory traversal
                                                  target_path = os.path.join(extract_temp, member.filename)
                                                  if not os.path.abspath(target_path).startswith(os.path.abspath(extract_temp)):
                                                      logging.warning("Plugin Installer: skipping suspicious zip entry (path traversal): %s", member.filename)
                                                      continue
                                                  z.extract(member, extract_temp)
                                          
                                          # Determine source content
                                          # If zip contains a single folder, use that inner folder as source
                                          items = os.listdir(extract_temp)
                                          if len(items) == 1 and os.path.isdir(os.path.join(extract_temp, items[0])):
                                               source = os.path.join(extract_temp, items[0])
                                          else:
                                               source = extract_temp
                                          
                                          # Backup settings.json
                                          settings_backup = None
                                          settings_path = os.path.join(target_dir, "settings.json")
                                          if os.path.exists(settings_path):
                                               try:
                                                   with open(settings_path, 'r', encoding='utf-8') as f:
                                                       settings_backup = f.read()
                                               except Exception as e:
                                                   logging.warning("Plugin Installer: failed to backup settings.json (folder): %s", e)

                                          # Overwrite using copytree with dirs_exist_ok=True
                                          # This overwrites existing files and adds new ones
                                          shutil.copytree(source, target_dir, dirs_exist_ok=True)
                                          
                                          # Restore settings.json
                                          if settings_backup:
                                               try:
                                                   with open(settings_path, 'w', encoding='utf-8') as f:
                                                       f.write(settings_backup)
                                                   logging.info("Plugin Installer: restored settings.json (folder)")
                                               except Exception as e:
                                                   logging.warning("Plugin Installer: failed to restore settings.json (folder): %s", e)
                                          
                                          did_manual_overwrite = True

                                      except Exception as e:
                                          logging.warning("Plugin Installer: manual folder overwrite failed: %s — falling back to manager", e)
                            
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
                                                logging.warning("Plugin Installer: failed to backup settings.json: %s", e)

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
                                            logging.info("Plugin Installer: restored settings.json")
                                        except Exception as e:
                                            logging.warning("Plugin Installer: failed to restore settings.json: %s", e)
                                
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
                                        except Exception as e:
                                            logging.warning("Plugin Installer: failed to remove old .py during transition: %s", e)

                                    # Case 2: Transition from Package (Folder) to Single File (.py)
                                    elif not new_is_zip and old_is_init:
                                        try:
                                            parent_dir = os.path.dirname(target_file)
                                            shutil.rmtree(parent_dir)
                                        except Exception as e:
                                            logging.warning("Plugin Installer: failed to remove old folder during transition: %s", e)
                            
                            # Success Message
                            # Success Message & Dependency Check
                            verb_past = "updated" if is_installed else "installed"
                            
                            if missing_deps:
                                # Reduce dialogs: Combine Success + Warning
                                QMessageBox.warning(self, "Installation Complete (Dependencies Missing)", 
                                                    f"Successfully {verb_past} '{plugin_name}'.\n\nHowever, you must install the missing dependencies for it to work.\n\nOpening details window...")
                                dialog = PluginDetailsDialog(self, plugin_name, author, version, description, dependencies, local_info, target_file)
                                dialog.exec()
                            elif not getattr(self, '_batch_updating', False):
                                QMessageBox.information(self, "Success", f"Successfully {verb_past} '{plugin_name}'.")
                            
                            # Reload plugins and immediately update Plugin menu + toolbar
                            self.main_window.plugin_manager.discover_plugins(self.main_window)
                            _refresh_plugin_menus(self.main_window)
                            # Refresh our list (reuse already-fetched remote_data — no extra network call)
                            self.populate_table()
                            
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



