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
    QHeaderView, QMessageBox, QAbstractItemView, QApplication, QCheckBox
)
from PyQt6.QtCore import Qt, QUrl
from PyQt6.QtGui import QDesktopServices, QColor, QIcon

# --- Metadata ---
PLUGIN_NAME = "Plugin Installer"
PLUGIN_VERSION = "2026.01.04"
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
        
        info_layout.addStretch()
        layout.addLayout(info_layout)
        
        # Add a small line
        line = QLabel()
        line.setFrameStyle(QLabel.Shape.HLine | QLabel.Shadow.Sunken)
        layout.addWidget(line)

        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(5) # Added "Action" column
        self.table.setHorizontalHeaderLabels(["Plugin Name", "Local Version", "Latest Version", "Status", "Action"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents) # Action
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
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
        close_btn.clicked.connect(self.accept)
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
            
            local_ver = "Unknown"
            remote_ver = "Unknown"
            status = "Unknown"
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
            
            if is_installed:
                if remote_info:
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
                            status = "Up to date (Newer?)"
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
            
            self.table.setItem(row, 0, QTableWidgetItem(name))
            self.table.setItem(row, 1, QTableWidgetItem(str(local_ver)))
            self.table.setItem(row, 2, QTableWidgetItem(str(remote_ver)))
            
            status_item = QTableWidgetItem(status)
            if color:
                status_item.setBackground(color)
            self.table.setItem(row, 3, status_item)
            
            # Action Button
            if remote_info and 'downloadUrl' in remote_info:
                label = ""
                if can_update:
                    if status == "Update Available":
                        label = "Update"
                    else:
                        label = "Reinstall"
                elif can_download:
                     label = "Download"
                
                if label:
                    btn_action = QPushButton(label)
                    btn_action.setProperty("plugin_name", name)
                    btn_action.setProperty("download_url", remote_info['downloadUrl'])
                    btn_action.setProperty("target_file", target_file)
                    # We might need to differentiate clean install vs overwrite in message?
                    # The on_update_clicked handles prompt.
                    btn_action.clicked.connect(self.on_update_clicked)
                    self.table.setCellWidget(row, 4, btn_action)
            else:
                self.table.setItem(row, 4, QTableWidgetItem(""))

    def on_update_clicked(self):
        btn = self.sender()
        if not btn: return
        
        plugin_name = btn.property("plugin_name")
        download_url = btn.property("download_url")
        target_file = btn.property("target_file")
        
        if not download_url or not target_file:
            QMessageBox.warning(self, "Error", "Cannot update: Missing URL or file path.")
            return

        # Confirm
        ret = QMessageBox.question(self, "Update Plugin", 
                                   f"Update '{plugin_name}'?\nThis will overwrite the local file.",
                                   QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if ret != QMessageBox.StandardButton.Yes:
            return

        # Resolve URL (handle relative)
        final_url = download_url
        if not download_url.startswith("http"):
            # Construct absolute URL relative to plugins.json location
            base_url = REMOTE_JSON_URL.rsplit('/', 1)[0]
            # Handle ../ parent traversal if needed, but urllib.parse.urljoin is best
            from urllib.parse import urljoin
            final_url = urljoin(base_url + "/", download_url)

        print(f"Downloading update from: {final_url}")
        
        # Download and Overwrite
        try:
            # Ensure target directory exists
            target_dir = os.path.dirname(target_file)
            if not os.path.exists(target_dir):
                os.makedirs(target_dir, exist_ok=True)
                
            with urllib.request.urlopen(final_url, timeout=10) as response:
                if response.status == 200:
                    content = response.read()
                    with open(target_file, 'wb') as f:
                        f.write(content)
                    
                    QMessageBox.information(self, "Success", f"Updated '{plugin_name}'.\nReloading plugins...")
                    
                    # Reload plugins
                    self.main_window.plugin_manager.discover_plugins(self.main_window)
                    # Refresh our list
                    self.check_updates()
                    
                else:
                    QMessageBox.warning(self, "Error", f"Download failed with status: {response.status}")
        except Exception as e:
             QMessageBox.warning(self, "Update Error", f"Failed to update plugin:\n{e}")
