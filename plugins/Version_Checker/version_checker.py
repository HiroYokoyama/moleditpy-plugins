#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Version Checker Plugin
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
PLUGIN_NAME = "Version Checker"
PLUGIN_VERSION = "2025.12.26"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Checks installed plugins against the official repository for updates."

REMOTE_JSON_URL = "https://raw.githubusercontent.com/HiroYokoyama/moleditpy-plugins/refs/heads/main/explorer/plugins.json"
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "version_checker.json")

def load_settings():
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, 'r') as f:
                return json.load(f)
        except:
            pass
    return {"check_at_startup": True}

def save_settings(settings):
    try:
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(settings, f)
    except Exception as e:
        print(f"Error saving version checker settings: {e}")

def initialize(context):
    """
    Called by PluginManager at startup.
    Checks for updates if enabled in settings.
    """
    # Register menu action
    settings = load_settings()
    if settings.get("check_at_startup", True):
        mw = context.get_main_window()
        if mw:
            # Delay check to ensure main window is visible
            from PyQt6.QtCore import QTimer
            # Delay 2 seconds
            QTimer.singleShot(2000, lambda: perform_startup_check(mw))

def perform_startup_check(mw):
    try:
        # Run silent check
        checker = VersionCheckerWindow(mw, auto_check=True)
        if checker.updates_found:
            # Ask user if they want to see updates or skip
            reply = QMessageBox.question(
                mw, "Version Checker", 
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
        print(f"Version Check Auto-run failed: {e}")

def run(main_window):
    """
    Entry point for the plugin.
    Display the Version Checker window.
    """
    dialog = VersionCheckerWindow(main_window)
    dialog.exec()

class VersionCheckerWindow(QDialog):
    def __init__(self, main_window, auto_check=False):
        super().__init__(main_window)
        self.setWindowTitle("Plugin Version Checker")
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
        self.chk_startup.setChecked(self.settings.get("check_at_startup", True))
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
             self.setWindowTitle("Plugin Version Checker (Updates Available!)")

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
        
        # Create map
        remote_map = {}
        for entry in self.remote_data:
            name = entry.get('name', '')
            remote_map[name] = entry
            if 'id' in entry:
                 remote_map[entry['id']] = entry

        installed_plugins = self.main_window.plugin_manager.plugins
        
        for plugin in installed_plugins:
            name = plugin.get('name', 'Unknown')
            local_ver = plugin.get('version', 'Unknown')
            
            remote_info = remote_map.get(name)
            remote_ver = "Unknown"
            status = "Unknown"
            color = None
            update_available = False
            
            if remote_info:
                remote_ver = remote_info.get('version', 'Unknown')
                if local_ver != 'Unknown' and remote_ver != 'Unknown':
                    if local_ver == remote_ver:
                        status = "Up to date"
                        color = QColor("#d4edda")
                    else:
                        if remote_ver > local_ver:
                            status = "Update Available"
                            color = QColor("#f8d7da")
                            update_available = True
                            self.updates_found = True
                        else:
                            status = "Up to date (Newer?)"
                            color = QColor("#cce5ff")
                else:
                    status = "Version Mismatch"
            else:
                status = "Not in Registry"
                color = QColor("#e2e3e5")

            if silent and not update_available:
                # If doing silent check, we might not populate everything? 
                # Actually we should populate so if user shows window it's ready.
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
            if update_available and remote_info and 'downloadUrl' in remote_info:
                btn_update = QPushButton("Update")
                # Store data in button property to retrieve later
                btn_update.setProperty("plugin_name", name)
                btn_update.setProperty("download_url", remote_info['downloadUrl'])
                btn_update.setProperty("target_file", plugin.get('filepath'))
                btn_update.clicked.connect(self.on_update_clicked)
                self.table.setCellWidget(row, 4, btn_update)
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
