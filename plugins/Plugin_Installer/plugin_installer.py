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
import urllib.parse
import urllib.request
import urllib.error
import ast
import logging
from PyQt6.QtCore import QThread, pyqtSignal, Qt, QTimer
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QHeaderView,
    QMessageBox,
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QLineEdit,
    QProgressBar,
    QProgressDialog,
)
from PyQt6.QtGui import QColor
import importlib.metadata
import hashlib
import re

import shutil
import tempfile

# --- Metadata ---
PLUGIN_NAME = "Plugin Installer"
PLUGIN_VERSION = "2026.06.26"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Checks for updates, installs new plugins, and allows manual reinstallation."
)

REMOTE_JSON_URL = (
    "https://hiroyokoyama.github.io/moleditpy-plugins/REGISTRY/plugins.json"
)
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "plugin_installer.json")

_startup_check_performed = False
_context = None


def _read_plugin_version_ast(filepath: str) -> str:
    """Read PLUGIN_VERSION from a .py file using AST — no import, no side effects."""
    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            source = f.read()
        tree = ast.parse(source, filename=filepath)
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id == "PLUGIN_VERSION":
                        if isinstance(node.value, ast.Constant):
                            return str(node.value.value)
    except Exception as e:
        logging.warning("Plugin Installer: AST version read failed for %s: %s", filepath, e)
    return "Unknown"


def _refresh_plugin_menus(mw):
    global _context
    if _context is not None and hasattr(_context, "rebuild_menus"):
        _context.rebuild_menus()
        return
    if hasattr(mw, "plugin_menu_manager"):
        try:
            mw.plugin_menu_manager.rebuild_plugin_menus()
        except Exception as e:
            logging.warning("Plugin Installer: menu rebuild error: %s", e)


def load_settings():
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            logging.warning("Plugin Installer: failed to load settings: %s", e)
    return {}


def save_settings(settings):
    try:
        with open(SETTINGS_FILE, "w", encoding="utf-8") as f:
            json.dump(settings, f)
    except Exception as e:
        logging.warning("Plugin Installer: failed to save settings: %s", e)


def initialize(context):
    global _startup_check_performed, _context
    _context = context

    def _open_installer():
        dlg = PluginInstallerWindow(context.get_main_window())
        dlg.exec()

    context.add_menu_action("Plugin/Plugin Installer...", _open_installer)

    mw = context.get_main_window()
    if not mw:
        return

    if not _startup_check_performed:
        _startup_check_performed = True
        settings = load_settings()
        if "check_at_startup" not in settings:
            QTimer.singleShot(3000, lambda: ask_user_permission(mw))
        elif settings.get("check_at_startup"):
            QTimer.singleShot(2000, lambda: perform_startup_check(mw))


def ask_user_permission(mw):
    reply = QMessageBox.question(
        mw,
        "Plugin Installer",
        "Do you want to enable automatic plugin update checks at startup?\n\n(You can change this later in the Plugin Installer)",
        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        QMessageBox.StandardButton.No,
    )
    settings = load_settings()
    if reply == QMessageBox.StandardButton.Yes:
        settings["check_at_startup"] = True
        save_settings(settings)
        perform_startup_check(mw)
    else:
        settings["check_at_startup"] = False
        save_settings(settings)


def perform_startup_check(mw):
    try:
        checker = PluginInstallerWindow(mw, auto_check=True)
        if checker.updates_found:
            reply = QMessageBox.question(
                mw,
                "Plugin Installer",
                "Updates are available for MoleditPy or installed plugins.\nDo you want to check them now?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.Yes,
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


def is_app_version_compatible(app_version: str, specifier: str) -> bool:
    if not specifier or specifier.strip() in ("*", ""):
        return True

    def parse_ver(v_str):
        parts = [int(x) for x in re.findall(r"\d+", v_str)]
        while len(parts) < 3:
            parts.append(0)
        return tuple(parts[:3])

    try:
        app_v = parse_ver(app_version)
    except Exception:
        return True

    specs = [s.strip() for s in specifier.split(",")]
    for spec in specs:
        if not spec:
            continue
        if "*" in spec:
            clean_spec = spec.replace("==", "").strip()
            prefix = clean_spec.split("*")[0].rstrip(".")
            try:
                prefix_parts = [int(x) for x in re.findall(r"\d+", prefix)]
                if app_v[: len(prefix_parts)] != tuple(prefix_parts):
                    return False
            except Exception:
                pass
            continue

        match = re.match(r"^(?P<op>>=|<=|~=|>=|>|<|==|!=)?\s*(?P<ver>[\d\.]+)", spec)
        if not match:
            continue
        op = match.group("op") or "=="
        spec_v_str = match.group("ver")
        try:
            spec_v = parse_ver(spec_v_str)
            if op == "==":
                spec_parts = [int(x) for x in re.findall(r"\d+", spec_v_str)]
                if app_v[: len(spec_parts)] != tuple(spec_parts):
                    return False
            elif op == "!=":
                spec_parts = [int(x) for x in re.findall(r"\d+", spec_v_str)]
                if app_v[: len(spec_parts)] == tuple(spec_parts):
                    return False
            elif op == ">=":
                if app_v < spec_v:
                    return False
            elif op == "<=":
                if app_v > spec_v:
                    return False
            elif op == ">":
                if app_v <= spec_v:
                    return False
            elif op == "<":
                if app_v >= spec_v:
                    return False
            elif op == "~=":
                spec_parts = [int(x) for x in re.findall(r"\d+", spec_v_str)]
                if not spec_parts:
                    return False
                if app_v < spec_v:
                    return False
                n_prefix = len(spec_parts) - 1
                if n_prefix > 0:
                    if app_v[:n_prefix] != tuple(spec_parts[:n_prefix]):
                        return False
                else:
                    if app_v[0] != spec_parts[0]:
                        return False
        except Exception:
            pass
    return True


def parse_dependency(dep_str: str) -> tuple[str, str]:
    if ":" in dep_str:
        return "", ""
    dep_str = dep_str.strip()
    if not dep_str:
        return "", ""
    match = re.match(r"^([a-zA-Z0-9_\-\.]+)(.*)$", dep_str)
    if not match:
        return dep_str, ""
    name = match.group(1)
    rest = match.group(2).strip()
    if rest.startswith("["):
        bracket_count = 0
        end_idx = -1
        for idx, char in enumerate(rest):
            if char == "[":
                bracket_count += 1
            elif char == "]":
                bracket_count -= 1
                if bracket_count == 0:
                    end_idx = idx
                    break
        if end_idx != -1:
            rest = rest[end_idx + 1 :].strip()
    if ";" in rest:
        rest = rest.split(";", 1)[0].strip()
    return name, rest


def check_dependency_satisfied(dep_str: str) -> bool:
    name, specifier = parse_dependency(dep_str)
    if not name:
        return False
    try:
        dist = importlib.metadata.distribution(name)
        installed_ver = dist.version
    except importlib.metadata.PackageNotFoundError:
        return False
    if not specifier:
        return True
    return is_app_version_compatible(installed_ver, specifier)


def sanitize_and_quote_dependency(dep_str: str) -> str:
    if ";" in dep_str:
        return ""
    name, specifier = parse_dependency(dep_str)
    if not name:
        return ""
    if not re.match(r"^[a-zA-Z0-9_\-\.]+$", name):
        return ""
    if specifier:
        if not re.match(r"^[a-zA-Z0-9_\-\.\*>=<!~,\s]+$", specifier):
            return ""
        return f'"{name}{specifier}"'
    return name


# ---------------------------------------------------------------------------
# Background worker: fetches PyPI version + remote registry without blocking
# ---------------------------------------------------------------------------


class _FetchWorker(QThread):
    """Fetches the PyPI app version and remote registry JSON in a background thread."""

    done = pyqtSignal(str, list)  # (pypi_version, remote_entries)

    def __init__(self, package_name: str, remote_url: str) -> None:
        super().__init__()
        self._pkg = package_name
        self._url = remote_url

    def run(self) -> None:
        pypi_ver = self._fetch_pypi()
        remote = self._fetch_remote()
        self.done.emit(pypi_ver, remote)

    def _fetch_pypi(self) -> str:
        try:
            url = f"https://pypi.org/pypi/{self._pkg}/json"
            with urllib.request.urlopen(url, timeout=5) as r:
                if r.status == 200:
                    return json.loads(r.read().decode("utf-8"))["info"]["version"]
        except Exception as e:
            logging.warning("Plugin Installer: PyPI fetch failed: %s", e)
        return "Unknown"

    def _fetch_remote(self) -> list:
        try:
            with urllib.request.urlopen(self._url, timeout=5) as r:
                if r.status == 200:
                    return json.loads(r.read().decode("utf-8-sig"))
        except Exception as e:
            logging.warning("Plugin Installer: registry fetch failed: %s", e)
        return []


# ---------------------------------------------------------------------------
# Details dialog
# ---------------------------------------------------------------------------


class PluginDetailsDialog(QDialog):
    def __init__(
        self,
        parent,
        name,
        author,
        version,
        description,
        dependencies,
        local_info,
        target_file,
        supported_version="Unknown",
    ):
        super().__init__(parent)
        self.setWindowTitle(f"{name} - Details")
        self.resize(400, 370)
        self.parent_installer = parent
        self.plugin_name = name
        self.target_file = target_file

        layout = QVBoxLayout(self)

        lbl_name = QLabel(f"<h2>{name}</h2>")
        lbl_author = QLabel(f"<b>Author:</b> {author}")
        lbl_ver = QLabel(f"<b>Version:</b> {version}")
        lbl_supported = QLabel(f"<b>Supported MoleditPy:</b> {supported_version}")
        lbl_desc = QLabel(f"<b>Description:</b><br>{description}")
        lbl_desc.setWordWrap(True)

        layout.addWidget(lbl_name)
        layout.addWidget(lbl_author)
        layout.addWidget(lbl_ver)
        layout.addWidget(lbl_supported)
        layout.addSpacing(10)
        layout.addWidget(lbl_desc)

        if dependencies:
            layout.addSpacing(10)
            lbl_dep_header = QLabel("<b>Dependencies:</b>")
            layout.addWidget(lbl_dep_header)

            installed_deps = []
            missing_deps = []
            for dep in dependencies:
                if check_dependency_satisfied(dep):
                    installed_deps.append(dep)
                else:
                    missing_deps.append(dep)

            dep_text_parts = []
            if installed_deps:
                text = ", ".join(installed_deps)
                dep_text_parts.append(
                    f"<span style='color:green'>{text} (Installed)</span>"
                )
            if missing_deps:
                text = ", ".join(missing_deps)
                dep_text_parts.append(
                    f"<span style='color:red'>{text} (Missing)</span>"
                )

            lbl_deps = QLabel("<br>".join(dep_text_parts))
            lbl_deps.setWordWrap(True)
            layout.addWidget(lbl_deps)

            if missing_deps:
                self.missing_deps = missing_deps
                btn_install_all = QPushButton(
                    f"Copy install command ({len(missing_deps)})"
                )
                btn_install_all.setStyleSheet("color: blue;")
                btn_install_all.clicked.connect(self.copy_install_all_command)
                layout.addWidget(btn_install_all)

        layout.addStretch()

        btn_layout = QHBoxLayout()
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
        if not getattr(self, "missing_deps", None):
            return
        safe_deps = []
        for d in self.missing_deps:
            sanitized = sanitize_and_quote_dependency(d)
            if sanitized:
                safe_deps.append(sanitized)
        if len(safe_deps) != len(self.missing_deps):
            QMessageBox.warning(
                self,
                "Security Warning",
                "Some dependencies contain invalid characters and were skipped.",
            )
        if not safe_deps:
            return
        deps_str = " ".join(safe_deps)
        cmd = f"pip install {deps_str}"
        QApplication.clipboard().setText(cmd)
        QMessageBox.information(
            self,
            "Install Command",
            f"Command copied to clipboard:\n\n{cmd}\n\nPlease close MoleditPy and run this terminal command.",
        )

    def on_delete(self):
        ret = QMessageBox.question(
            self,
            "Delete Plugin",
            f"Are you sure you want to delete '{self.plugin_name}'?\nThis cannot be undone.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if ret != QMessageBox.StandardButton.Yes:
            return
        try:
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
            self.accept()

            mw = self.parent_installer.main_window
            mw.plugin_manager.discover_plugins(mw)
            _refresh_plugin_menus(mw)
            self.parent_installer.check_updates()
        except Exception as e:
            QMessageBox.warning(self, "Delete Error", f"Failed to delete plugin:\n{e}")


# ---------------------------------------------------------------------------
# Main installer window
# ---------------------------------------------------------------------------


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
        self._fetch_worker = None
        self._batch_updating = False
        self._batch_progress = None
        # Plugins installed this session but not yet reloaded into the host.
        # Keyed by plugin name; holds a minimal dict for table display so the
        # table shows correct status immediately without calling discover_plugins.
        # discover_plugins is called exactly once, in closeEvent.
        self._pending_installs: dict = {}
        self._needs_plugin_reload = False

        self.init_ui()
        self.finished.connect(self._on_finished)

        if auto_check:
            self.check_updates_silent()
        else:
            self.check_updates()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def init_ui(self):
        self.latest_app_version = "Checking..."
        layout = QVBoxLayout(self)

        app_ver = self.get_app_version()
        pkg_name = self._get_package_name()

        info_layout = QHBoxLayout()

        self.lbl_current_version = QLabel(f"<b>{pkg_name} Version:</b> {app_ver}")
        self.lbl_current_version.setStyleSheet("font-size: 14px;")
        info_layout.addWidget(self.lbl_current_version)

        info_layout.addSpacing(20)

        self.lbl_latest_version = QLabel(
            f"<b>Latest ({pkg_name} on PyPI):</b> {self.latest_app_version}"
        )
        self.lbl_latest_version.setStyleSheet("font-size: 14px; color: gray;")
        info_layout.addWidget(self.lbl_latest_version)

        self.btn_upgrade_app = QPushButton("Copy upgrade command")
        self.btn_upgrade_app.setVisible(False)
        self.btn_upgrade_app.setStyleSheet("color: blue; font-weight: bold;")
        self.btn_upgrade_app.clicked.connect(self.copy_upgrade_command)
        info_layout.addWidget(self.btn_upgrade_app)

        info_layout.addStretch()
        layout.addLayout(info_layout)

        # Thin indeterminate progress bar shown while fetching
        self._progress_bar = QProgressBar()
        self._progress_bar.setRange(0, 0)  # indeterminate
        self._progress_bar.setMaximumHeight(5)
        self._progress_bar.setTextVisible(False)
        self._progress_bar.setVisible(False)
        layout.addWidget(self._progress_bar)

        line = QLabel()
        line.setFrameStyle(QLabel.Shape.HLine | QLabel.Shadow.Sunken)
        layout.addWidget(line)

        search_layout = QHBoxLayout()
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search plugins by name or author...")
        self.search_input.textChanged.connect(self.filter_plugins)
        search_layout.addWidget(QLabel("Search:"))
        search_layout.addWidget(self.search_input)
        layout.addLayout(search_layout)

        self.table = QTableWidget()
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels(
            ["Plugin Name", "Author", "Local Ver", "Latest Ver", "Status", "Action"]
        )
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.cellDoubleClicked.connect(self.show_plugin_details)
        layout.addWidget(self.table)

        self.chk_startup = QCheckBox("Check for updates at startup")
        self.chk_startup.setChecked(self.settings.get("check_at_startup", False))
        self.chk_startup.toggled.connect(self.save_settings_ui)
        layout.addWidget(self.chk_startup)

        self._status_label = QLabel("")
        self._status_label.setStyleSheet("color: #856404; font-weight: bold;")
        layout.addWidget(self._status_label)

        btn_layout = QHBoxLayout()

        self._btn_refresh = QPushButton("Refresh")
        self._btn_refresh.clicked.connect(self.check_updates)
        btn_layout.addWidget(self._btn_refresh)

        self.btn_update_all = QPushButton("Update All")
        self.btn_update_all.clicked.connect(self.update_all_plugins)
        self.btn_update_all.setEnabled(False)
        btn_layout.addWidget(self.btn_update_all)

        btn_layout.addStretch()

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.on_close_clicked)
        btn_layout.addWidget(close_btn)

        layout.addLayout(btn_layout)

    # ------------------------------------------------------------------
    # Fetch: non-blocking path (visible dialog)
    # ------------------------------------------------------------------

    def _set_ui_busy(self, busy: bool) -> None:
        self._progress_bar.setVisible(busy)
        self._btn_refresh.setEnabled(not busy)
        if busy:
            self.btn_update_all.setEnabled(False)

    def check_updates(self):
        pkg_name = self._get_package_name()
        self.lbl_latest_version.setText(
            f"<b>Latest ({pkg_name} on PyPI):</b> Fetching..."
        )
        self._set_ui_busy(True)
        QApplication.processEvents()  # paint the bar before the thread starts

        # Cancel any in-flight fetch
        if self._fetch_worker and self._fetch_worker.isRunning():
            try:
                self._fetch_worker.done.disconnect()
            except Exception:
                pass
            self._fetch_worker.quit()
            self._fetch_worker.wait(300)

        self._fetch_worker = _FetchWorker(pkg_name, REMOTE_JSON_URL)
        self._fetch_worker.done.connect(self._on_fetch_done)
        self._fetch_worker.start()

    def _on_fetch_done(self, pypi_ver: str, remote_data: list) -> None:
        self._set_ui_busy(False)
        self.latest_app_version = pypi_ver
        self.remote_data = remote_data

        pkg_name = self._get_package_name()
        app_ver = self.get_app_version()

        color = "gray"
        if pypi_ver != "Unknown":
            comp = self.compare_versions(pypi_ver, app_ver)
            if comp > 0:
                color = "#dc3545"
                status_text = f"{pypi_ver} (Update Available)"
            elif comp == 0:
                color = "green"
                status_text = f"{pypi_ver} (Up to date)"
            else:
                color = "blue"
                status_text = f"{pypi_ver} (Dev/Newer)"
        else:
            status_text = "Unknown"

        self.lbl_latest_version.setText(
            f"<b>Latest ({pkg_name} on PyPI):</b> "
            f"<span style='color:{color}'>{status_text}</span>"
        )
        self.populate_table()

    # ------------------------------------------------------------------
    # Fetch: synchronous path (auto_check mode only — dialog not yet shown)
    # ------------------------------------------------------------------

    def check_updates_silent(self):
        latest = self.fetch_pypi_version()
        self.latest_app_version = latest
        app_ver = self.get_app_version()

        if latest != "Unknown" and self.compare_versions(latest, app_ver) > 0:
            self.updates_found = True

        self.remote_data = self.fetch_remote_data()
        self.populate_table(silent=True)

        if self.updates_found:
            self.setWindowTitle("Plugin Installer (Updates Available!)")

    def fetch_pypi_version(self):
        package_name = self._get_package_name()
        try:
            url = f"https://pypi.org/pypi/{package_name}/json"
            with urllib.request.urlopen(url, timeout=5) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode("utf-8"))
                    return data["info"]["version"]
        except Exception as e:
            logging.warning(
                "Plugin Installer: failed to fetch PyPI version for %s: %s",
                package_name,
                e,
            )
        return "Unknown"

    def fetch_remote_data(self):
        try:
            with urllib.request.urlopen(REMOTE_JSON_URL, timeout=5) as response:
                if response.status == 200:
                    return json.loads(response.read().decode("utf-8-sig"))
        except Exception as e:
            if not self.auto_check_mode:
                QMessageBox.warning(
                    self, "Network Error", f"Failed to fetch plugin data:\n{e}"
                )
        return []

    # ------------------------------------------------------------------
    # Download helper — chunked reads, callback per chunk
    # ------------------------------------------------------------------

    def _download_chunked(self, url: str, dest_path: str, on_progress=None) -> bool:
        """
        Download *url* to *dest_path* in 8 KB chunks.

        *on_progress(received: int, total: int) -> bool*
        Called after each chunk; return False to abort the download.
        If omitted the download still runs without cancellation support.

        Returns True on success, False if aborted or on network error.
        The caller is responsible for showing any progress UI — this method
        only drives the data transfer and reports progress via the callback.
        """
        try:
            with urllib.request.urlopen(url, timeout=15) as response:
                total = int(response.headers.get("Content-Length") or 0)
                received = 0
                with open(dest_path, "wb") as f:
                    while True:
                        chunk = response.read(8192)
                        if not chunk:
                            break
                        f.write(chunk)
                        received += len(chunk)
                        if on_progress is not None:
                            if not on_progress(received, total):
                                return False
            return True
        except Exception as e:
            logging.exception("Plugin Installer: download failed from %s: %s", url, e)
            QMessageBox.warning(self, "Download Error", f"Download failed:\n{e}")
            return False

    # ------------------------------------------------------------------
    # Table population
    # ------------------------------------------------------------------

    def populate_table(self, silent=False):
        self.table.setRowCount(0)
        self.updates_found = False
        plugin_update_count = 0

        app_ver = self.get_app_version()
        if (
            self.latest_app_version != "Checking..."
            and self.latest_app_version != "Unknown"
        ):
            try:
                if self.compare_versions(self.latest_app_version, app_ver) > 0:
                    self.updates_found = True
                    self.btn_upgrade_app.setVisible(True)
                    self.btn_upgrade_app.setText(
                        f"Copy upgrade command (v{self.latest_app_version})"
                    )
                else:
                    self.btn_upgrade_app.setVisible(False)
            except Exception as e:
                logging.warning(
                    "Plugin Installer: upgrade button version check failed: %s", e
                )
                self.btn_upgrade_app.setVisible(False)

        remote_map = {entry.get("name", ""): entry for entry in self.remote_data}
        installed_plugins = self.main_window.plugin_manager.plugins
        installed_map = {p.get("name"): p for p in installed_plugins}
        # Merge plugins installed this session so the table shows correct status
        # without waiting for discover_plugins (which only runs on close).
        for name, info in self._pending_installs.items():
            if name not in installed_map:
                installed_map[name] = info
        all_names = sorted(set(remote_map.keys()) | set(installed_map.keys()))

        # Collect row descriptors, then sort before rendering
        rows = []
        for name in all_names:
            remote_info = remote_map.get(name)
            local_info = installed_map.get(name)

            local_ver = "-"
            remote_ver = "Unknown"
            status = "Unknown"
            author = "Unknown"
            color = None
            is_installed = local_info is not None

            is_compatible = True
            supported_ver_str = ""
            if remote_info:
                supported_ver_str = remote_info.get("supported_moleditpy_version", "")
                if supported_ver_str:
                    is_compatible = is_app_version_compatible(
                        app_ver, supported_ver_str
                    )

            # Skip non-visible, not-yet-installed plugins
            if remote_info and not is_installed:
                if not remote_info.get("visible", True):
                    continue

            can_update = False
            can_download = False
            target_file = None

            if is_installed:
                target_file = local_info.get("filepath", None)
                if target_file and os.path.isfile(target_file):
                    local_ver = _read_plugin_version_ast(target_file)
                else:
                    local_ver = local_info.get("version", "Unknown")

            if remote_info:
                remote_ver = remote_info.get("version", "Unknown")
                author = remote_info.get("author", "Unknown")

            if author == "Unknown" and is_installed:
                author = local_info.get("author", "Unknown")

            if is_installed:
                if remote_info:
                    if local_ver != "Unknown" and remote_ver != "Unknown":
                        comp = self.compare_versions(remote_ver, local_ver)
                        if comp > 0:
                            status = "Update Available"
                            color = QColor("#f8d7da")
                            can_update = True
                            self.updates_found = True
                            plugin_update_count += 1
                        elif comp == 0:
                            status = "Up to date"
                            color = QColor("#d4edda")
                            can_update = True
                        else:
                            status = "Newer"
                            color = QColor("#cce5ff")
                            can_update = True
                    else:
                        status = "Version Mismatch"
                        can_update = True
                else:
                    status = "Not in Registry"
                    color = QColor("#e2e3e5")
            else:
                if not is_compatible:
                    status = "Incompatible"
                    color = QColor("#fff3cd")
                else:
                    status = "Not Installed"
                can_download = True
                if remote_info and "downloadUrl" in remote_info:
                    d_url = remote_info["downloadUrl"]
                    filename = os.path.basename(d_url)
                    base_dir = os.path.dirname(
                        os.path.dirname(os.path.abspath(__file__))
                    )
                    target_file = os.path.join(base_dir, filename)

            if not is_compatible:
                color = QColor("#fff3cd")

            rows.append(
                dict(
                    name=name,
                    author=author,
                    local_ver=local_ver,
                    remote_ver=remote_ver,
                    status=status,
                    color=color,
                    is_installed=is_installed,
                    is_compatible=is_compatible,
                    supported_ver_str=supported_ver_str,
                    can_update=can_update,
                    can_download=can_download,
                    target_file=target_file,
                    remote_info=remote_info,
                    local_info=local_info,
                )
            )

        rows.sort(key=lambda r: r["name"].lower())

        # Disable repaints during bulk insertion so the table doesn't redraw
        # after every row — avoids visible stutter on 40+ plugins.
        self.table.setUpdatesEnabled(False)
        for row_data in rows:
            row = self.table.rowCount()
            self.table.insertRow(row)
            name = row_data["name"]
            remote_info = row_data["remote_info"]

            self.table.setItem(row, 0, QTableWidgetItem(name))
            self.table.setItem(row, 1, QTableWidgetItem(row_data["author"]))
            self.table.setItem(row, 2, QTableWidgetItem(str(row_data["local_ver"])))
            self.table.setItem(row, 3, QTableWidgetItem(str(row_data["remote_ver"])))

            status_item = QTableWidgetItem(row_data["status"])
            if row_data["color"]:
                status_item.setBackground(row_data["color"])
            if not row_data["is_compatible"] and row_data["supported_ver_str"]:
                status_item.setToolTip(
                    f"Requires MoleditPy {row_data['supported_ver_str']} "
                    f"(Current: {app_ver})"
                )
            self.table.setItem(row, 4, status_item)

            if remote_info and "downloadUrl" in remote_info:
                label = ""
                if row_data["can_update"]:
                    label = (
                        "Update"
                        if row_data["status"] == "Update Available"
                        else "Reinstall"
                    )
                elif row_data["can_download"]:
                    label = "Install"

                if label:
                    btn_action = QPushButton(label)
                    btn_action.setProperty("plugin_name", name)
                    btn_action.setProperty("download_url", remote_info["downloadUrl"])
                    btn_action.setProperty("target_file", row_data["target_file"])
                    btn_action.setProperty(
                        "dependencies", remote_info.get("dependencies", [])
                    )
                    if not row_data["is_installed"] and not row_data["is_compatible"]:
                        btn_action.setToolTip(
                            f"Warning: Incompatible with MoleditPy "
                            f"(Requires: {row_data['supported_ver_str']})"
                        )
                    btn_action.clicked.connect(self.on_update_clicked)
                    self.table.setCellWidget(row, 5, btn_action)
            else:
                self.table.setItem(row, 5, QTableWidgetItem(""))

        self.table.setUpdatesEnabled(True)

        if getattr(self, "btn_update_all", None) is not None:
            self.btn_update_all.setEnabled(plugin_update_count > 0)

        self._update_status_label(plugin_update_count)

        if getattr(self, "search_input", None) is not None:
            self.filter_plugins()

    def _update_status_label(self, plugin_update_count: int) -> None:
        if getattr(self, "_status_label", None) is None:
            return
        if plugin_update_count == 1:
            self._status_label.setStyleSheet("color: #856404; font-weight: bold;")
            self._status_label.setText("1 plugin has an update available")
        elif plugin_update_count > 1:
            self._status_label.setStyleSheet("color: #856404; font-weight: bold;")
            self._status_label.setText(
                f"{plugin_update_count} plugins have updates available"
            )
        elif self.table.rowCount() > 0:
            self._status_label.setStyleSheet("color: #155724; font-weight: bold;")
            self._status_label.setText("All plugins are up to date")
        else:
            self._status_label.setText("")

    # ------------------------------------------------------------------
    # Misc helpers
    # ------------------------------------------------------------------

    def save_settings_ui(self):
        self.settings["check_at_startup"] = self.chk_startup.isChecked()
        save_settings(self.settings)

    def calculate_sha256(self, filepath):
        sha256_hash = hashlib.sha256()
        with open(filepath, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()

    def compare_versions(self, v1, v2):
        try:

            def parse(v):
                return [int(x) for x in re.findall(r"\d+", v)]

            p1 = parse(v1)
            p2 = parse(v2)
            length = max(len(p1), len(p2))
            p1 += [0] * (length - len(p1))
            p2 += [0] * (length - len(p2))
            if p1 > p2:
                return 1
            if p1 < p2:
                return -1
            return 0
        except Exception as e:
            logging.warning(
                "Plugin Installer: version comparison failed (%r vs %r): %s", v1, v2, e
            )
            if v1 > v2:
                return 1
            if v1 < v2:
                return -1
            return 0

    def _get_package_name(self):
        try:
            from moleditpy.utils.constants import VERSION as _ver  # noqa: F401

            return "moleditpy"
        except ImportError:
            pass
        try:
            from moleditpy_linux.utils.constants import VERSION as _ver2  # noqa: F401

            return "moleditpy-linux"
        except ImportError:
            pass
        return "moleditpy"

    def get_app_version(self):
        try:
            from moleditpy.utils.constants import VERSION as APP_VERSION

            return APP_VERSION
        except ImportError:
            try:
                from moleditpy_linux.utils.constants import VERSION as APP_VERSION

                return APP_VERSION
            except ImportError:
                try:
                    main_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
                    if main_dir not in sys.path:
                        sys.path.append(main_dir)
                    try:
                        from moleditpy.utils.constants import VERSION as APP_VERSION
                    except ImportError:
                        from moleditpy_linux.utils.constants import (
                            VERSION as APP_VERSION,
                        )
                    return APP_VERSION
                except Exception as e:
                    logging.warning("Plugin Installer: failed to detect version: %s", e)
                    return "0.0.0"

    def filter_plugins(self):
        text = self.search_input.text().lower()
        for r in range(self.table.rowCount()):
            name_item = self.table.item(r, 0)
            author_item = self.table.item(r, 1)
            name = name_item.text().lower() if name_item else ""
            author = author_item.text().lower() if author_item else ""
            self.table.setRowHidden(r, not ((text in name) or (text in author)))

    # ------------------------------------------------------------------
    # Batch update with a shared progress dialog
    # ------------------------------------------------------------------

    def update_all_plugins(self):
        rows_to_update = []
        for row in range(self.table.rowCount()):
            status_item = self.table.item(row, 4)
            if status_item and status_item.text() == "Update Available":
                btn = self.table.cellWidget(row, 5)
                if btn and isinstance(btn, QPushButton):
                    rows_to_update.append(btn)

        if not rows_to_update:
            QMessageBox.information(
                self, "Update All", "No plugins with available updates found."
            )
            return

        names = [btn.property("plugin_name") for btn in rows_to_update]
        name_list = "\n".join(f"  - {n}" for n in names)
        ret = QMessageBox.question(
            self,
            "Update All",
            f"Update the following {len(names)} plugin(s)?\n\n{name_list}",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes,
        )
        if ret != QMessageBox.StandardButton.Yes:
            return

        progress = QProgressDialog(
            "Preparing updates...", "Cancel", 0, len(rows_to_update), self
        )
        progress.setWindowTitle("Updating Plugins")
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        QApplication.processEvents()

        self._batch_updating = True
        self._batch_progress = progress
        succeeded = []
        failed = []

        try:
            for i, btn in enumerate(rows_to_update):
                if progress.wasCanceled():
                    break
                name = btn.property("plugin_name")
                progress.setLabelText(
                    f"Updating {name}  ({i + 1} / {len(rows_to_update)})"
                )
                progress.setValue(i)
                QApplication.processEvents()

                try:
                    self._last_install_succeeded = False
                    btn.click()
                    if self._last_install_succeeded:
                        succeeded.append(name)
                    else:
                        failed.append(f"{name} (Skipped/Cancelled)")
                except Exception as e:
                    failed.append(f"{name}: {e}")

            progress.setValue(len(rows_to_update))
        finally:
            self._batch_updating = False
            self._batch_progress = None
            progress.close()
            # discover_plugins / _refresh_plugin_menus run once in closeEvent.
            # _pending_installs was populated per-plugin inside on_update_clicked.
            QApplication.processEvents()
            self.populate_table()

        lines = [f"Updated {len(succeeded)} of {len(names)} plugin(s)."]
        if succeeded:
            lines.append("\nSucceeded:\n" + "\n".join(f"  - {n}" for n in succeeded))
        if failed:
            lines.append(
                "\nFailed / Skipped:\n" + "\n".join(f"  - {n}" for n in failed)
            )
        QMessageBox.information(self, "Update All Complete", "\n".join(lines))

    # ------------------------------------------------------------------
    # Details dialog
    # ------------------------------------------------------------------

    def show_plugin_details(self, row, col):
        name_item = self.table.item(row, 0)
        if not name_item:
            return
        name = name_item.text()

        remote_info = next((e for e in self.remote_data if e.get("name") == name), None)
        local_info = next(
            (
                p
                for p in self.main_window.plugin_manager.plugins
                if p.get("name") == name
            ),
            None,
        )

        description = "No description available."
        author = "Unknown"
        version = "Unknown"
        dependencies = []
        target_file = None

        if remote_info:
            description = remote_info.get("description", description)
            author = remote_info.get("author", author)
            version = remote_info.get("version", version)
            dependencies = remote_info.get("dependencies", [])

        if local_info:
            if description == "No description available.":
                description = local_info.get("description", description)
            if author == "Unknown":
                author = local_info.get("author", author)
            local_ver = local_info.get("version", "Unknown")
            target_file = local_info.get("filepath", None)
            if remote_info:
                version = f"Installed: {local_ver}\nLatest: {version}"
            else:
                version = f"Installed: {local_ver}"

        supported_version = (
            remote_info.get("supported_moleditpy_version", "Unknown")
            if remote_info
            else "Unknown"
        )

        dialog = PluginDetailsDialog(
            self,
            name,
            author,
            version,
            description,
            dependencies,
            local_info,
            target_file,
            supported_version,
        )
        dialog.exec()

    # ------------------------------------------------------------------
    # Install / update click handler
    # ------------------------------------------------------------------

    def on_update_clicked(self):
        btn = self.sender()
        if not btn:
            return

        plugin_name = btn.property("plugin_name")
        download_url = btn.property("download_url")
        target_file = btn.property("target_file")
        dependencies = btn.property("dependencies") or []
        missing_deps = []
        user_confirmed_intent = False
        self._last_install_succeeded = False

        is_installed = bool(target_file and os.path.exists(target_file))

        remote_info = next(
            (e for e in self.remote_data if e.get("name") == plugin_name), None
        )

        # App version compatibility check
        if remote_info:
            supported_ver_str = remote_info.get("supported_moleditpy_version", "")
            app_ver = self.get_app_version()
            if not is_app_version_compatible(app_ver, supported_ver_str):
                action_word = "update" if is_installed else "install"
                ret = QMessageBox.warning(
                    self,
                    "Incompatible Version Warning",
                    f"Warning: '{plugin_name}' requires MoleditPy version "
                    f"'{supported_ver_str}', but you are running version "
                    f"'{app_ver}'.\n\nThis plugin may fail to function correctly.\n\n"
                    f"Do you want to proceed with the {action_word} anyway?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.No,
                )
                if ret != QMessageBox.StandardButton.Yes:
                    return
                user_confirmed_intent = True

        # Dependency check
        if dependencies:
            missing_deps = [
                d for d in dependencies if not check_dependency_satisfied(d)
            ]
            if missing_deps:
                local_info = next(
                    (
                        p
                        for p in self.main_window.plugin_manager.plugins
                        if p.get("name") == plugin_name
                    ),
                    None,
                )
                description = remote_info.get("description", "") if remote_info else ""
                author = (
                    remote_info.get("author", "Unknown") if remote_info else "Unknown"
                )
                version = (
                    remote_info.get("version", "Unknown") if remote_info else "Unknown"
                )
                if local_info:
                    local_ver = local_info.get("version", "Unknown")
                    version = f"Installed: {local_ver}\nLatest: {version}"

                ret = QMessageBox.question(
                    self,
                    "Missing Dependencies",
                    f"The following dependencies are missing: {', '.join(missing_deps)}.\n\n"
                    f"The plugin may not work without them.\nDo you want to install it anyway?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                )
                if ret == QMessageBox.StandardButton.No:
                    supported_version = (
                        remote_info.get("supported_moleditpy_version", "Unknown")
                        if remote_info
                        else "Unknown"
                    )
                    dialog = PluginDetailsDialog(
                        self,
                        plugin_name,
                        author,
                        version,
                        description,
                        dependencies,
                        local_info,
                        target_file,
                        supported_version,
                    )
                    dialog.exec()
                    return
                user_confirmed_intent = True

        if not download_url:
            QMessageBox.warning(self, "Error", "Cannot update: Missing download URL.")
            return

        is_installed = bool(target_file and os.path.exists(target_file))
        action_verb = "Update" if is_installed else "Install"

        if not user_confirmed_intent and not self._batch_updating:
            ret = QMessageBox.question(
                self,
                f"{action_verb} Plugin",
                f"{action_verb} '{plugin_name}'?\n"
                f"This will download and install the plugin.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if ret != QMessageBox.StandardButton.Yes:
                return

        # Resolve URL
        final_url = download_url
        if not download_url.startswith("http"):
            base_url = REMOTE_JSON_URL.rsplit("/", 1)[0]
            final_url = urllib.parse.urljoin(base_url + "/", download_url)

        logging.info(
            "Plugin Installer: %sing %s from: %s", action_verb, plugin_name, final_url
        )

        temp_dir = tempfile.mkdtemp()
        try:
            # Determine filename
            filename = None
            if target_file and os.path.exists(target_file):
                file_name_existing = os.path.basename(target_file)
                if (
                    file_name_existing != "__init__.py"
                    and file_name_existing.endswith(".py")
                    and not final_url.lower().endswith(".zip")
                ):
                    filename = file_name_existing

            if not filename:
                path = urllib.parse.urlparse(final_url).path
                filename = os.path.basename(path)

            if not filename:
                filename = (
                    "plugin_update.zip"
                    if final_url.lower().endswith(".zip")
                    else "plugin_update.py"
                )

            download_path = os.path.join(temp_dir, filename)
            dl_label = (
                f"{'Updating' if is_installed else 'Installing'} {plugin_name}..."
            )

            # ---- Download with progress ----
            if self._batch_updating:
                # Batch mode: update the shared dialog label; no own dialog created,
                # so the batch dialog's range (0..N plugins) is never touched.
                batch_prog = self._batch_progress

                def _batch_cb(recv, total):
                    if batch_prog is None or batch_prog.wasCanceled():
                        return False
                    kb = recv // 1024
                    suffix = f" / {total // 1024} KB" if total else " KB"
                    batch_prog.setLabelText(
                        f"Downloading {plugin_name}...\n{kb}{suffix}"
                    )
                    QApplication.processEvents()
                    return not batch_prog.wasCanceled()

                ok = self._download_chunked(final_url, download_path, _batch_cb)
            else:
                # Single mode: own QProgressDialog with byte-level progress
                prog = QProgressDialog(dl_label, "Cancel", 0, 0, self)
                prog.setWindowTitle("Downloading")
                prog.setWindowModality(Qt.WindowModality.WindowModal)
                prog.setMinimumDuration(0)
                prog.setValue(0)
                QApplication.processEvents()

                def _single_cb(recv, total):
                    if total > 0:
                        prog.setRange(0, total)
                        prog.setValue(recv)
                    kb = recv // 1024
                    suffix = f" / {total // 1024} KB" if total else " KB"
                    prog.setLabelText(f"{dl_label}\n{kb}{suffix}")
                    QApplication.processEvents()
                    return not prog.wasCanceled()

                ok = self._download_chunked(final_url, download_path, _single_cb)
                prog.close()

            if not ok:
                return  # Cancelled or network error (already shown to user)

            QApplication.processEvents()

            # ---- SHA256 verification ----
            remote_sha256 = remote_info.get("sha256") if remote_info else None

            if not remote_sha256:
                ret = QMessageBox.warning(
                    self,
                    "Security Warning",
                    f"SHA256 checksum is not available for '{plugin_name}'.\n\n"
                    f"The plugin file's integrity cannot be verified.\n\n"
                    f"Do you want to proceed anyway?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.No,
                )
                if ret != QMessageBox.StandardButton.Yes:
                    return
            else:
                calculated_sha256 = self.calculate_sha256(download_path)
                QApplication.processEvents()
                if calculated_sha256 != remote_sha256:
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

            # ---- Install / update ----
            try:
                did_manual_overwrite = False

                if target_file and os.path.exists(target_file):
                    file_name_existing = os.path.basename(target_file)
                    if (
                        file_name_existing != "__init__.py"
                        and file_name_existing.endswith(".py")
                        and not final_url.lower().endswith(".zip")
                    ):
                        try:
                            logging.info(
                                "Plugin Installer: overwriting %s", target_file
                            )
                            shutil.copy2(download_path, target_file)
                            did_manual_overwrite = True
                        except Exception as e:
                            logging.warning(
                                "Plugin Installer: manual overwrite failed: %s "
                                "— falling back to manager",
                                e,
                            )

                    elif (
                        file_name_existing == "__init__.py"
                        and final_url.lower().endswith(".zip")
                    ):
                        try:
                            target_dir = os.path.dirname(target_file)
                            logging.info(
                                "Plugin Installer: overwriting folder plugin %s",
                                target_dir,
                            )
                            extract_temp = os.path.join(temp_dir, "extracted")
                            os.makedirs(extract_temp, exist_ok=True)
                            with zipfile.ZipFile(download_path, "r") as z:
                                for member in z.infolist():
                                    target_path = os.path.join(
                                        extract_temp, member.filename
                                    )
                                    if not os.path.abspath(target_path).startswith(
                                        os.path.abspath(extract_temp)
                                    ):
                                        logging.warning(
                                            "Plugin Installer: skipping suspicious "
                                            "zip entry: %s",
                                            member.filename,
                                        )
                                        continue
                                    z.extract(member, extract_temp)

                            items = os.listdir(extract_temp)
                            if len(items) == 1 and os.path.isdir(
                                os.path.join(extract_temp, items[0])
                            ):
                                source = os.path.join(extract_temp, items[0])
                            else:
                                source = extract_temp

                            settings_backup = None
                            settings_path = os.path.join(target_dir, "settings.json")
                            if os.path.exists(settings_path):
                                try:
                                    with open(
                                        settings_path, "r", encoding="utf-8"
                                    ) as f:
                                        settings_backup = f.read()
                                except Exception as e:
                                    logging.warning(
                                        "Plugin Installer: failed to backup "
                                        "settings.json: %s",
                                        e,
                                    )

                            shutil.copytree(source, target_dir, dirs_exist_ok=True)
                            QApplication.processEvents()

                            if settings_backup:
                                try:
                                    with open(
                                        settings_path, "w", encoding="utf-8"
                                    ) as f:
                                        f.write(settings_backup)
                                    logging.info(
                                        "Plugin Installer: restored settings.json"
                                    )
                                except Exception as e:
                                    logging.warning(
                                        "Plugin Installer: failed to restore "
                                        "settings.json: %s",
                                        e,
                                    )

                            did_manual_overwrite = True
                        except Exception as e:
                            logging.warning(
                                "Plugin Installer: manual folder overwrite failed: %s"
                                " — falling back to manager",
                                e,
                            )

                if not did_manual_overwrite:
                    saved_settings_content = None
                    target_dir_for_restore = None

                    if target_file and os.path.exists(target_file):
                        if os.path.basename(target_file) == "__init__.py":
                            folder_dir = os.path.dirname(target_file)
                            settings_json = os.path.join(folder_dir, "settings.json")
                            if os.path.exists(settings_json):
                                try:
                                    with open(
                                        settings_json, "r", encoding="utf-8"
                                    ) as f:
                                        saved_settings_content = f.read()
                                    target_dir_for_restore = folder_dir
                                except Exception as e:
                                    logging.warning(
                                        "Plugin Installer: failed to backup "
                                        "settings.json: %s",
                                        e,
                                    )

                    QApplication.processEvents()
                    self.main_window.plugin_manager.install_plugin(download_path)
                    QApplication.processEvents()

                    if saved_settings_content and target_dir_for_restore:
                        if os.path.exists(target_dir_for_restore):
                            settings_json = os.path.join(
                                target_dir_for_restore, "settings.json"
                            )
                            try:
                                with open(settings_json, "w", encoding="utf-8") as f:
                                    f.write(saved_settings_content)
                                logging.info("Plugin Installer: restored settings.json")
                            except Exception as e:
                                logging.warning(
                                    "Plugin Installer: failed to restore "
                                    "settings.json: %s",
                                    e,
                                )

                    if target_file and os.path.exists(target_file):
                        _, new_ext = os.path.splitext(filename)
                        new_is_zip = new_ext.lower() == ".zip"
                        old_is_py = target_file.lower().endswith(".py")
                        old_is_init = os.path.basename(target_file) == "__init__.py"

                        if new_is_zip and old_is_py and not old_is_init:
                            try:
                                os.remove(target_file)
                            except Exception as e:
                                logging.warning(
                                    "Plugin Installer: failed to remove old .py: %s",
                                    e,
                                )
                        elif not new_is_zip and old_is_init:
                            try:
                                parent_dir = os.path.dirname(target_file)
                                shutil.rmtree(parent_dir)
                            except Exception as e:
                                logging.warning(
                                    "Plugin Installer: failed to remove old folder: %s",
                                    e,
                                )

                verb_past = "updated" if is_installed else "installed"
                self._last_install_succeeded = True

                # Record installed plugin for immediate table display.
                # discover_plugins / _refresh_plugin_menus run once in closeEvent.
                installed_ver = (
                    remote_info.get("version", "Unknown") if remote_info else "Unknown"
                )
                installed_author = (
                    remote_info.get("author", "Unknown") if remote_info else "Unknown"
                )
                self._pending_installs[plugin_name] = {
                    "name": plugin_name,
                    "version": installed_ver,
                    "filepath": target_file or "",
                    "author": installed_author,
                }
                self._needs_plugin_reload = True

                if missing_deps:
                    QMessageBox.warning(
                        self,
                        "Installation Complete (Dependencies Missing)",
                        f"Successfully {verb_past} '{plugin_name}'.\n\n"
                        f"However, you must install the missing dependencies for it "
                        f"to work.\n\nOpening details window...",
                    )
                    local_info = next(
                        (
                            p
                            for p in self.main_window.plugin_manager.plugins
                            if p.get("name") == plugin_name
                        ),
                        None,
                    )
                    description = (
                        remote_info.get("description", "") if remote_info else ""
                    )
                    author = (
                        remote_info.get("author", "Unknown")
                        if remote_info
                        else "Unknown"
                    )
                    version = (
                        remote_info.get("version", "Unknown")
                        if remote_info
                        else "Unknown"
                    )
                    supported_version = (
                        remote_info.get("supported_moleditpy_version", "Unknown")
                        if remote_info
                        else "Unknown"
                    )
                    dialog = PluginDetailsDialog(
                        self,
                        plugin_name,
                        author,
                        version,
                        description,
                        dependencies,
                        local_info,
                        target_file,
                        supported_version,
                    )
                    dialog.exec()
                elif not self._batch_updating:
                    QMessageBox.information(
                        self,
                        "Success",
                        f"Successfully {verb_past} '{plugin_name}'.",
                    )

                if not self._batch_updating:
                    self.populate_table()

            except Exception as e:
                logging.exception(
                    "Plugin Installer: installation error for %s: %s", plugin_name, e
                )
                QMessageBox.warning(
                    self,
                    "Installation Error",
                    f"Failed to install plugin:\n{e}",
                )

        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Close
    # ------------------------------------------------------------------

    def on_close_clicked(self):
        self.close()

    def _on_finished(self, _result: int) -> None:
        """Runs after the dialog closes regardless of how it was closed."""
        if not (
            self._needs_plugin_reload
            and self.main_window
            and hasattr(self.main_window, "plugin_manager")
        ):
            return

        try:
            self.main_window.plugin_manager.discover_plugins(self.main_window)
            _refresh_plugin_menus(self.main_window)
        except Exception as e:
            logging.warning("Plugin Installer: reload on close failed: %s", e)
        finally:
            self._pending_installs.clear()
            self._needs_plugin_reload = False

        for widget in QApplication.topLevelWidgets():
            if type(widget).__name__ == "PluginManagerWindow":
                if hasattr(widget, "refresh_plugin_list"):
                    try:
                        widget.refresh_plugin_list()
                    except Exception as e:
                        logging.warning(
                            "Plugin Installer: failed to sync PluginManagerWindow: %s",
                            e,
                        )

    def closeEvent(self, event):
        # Stop any in-flight fetch thread before closing
        if self._fetch_worker and self._fetch_worker.isRunning():
            try:
                self._fetch_worker.done.disconnect()
            except Exception:
                pass
            self._fetch_worker.quit()
            self._fetch_worker.wait(500)

        super().closeEvent(event)

    def copy_upgrade_command(self):
        package_name = self._get_package_name()
        cmd = f"pip install --upgrade {package_name}"
        QApplication.clipboard().setText(cmd)
        QMessageBox.information(
            self,
            "Upgrade Command",
            f"Command copied to clipboard:\n\n{cmd}\n\n"
            f"Please close MoleditPy and run this command in your terminal.",
        )
