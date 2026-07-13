"""
Metadata Saver — MoleditPy Plugin
==================================
Captures debug/trace metadata (save timestamp, file path, username, OS, …)
and embeds it into the .pmeprj project file via the save hook.

⚠  This data is personal. Do NOT share project files containing
   this metadata with anonymous users.

Author : HiroYokoyama
License: GPL-3.0
"""

from __future__ import annotations

import getpass
import json
import logging
import os
import platform
import socket
from datetime import datetime, timezone

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

# ---------------------------------------------------------------------------
# Plugin metadata
# ---------------------------------------------------------------------------

PLUGIN_NAME = "Metadata Saver"
PLUGIN_VERSION = "2026.07.13"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Saves debug/trace metadata (timestamp, file path, username, OS info, …) "
    "into the .pmeprj project file. For debugging and tracing purposes only — "
    "do not share project files containing this data with anonymous users."
)
PLUGIN_CATEGORY = "Utility"
PLUGIN_TAGS = ["Debug", "Tracing", "Utility"]
PLUGIN_DEPENDENCIES = []
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"

# ---------------------------------------------------------------------------
# Module-level state
# ---------------------------------------------------------------------------

_CONTEXT = None
_LOADED_METADATA: dict = {}  # metadata unpacked from the last loaded .pmeprj

# ---------------------------------------------------------------------------
# Field registry — ordered list of all capturable fields
# ---------------------------------------------------------------------------

ALL_FIELDS: list[dict] = [
    {
        "key": "saved_at",
        "label": "Save Timestamp",
        "description": "ISO 8601 datetime when the project was saved (UTC).",
        "default": True,
    },
    {
        "key": "saved_path",
        "label": "Save File Path",
        "description": "Full absolute path of the saved project file.",
        "default": True,
    },
    {
        "key": "username",
        "label": "Username",
        "description": "OS login name of the user who saved the file.",
        "default": True,
    },
    {
        "key": "hostname",
        "label": "Hostname",
        "description": "Network hostname of the machine.",
        "default": True,
    },
    {
        "key": "os_name",
        "label": "OS Name",
        "description": "Operating system name (e.g. Windows, Linux, Darwin).",
        "default": True,
    },
    {
        "key": "os_version",
        "label": "OS Version",
        "description": "Detailed OS version string.",
        "default": True,
    },
    {
        "key": "os_release",
        "label": "OS Release",
        "description": "OS release identifier (e.g. 10, 11, 22.04).",
        "default": True,
    },
    {
        "key": "machine",
        "label": "Machine Architecture",
        "description": "CPU architecture (e.g. AMD64, x86_64, arm64).",
        "default": True,
    },
    {
        "key": "platform_full",
        "label": "Full Platform String",
        "description": "Complete platform identification string.",
        "default": True,
    },
    {
        "key": "python_version",
        "label": "Python Version",
        "description": "Python interpreter version.",
        "default": True,
    },
    {
        "key": "app_version",
        "label": "App Version",
        "description": "MoleditPy application version (if available).",
        "default": True,
    },
    {
        "key": "installed_plugins",
        "label": "Installed Plugins",
        "description": "List of currently installed plugins and their versions.",
        "default": False,
    },
    {
        "key": "note",
        "label": "Custom Note",
        "description": "User-defined free-text note attached to each save.",
        "default": False,
    },
]

# ---------------------------------------------------------------------------
# Config persistence  (companion JSON, survives plugin updates)
# ---------------------------------------------------------------------------

_CONFIG_FILENAME = "metadata_saver.json"


def _config_path() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), _CONFIG_FILENAME)


def _default_config() -> dict:
    return {
        "enabled": True,
        "silence_notice": False,
        "enabled_fields": {f["key"]: f["default"] for f in ALL_FIELDS},
        "custom_note": "",
    }


def load_config() -> dict:
    """Load plugin config from companion JSON, merging with defaults."""
    path = _config_path()
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
            cfg = _default_config()
            cfg["enabled"] = bool(data.get("enabled", True))
            cfg["silence_notice"] = bool(data.get("silence_notice", False))
            if "enabled_fields" in data and isinstance(data["enabled_fields"], dict):
                cfg["enabled_fields"].update(data["enabled_fields"])
            cfg["custom_note"] = data.get("custom_note", "")
            return cfg
        except Exception:
            logging.debug(
                "Metadata Saver: failed to load config — using defaults", exc_info=True
            )
    return _default_config()


def save_config(cfg: dict) -> None:
    """Persist plugin config to companion JSON."""
    try:
        with open(_config_path(), "w", encoding="utf-8") as fh:
            json.dump(cfg, fh, indent=4, ensure_ascii=False)
    except Exception:
        logging.warning("Metadata Saver: failed to save config", exc_info=True)


# ---------------------------------------------------------------------------
# Metadata collection helpers
# ---------------------------------------------------------------------------

def _try_getattr_chain(obj, *attr_paths: str) -> str | None:
    """
    Try a list of dotted attribute paths on *obj* in order.
    Returns the first non-empty string value found, or None.
    """
    for attr_path in attr_paths:
        try:
            current = obj
            for part in attr_path.split("."):
                current = getattr(current, part, None)
                if current is None:
                    break
            if current and isinstance(current, str) and current.strip():
                return current
        except Exception:
            pass
    return None


def _collect_field(key: str, mw) -> str:
    """Collect the value for a single metadata field. Never raises."""
    try:
        if key == "saved_at":
            return datetime.now(timezone.utc).isoformat()

        if key == "saved_path":
            if mw is not None:
                val = None
                if hasattr(mw, "get_current_file_path"):
                    val = mw.get_current_file_path()
                if not val and hasattr(mw, "init_manager") and hasattr(mw.init_manager, "current_file_path"):
                    val = mw.init_manager.current_file_path
                
                if val:
                    return os.path.abspath(val)
            return "(unknown — save path not yet available)"

        if key == "username":
            return getpass.getuser()

        if key == "hostname":
            return socket.gethostname()

        if key == "os_name":
            return platform.system()

        if key == "os_version":
            return platform.version()

        if key == "os_release":
            return platform.release()

        if key == "machine":
            return platform.machine()

        if key == "platform_full":
            return platform.platform()

        if key == "python_version":
            return platform.python_version()

        if key == "app_version":
            try:
                from moleditpy.utils.constants import VERSION
                return str(VERSION)
            except ImportError:
                pass
            try:
                from moleditpy_linux.utils.constants import VERSION
                return str(VERSION)
            except ImportError:
                pass
            return "(unknown)"

        if key == "installed_plugins":
            if mw is not None and hasattr(mw, "plugin_manager"):
                plugins = getattr(mw.plugin_manager, "plugins", [])
                if plugins:
                    return ", ".join(f"{p.get('name', 'Unknown')} v{p.get('version', 'Unknown')}" for p in plugins)
            return "(none/unknown)"

    except Exception:
        logging.debug(
            "Metadata Saver: error collecting field '%s'", key, exc_info=True
        )

    return "(error)"


def collect_metadata(cfg: dict, mw) -> dict:
    """Return a dict of all enabled metadata fields with their current values."""
    enabled: dict = cfg.get("enabled_fields", {})
    result: dict = {}

    for field in ALL_FIELDS:
        key = field["key"]
        if not enabled.get(key, field["default"]):
            continue

        if key == "note":
            note_text = cfg.get("custom_note", "").strip()
            if note_text:
                result[key] = note_text
            # skip if note is empty
        else:
            result[key] = _collect_field(key, mw)

    return result


# ---------------------------------------------------------------------------
# Plugin entry point
# ---------------------------------------------------------------------------

def initialize(context) -> None:
    """Register menu action and project persistence hooks."""
    global _CONTEXT
    _CONTEXT = context

    context.add_menu_action(
        "Settings/Metadata Saver...",
        lambda: _open_dialog(context),
    )
    context.register_save_handler(on_save_project)
    context.register_load_handler(on_load_project)
    context.register_document_reset_handler(on_document_reset)

    cfg = load_config()
    if not cfg.get("silence_notice", False):
        logging.info("Metadata Saver: initialized (v%s).", PLUGIN_VERSION)


# ---------------------------------------------------------------------------
# Save / load / reset handlers
# ---------------------------------------------------------------------------

def on_save_project() -> dict | None:
    """
    Called by the host when the user saves a .pmeprj project.
    Returns a dict that is embedded in the project file under this plugin's key.

    ⚠  The returned data contains personal information (username, file path, …).
       The resulting .pmeprj file should NOT be shared with anonymous users.
    """
    try:
        cfg = load_config()

        # Master switch — if disabled, do not embed anything
        if not cfg.get("enabled", True):
            return None

        mw = _CONTEXT.get_main_window() if _CONTEXT is not None else None
        metadata = collect_metadata(cfg, mw)

        if not metadata:
            return None

        if not cfg.get("silence_notice", False):
            logging.info(
                "Metadata Saver: embedding %d metadata field(s) in project. "
                "NOTE — this data is personal; do not share with anonymous users.",
                len(metadata),
            )
        return {
            "metadata": metadata,
            "plugin_version": PLUGIN_VERSION,
        }
    except Exception:
        logging.warning("Metadata Saver: on_save_project failed", exc_info=True)
        return None


def on_load_project(data: dict) -> None:
    """Restore metadata from a loaded .pmeprj file."""
    global _LOADED_METADATA
    try:
        if isinstance(data, dict) and "metadata" in data:
            _LOADED_METADATA = dict(data["metadata"])
            cfg = load_config()
            if not cfg.get("silence_notice", False):
                logging.info(
                    "Metadata Saver: loaded %d metadata field(s) from project.",
                    len(_LOADED_METADATA),
                )
        else:
            _LOADED_METADATA = {}
    except Exception:
        logging.debug("Metadata Saver: on_load_project failed", exc_info=True)
        _LOADED_METADATA = {}


def on_document_reset() -> None:
    """Clear loaded metadata when a new document is created (File > New)."""
    global _LOADED_METADATA
    _LOADED_METADATA = {}


# ---------------------------------------------------------------------------
# Dialog helpers
# ---------------------------------------------------------------------------

def _open_dialog(context) -> None:
    """Open the settings dialog (singleton pattern)."""
    win = context.get_window("settings_dialog")
    if win is not None and not win.isHidden():
        win.raise_()
        win.activateWindow()
        return

    mw = context.get_main_window()
    dlg = MetadataSaverDialog(context, mw)
    context.register_window("settings_dialog", dlg)
    dlg.show()


# ---------------------------------------------------------------------------
# Settings dialog
# ---------------------------------------------------------------------------

class MetadataSaverDialog(QDialog):
    """
    Settings dialog for the Metadata Saver plugin.

    Left column  : checkboxes to select which fields are captured on save.
    Right area   : read-only viewer showing metadata embedded in the current project,
                   plus a "Preview" button to see what would be saved right now.
    """

    def __init__(self, context, parent=None):
        super().__init__(parent)
        self.context = context
        self.setWindowTitle("Metadata Saver — Settings")
        self.setMinimumSize(560, 480)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, False)

        self._cfg: dict = load_config()
        self._checkboxes: dict[str, QCheckBox] = {}
        self._note_edit: QTextEdit | None = None

        self._build_ui()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(10)
        root.setContentsMargins(12, 12, 12, 12)

        # --- Master enable toggle (very top) ---
        self._chk_master = QCheckBox("Enable Metadata Saver (capture metadata on every Save)")
        self._chk_master.setChecked(self._cfg.get("enabled", True))
        self._chk_master.setStyleSheet("font-weight: bold; font-size: 12px;")
        self._chk_master.toggled.connect(self._on_master_toggled)
        root.addWidget(self._chk_master)

        self._chk_silence = QCheckBox("Silence privacy notice in console on save")
        self._chk_silence.setChecked(self._cfg.get("silence_notice", False))
        root.addWidget(self._chk_silence)

        # Privacy banner
        banner = QLabel(
            "⚠  This data is personal (username, hostname, file path).\n"
            "Do NOT share project files containing this metadata with anonymous users."
        )
        banner.setWordWrap(True)
        banner.setStyleSheet(
            "color: #e8a000; font-size: 11px; padding: 6px 8px;"
            "border: 1px solid #7a5500; border-radius: 4px;"
            "background-color: rgba(232, 160, 0, 0.08);"
        )
        root.addWidget(banner)

        # --- Fields group -----------------------------------------------
        fields_group = QGroupBox("Fields to capture on each Save")
        fields_layout = QVBoxLayout(fields_group)
        fields_layout.setSpacing(5)

        enabled_map: dict = self._cfg.get("enabled_fields", {})

        for field in ALL_FIELDS:
            key = field["key"]
            chk = QCheckBox(field["label"])
            chk.setToolTip(field["description"])
            chk.setChecked(enabled_map.get(key, field["default"]))
            self._checkboxes[key] = chk

            if key == "note":
                note_row = QHBoxLayout()
                note_row.setSpacing(8)
                note_row.addWidget(chk)

                self._note_edit = QTextEdit()
                self._note_edit.setPlaceholderText("Enter custom note text here…")
                self._note_edit.setMaximumHeight(48)
                self._note_edit.setSizePolicy(
                    QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
                )
                self._note_edit.setText(self._cfg.get("custom_note", ""))
                note_row.addWidget(self._note_edit)
                fields_layout.addLayout(note_row)
            else:
                fields_layout.addWidget(chk)

        # Wrap fields in a container that can be disabled en masse
        self._fields_container = QWidget()
        self._fields_container.setLayout(fields_layout)
        self._fields_container.setEnabled(self._cfg.get("enabled", True))

        scroll_widget = QWidget()
        container_outer = QVBoxLayout(scroll_widget)
        container_outer.setContentsMargins(0, 0, 0, 0)
        container_outer.addWidget(self._fields_container)
        scroll_area = QScrollArea()
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll_area.setEnabled(self._cfg.get("enabled", True))

        root.addWidget(scroll_area, stretch=2)

        # --- Metadata viewer group -------------------------------------
        viewer_group = QGroupBox("Metadata in current project")
        viewer_layout = QVBoxLayout(viewer_group)
        viewer_layout.setSpacing(6)

        self._metadata_view = QTextEdit()
        self._metadata_view.setReadOnly(True)
        self._metadata_view.setMaximumHeight(160)
        self._metadata_view.setFontFamily("Consolas")
        self._metadata_view.setFontPointSize(9)
        self._metadata_view.setPlaceholderText(
            "(no metadata in current project — save the project first)"
        )
        viewer_layout.addWidget(self._metadata_view)

        btn_row = QHBoxLayout()
        btn_preview = QPushButton("Preview current capture")
        btn_preview.setToolTip(
            "Show what would be saved right now based on current settings.\n"
            "Does not modify the project file."
        )
        btn_preview.clicked.connect(self._on_preview)
        btn_row.addWidget(btn_preview)
        btn_row.addStretch()
        viewer_layout.addLayout(btn_row)

        root.addWidget(viewer_group, stretch=1)

        # --- Dialog buttons -------------------------------------------
        btn_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save
            | QDialogButtonBox.StandardButton.Close
        )
        btn_box.button(QDialogButtonBox.StandardButton.Save).setText("Save Settings")
        btn_box.accepted.connect(self._on_save_settings)
        btn_box.rejected.connect(self.close)
        root.addWidget(btn_box)

        # Populate viewer with metadata already loaded from project
        self._refresh_viewer_from_loaded()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _refresh_viewer_from_loaded(self) -> None:
        """Fill the viewer with metadata loaded from the current .pmeprj."""
        if _LOADED_METADATA:
            lines = [f"{k}: {v}" for k, v in _LOADED_METADATA.items()]
            self._metadata_view.setPlainText("\n".join(lines))
        else:
            self._metadata_view.setPlainText("")

    def _on_master_toggled(self, checked: bool) -> None:
        """Grey out / restore the fields section when master toggle changes."""
        self._fields_container.setEnabled(checked)

    def _build_cfg_from_ui(self) -> dict:
        """Assemble a config dict from the current UI state."""
        enabled = {k: chk.isChecked() for k, chk in self._checkboxes.items()}
        custom_note = (
            self._note_edit.toPlainText().strip()
            if self._note_edit is not None
            else ""
        )
        return {
            "enabled": self._chk_master.isChecked(),
            "silence_notice": self._chk_silence.isChecked(),
            "enabled_fields": enabled,
            "custom_note": custom_note,
        }

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _on_preview(self) -> None:
        """Show a live preview of what would be saved right now."""
        cfg = self._build_cfg_from_ui()
        mw = self.context.get_main_window() if self.context else None
        meta = collect_metadata(cfg, mw)
        if not meta:
            self._metadata_view.setPlainText("(no fields selected)")
            return
        lines = [f"{k}: {v}" for k, v in meta.items()]
        self._metadata_view.setPlainText("\n".join(lines))

    def _on_save_settings(self) -> None:
        """Persist the current UI state to the companion JSON."""
        cfg = self._build_cfg_from_ui()
        save_config(cfg)
        self._cfg = cfg
        self.context.show_status_message(
            "Metadata Saver: settings saved.", 3000
        )
        QMessageBox.information(
            self,
            "Settings Saved",
            "Settings saved.\n\n"
            "Metadata will be captured on the next project save (File > Save).",
        )
