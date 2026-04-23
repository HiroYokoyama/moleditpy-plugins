"""
Unit tests for the Plugin Installer plugin.

All tests run headlessly — PyQt6 and network are fully mocked.
"""

import json
import sys
import types
import importlib
from pathlib import Path
from unittest.mock import MagicMock, patch, call
import pytest

# ---------------------------------------------------------------------------
# Helpers — load the module without executing Qt or network code at import
# ---------------------------------------------------------------------------

PLUGIN_PATH = (
    Path(__file__).resolve().parents[1]
    / "plugins" / "Plugin_Installer" / "plugin_installer.py"
)


def _load_module():
    """Import plugin_installer with PyQt6 stubbed out."""
    # Provide a minimal PyQt6 stub so the top-level import succeeds
    if "PyQt6" not in sys.modules:
        pyqt6 = types.ModuleType("PyQt6")
        pyqt6.QtWidgets = types.ModuleType("PyQt6.QtWidgets")
        pyqt6.QtGui = types.ModuleType("PyQt6.QtGui")
        for name in [
            "QDialog", "QVBoxLayout", "QHBoxLayout", "QLabel",
            "QTableWidget", "QTableWidgetItem", "QPushButton",
            "QHeaderView", "QMessageBox", "QAbstractItemView",
            "QApplication", "QCheckBox", "QLineEdit",
        ]:
            setattr(pyqt6.QtWidgets, name, MagicMock())
        setattr(pyqt6.QtGui, "QColor", MagicMock())
        sys.modules["PyQt6"] = pyqt6
        sys.modules["PyQt6.QtWidgets"] = pyqt6.QtWidgets
        sys.modules["PyQt6.QtGui"] = pyqt6.QtGui

    spec = importlib.util.spec_from_file_location("plugin_installer", PLUGIN_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PI = _load_module()


# ---------------------------------------------------------------------------
# compare_versions — mirrors PluginInstallerWindow.compare_versions exactly
# (tested standalone because the class inherits from a mocked QDialog)
# ---------------------------------------------------------------------------

import re as _re


def _compare_versions(v1, v2):
    """Mirrors PluginInstallerWindow.compare_versions for unit testing."""
    try:
        def parse(v):
            return [int(x) for x in _re.findall(r'\d+', v)]
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
    except Exception:
        if v1 > v2:
            return 1
        if v1 < v2:
            return -1
        return 0


class TestCompareVersions:
    def test_greater(self):
        assert _compare_versions("2.10.0", "2.5.0") == 1

    def test_less(self):
        assert _compare_versions("1.0.0", "1.0.1") == -1

    def test_equal(self):
        assert _compare_versions("2026.04.13", "2026.04.13") == 0

    def test_different_length(self):
        assert _compare_versions("2.5.1", "2.5") == 1
        assert _compare_versions("2.5", "2.5.1") == -1

    def test_date_versions(self):
        assert _compare_versions("2026.04.23", "2026.04.13") == 1
        assert _compare_versions("2026.04.13", "2026.04.23") == -1

    def test_equal_date_versions(self):
        assert _compare_versions("2026.04.23", "2026.04.23") == 0


# ---------------------------------------------------------------------------
# load_settings / save_settings
# ---------------------------------------------------------------------------

class TestSettings:
    def test_round_trip(self, tmp_path, monkeypatch):
        settings_file = tmp_path / "plugin_installer.json"
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(settings_file))

        PI.save_settings({"check_at_startup": True, "foo": "bar"})
        loaded = PI.load_settings()
        assert loaded["check_at_startup"] is True
        assert loaded["foo"] == "bar"

    def test_missing_file_returns_empty(self, tmp_path, monkeypatch):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "nonexistent.json"))
        result = PI.load_settings()
        assert result == {}

    def test_corrupt_file_returns_empty(self, tmp_path, monkeypatch):
        f = tmp_path / "plugin_installer.json"
        f.write_text("not valid json")
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(f))
        result = PI.load_settings()
        assert result == {}


# ---------------------------------------------------------------------------
# initialize — startup check scheduling
# ---------------------------------------------------------------------------

class TestInitialize:
    def _make_context(self, mw=None):
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw or MagicMock()
        return ctx

    def setup_method(self):
        # Reset the per-session flag before each test
        PI._startup_check_performed = False

    def test_no_main_window_does_not_crash(self):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        PI.initialize(ctx)  # must not raise

    def test_check_at_startup_true_schedules_perform_check(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        PI.save_settings({"check_at_startup": True})

        timer_mock = MagicMock()
        qtcore = types.ModuleType("PyQt6.QtCore")
        qtcore.QTimer = MagicMock()
        qtcore.QTimer.singleShot = timer_mock
        sys.modules["PyQt6.QtCore"] = qtcore

        ctx = self._make_context()
        PI.initialize(ctx)

        timer_mock.assert_called_once()
        delay, fn = timer_mock.call_args[0]
        assert delay == 2000

    def test_first_run_no_key_schedules_ask_permission(self, monkeypatch, tmp_path):
        """When check_at_startup key is absent (first run), ask_user_permission is scheduled."""
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        # No settings file → empty dict → key absent
        (tmp_path / "s.json").unlink(missing_ok=True)

        timer_mock = MagicMock()
        qtcore = types.ModuleType("PyQt6.QtCore")
        qtcore.QTimer = MagicMock()
        qtcore.QTimer.singleShot = timer_mock
        sys.modules["PyQt6.QtCore"] = qtcore

        ctx = self._make_context()
        PI.initialize(ctx)

        timer_mock.assert_called_once()
        delay, fn = timer_mock.call_args[0]
        assert delay == 3000  # first-run uses 3 s delay

    def test_check_at_startup_false_schedules_nothing(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        PI.save_settings({"check_at_startup": False})

        timer_mock = MagicMock()
        qtcore = types.ModuleType("PyQt6.QtCore")
        qtcore.QTimer = MagicMock()
        qtcore.QTimer.singleShot = timer_mock
        sys.modules["PyQt6.QtCore"] = qtcore

        ctx = self._make_context()
        PI.initialize(ctx)

        timer_mock.assert_not_called()

    def test_startup_check_runs_only_once(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        PI.save_settings({"check_at_startup": True})

        timer_mock = MagicMock()
        qtcore = types.ModuleType("PyQt6.QtCore")
        qtcore.QTimer = MagicMock()
        qtcore.QTimer.singleShot = timer_mock
        sys.modules["PyQt6.QtCore"] = qtcore

        ctx = self._make_context()
        PI.initialize(ctx)
        PI.initialize(ctx)  # second call — should do nothing

        assert timer_mock.call_count == 1


# ---------------------------------------------------------------------------
# Version metadata check
# ---------------------------------------------------------------------------

def test_plugin_version_constant_present():
    assert hasattr(PI, "PLUGIN_VERSION")
    assert PI.PLUGIN_VERSION  # non-empty


def test_plugin_version_matches_registry():
    """PLUGIN_VERSION in source must match the registry entry."""
    import re
    registry = json.loads(
        (Path(__file__).resolve().parents[1] / "REGISTRY" / "plugins.json")
        .read_text(encoding="utf-8-sig")
    )
    entry = next((p for p in registry if p.get("id") == "plugin_installer"), None)
    assert entry is not None, "plugin_installer not found in registry"
    assert PI.PLUGIN_VERSION == entry["version"], (
        f"Source version {PI.PLUGIN_VERSION!r} != registry {entry['version']!r}"
    )
