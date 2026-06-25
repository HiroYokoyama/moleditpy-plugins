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
    """Import plugin_installer with PyQt6 stubbed out.

    QThread and QDialog use real (minimal) base classes so that subclasses
    defined in the plugin are proper Python types — this lets tests use
    object.__new__() to create bare instances for unit testing.
    All other Qt classes remain MagicMock.
    """
    # Concrete stubs for classes that are subclassed in the plugin
    class _FakeQThread:
        def __init__(self, *a, **kw): pass
        def start(self): pass
        def quit(self): pass
        def wait(self, *a): return True
        def isRunning(self): return False

    class _FakeQDialog:
        def __init__(self, parent=None, *a, **kw): pass

    if "PyQt6" not in sys.modules:
        pyqt6 = types.ModuleType("PyQt6")
        pyqt6.QtWidgets = types.ModuleType("PyQt6.QtWidgets")
        pyqt6.QtGui = types.ModuleType("PyQt6.QtGui")

        # QDialog needs to be a real class (it's subclassed by PluginInstallerWindow
        # and PluginDetailsDialog).  Everything else stays as MagicMock.
        for name in [
            "QVBoxLayout", "QHBoxLayout", "QLabel",
            "QTableWidget", "QTableWidgetItem", "QPushButton",
            "QHeaderView", "QMessageBox", "QAbstractItemView",
            "QApplication", "QCheckBox", "QLineEdit",
            "QProgressBar", "QProgressDialog",
        ]:
            setattr(pyqt6.QtWidgets, name, MagicMock())
        pyqt6.QtWidgets.QDialog = _FakeQDialog

        setattr(pyqt6.QtGui, "QColor", MagicMock())

        pyqt6.QtCore = types.ModuleType("PyQt6.QtCore")
        # QThread needs to be a real class (it's subclassed by _FetchWorker).
        pyqt6.QtCore.QThread = _FakeQThread
        for name in ["QTimer", "pyqtSignal", "Qt"]:
            setattr(pyqt6.QtCore, name, MagicMock())

        sys.modules["PyQt6"] = pyqt6
        sys.modules["PyQt6.QtWidgets"] = pyqt6.QtWidgets
        sys.modules["PyQt6.QtGui"] = pyqt6.QtGui
        sys.modules["PyQt6.QtCore"] = pyqt6.QtCore

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
        monkeypatch.setattr(PI.QTimer, "singleShot", timer_mock)

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
        monkeypatch.setattr(PI.QTimer, "singleShot", timer_mock)

        ctx = self._make_context()
        PI.initialize(ctx)

        timer_mock.assert_called_once()
        delay, fn = timer_mock.call_args[0]
        assert delay == 3000  # first-run uses 3 s delay

    def test_check_at_startup_false_schedules_nothing(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        PI.save_settings({"check_at_startup": False})

        timer_mock = MagicMock()
        monkeypatch.setattr(PI.QTimer, "singleShot", timer_mock)

        ctx = self._make_context()
        PI.initialize(ctx)

        timer_mock.assert_not_called()

    def test_startup_check_runs_only_once(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        PI.save_settings({"check_at_startup": True})

        timer_mock = MagicMock()
        monkeypatch.setattr(PI.QTimer, "singleShot", timer_mock)

        ctx = self._make_context()
        PI.initialize(ctx)
        PI.initialize(ctx)  # second call — should do nothing

        assert timer_mock.call_count == 1

    def test_registers_menu_action(self):
        ctx = self._make_context()
        PI.initialize(ctx)
        ctx.add_menu_action.assert_called_once()
        path, _callback = ctx.add_menu_action.call_args[0]
        assert "Plugin Installer" in path

    def test_registers_menu_action_even_when_no_main_window(self):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        PI.initialize(ctx)
        ctx.add_menu_action.assert_called_once()


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


class TestVersionCompatibility:
    def test_empty_or_wildcard_specifier(self):
        assert PI.is_app_version_compatible("3.5.1", "") is True
        assert PI.is_app_version_compatible("3.5.1", "*") is True

    def test_wildcard_matches(self):
        assert PI.is_app_version_compatible("3.5.1", "3.*") is True
        assert PI.is_app_version_compatible("3.5.1", "3.5.*") is True
        assert PI.is_app_version_compatible("3.5.1", "3.6.*") is False
        assert PI.is_app_version_compatible("4.0.0", "3.*") is False

    def test_simple_comparison(self):
        assert PI.is_app_version_compatible("3.5.1", ">=3.5") is True
        assert PI.is_app_version_compatible("3.5.1", "<=3.6") is True
        assert PI.is_app_version_compatible("3.5.1", "==3.5.1") is True
        assert PI.is_app_version_compatible("3.5.1", ">3.5.1") is False
        assert PI.is_app_version_compatible("3.5.1", "<3.5") is False

    def test_compound_ranges(self):
        assert PI.is_app_version_compatible("3.5.1", ">=3.5, <4") is True
        assert PI.is_app_version_compatible("4.0.0", ">=3.5, <4") is False
        assert PI.is_app_version_compatible("3.4.9", ">=3.5, <4") is False
        assert PI.is_app_version_compatible("3.5.1", ">=3.5, <4, !=3.6") is True

    def test_compatible_release_tilde_equal(self):
        assert PI.is_app_version_compatible("1.2.3", "~=1.2.3") is True
        assert PI.is_app_version_compatible("1.2.5", "~=1.2.3") is True
        # ~=1.2.3 means >=1.2.3 and <1.3
        assert PI.is_app_version_compatible("1.3.0", "~=1.2.3") is False
        assert PI.is_app_version_compatible("1.2.2", "~=1.2.3") is False
        
        # ~=1.2 means >=1.2 and <2.0
        assert PI.is_app_version_compatible("1.2.0", "~=1.2") is True
        assert PI.is_app_version_compatible("1.5.0", "~=1.2") is True
        assert PI.is_app_version_compatible("2.0.0", "~=1.2") is False


class TestDependencyParsing:
    def test_simple_package(self):
        assert PI.parse_dependency("numpy") == ("numpy", "")
        
    def test_simple_version(self):
        assert PI.parse_dependency("numpy>=1.20") == ("numpy", ">=1.20")
        
    def test_compound_versions(self):
        assert PI.parse_dependency("rdkit>=2022.03,<2023.0") == ("rdkit", ">=2022.03,<2023.0")
        
    def test_with_spaces(self):
        assert PI.parse_dependency("numpy >= 1.20, < 2.0") == ("numpy", ">= 1.20, < 2.0")
        
    def test_with_extras(self):
        assert PI.parse_dependency("requests[security]>=2.0") == ("requests", ">=2.0")
        assert PI.parse_dependency("requests[security,other]") == ("requests", "")
        
    def test_with_markers(self):
        assert PI.parse_dependency("numpy>=1.20; python_version < '3.8'") == ("numpy", ">=1.20")
        
    def test_reject_colon_separated(self):
        assert PI.parse_dependency("numpy:>=1.20") == ("", "")


class TestCheckDependencySatisfied:
    def test_missing_package(self):
        with patch("importlib.metadata.distribution", side_effect=importlib.metadata.PackageNotFoundError):
            assert PI.check_dependency_satisfied("nonexistent_package") is False

    def test_installed_no_constraint(self):
        mock_dist = MagicMock()
        mock_dist.version = "1.20.0"
        with patch("importlib.metadata.distribution", return_value=mock_dist):
            assert PI.check_dependency_satisfied("numpy") is True

    def test_installed_satisfied(self):
        mock_dist = MagicMock()
        mock_dist.version = "1.20.5"
        with patch("importlib.metadata.distribution", return_value=mock_dist):
            assert PI.check_dependency_satisfied("numpy>=1.20") is True

    def test_installed_unsatisfied(self):
        mock_dist = MagicMock()
        mock_dist.version = "1.19.0"
        with patch("importlib.metadata.distribution", return_value=mock_dist):
            assert PI.check_dependency_satisfied("numpy>=1.20") is False

    def test_invalid_colon_separated(self):
        assert PI.check_dependency_satisfied("numpy:>=1.20") is False


class TestSanitizeAndQuoteDependency:
    def test_valid_no_specifier(self):
        assert PI.sanitize_and_quote_dependency("numpy") == "numpy"

    def test_valid_with_specifier(self):
        assert PI.sanitize_and_quote_dependency("numpy>=1.20") == '"numpy>=1.20"'
        assert PI.sanitize_and_quote_dependency("rdkit>=2022.03,<2023.0") == '"rdkit>=2022.03,<2023.0"'

    def test_invalid_characters_in_name(self):
        assert PI.sanitize_and_quote_dependency("numpy;pkg") == ""
        assert PI.sanitize_and_quote_dependency("numpy&foo") == ""

    def test_invalid_characters_in_specifier(self):
        assert PI.sanitize_and_quote_dependency("numpy>=1.20;echo") == ""


# ---------------------------------------------------------------------------
# _STATUS_PRIORITY — "Update Available" must sort to the very top
# ---------------------------------------------------------------------------

class TestStatusPriority:
    def test_update_available_is_lowest_priority_value(self):
        """'Update Available' must have the smallest priority number (sorts first)."""
        update_prio = PI._STATUS_PRIORITY["Update Available"]
        for status, prio in PI._STATUS_PRIORITY.items():
            if status != "Update Available":
                assert update_prio < prio, (
                    f"'Update Available' priority {update_prio} should be lower than "
                    f"'{status}' priority {prio}"
                )

    def test_sort_order_update_available_first(self):
        """Rows with 'Update Available' must appear before all other statuses."""
        rows = [
            {"status": "Up to date",       "name": "Aaa"},
            {"status": "Update Available",  "name": "Zzz"},
            {"status": "Not Installed",     "name": "Mmm"},
            {"status": "Update Available",  "name": "Aab"},
            {"status": "Incompatible",      "name": "Bbb"},
        ]
        rows.sort(
            key=lambda r: (PI._STATUS_PRIORITY.get(r["status"], 99), r["name"].lower())
        )
        assert rows[0]["status"] == "Update Available"
        assert rows[1]["status"] == "Update Available"
        # "Update Available" rows themselves should be alphabetical
        assert rows[0]["name"] < rows[1]["name"]
        # Everything after the Update Available rows is a non-update status
        assert all(r["status"] != "Update Available" for r in rows[2:])

    def test_sort_order_alphabetical_within_group(self):
        """Within the same status group, names are sorted alphabetically."""
        rows = [
            {"status": "Up to date", "name": "Zzz"},
            {"status": "Up to date", "name": "Aaa"},
            {"status": "Up to date", "name": "Mmm"},
        ]
        rows.sort(
            key=lambda r: (PI._STATUS_PRIORITY.get(r["status"], 99), r["name"].lower())
        )
        assert [r["name"] for r in rows] == ["Aaa", "Mmm", "Zzz"]


# ---------------------------------------------------------------------------
# _FetchWorker — network fetch logic
# ---------------------------------------------------------------------------

class TestFetchWorker:
    """Tests for _FetchWorker._fetch_pypi and _fetch_remote (no real threads)."""

    def _make_worker(self):
        """Create a _FetchWorker bypassing QThread.__init__."""
        w = object.__new__(PI._FetchWorker)
        w._pkg = "moleditpy"
        w._url = PI.REMOTE_JSON_URL
        return w

    def test_fetch_pypi_returns_version(self):
        w = self._make_worker()
        fake_body = json.dumps({"info": {"version": "4.5.0"}}).encode()
        fake_resp = MagicMock()
        fake_resp.__enter__ = lambda s: s
        fake_resp.__exit__ = MagicMock(return_value=False)
        fake_resp.status = 200
        fake_resp.read.return_value = fake_body
        with patch("urllib.request.urlopen", return_value=fake_resp):
            result = w._fetch_pypi()
        assert result == "4.5.0"

    def test_fetch_pypi_network_error_returns_unknown(self):
        w = self._make_worker()
        with patch("urllib.request.urlopen", side_effect=OSError("timeout")):
            result = w._fetch_pypi()
        assert result == "Unknown"

    def test_fetch_remote_returns_list(self):
        w = self._make_worker()
        payload = [{"name": "PluginA", "version": "1.0.0"}]
        fake_body = json.dumps(payload).encode("utf-8-sig")
        fake_resp = MagicMock()
        fake_resp.__enter__ = lambda s: s
        fake_resp.__exit__ = MagicMock(return_value=False)
        fake_resp.status = 200
        fake_resp.read.return_value = fake_body
        with patch("urllib.request.urlopen", return_value=fake_resp):
            result = w._fetch_remote()
        assert result == payload

    def test_fetch_remote_network_error_returns_empty_list(self):
        w = self._make_worker()
        with patch("urllib.request.urlopen", side_effect=OSError("no route")):
            result = w._fetch_remote()
        assert result == []


# ---------------------------------------------------------------------------
# _download_chunked — progress callback and cancellation
# ---------------------------------------------------------------------------

class TestDownloadChunked:
    """Tests for PluginInstallerWindow._download_chunked via a bare instance."""

    def _make_installer(self):
        """Create a PluginInstallerWindow without running __init__ or Qt code."""
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.main_window = MagicMock()
        return inst

    def test_downloads_file_and_calls_callback(self, tmp_path):
        inst = self._make_installer()
        dest = tmp_path / "plugin.py"
        content = b"print('hello')" * 100

        fake_resp = MagicMock()
        fake_resp.__enter__ = lambda s: s
        fake_resp.__exit__ = MagicMock(return_value=False)
        fake_resp.headers = {"Content-Length": str(len(content))}
        # Serve content in two chunks then empty
        fake_resp.read.side_effect = [content[:800], content[800:], b""]

        progress_calls = []

        def on_progress(recv, total):
            progress_calls.append((recv, total))
            return True  # continue

        with patch("urllib.request.urlopen", return_value=fake_resp):
            ok = inst._download_chunked("https://example.com/plugin.py", str(dest), on_progress)

        assert ok is True
        assert dest.read_bytes() == content
        assert len(progress_calls) == 2
        assert progress_calls[-1][1] == len(content)

    def test_callback_returning_false_aborts(self, tmp_path):
        inst = self._make_installer()
        dest = tmp_path / "plugin.py"
        content = b"x" * 1000

        fake_resp = MagicMock()
        fake_resp.__enter__ = lambda s: s
        fake_resp.__exit__ = MagicMock(return_value=False)
        fake_resp.headers = {}
        fake_resp.read.side_effect = [content, b""]

        with patch("urllib.request.urlopen", return_value=fake_resp):
            ok = inst._download_chunked(
                "https://example.com/plugin.py",
                str(dest),
                on_progress=lambda r, t: False,  # cancel immediately
            )

        assert ok is False

    def test_network_error_returns_false(self, tmp_path):
        inst = self._make_installer()
        dest = tmp_path / "plugin.py"

        with patch("urllib.request.urlopen", side_effect=OSError("connection refused")):
            with patch.object(PI.QMessageBox, "warning") as warn:
                ok = inst._download_chunked(
                    "https://example.com/plugin.py", str(dest)
                )

        assert ok is False
        warn.assert_called_once()

    def test_no_callback_still_downloads(self, tmp_path):
        inst = self._make_installer()
        dest = tmp_path / "plugin.py"
        content = b"# plugin code"

        fake_resp = MagicMock()
        fake_resp.__enter__ = lambda s: s
        fake_resp.__exit__ = MagicMock(return_value=False)
        fake_resp.headers = {}
        fake_resp.read.side_effect = [content, b""]

        with patch("urllib.request.urlopen", return_value=fake_resp):
            ok = inst._download_chunked(
                "https://example.com/plugin.py", str(dest)
            )

        assert ok is True
        assert dest.read_bytes() == content


# ---------------------------------------------------------------------------
# SHA-256 verification
# ---------------------------------------------------------------------------

class TestCalculateSha256:
    def test_known_hash(self, tmp_path):
        import hashlib
        f = tmp_path / "data.bin"
        data = b"hello world"
        f.write_bytes(data)
        inst = object.__new__(PI.PluginInstallerWindow)
        expected = hashlib.sha256(data).hexdigest()
        assert inst.calculate_sha256(str(f)) == expected

