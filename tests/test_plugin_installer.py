"""
Unit tests for the Plugin Installer plugin.

All tests run headlessly — PySide6 and network are fully mocked.
"""

import json
import os
import sys
import types
import importlib
import pytest
from pathlib import Path
import re as _re
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# Helpers — load the module without executing Qt or network code at import
# ---------------------------------------------------------------------------

PLUGIN_PATH = (
    Path(__file__).resolve().parents[1]
    / "plugins"
    / "Plugin_Installer"
    / "plugin_installer.py"
)


def _load_module():
    """Import plugin_installer with PySide6 stubbed out.

    QThread and QDialog use real (minimal) base classes so that subclasses
    defined in the plugin are proper Python types — this lets tests use
    object.__new__() to create bare instances for unit testing.
    All other Qt classes remain MagicMock.
    """

    # Concrete stubs for classes that are subclassed in the plugin
    class _FakeQThread:
        def __init__(self, *a, **kw):
            pass

        def start(self):
            pass

        def quit(self):
            pass

        def wait(self, *a):
            return True

        def isRunning(self):
            return False

    class _FakeQDialog:
        def __init__(self, parent=None, *a, **kw):
            pass

    # Always install the stub PySide6 while loading the plugin, even if pytest-qt
    # has already imported the real PySide6 during collection. Otherwise the plugin
    # subclasses the real QDialog and object.__new__() on a sip type is unsafe
    # (this is why the suite passed in CI, where PySide6 is absent, but failed
    # locally with PySide6 installed). The real modules are restored afterwards.
    _saved_pyside6 = {
        key: sys.modules.get(key)
        for key in ("PySide6", "PySide6.QtWidgets", "PySide6.QtGui", "PySide6.QtCore")
    }
    if True:  # force the stub regardless of any real PySide6 already imported
        pyside6 = types.ModuleType("PySide6")
        pyside6.QtWidgets = types.ModuleType("PySide6.QtWidgets")
        pyside6.QtGui = types.ModuleType("PySide6.QtGui")

        # QDialog needs to be a real class (it's subclassed by PluginInstallerWindow
        # and PluginDetailsDialog).  Everything else stays as MagicMock.
        for name in [
            "QVBoxLayout",
            "QHBoxLayout",
            "QLabel",
            "QTableWidget",
            "QTableWidgetItem",
            "QPushButton",
            "QHeaderView",
            "QMessageBox",
            "QAbstractItemView",
            "QApplication",
            "QCheckBox",
            "QLineEdit",
            "QProgressBar",
            "QProgressDialog",
        ]:
            setattr(pyside6.QtWidgets, name, MagicMock())
        pyside6.QtWidgets.QDialog = _FakeQDialog

        for name in ["QColor", "QDesktopServices"]:
            setattr(pyside6.QtGui, name, MagicMock())

        pyside6.QtCore = types.ModuleType("PySide6.QtCore")
        # QThread needs to be a real class (it's subclassed by _FetchWorker).
        pyside6.QtCore.QThread = _FakeQThread
        for name in ["QTimer", "Signal", "Qt", "QUrl"]:
            setattr(pyside6.QtCore, name, MagicMock())

        sys.modules["PySide6"] = pyside6
        sys.modules["PySide6.QtWidgets"] = pyside6.QtWidgets
        sys.modules["PySide6.QtGui"] = pyside6.QtGui
        sys.modules["PySide6.QtCore"] = pyside6.QtCore

    spec = importlib.util.spec_from_file_location("plugin_installer", PLUGIN_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PI = _load_module()


# ---------------------------------------------------------------------------
# compare_versions — mirrors PluginInstallerWindow.compare_versions exactly
# (tested standalone because the class inherits from a mocked QDialog)
# ---------------------------------------------------------------------------


def _compare_versions(v1, v2):
    """Mirrors PluginInstallerWindow.compare_versions for unit testing."""
    try:

        def parse(v):
            return [int(x) for x in _re.findall(r"\d+", v)]

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
        PI._startup_check_performed = False
        app_mock = PI.QApplication.instance.return_value
        app_mock.reset_mock()
        app_mock.property.return_value = False

    def test_no_main_window_does_not_crash(self):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        PI.initialize(ctx)  # must not raise

    def test_check_at_startup_true_schedules_perform_check(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        PI.save_settings({"check_at_startup": True})

        timer_mock = MagicMock()
        monkeypatch.setattr(PI, "QTimer", MagicMock(singleShot=timer_mock))

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
        monkeypatch.setattr(PI, "QTimer", MagicMock(singleShot=timer_mock))

        ctx = self._make_context()
        PI.initialize(ctx)

        timer_mock.assert_called_once()
        delay, fn = timer_mock.call_args[0]
        assert delay == 3000  # first-run uses 3 s delay

    def test_check_at_startup_false_schedules_nothing(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        PI.save_settings({"check_at_startup": False})

        timer_mock = MagicMock()
        monkeypatch.setattr(PI, "QTimer", MagicMock(singleShot=timer_mock))

        ctx = self._make_context()
        PI.initialize(ctx)

        timer_mock.assert_not_called()

    def test_startup_check_runs_only_once(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        PI.save_settings({"check_at_startup": True})

        timer_mock = MagicMock()
        monkeypatch.setattr(PI, "QTimer", MagicMock(singleShot=timer_mock))

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

    def test_asks_for_the_header_slot(self):
        """MoleditPy >= 4.8.1 is told where this entry belongs."""
        ctx = self._make_context()
        PI.initialize(ctx)
        assert ctx.add_menu_action.call_args.kwargs.get("pin") == "header"

    def test_falls_back_when_pin_is_unsupported(self):
        """An older host has no `pin` parameter and must not lose the plugin.

        Letting the TypeError escape initialize() would drop the installer
        entirely; those releases place it by name instead.
        """
        ctx = self._make_context()
        ctx.add_menu_action.side_effect = [TypeError("unexpected keyword 'pin'"), None]

        PI.initialize(ctx)

        assert ctx.add_menu_action.call_count == 2
        retry_args, retry_kwargs = ctx.add_menu_action.call_args
        assert "pin" not in retry_kwargs
        assert "Plugin Installer" in retry_args[0]

    def test_registers_menu_action_even_when_no_main_window(self):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        PI.initialize(ctx)
        ctx.add_menu_action.assert_called_once()


# ---------------------------------------------------------------------------
# Startup-check guard — module-reload resilience
# ---------------------------------------------------------------------------


class TestStartupCheckGuard:
    """Verifies the module-reload-safe startup check guard."""

    def setup_method(self):
        PI._startup_check_performed = False
        app_mock = PI.QApplication.instance.return_value
        app_mock.reset_mock()
        app_mock.property.return_value = False

    def test_not_performed_initially(self):
        assert PI._is_startup_check_performed() is False

    def test_mark_sets_module_level_flag(self):
        PI._mark_startup_check_performed()
        assert PI._startup_check_performed is True

    def test_mark_sets_qapplication_property(self):
        PI._mark_startup_check_performed()
        app_mock = PI.QApplication.instance.return_value
        app_mock.setProperty.assert_called_once_with(PI._STARTUP_CHECK_PROP, True)

    def test_is_performed_reads_qapplication_property(self):
        # Module-level flag is False (simulates post-reload state); QApp prop is True.
        PI._startup_check_performed = False
        app_mock = PI.QApplication.instance.return_value
        app_mock.property.return_value = True
        assert PI._is_startup_check_performed() is True

    def test_module_reload_does_not_retrigger_startup_check(
        self, monkeypatch, tmp_path
    ):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        PI.save_settings({"check_at_startup": True})
        timer_mock = MagicMock()
        monkeypatch.setattr(PI, "QTimer", MagicMock(singleShot=timer_mock))

        app_mock = PI.QApplication.instance.return_value
        app_mock.property.return_value = False
        ctx = MagicMock()
        ctx.get_main_window.return_value = MagicMock()
        PI.initialize(ctx)
        assert timer_mock.call_count == 1

        # Simulate reload: module-level flag reset, QApp property survives.
        PI._startup_check_performed = False
        app_mock.property.return_value = True
        PI.initialize(ctx)
        assert timer_mock.call_count == 1

    def test_no_qapplication_falls_back_to_module_level(self, monkeypatch):
        monkeypatch.setattr(PI.QApplication, "instance", lambda: None)
        PI._startup_check_performed = False
        assert PI._is_startup_check_performed() is False
        PI._mark_startup_check_performed()
        assert PI._startup_check_performed is True
        assert PI._is_startup_check_performed() is True


# ---------------------------------------------------------------------------
# ask_user_permission — pre-dialog settings re-check
# ---------------------------------------------------------------------------


class TestAskUserPermission:
    """ask_user_permission() re-reads settings before showing the dialog."""

    def setup_method(self):
        # Reset the shared QMessageBox mock so call counts are fresh per test.
        PI.QMessageBox.question.reset_mock()

    def _setup(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI, "SETTINGS_FILE", str(tmp_path / "s.json"))
        return MagicMock()

    def test_shows_dialog_when_key_absent(self, monkeypatch, tmp_path):
        mw = self._setup(monkeypatch, tmp_path)
        PI.QMessageBox.question.return_value = PI.QMessageBox.StandardButton.No
        PI.ask_user_permission(mw)
        PI.QMessageBox.question.assert_called_once()

    def test_skips_dialog_when_key_already_false(self, monkeypatch, tmp_path):
        mw = self._setup(monkeypatch, tmp_path)
        PI.save_settings({"check_at_startup": False})
        perform_mock = MagicMock()
        monkeypatch.setattr(PI, "perform_startup_check", perform_mock)
        PI.ask_user_permission(mw)
        PI.QMessageBox.question.assert_not_called()
        perform_mock.assert_not_called()

    def test_skips_dialog_and_runs_check_when_key_already_true(
        self, monkeypatch, tmp_path
    ):
        mw = self._setup(monkeypatch, tmp_path)
        PI.save_settings({"check_at_startup": True})
        perform_mock = MagicMock()
        monkeypatch.setattr(PI, "perform_startup_check", perform_mock)
        PI.ask_user_permission(mw)
        PI.QMessageBox.question.assert_not_called()
        perform_mock.assert_called_once_with(mw)

    def test_yes_saves_true_and_runs_check(self, monkeypatch, tmp_path):
        mw = self._setup(monkeypatch, tmp_path)
        PI.QMessageBox.question.return_value = PI.QMessageBox.StandardButton.Yes
        perform_mock = MagicMock()
        monkeypatch.setattr(PI, "perform_startup_check", perform_mock)
        PI.ask_user_permission(mw)
        assert PI.load_settings().get("check_at_startup") is True
        perform_mock.assert_called_once_with(mw)

    def test_no_saves_false_and_does_not_run_check(self, monkeypatch, tmp_path):
        mw = self._setup(monkeypatch, tmp_path)
        PI.QMessageBox.question.return_value = PI.QMessageBox.StandardButton.No
        perform_mock = MagicMock()
        monkeypatch.setattr(PI, "perform_startup_check", perform_mock)
        PI.ask_user_permission(mw)
        assert PI.load_settings().get("check_at_startup") is False
        perform_mock.assert_not_called()


# ---------------------------------------------------------------------------
# Version metadata check
# ---------------------------------------------------------------------------


def test_plugin_version_constant_present():
    assert hasattr(PI, "PLUGIN_VERSION")
    assert PI.PLUGIN_VERSION  # non-empty


def test_plugin_version_matches_registry():
    """PLUGIN_VERSION in source must match the registry entry."""
    registry = json.loads(
        (Path(__file__).resolve().parents[1] / "REGISTRY" / "plugins.json").read_text(
            encoding="utf-8-sig"
        )
    )
    entry = next((p for p in registry if p.get("id") == "plugin_installer"), None)
    assert entry is not None, "plugin_installer not found in registry"
    assert PI.PLUGIN_VERSION == entry["version"], (
        f"Source version {PI.PLUGIN_VERSION!r} != registry {entry['version']!r}"
    )


class TestOSDetection:
    def test_windows(self, monkeypatch):
        monkeypatch.setattr(PI.sys, "platform", "win32")
        assert PI.get_current_os() == "Windows"

    def test_macos(self, monkeypatch):
        monkeypatch.setattr(PI.sys, "platform", "darwin")
        assert PI.get_current_os() == "macOS"

    def test_native_linux(self, monkeypatch):
        monkeypatch.setattr(PI.sys, "platform", "linux")
        monkeypatch.setattr(PI, "_running_under_wsl", lambda: False)
        assert PI.get_current_os() == "Linux"

    def test_wsl_detected_via_env(self, monkeypatch):
        monkeypatch.setattr(PI.sys, "platform", "linux")
        monkeypatch.setenv("WSL_DISTRO_NAME", "Ubuntu")
        assert PI.get_current_os() == "WSL"

    def test_wsl_detected_via_proc_version(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI.sys, "platform", "linux")
        monkeypatch.delenv("WSL_DISTRO_NAME", raising=False)
        proc = tmp_path / "version"
        proc.write_text("Linux version 5.15.0-microsoft-standard-WSL2")
        real_open = open
        monkeypatch.setattr(
            PI, "open", lambda p, *a, **kw: real_open(proc, *a, **kw), raising=False
        )
        assert PI.get_current_os() == "WSL"

    def test_proc_version_without_microsoft_is_plain_linux(self, monkeypatch, tmp_path):
        monkeypatch.setattr(PI.sys, "platform", "linux")
        monkeypatch.delenv("WSL_DISTRO_NAME", raising=False)
        proc = tmp_path / "version"
        proc.write_text("Linux version 6.1.0-generic")
        real_open = open
        monkeypatch.setattr(
            PI, "open", lambda p, *a, **kw: real_open(proc, *a, **kw), raising=False
        )
        assert PI.get_current_os() == "Linux"

    def test_unknown_platform(self, monkeypatch):
        assert PI.get_current_os.__module__ is not None
        monkeypatch.setattr(PI.sys, "platform", "sunos5")
        assert PI.get_current_os() == "Unknown"


class TestOSCompatibility:
    def test_unspecified_is_compatible(self):
        # A missing/empty list means "unstated", never a warning.
        assert PI.is_os_compatible("Windows", None) is True
        assert PI.is_os_compatible("Windows", []) is True
        assert PI.is_os_compatible("Windows", "not-a-list") is True

    def test_listed_os_is_compatible(self):
        assert PI.is_os_compatible("Windows", ["Windows", "macOS"]) is True
        assert PI.is_os_compatible("macOS", ["macOS", "Linux", "WSL"]) is True

    def test_unlisted_os_is_incompatible(self):
        assert PI.is_os_compatible("Windows", ["macOS", "Linux", "WSL"]) is False

    def test_wsl_satisfies_plain_linux(self):
        # WSL is a Linux userland, so a "Linux"-only plugin runs there.
        assert PI.is_os_compatible("WSL", ["macOS", "Linux"]) is True

    def test_wsl_satisfies_explicit_wsl(self):
        assert PI.is_os_compatible("WSL", ["macOS", "Linux", "WSL"]) is True

    def test_linux_is_not_satisfied_by_wsl_only_entry(self):
        assert PI.is_os_compatible("Linux", ["Windows", "WSL"]) is False

    def test_matching_is_case_insensitive(self):
        assert PI.is_os_compatible("macOS", ["macos"]) is True
        assert PI.is_os_compatible("Windows", ["WINDOWS"]) is True

    def test_unknown_current_os_never_warns(self):
        assert PI.is_os_compatible("Unknown", ["Windows"]) is True


class TestFormatSupportedOS:
    def test_list(self):
        assert PI.format_supported_os(["macOS", "Linux"]) == "macOS, Linux"

    def test_string_passthrough(self):
        assert PI.format_supported_os("Linux") == "Linux"

    def test_empty_is_unknown(self):
        assert PI.format_supported_os(None) == "Unknown"
        assert PI.format_supported_os([]) == "Unknown"
        assert PI.format_supported_os("  ") == "Unknown"


class TestOSIncompatibilityHint:
    def test_windows_user_told_about_wsl(self):
        hint = PI.os_incompatibility_hint("Windows", ["macOS", "Linux", "WSL"])
        assert "WSL" in hint

    def test_windows_user_told_about_wsl_when_only_linux_listed(self):
        hint = PI.os_incompatibility_hint("Windows", ["macOS", "Linux"])
        assert "WSL" in hint

    def test_no_hint_for_non_windows(self):
        assert PI.os_incompatibility_hint("macOS", ["Windows"]) == ""

    def test_no_hint_when_wsl_route_absent(self):
        assert PI.os_incompatibility_hint("Windows", ["macOS"]) == ""


class TestVersionCompatibility:
    @pytest.mark.parametrize(
        "host, expected",
        [("5.0.0a1", True), ("5.0.0", True), ("5.2.0", True),
         ("4.13.0", False), ("6.0.0a1", False)],
    )
    def test_stable_host_range_also_accepts_its_alpha(self, host, expected):
        assert PI.is_app_version_compatible(host, ">=5.0.0, <6.0.0") is expected

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
        assert PI.parse_dependency("rdkit>=2022.03,<2023.0") == (
            "rdkit",
            ">=2022.03,<2023.0",
        )

    def test_with_spaces(self):
        assert PI.parse_dependency("numpy >= 1.20, < 2.0") == (
            "numpy",
            ">= 1.20, < 2.0",
        )

    def test_with_extras(self):
        assert PI.parse_dependency("requests[security]>=2.0") == ("requests", ">=2.0")
        assert PI.parse_dependency("requests[security,other]") == ("requests", "")

    def test_with_markers(self):
        assert PI.parse_dependency("numpy>=1.20; python_version < '3.8'") == (
            "numpy",
            ">=1.20",
        )

    def test_reject_colon_separated(self):
        assert PI.parse_dependency("numpy:>=1.20") == ("", "")


class TestCheckDependencySatisfied:
    def test_missing_package(self):
        with patch(
            "importlib.metadata.distribution",
            side_effect=importlib.metadata.PackageNotFoundError,
        ):
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
        assert (
            PI.sanitize_and_quote_dependency("rdkit>=2022.03,<2023.0")
            == '"rdkit>=2022.03,<2023.0"'
        )

    def test_invalid_characters_in_name(self):
        assert PI.sanitize_and_quote_dependency("numpy;pkg") == ""
        assert PI.sanitize_and_quote_dependency("numpy&foo") == ""

    def test_invalid_characters_in_specifier(self):
        assert PI.sanitize_and_quote_dependency("numpy>=1.20;echo") == ""


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
            ok = inst._download_chunked(
                "https://example.com/plugin.py", str(dest), on_progress
            )

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

    def test_network_error_retries_3_times_then_shows_dialog(self, tmp_path):
        """URLError retried up to 3 times; dialog shown only after all attempts fail."""
        import urllib.error

        inst = self._make_installer()
        dest = tmp_path / "plugin.py"

        with patch(
            "urllib.request.urlopen",
            side_effect=urllib.error.URLError("timeout"),
        ) as mock_open:
            with patch.object(PI.QMessageBox, "warning") as warn:
                with patch("time.sleep") as mock_sleep:
                    ok = inst._download_chunked(
                        "https://example.com/plugin.py", str(dest)
                    )

        assert ok is False
        assert mock_open.call_count == 3
        warn.assert_called_once()
        # sleep(2) called between retries: after attempt 1 and 2, not after 3
        assert mock_sleep.call_count == 2
        mock_sleep.assert_called_with(2)

    def test_transient_network_error_recovers_on_second_attempt(self, tmp_path):
        """Fails once then succeeds — returns True without showing a dialog."""
        import urllib.error

        inst = self._make_installer()
        dest = tmp_path / "plugin.py"
        content = b"# plugin"

        fake_resp = MagicMock()
        fake_resp.__enter__ = lambda s: s
        fake_resp.__exit__ = MagicMock(return_value=False)
        fake_resp.headers = {}
        fake_resp.read.side_effect = [content, b""]

        with patch(
            "urllib.request.urlopen",
            side_effect=[urllib.error.URLError("timeout"), fake_resp],
        ) as mock_open:
            with patch.object(PI.QMessageBox, "warning") as warn:
                with patch("time.sleep"):
                    ok = inst._download_chunked(
                        "https://example.com/plugin.py", str(dest)
                    )

        assert ok is True
        assert dest.read_bytes() == content
        assert mock_open.call_count == 2
        warn.assert_not_called()

    def test_non_network_error_does_not_retry(self, tmp_path):
        """Non-retryable exception shows dialog immediately without retrying."""
        inst = self._make_installer()
        dest = tmp_path / "plugin.py"

        with patch(
            "urllib.request.urlopen",
            side_effect=ValueError("unexpected"),
        ) as mock_open:
            with patch.object(PI.QMessageBox, "warning") as warn:
                with patch("time.sleep") as mock_sleep:
                    ok = inst._download_chunked(
                        "https://example.com/plugin.py", str(dest)
                    )

        assert ok is False
        assert mock_open.call_count == 1
        warn.assert_called_once()
        mock_sleep.assert_not_called()

    def test_network_error_returns_false(self, tmp_path):
        inst = self._make_installer()
        dest = tmp_path / "plugin.py"

        with patch("urllib.request.urlopen", side_effect=OSError("connection refused")):
            with patch.object(PI.QMessageBox, "warning") as warn:
                with patch("time.sleep"):
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
            ok = inst._download_chunked("https://example.com/plugin.py", str(dest))

        assert ok is True
        assert dest.read_bytes() == content


# ---------------------------------------------------------------------------
# _read_plugin_version_ast — version extraction without import
# ---------------------------------------------------------------------------


class TestReadPluginVersionAst:
    def test_reads_version(self, tmp_path):
        f = tmp_path / "plugin.py"
        f.write_text('PLUGIN_VERSION = "2026.06.26"\n')
        assert PI._read_plugin_version_ast(str(f)) == "2026.06.26"

    def test_missing_constant_returns_unknown(self, tmp_path):
        f = tmp_path / "plugin.py"
        f.write_text("x = 1\n")
        assert PI._read_plugin_version_ast(str(f)) == "Unknown"

    def test_nonexistent_file_returns_unknown(self, tmp_path):
        assert PI._read_plugin_version_ast(str(tmp_path / "missing.py")) == "Unknown"

    def test_syntax_error_returns_unknown(self, tmp_path):
        f = tmp_path / "plugin.py"
        f.write_text("def (:\n")
        assert PI._read_plugin_version_ast(str(f)) == "Unknown"

    def test_reads_version_from_package_init(self, tmp_path):
        pkg = tmp_path / "My_Plugin"
        pkg.mkdir()
        init = pkg / "__init__.py"
        init.write_text('PLUGIN_VERSION = "2026.01.01"\nPLUGIN_NAME = "My Plugin"\n')
        assert PI._read_plugin_version_ast(str(init)) == "2026.01.01"


# ---------------------------------------------------------------------------
# _update_status_label — update-count label text and colour
# ---------------------------------------------------------------------------


class TestUpdateStatusLabel:
    def _make_installer(self, row_count=5):
        inst = object.__new__(PI.PluginInstallerWindow)
        inst._status_label = MagicMock()
        inst.table = MagicMock()
        inst.table.rowCount.return_value = row_count
        return inst

    def test_singular_update(self):
        inst = self._make_installer()
        inst._update_status_label(1)
        inst._status_label.setText.assert_called_once_with(
            "1 plugin has an update available"
        )

    def test_plural_updates(self):
        inst = self._make_installer()
        inst._update_status_label(3)
        inst._status_label.setText.assert_called_once_with(
            "3 plugins have updates available"
        )

    def test_all_up_to_date(self):
        inst = self._make_installer()
        inst._update_status_label(0)
        inst._status_label.setText.assert_called_once_with("All plugins are up to date")

    def test_no_rows_clears_label(self):
        inst = self._make_installer(row_count=0)
        inst._update_status_label(0)
        inst._status_label.setText.assert_called_once_with("")

    def test_no_label_does_not_crash(self):
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.table = MagicMock()
        inst.table.rowCount.return_value = 5
        inst._update_status_label(2)  # must not raise


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


# ---------------------------------------------------------------------------
# _overwrite_folder_plugin — wipe-and-replace with settings.json preservation
# ---------------------------------------------------------------------------


class TestOverwriteFolderPlugin:
    """Updating a folder plugin must remove orphaned files while keeping settings."""

    def _make_source(self, tmp_path, files):
        src = tmp_path / "source"
        src.mkdir()
        for name, content in files.items():
            (src / name).write_text(content, encoding="utf-8")
        return src

    def test_copy_failure_preserves_installed_plugin(self, tmp_path):
        target = self._make_installed(tmp_path, {
            "__init__.py": "old", "settings.json": "user settings"
        })
        source = self._make_source(tmp_path, {"__init__.py": "new"})
        with patch.object(PI.shutil, "copytree", side_effect=OSError("disk full")):
            import pytest
            with pytest.raises(OSError, match="disk full"):
                PI.PluginInstallerWindow._overwrite_folder_plugin(str(source), str(target))
        assert (target / "__init__.py").read_text() == "old"
        assert (target / "settings.json").read_text() == "user settings"
        assert not list(tmp_path.glob(".plugin-update-*"))

    def test_promotion_failure_rolls_back(self, tmp_path):
        target = self._make_installed(tmp_path, {"__init__.py": "old"})
        source = self._make_source(tmp_path, {"__init__.py": "new"})
        real_replace = PI.os.replace
        def fail_promotion(src, dst):
            if Path(src).name == "replacement":
                raise OSError("promotion failed")
            return real_replace(src, dst)
        with patch.object(PI.os, "replace", side_effect=fail_promotion):
            import pytest
            with pytest.raises(OSError, match="promotion failed"):
                PI.PluginInstallerWindow._overwrite_folder_plugin(str(source), str(target))
        assert (target / "__init__.py").read_text() == "old"
        assert not list(tmp_path.glob(".plugin-update-*"))

    def test_failed_rollback_retains_recovery_copy(self, tmp_path):
        target = self._make_installed(tmp_path, {"__init__.py": "old"})
        source = self._make_source(tmp_path, {"__init__.py": "new"})
        real_replace = PI.os.replace
        def fail_after_backup(src, dst):
            if Path(src).name in ("replacement", "backup"):
                raise OSError("rename failed")
            return real_replace(src, dst)
        with patch.object(PI.os, "replace", side_effect=fail_after_backup):
            import pytest
            with pytest.raises(OSError, match="rename failed"):
                PI.PluginInstallerWindow._overwrite_folder_plugin(str(source), str(target))
        backups = list(tmp_path.glob(".plugin-update-*/backup/__init__.py"))
        assert len(backups) == 1
        assert backups[0].read_text() == "old"

    def _make_installed(self, tmp_path, files):
        tgt = tmp_path / "installed_plugin"
        tgt.mkdir()
        for name, content in files.items():
            (tgt / name).write_text(content, encoding="utf-8")
        return tgt

    def test_removes_orphaned_file_on_update(self, tmp_path):
        # Installed copy has a stale LICENSE.txt no longer shipped upstream.
        target = self._make_installed(
            tmp_path,
            {"__init__.py": "old", "LICENSE.txt": "stale", "LICENSE": "gpl"},
        )
        source = self._make_source(tmp_path, {"__init__.py": "new", "LICENSE": "gpl"})

        PI.PluginInstallerWindow._overwrite_folder_plugin(str(source), str(target))

        assert not (target / "LICENSE.txt").exists()  # orphan removed
        assert (target / "LICENSE").exists()
        assert (target / "__init__.py").read_text(encoding="utf-8") == "new"

    def test_preserves_user_settings_json(self, tmp_path):
        target = self._make_installed(
            tmp_path,
            {"__init__.py": "old", "settings.json": '{"user": "keepme"}'},
        )
        source = self._make_source(tmp_path, {"__init__.py": "new"})

        PI.PluginInstallerWindow._overwrite_folder_plugin(str(source), str(target))

        assert (target / "settings.json").read_text(
            encoding="utf-8"
        ) == '{"user": "keepme"}'
        assert (target / "__init__.py").read_text(encoding="utf-8") == "new"

    def test_no_settings_json_is_fine(self, tmp_path):
        target = self._make_installed(tmp_path, {"__init__.py": "old"})
        source = self._make_source(
            tmp_path, {"__init__.py": "new", "helper.py": "x = 1"}
        )

        PI.PluginInstallerWindow._overwrite_folder_plugin(str(source), str(target))

        assert not (target / "settings.json").exists()
        assert (target / "helper.py").read_text(encoding="utf-8") == "x = 1"

    def test_nested_orphan_directory_removed(self, tmp_path):
        target = self._make_installed(tmp_path, {"__init__.py": "old"})
        (target / "old_subdir").mkdir()
        (target / "old_subdir" / "data.txt").write_text("gone", encoding="utf-8")
        source = self._make_source(tmp_path, {"__init__.py": "new"})

        PI.PluginInstallerWindow._overwrite_folder_plugin(str(source), str(target))

        assert not (target / "old_subdir").exists()


# ---------------------------------------------------------------------------
# _is_within_directory — zip-slip / path-traversal guard
# ---------------------------------------------------------------------------


class TestIsWithinDirectory:
    def test_file_directly_inside(self, tmp_path):
        base = tmp_path / "extracted"
        assert PI._is_within_directory(str(base), str(base / "a.py")) is True

    def test_nested_inside(self, tmp_path):
        base = tmp_path / "extracted"
        assert PI._is_within_directory(str(base), str(base / "sub" / "a.py")) is True

    def test_sibling_prefix_rejected(self, tmp_path):
        # The bug a plain startswith() has: /x/extracted_evil is NOT inside /x/extracted
        base = tmp_path / "extracted"
        sibling = tmp_path / "extracted_evil" / "a.py"
        assert PI._is_within_directory(str(base), str(sibling)) is False

    def test_parent_traversal_rejected(self, tmp_path):
        base = tmp_path / "extracted"
        escaped = os.path.join(str(base), "..", "..", "evil.py")
        assert PI._is_within_directory(str(base), escaped) is False


# ---------------------------------------------------------------------------
# _is_allowed_download_url — scheme allowlist
# ---------------------------------------------------------------------------


class TestIsAllowedDownloadUrl:
    def test_https_allowed(self):
        assert PI._is_allowed_download_url("https://github.com/x/y/z.zip") is True

    def test_http_allowed(self):
        assert PI._is_allowed_download_url("http://example.com/z.py") is True

    def test_file_scheme_rejected(self):
        assert PI._is_allowed_download_url("file:///etc/passwd") is False

    def test_ftp_scheme_rejected(self):
        assert PI._is_allowed_download_url("ftp://example.com/z.zip") is False

    def test_data_scheme_rejected(self):
        assert PI._is_allowed_download_url("data:text/plain;base64,QQ==") is False

    def test_empty_rejected(self):
        assert PI._is_allowed_download_url("") is False


# ---------------------------------------------------------------------------
# on_update_clicked — dependency warning suppressed when already installed
# ---------------------------------------------------------------------------


class TestUpdateSkipsDependencyWarning:
    """Updating an installed plugin must not nag about missing libraries."""

    def _make_installer(self):
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.main_window = MagicMock()
        inst.main_window.plugin_manager.plugins = []
        inst.main_window.plugin_manager.install_plugin.return_value = (True, "Installed")
        inst.remote_data = []  # -> remote_info None, skips app-compat block
        inst._batch_updating = False
        inst._last_install_succeeded = False
        return inst

    def _btn(self, target_file, deps, optional_deps=None):
        props = {
            "plugin_name": "Demo",
            "download_url": "https://example.com/demo.zip",
            "target_file": target_file,
            "dependencies": deps,
            "optional_dependencies": optional_deps or [],
        }
        btn = MagicMock()
        btn.property.side_effect = lambda k: props.get(k)
        return btn

    @staticmethod
    def _question_titles(mock_question):
        titles = []
        for call in mock_question.call_args_list:
            if len(call.args) > 1:
                titles.append(call.args[1])
        return titles

    def test_no_dependency_warning_when_already_installed(self, tmp_path, monkeypatch):
        installed = tmp_path / "demo.py"
        installed.write_text("x = 1", encoding="utf-8")
        inst = self._make_installer()
        btn = self._btn(str(installed), ["definitely_missing_pkg_xyz"])
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)

        with (
            patch.object(PI.QMessageBox, "question") as q,
            patch.object(PI.QMessageBox, "warning"),
        ):
            inst.on_update_clicked()  # returns at the confirm dialog

        assert not any("Missing Dependencies" in t for t in self._question_titles(q))

    def test_dependency_warning_shown_on_fresh_install(self, tmp_path, monkeypatch):
        inst = self._make_installer()
        # target_file does not exist -> is_installed False -> warning applies
        btn = self._btn(str(tmp_path / "missing.py"), ["definitely_missing_pkg_xyz"])
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)

        with (
            patch.object(
                PI.QMessageBox,
                "question",
                return_value=PI.QMessageBox.StandardButton.No,
            ) as q,
            patch.object(PI, "PluginDetailsDialog"),
        ):
            inst.on_update_clicked()

        assert any("Missing Dependencies" in t for t in self._question_titles(q))

    def test_missing_optional_dependency_does_not_warn(self, tmp_path, monkeypatch):
        """An optional package is never a reason to interrupt a fresh install."""
        inst = self._make_installer()
        btn = self._btn(
            str(tmp_path / "missing.py"),
            [],
            optional_deps=["definitely_missing_pkg_xyz"],
        )
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)

        with (
            patch.object(
                PI.QMessageBox,
                "question",
                return_value=PI.QMessageBox.StandardButton.No,
            ) as q,
            patch.object(PI, "PluginDetailsDialog"),
        ):
            inst.on_update_clicked()

        assert not any("Missing Dependencies" in t for t in self._question_titles(q))

    def test_optional_dependencies_reach_the_details_dialog(
        self, tmp_path, monkeypatch
    ):
        inst = self._make_installer()
        btn = self._btn(
            str(tmp_path / "missing.py"),
            ["definitely_missing_pkg_xyz"],
            optional_deps=["optional_pkg_xyz"],
        )
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)

        with (
            patch.object(
                PI.QMessageBox,
                "question",
                return_value=PI.QMessageBox.StandardButton.No,
            ),
            patch.object(PI, "PluginDetailsDialog") as dlg,
        ):
            inst.on_update_clicked()

        assert dlg.call_args.args[-1] == ["optional_pkg_xyz"]


class _Item:
    """Recorder stand-in for QTableWidgetItem."""

    def __init__(self, text=""):
        self._text = str(text)
        self.tooltip = None
        self.background = None
        self.data_by_role = {}

    def setData(self, role, value):
        self.data_by_role[role] = value

    def data(self, role):
        return self.data_by_role.get(role)

    def setBackground(self, color):
        self.background = color

    def setToolTip(self, tip):
        self.tooltip = tip

    def text(self):
        return self._text


class _FakeTable:
    """Recorder stand-in for the installer's QTableWidget."""

    def __init__(self):
        self.rows = 0
        self.items = {}  # (row, col) -> _Item
        self.cell_widgets = {}  # (row, col) -> widget
        self.hidden = {}  # row -> bool

    def setRowCount(self, n):
        self.rows = n
        if n == 0:
            self.items.clear()
            self.cell_widgets.clear()

    def rowCount(self):
        return self.rows

    def insertRow(self, row):
        self.rows += 1

    def setItem(self, row, col, item):
        self.items[(row, col)] = item

    def item(self, row, col):
        return self.items.get((row, col))

    def setCellWidget(self, row, col, widget):
        self.cell_widgets[(row, col)] = widget

    def setUpdatesEnabled(self, flag):
        pass

    def setRowHidden(self, row, hidden):
        self.hidden[row] = hidden


class TestPopulateTable:
    def _make_installer(self, remote_data, installed, app_version="5.0.0a1"):
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.table = _FakeTable()
        inst.latest_app_version = "Unknown"  # skips the app-upgrade block
        inst.remote_data = remote_data
        inst.main_window = MagicMock()
        inst.main_window.plugin_manager.plugins = installed
        inst._pending_installs = {}
        inst.btn_upgrade_app = MagicMock()
        inst.btn_update_all = None
        inst.get_app_version = lambda: app_version
        inst._update_status_label = MagicMock()
        return inst

    def _statuses(self, inst):
        """Return {plugin name: status text} from the fake table."""
        out = {}
        for (row, col), item in inst.table.items.items():
            if col == 4:
                name_item = inst.table.items.get((row, 0))
                out[name_item.text()] = item.text()
        return out

    def _run(self, inst, monkeypatch):
        monkeypatch.setattr(PI, "QTableWidgetItem", _Item)
        inst.populate_table(silent=True)

    def test_invisible_uninstalled_plugin_skipped(self, monkeypatch):
        remote = [
            {
                "name": "Hidden",
                "version": "1.0.0",
                "visible": False,
                "downloadUrl": "https://x/h.py",
            },
        ]
        inst = self._make_installer(remote, installed=[])
        self._run(inst, monkeypatch)
        assert inst.table.rows == 0
        inst._update_status_label.assert_called_once_with(0)

    def test_incompatible_uninstalled_shows_incompatible_status(self, monkeypatch):
        remote = [
            {
                "name": "Future",
                "version": "9.0.0",
                "visible": True,
                "supported_moleditpy_version": ">=5.0.0",
                "downloadUrl": "https://x/f.py",
            },
        ]
        inst = self._make_installer(remote, installed=[], app_version="4.2.0")
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Future": "Incompatible"}
        status_item = inst.table.items[(0, 4)]
        assert ">=5.0.0" in status_item.tooltip
        assert "4.2.0" in status_item.tooltip

    def test_compatible_uninstalled_gets_install_button(self, monkeypatch):
        remote = [
            {
                "name": "Fresh",
                "version": "1.0.0",
                "visible": True,
                "supported_moleditpy_version": ">=5.0.0, <6.0.0",
                "downloadUrl": "https://x/fresh.py",
            },
        ]
        inst = self._make_installer(remote, installed=[])
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Fresh": "Not Installed"}
        assert (0, 5) in inst.table.cell_widgets  # Install button placed

    def test_update_available_counted(self, monkeypatch):
        remote = [
            {
                "name": "Old",
                "version": "2.0.0",
                "visible": True,
                "downloadUrl": "https://x/old.py",
            },
        ]
        installed = [{"name": "Old", "version": "1.0.0", "filepath": None}]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Old": "Update Available"}
        assert inst.updates_found is True
        inst._update_status_label.assert_called_once_with(1)

    def test_up_to_date_not_counted(self, monkeypatch):
        remote = [
            {
                "name": "Same",
                "version": "1.0.0",
                "visible": True,
                "downloadUrl": "https://x/same.py",
            },
        ]
        installed = [{"name": "Same", "version": "1.0.0", "filepath": None}]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Same": "Up to date"}
        inst._update_status_label.assert_called_once_with(0)

    def test_local_newer_than_registry(self, monkeypatch):
        remote = [
            {
                "name": "Dev",
                "version": "1.0.0",
                "visible": True,
                "downloadUrl": "https://x/dev.py",
            },
        ]
        installed = [{"name": "Dev", "version": "2.0.0", "filepath": None}]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Dev": "Newer"}

    def test_installed_but_not_in_registry(self, monkeypatch):
        installed = [{"name": "Private", "version": "0.1.0", "filepath": None}]
        inst = self._make_installer([], installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Private": "Not in Registry"}
        assert (0, 5) not in inst.table.cell_widgets  # no action button

    def test_local_version_read_from_file(self, monkeypatch, tmp_path):
        f = tmp_path / "p.py"
        f.write_text('PLUGIN_VERSION = "3.3.3"\n', encoding="utf-8")
        remote = [
            {
                "name": "OnDisk",
                "version": "3.3.3",
                "visible": True,
                "downloadUrl": "https://x/p.py",
            },
        ]
        installed = [{"name": "OnDisk", "version": "0.0.1", "filepath": str(f)}]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        # AST-read version (3.3.3) wins over the stale metadata (0.0.1)
        assert self._statuses(inst) == {"OnDisk": "Up to date"}
        assert inst.table.items[(0, 2)].text() == "3.3.3"

    def test_disabled_plugin_shows_disabled(self, monkeypatch):
        remote = [
            {
                "name": "Off",
                "version": "1.0.0",
                "visible": True,
                "downloadUrl": "https://x/off.py",
            },
        ]
        installed = [
            {
                "name": "Off",
                "version": "1.0.0",
                "filepath": None,
                "status": "Disabled",
                "disabled": True,
            }
        ]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Off": "Disabled"}
        item = inst.table.items[(0, 4)]
        assert item.data(PI.Qt.ItemDataRole.UserRole) == "Up to date"
        assert "Plugin Manager" in item.tooltip
        inst._update_status_label.assert_called_once_with(0)

    def test_disabled_plugin_with_update_still_counted(self, monkeypatch):
        remote = [
            {
                "name": "Off",
                "version": "2.0.0",
                "visible": True,
                "downloadUrl": "https://x/off.py",
            },
        ]
        installed = [
            {"name": "Off", "version": "1.0.0", "filepath": None, "disabled": True}
        ]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Off": "Update Available (Disabled)"}
        assert inst.updates_found is True
        inst._update_status_label.assert_called_once_with(1)


class TestProjectPageAndDisabledHelpers:
    def test_project_url_github(self):
        url = "https://github.com/HiroYokoyama/moleditpy_pyscf-calculator"
        assert PI._project_page_url({"projectUrl": url}) == url
        assert PI._project_page_label(url) == "Open GitHub Page"

    def test_project_url_rejects_unsafe_scheme(self):
        assert PI._project_page_url({"projectUrl": "file:///C:/x"}) is None
        assert PI._project_page_url({"projectUrl": ""}) is None
        assert PI._project_page_url(None) is None

    def test_project_label_other_host(self):
        assert PI._project_page_label("https://gitlab.com/x") == "Open Project Page"

    def test_is_plugin_disabled(self):
        assert PI._is_plugin_disabled({"disabled": True}) is True
        assert PI._is_plugin_disabled({"status": "Disabled"}) is True
        assert PI._is_plugin_disabled({"status": "Loaded"}) is False
        assert PI._is_plugin_disabled(None) is False


class TestFilterPlugins:
    def _make_installer(self, rows):
        """rows: list of (name, author)."""
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.table = _FakeTable()
        inst.table.rows = len(rows)
        for r, (name, author) in enumerate(rows):
            inst.table.items[(r, 0)] = _Item(name)
            inst.table.items[(r, 1)] = _Item(author)
        inst.search_input = MagicMock()
        return inst

    def test_empty_filter_shows_all(self):
        inst = self._make_installer([("Alpha", "Ann"), ("Beta", "Bob")])
        inst.search_input.text.return_value = ""
        inst.filter_plugins()
        assert inst.table.hidden == {0: False, 1: False}

    def test_matches_name_case_insensitive(self):
        inst = self._make_installer([("Alpha", "Ann"), ("Beta", "Bob")])
        inst.search_input.text.return_value = "ALPHA"
        inst.filter_plugins()
        assert inst.table.hidden == {0: False, 1: True}

    def test_matches_author(self):
        inst = self._make_installer([("Alpha", "Ann"), ("Beta", "Bob")])
        inst.search_input.text.return_value = "bob"
        inst.filter_plugins()
        assert inst.table.hidden == {0: True, 1: False}

    def test_no_match_hides_all(self):
        inst = self._make_installer([("Alpha", "Ann"), ("Beta", "Bob")])
        inst.search_input.text.return_value = "zzz"
        inst.filter_plugins()
        assert inst.table.hidden == {0: True, 1: True}

    def test_missing_items_treated_as_empty(self):
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.table = _FakeTable()
        inst.table.rows = 1  # row exists but has no items
        inst.search_input = MagicMock()
        inst.search_input.text.return_value = "x"
        inst.filter_plugins()  # must not raise
        assert inst.table.hidden == {0: True}


def test_dependency_commands_reject_pip_options():
    for dependency in ("--upgrade", "--no-deps", "--pre", "-r", ".", ".."):
        assert PI.sanitize_and_quote_dependency(dependency) == ""
    assert PI.sanitize_and_quote_dependency("numpy") == "numpy"
    assert PI.sanitize_and_quote_dependency("scikit-learn>=1.0") == '"scikit-learn>=1.0"'


class TestOnUpdateClickedSecurity:
    """Exercises the SHA256-verification / download-URL branches of
    on_update_clicked (the actual install path) rather than the earlier
    dependency-warning gate."""

    def _make_installer(self):
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.main_window = MagicMock()
        inst.main_window.plugin_manager.plugins = []
        inst.main_window.plugin_manager.install_plugin.return_value = (True, 'Installed')
        inst._batch_updating = False
        inst._last_install_succeeded = False
        inst._pending_installs = {}
        inst.populate_table = MagicMock()
        return inst

    def test_host_failure_keeps_old_plugin_and_reports_failure(self, tmp_path, monkeypatch):
        import hashlib
        content = b"zip payload"
        target = tmp_path / "old.py"
        target.write_text("old plugin")
        btn, entry = self._btn("Demo", "https://example.com/demo.zip", str(target),
                               sha256=hashlib.sha256(content).hexdigest())
        inst = self._make_installer()
        inst.remote_data = [entry]
        inst.main_window.plugin_manager.install_plugin.return_value = (False, "disk full")
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)
        def download(url, path, cb=None):
            Path(path).write_bytes(content)
            return True
        monkeypatch.setattr(inst, "_download_chunked", download)
        with patch.object(PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes), \
             patch.object(PI.QMessageBox, "warning") as warning, \
             patch.object(PI.QMessageBox, "information") as information:
            inst.on_update_clicked()
        assert target.read_text() == "old plugin"
        assert inst._last_install_succeeded is False
        assert inst._pending_installs == {}
        warning.assert_called_once()
        assert "disk full" in warning.call_args[0][2]
        information.assert_not_called()

    def _btn(self, plugin_name, download_url, target_file, sha256=None):
        props = {
            "plugin_name": plugin_name,
            "download_url": download_url,
            "target_file": target_file,
            "dependencies": [],
        }
        btn = MagicMock()
        btn.property.side_effect = lambda k: props.get(k)
        entry = {
            "name": plugin_name,
            "downloadUrl": download_url,
            "version": "1.0.0",
            "author": "Someone",
            "description": "desc",
        }
        if sha256 is not None:
            entry["sha256"] = sha256
        return btn, entry

    def test_sha256_mismatch_blocks_install(self, tmp_path, monkeypatch):
        content = b"actual downloaded bytes"
        wrong_hash = "0" * 64
        btn, entry = self._btn(
            "Demo",
            "https://example.com/demo.py",
            str(tmp_path / "missing.py"),
            sha256=wrong_hash,
        )
        inst = self._make_installer()
        inst.remote_data = [entry]
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)
        monkeypatch.setattr(
            inst,
            "_download_chunked",
            lambda url, path, cb=None: (open(path, "wb").write(content) or True)
            and True,
            raising=False,
        )

        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes
        ):
            with patch.object(PI.QMessageBox, "critical") as crit:
                inst.on_update_clicked()

        crit.assert_called_once()
        title = crit.call_args[0][1]
        assert "Security" in title
        inst.main_window.plugin_manager.install_plugin.assert_not_called()
        assert inst._last_install_succeeded is False

    def test_sha256_match_allows_install(self, tmp_path, monkeypatch):
        import hashlib

        content = b"actual downloaded bytes"
        correct_hash = hashlib.sha256(content).hexdigest()
        btn, entry = self._btn(
            "Demo",
            "https://example.com/demo.py",
            str(tmp_path / "missing.py"),
            sha256=correct_hash,
        )
        inst = self._make_installer()
        inst.remote_data = [entry]
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)
        monkeypatch.setattr(
            inst,
            "_download_chunked",
            lambda url, path, cb=None: open(path, "wb").write(content) or True,
            raising=False,
        )

        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes
        ):
            with patch.object(PI.QMessageBox, "information"):
                inst.on_update_clicked()

        inst.main_window.plugin_manager.install_plugin.assert_called_once()
        assert inst._last_install_succeeded is True

    def test_missing_sha256_prompt_no_aborts(self, tmp_path, monkeypatch):
        btn, entry = self._btn(
            "Demo", "https://example.com/demo.py", str(tmp_path / "missing.py")
        )  # no sha256 key
        inst = self._make_installer()
        inst.remote_data = [entry]
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)
        monkeypatch.setattr(
            inst,
            "_download_chunked",
            lambda url, path, cb=None: open(path, "wb").write(b"x") or True,
            raising=False,
        )

        # First question() = confirm install -> Yes.
        # Second question() would be missing-sha warning is actually QMessageBox.warning,
        # not question -- patch warning to return No.
        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes
        ):
            with patch.object(
                PI.QMessageBox, "warning", return_value=PI.QMessageBox.StandardButton.No
            ) as warn:
                inst.on_update_clicked()

        warn.assert_called_once()
        assert "checksum" in warn.call_args[0][2].lower()
        inst.main_window.plugin_manager.install_plugin.assert_not_called()

    def test_missing_sha256_prompt_yes_proceeds(self, tmp_path, monkeypatch):
        btn, entry = self._btn(
            "Demo", "https://example.com/demo.py", str(tmp_path / "missing.py")
        )  # no sha256 key
        inst = self._make_installer()
        inst.remote_data = [entry]
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)
        monkeypatch.setattr(
            inst,
            "_download_chunked",
            lambda url, path, cb=None: open(path, "wb").write(b"x") or True,
            raising=False,
        )

        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes
        ):
            with patch.object(
                PI.QMessageBox,
                "warning",
                return_value=PI.QMessageBox.StandardButton.Yes,
            ):
                with patch.object(PI.QMessageBox, "information"):
                    inst.on_update_clicked()

        inst.main_window.plugin_manager.install_plugin.assert_called_once()
        assert inst._last_install_succeeded is True

    def test_download_failure_skips_sha_and_install(self, tmp_path, monkeypatch):
        btn, entry = self._btn(
            "Demo",
            "https://example.com/demo.py",
            str(tmp_path / "missing.py"),
            sha256="deadbeef",
        )
        inst = self._make_installer()
        inst.remote_data = [entry]
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)
        monkeypatch.setattr(
            inst, "_download_chunked", lambda *a, **k: False, raising=False
        )
        calc_mock = MagicMock()
        monkeypatch.setattr(inst, "calculate_sha256", calc_mock, raising=False)

        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes
        ):
            inst.on_update_clicked()

        calc_mock.assert_not_called()
        inst.main_window.plugin_manager.install_plugin.assert_not_called()

    def test_disallowed_url_scheme_rejected_before_download(
        self, tmp_path, monkeypatch
    ):
        # Starts with "http" (passes the naive startswith check) but the scheme
        # itself ("httpevil") is not in the http/https allowlist.
        btn, entry = self._btn(
            "Demo",
            "httpevil://example.com/demo.py",
            str(tmp_path / "missing.py"),
            sha256="deadbeef",
        )
        inst = self._make_installer()
        inst.remote_data = [entry]
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)
        dl_mock = MagicMock()
        monkeypatch.setattr(inst, "_download_chunked", dl_mock, raising=False)

        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes
        ):
            with patch.object(PI.QMessageBox, "warning") as warn:
                inst.on_update_clicked()

        dl_mock.assert_not_called()
        warn.assert_called_once()
        assert (
            "unsafe" in warn.call_args[0][2].lower()
            or "scheme" in warn.call_args[0][2].lower()
        )

    def test_no_sender_returns_immediately(self):
        inst = self._make_installer()
        inst.remote_data = []
        with patch.object(inst, "sender", return_value=None, create=True):
            inst.on_update_clicked()  # must not raise
        inst.main_window.plugin_manager.install_plugin.assert_not_called()

    def test_no_download_url_warns_and_returns(self, tmp_path, monkeypatch):
        btn, entry = self._btn("Demo", "", str(tmp_path / "missing.py"))
        inst = self._make_installer()
        inst.remote_data = [entry]
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)

        with patch.object(PI.QMessageBox, "warning") as warn:
            inst.on_update_clicked()

        warn.assert_called_once()
        inst.main_window.plugin_manager.install_plugin.assert_not_called()

    def test_incompatible_version_no_aborts(self, tmp_path, monkeypatch):
        btn, entry = self._btn(
            "Demo",
            "https://example.com/demo.py",
            str(tmp_path / "missing.py"),
            sha256="deadbeef",
        )
        entry["supported_moleditpy_version"] = ">=99.0.0"
        inst = self._make_installer()
        inst.remote_data = [entry]
        inst.get_app_version = lambda: "4.2.0"
        monkeypatch.setattr(inst, "sender", lambda: btn, raising=False)
        dl_mock = MagicMock()
        monkeypatch.setattr(inst, "_download_chunked", dl_mock, raising=False)

        with patch.object(
            PI.QMessageBox, "warning", return_value=PI.QMessageBox.StandardButton.No
        ):
            inst.on_update_clicked()

        dl_mock.assert_not_called()
        inst.main_window.plugin_manager.install_plugin.assert_not_called()


class TestPluginDetailsDialogActions:
    def _dialog(self, missing_deps=None, target_file=None, plugin_name="Demo"):
        d = object.__new__(PI.PluginDetailsDialog)
        d.plugin_name = plugin_name
        d.target_file = target_file
        d.parent_installer = MagicMock()
        if missing_deps is not None:
            d.missing_deps = missing_deps
        return d

    def test_copy_install_all_no_missing_deps_is_noop(self):
        d = self._dialog(missing_deps=None)
        with patch.object(PI.QApplication, "clipboard") as clip:
            d.copy_install_all_command()
        clip.assert_not_called()

    def test_copy_install_all_builds_pip_command(self):
        d = self._dialog(missing_deps=["numpy>=1.20", "rdkit"])
        clipboard = MagicMock()
        with patch.object(PI.QApplication, "clipboard", return_value=clipboard):
            with patch.object(PI.QMessageBox, "information") as info:
                d.copy_install_all_command()
        clipboard.setText.assert_called_once()
        cmd = clipboard.setText.call_args[0][0]
        assert cmd == 'pip install "numpy>=1.20" rdkit'
        info.assert_called_once()

    def test_copy_install_all_skips_unsafe_dep_with_warning(self):
        d = self._dialog(missing_deps=["numpy;pkg", "rdkit"])
        clipboard = MagicMock()
        with patch.object(PI.QApplication, "clipboard", return_value=clipboard):
            with patch.object(PI.QMessageBox, "warning") as warn:
                with patch.object(PI.QMessageBox, "information"):
                    d.copy_install_all_command()
        warn.assert_called_once()
        cmd = clipboard.setText.call_args[0][0]
        assert cmd == "pip install rdkit"

    def test_copy_optional_install_command_uses_optional_list(self):
        d = self._dialog(missing_deps=["numpy"])
        d.missing_optional_deps = ["matplotlib>=3.5"]
        clipboard = MagicMock()
        with patch.object(PI.QApplication, "clipboard", return_value=clipboard):
            with patch.object(PI.QMessageBox, "information"):
                d.copy_optional_install_command()
        assert clipboard.setText.call_args[0][0] == 'pip install "matplotlib>=3.5"'

    def test_copy_optional_install_command_without_optional_deps_is_noop(self):
        d = self._dialog(missing_deps=["numpy"])
        with patch.object(PI.QApplication, "clipboard") as clip:
            d.copy_optional_install_command()
        clip.assert_not_called()

    def test_copy_install_all_all_unsafe_sets_nothing(self):
        d = self._dialog(missing_deps=["numpy;pkg", "rdkit&x"])
        clipboard = MagicMock()
        with patch.object(PI.QApplication, "clipboard", return_value=clipboard):
            with patch.object(PI.QMessageBox, "warning"):
                d.copy_install_all_command()
        clipboard.setText.assert_not_called()

    def test_on_delete_user_cancels_does_not_delete(self, tmp_path):
        f = tmp_path / "demo.py"
        f.write_text("x = 1")
        d = self._dialog(target_file=str(f))
        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.No
        ):
            d.on_delete()
        assert f.exists()

    def test_on_delete_removes_single_file(self, tmp_path):
        f = tmp_path / "demo.py"
        f.write_text("x = 1")
        d = self._dialog(target_file=str(f))
        d.accept = MagicMock()
        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes
        ):
            with patch.object(PI.QMessageBox, "information"):
                d.on_delete()
        assert not f.exists()
        d.accept.assert_called_once()
        d.parent_installer.main_window.plugin_manager.discover_plugins.assert_called_once()
        d.parent_installer.check_updates.assert_called_once()

    def test_on_delete_removes_folder_plugin(self, tmp_path):
        pkg = tmp_path / "My_Plugin"
        pkg.mkdir()
        init = pkg / "__init__.py"
        init.write_text("x = 1")
        d = self._dialog(target_file=str(init))
        d.accept = MagicMock()
        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes
        ):
            with patch.object(PI.QMessageBox, "information"):
                d.on_delete()
        assert not pkg.exists()

    def test_on_delete_missing_file_warns(self, tmp_path):
        d = self._dialog(target_file=str(tmp_path / "gone.py"))
        with patch.object(
            PI.QMessageBox, "question", return_value=PI.QMessageBox.StandardButton.Yes
        ):
            with patch.object(PI.QMessageBox, "warning") as warn:
                d.on_delete()
        warn.assert_called_once()
        d.parent_installer.check_updates.assert_not_called()


def _fake_pkg(monkeypatch, top, version):
    """Install a fake <top>.utils.constants module tree with VERSION."""
    pkg = types.ModuleType(top)
    pkg.__path__ = []
    utils = types.ModuleType(f"{top}.utils")
    utils.__path__ = []
    constants = types.ModuleType(f"{top}.utils.constants")
    constants.VERSION = version
    pkg.utils = utils
    utils.constants = constants
    monkeypatch.setitem(sys.modules, top, pkg)
    monkeypatch.setitem(sys.modules, f"{top}.utils", utils)
    monkeypatch.setitem(sys.modules, f"{top}.utils.constants", constants)


def _block_pkg(monkeypatch, top):
    """Make ``import <top>`` raise ImportError even if really installed."""
    monkeypatch.setitem(sys.modules, top, None)


class TestAppVersionDetection:
    def _inst(self):
        return object.__new__(PI.PluginInstallerWindow)

    def test_get_app_version_from_moleditpy(self, monkeypatch):
        _fake_pkg(monkeypatch, "moleditpy", "4.9.9")
        assert self._inst().get_app_version() == "4.9.9"

    def test_get_app_version_from_linux_variant(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _fake_pkg(monkeypatch, "moleditpy_linux", "4.8.0")
        assert self._inst().get_app_version() == "4.8.0"

    def test_get_app_version_from_main_module(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        fake_main = types.ModuleType("__main__")
        fake_main.VERSION = "5.0.0"
        monkeypatch.setitem(sys.modules, "__main__", fake_main)
        assert self._inst().get_app_version() == "5.0.0"

    def test_get_app_version_ignores_unknown_version(self, monkeypatch):
        _fake_pkg(monkeypatch, "moleditpy", "Unknown")
        fake_main = types.ModuleType("__main__")
        fake_main.VERSION = "Unknown"
        monkeypatch.setitem(sys.modules, "__main__", fake_main)
        assert self._inst().get_app_version() != "Unknown"

    def test_get_app_version_parses_pyproject_toml(self, monkeypatch, tmp_path):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        fake_main = types.ModuleType("__main__")
        fake_main.VERSION = "Unknown"
        monkeypatch.setitem(sys.modules, "__main__", fake_main)

        proj_dir = tmp_path / "myproj"
        proj_dir.mkdir()
        toml_file = proj_dir / "pyproject.toml"
        toml_file.write_text(
            '[project]\nname = "test"\nversion = "3.2.1"\n', encoding="utf-8"
        )

        script_file = proj_dir / "src" / "main.py"
        script_file.parent.mkdir(parents=True)
        script_file.write_text("# main", encoding="utf-8")

        monkeypatch.setattr(sys, "argv", [str(script_file)])
        assert self._inst().get_app_version() == "3.2.1"

    def test_get_app_version_fallback_zero(self, monkeypatch, tmp_path):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        fake_main = types.ModuleType("__main__")
        fake_main.VERSION = "Unknown"
        monkeypatch.setitem(sys.modules, "__main__", fake_main)
        monkeypatch.setattr(sys, "argv", [])
        monkeypatch.chdir(tmp_path)
        assert self._inst().get_app_version() == "0.0.0"

    def test_package_name_moleditpy(self, monkeypatch):
        _fake_pkg(monkeypatch, "moleditpy", "4.9.9")
        assert self._inst()._get_package_name() == "moleditpy"

    def test_package_name_linux(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _fake_pkg(monkeypatch, "moleditpy_linux", "4.8.0")
        assert self._inst()._get_package_name() == "moleditpy-linux"

    def test_package_name_default_when_neither(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        assert self._inst()._get_package_name() == "moleditpy"

    def test_get_app_version_fallback_appends_sys_path(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        fake_argv = os.path.abspath("/fake/src/moleditpy_linux/main.py")
        monkeypatch.setattr(sys, "argv", [fake_argv])
        fake_path = []
        monkeypatch.setattr(sys, "path", fake_path)

        self._inst().get_app_version()

        normalized_paths = [os.path.normpath(p) for p in fake_path]
        assert os.path.dirname(fake_argv) in normalized_paths
        assert os.path.dirname(os.path.dirname(fake_argv)) in normalized_paths

    def test_get_package_name_fallback_appends_sys_path(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        fake_argv = os.path.abspath("/fake/src/moleditpy_linux/main.py")
        monkeypatch.setattr(sys, "argv", [fake_argv])
        fake_path = []
        monkeypatch.setattr(sys, "path", fake_path)

        self._inst()._get_package_name()

        normalized_paths = [os.path.normpath(p) for p in fake_path]
        assert os.path.dirname(fake_argv) in normalized_paths
        assert os.path.dirname(os.path.dirname(fake_argv)) in normalized_paths


class TestPythonVersionCheck:
    def test_get_python_version_matches_interpreter(self):
        assert PI.get_python_version() == "%d.%d.%d" % sys.version_info[:3]

    def test_get_python_version_format(self):
        assert _re.fullmatch(r"\d+\.\d+\.\d+", PI.get_python_version())

    def test_constant_present(self):
        assert PI.is_app_version_compatible("3.9.0", PI.PLUGIN_SUPPORTED_PYTHON_VERSION)
        assert not PI.is_app_version_compatible("3.8.0", PI.PLUGIN_SUPPORTED_PYTHON_VERSION)

    def test_constant_matches_registry(self):
        registry = json.loads(
            (
                Path(__file__).resolve().parents[1] / "REGISTRY" / "plugins.json"
            ).read_text(encoding="utf-8-sig")
        )
        entry = next((p for p in registry if p.get("id") == "plugin_installer"), None)
        assert entry is not None
        assert (
            entry.get("supported_python_version") == PI.PLUGIN_SUPPORTED_PYTHON_VERSION
        )

    def test_default_spec_boundaries(self):
        spec = PI.PLUGIN_SUPPORTED_PYTHON_VERSION
        assert PI.is_app_version_compatible("3.9.0", spec) is True
        assert PI.is_app_version_compatible("3.12.4", spec) is True
        assert PI.is_app_version_compatible("3.14.2", spec) is True
        assert PI.is_app_version_compatible("3.8.19", spec) is False
        assert PI.is_app_version_compatible("3.15.0", spec) is False

    def test_running_interpreter_compatible_with_default(self):
        assert (
            PI.is_app_version_compatible(PI.get_python_version(), PI.PLUGIN_SUPPORTED_PYTHON_VERSION)
            is True
        )

    def test_missing_spec_is_compatible(self):
        assert PI.is_app_version_compatible(PI.get_python_version(), "") is True


# ---------------------------------------------------------------------------
# pyproject.toml parsing (source-checkout version fallback)
# ---------------------------------------------------------------------------


class TestVersionFromPyprojectText:
    def test_reads_project_version(self):
        text = '[project]\nname = "moleditpy"\nversion = "4.2.0"\n'
        assert PI._version_from_pyproject_text(text) == "4.2.0"

    def test_ignores_version_in_other_tables(self):
        """A [tool.*] version belongs to that tool, not to the application.

        Regression: the section filter was dead code, so the first `version =`
        line anywhere in the file won — here 9.9.9 instead of the real 4.2.0.
        """
        text = (
            '[build-system]\nrequires = ["setuptools"]\n\n'
            '[tool.commitizen]\nversion = "9.9.9"\n\n'
            '[project]\nname = "moleditpy"\nversion = "4.2.0"\n'
        )
        assert PI._version_from_pyproject_text(text) == "4.2.0"

    def test_ignores_version_after_the_project_table_ends(self):
        text = '[project]\nname = "x"\n\n[tool.poetry]\nversion = "9.9.9"\n'
        assert PI._version_from_pyproject_text(text) is None

    def test_ignores_preamble_before_any_table(self):
        text = 'version = "9.9.9"\n\n[project]\nname = "x"\n'
        assert PI._version_from_pyproject_text(text) is None

    def test_accepts_single_quotes_and_no_spaces(self):
        assert (
            PI._version_from_pyproject_text("[project]\nversion='1.2.3'\n") == "1.2.3"
        )

    def test_ignores_keys_that_merely_start_with_version(self):
        text = '[project]\nversion_scheme = "calver"\nversion = "1.0.0"\n'
        assert PI._version_from_pyproject_text(text) == "1.0.0"

    def test_rejects_unknown_placeholder(self):
        assert (
            PI._version_from_pyproject_text('[project]\nversion = "Unknown"\n') is None
        )

    def test_rejects_empty_value(self):
        assert PI._version_from_pyproject_text('[project]\nversion = ""\n') is None

    def test_returns_none_without_a_project_table(self):
        assert (
            PI._version_from_pyproject_text("[tool.black]\nline-length = 88\n") is None
        )

    def test_returns_none_for_empty_text(self):
        assert PI._version_from_pyproject_text("") is None


class TestVersionFromPyproject:
    def _script(self, tmp_path, rel="app/main.py"):
        script = tmp_path / rel
        script.parent.mkdir(parents=True, exist_ok=True)
        script.write_text("# main", encoding="utf-8")
        return script

    def test_finds_pyproject_next_to_the_script(self, monkeypatch, tmp_path):
        script = self._script(tmp_path, "main.py")
        (tmp_path / "pyproject.toml").write_text(
            '[project]\nversion = "4.2.0"\n', encoding="utf-8"
        )
        monkeypatch.setattr(sys, "argv", [str(script)])
        assert PI._version_from_pyproject() == "4.2.0"

    def test_walks_up_to_a_parent_directory(self, monkeypatch, tmp_path):
        script = self._script(tmp_path, "src/app/main.py")
        (tmp_path / "pyproject.toml").write_text(
            '[project]\nversion = "4.3.0"\n', encoding="utf-8"
        )
        monkeypatch.setattr(sys, "argv", [str(script)])
        assert PI._version_from_pyproject() == "4.3.0"

    def test_nearest_pyproject_wins(self, monkeypatch, tmp_path):
        script = self._script(tmp_path, "src/app/main.py")
        (tmp_path / "pyproject.toml").write_text(
            '[project]\nversion = "1.0.0"\n', encoding="utf-8"
        )
        (tmp_path / "src" / "pyproject.toml").write_text(
            '[project]\nversion = "2.0.0"\n', encoding="utf-8"
        )
        monkeypatch.setattr(sys, "argv", [str(script)])
        assert PI._version_from_pyproject() == "2.0.0"

    def test_stops_after_the_search_depth(self, monkeypatch, tmp_path):
        deep = tmp_path.joinpath(
            *[f"d{i}" for i in range(PI._PYPROJECT_SEARCH_DEPTH + 2)]
        )
        deep.mkdir(parents=True)
        script = deep / "main.py"
        script.write_text("# main", encoding="utf-8")
        (tmp_path / "pyproject.toml").write_text(
            '[project]\nversion = "4.2.0"\n', encoding="utf-8"
        )
        monkeypatch.setattr(sys, "argv", [str(script)])
        monkeypatch.chdir(deep)
        assert PI._version_from_pyproject() is None

    def test_falls_back_to_the_working_directory(self, monkeypatch, tmp_path):
        (tmp_path / "pyproject.toml").write_text(
            '[project]\nversion = "5.1.0"\n', encoding="utf-8"
        )
        monkeypatch.setattr(sys, "argv", [])
        monkeypatch.chdir(tmp_path)
        assert PI._version_from_pyproject() == "5.1.0"

    def test_unreadable_pyproject_is_skipped(self, monkeypatch, tmp_path):
        import builtins

        (tmp_path / "pyproject.toml").write_text(
            '[project]\nversion = "4.2.0"\n', encoding="utf-8"
        )
        monkeypatch.setattr(sys, "argv", [])
        monkeypatch.chdir(tmp_path)

        def _boom(*a, **kw):
            raise OSError("locked")

        monkeypatch.setattr(builtins, "open", _boom)
        assert PI._version_from_pyproject() is None

    def test_malformed_pyproject_is_skipped(self, monkeypatch, tmp_path):
        (tmp_path / "pyproject.toml").write_text("not a toml file", encoding="utf-8")
        monkeypatch.setattr(sys, "argv", [])
        monkeypatch.chdir(tmp_path)
        assert PI._version_from_pyproject() is None

    def test_directory_argv_is_handled(self, monkeypatch, tmp_path):
        (tmp_path / "pyproject.toml").write_text(
            '[project]\nversion = "6.0.0"\n', encoding="utf-8"
        )
        monkeypatch.setattr(sys, "argv", [str(tmp_path)])  # a dir, not a file
        assert PI._version_from_pyproject() == "6.0.0"


# ---------------------------------------------------------------------------
# Host-package resolution
# ---------------------------------------------------------------------------


class TestResolveAppPackage:
    def test_prefers_an_already_imported_package(self, monkeypatch):
        _fake_pkg(monkeypatch, "moleditpy", "4.9.9")
        assert PI._resolve_app_package() == "moleditpy"

    def test_finds_the_linux_build(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _fake_pkg(monkeypatch, "moleditpy_linux", "4.8.0")
        assert PI._resolve_app_package() == "moleditpy_linux"

    def test_returns_none_when_nothing_resolves(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        monkeypatch.setattr(sys, "argv", [])
        assert PI._resolve_app_package() is None

    def test_adds_source_checkout_dirs_to_sys_path(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        argv0 = os.path.abspath(
            os.path.join("fake", "src", "moleditpy_linux", "main.py")
        )
        monkeypatch.setattr(sys, "argv", [argv0])
        fake_path = []
        monkeypatch.setattr(sys, "path", fake_path)

        PI._resolve_app_package()

        normalized = [os.path.normpath(p) for p in fake_path]
        assert os.path.dirname(argv0) in normalized
        assert os.path.dirname(os.path.dirname(argv0)) in normalized

    def test_source_checkout_dirs_are_not_duplicated(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        argv0 = os.path.abspath(os.path.join("fake", "app", "main.py"))
        monkeypatch.setattr(sys, "argv", [argv0])
        fake_path = [os.path.dirname(argv0)]
        monkeypatch.setattr(sys, "path", fake_path)

        PI._resolve_app_package()

        assert fake_path.count(os.path.dirname(argv0)) == 1

    def test_no_argv_yields_no_source_dirs(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", [])
        assert PI._source_checkout_dirs() == []


class TestPackageAndVersionAgree:
    """The distribution name and the version must describe the same build."""

    def _inst(self):
        return object.__new__(PI.PluginInstallerWindow)

    def test_linux_build_reports_linux_name_and_its_version(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _fake_pkg(monkeypatch, "moleditpy_linux", "4.8.0")
        monkeypatch.delitem(sys.modules, "__main__", raising=False)
        inst = self._inst()
        assert inst._get_package_name() == "moleditpy-linux"
        assert inst.get_app_version() == "4.8.0"

    def test_pypi_build_reports_pypi_name_and_its_version(self, monkeypatch):
        _fake_pkg(monkeypatch, "moleditpy", "4.9.9")
        monkeypatch.delitem(sys.modules, "__main__", raising=False)
        inst = self._inst()
        assert inst._get_package_name() == "moleditpy"
        assert inst.get_app_version() == "4.9.9"

    def test_both_present_resolve_to_the_same_package(self, monkeypatch):
        _fake_pkg(monkeypatch, "moleditpy", "4.9.9")
        _fake_pkg(monkeypatch, "moleditpy_linux", "4.8.0")
        monkeypatch.delitem(sys.modules, "__main__", raising=False)
        inst = self._inst()
        assert inst._get_package_name() == "moleditpy"
        assert inst.get_app_version() == "4.9.9"

    def test_package_name_defaults_when_nothing_resolves(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        monkeypatch.setattr(sys, "argv", [])
        assert self._inst()._get_package_name() == "moleditpy"

    def test_pyproject_wins_over_a_tool_table_end_to_end(self, monkeypatch, tmp_path):
        """The whole fallback chain, not just the text parser."""
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        fake_main = types.ModuleType("__main__")
        fake_main.VERSION = "Unknown"
        monkeypatch.setitem(sys.modules, "__main__", fake_main)

        (tmp_path / "pyproject.toml").write_text(
            '[tool.commitizen]\nversion = "9.9.9"\n\n[project]\nversion = "4.2.0"\n',
            encoding="utf-8",
        )
        script = tmp_path / "main.py"
        script.write_text("# main", encoding="utf-8")
        monkeypatch.setattr(sys, "argv", [str(script)])

        assert self._inst().get_app_version() == "4.2.0"


class TestIsUsableVersion:
    def test_accepts_a_real_version(self):
        assert PI._is_usable_version("4.2.0")

    def test_rejects_unknown_none_and_empty(self):
        assert not PI._is_usable_version("Unknown")
        assert not PI._is_usable_version(None)
        assert not PI._is_usable_version("")


class TestRegistryRetirementMetadata:
    @pytest.mark.parametrize("entry,expected", [
        (None, False), ({}, False), ({"visible": True}, False),
        ({"visible": False}, True), ({"visible": "false"}, False),
        ({"visible": 0}, False), ([], False),
    ])
    def test_only_explicit_hidden_flag_is_stale(self, entry, expected):
        assert PI._is_registry_hidden(entry) is expected

    @pytest.mark.parametrize("entry", [None, {}, {"migration_to": None},
        {"migration_to": []}, {"migration_to": " "}, {"migration_to": 12}])
    def test_missing_or_malformed_migration_is_ignored(self, entry):
        assert PI._migration_target_name(entry, []) == ""

    def test_target_is_resolved_by_id_not_old_display_name(self):
        assert PI._migration_target_name({"migration_to": "new"}, [
            {"id": "new", "name": "Renamed Replacement"}
        ]) == "Renamed Replacement"

    def test_missing_target_keeps_reference_as_plain_text(self):
        assert PI._migration_target_name({"migration_to": " missing "}, []) == "missing"

    def test_malformed_target_name_uses_reference(self):
        assert PI._migration_target_name({"migration_to": "new"}, [None,
            {"id": "new", "name": []}]) == "new"


# ===========================================================================
# Qt-free core shared by the dialog and the command line
# ===========================================================================

import hashlib as _hashlib
import io as _io
import zipfile as _zipfile


def _sha(data: bytes) -> str:
    return _hashlib.sha256(data).hexdigest()


def _plugin_source(version, name="Alpha Tool"):
    return (
        f'PLUGIN_NAME = "{name}"\nPLUGIN_VERSION = "{version}"\n'
        "def initialize(context):\n    pass\n"
    ).encode()


def _fake_response(content: bytes, length=True):
    resp = MagicMock()
    resp.__enter__ = lambda s: s
    resp.__exit__ = MagicMock(return_value=False)
    resp.headers = {"Content-Length": str(len(content))} if length else {}
    resp.read.side_effect = [content, b""]
    resp.status = 200
    return resp


class TestModuleLevelHelpers:
    def test_compare_versions_matches_the_method(self):
        inst = object.__new__(PI.PluginInstallerWindow)
        for a, b in [("2.10", "2.5"), ("1.0", "1.0.0"), ("1", "2"), ("x", "y")]:
            assert PI.compare_versions(a, b) == inst.compare_versions(a, b)
        assert PI.compare_versions("2026.11.22", "2026.10.30") == 1

    def test_compare_versions_falls_back_to_string_order(self):
        with patch.object(PI.re, "findall", side_effect=TypeError("boom")):
            assert PI.compare_versions("b", "a") == 1
            assert PI.compare_versions("a", "b") == -1

    def test_method_wrappers_delegate(self, monkeypatch, tmp_path):
        inst = object.__new__(PI.PluginInstallerWindow)
        monkeypatch.setattr(PI, "_resolve_app_package", lambda: "moleditpy_linux")
        assert inst._get_package_name() == "moleditpy-linux"
        f = tmp_path / "x.bin"
        f.write_bytes(b"abc")
        assert inst.calculate_sha256(str(f)) == _sha(b"abc")

    @pytest.mark.parametrize(
        ("url", "expected"),
        [
            ("https://example.com/p.py", "https://example.com/p.py"),
            (
                "../plugins/A/a.py",
                "https://hiroyokoyama.github.io/moleditpy-plugins/plugins/A/a.py",
            ),
        ],
    )
    def test_resolve_download_url(self, url, expected):
        assert PI._resolve_download_url(url) == expected

    def test_download_filename_keeps_installed_single_file_name(self, tmp_path):
        existing = tmp_path / "my_local_name.py"
        existing.write_text("x")
        url = "https://example.com/plugins/remote_name.py"
        assert PI._download_filename(url, str(existing)) == "my_local_name.py"

    def test_download_filename_uses_url_for_zip_and_fresh_installs(self, tmp_path):
        existing = tmp_path / "old.py"
        existing.write_text("x")
        assert PI._download_filename("https://e.com/p.zip", str(existing)) == "p.zip"
        assert PI._download_filename("https://e.com/new.py", None) == "new.py"
        init = tmp_path / "pkg" / "__init__.py"
        init.parent.mkdir()
        init.write_text("")
        assert PI._download_filename("https://e.com/new.py", str(init)) == "new.py"

    def test_download_filename_defaults_when_url_has_no_name(self):
        assert PI._download_filename("https://e.com/", None) == "plugin_update.py"

    def test_fetch_registry_retries_then_gives_up(self):
        import urllib.error

        with patch(
            "urllib.request.urlopen", side_effect=urllib.error.URLError("down")
        ) as op, patch("time.sleep") as sleep:
            assert PI._fetch_registry("https://example.com/r.json") == []
        assert op.call_count == PI._FETCH_RETRIES
        assert sleep.call_count == PI._FETCH_RETRIES - 1

    def test_fetch_registry_reads_json_with_bom(self):
        body = "﻿[{\"name\": \"A\"}]".encode("utf-8")
        with patch("urllib.request.urlopen", return_value=_fake_response(body)):
            assert PI._fetch_registry() == [{"name": "A"}]

    def test_fetch_pypi_version(self):
        body = json.dumps({"info": {"version": "9.9.9"}}).encode()
        with patch("urllib.request.urlopen", return_value=_fake_response(body)):
            assert PI._fetch_pypi_version("moleditpy") == "9.9.9"
        with patch("urllib.request.urlopen", side_effect=OSError("x")), patch(
            "time.sleep"
        ):
            assert PI._fetch_pypi_version("moleditpy") == "Unknown"

    def test_worker_uses_the_shared_fetchers(self, monkeypatch):
        monkeypatch.setattr(PI, "_fetch_pypi_version", lambda pkg: f"v-{pkg}")
        monkeypatch.setattr(PI, "_fetch_registry", lambda url: [url])
        worker = PI._FetchWorker("pkg", "https://r")
        assert worker._fetch_pypi() == "v-pkg"
        assert worker._fetch_remote() == ["https://r"]


class TestDownloadFile:
    def test_writes_the_file(self, tmp_path):
        dest = tmp_path / "a.py"
        with patch("urllib.request.urlopen", return_value=_fake_response(b"hello")):
            assert PI._download_file("https://e.com/a.py", str(dest)) is True
        assert dest.read_bytes() == b"hello"

    def test_abort_from_progress_returns_false(self, tmp_path):
        with patch("urllib.request.urlopen", return_value=_fake_response(b"hi")):
            assert (
                PI._download_file(
                    "https://e.com/a.py", str(tmp_path / "a"), lambda r, t: False
                )
                is False
            )

    def test_oversized_download_raises_without_retry(self, tmp_path):
        resp = _fake_response(b"x")
        resp.headers = {"Content-Length": str(PI.MAX_DOWNLOAD_BYTES + 1)}
        with patch("urllib.request.urlopen", return_value=resp) as op:
            with pytest.raises(PI.DownloadError, match="100 MiB"):
                PI._download_file("https://e.com/a.py", str(tmp_path / "a"))
        assert op.call_count == 1

    def test_network_failure_raises_after_retries(self, tmp_path):
        import urllib.error

        with patch(
            "urllib.request.urlopen", side_effect=urllib.error.URLError("down")
        ), patch("time.sleep"):
            with pytest.raises(PI.DownloadError, match="after 3 attempts"):
                PI._download_file("https://e.com/a.py", str(tmp_path / "a"))

    def test_incomplete_download_is_retried(self, tmp_path):
        short = _fake_response(b"abc")
        short.headers = {"Content-Length": "10"}
        with patch(
            "urllib.request.urlopen", side_effect=[short, _fake_response(b"abc")]
        ), patch("time.sleep"):
            assert PI._download_file("https://e.com/a.py", str(tmp_path / "a"))


def _zip_bytes(files, top="pkg"):
    buf = _io.BytesIO()
    with _zipfile.ZipFile(buf, "w") as z:
        for rel, data in files.items():
            z.writestr(f"{top}/{rel}" if top else rel, data)
    return buf.getvalue()


class TestInstallDownload:
    def test_single_file_is_overwritten_in_place(self, tmp_path):
        target = tmp_path / "plugins" / "tool.py"
        target.parent.mkdir()
        target.write_bytes(b"old")
        download = tmp_path / "dl" / "tool.py"
        download.parent.mkdir()
        download.write_bytes(b"new")
        host = MagicMock()
        PI._install_download(str(download), "tool.py", str(target), str(tmp_path), host)
        assert target.read_bytes() == b"new"
        host.assert_not_called()

    def test_single_file_falls_back_to_host_when_copy_fails(self, tmp_path):
        target = tmp_path / "tool.py"
        target.write_bytes(b"old")
        download = tmp_path / "dl.py"
        download.write_bytes(b"new")
        host = MagicMock(return_value=(True, "ok"))
        with patch.object(PI.shutil, "copy2", side_effect=OSError("locked")):
            PI._install_download(
                str(download), "tool.py", str(target), str(tmp_path), host
            )
        host.assert_called_once_with(str(download))

    def test_folder_plugin_is_replaced_and_keeps_settings(self, tmp_path):
        pkg = tmp_path / "plugins" / "pkg"
        pkg.mkdir(parents=True)
        (pkg / "__init__.py").write_text("old")
        (pkg / "stale.py").write_text("gone")
        (pkg / "settings.json").write_text('{"keep": 1}')
        download = tmp_path / "pkg.zip"
        download.write_bytes(_zip_bytes({"__init__.py": "new", "extra.py": "x"}))
        work = tmp_path / "work"
        work.mkdir()
        PI._install_download(
            str(download), "pkg.zip", str(pkg / "__init__.py"), str(work), MagicMock()
        )
        assert (pkg / "__init__.py").read_text() == "new"
        assert (pkg / "extra.py").exists()
        assert not (pkg / "stale.py").exists()
        assert (pkg / "settings.json").read_text() == '{"keep": 1}'

    def test_folder_plugin_zip_without_top_folder(self, tmp_path):
        pkg = tmp_path / "pkg"
        pkg.mkdir()
        (pkg / "__init__.py").write_text("old")
        download = tmp_path / "pkg.zip"
        download.write_bytes(_zip_bytes({"__init__.py": "new"}, top=None))
        work = tmp_path / "work"
        work.mkdir()
        overwrite = MagicMock()
        PI._install_download(
            str(download), "pkg.zip", str(pkg / "__init__.py"), str(work), MagicMock(),
            overwrite_folder=overwrite,
        )
        source, target = overwrite.call_args.args
        assert target == str(pkg)
        assert os.path.basename(source) == "extracted"

    def test_unsafe_zip_entry_is_rejected(self, tmp_path):
        pkg = tmp_path / "pkg"
        pkg.mkdir()
        (pkg / "__init__.py").write_text("old")
        download = tmp_path / "evil.zip"
        download.write_bytes(_zip_bytes({"../../evil.py": "x"}, top=None))
        work = tmp_path / "work"
        work.mkdir()
        with pytest.raises(ValueError, match="Unsafe zip entry"):
            PI._install_download(
                str(download), "evil.zip", str(pkg / "__init__.py"), str(work),
                MagicMock(),
            )
        assert (pkg / "__init__.py").read_text() == "old"

    def test_fresh_install_goes_through_the_host(self, tmp_path):
        download = tmp_path / "new.py"
        download.write_bytes(b"x")
        host = MagicMock(return_value=(True, "Installed"))
        pumped = []
        PI._install_download(
            str(download), "new.py", None, str(tmp_path), host,
            pump=lambda: pumped.append(1),
        )
        host.assert_called_once_with(str(download))
        assert pumped

    def test_host_failure_raises(self, tmp_path):
        download = tmp_path / "new.py"
        download.write_bytes(b"x")
        with pytest.raises(RuntimeError, match="nope"):
            PI._install_download(
                str(download), "new.py", None, str(tmp_path),
                MagicMock(return_value=(False, "nope")),
            )
        with pytest.raises(RuntimeError, match="Host plugin installation failed"):
            PI._install_download(
                str(download), "new.py", None, str(tmp_path),
                MagicMock(return_value=(False, "")),
            )

    def test_single_file_turning_into_a_package_removes_the_old_file(self, tmp_path):
        old = tmp_path / "tool.py"
        old.write_text("old")
        download = tmp_path / "tool.zip"
        download.write_bytes(b"zip")
        PI._install_download(
            str(download), "tool.zip", str(old), str(tmp_path),
            MagicMock(return_value=(True, "ok")),
        )
        assert not old.exists()

    def test_package_turning_into_a_single_file_removes_the_old_folder(
        self, tmp_path
    ):
        pkg = tmp_path / "tool"
        pkg.mkdir()
        (pkg / "__init__.py").write_text("old")
        (pkg / "settings.json").write_text("{}")
        download = tmp_path / "tool.py"
        download.write_bytes(b"new")
        PI._install_download(
            str(download), "tool.py", str(pkg / "__init__.py"), str(tmp_path),
            MagicMock(return_value=(True, "ok")),
        )
        assert not pkg.exists()


# ===========================================================================
# Status classification
# ===========================================================================


def _remote(name="Alpha Tool", version="2.0", **extra):
    entry = {
        "id": name.lower().replace(" ", "_"),
        "name": name,
        "version": version,
        "author": "Someone",
        "downloadUrl": f"https://example.com/{name.lower().replace(' ', '_')}.py",
        "sha256": "0" * 64,
    }
    entry.update(extra)
    return entry


def _local(path, name="Alpha Tool", **extra):
    info = {"name": name, "version": "Unknown", "filepath": str(path)}
    info.update(extra)
    return info


class TestClassifyPlugin:
    def _installed(self, tmp_path, version):
        path = tmp_path / "alpha.py"
        path.write_bytes(_plugin_source(version))
        return path

    @pytest.mark.parametrize(
        ("local", "remote", "status"),
        [
            ("1.0", "2.0", "Update Available"),
            ("2.0", "2.0", "Up to date"),
            ("3.0", "2.0", "Newer"),
            ("2.0", "Unknown", "Version Mismatch"),
        ],
    )
    def test_installed_statuses(self, tmp_path, local, remote, status):
        path = self._installed(tmp_path, local)
        row = PI._classify_plugin(
            "Alpha Tool", _remote(version=remote), _local(path), "5.3.2", []
        )
        assert row["status"] == status
        assert row["local_ver"] == local
        assert row["can_update"] is True

    def test_version_comes_from_the_file_not_the_cached_info(self, tmp_path):
        path = self._installed(tmp_path, "1.5")
        row = PI._classify_plugin(
            "Alpha Tool", _remote(), _local(path, version="9.9"), "5.3.2", []
        )
        assert row["local_ver"] == "1.5"

    def test_package_version_comes_from_local_info(self, tmp_path):
        row = PI._classify_plugin(
            "Alpha Tool", _remote(), _local(tmp_path / "missing", version="1.0"),
            "5.3.2", [],
        )
        assert row["local_ver"] == "1.0"
        assert row["status"] == "Update Available"

    def test_not_in_registry(self, tmp_path):
        path = self._installed(tmp_path, "1.0")
        row = PI._classify_plugin("Alpha Tool", None, _local(path), "5.3.2", [])
        assert row["status"] == "Not in Registry"
        assert row["can_update"] is False

    def test_stale_entry_with_migration(self, tmp_path):
        path = self._installed(tmp_path, "1.0")
        registry = [{"id": "beta", "name": "Beta Tool"}]
        row = PI._classify_plugin(
            "Alpha Tool",
            _remote(visible=False, migration_to="beta"),
            _local(path),
            "5.3.2",
            registry,
        )
        assert row["status"] == "Stale"
        assert row["migration_target"] == "Beta Tool"

    def test_hidden_entry_that_is_not_installed_is_skipped(self):
        assert (
            PI._classify_plugin("Alpha Tool", _remote(visible=False), None, "5.3.2", [])
            is None
        )

    def test_not_installed_and_incompatible(self):
        row = PI._classify_plugin("Alpha Tool", _remote(), None, "5.3.2", [])
        assert row["status"] == "Not Installed"
        assert row["can_download"] is True
        assert row["target_file"].endswith("alpha_tool.py")
        row = PI._classify_plugin(
            "Alpha Tool",
            _remote(supported_moleditpy_version=">=6.0"),
            None,
            "5.3.2",
            [],
        )
        assert row["status"] == "Incompatible"
        assert row["is_compatible"] is False

    def test_python_and_os_requirements(self, monkeypatch):
        monkeypatch.setattr(PI, "get_current_os", lambda: "Windows")
        row = PI._classify_plugin(
            "Alpha Tool",
            _remote(supported_python_version=">=99", supported_os=["Linux"]),
            None,
            "5.3.2",
            [],
        )
        assert row["py_compatible"] is False
        assert row["os_compatible"] is False
        assert row["is_compatible"] is False

    @pytest.mark.parametrize(
        ("local", "visible", "text"),
        [
            ("1.0", True, "Update Available (Disabled)"),
            ("2.0", True, "Disabled"),
            ("1.0", False, "Stale (Disabled)"),
        ],
    )
    def test_disabled_status_text(self, tmp_path, local, visible, text):
        path = self._installed(tmp_path, local)
        row = PI._classify_plugin(
            "Alpha Tool",
            _remote(visible=visible),
            _local(path, disabled=True),
            "5.3.2",
            [],
        )
        assert row["status_text"] == text

    def test_author_falls_back_to_local_info(self, tmp_path):
        path = self._installed(tmp_path, "1.0")
        row = PI._classify_plugin(
            "Alpha Tool", None, _local(path, author="Me"), "5.3.2", []
        )
        assert row["author"] == "Me"


class TestRowColor:
    def _row(self, status, compatible=True, disabled=False, stale=False):
        return {
            "status": status,
            "is_compatible": compatible,
            "is_disabled": disabled,
            "is_stale": stale,
        }

    @pytest.mark.parametrize(
        ("status", "color"),
        [
            ("Update Available", "#f8d7da"),
            ("Up to date", "#d4edda"),
            ("Newer", "#cce5ff"),
            ("Not in Registry", "#e2e3e5"),
            ("Stale", "#fff3cd"),
            ("Version Mismatch", None),
            ("Not Installed", None),
        ],
    )
    def test_status_colors(self, status, color):
        assert PI._row_color(self._row(status)) == color

    def test_incompatible_and_disabled_override(self):
        assert PI._row_color(self._row("Up to date", compatible=False)) == "#fff3cd"
        assert PI._row_color(self._row("Up to date", disabled=True)) == "#d6d8db"
        assert (
            PI._row_color(self._row("Update Available", disabled=True)) == "#f8d7da"
        )
        assert PI._row_color(self._row("Stale", disabled=True, stale=True)) == "#fff3cd"


# ===========================================================================
# Command line
# ===========================================================================


class _Cli:
    """A plugin folder, a registry and a fake host, wired into the CLI."""

    def __init__(self, tmp_path, monkeypatch):
        self.dir = tmp_path / "plugins"
        self.dir.mkdir()
        self.registry = []
        self.installed = {}
        self.downloads = {}
        self.monkeypatch = monkeypatch
        self.app_version = "5.3.2"
        self.pypi_version = "5.3.2"
        monkeypatch.setattr(PI, "_fetch_registry", lambda url=None: self.registry)
        monkeypatch.setattr(PI, "_installed_plugin_infos", lambda: self.installed)
        monkeypatch.setattr(PI, "get_app_version", lambda: self.app_version)
        monkeypatch.setattr(PI, "get_package_name", lambda: "moleditpy")
        monkeypatch.setattr(PI, "_fetch_pypi_version", lambda pkg: self.pypi_version)
        monkeypatch.setattr(PI, "_download_file", self._download)
        monkeypatch.setattr(PI, "_host_install_plugin", self._host_install)
        self.host_installs = []

    def _download(self, url, dest, on_progress=None):
        if url not in self.downloads:
            raise PI.DownloadError(f"Download failed:\n404 for {url}")
        with open(dest, "wb") as f:
            f.write(self.downloads[url])
        return True

    def _host_install(self, path):
        self.host_installs.append(path)
        dest = self.dir / os.path.basename(path)
        dest.write_bytes(Path(path).read_bytes())
        return True, f"Installed {dest.name}"

    def publish(self, name, version, sha=True, **extra):
        data = _plugin_source(version, name)
        entry = _remote(name=name, version=version, **extra)
        entry["sha256"] = _sha(data) if sha is True else sha
        if sha is None:
            del entry["sha256"]
        self.registry.append(entry)
        self.downloads[entry["downloadUrl"]] = data
        return entry

    def install(self, name, version, **extra):
        path = self.dir / (name.lower().replace(" ", "_") + ".py")
        path.write_bytes(_plugin_source(version, name))
        self.installed[name] = _local(path, name=name, version=version, **extra)
        return path

    def run(self, *argv):
        out, err = _io.StringIO(), _io.StringIO()
        code = PI.cli_main(list(argv), out, err)
        return code, out.getvalue(), err.getvalue()


@pytest.fixture
def cli(tmp_path, monkeypatch):
    return _Cli(tmp_path, monkeypatch)


def _installed_version(path):
    return PI._read_plugin_version_ast(str(path))


class TestCliBasics:
    def test_no_command_prints_help(self, cli):
        code, out, _ = cli.run()
        assert code == PI.EXIT_USAGE
        assert "list" in out and "update" in out and "install" in out

    def test_help_exits_zero(self, cli, capsys):
        assert PI.cli_main(["--help"]) == 0
        assert "plugin_installer" in capsys.readouterr().out

    def test_unknown_command_is_a_usage_error(self, cli):
        assert PI.cli_main(["frobnicate"]) == PI.EXIT_USAGE

    def test_unavailable_registry(self, cli):
        cli.registry = []
        code, _, err = cli.run("list")
        assert code == PI.EXIT_UNAVAILABLE
        assert "Could not fetch the plugin registry" in err

    def test_non_dict_registry_entries_are_ignored(self, cli):
        cli.publish("Alpha Tool", "2.0")
        cli.registry.append("garbage")
        code, out, _ = cli.run("list", "--all")
        assert code == 0
        assert "Alpha Tool" in out

    def test_keyboard_interrupt(self, cli, monkeypatch):
        def boom(args, out, err):
            raise KeyboardInterrupt

        monkeypatch.setattr(PI, "_cmd_list", boom)
        code, _, err = cli.run("list")
        assert code == PI.EXIT_FAILED
        assert "Interrupted" in err

    def test_main_guard_is_top_level(self):
        import ast

        tree = ast.parse(PLUGIN_PATH.read_text(encoding="utf-8"))
        guards = [
            node
            for node in tree.body
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.Compare)
            and getattr(node.test.left, "id", None) == "__name__"
        ]
        assert len(guards) == 1

    def test_the_main_apps_runner_sees_the_installer_as_runnable(self):
        src = (
            Path(__file__).resolve().parents[2]
            / "python_molecular_editor"
            / "moleditpy"
            / "src"
        )
        runner = src / "moleditpy" / "plugins" / "plugin_runner.py"
        if not runner.exists() or "def has_main_guard" not in runner.read_text(
            encoding="utf-8"
        ):
            pytest.skip("main app with single-file --run-plugin support not present")
        import ast

        ns = {"ast": ast, "logging": __import__("logging")}
        tree = ast.parse(runner.read_text(encoding="utf-8"))
        wanted = [
            n
            for n in tree.body
            if isinstance(n, ast.FunctionDef)
            and n.name in ("_is_main_guard", "has_main_guard")
        ]
        exec(  # noqa: S102 - two pure functions from the main app's source
            compile(ast.Module(body=wanted, type_ignores=[]), str(runner), "exec"),
            ns,
        )
        assert ns["has_main_guard"](str(PLUGIN_PATH)) is True


class TestCliInstalledPluginInfos:
    def test_maps_the_hosts_plugin_list(self, monkeypatch, tmp_path):
        from dataclasses import dataclass

        @dataclass
        class Fake:
            name: str
            path: str
            display_name: str
            version: str
            package: bool
            enabled: bool

        plugins = [
            Fake("a", str(tmp_path / "a.py"), "Alpha", "1.0", False, True),
            Fake("b", str(tmp_path / "b"), "Beta", "2.0", True, False),
            Fake("a2", str(tmp_path / "x" / "a.py"), "Alpha", "0.1", False, True),
        ]
        runner = types.ModuleType("moleditpy.plugins.plugin_runner")
        runner.installed_plugins = lambda: plugins
        monkeypatch.setitem(sys.modules, "moleditpy.plugins.plugin_runner", runner)
        infos = PI._installed_plugin_infos()
        assert infos["Alpha"]["filepath"] == str(tmp_path / "a.py")
        assert infos["Alpha"]["version"] == "1.0"
        assert infos["Beta"]["filepath"] == os.path.join(str(tmp_path / "b"), "__init__.py")
        assert infos["Beta"]["disabled"] is True

    def test_old_host_is_reported(self, monkeypatch):
        monkeypatch.setitem(sys.modules, "moleditpy.plugins.plugin_runner", None)
        with pytest.raises(PI._CliError, match="5.3.2 or later") as info:
            PI._installed_plugin_infos()
        assert info.value.code == PI.EXIT_UNAVAILABLE

    def test_old_host_exit_status(self, tmp_path, monkeypatch):
        monkeypatch.setitem(sys.modules, "moleditpy.plugins.plugin_runner", None)
        monkeypatch.setattr(PI, "_fetch_registry", lambda url=None: [_remote()])
        err = _io.StringIO()
        assert PI.cli_main(["list"], _io.StringIO(), err) == PI.EXIT_UNAVAILABLE
        assert "5.3.2" in err.getvalue()

    def test_host_install_uses_the_plugin_manager(self, monkeypatch):
        manager = MagicMock()
        manager.return_value.install_plugin.return_value = (True, "ok")
        pm = types.ModuleType("moleditpy.plugins.plugin_manager")
        pm.PluginManager = manager
        monkeypatch.setitem(sys.modules, "moleditpy.plugins.plugin_manager", pm)
        assert PI._host_install_plugin("/tmp/x.py") == (True, "ok")
        manager.return_value.install_plugin.assert_called_once_with("/tmp/x.py")


class TestCliList:
    def test_lists_installed_plugins_only_by_default(self, cli):
        cli.publish("Alpha Tool", "2.0")
        cli.publish("Beta Tool", "1.0")
        cli.install("Alpha Tool", "1.0")
        code, out, _ = cli.run("list")
        assert code == 0
        assert "Alpha Tool" in out
        assert "Update Available" in out
        assert "Beta Tool" not in out

    def test_all_includes_plugins_that_are_not_installed(self, cli):
        cli.publish("Beta Tool", "1.0")
        code, out, _ = cli.run("list", "--all")
        assert "Beta Tool" in out and "Not Installed" in out
        assert "beta_tool" in out

    def test_empty(self, cli):
        cli.publish("Beta Tool", "1.0")
        code, out, _ = cli.run("list")
        assert code == 0
        assert "No plugins to show." in out

    def test_json(self, cli):
        cli.publish("Alpha Tool", "2.0")
        path = cli.install("Alpha Tool", "1.0", disabled=True)
        code, out, _ = cli.run("list", "--json")
        data = json.loads(out)
        assert data == [
            {
                "name": "Alpha Tool",
                "id": "alpha_tool",
                "installed_version": "1.0",
                "latest_version": "2.0",
                "status": "Update Available",
                "disabled": True,
                "compatible": True,
                "path": str(path),
            }
        ]

    def test_json_for_a_plugin_outside_the_registry(self, cli):
        cli.publish("Beta Tool", "1.0")
        cli.install("Local Only", "0.1")
        data = json.loads(cli.run("list", "--json")[1])
        assert data[0]["id"] is None
        assert data[0]["latest_version"] is None
        assert data[0]["status"] == "Not in Registry"


class TestCliCheck:
    def test_up_to_date(self, cli):
        cli.publish("Alpha Tool", "1.0")
        cli.install("Alpha Tool", "1.0")
        code, out, _ = cli.run("check")
        assert code == 0
        assert "MoleditPy 5.3.2 is up to date." in out
        assert "All installed plugins are up to date." in out

    def test_updates_and_app_update(self, cli):
        cli.publish("Alpha Tool", "2.0")
        cli.install("Alpha Tool", "1.0")
        cli.pypi_version = "5.4.0"
        code, out, _ = cli.run("check")
        assert code == 0
        assert "MoleditPy 5.4.0 is available" in out
        assert "pip install --upgrade moleditpy" in out
        assert "1 plugin update(s) available" in out
        assert "update --all" in out

    def test_exit_code_flag(self, cli):
        cli.publish("Alpha Tool", "2.0")
        cli.install("Alpha Tool", "1.0")
        assert cli.run("check", "--exit-code")[0] == PI.EXIT_FAILED
        cli.installed.clear()
        cli.install("Alpha Tool", "2.0")
        assert cli.run("check", "--exit-code")[0] == 0
        cli.pypi_version = "9.0"
        assert cli.run("check", "--exit-code")[0] == PI.EXIT_FAILED

    def test_pypi_unreachable(self, cli):
        cli.publish("Alpha Tool", "1.0")
        cli.pypi_version = "Unknown"
        out = cli.run("check")[1]
        assert "could not check PyPI" in out

    def test_stale_plugins_are_reported_with_migration(self, cli):
        cli.publish("Beta Tool", "1.0")
        cli.publish("Alpha Tool", "1.0", visible=False, migration_to="beta_tool")
        cli.install("Alpha Tool", "1.0")
        out = cli.run("check")[1]
        assert "Alpha Tool is no longer offered by the registry" in out
        assert "moved to Beta Tool" in out

    def test_json(self, cli):
        cli.publish("Alpha Tool", "2.0")
        cli.publish("Old Tool", "1.0", visible=False)
        cli.install("Alpha Tool", "1.0")
        cli.install("Old Tool", "1.0")
        cli.pypi_version = "Unknown"
        data = json.loads(cli.run("check", "--json")[1])
        assert data["app"] == {
            "package": "moleditpy",
            "installed_version": "5.3.2",
            "latest_version": None,
            "update_available": False,
        }
        assert [u["name"] for u in data["updates"]] == ["Alpha Tool"]
        assert [s["name"] for s in data["stale"]] == ["Old Tool"]


class TestCliUpdate:
    def test_needs_names_or_all(self, cli):
        code, _, err = cli.run("update")
        assert code == PI.EXIT_USAGE
        assert "--all" in err
        assert cli.run("update", "x", "--all")[0] == PI.EXIT_USAGE

    def test_update_all(self, cli):
        cli.publish("Alpha Tool", "2.0")
        cli.publish("Beta Tool", "1.0")
        alpha = cli.install("Alpha Tool", "1.0")
        beta = cli.install("Beta Tool", "1.0")
        code, out, err = cli.run("update", "--all", "--yes")
        assert (code, err) == (0, "")
        assert _installed_version(alpha) == "2.0"
        assert _installed_version(beta) == "1.0"
        assert "Alpha Tool  1.0 -> 2.0" in out
        assert "Alpha Tool: updated to 2.0." in out
        assert "Restart MoleditPy" in out
        assert cli.host_installs == []  # single files are overwritten in place

    def test_update_all_when_nothing_to_do(self, cli):
        cli.publish("Alpha Tool", "1.0")
        cli.install("Alpha Tool", "1.0")
        code, out, _ = cli.run("update", "--all")
        assert code == 0
        assert "All installed plugins are up to date." in out

    def test_update_all_skips_stale_and_includes_disabled(self, cli):
        cli.publish("Alpha Tool", "2.0", visible=False)
        cli.publish("Beta Tool", "2.0")
        alpha = cli.install("Alpha Tool", "1.0")
        beta = cli.install("Beta Tool", "1.0", disabled=True)
        assert cli.run("update", "--all", "--yes")[0] == 0
        assert _installed_version(alpha) == "1.0"
        assert _installed_version(beta) == "2.0"

    def test_dry_run_changes_nothing(self, cli):
        cli.publish("Alpha Tool", "2.0")
        alpha = cli.install("Alpha Tool", "1.0")
        code, out, _ = cli.run("update", "--all", "--dry-run")
        assert code == 0
        assert "Dry run" in out
        assert _installed_version(alpha) == "1.0"

    def test_refuses_to_prompt_without_a_terminal(self, cli, monkeypatch):
        cli.publish("Alpha Tool", "2.0")
        alpha = cli.install("Alpha Tool", "1.0")
        monkeypatch.setattr(sys, "stdin", _io.StringIO("y\n"))
        code, _, err = cli.run("update", "--all")
        assert code == PI.EXIT_USAGE
        assert "--yes" in err
        assert _installed_version(alpha) == "1.0"

    @pytest.mark.parametrize(("answer", "updated"), [("y\n", True), ("n\n", False)])
    def test_prompt_on_a_terminal(self, cli, monkeypatch, answer, updated):
        cli.publish("Alpha Tool", "2.0")
        alpha = cli.install("Alpha Tool", "1.0")
        stdin = _io.StringIO(answer)
        stdin.isatty = lambda: True
        monkeypatch.setattr(sys, "stdin", stdin)
        code, out, _ = cli.run("update", "--all")
        assert "Proceed? [y/N]" in out
        assert _installed_version(alpha) == ("2.0" if updated else "1.0")
        assert code == (0 if updated else PI.EXIT_FAILED)
        if not updated:
            assert "Cancelled." in out

    def test_prompt_without_stdin(self, cli, monkeypatch):
        cli.publish("Alpha Tool", "2.0")
        cli.install("Alpha Tool", "1.0")
        monkeypatch.setattr(sys, "stdin", None)
        assert cli.run("update", "--all")[0] == PI.EXIT_USAGE

    def test_checksum_mismatch_stops_the_update(self, cli):
        cli.publish("Alpha Tool", "2.0", sha="f" * 64)
        alpha = cli.install("Alpha Tool", "1.0")
        code, _, err = cli.run("update", "--all", "-y")
        assert code == PI.EXIT_FAILED
        assert "SHA256 checksum mismatch" in err
        assert "f" * 64 in err
        assert _installed_version(alpha) == "1.0"

    def test_missing_checksum_needs_allow_unverified(self, cli):
        cli.publish("Alpha Tool", "2.0", sha=None)
        alpha = cli.install("Alpha Tool", "1.0")
        code, _, err = cli.run("update", "--all", "-y")
        assert code == PI.EXIT_FAILED
        assert "--allow-unverified" in err
        assert _installed_version(alpha) == "1.0"
        assert cli.run("update", "--all", "-y", "--allow-unverified")[0] == 0
        assert _installed_version(alpha) == "2.0"

    def test_incompatible_plugin_needs_force(self, cli):
        cli.publish("Alpha Tool", "2.0", supported_moleditpy_version=">=6.0")
        alpha = cli.install("Alpha Tool", "1.0")
        code, _, err = cli.run("update", "--all", "-y")
        assert code == PI.EXIT_FAILED
        assert "requires MoleditPy >=6.0 (running 5.3.2)" in err
        assert "--force" in err
        assert _installed_version(alpha) == "1.0"
        assert cli.run("update", "--all", "-y", "--force")[0] == 0
        assert _installed_version(alpha) == "2.0"

    def test_incompatible_python_and_os_are_explained(self, cli, monkeypatch):
        monkeypatch.setattr(PI, "get_current_os", lambda: "Windows")
        cli.publish(
            "Alpha Tool",
            "2.0",
            supported_python_version=">=99",
            supported_os=["Linux"],
        )
        cli.install("Alpha Tool", "1.0")
        err = cli.run("update", "--all", "-y")[2]
        assert "Python >=99" in err
        assert "OS Linux (running Windows)" in err

    def test_one_failure_does_not_stop_the_rest(self, cli):
        cli.publish("Alpha Tool", "2.0", sha="0" * 64)
        cli.publish("Beta Tool", "2.0")
        alpha = cli.install("Alpha Tool", "1.0")
        beta = cli.install("Beta Tool", "1.0")
        code, out, _ = cli.run("update", "--all", "-y")
        assert code == PI.EXIT_FAILED
        assert _installed_version(alpha) == "1.0"
        assert _installed_version(beta) == "2.0"
        assert "Restart MoleditPy" in out

    def test_named_update_by_id_or_any_case_name(self, cli):
        cli.publish("Alpha Tool", "2.0")
        cli.publish("Beta Tool", "2.0")
        alpha = cli.install("Alpha Tool", "1.0")
        beta = cli.install("Beta Tool", "1.0")
        assert cli.run("update", "ALPHA_TOOL", "-y")[0] == 0
        assert _installed_version(alpha) == "2.0"
        assert _installed_version(beta) == "1.0"
        assert cli.run("update", "beta tool", "-y")[0] == 0
        assert _installed_version(beta) == "2.0"

    def test_named_up_to_date_plugin_needs_reinstall(self, cli):
        cli.publish("Alpha Tool", "2.0")
        alpha = cli.install("Alpha Tool", "2.0")
        alpha.write_bytes(_plugin_source("2.0") + b"# edited locally\n")
        code, out, _ = cli.run("update", "alpha_tool", "-y")
        assert code == 0
        assert "up to date (2.0)" in out
        assert "--reinstall" in out
        assert b"edited locally" in alpha.read_bytes()
        assert cli.run("update", "alpha_tool", "-y", "--reinstall")[0] == 0
        assert b"edited locally" not in alpha.read_bytes()

    def test_named_stale_plugin_can_be_reinstalled(self, cli):
        cli.publish("Alpha Tool", "2.0", visible=False)
        alpha = cli.install("Alpha Tool", "1.0")
        assert cli.run("update", "alpha_tool", "-y")[0] == 0
        assert _installed_version(alpha) == "1.0"
        assert cli.run("update", "alpha_tool", "-y", "--reinstall")[0] == 0
        assert _installed_version(alpha) == "2.0"

    def test_named_plugin_that_is_not_installed(self, cli):
        cli.publish("Alpha Tool", "2.0")
        code, _, err = cli.run("update", "alpha_tool")
        assert code == PI.EXIT_USAGE
        assert "use 'install'" in err

    def test_named_plugin_outside_the_registry(self, cli):
        cli.publish("Beta Tool", "1.0")
        cli.install("Local Only", "1.0")
        code, _, err = cli.run("update", "local only", "-y")
        assert code == 0
        assert "not in the registry" in err

    def test_unknown_name(self, cli):
        cli.publish("Alpha Tool", "2.0")
        code, _, err = cli.run("update", "nothing")
        assert code == PI.EXIT_USAGE
        assert "No plugin named 'nothing'" in err

    def test_blank_name_matches_nothing(self, cli):
        cli.publish("Alpha Tool", "2.0")
        cli.registry[0]["id"] = ""
        assert cli.run("update", " ")[0] == PI.EXIT_USAGE

    def test_download_failure(self, cli):
        entry = cli.publish("Alpha Tool", "2.0")
        del cli.downloads[entry["downloadUrl"]]
        cli.install("Alpha Tool", "1.0")
        code, _, err = cli.run("update", "--all", "-y")
        assert code == PI.EXIT_FAILED
        assert "404" in err

    def test_unsafe_url_is_refused(self, cli):
        cli.publish("Alpha Tool", "2.0", downloadUrl="file:///etc/passwd")
        cli.install("Alpha Tool", "1.0")
        code, _, err = cli.run("update", "--all", "-y")
        assert code == PI.EXIT_FAILED
        assert "unsafe URL" in err

    def test_entry_without_download_url(self, cli):
        cli.publish("Alpha Tool", "2.0")
        del cli.registry[0]["downloadUrl"]
        cli.install("Alpha Tool", "1.0")
        code, _, err = cli.run("update", "--all", "-y")
        assert code == PI.EXIT_FAILED
        assert "no download URL" in err

    def test_install_step_failure_is_reported(self, cli, monkeypatch):
        cli.publish("Alpha Tool", "2.0")
        cli.install("Alpha Tool", "1.0")

        def broken(*a, **k):
            raise RuntimeError("disk full")

        monkeypatch.setattr(PI, "_install_download", broken)
        code, _, err = cli.run("update", "--all", "-y")
        assert code == PI.EXIT_FAILED
        assert "installation failed: disk full" in err

    def test_relative_download_url_is_resolved(self, cli, monkeypatch):
        entry = cli.publish("Alpha Tool", "2.0")
        data = cli.downloads.pop(entry["downloadUrl"])
        entry["downloadUrl"] = "../plugins/Alpha/alpha_tool.py"
        cli.downloads[PI._resolve_download_url(entry["downloadUrl"])] = data
        alpha = cli.install("Alpha Tool", "1.0")
        assert cli.run("update", "--all", "-y")[0] == 0
        assert _installed_version(alpha) == "2.0"


class TestCliInstall:
    def test_installs_through_the_host(self, cli):
        cli.publish("Beta Tool", "1.0")
        code, out, err = cli.run("install", "beta_tool", "-y")
        assert (code, err) == (0, "")
        assert len(cli.host_installs) == 1
        assert _installed_version(cli.dir / "beta_tool.py") == "1.0"
        assert "Beta Tool: installed to 1.0." in out
        assert "not installed -> 1.0" in out

    def test_already_installed(self, cli):
        cli.publish("Alpha Tool", "2.0")
        cli.install("Alpha Tool", "1.0")
        code, out, _ = cli.run("install", "Alpha Tool", "-y")
        assert code == 0
        assert "already installed (1.0)" in out
        assert cli.host_installs == []

    def test_missing_dependencies_are_listed_with_a_pip_command(
        self, cli, monkeypatch
    ):
        cli.publish("Beta Tool", "1.0", dependencies=["numpy>=1.20", "evil; rm"])
        monkeypatch.setattr(PI, "check_dependency_satisfied", lambda d: False)
        code, out, _ = cli.run("install", "beta_tool", "--dry-run")
        assert code == 0
        assert "missing dependencies: numpy>=1.20, evil; rm" in out
        assert 'pip install "numpy>=1.20"' in out
        assert "rm" not in out.split("pip install", 1)[1].splitlines()[0]

    def test_missing_dependency_without_a_safe_pip_command(self, cli, monkeypatch):
        cli.publish("Beta Tool", "1.0", dependencies=["evil; rm"])
        monkeypatch.setattr(PI, "check_dependency_satisfied", lambda d: False)
        out = cli.run("install", "beta_tool", "--dry-run")[1]
        assert "missing dependencies: evil; rm" in out
        assert "pip install" not in out

    def test_satisfied_dependencies_are_not_mentioned(self, cli, monkeypatch):
        cli.publish("Beta Tool", "1.0", dependencies=["numpy"])
        monkeypatch.setattr(PI, "check_dependency_satisfied", lambda d: True)
        out = cli.run("install", "beta_tool", "--dry-run")[1]
        assert "missing dependencies" not in out

    def test_incompatible_install_is_skipped(self, cli):
        cli.publish("Beta Tool", "1.0", supported_moleditpy_version="<5")
        code, _, err = cli.run("install", "beta_tool", "-y")
        assert code == PI.EXIT_FAILED
        assert "use --force to install anyway" in err
        assert cli.host_installs == []

    def test_install_needs_a_name(self, cli, capsys):
        assert cli.run("install")[0] == PI.EXIT_USAGE


class TestCoreEdgeCases:
    def test_stream_larger_than_the_limit_without_a_length(self, tmp_path, monkeypatch):
        monkeypatch.setattr(PI, "MAX_DOWNLOAD_BYTES", 4)
        resp = _fake_response(b"123456", length=False)
        with patch("urllib.request.urlopen", return_value=resp):
            with pytest.raises(PI.DownloadError, match="100 MiB"):
                PI._download_file("https://e.com/a.py", str(tmp_path / "a"))

    def test_failed_cleanup_of_the_old_form_is_only_logged(self, tmp_path, caplog):
        old = tmp_path / "tool.py"
        old.write_text("old")
        download = tmp_path / "tool.zip"
        download.write_bytes(b"zip")
        with patch.object(PI.os, "remove", side_effect=OSError("busy")):
            PI._install_download(
                str(download), "tool.zip", str(old), str(tmp_path),
                MagicMock(return_value=(True, "ok")),
            )
        assert "failed to remove old .py" in caplog.text

        pkg = tmp_path / "pkg"
        pkg.mkdir()
        (pkg / "__init__.py").write_text("")
        single = tmp_path / "pkg.py"
        single.write_text("new")
        with patch.object(PI.shutil, "rmtree", side_effect=OSError("busy")):
            PI._install_download(
                str(single), "pkg.py", str(pkg / "__init__.py"), str(tmp_path),
                MagicMock(return_value=(True, "ok")),
            )
        assert "failed to remove old folder" in caplog.text

    def test_settings_are_restored_after_a_host_install_into_the_folder(
        self, tmp_path
    ):
        pkg = tmp_path / "pkg"
        pkg.mkdir()
        (pkg / "__init__.py").write_text("")
        (pkg / "settings.json").write_text('{"a": 1}')

        def host(path):
            (pkg / "settings.json").write_text("{}")
            return True, "ok"

        download = tmp_path / "pkg.py"
        download.write_text("new")
        with patch.object(PI.shutil, "rmtree"):
            PI._install_download(
                str(download), "pkg.py", str(pkg / "__init__.py"), str(tmp_path), host
            )
        assert (pkg / "settings.json").read_text() == '{"a": 1}'

    def test_a_stdin_read_error_counts_as_no(self, cli, monkeypatch):
        cli.publish("Alpha Tool", "2.0")
        alpha = cli.install("Alpha Tool", "1.0")
        stdin = MagicMock()
        stdin.isatty.return_value = True
        stdin.readline.side_effect = OSError("closed")
        monkeypatch.setattr(sys, "stdin", stdin)
        assert cli.run("update", "--all")[0] == PI.EXIT_FAILED
        assert _installed_version(alpha) == "1.0"
