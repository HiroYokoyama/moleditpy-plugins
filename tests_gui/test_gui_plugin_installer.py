"""
Headless GUI tests for the Plugin Installer plugin.

Covers: PluginDetailsDialog, PluginInstallerWindow, and pure-function unit
tests (is_app_version_compatible, parse_dependency, check_dependency_satisfied).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PLUGIN_INSTALLER_PATH = PLUGINS_DIR / "Plugin_Installer" / "plugin_installer.py"

with mock_chemistry_imports():
    _plugin_installer = load_plugin_for_gui(PLUGIN_INSTALLER_PATH)


# ===========================================================================
# PluginDetailsDialog  (inside Plugin Installer)
# ===========================================================================


def _make_details_dialog(*, dependencies=None, local_info=None) -> object:
    """Create a PluginDetailsDialog with safe defaults."""
    return _plugin_installer.PluginDetailsDialog(
        None,
        "Test Plugin",
        "Test Author",
        "2026.01.01",
        "A headless test plugin.",
        dependencies if dependencies is not None else [],
        local_info,
        None,
        ">=4.0.0",
    )


class TestPluginDetailsDialog:
    """PluginDetailsDialog — lightweight details popup, no network calls."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _make_details_dialog()
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title_includes_plugin_name(self, dlg):
        assert "Test Plugin" in dlg.windowTitle()

    def test_no_error_with_empty_dependencies(self, qapp):
        d = _make_details_dialog(dependencies=[])
        assert d is not None
        d.destroy()

    def test_no_error_with_missing_dependency(self, qapp):
        d = _make_details_dialog(dependencies=["fake-nonexistent-pkg-xyz"])
        assert d is not None
        d.destroy()

    def test_no_error_with_installed_dependency(self, qapp):
        d = _make_details_dialog(dependencies=["pytest"])
        assert d is not None
        d.destroy()

    def test_target_file_stored(self, qapp):
        d = _plugin_installer.PluginDetailsDialog(
            None, "X", "A", "1.0", "Desc", [], None, "/some/path.py"
        )
        assert d.target_file == "/some/path.py"
        d.destroy()


# ===========================================================================
# PluginInstallerWindow  (visible plugin: "Plugin Installer")
# ===========================================================================


class TestPluginInstallerWindow:
    """PluginInstallerWindow — main window of the Plugin Installer.

    check_updates() is patched to prevent any network calls during __init__.
    """

    @pytest.fixture
    def installer(self, qapp):
        from unittest.mock import patch

        with patch.object(
            _plugin_installer.PluginInstallerWindow,
            "check_updates",
            lambda self: None,
        ):
            d = _plugin_installer.PluginInstallerWindow(None, auto_check=False)
        yield d
        d.destroy()

    def test_creates_without_error(self, installer):
        assert installer is not None

    def test_window_title(self, installer):
        assert installer.windowTitle() == "Plugin Installer"

    def test_table_has_six_columns(self, installer):
        assert installer.table.columnCount() == 6

    def test_search_input_exists(self, installer):
        assert installer.search_input.placeholderText() != ""

    def test_startup_checkbox_exists(self, installer):
        assert installer.chk_startup is not None


# ===========================================================================
# Plugin Installer — pure-function unit tests (no Qt needed)
# ===========================================================================


class TestPluginInstallerVersionCompat:
    """is_app_version_compatible(): pure Python, no Qt required."""

    def test_within_range(self):
        assert _plugin_installer.is_app_version_compatible("4.2.0", ">=4.0.0,<5.0.0")

    def test_equal_lower_bound(self):
        assert _plugin_installer.is_app_version_compatible("4.0.0", ">=4.0.0")

    def test_below_lower_bound(self):
        assert not _plugin_installer.is_app_version_compatible("3.9.9", ">=4.0.0")

    def test_equal_upper_bound_exclusive(self):
        assert not _plugin_installer.is_app_version_compatible("5.0.0", "<5.0.0")

    def test_wildcard_always_passes(self):
        assert _plugin_installer.is_app_version_compatible("1.0.0", "*")

    def test_empty_specifier_always_passes(self):
        assert _plugin_installer.is_app_version_compatible("1.0.0", "")

    def test_not_equal(self):
        assert not _plugin_installer.is_app_version_compatible("4.0.0", "!=4.0.0")

    def test_compatible_release_same_major(self):
        assert _plugin_installer.is_app_version_compatible("4.2.0", "~=4.0")

    def test_compatible_release_different_major(self):
        assert not _plugin_installer.is_app_version_compatible("5.0.0", "~=4.0")


class TestPluginInstallerParseDependency:
    """parse_dependency(): pure Python, no Qt required."""

    def test_name_only(self):
        name, spec = _plugin_installer.parse_dependency("numpy")
        assert name == "numpy"
        assert spec == ""

    def test_name_with_specifier(self):
        name, spec = _plugin_installer.parse_dependency("rdkit>=2022.03.5")
        assert name == "rdkit"
        assert "2022" in spec

    def test_empty_string_returns_empty_pair(self):
        name, spec = _plugin_installer.parse_dependency("")
        assert name == ""
        assert spec == ""

    def test_url_with_colon_is_rejected(self):
        name, spec = _plugin_installer.parse_dependency("http://example.com/package")
        assert name == ""

    def test_dash_in_name(self):
        name, spec = _plugin_installer.parse_dependency("my-package>=1.0")
        assert name == "my-package"


class TestPluginInstallerCheckDependency:
    """check_dependency_satisfied(): queries importlib.metadata, pure Python."""

    def test_installed_package_returns_true(self):
        assert _plugin_installer.check_dependency_satisfied("pytest")

    def test_nonexistent_package_returns_false(self):
        assert not _plugin_installer.check_dependency_satisfied(
            "fake-nonexistent-dep-xyz123"
        )

    def test_empty_string_returns_false(self):
        assert not _plugin_installer.check_dependency_satisfied("")


# ===========================================================================
# Plugin Installer — PluginInstallerWindow more widget tests
# ===========================================================================


class TestPluginInstallerWindowWidgets:
    """Additional widget-state tests for PluginInstallerWindow."""

    @pytest.fixture
    def installer(self, qapp):
        from unittest.mock import patch

        with patch.object(
            _plugin_installer.PluginInstallerWindow,
            "check_updates",
            lambda self: None,
        ):
            d = _plugin_installer.PluginInstallerWindow(None, auto_check=False)
        yield d
        d.destroy()

    def test_table_header_first_column(self, installer):
        assert installer.table.horizontalHeaderItem(0).text() == "Plugin Name"

    def test_table_header_action_column(self, installer):
        assert installer.table.horizontalHeaderItem(5).text() == "Action"

    def test_progress_bar_initially_hidden(self, installer):
        assert installer._progress_bar.isHidden()

    def test_upgrade_button_initially_hidden(self, installer):
        assert installer.btn_upgrade_app.isHidden()

    def test_update_all_initially_disabled(self, installer):
        assert not installer.btn_update_all.isEnabled()

    def test_table_empty_before_fetch(self, installer):
        assert installer.table.rowCount() == 0

    def test_remote_data_initially_empty(self, installer):
        assert installer.remote_data == []
