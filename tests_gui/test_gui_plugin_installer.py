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


# ===========================================================================
# Plugin Installer — supported_python_version compatibility (GUI)
# ===========================================================================

PY_OK_SPEC = ">=3.9, <3.15"
PY_BAD_SPEC = ">=99.0"
APP_OK_SPEC = ">=4.0.0, <5.0.0"
APP_BAD_SPEC = ">=99.0"


def _remote_entry(name, *, moleditpy_spec=APP_OK_SPEC, python_spec=None):
    entry = {
        "name": name,
        "visible": True,
        "version": "1.0.0",
        "author": "Test Author",
        "description": "d",
        "supported_moleditpy_version": moleditpy_spec,
        "downloadUrl": "../plugins/Fake/fake.py",
    }
    if python_spec is not None:
        entry["supported_python_version"] = python_spec
    return entry


class TestPythonVersionCompatGui:
    """populate_table() + install path honour supported_python_version."""

    @pytest.fixture
    def installer(self, qapp):
        from unittest.mock import MagicMock, patch

        mw = MagicMock()
        mw.plugin_manager.plugins = []
        with patch.object(
            _plugin_installer.PluginInstallerWindow,
            "check_updates",
            lambda self: None,
        ):
            d = _plugin_installer.PluginInstallerWindow(None, auto_check=False)
        d.main_window = mw
        d.get_app_version = lambda: "4.2.0"
        d.latest_app_version = "Checking..."
        yield d
        d.destroy()

    @staticmethod
    def _row_of(installer, name):
        for row in range(installer.table.rowCount()):
            if installer.table.item(row, 0).text() == name:
                return row
        raise AssertionError(f"row for {name!r} not found")

    def _populate(self, installer, entries):
        installer.remote_data = entries
        installer.populate_table(silent=True)

    def test_python_incompatible_marks_incompatible(self, installer):
        self._populate(
            installer, [_remote_entry("Py Bad", python_spec=PY_BAD_SPEC)]
        )
        row = self._row_of(installer, "Py Bad")
        status = installer.table.item(row, 4)
        assert status.text() == "Incompatible"
        assert "Python" in status.toolTip()
        assert PY_BAD_SPEC in status.toolTip()

    def test_python_compatible_is_not_installed(self, installer):
        self._populate(
            installer, [_remote_entry("Py Ok", python_spec=PY_OK_SPEC)]
        )
        row = self._row_of(installer, "Py Ok")
        status = installer.table.item(row, 4)
        assert status.text() == "Not Installed"
        assert "Python" not in status.toolTip()

    def test_missing_python_spec_is_compatible(self, installer):
        self._populate(installer, [_remote_entry("No Spec")])
        row = self._row_of(installer, "No Spec")
        assert installer.table.item(row, 4).text() == "Not Installed"

    def test_both_incompatible_tooltip_has_both_lines(self, installer):
        self._populate(
            installer,
            [
                _remote_entry(
                    "Both Bad",
                    moleditpy_spec=APP_BAD_SPEC,
                    python_spec=PY_BAD_SPEC,
                )
            ],
        )
        row = self._row_of(installer, "Both Bad")
        tip = installer.table.item(row, 4).toolTip()
        assert "Requires MoleditPy" in tip
        assert "Requires Python" in tip

    def test_install_button_tooltip_mentions_python(self, installer):
        self._populate(
            installer, [_remote_entry("Py Bad", python_spec=PY_BAD_SPEC)]
        )
        row = self._row_of(installer, "Py Bad")
        btn = installer.table.cellWidget(row, 5)
        assert btn is not None
        assert "Python" in btn.toolTip()
        assert PY_BAD_SPEC in btn.toolTip()

    def test_install_button_no_python_tooltip_when_compatible(self, installer):
        self._populate(
            installer, [_remote_entry("Py Ok", python_spec=PY_OK_SPEC)]
        )
        row = self._row_of(installer, "Py Ok")
        btn = installer.table.cellWidget(row, 5)
        assert btn is not None
        assert "Python" not in btn.toolTip()

    # -- install-time warning path ------------------------------------------

    def _click_install(self, installer, name):
        """Fire on_update_clicked via a real button so self.sender() works.

        download_url is left unset: after the compatibility gates the handler
        aborts at 'Missing download URL', so no network is ever touched.
        """
        from PyQt6.QtWidgets import QPushButton

        btn = QPushButton(parent=installer)
        btn.setProperty("plugin_name", name)
        btn.setProperty("target_file", None)
        btn.setProperty("dependencies", [])
        btn.clicked.connect(installer.on_update_clicked)
        btn.click()

    def test_install_warning_no_aborts(self, installer, qapp):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        installer.remote_data = [
            _remote_entry("Py Bad", python_spec=PY_BAD_SPEC)
        ]
        with patch.object(
            _plugin_installer.QMessageBox,
            "warning",
            return_value=QMessageBox.StandardButton.No,
        ) as warn:
            self._click_install(installer, "Py Bad")
        assert warn.call_count == 1
        assert "Incompatible Python Version Warning" in warn.call_args[0]

    def test_install_warning_yes_proceeds_past_check(self, installer, qapp):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        installer.remote_data = [
            _remote_entry("Py Bad", python_spec=PY_BAD_SPEC)
        ]
        with patch.object(
            _plugin_installer.QMessageBox,
            "warning",
            return_value=QMessageBox.StandardButton.Yes,
        ) as warn:
            self._click_install(installer, "Py Bad")
        # Warning #1: python check (answered Yes) → proceeds; warning #2:
        # the missing-download-URL abort, proving the gate was passed.
        assert warn.call_count == 2
        assert "Incompatible Python Version Warning" in warn.call_args_list[0][0]
        assert any(
            "Missing download URL" in str(a) for a in warn.call_args_list[1][0]
        )

    def test_install_compatible_no_python_warning(self, installer, qapp):
        from unittest.mock import patch
        from PyQt6.QtWidgets import QMessageBox

        installer.remote_data = [
            _remote_entry("Py Ok", python_spec=PY_OK_SPEC)
        ]
        with patch.object(
            _plugin_installer.QMessageBox,
            "warning",
            return_value=QMessageBox.StandardButton.No,
        ) as warn:
            self._click_install(installer, "Py Ok")
        # Only the missing-download-URL abort fires — no compatibility warning.
        assert warn.call_count == 1
        assert any("Missing download URL" in str(a) for a in warn.call_args[0])


class TestDetailsDialogSupportedPython:
    """PluginDetailsDialog shows the Supported Python line."""

    def _find_labels(self, dlg):
        from PyQt6.QtWidgets import QLabel

        return [w.text() for w in dlg.findChildren(QLabel)]

    def test_shows_passed_spec(self, qapp):
        d = _plugin_installer.PluginDetailsDialog(
            None, "X", "A", "1.0", "Desc", [], None, None,
            ">=4.0.0", PY_OK_SPEC,
        )
        import html as _html

        texts = self._find_labels(d)
        escaped = _html.escape(PY_OK_SPEC)
        assert any(
            "Supported Python:" in t and escaped in t for t in texts
        )
        d.destroy()

    def test_defaults_to_unknown(self, qapp):
        d = _make_details_dialog()
        texts = self._find_labels(d)
        assert any(
            "Supported Python:" in t and "Unknown" in t for t in texts
        )
        d.destroy()
