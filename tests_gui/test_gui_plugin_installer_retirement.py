"""Installed hidden entries remain visible, with explicit migration guidance."""

from unittest.mock import MagicMock, patch

import pytest
from test_gui_plugin_installer import _plugin_installer as installer_module


@pytest.fixture
def bare_installer(qapp):
    host = MagicMock()
    host.plugin_manager.plugins = []
    with patch.object(
        installer_module.PluginInstallerWindow, "check_updates", lambda self: None
    ):
        window = installer_module.PluginInstallerWindow(None, auto_check=False)
    window.main_window = host
    window.get_app_version = lambda: "5.0.0"
    window.latest_app_version = "Checking..."
    yield window
    window.destroy()


def _row(window, name):
    for index in range(window.table.rowCount()):
        if window.table.item(index, 0).text() == name:
            return index
    raise AssertionError(f"Missing plugin row: {name}")


def _populate(
    window,
    *,
    local_version="1.0",
    hidden=True,
    disabled=False,
    migration="new",
    replacement_name="Replacement Pro",
):
    old = {
        "id": "old",
        "name": "Old",
        "visible": not hidden,
        "version": "9.0",
        "author": "A",
        "description": "Old plugin",
        "downloadUrl": "https://example.test/old.py",
        "supported_moleditpy_version": ">=5.0.0, <6.0.0",
    }
    if migration is not None:
        old["migration_to"] = migration
    new = {
        "id": "new",
        "name": replacement_name,
        "visible": True,
        "version": "1.0",
        "author": "A",
    }
    window.main_window.plugin_manager.plugins = [
        {
            "name": "Old",
            "version": local_version,
            "filepath": None,
            "disabled": disabled,
        }
    ]
    window.remote_data = [old, new]
    window.populate_table(silent=True)
    return _row(window, "Old")


@pytest.mark.parametrize("version", ["1.0", "9.0", "10.0", "Unknown"])
def test_hidden_installed_status_overrides_version_comparison(bare_installer, version):
    from PySide6.QtCore import Qt

    row = _populate(bare_installer, local_version=version)
    status = bare_installer.table.item(row, 4)
    assert status.text() == "Stale"
    assert status.data(Qt.ItemDataRole.UserRole) == "Stale"
    assert bare_installer.table.cellWidget(row, 5).text() == "Reinstall"
    assert not bare_installer.btn_update_all.isEnabled()
    assert not bare_installer.updates_found
    assert "1 installed plugin is stale" == bare_installer._status_label.text()


def test_tooltips_resolve_migration_id_to_current_display_name(bare_installer):
    row = _populate(bare_installer, replacement_name="Renamed Pro")
    for column in (0, 4):
        tooltip = bare_installer.table.item(row, column).toolTip()
        assert "Hidden in the plugin registry" in tooltip
        assert "Migration to: Renamed Pro" in tooltip


def test_hidden_uninstalled_entry_stays_out_of_table(bare_installer):
    _populate(bare_installer)
    bare_installer.main_window.plugin_manager.plugins = []
    bare_installer.populate_table(silent=True)
    assert bare_installer.table.rowCount() == 1
    assert bare_installer.table.item(0, 0).text() == "Replacement Pro"


def test_hidden_disabled_keeps_both_statuses_and_explanations(bare_installer):
    row = _populate(bare_installer, disabled=True)
    status = bare_installer.table.item(row, 4)
    assert status.text() == "Stale (Disabled)"
    assert "not loaded" in status.toolTip()
    assert "Migration to: Replacement Pro" in status.toolTip()


@pytest.mark.parametrize("migration", [None, "", [], 12])
def test_missing_or_invalid_migration_shows_only_stale(bare_installer, migration):
    row = _populate(bare_installer, migration=migration)
    status = bare_installer.table.item(row, 4)
    assert status.text() == "Stale"
    assert "Migration to:" not in status.toolTip()


def test_unresolved_migration_is_displayed_as_plain_reference(bare_installer):
    row = _populate(bare_installer, migration="missing_id")
    assert "Migration to: missing_id" in bare_installer.table.item(row, 4).toolTip()


def test_incompatible_stale_entry_keeps_compatibility_warning(bare_installer):
    _populate(bare_installer)
    bare_installer.remote_data[0]["supported_moleditpy_version"] = ">=99.0"
    bare_installer.remote_data[0]["supported_python_version"] = ">=99.0"
    bare_installer.remote_data[0]["supported_os"] = ["unsupported"]
    bare_installer.populate_table(silent=True)
    status = bare_installer.table.item(_row(bare_installer, "Old"), 4)
    assert status.text() == "Stale"
    assert "Requires MoleditPy" in status.toolTip()
    assert "Requires Python" in status.toolTip()
    assert "Requires OS" in status.toolTip()
    assert "Migration to:" in status.toolTip()


def test_missing_registry_entry_is_not_classified_as_stale(bare_installer):
    _populate(bare_installer)
    bare_installer.remote_data = []
    bare_installer.populate_table(silent=True)
    assert bare_installer.table.item(0, 4).text() == "Not in Registry"


def test_visible_entry_keeps_update_behavior_and_ignores_retirement_notice(
    bare_installer,
):
    row = _populate(bare_installer, hidden=False)
    assert bare_installer.table.item(row, 4).text() == "Update Available"
    assert bare_installer.table.cellWidget(row, 5).text() == "Update"
    assert bare_installer.btn_update_all.isEnabled()
    assert "Migration to:" not in bare_installer.table.item(row, 4).toolTip()


def test_refresh_to_hidden_clears_update_recommendation(bare_installer):
    _populate(bare_installer, hidden=False)
    bare_installer.remote_data[0]["visible"] = False
    bare_installer.populate_table(silent=True)
    assert bare_installer.table.item(_row(bare_installer, "Old"), 4).text() == "Stale"
    assert not bare_installer.updates_found
    assert not bare_installer.btn_update_all.isEnabled()


def test_pending_install_is_recognized_when_entry_becomes_hidden(bare_installer):
    _populate(bare_installer)
    bare_installer._pending_installs["Old"] = (
        bare_installer.main_window.plugin_manager.plugins[0]
    )
    bare_installer.main_window.plugin_manager.plugins = []
    bare_installer.populate_table(silent=True)
    assert bare_installer.table.item(_row(bare_installer, "Old"), 4).text() == "Stale"


def test_update_all_excludes_stale_but_includes_active_update(
    bare_installer, monkeypatch
):
    from PySide6.QtWidgets import QMessageBox

    _populate(bare_installer)
    bare_installer.remote_data[1]["downloadUrl"] = "https://example.test/new.py"
    bare_installer.main_window.plugin_manager.plugins.append(
        {"name": "Replacement Pro", "version": "0.1", "filepath": None}
    )
    bare_installer.populate_table(silent=True)
    confirm = MagicMock(return_value=QMessageBox.StandardButton.No)
    monkeypatch.setattr(installer_module.QMessageBox, "question", confirm)
    bare_installer.update_all_plugins()
    text = confirm.call_args.args[2]
    assert "Replacement Pro" in text
    assert "Old" not in text
    assert (
        "1 plugin has an update available; 1 installed plugin is stale"
        == bare_installer._status_label.text()
    )


def test_details_receive_stale_flag_and_resolved_migration(bare_installer):
    row = _populate(bare_installer)
    with patch.object(installer_module, "PluginDetailsDialog") as dialog:
        bare_installer.show_plugin_details(row, 0)
    assert dialog.call_args.kwargs["registry_hidden"] is True
    assert dialog.call_args.kwargs["migration_target"] == "Replacement Pro"
    dialog.return_value.exec.assert_called_once()


@pytest.mark.parametrize("disabled", [False, True])
def test_real_details_show_stale_migration_and_optional_disabled(qapp, disabled):
    from PySide6.QtWidgets import QLabel

    dialog = installer_module.PluginDetailsDialog(
        None,
        "Old",
        "A",
        "1.0",
        "Desc",
        [],
        {"disabled": disabled},
        None,
        registry_hidden=True,
        migration_target="New <Pro> & Tools",
    )
    try:
        texts = [label.text() for label in dialog.findChildren(QLabel)]
        assert any("Stale" in text for text in texts)
        assert any(
            "Migration to:" in text and "New &lt;Pro&gt; &amp; Tools" in text
            for text in texts
        )
        assert any("Disabled" in text for text in texts) is disabled
    finally:
        dialog.destroy()


def test_details_without_target_omit_migration_line(qapp):
    from PySide6.QtWidgets import QLabel

    dialog = installer_module.PluginDetailsDialog(
        None, "Old", "A", "1.0", "Desc", [], {}, None, registry_hidden=True
    )
    try:
        texts = [label.text() for label in dialog.findChildren(QLabel)]
        assert any("Stale" in text for text in texts)
        assert not any("Migration to:" in text for text in texts)
    finally:
        dialog.destroy()
