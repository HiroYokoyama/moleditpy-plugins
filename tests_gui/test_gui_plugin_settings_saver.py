"""
Headless GUI tests for the Settings Saver plugin.

Covers: SettingsSaverDialog.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

SETTINGS_PATH = PLUGINS_DIR / "Settings_Saver" / "settings_saver.py"

with mock_chemistry_imports():
    _settings = load_plugin_for_gui(SETTINGS_PATH)

# Alias used inside nested stub classes, where `_settings` would shadow a
# local parameter name.
_settings_mod = _settings

from PyQt6.QtCore import Qt  # noqa: E402
from PyQt6.QtGui import QColor  # noqa: E402


# ===========================================================================
# SettingsSaverDialog  (Settings Saver)
# ===========================================================================


def _settings_context() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


class TestSettingsSaverDialog:
    """SettingsSaverDialog with MagicMock context and parent=None."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _settings_context()
        d = _settings.SettingsSaverDialog(context=ctx, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Settings Saver Manager"

    def test_preset_list_initially_empty(self, dlg):
        assert dlg.preset_list.count() == 0

    def test_load_button_exists(self, dlg):
        assert dlg.btn_load.text() == "Load Preset"

    def test_save_button_exists(self, dlg):
        assert dlg.btn_save.text() == "Save New..."

    def test_delete_button_exists(self, dlg):
        assert dlg.btn_delete.text() == "Delete"

    def test_embed_checkbox_unchecked_by_default(self, dlg):
        assert not dlg.chk_embed.isChecked()

    def test_global_default_button_disabled_when_embed_off(self, dlg):
        assert not dlg.btn_set_global.isEnabled()

    def test_export_button_has_menu(self, dlg):
        assert dlg.btn_export.menu() is not None

    def test_close_button_exists(self, dlg):
        assert dlg.btn_close.text() == "Close"


# ===========================================================================
# SettingsSaverDialog: real event-handler methods driven through a real
# QDialog, with a fake main window and patched QMessageBox/QInputDialog/
# QFileDialog static methods so nothing modal actually appears.
# ===========================================================================


class _FakeInitManager2:
    def __init__(self, settings=None):
        self.settings = settings if settings is not None else {}
        self.settings_dirty = None
        self.initial_settings = {}
        self.saved = False

    def save_settings(self):
        self.saved = True


class _FakeEditActionsManager:
    def __init__(self):
        self.settings_dirty = None


from PyQt6.QtWidgets import QWidget  # noqa: E402


class _FakeMW2(QWidget):
    """A real QWidget (so it can be used as a QDialog parent) with fake manager attrs."""

    def __init__(self, settings=None):
        super().__init__()
        self.init_manager = _FakeInitManager2(settings)
        self.edit_actions_manager = _FakeEditActionsManager()
        # The real main window owns initial_settings; init_manager owns
        # settings_dirty (see moleditpy/ui/main_window.py).
        self.initial_settings = None


class _BareMW(QWidget):
    """QWidget whose init_manager lacks 'settings'/'save_settings' attributes."""

    def __init__(self):
        super().__init__()
        self.init_manager = SimpleNamespace()


def _make_dialog(monkeypatch, tmp_path, mw=None, library_data=None, embed_enabled=False):
    data_path = str(tmp_path / "settings_saver.json")
    monkeypatch.setattr(_settings, "get_plugin_data_path", lambda: data_path)
    if library_data is not None:
        Path(data_path).write_text(json.dumps(library_data), encoding="utf-8")
    _settings.EMBED_SETTINGS["enabled"] = embed_enabled
    ctx = _settings_context()
    d = _settings.SettingsSaverDialog(context=ctx, parent=mw)
    return d, data_path


class TestSettingsSaverDialogRealMethods:
    """Exercise the dialog's actual event-handler bodies (not extracted)."""

    @pytest.fixture(autouse=True)
    def _reset_globals(self):
        yield
        _settings.EMBED_SETTINGS["enabled"] = False
        _settings.PROJECT_PRESETS.clear()
        _settings.PLUGIN_CONFIG["always_save_to_project"] = False

    def test_init_with_project_mode_resets_dirty_and_syncs_initial_settings(
        self, qapp, tmp_path, monkeypatch
    ):
        mw = _FakeMW2(settings={"a": 1})
        mw.init_manager.settings_dirty = True
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw, embed_enabled=True)
        try:
            assert mw.init_manager.settings_dirty is False
            assert mw.initial_settings == {"a": 1}
        finally:
            d.destroy()

    def test_refresh_list_populates_library_and_project_items(
        self, qapp, tmp_path, monkeypatch
    ):
        library = {"Alpha": {"x": 1}, "Zeta": {"y": 2}}
        _settings.PROJECT_PRESETS["Project Settings"] = {"z": 3}
        d, _ = _make_dialog(monkeypatch, tmp_path, library_data=library)
        try:
            assert d.preset_list.count() == 3
        finally:
            d.destroy()

    def test_on_always_save_toggled_persists_config(self, qapp, tmp_path, monkeypatch):
        d, data_path = _make_dialog(monkeypatch, tmp_path)
        try:
            d.chk_always_save.setChecked(True)
            assert _settings.PLUGIN_CONFIG["always_save_to_project"] is True
            saved = json.loads(Path(data_path).read_text(encoding="utf-8"))
            assert saved["_PLUGIN_CONFIG"]["always_save_to_project"] is True
        finally:
            d.destroy()

    def test_on_embed_toggled_true_enables_project_mode(self, qapp, tmp_path, monkeypatch):
        mw = _FakeMW2()
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw)
        try:
            d.chk_embed.setChecked(True)
            assert _settings.EMBED_SETTINGS["enabled"] is True
            assert d.btn_set_global.isEnabled()
        finally:
            d.destroy()

    def test_on_embed_toggled_false_disables_and_marks_dirty(self, qapp, tmp_path, monkeypatch):
        mw = _FakeMW2()
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw, embed_enabled=True)
        try:
            d.chk_embed.setChecked(False)
            assert _settings.EMBED_SETTINGS["enabled"] is False
            assert mw.init_manager.settings_dirty is True
        finally:
            d.destroy()

    def test_on_load_no_selection_warns(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            with patch.object(_settings.QMessageBox, "warning") as warn:
                d.on_load()
            warn.assert_called_once()
        finally:
            d.destroy()

    def test_on_load_applies_library_preset(self, qapp, tmp_path, monkeypatch):
        mw = _FakeMW2(settings={})
        library = {"P1": {"x": 42}}
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw, library_data=library)
        try:
            d.preset_list.setCurrentRow(0)
            d.on_load()
            assert mw.init_manager.settings.get("x") == 42
        finally:
            d.destroy()

    def test_on_load_applies_project_preset(self, qapp, tmp_path, monkeypatch):
        mw = _FakeMW2(settings={})
        _settings.PROJECT_PRESETS["Project Settings"] = {"z": 9}
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw)
        try:
            d.preset_list.setCurrentRow(0)
            d.on_load()
            assert mw.init_manager.settings.get("z") == 9
        finally:
            d.destroy()

    def test_on_save_adds_new_preset(self, qapp, tmp_path, monkeypatch):
        mw = _FakeMW2(settings={"k": "v"})
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw)
        try:
            with patch.object(
                _settings.QInputDialog, "getText", return_value=("NewPreset", True)
            ), patch.object(_settings.QMessageBox, "information"):
                d.on_save()
            assert d.library.get("NewPreset") == {"k": "v"}
            assert d.preset_list.count() == 1
        finally:
            d.destroy()

    def test_on_save_rejects_underscore_name(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            with patch.object(
                _settings.QInputDialog, "getText", return_value=("_bad", True)
            ), patch.object(_settings.QMessageBox, "warning") as warn:
                d.on_save()
            warn.assert_called_once()
            assert "_bad" not in d.library
        finally:
            d.destroy()

    def test_on_save_cancelled_noop(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            with patch.object(_settings.QInputDialog, "getText", return_value=("", False)):
                d.on_save()
            assert d.library == {}
        finally:
            d.destroy()

    def test_on_save_overwrite_confirmed(self, qapp, tmp_path, monkeypatch):
        mw = _FakeMW2(settings={"k": 2})
        library = {"Existing": {"k": 1}}
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw, library_data=library)
        try:
            with patch.object(
                _settings.QInputDialog, "getText", return_value=("Existing", True)
            ), patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.Yes,
            ), patch.object(_settings.QMessageBox, "information"):
                d.on_save()
            assert d.library["Existing"] == {"k": 2}
        finally:
            d.destroy()

    def test_on_save_overwrite_declined(self, qapp, tmp_path, monkeypatch):
        mw = _FakeMW2(settings={"k": 2})
        library = {"Existing": {"k": 1}}
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw, library_data=library)
        try:
            with patch.object(
                _settings.QInputDialog, "getText", return_value=("Existing", True)
            ), patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.No,
            ):
                d.on_save()
            assert d.library["Existing"] == {"k": 1}
        finally:
            d.destroy()

    def test_on_save_missing_settings_warns(self, qapp, tmp_path, monkeypatch):
        mw = _BareMW()  # init_manager has no "settings" attr
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw)
        try:
            with patch.object(
                _settings.QInputDialog, "getText", return_value=("Name", True)
            ), patch.object(_settings.QMessageBox, "warning") as warn:
                d.on_save()
            warn.assert_called_once()
        finally:
            d.destroy()

    def test_on_delete_removes_library_preset(self, qapp, tmp_path, monkeypatch):
        library = {"ToDelete": {"a": 1}}
        d, _ = _make_dialog(monkeypatch, tmp_path, library_data=library)
        try:
            d.preset_list.setCurrentRow(0)
            with patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.Yes,
            ):
                d.on_delete()
            assert "ToDelete" not in d.library
        finally:
            d.destroy()

    def test_on_delete_project_preset_blocked(self, qapp, tmp_path, monkeypatch):
        _settings.PROJECT_PRESETS["Project Settings"] = {"a": 1}
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            d.preset_list.setCurrentRow(0)
            with patch.object(_settings.QMessageBox, "information") as info:
                d.on_delete()
            info.assert_called_once()
        finally:
            d.destroy()

    def test_on_export_single_writes_selected_preset(self, qapp, tmp_path, monkeypatch):
        library = {"P1": {"x": 1}}
        d, _ = _make_dialog(monkeypatch, tmp_path, library_data=library)
        out = tmp_path / "out.json"
        try:
            d.preset_list.setCurrentRow(0)
            with patch.object(
                _settings.QFileDialog, "getSaveFileName", return_value=(str(out), "")
            ), patch.object(_settings.QMessageBox, "information"):
                d.on_export_single()
            assert json.loads(out.read_text(encoding="utf-8")) == {"x": 1}
        finally:
            d.destroy()

    def test_on_export_all_writes_all_presets(self, qapp, tmp_path, monkeypatch):
        library = {"P1": {"x": 1}, "P2": {"y": 2}}
        d, _ = _make_dialog(monkeypatch, tmp_path, library_data=library)
        out = tmp_path / "all.json"
        try:
            with patch.object(
                _settings.QFileDialog, "getSaveFileName", return_value=(str(out), "")
            ), patch.object(_settings.QMessageBox, "information"):
                d.on_export_all()
            data = json.loads(out.read_text(encoding="utf-8"))
            assert set(data.keys()) == {"P1", "P2"}
        finally:
            d.destroy()

    def test_on_export_all_no_presets_shows_info(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            with patch.object(_settings.QMessageBox, "information") as info:
                d.on_export_all()
            info.assert_called_once()
        finally:
            d.destroy()

    def test_on_import_library_format(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        src = tmp_path / "in.json"
        src.write_text(json.dumps({"P1": {"x": 1}, "_PLUGIN_CONFIG": {"a": 1}}), encoding="utf-8")
        try:
            with patch.object(
                _settings.QFileDialog, "getOpenFileName", return_value=(str(src), "")
            ), patch.object(_settings.QMessageBox, "information"):
                d.on_import()
            assert d.library.get("P1") == {"x": 1}
            assert "_PLUGIN_CONFIG" not in d.library
        finally:
            d.destroy()

    def test_on_import_no_path_noop(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            with patch.object(
                _settings.QFileDialog, "getOpenFileName", return_value=("", "")
            ):
                d.on_import()  # no path -> early return
            assert d.library == {}
        finally:
            d.destroy()

    def test_on_import_invalid_json_shows_critical(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        src = tmp_path / "bad.json"
        src.write_text("not valid json", encoding="utf-8")
        try:
            with patch.object(
                _settings.QFileDialog, "getOpenFileName", return_value=(str(src), "")
            ), patch.object(_settings.QMessageBox, "critical") as crit:
                d.on_import()
            crit.assert_called_once()
        finally:
            d.destroy()

    def test_on_import_non_dict_shows_critical(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        src = tmp_path / "list.json"
        src.write_text(json.dumps([1, 2, 3]), encoding="utf-8")
        try:
            with patch.object(
                _settings.QFileDialog, "getOpenFileName", return_value=(str(src), "")
            ), patch.object(_settings.QMessageBox, "critical") as crit:
                d.on_import()
            crit.assert_called_once()
        finally:
            d.destroy()

    def test_on_import_single_preset_format(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        src = tmp_path / "single.json"
        src.write_text(json.dumps({"opacity": 0.5}), encoding="utf-8")
        try:
            with patch.object(
                _settings.QFileDialog, "getOpenFileName", return_value=(str(src), "")
            ), patch.object(
                _settings.QInputDialog, "getText", return_value=("Imported1", True)
            ), patch.object(_settings.QMessageBox, "information"):
                d.on_import()
            assert d.library.get("Imported1") == {"opacity": 0.5}
        finally:
            d.destroy()

    def test_on_import_single_underscore_name_rejected(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        src = tmp_path / "single.json"
        src.write_text(json.dumps({"opacity": 0.5}), encoding="utf-8")
        try:
            with patch.object(
                _settings.QFileDialog, "getOpenFileName", return_value=(str(src), "")
            ), patch.object(
                _settings.QInputDialog, "getText", return_value=("_bad", True)
            ), patch.object(_settings.QMessageBox, "warning") as warn:
                d.on_import()
            warn.assert_called_once()
            assert "_bad" not in d.library
        finally:
            d.destroy()

    def test_on_import_single_cancelled_noop(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        src = tmp_path / "single.json"
        src.write_text(json.dumps({"opacity": 0.5}), encoding="utf-8")
        try:
            with patch.object(
                _settings.QFileDialog, "getOpenFileName", return_value=(str(src), "")
            ), patch.object(_settings.QInputDialog, "getText", return_value=("", False)):
                d.on_import()
            assert d.library == {}
        finally:
            d.destroy()

    def test_on_save_as_global_default_cancelled_noop(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            with patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.No,
            ):
                d.on_save_as_global_default()  # must not raise
        finally:
            d.destroy()

    def test_on_save_as_global_default_no_main_window_warns(
        self, qapp, tmp_path, monkeypatch
    ):
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=None)
        try:
            with patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.Yes,
            ), patch.object(_settings.QMessageBox, "warning") as warn:
                d.on_save_as_global_default()
            warn.assert_called_once()
        finally:
            d.destroy()

    def test_on_save_as_global_default_executes_save_settings(
        self, qapp, tmp_path, monkeypatch
    ):
        mw = _FakeMW2(settings={"a": 1})
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw)
        try:
            with patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.Yes,
            ), patch.object(_settings.QMessageBox, "information") as info:
                d.on_save_as_global_default()
            assert mw.init_manager.saved is True
            info.assert_called_once()
        finally:
            d.destroy()

    def test_on_save_as_global_default_uses_original_save_settings(
        self, qapp, tmp_path, monkeypatch
    ):
        mw = _FakeMW2(settings={"a": 1})
        orig = MagicMock()
        mw._original_save_settings = orig
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw)
        try:
            with patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.Yes,
            ), patch.object(_settings.QMessageBox, "information") as info:
                d.on_save_as_global_default()
            orig.assert_called_once()
            info.assert_called_once()
        finally:
            d.destroy()

    def test_on_save_as_global_default_no_save_func_warns(self, qapp, tmp_path, monkeypatch):
        mw = _BareMW()  # init_manager has no save_settings attr
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw)
        try:
            with patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.Yes,
            ), patch.object(_settings.QMessageBox, "warning") as warn:
                d.on_save_as_global_default()
            warn.assert_called_once()
        finally:
            d.destroy()

    def test_on_save_as_global_default_exception_shows_critical(
        self, qapp, tmp_path, monkeypatch
    ):
        mw = _FakeMW2(settings={"a": 1})
        mw.init_manager.save_settings = MagicMock(side_effect=RuntimeError("boom"))
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw)
        try:
            with patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.Yes,
            ), patch.object(_settings.QMessageBox, "critical") as crit:
                d.on_save_as_global_default()
            crit.assert_called_once()
        finally:
            d.destroy()


# ===========================================================================
# PresetEditorDialog
# ===========================================================================


def _editor(settings, skipped=None, name="P1", existing=()):
    return _settings.PresetEditorDialog(
        name, settings, skipped=skipped, existing_names=existing, parent=None
    )


def _row_of(dlg, key):
    for row in range(dlg.table.rowCount()):
        if dlg.table.item(row, dlg.COL_KEY).text() == key:
            return row
    raise AssertionError(f"row for {key!r} not found")


def _set_value(dlg, key, text):
    dlg.table.item(_row_of(dlg, key), dlg.COL_VALUE).setText(text)


def _set_checked(dlg, key, checked):
    state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
    dlg.table.item(_row_of(dlg, key), dlg.COL_INCLUDE).setCheckState(state)


class TestPresetEditorDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _editor({"bg": "#FFF", "size": 12, "flag": True})
        yield d
        d.destroy()

    def test_title_names_the_preset(self, dlg):
        assert "P1" in dlg.windowTitle()

    def test_name_field_prefilled(self, dlg):
        assert dlg.edit_name.text() == "P1"

    def test_one_row_per_key(self, dlg):
        assert dlg.table.rowCount() == 3

    def test_rows_sorted_by_key(self, dlg):
        keys = [
            dlg.table.item(r, dlg.COL_KEY).text() for r in range(dlg.table.rowCount())
        ]
        assert keys == sorted(keys)

    def test_all_keys_checked_by_default(self, qapp, dlg):
        for row in range(dlg.table.rowCount()):
            assert (
                dlg.table.item(row, dlg.COL_INCLUDE).checkState()
                == Qt.CheckState.Checked
            )

    def test_skipped_keys_start_unchecked(self, qapp):
        d = _editor({"bg": "#FFF"}, skipped={"size": 12})
        try:
            assert d.table.rowCount() == 2
            assert (
                d.table.item(_row_of(d, "size"), d.COL_INCLUDE).checkState()
                == Qt.CheckState.Unchecked
            )
        finally:
            d.destroy()

    def test_settings_win_over_a_stale_skipped_entry(self, qapp):
        d = _editor({"bg": "#FFF"}, skipped={"bg": "#000"})
        try:
            assert d.table.rowCount() == 1
            assert d.table.item(0, d.COL_VALUE).text() == "#FFF"
        finally:
            d.destroy()

    def test_key_column_is_not_editable(self, dlg):
        item = dlg.table.item(0, dlg.COL_KEY)
        assert not (item.flags() & Qt.ItemFlag.ItemIsEditable)

    # -- collect() --------------------------------------------------------

    def test_collect_returns_unchanged_values(self, dlg):
        settings, skipped, errors = dlg.collect()
        assert settings == {"bg": "#FFF", "size": 12, "flag": True}
        assert skipped == {}
        assert errors == []

    def test_collect_picks_up_an_edited_value(self, dlg):
        _set_value(dlg, "size", "20")
        settings, _, errors = dlg.collect()
        assert errors == []
        assert settings["size"] == 20

    def test_collect_keeps_string_type_for_string_keys(self, dlg):
        _set_value(dlg, "bg", "123")
        settings, _, _ = dlg.collect()
        assert settings["bg"] == "123"

    def test_collect_moves_unchecked_key_to_skipped(self, dlg):
        _set_checked(dlg, "size", False)
        settings, skipped, _ = dlg.collect()
        assert "size" not in settings
        assert skipped["size"] == 12

    def test_collect_skips_empty_value_even_when_checked(self, dlg):
        _set_value(dlg, "bg", "")
        settings, skipped, errors = dlg.collect()
        assert "bg" not in settings
        assert "bg" in skipped
        assert errors == []

    def test_collect_skips_whitespace_only_value(self, dlg):
        _set_value(dlg, "bg", "   ")
        settings, skipped, _ = dlg.collect()
        assert "bg" not in settings and "bg" in skipped

    def test_empty_value_stash_keeps_the_previous_value(self, dlg):
        _set_value(dlg, "size", "")
        _, skipped, _ = dlg.collect()
        assert skipped["size"] == 12

    def test_collect_reports_unreadable_values(self, dlg):
        _set_value(dlg, "size", "huge")
        settings, _, errors = dlg.collect()
        assert "size" not in settings
        assert any("size" in e for e in errors)

    # -- row actions ------------------------------------------------------

    def test_add_key_appends_a_checked_row(self, dlg):
        with patch.object(
            _settings.QInputDialog, "getText", return_value=("theme", True)
        ):
            dlg.on_add_key()
        assert "theme" in dlg.existing_keys()

    def test_add_key_rejects_duplicates(self, dlg):
        before = dlg.table.rowCount()
        with patch.object(
            _settings.QInputDialog, "getText", return_value=("bg", True)
        ), patch.object(_settings.QMessageBox, "warning") as warn:
            dlg.on_add_key()
        warn.assert_called_once()
        assert dlg.table.rowCount() == before

    def test_add_key_cancelled_changes_nothing(self, dlg):
        before = dlg.table.rowCount()
        with patch.object(_settings.QInputDialog, "getText", return_value=("", False)):
            dlg.on_add_key()
        assert dlg.table.rowCount() == before

    def test_add_key_asks_for_a_value(self, dlg):
        with patch.object(
            _settings.QInputDialog,
            "getText",
            side_effect=[("theme", True), ("dark", True)],
        ):
            dlg.on_add_key()
        settings, _, errors = dlg.collect()
        assert errors == []
        assert settings["theme"] == "dark"

    def test_add_key_cancelled_at_the_value_prompt_adds_nothing(self, dlg):
        before = dlg.table.rowCount()
        with patch.object(
            _settings.QInputDialog,
            "getText",
            side_effect=[("theme", True), ("", False)],
        ):
            dlg.on_add_key()
        assert dlg.table.rowCount() == before

    def test_added_key_with_empty_value_is_skipped(self, dlg):
        with patch.object(
            _settings.QInputDialog,
            "getText",
            side_effect=[("theme", True), ("", True)],
        ):
            dlg.on_add_key()
        settings, skipped, errors = dlg.collect()
        assert errors == []
        assert "theme" not in settings
        assert "theme" in skipped

    def test_remove_key_drops_the_row(self, dlg):
        dlg.table.setCurrentCell(_row_of(dlg, "size"), dlg.COL_VALUE)
        with patch.object(
            _settings.QMessageBox,
            "question",
            return_value=_settings.QMessageBox.StandardButton.Yes,
        ):
            dlg.on_remove_key()
        assert "size" not in dlg.existing_keys()

    def test_remove_key_declined_keeps_the_row(self, dlg):
        dlg.table.setCurrentCell(_row_of(dlg, "size"), dlg.COL_VALUE)
        with patch.object(
            _settings.QMessageBox,
            "question",
            return_value=_settings.QMessageBox.StandardButton.No,
        ):
            dlg.on_remove_key()
        assert "size" in dlg.existing_keys()

    def test_remove_key_without_selection_informs(self, dlg):
        dlg.table.setCurrentCell(-1, -1)
        with patch.object(_settings.QMessageBox, "information") as info:
            dlg.on_remove_key()
        info.assert_called_once()

    def test_skip_all_then_collect_yields_no_settings(self, dlg):
        dlg.set_all_checked(False)
        settings, skipped, _ = dlg.collect()
        assert settings == {}
        assert set(skipped) == {"bg", "size", "flag"}

    def test_apply_all_rechecks_every_row(self, dlg):
        dlg.set_all_checked(False)
        dlg.set_all_checked(True)
        settings, skipped, _ = dlg.collect()
        assert skipped == {}
        assert len(settings) == 3

    def test_filter_hides_non_matching_rows(self, dlg):
        dlg.edit_filter.setText("si")
        hidden = [
            dlg.table.isRowHidden(r) for r in range(dlg.table.rowCount())
        ]
        assert hidden.count(False) == 1
        assert not dlg.table.isRowHidden(_row_of(dlg, "size"))

    def test_clearing_the_filter_shows_everything(self, dlg):
        dlg.edit_filter.setText("si")
        dlg.edit_filter.setText("")
        assert not any(
            dlg.table.isRowHidden(r) for r in range(dlg.table.rowCount())
        )

    def test_check_all_ignores_filtered_out_rows(self, dlg):
        dlg.edit_filter.setText("si")
        dlg.set_all_checked(False)
        assert (
            dlg.table.item(_row_of(dlg, "bg"), dlg.COL_INCLUDE).checkState()
            == Qt.CheckState.Checked
        )

    # -- on_accept --------------------------------------------------------

    def test_accept_publishes_results(self, dlg):
        _set_checked(dlg, "flag", False)
        dlg.on_accept()
        assert dlg.result_name == "P1"
        assert dlg.result_settings == {"bg": "#FFF", "size": 12}
        assert dlg.result_skipped == {"flag": True}

    def test_accept_carries_a_rename(self, dlg):
        dlg.edit_name.setText("Renamed")
        dlg.on_accept()
        assert dlg.result_name == "Renamed"

    def test_accept_rejects_an_empty_name(self, dlg):
        dlg.edit_name.setText("   ")
        with patch.object(_settings.QMessageBox, "warning") as warn:
            dlg.on_accept()
        warn.assert_called_once()
        assert dlg.result() != _settings.QDialog.DialogCode.Accepted

    def test_accept_rejects_underscore_name(self, dlg):
        dlg.edit_name.setText("_hidden")
        with patch.object(_settings.QMessageBox, "warning") as warn:
            dlg.on_accept()
        warn.assert_called_once()

    def test_accept_blocks_on_unreadable_value(self, dlg):
        _set_value(dlg, "size", "huge")
        with patch.object(_settings.QMessageBox, "warning") as warn:
            dlg.on_accept()
        warn.assert_called_once()
        assert dlg.result_settings == {}

    def test_rename_onto_existing_name_asks_first(self, qapp):
        d = _editor({"bg": "#FFF"}, existing=("P1", "Other"))
        try:
            d.edit_name.setText("Other")
            with patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.No,
            ) as q:
                d.on_accept()
            q.assert_called_once()
            assert d.result_name == "P1"
        finally:
            d.destroy()

    def test_keeping_own_name_does_not_ask(self, qapp):
        d = _editor({"bg": "#FFF"}, existing=("P1", "Other"))
        try:
            with patch.object(_settings.QMessageBox, "question") as q:
                d.on_accept()
            q.assert_not_called()
            assert d.result_name == "P1"
        finally:
            d.destroy()


# ===========================================================================
# SettingsSaverDialog.on_edit
# ===========================================================================


class TestSettingsSaverDialogEdit:
    @pytest.fixture(autouse=True)
    def _reset_globals(self):
        yield
        _settings.EMBED_SETTINGS["enabled"] = False
        _settings.PROJECT_PRESETS.clear()

    def _accept_editor(self, name=None, settings=None, skipped=None):
        """Patch PresetEditorDialog with a stub that accepts immediately."""
        want_name, want_settings, want_skipped = name, settings, skipped

        class _StubEditor:
            def __init__(self, in_name, in_settings, skipped=None,
                         existing_names=(), parent=None):
                self.result_name = want_name if want_name is not None else in_name
                self.result_settings = (
                    dict(want_settings)
                    if want_settings is not None
                    else dict(in_settings)
                )
                self.result_skipped = dict(want_skipped or {})

            def exec(self):
                return _settings_mod.QDialog.DialogCode.Accepted

        return _StubEditor

    def test_edit_button_exists(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            assert d.btn_edit.text() == "Edit..."
        finally:
            d.destroy()

    def test_on_edit_without_selection_warns(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            with patch.object(_settings.QMessageBox, "warning") as warn:
                d.on_edit()
            warn.assert_called_once()
        finally:
            d.destroy()

    def test_on_edit_writes_changes_to_the_library(self, qapp, tmp_path, monkeypatch):
        d, data_path = _make_dialog(
            monkeypatch, tmp_path, library_data={"P1": {"x": 1, "y": 2}}
        )
        try:
            d.preset_list.setCurrentRow(0)
            monkeypatch.setattr(
                _settings,
                "PresetEditorDialog",
                self._accept_editor(settings={"x": 9}, skipped={"y": 2}),
            )
            d.on_edit()
            assert d.library["P1"] == {"x": 9}
            saved = json.loads(Path(data_path).read_text(encoding="utf-8"))
            assert saved["P1"] == {"x": 9}
            assert saved[_settings.SKIPPED_STORE_KEY]["P1"] == {"y": 2}
        finally:
            d.destroy()

    def test_on_edit_renames_and_drops_the_old_entry(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path, library_data={"P1": {"x": 1}})
        try:
            d.preset_list.setCurrentRow(0)
            monkeypatch.setattr(
                _settings, "PresetEditorDialog", self._accept_editor(name="P2")
            )
            d.on_edit()
            assert "P1" not in d.library
            assert d.library["P2"] == {"x": 1}
        finally:
            d.destroy()

    def test_rename_moves_the_skipped_stash(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(
            monkeypatch,
            tmp_path,
            library_data={
                "P1": {"x": 1},
                _settings.SKIPPED_STORE_KEY: {"P1": {"y": 2}},
            },
        )
        try:
            d.preset_list.setCurrentRow(0)
            monkeypatch.setattr(
                _settings,
                "PresetEditorDialog",
                self._accept_editor(name="P2", skipped={"y": 2}),
            )
            d.on_edit()
            assert _settings.get_skipped_map(d.library, "P1") == {}
            assert _settings.get_skipped_map(d.library, "P2") == {"y": 2}
        finally:
            d.destroy()

    def test_on_edit_passes_the_stored_skipped_map_to_the_editor(
        self, qapp, tmp_path, monkeypatch
    ):
        d, _ = _make_dialog(
            monkeypatch,
            tmp_path,
            library_data={
                "P1": {"x": 1},
                _settings.SKIPPED_STORE_KEY: {"P1": {"y": 2}},
            },
        )
        try:
            d.preset_list.setCurrentRow(0)
            seen = {}

            class _Capture:
                def __init__(self, name, settings, skipped=None, existing_names=(),
                             parent=None):
                    seen["skipped"] = dict(skipped or {})
                    seen["settings"] = dict(settings)
                    self.result_name = name
                    self.result_settings = dict(settings)
                    self.result_skipped = dict(skipped or {})

                def exec(self):
                    return _settings_mod.QDialog.DialogCode.Rejected

            monkeypatch.setattr(_settings, "PresetEditorDialog", _Capture)
            d.on_edit()
            assert seen["skipped"] == {"y": 2}
            assert seen["settings"] == {"x": 1}
        finally:
            d.destroy()

    def test_cancelled_editor_leaves_the_library_alone(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path, library_data={"P1": {"x": 1}})
        try:
            d.preset_list.setCurrentRow(0)

            class _Rejecting:
                def __init__(self, *a, **k):
                    self.result_name = "ignored"
                    self.result_settings = {}
                    self.result_skipped = {}

                def exec(self):
                    return _settings_mod.QDialog.DialogCode.Rejected

            monkeypatch.setattr(_settings, "PresetEditorDialog", _Rejecting)
            d.on_edit()
            assert d.library["P1"] == {"x": 1}
        finally:
            d.destroy()

    def test_editing_a_project_preset_updates_project_presets_only(
        self, qapp, tmp_path, monkeypatch
    ):
        _settings.PROJECT_PRESETS["Project Settings"] = {"z": 3}
        d, _ = _make_dialog(monkeypatch, tmp_path)
        try:
            d.preset_list.setCurrentRow(0)
            monkeypatch.setattr(
                _settings, "PresetEditorDialog", self._accept_editor(settings={"z": 4})
            )
            d.on_edit()
            assert _settings.PROJECT_PRESETS["Project Settings"] == {"z": 4}
            assert d.library == {}
        finally:
            d.destroy()

    def test_edit_of_a_preset_without_settings_warns(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(monkeypatch, tmp_path, library_data={"P1": "not-a-dict"})
        try:
            d.preset_list.setCurrentRow(0)
            with patch.object(_settings.QMessageBox, "warning") as warn:
                d.on_edit()
            warn.assert_called_once()
        finally:
            d.destroy()

    def test_deleting_a_preset_clears_its_skipped_stash(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(
            monkeypatch,
            tmp_path,
            library_data={
                "P1": {"x": 1},
                _settings.SKIPPED_STORE_KEY: {"P1": {"y": 2}},
            },
        )
        try:
            d.preset_list.setCurrentRow(0)
            with patch.object(
                _settings.QMessageBox,
                "question",
                return_value=_settings.QMessageBox.StandardButton.Yes,
            ):
                d.on_delete()
            assert "P1" not in d.library
            assert _settings.get_skipped_map(d.library, "P1") == {}
        finally:
            d.destroy()

    def test_skipped_store_is_not_listed_as_a_preset(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(
            monkeypatch,
            tmp_path,
            library_data={
                "P1": {"x": 1},
                _settings.SKIPPED_STORE_KEY: {"P1": {"y": 2}},
            },
        )
        try:
            assert d.preset_list.count() == 1
            assert d.preset_list.item(0).text() == "P1"
        finally:
            d.destroy()


class TestMissingValues:
    """Keys that arrive without a value: skipped, greyed, locked."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _editor({"bg": "#FFF", "gone": None, "blank": "   "})
        yield d
        d.destroy()

    def test_null_value_row_is_flagged_missing(self, dlg):
        assert dlg.is_missing_row(_row_of(dlg, "gone"))

    def test_blank_string_row_is_flagged_missing(self, dlg):
        assert dlg.is_missing_row(_row_of(dlg, "blank"))

    def test_normal_row_is_not_flagged(self, dlg):
        assert not dlg.is_missing_row(_row_of(dlg, "bg"))

    def test_missing_row_value_is_not_editable(self, dlg):
        item = dlg.table.item(_row_of(dlg, "gone"), dlg.COL_VALUE)
        assert not (item.flags() & Qt.ItemFlag.ItemIsEditable)

    def test_missing_row_checkbox_is_locked(self, dlg):
        item = dlg.table.item(_row_of(dlg, "gone"), dlg.COL_INCLUDE)
        assert not (item.flags() & Qt.ItemFlag.ItemIsUserCheckable)
        assert item.checkState() == Qt.CheckState.Unchecked

    def test_missing_row_shows_a_placeholder(self, dlg):
        item = dlg.table.item(_row_of(dlg, "gone"), dlg.COL_VALUE)
        assert item.text() == _settings.MISSING_VALUE_PLACEHOLDER

    def test_missing_row_is_greyed_out(self, dlg):
        row = _row_of(dlg, "gone")
        grey = QColor("gray").name()
        assert dlg.table.item(row, dlg.COL_VALUE).foreground().color().name() == grey
        assert dlg.table.item(row, dlg.COL_KEY).foreground().color().name() == grey

    def test_missing_row_has_a_tooltip(self, dlg):
        row = _row_of(dlg, "gone")
        tip = _settings.MISSING_VALUE_TIP
        assert dlg.table.item(row, dlg.COL_VALUE).toolTip() == tip
        assert dlg.table.item(row, dlg.COL_KEY).toolTip() == tip

    def test_collect_skips_missing_rows_with_their_original_value(self, dlg):
        settings, skipped, errors = dlg.collect()
        assert errors == []
        assert settings == {"bg": "#FFF"}
        assert skipped == {"gone": None, "blank": "   "}

    def test_placeholder_is_never_stored_as_a_value(self, dlg):
        _, skipped, _ = dlg.collect()
        assert _settings.MISSING_VALUE_PLACEHOLDER not in skipped.values()

    def test_apply_all_cannot_revive_a_missing_row(self, dlg):
        dlg.set_all_checked(True)
        settings, _, _ = dlg.collect()
        assert "gone" not in settings

    def test_missing_row_can_still_be_removed(self, dlg):
        dlg.table.setCurrentCell(_row_of(dlg, "gone"), dlg.COL_KEY)
        with patch.object(
            _settings.QMessageBox,
            "question",
            return_value=_settings.QMessageBox.StandardButton.Yes,
        ):
            dlg.on_remove_key()
        assert "gone" not in dlg.existing_keys()

    def test_accept_is_not_blocked_by_missing_rows(self, dlg):
        dlg.on_accept()
        assert dlg.result_settings == {"bg": "#FFF"}


class TestEmptiedValueFeedback:
    """Greying/tooltip applied while the user edits a value down to empty."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _editor({"bg": "#FFF"})
        yield d
        d.destroy()

    def test_emptying_a_value_greys_it(self, dlg):
        _set_value(dlg, "bg", "")
        item = dlg.table.item(_row_of(dlg, "bg"), dlg.COL_VALUE)
        assert item.foreground().color().name() == QColor("gray").name()

    def test_emptying_a_value_adds_a_tooltip(self, dlg):
        _set_value(dlg, "bg", "")
        item = dlg.table.item(_row_of(dlg, "bg"), dlg.COL_VALUE)
        assert item.toolTip() == _settings.EMPTIED_VALUE_TIP

    def test_emptying_a_value_unchecks_the_row(self, dlg):
        _set_value(dlg, "bg", "")
        item = dlg.table.item(_row_of(dlg, "bg"), dlg.COL_INCLUDE)
        assert item.checkState() == Qt.CheckState.Unchecked

    def test_refilling_a_value_clears_the_tooltip(self, dlg):
        _set_value(dlg, "bg", "")
        _set_value(dlg, "bg", "#000")
        item = dlg.table.item(_row_of(dlg, "bg"), dlg.COL_VALUE)
        assert item.toolTip() == ""

    def test_refilled_value_is_collected(self, dlg):
        _set_value(dlg, "bg", "")
        _set_value(dlg, "bg", "#000")
        _set_checked(dlg, "bg", True)
        settings, skipped, _ = dlg.collect()
        assert settings == {"bg": "#000"}
        assert skipped == {}

    def test_locked_row_stays_locked_when_its_text_changes(self, qapp):
        d = _editor({"gone": None})
        try:
            d.table.item(0, d.COL_VALUE).setText("x")
            assert d.is_missing_row(0)
            _, skipped, _ = d.collect()
            assert skipped == {"gone": None}
        finally:
            d.destroy()


class TestLegacyPresetSupport:
    """Presets written before per-key skipping existed must keep working."""

    @pytest.fixture(autouse=True)
    def _reset_globals(self):
        yield
        _settings.EMBED_SETTINGS["enabled"] = False
        _settings.PROJECT_PRESETS.clear()
        _settings.PLUGIN_CONFIG["always_save_to_project"] = False

    def test_legacy_library_lists_its_presets(self, qapp, tmp_path, monkeypatch):
        d, _ = _make_dialog(
            monkeypatch,
            tmp_path,
            library_data={
                "Old": {"a": 1},
                "_PLUGIN_CONFIG": {"always_save_to_project": False},
            },
        )
        try:
            assert d.preset_list.count() == 1
        finally:
            d.destroy()

    def test_legacy_preset_opens_with_every_key_applied(self, qapp):
        d = _editor({"a": 1, "b": "two"})
        try:
            settings, skipped, errors = d.collect()
            assert settings == {"a": 1, "b": "two"}
            assert skipped == {} and errors == []
        finally:
            d.destroy()

    def test_editing_a_legacy_preset_keeps_plugin_config(
        self, qapp, tmp_path, monkeypatch
    ):
        d, data_path = _make_dialog(
            monkeypatch,
            tmp_path,
            library_data={
                "Old": {"a": 1},
                "_PLUGIN_CONFIG": {"always_save_to_project": True},
            },
        )
        try:
            d.preset_list.setCurrentRow(0)

            class _StubEditor:
                def __init__(self, name, settings, skipped=None, existing_names=(),
                             parent=None):
                    self.result_name = name
                    self.result_settings = {"a": 2}
                    self.result_skipped = {}

                def exec(self):
                    return _settings_mod.QDialog.DialogCode.Accepted

            monkeypatch.setattr(_settings, "PresetEditorDialog", _StubEditor)
            d.on_edit()
            saved = json.loads(Path(data_path).read_text(encoding="utf-8"))
            assert saved["Old"] == {"a": 2}
            assert saved["_PLUGIN_CONFIG"] == {"always_save_to_project": True}
        finally:
            d.destroy()

    def test_legacy_preset_with_null_values_is_readable(self, qapp):
        d = _editor({"a": 1, "dead": None})
        try:
            settings, skipped, errors = d.collect()
            assert settings == {"a": 1}
            assert skipped == {"dead": None}
            assert errors == []
        finally:
            d.destroy()

    def test_saved_preset_stays_a_flat_dict(self, qapp, tmp_path, monkeypatch):
        """The on-disk preset shape is unchanged, so older builds can read it."""
        d, data_path = _make_dialog(
            monkeypatch, tmp_path, library_data={"P1": {"x": 1, "y": 2}}
        )
        try:
            d.preset_list.setCurrentRow(0)

            class _StubEditor:
                def __init__(self, name, settings, skipped=None, existing_names=(),
                             parent=None):
                    self.result_name = name
                    self.result_settings = {"x": 1}
                    self.result_skipped = {"y": 2}

                def exec(self):
                    return _settings_mod.QDialog.DialogCode.Accepted

            monkeypatch.setattr(_settings, "PresetEditorDialog", _StubEditor)
            d.on_edit()
            saved = json.loads(Path(data_path).read_text(encoding="utf-8"))
            assert saved["P1"] == {"x": 1}
        finally:
            d.destroy()
