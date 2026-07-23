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
        mw.edit_actions_manager.settings_dirty = True
        d, _ = _make_dialog(monkeypatch, tmp_path, mw=mw, embed_enabled=True)
        try:
            assert mw.edit_actions_manager.settings_dirty is False
            assert mw.init_manager.initial_settings == {"a": 1}
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
            assert mw.edit_actions_manager.settings_dirty is True
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
