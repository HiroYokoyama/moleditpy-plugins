"""
Unit tests for the Settings Saver plugin.

Covers functions not tested by test_save_load.py:
  - get_plugin_data_path()
  - _get_live_settings(mw)
  - _sync_legacy_settings_alias(mw, settings)
  - on_save_project() / on_load_project(data)
  - enable_project_mode(mw) / disable_project_mode(mw)
  - on_document_reset()
  - load_library() / save_library()
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGIN_PATH = (
    Path(__file__).resolve().parents[1]
    / "plugins" / "Settings_Saver" / "settings_saver.py"
)

with mock_optional_imports():
    MOD = load_plugin(PLUGIN_PATH)


# ---------------------------------------------------------------------------
# get_plugin_data_path
# ---------------------------------------------------------------------------

class TestGetPluginDataPath:
    def test_returns_string(self):
        assert isinstance(MOD.get_plugin_data_path(), str)

    def test_ends_with_settings_saver_json(self):
        assert MOD.get_plugin_data_path().endswith("settings_saver.json")

    def test_parent_directory_exists(self):
        parent = Path(MOD.get_plugin_data_path()).parent
        assert parent.exists()


# ---------------------------------------------------------------------------
# _get_live_settings
# ---------------------------------------------------------------------------

class TestGetLiveSettings:
    def test_returns_settings_dict_when_present(self):
        mw = MagicMock(spec=[])
        mw.init_manager = MagicMock(spec=[])
        expected = {"theme": "dark"}
        mw.init_manager.settings = expected
        assert MOD._get_live_settings(mw) is expected

    def test_returns_none_when_no_init_manager(self):
        mw = MagicMock(spec=[])  # no init_manager attribute
        assert MOD._get_live_settings(mw) is None

    def test_returns_none_when_init_manager_has_no_settings(self):
        mw = MagicMock(spec=[])
        mw.init_manager = MagicMock(spec=[])  # no settings attribute
        assert MOD._get_live_settings(mw) is None


# ---------------------------------------------------------------------------
# _sync_legacy_settings_alias
# ---------------------------------------------------------------------------

class TestSyncLegacySettingsAlias:
    def _make_mw_with_dict_settings(self, initial: dict) -> MagicMock:
        mw = MagicMock(spec=[])
        mw.init_manager = MagicMock(spec=[])
        mw.init_manager.settings = dict(initial)
        return mw

    def test_clears_and_updates_when_different_dict(self):
        mw = self._make_mw_with_dict_settings({"old": 1})
        new_settings = {"new": 2, "extra": 3}
        MOD._sync_legacy_settings_alias(mw, new_settings)
        assert mw.init_manager.settings == new_settings
        assert "old" not in mw.init_manager.settings

    def test_noop_when_same_object(self):
        mw = self._make_mw_with_dict_settings({"key": "val"})
        same = mw.init_manager.settings
        MOD._sync_legacy_settings_alias(mw, same)
        # No mutation: still the same content
        assert mw.init_manager.settings == {"key": "val"}

    def test_noop_when_no_init_manager(self):
        mw = MagicMock(spec=[])
        MOD._sync_legacy_settings_alias(mw, {"x": 1})  # must not raise

    def test_noop_when_settings_is_not_dict(self):
        mw = MagicMock(spec=[])
        mw.init_manager = MagicMock(spec=[])
        mw.init_manager.settings = MagicMock()  # not a dict
        MOD._sync_legacy_settings_alias(mw, {"x": 1})  # must not raise

    def test_exception_during_update_is_swallowed(self):
        class BadDict(dict):
            def update(self, *args, **kwargs):
                raise RuntimeError("oops")

        mw = MagicMock(spec=[])
        mw.init_manager = MagicMock(spec=[])
        mw.init_manager.settings = BadDict({"a": 1})
        MOD._sync_legacy_settings_alias(mw, {"b": 2})  # must not raise


# ---------------------------------------------------------------------------
# on_save_project / on_load_project
# ---------------------------------------------------------------------------

class TestSaveLoadProject:
    def setup_method(self):
        self._orig_ctx = MOD.PLUGIN_CONTEXT
        self._orig_embed = dict(MOD.EMBED_SETTINGS)
        self._orig_presets = dict(MOD.PROJECT_PRESETS)
        self._orig_original = MOD.ORIGINAL_SETTINGS

    def teardown_method(self):
        MOD.PLUGIN_CONTEXT = self._orig_ctx
        MOD.EMBED_SETTINGS.clear()
        MOD.EMBED_SETTINGS.update(self._orig_embed)
        MOD.PROJECT_PRESETS.clear()
        MOD.PROJECT_PRESETS.update(self._orig_presets)
        MOD.ORIGINAL_SETTINGS = self._orig_original

    def test_on_save_project_returns_none_when_not_enabled(self):
        MOD.EMBED_SETTINGS["enabled"] = False
        result = MOD.on_save_project()
        assert result is None

    def test_on_save_project_returns_dict_when_enabled(self):
        ctx = make_context()
        mw = MagicMock(spec=[])
        mw.init_manager = MagicMock(spec=[])
        mw.init_manager.settings = {"theme": "dark"}
        ctx.get_main_window.return_value = mw

        MOD.PLUGIN_CONTEXT = ctx
        MOD.EMBED_SETTINGS["enabled"] = True

        result = MOD.on_save_project()
        assert isinstance(result, dict)
        assert "settings" in result

    def test_on_load_project_empty_dict_does_not_raise(self):
        MOD.PLUGIN_CONTEXT = make_context()
        MOD.on_load_project({})  # no "settings" key → early return

    def test_on_load_project_with_settings_key_does_not_raise(self):
        MOD.PLUGIN_CONTEXT = make_context()
        MOD.on_load_project({"settings": {"theme": "light"}})

    def test_on_load_project_stores_preset(self):
        MOD.PLUGIN_CONTEXT = None  # no context, so only stores
        MOD.on_load_project({"settings": {"color": "blue"}})
        assert "Project Settings" in MOD.PROJECT_PRESETS
        assert MOD.PROJECT_PRESETS["Project Settings"] == {"color": "blue"}

    def test_on_load_project_non_dict_does_not_raise(self):
        MOD.PLUGIN_CONTEXT = make_context()
        MOD.on_load_project(None)  # not a dict → early return
        MOD.on_load_project("bad")

    def test_on_save_project_no_context_uses_active_window(self):
        MOD.PLUGIN_CONTEXT = None
        MOD.EMBED_SETTINGS["enabled"] = True
        result = MOD.on_save_project()  # falls back to QApplication.activeWindow()
        assert result is None

    def test_on_load_project_applies_settings_full_path(self):
        mw = SimpleNamespace(init_manager=SimpleNamespace(settings={"orig": 1}))
        ctx = make_context()
        ctx.get_main_window.return_value = mw
        MOD.PLUGIN_CONTEXT = ctx
        MOD.ORIGINAL_SETTINGS = None
        MOD.on_load_project({"settings": {"theme": "dark"}})
        assert MOD.ORIGINAL_SETTINGS == {"orig": 1}
        assert mw.init_manager.settings["theme"] == "dark"
        assert MOD.EMBED_SETTINGS["enabled"] is True
        ctx.show_status_message.assert_called_once()

    def test_on_load_project_exception_swallowed(self):
        class BadDict(dict):
            def update(self, *a, **k):
                raise RuntimeError("boom")

        mw = SimpleNamespace(init_manager=SimpleNamespace(settings=BadDict({"a": 1})))
        ctx = make_context()
        ctx.get_main_window.return_value = mw
        MOD.PLUGIN_CONTEXT = ctx
        MOD.on_load_project({"settings": {"x": 1}})  # must not raise


# ---------------------------------------------------------------------------
# enable_project_mode / disable_project_mode
# ---------------------------------------------------------------------------

class TestProjectMode:
    def setup_method(self):
        self._orig_embed = dict(MOD.EMBED_SETTINGS)
        self._orig_original = MOD.ORIGINAL_SETTINGS

    def teardown_method(self):
        MOD.EMBED_SETTINGS.clear()
        MOD.EMBED_SETTINGS.update(self._orig_embed)
        MOD.ORIGINAL_SETTINGS = self._orig_original

    def _bare_mw(self):
        """MagicMock with no spec so attribute access never raises."""
        return MagicMock()

    def test_enable_sets_embed_enabled_true(self):
        MOD.EMBED_SETTINGS["enabled"] = False
        MOD.enable_project_mode(self._bare_mw())
        assert MOD.EMBED_SETTINGS["enabled"] is True

    def test_enable_does_not_raise_with_mock_mw(self):
        MOD.enable_project_mode(self._bare_mw())  # must not raise

    def test_disable_sets_embed_enabled_false(self):
        MOD.EMBED_SETTINGS["enabled"] = True
        MOD.disable_project_mode(self._bare_mw(), restore_content=False)
        assert MOD.EMBED_SETTINGS["enabled"] is False

    def test_disable_does_not_raise_with_mock_mw(self):
        MOD.disable_project_mode(self._bare_mw(), restore_content=False)

    def test_enable_then_disable_round_trip(self):
        mw = self._bare_mw()
        MOD.enable_project_mode(mw)
        assert MOD.EMBED_SETTINGS["enabled"] is True
        MOD.disable_project_mode(mw, restore_content=False)
        assert MOD.EMBED_SETTINGS["enabled"] is False

    def test_enable_project_mode_full_paths(self):
        """Exercise dirty-flag reset, initial_settings alias, and save_settings monkeypatch."""
        mw = SimpleNamespace(
            settings_dirty=None,
            initial_settings=None,
            save_settings=lambda: None,
            init_manager=SimpleNamespace(
                settings={"a": 1}, settings_dirty=None, save_settings=lambda: "orig"
            ),
        )
        MOD.enable_project_mode(mw)
        assert mw.init_manager.settings_dirty is False
        assert mw.initial_settings == {"a": 1}
        assert hasattr(mw, "_original_save_settings")
        assert mw.init_manager.save_settings() is None  # proxy blocks the call

    def test_disable_project_mode_restores_original(self):
        """Exercise save_settings restoration and content restoration from ORIGINAL_SETTINGS."""
        mw = SimpleNamespace(
            initial_settings={"old": True},
            init_manager=SimpleNamespace(settings={"live": 1}, save_settings=lambda: "orig"),
        )
        mw._original_save_settings = lambda: "orig"
        MOD.ORIGINAL_SETTINGS = {"restored": 2}
        MOD.disable_project_mode(mw, restore_content=True)
        assert not hasattr(mw, "_original_save_settings")
        assert mw.init_manager.save_settings() == "orig"
        assert mw.initial_settings == {"restored": 2}
        assert mw.init_manager.settings == {"restored": 2}
        assert MOD.ORIGINAL_SETTINGS is None

    def test_disable_project_mode_restore_content_exception_swallowed(self):
        class BadDict(dict):
            def clear(self):
                raise RuntimeError("boom")

        mw = SimpleNamespace(init_manager=SimpleNamespace(settings=BadDict({"a": 1})))
        MOD.ORIGINAL_SETTINGS = {"x": 1}
        MOD.disable_project_mode(mw, restore_content=True)  # must not raise
        assert MOD.ORIGINAL_SETTINGS is None


# ---------------------------------------------------------------------------
# on_document_reset
# ---------------------------------------------------------------------------

class TestOnDocumentReset:
    def setup_method(self):
        self._orig_ctx = MOD.PLUGIN_CONTEXT
        self._orig_embed = dict(MOD.EMBED_SETTINGS)
        self._orig_presets = dict(MOD.PROJECT_PRESETS)
        self._orig_config = dict(MOD.PLUGIN_CONFIG)
        self._orig_original = MOD.ORIGINAL_SETTINGS

    def teardown_method(self):
        MOD.PLUGIN_CONTEXT = self._orig_ctx
        MOD.EMBED_SETTINGS.clear()
        MOD.EMBED_SETTINGS.update(self._orig_embed)
        MOD.PROJECT_PRESETS.clear()
        MOD.PROJECT_PRESETS.update(self._orig_presets)
        MOD.PLUGIN_CONFIG.clear()
        MOD.PLUGIN_CONFIG.update(self._orig_config)
        MOD.ORIGINAL_SETTINGS = self._orig_original

    def test_does_not_raise_with_no_context(self):
        MOD.PLUGIN_CONTEXT = None
        MOD.on_document_reset()

    def test_does_not_raise_with_mock_context(self):
        MOD.PLUGIN_CONTEXT = make_context()
        MOD.on_document_reset()

    def test_clears_project_presets(self):
        MOD.PLUGIN_CONTEXT = None
        MOD.PROJECT_PRESETS["Project Settings"] = {"x": 1}
        MOD.on_document_reset()
        assert MOD.PROJECT_PRESETS == {}

    def test_disables_embed_when_always_save_false(self):
        MOD.PLUGIN_CONTEXT = None
        MOD.EMBED_SETTINGS["enabled"] = True
        MOD.PLUGIN_CONFIG["always_save_to_project"] = False
        MOD.on_document_reset()
        assert MOD.EMBED_SETTINGS["enabled"] is False

    def test_reenables_project_mode_when_always_save_true(self):
        mw = SimpleNamespace(init_manager=SimpleNamespace(settings={}))
        ctx = make_context()
        ctx.get_main_window.return_value = mw
        MOD.PLUGIN_CONTEXT = ctx
        MOD.PLUGIN_CONFIG["always_save_to_project"] = True
        MOD.on_document_reset()
        assert MOD.EMBED_SETTINGS["enabled"] is True


class TestSettingsSaver:
    def _patch_path(self, tmp_path, monkeypatch) -> str:
        data_path = str(tmp_path / "settings_saver.json")
        monkeypatch.setattr(MOD, "get_plugin_data_path", lambda: data_path)
        return data_path

    def test_round_trip(self, tmp_path, monkeypatch):
        self._patch_path(tmp_path, monkeypatch)
        payload = {
            "preset_A": {"opacity": 0.5},
            "_PLUGIN_CONFIG": {"always_save_to_project": True},
        }
        MOD.save_library(payload)
        loaded = MOD.load_library()
        assert loaded["preset_A"]["opacity"] == pytest.approx(0.5)
        assert loaded["_PLUGIN_CONFIG"]["always_save_to_project"] is True

    def test_load_missing_returns_empty_dict(self, tmp_path, monkeypatch):
        self._patch_path(tmp_path, monkeypatch)
        assert MOD.load_library() == {}

    def test_save_overwrites_existing(self, tmp_path, monkeypatch):
        self._patch_path(tmp_path, monkeypatch)
        MOD.save_library({"old": 1})
        MOD.save_library({"new": 2})
        assert MOD.load_library() == {"new": 2}

    def test_get_plugin_data_path_returns_json(self):
        path = MOD.get_plugin_data_path()
        assert isinstance(path, str)
        assert path.endswith(".json")
        assert os.path.basename(path) == MOD.SETTINGS_FILENAME

    def test_load_corrupted_json_returns_empty_dict(self, tmp_path, monkeypatch):
        data_path = self._patch_path(tmp_path, monkeypatch)
        Path(data_path).write_text("{not valid json", encoding="utf-8")
        assert MOD.load_library() == {}

    def test_save_library_write_error_shows_critical(self, tmp_path, monkeypatch):
        bad_path = str(tmp_path / "missing_dir" / "file.json")
        monkeypatch.setattr(MOD, "get_plugin_data_path", lambda: bad_path)
        parent = MagicMock()
        with patch.object(MOD.QMessageBox, "critical") as crit:
            MOD.save_library({"a": 1}, parent)
        crit.assert_called_once()

    def test_save_library_write_error_no_parent_silent(self, tmp_path, monkeypatch):
        bad_path = str(tmp_path / "missing_dir2" / "file.json")
        monkeypatch.setattr(MOD, "get_plugin_data_path", lambda: bad_path)
        MOD.save_library({"a": 1})  # no parent_window -> silently swallowed


# ---------------------------------------------------------------------------
# initialize
# ---------------------------------------------------------------------------


class TestInitialize:
    def setup_method(self):
        self._orig_ctx = MOD.PLUGIN_CONTEXT
        self._orig_config = dict(MOD.PLUGIN_CONFIG)
        self._orig_embed = dict(MOD.EMBED_SETTINGS)

    def teardown_method(self):
        MOD.PLUGIN_CONTEXT = self._orig_ctx
        MOD.PLUGIN_CONFIG.clear()
        MOD.PLUGIN_CONFIG.update(self._orig_config)
        MOD.EMBED_SETTINGS.clear()
        MOD.EMBED_SETTINGS.update(self._orig_embed)

    def test_registers_handlers_and_menu(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            MOD, "get_plugin_data_path", lambda: str(tmp_path / "settings_saver.json")
        )
        ctx = make_context()
        MOD.initialize(ctx)
        ctx.add_menu_action.assert_called_once()
        ctx.register_save_handler.assert_called_once_with(MOD.on_save_project)
        ctx.register_load_handler.assert_called_once_with(MOD.on_load_project)
        ctx.register_document_reset_handler.assert_called_once_with(MOD.on_document_reset)

    def test_loads_config_and_enables_project_mode(self, tmp_path, monkeypatch):
        data_path = str(tmp_path / "settings_saver.json")
        monkeypatch.setattr(MOD, "get_plugin_data_path", lambda: data_path)
        Path(data_path).write_text(
            json.dumps({"_PLUGIN_CONFIG": {"always_save_to_project": True}}),
            encoding="utf-8",
        )
        ctx = make_context()
        MOD.initialize(ctx)
        assert MOD.PLUGIN_CONFIG["always_save_to_project"] is True
        assert MOD.EMBED_SETTINGS["enabled"] is True


# ---------------------------------------------------------------------------
# apply_settings_hot / refresh_loaded_scene
# ---------------------------------------------------------------------------


class _FakeInitManager:
    def __init__(self, settings):
        self.settings = settings
        self.settings_dirty = None
        self.cpk_calls = 0

    def update_cpk_colors_from_settings(self):
        self.cpk_calls += 1


class _FakeView3D:
    def __init__(self):
        self.applied = False

    def apply_3d_settings(self):
        self.applied = True


class _FakeScene:
    def __init__(self):
        self.background_calls = []
        self._items = []

    def setBackgroundBrush(self, color):
        self.background_calls.append(color)

    def items(self):
        return self._items


class TestApplySettingsHot:
    def setup_method(self):
        self._orig_ctx = MOD.PLUGIN_CONTEXT
        self._orig_embed = dict(MOD.EMBED_SETTINGS)

    def teardown_method(self):
        MOD.PLUGIN_CONTEXT = self._orig_ctx
        MOD.EMBED_SETTINGS.clear()
        MOD.EMBED_SETTINGS.update(self._orig_embed)

    def _mw(self):
        from types import SimpleNamespace

        return SimpleNamespace(
            init_manager=_FakeInitManager({"background_color_2d": "#112233"}),
            view_3d_manager=_FakeView3D(),
            scene=_FakeScene(),
        )

    def test_no_settings_returns_early_without_raising(self):
        mw = MagicMock(spec=[])
        MOD.PLUGIN_CONTEXT = None
        MOD.apply_settings_hot(mw)  # must not raise

    def test_applies_3d_settings_and_cpk_colors(self):
        MOD.PLUGIN_CONTEXT = None
        MOD.EMBED_SETTINGS["enabled"] = False
        mw = self._mw()
        MOD.apply_settings_hot(mw)
        assert mw.view_3d_manager.applied is True
        assert mw.init_manager.cpk_calls == 1

    def test_dirty_flag_set_true_when_project_mode_disabled(self):
        MOD.PLUGIN_CONTEXT = None
        MOD.EMBED_SETTINGS["enabled"] = False
        mw = self._mw()
        MOD.apply_settings_hot(mw)
        assert mw.init_manager.settings_dirty is True

    def test_dirty_flag_cleared_when_project_mode_enabled(self):
        MOD.PLUGIN_CONTEXT = None
        MOD.EMBED_SETTINGS["enabled"] = True
        mw = self._mw()
        MOD.apply_settings_hot(mw)
        assert mw.init_manager.settings_dirty is False

    def test_updates_2d_scene_background(self):
        MOD.PLUGIN_CONTEXT = None
        MOD.EMBED_SETTINGS["enabled"] = False
        mw = self._mw()
        MOD.apply_settings_hot(mw)
        assert len(mw.scene.background_calls) == 1

    def test_dirty_flag_assignment_exception_swallowed(self):
        class _RaisingIM:
            def __init__(self, settings):
                self.settings = settings

            @property
            def settings_dirty(self):
                return None

            @settings_dirty.setter
            def settings_dirty(self, value):
                raise RuntimeError("boom")

        MOD.PLUGIN_CONTEXT = None
        MOD.EMBED_SETTINGS["enabled"] = False
        mw = SimpleNamespace(init_manager=_RaisingIM({"a": 1}))
        MOD.apply_settings_hot(mw)  # must not raise

    def test_3d_settings_exception_swallowed(self):
        view3d = MagicMock()
        view3d.apply_3d_settings.side_effect = RuntimeError("boom")
        MOD.PLUGIN_CONTEXT = None
        mw = SimpleNamespace(init_manager=_FakeInitManager({"a": 1}), view_3d_manager=view3d)
        MOD.apply_settings_hot(mw)  # must not raise

    def test_cpk_exception_swallowed(self):
        im = _FakeInitManager({"a": 1})
        im.update_cpk_colors_from_settings = MagicMock(side_effect=RuntimeError("boom"))
        MOD.PLUGIN_CONTEXT = None
        mw = SimpleNamespace(init_manager=im)
        MOD.apply_settings_hot(mw)  # must not raise

    def test_refresh_ui_inner_exception_swallowed(self):
        widget = MagicMock()
        widget.refresh_ui.side_effect = RuntimeError("boom")
        MOD.PLUGIN_CONTEXT = None
        mw = self._mw()
        with patch.object(MOD.QApplication, "topLevelWidgets", return_value=[widget]):
            MOD.apply_settings_hot(mw)  # must not raise

    def test_top_level_widgets_outer_exception_swallowed(self):
        MOD.PLUGIN_CONTEXT = None
        mw = self._mw()
        with patch.object(
            MOD.QApplication, "topLevelWidgets", side_effect=RuntimeError("boom")
        ):
            MOD.apply_settings_hot(mw)  # must not raise

    def test_draw_molecule_exception_swallowed(self):
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        ctx.draw_molecule_3d.side_effect = RuntimeError("boom")
        MOD.PLUGIN_CONTEXT = ctx
        mw = self._mw()
        MOD.apply_settings_hot(mw)  # must not raise

    def test_2d_scene_exception_swallowed(self):
        scene = MagicMock()
        scene.setBackgroundBrush.side_effect = RuntimeError("boom")
        MOD.PLUGIN_CONTEXT = None
        mw = SimpleNamespace(init_manager=_FakeInitManager({"a": 1}), scene=scene)
        MOD.apply_settings_hot(mw)  # must not raise

    def test_undo_checkpoint_exception_swallowed(self):
        ctx = make_context()
        ctx.push_undo_checkpoint.side_effect = RuntimeError("boom")
        MOD.PLUGIN_CONTEXT = ctx
        mw = self._mw()
        MOD.apply_settings_hot(mw)  # must not raise


class _FakeScene2:
    def __init__(self):
        self.updated = False
        self._views = []

    def update(self):
        self.updated = True

    def views(self):
        return self._views


class TestRefreshLoadedScene:
    def test_immediate_refresh_updates_scene(self):
        from types import SimpleNamespace

        mw = SimpleNamespace(scene=_FakeScene2())
        MOD.refresh_loaded_scene(mw, defer=False)
        assert mw.scene.updated is True

    def test_deferred_refresh_uses_qtimer_singleshot(self):
        from types import SimpleNamespace

        mw = SimpleNamespace(scene=_FakeScene2())
        MOD.QTimer.singleShot.reset_mock()
        MOD.refresh_loaded_scene(mw, defer=True)
        MOD.QTimer.singleShot.assert_called_once()
        args = MOD.QTimer.singleShot.call_args[0]
        assert args[0] == 0
        assert callable(args[1])
        # scene not updated synchronously (deferred via QTimer)
        assert mw.scene.updated is False

    def test_views_error_swallowed(self):
        class _SceneNoViewsError(_FakeScene2):
            def views(self):
                raise RuntimeError("boom")

        mw = SimpleNamespace(scene=_SceneNoViewsError())
        MOD.refresh_loaded_scene(mw, defer=False)  # must not raise

    def test_view2d_viewport_error_swallowed(self):
        view2d = MagicMock()
        view2d.viewport.side_effect = RuntimeError("boom")
        mw = SimpleNamespace(init_manager=SimpleNamespace(view_2d=view2d))
        MOD.refresh_loaded_scene(mw, defer=False)  # must not raise

    def test_edit_3d_manager_error_swallowed(self):
        e3d = MagicMock()
        e3d.update_2d_measurement_labels.side_effect = RuntimeError("boom")
        mw = SimpleNamespace(edit_3d_manager=e3d)
        MOD.refresh_loaded_scene(mw, defer=False)  # must not raise

    def test_draw_molecule_error_swallowed(self):
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        ctx.draw_molecule_3d.side_effect = RuntimeError("boom")
        MOD.PLUGIN_CONTEXT = ctx
        try:
            mw = SimpleNamespace()
            MOD.refresh_loaded_scene(mw, defer=False)  # must not raise
        finally:
            MOD.PLUGIN_CONTEXT = None

    def test_outer_exception_swallowed(self):
        class _RaisingScene:
            def update(self):
                raise RuntimeError("boom")

        mw = SimpleNamespace(scene=_RaisingScene())
        MOD.refresh_loaded_scene(mw, defer=False)  # must not raise


# ---------------------------------------------------------------------------
# SettingsSaverDialog: real-instance behavior tests (QDialog base is mocked)
# ---------------------------------------------------------------------------


class _FakeListItem:
    def __init__(self, text, is_project=False):
        self._text = text
        self._is_project = is_project

    def text(self):
        return self._text

    def data(self, role):
        return "project" if self._is_project else None


class _FakePresetList:
    def __init__(self, items=None):
        self._selected = items or []
        self.cleared = False
        self.added = []
        self.find_items_result = []
        self.set_current = []

    def selectedItems(self):
        return self._selected

    def clear(self):
        self.cleared = True
        self.added = []

    def addItem(self, item):
        self.added.append(item)

    def findItems(self, name, flag):
        return self.find_items_result

    def setCurrentItem(self, item):
        self.set_current.append(item)


def _extract_dialog_method(method_name):
    """
    Extract a SettingsSaverDialog method via AST and exec it with MOD's real
    module namespace as globals (so Qt/QMessageBox/save_library/etc. resolve
    to the same mocked objects the module itself uses).
    """
    return extract_function(
        PLUGIN_PATH, "SettingsSaverDialog", method_name, extra_globals=dict(vars(MOD))
    )


def _fake_dialog_self(mw=None, library=None):
    return SimpleNamespace(
        library=library if library is not None else {},
        preset_list=_FakePresetList(),
        main_window=mw if mw is not None else MagicMock(),
        context=make_context(),
        btn_set_global=MagicMock(),
        refresh_list=MagicMock(),
        is_project_preset=lambda item: item.data(None) == "project",
    )


class TestSettingsSaverDialogMethods:
    @pytest.fixture(autouse=True)
    def _patch_data_path(self, tmp_path, monkeypatch):
        # Every dialog method here may transitively call save_library()/
        # load_library(); always redirect them away from the real plugin
        # directory to avoid polluting it with test artifacts.
        data_path = str(tmp_path / "settings_saver.json")
        monkeypatch.setattr(MOD, "get_plugin_data_path", lambda: data_path)

    def test_get_selected_name_none_when_nothing_selected(self):
        fn = _extract_dialog_method("get_selected_name")
        self_ = _fake_dialog_self()
        assert fn(self_) is None

    def test_get_selected_name_returns_text(self):
        fn = _extract_dialog_method("get_selected_name")
        self_ = _fake_dialog_self()
        self_.preset_list = _FakePresetList([_FakeListItem("Preset A")])
        assert fn(self_) == "Preset A"

    def test_is_project_preset_true_for_project_item(self):
        fn = _extract_dialog_method("is_project_preset")
        self_ = _fake_dialog_self()
        assert fn(self_, _FakeListItem("x", is_project=True)) is True

    def test_is_project_preset_false_for_library_item(self):
        fn = _extract_dialog_method("is_project_preset")
        self_ = _fake_dialog_self()
        assert fn(self_, _FakeListItem("x", is_project=False)) is False

    def test_on_always_save_toggled_persists_config(self):
        fn = _extract_dialog_method("on_always_save_toggled")
        self_ = _fake_dialog_self()
        try:
            fn(self_, True)
            assert MOD.PLUGIN_CONFIG["always_save_to_project"] is True
            assert MOD.load_library()["_PLUGIN_CONFIG"]["always_save_to_project"] is True
        finally:
            MOD.PLUGIN_CONFIG["always_save_to_project"] = False

    def test_on_embed_toggled_true_enables_project_mode(self):
        fn = _extract_dialog_method("on_embed_toggled")
        mw = MagicMock()
        self_ = _fake_dialog_self(mw=mw)
        try:
            fn(self_, True)
            assert MOD.EMBED_SETTINGS["enabled"] is True
            self_.btn_set_global.setEnabled.assert_called_with(True)
        finally:
            MOD.EMBED_SETTINGS["enabled"] = False

    def test_on_embed_toggled_false_disables_and_marks_dirty(self):
        fn = _extract_dialog_method("on_embed_toggled")
        mw = MagicMock()
        self_ = _fake_dialog_self(mw=mw)
        fn(self_, False)
        assert MOD.EMBED_SETTINGS["enabled"] is False
        assert mw.edit_actions_manager.settings_dirty is True

    def test_on_save_adds_new_preset(self):
        fn = _extract_dialog_method("on_save")
        mw = MagicMock()
        mw.init_manager.settings = {"a": 1}
        self_ = _fake_dialog_self(mw=mw)
        with patch.object(MOD.QInputDialog, "getText", return_value=("MyPreset", True)):
            fn(self_)
        assert self_.library.get("MyPreset") == {"a": 1}
        self_.refresh_list.assert_called_once()

    def test_on_save_rejects_underscore_prefixed_name(self):
        fn = _extract_dialog_method("on_save")
        self_ = _fake_dialog_self()
        with patch.object(
            MOD.QInputDialog, "getText", return_value=("_bad", True)
        ), patch.object(MOD.QMessageBox, "warning") as warn:
            fn(self_)
        warn.assert_called_once()
        assert "_bad" not in self_.library

    def test_on_save_cancelled_does_nothing(self):
        fn = _extract_dialog_method("on_save")
        self_ = _fake_dialog_self()
        with patch.object(MOD.QInputDialog, "getText", return_value=("", False)):
            fn(self_)
        assert self_.library == {}

    def test_on_delete_removes_library_preset(self):
        fn = _extract_dialog_method("on_delete")
        self_ = _fake_dialog_self(library={"MyPreset": {"a": 1}})
        self_.preset_list = _FakePresetList([_FakeListItem("MyPreset")])
        with patch.object(
            MOD.QMessageBox, "question", return_value=MOD.QMessageBox.StandardButton.Yes
        ):
            fn(self_)
        assert "MyPreset" not in self_.library

    def test_on_delete_project_preset_blocked(self):
        fn = _extract_dialog_method("on_delete")
        self_ = _fake_dialog_self()
        self_.preset_list = _FakePresetList(
            [_FakeListItem("Project Settings", is_project=True)]
        )
        with patch.object(MOD.QMessageBox, "information") as info:
            fn(self_)
        info.assert_called_once()

    def test_on_load_applies_library_preset(self):
        fn = _extract_dialog_method("on_load")
        mw = MagicMock()
        self_ = _fake_dialog_self(mw=mw, library={"P1": {"x": 1}})
        self_.preset_list = _FakePresetList([_FakeListItem("P1")])
        fn(self_)
        mw.init_manager.settings.update.assert_called_once_with({"x": 1})

    def test_on_load_missing_preset_warns(self):
        fn = _extract_dialog_method("on_load")
        self_ = _fake_dialog_self()
        self_.preset_list = _FakePresetList([_FakeListItem("DoesNotExist")])
        with patch.object(MOD.QMessageBox, "warning") as warn:
            fn(self_)
        warn.assert_called_once()

    def test_on_export_single_writes_selected_preset(self, tmp_path):
        fn = _extract_dialog_method("on_export_single")
        self_ = _fake_dialog_self(library={"P1": {"x": 1}})
        self_.preset_list = _FakePresetList([_FakeListItem("P1")])
        out = tmp_path / "out.json"
        with patch.object(
            MOD.QFileDialog, "getSaveFileName", return_value=(str(out), "")
        ):
            fn(self_)
        assert json.loads(out.read_text(encoding="utf-8")) == {"x": 1}

    def test_on_export_all_excludes_underscore_keys(self, tmp_path):
        fn = _extract_dialog_method("on_export_all")
        self_ = _fake_dialog_self(
            library={
                "P1": {"x": 1},
                "_PLUGIN_CONFIG": {"always_save_to_project": True},
            }
        )
        out = tmp_path / "all.json"
        with patch.object(
            MOD.QFileDialog, "getSaveFileName", return_value=(str(out), "")
        ):
            fn(self_)
        data = json.loads(out.read_text(encoding="utf-8"))
        assert "P1" in data
        assert "_PLUGIN_CONFIG" not in data

    def test_on_import_library_format_skips_underscore_keys(self, tmp_path):
        fn = _extract_dialog_method("on_import")
        src = tmp_path / "in.json"
        src.write_text(
            json.dumps({"P1": {"x": 1}, "_PLUGIN_CONFIG": {"a": 1}}),
            encoding="utf-8",
        )
        self_ = _fake_dialog_self()
        with patch.object(
            MOD.QFileDialog, "getOpenFileName", return_value=(str(src), "")
        ):
            fn(self_)
        assert self_.library.get("P1") == {"x": 1}
        assert "_PLUGIN_CONFIG" not in self_.library

    def test_on_import_single_preset_format(self, tmp_path):
        fn = _extract_dialog_method("on_import")
        src = tmp_path / "single.json"
        src.write_text(json.dumps({"opacity": 0.5}), encoding="utf-8")
        self_ = _fake_dialog_self()
        with patch.object(
            MOD.QFileDialog, "getOpenFileName", return_value=(str(src), "")
        ), patch.object(MOD.QInputDialog, "getText", return_value=("Imported1", True)):
            fn(self_)
        assert self_.library.get("Imported1") == {"opacity": 0.5}

    def test_on_import_invalid_json_shows_critical(self, tmp_path):
        fn = _extract_dialog_method("on_import")
        src = tmp_path / "bad.json"
        src.write_text("not valid json", encoding="utf-8")
        self_ = _fake_dialog_self()
        with patch.object(
            MOD.QFileDialog, "getOpenFileName", return_value=(str(src), "")
        ), patch.object(MOD.QMessageBox, "critical") as crit:
            fn(self_)
        crit.assert_called_once()

    def test_on_save_as_global_default_cancelled_noop(self):
        fn = _extract_dialog_method("on_save_as_global_default")
        self_ = _fake_dialog_self()
        with patch.object(
            MOD.QMessageBox, "question", return_value=MOD.QMessageBox.StandardButton.No
        ):
            fn(self_)  # must not raise

    def test_on_save_as_global_default_no_save_func_warns(self):
        fn = _extract_dialog_method("on_save_as_global_default")
        mw = MagicMock(spec=["init_manager"])
        mw.init_manager = MagicMock(spec=["settings"])
        self_ = _fake_dialog_self(mw=mw)
        with patch.object(
            MOD.QMessageBox, "question", return_value=MOD.QMessageBox.StandardButton.Yes
        ), patch.object(MOD.QMessageBox, "warning") as warn:
            fn(self_)
        warn.assert_called_once()

    def test_on_save_as_global_default_executes_save_func(self):
        fn = _extract_dialog_method("on_save_as_global_default")
        # spec-limited so hasattr(mw, "_original_save_settings") is False,
        # forcing the mw.init_manager.save_settings branch.
        mw = MagicMock(spec=["init_manager", "initial_settings"])
        mw.init_manager = MagicMock(spec=["settings", "save_settings"])
        mw.init_manager.settings = {"a": 1}
        self_ = _fake_dialog_self(mw=mw)
        with patch.object(
            MOD.QMessageBox, "question", return_value=MOD.QMessageBox.StandardButton.Yes
        ), patch.object(MOD.QMessageBox, "information") as info:
            fn(self_)
        mw.init_manager.save_settings.assert_called_once()
        info.assert_called_once()

    def test_refresh_list_sorts_and_marks_project_presets(self):
        fn = _extract_dialog_method("refresh_list")
        self_ = _fake_dialog_self(library={"Zeta": {}, "Alpha": {}})
        self_.preset_list = _FakePresetList()
        MOD.PROJECT_PRESETS["Project Settings"] = {}
        try:
            fn(self_)
            assert self_.preset_list.cleared is True
            assert len(self_.preset_list.added) == 3
        finally:
            MOD.PROJECT_PRESETS.clear()


class TestOpenManager:
    def test_open_manager_does_not_raise(self, tmp_path, monkeypatch):
        data_path = str(tmp_path / "settings_saver.json")
        monkeypatch.setattr(MOD, "get_plugin_data_path", lambda: data_path)
        ctx = make_context()
        MOD.open_manager(ctx)  # must not raise
