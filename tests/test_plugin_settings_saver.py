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

import os
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

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
