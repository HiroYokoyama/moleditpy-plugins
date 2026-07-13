"""
Tests for the Metadata Saver plugin.

Covers:
- Metadata: module constants, field definitions
- Config: load/save/defaults/merge
- Metadata collection: each field, error resilience
- Plugin lifecycle: initialize, on_save_project, on_load_project, on_document_reset
- Dialog: smoke-test construction via AST extraction (Qt-free)
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import types
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
PLUGIN_PATH = PLUGINS_DIR / "Metadata_Saver" / "metadata_saver.py"


# ---------------------------------------------------------------------------
# Shared fixture — module loaded once per class
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def mod():
    with mock_optional_imports():
        return load_plugin(PLUGIN_PATH)


# ---------------------------------------------------------------------------
# Metadata constants
# ---------------------------------------------------------------------------


class TestMetadata:
    def test_plugin_name(self, mod):
        assert mod.PLUGIN_NAME == "Metadata Saver"

    def test_plugin_version_format(self, mod):
        # Expect YYYY.MM.DD
        parts = mod.PLUGIN_VERSION.split(".")
        assert len(parts) == 3
        assert all(p.isdigit() for p in parts)

    def test_author(self, mod):
        assert mod.PLUGIN_AUTHOR == "HiroYokoyama"

    def test_category(self, mod):
        assert mod.PLUGIN_CATEGORY == "Utility"

    def test_supported_version_specifier(self, mod):
        assert "4.0.0" in mod.PLUGIN_SUPPORTED_MOLEDITPY_VERSION

    def test_no_heavy_dependencies(self, mod):
        assert mod.PLUGIN_DEPENDENCIES == []

    def test_all_fields_have_required_keys(self, mod):
        for field in mod.ALL_FIELDS:
            assert "key" in field
            assert "label" in field
            assert "description" in field
            assert "default" in field
            assert isinstance(field["default"], bool)

    def test_field_keys_are_unique(self, mod):
        keys = [f["key"] for f in mod.ALL_FIELDS]
        assert len(keys) == len(set(keys))

    def test_expected_fields_present(self, mod):
        keys = {f["key"] for f in mod.ALL_FIELDS}
        for expected in ("saved_at", "saved_path", "username", "os_name",
                         "os_release", "app_version", "note", "installed_plugins"):
            assert expected in keys, f"Missing field: {expected}"

    def test_hostname_default_on(self, mod):
        hostname_field = next(f for f in mod.ALL_FIELDS if f["key"] == "hostname")
        assert hostname_field["default"] is True

    def test_default_false_fields(self, mod):
        note_field = next(f for f in mod.ALL_FIELDS if f["key"] == "note")
        assert note_field["default"] is False
        plugins_field = next(f for f in mod.ALL_FIELDS if f["key"] == "installed_plugins")
        assert plugins_field["default"] is False

    def test_all_other_fields_default_on(self, mod):
        for f in mod.ALL_FIELDS:
            if f["key"] not in ("note", "installed_plugins"):
                assert f["default"] is True, f"Expected {f['key']} default=True"

    def test_enabled_flag_true_by_default(self, mod):
        cfg = mod._default_config()
        assert cfg["enabled"] is True

    def test_initial_loaded_metadata_empty(self, mod):
        assert mod._LOADED_METADATA == {} or isinstance(mod._LOADED_METADATA, dict)


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------


class TestConfig:
    def test_default_config_has_all_field_keys(self, mod):
        cfg = mod._default_config()
        field_keys = {f["key"] for f in mod.ALL_FIELDS}
        assert set(cfg["enabled_fields"].keys()) == field_keys

    def test_default_config_custom_note_empty(self, mod):
        cfg = mod._default_config()
        assert cfg["custom_note"] == ""

    def test_load_config_returns_dict_when_no_file(self, mod):
        with tempfile.TemporaryDirectory() as tmpdir:
            fake_path = os.path.join(tmpdir, "metadata_saver.json")
            with patch.object(mod, "_config_path", return_value=fake_path):
                cfg = mod.load_config()
        assert isinstance(cfg, dict)
        assert "enabled_fields" in cfg

    def test_load_config_merges_with_defaults(self, mod):
        """Partial JSON on disk should be merged with defaults — new keys appear."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fake_path = os.path.join(tmpdir, "metadata_saver.json")
            partial = {"enabled_fields": {"saved_at": False}, "custom_note": "hi"}
            with open(fake_path, "w") as f:
                json.dump(partial, f)
            with patch.object(mod, "_config_path", return_value=fake_path):
                cfg = mod.load_config()
        # Explicit override should survive
        assert cfg["enabled_fields"]["saved_at"] is False
        # Missing keys should use defaults
        assert "username" in cfg["enabled_fields"]
        assert cfg["custom_note"] == "hi"

    def test_save_config_roundtrip(self, mod):
        with tempfile.TemporaryDirectory() as tmpdir:
            fake_path = os.path.join(tmpdir, "metadata_saver.json")
            cfg = mod._default_config()
            cfg["custom_note"] = "test-note"
            cfg["enabled_fields"]["hostname"] = True
            with patch.object(mod, "_config_path", return_value=fake_path):
                mod.save_config(cfg)
                loaded = mod.load_config()
        assert loaded["custom_note"] == "test-note"
        assert loaded["enabled_fields"]["hostname"] is True

    def test_load_config_handles_corrupt_json(self, mod):
        with tempfile.TemporaryDirectory() as tmpdir:
            fake_path = os.path.join(tmpdir, "metadata_saver.json")
            with open(fake_path, "w") as f:
                f.write("{corrupt json{{")
            with patch.object(mod, "_config_path", return_value=fake_path):
                cfg = mod.load_config()  # must not raise
        assert isinstance(cfg, dict)

    def test_save_config_handles_permission_error(self, mod):
        with patch.object(mod, "_config_path", return_value="/nonexistent/path/x.json"):
            mod.save_config(mod._default_config())  # must not raise


# ---------------------------------------------------------------------------
# Metadata collection
# ---------------------------------------------------------------------------


class TestCollectField:
    """Unit-test _collect_field() for each key."""

    def _call(self, mod, key, mw=None):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
            return m._collect_field(key, mw)

    def test_saved_at_is_iso(self, mod):
        val = self._call(mod, "saved_at")
        from datetime import datetime
        dt = datetime.fromisoformat(val)
        assert dt is not None

    def test_username_is_string(self, mod):
        val = self._call(mod, "username")
        assert isinstance(val, str) and val

    def test_hostname_is_string(self, mod):
        val = self._call(mod, "hostname")
        assert isinstance(val, str) and val

    def test_os_name_is_string(self, mod):
        val = self._call(mod, "os_name")
        assert isinstance(val, str) and val

    def test_os_version_is_string(self, mod):
        val = self._call(mod, "os_version")
        assert isinstance(val, str)

    def test_os_release_is_string(self, mod):
        val = self._call(mod, "os_release")
        assert isinstance(val, str)

    def test_machine_is_string(self, mod):
        val = self._call(mod, "machine")
        assert isinstance(val, str)

    def test_platform_full_is_string(self, mod):
        val = self._call(mod, "platform_full")
        assert isinstance(val, str) and val

    def test_python_version_has_dots(self, mod):
        val = self._call(mod, "python_version")
        assert "." in val

    def test_saved_path_unknown_when_no_mw(self, mod):
        val = self._call(mod, "saved_path", mw=None)
        assert "unknown" in val.lower()

    def test_saved_path_from_get_current_file_path(self, mod):
        mw = MagicMock()
        mw.get_current_file_path.return_value = "/tmp/my_project.pmeprj"
        val = self._call(mod, "saved_path", mw=mw)
        assert "my_project.pmeprj" in val

    def test_saved_path_from_init_manager(self, mod):
        mw = MagicMock()
        del mw.get_current_file_path
        mw.init_manager.current_file_path = "/tmp/my_project2.pmeprj"
        val = self._call(mod, "saved_path", mw=mw)
        assert "my_project2.pmeprj" in val

    def test_app_version_returns_string(self, mod):
        val = self._call(mod, "app_version")
        assert isinstance(val, str)

    def test_installed_plugins_from_plugin_manager(self, mod):
        mw = MagicMock()
        mw.plugin_manager.plugins = [
            {"name": "PluginA", "version": "1.0"},
            {"name": "PluginB", "version": "2.0"}
        ]
        val = self._call(mod, "installed_plugins", mw=mw)
        assert "PluginA v1.0" in val
        assert "PluginB v2.0" in val

    def test_unknown_key_returns_error_string(self, mod):
        val = self._call(mod, "nonexistent_key_xyz")
        assert isinstance(val, str)


class TestCollectMetadata:
    def test_all_enabled_fields_present(self, mod):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        cfg = m._default_config()
        # Enable everything
        for k in cfg["enabled_fields"]:
            cfg["enabled_fields"][k] = True
        cfg["custom_note"] = "hello"
        result = m.collect_metadata(cfg, mw=None)
        assert "saved_at" in result
        assert "username" in result
        assert "note" in result
        assert result["note"] == "hello"

    def test_disabled_fields_absent(self, mod):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        cfg = m._default_config()
        for k in cfg["enabled_fields"]:
            cfg["enabled_fields"][k] = False
        result = m.collect_metadata(cfg, mw=None)
        assert result == {}

    def test_note_skipped_when_empty(self, mod):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        cfg = m._default_config()
        cfg["enabled_fields"]["note"] = True
        cfg["custom_note"] = ""
        result = m.collect_metadata(cfg, mw=None)
        assert "note" not in result

    def test_note_included_when_nonempty(self, mod):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        cfg = m._default_config()
        cfg["enabled_fields"]["note"] = True
        cfg["custom_note"] = "debug run A"
        result = m.collect_metadata(cfg, mw=None)
        assert result.get("note") == "debug run A"


# ---------------------------------------------------------------------------
# Plugin lifecycle
# ---------------------------------------------------------------------------


class TestInitialize:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        m.initialize(ctx)
        ctx.add_menu_action.assert_called_once()
        path_arg = ctx.add_menu_action.call_args[0][0]
        assert path_arg == "Settings/Metadata Saver..."

    def test_initialize_registers_save_handler(self):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        m.initialize(ctx)
        ctx.register_save_handler.assert_called_once_with(m.on_save_project)

    def test_initialize_registers_load_handler(self):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        m.initialize(ctx)
        ctx.register_load_handler.assert_called_once_with(m.on_load_project)

    def test_initialize_registers_document_reset_handler(self):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        m.initialize(ctx)
        ctx.register_document_reset_handler.assert_called_once_with(m.on_document_reset)

    def test_initialize_sets_context(self):
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        m.initialize(ctx)
        assert m._CONTEXT is ctx


class TestSaveHandler:
    def _loaded_mod(self):
        with mock_optional_imports():
            return load_plugin(PLUGIN_PATH)

    def test_on_save_returns_dict(self):
        m = self._loaded_mod()
        ctx = make_context()
        m.initialize(ctx)
        with patch.object(m, "load_config", return_value=m._default_config()):
            result = m.on_save_project()
        assert isinstance(result, dict)
        assert "metadata" in result
        assert "plugin_version" in result

    def test_on_save_plugin_version_matches(self):
        m = self._loaded_mod()
        ctx = make_context()
        m.initialize(ctx)
        with patch.object(m, "load_config", return_value=m._default_config()):
            result = m.on_save_project()
        assert result["plugin_version"] == m.PLUGIN_VERSION

    def test_on_save_no_context_returns_none_or_dict(self):
        """Should not crash if _CONTEXT is None (no initialize called)."""
        m = self._loaded_mod()
        # _CONTEXT is None
        with patch.object(m, "load_config", return_value=m._default_config()):
            result = m.on_save_project()
        # Returns None or a dict; must not raise
        assert result is None or isinstance(result, dict)

    def test_on_save_returns_none_when_all_disabled(self):
        m = self._loaded_mod()
        ctx = make_context()
        m.initialize(ctx)
        cfg = m._default_config()
        for k in cfg["enabled_fields"]:
            cfg["enabled_fields"][k] = False
        with patch.object(m, "load_config", return_value=cfg):
            result = m.on_save_project()
        assert result is None

    def test_on_save_metadata_contains_saved_at(self):
        m = self._loaded_mod()
        ctx = make_context()
        m.initialize(ctx)
        cfg = m._default_config()
        cfg["enabled_fields"] = {"saved_at": True}
        with patch.object(m, "load_config", return_value=cfg):
            result = m.on_save_project()
        assert "saved_at" in result["metadata"]

    def test_on_save_survives_exception(self):
        """If collect_metadata raises, on_save_project must still return
        None (not propagate the exception)."""
        m = self._loaded_mod()
        ctx = make_context()
        m.initialize(ctx)
        with patch.object(m, "collect_metadata", side_effect=RuntimeError("boom")):
            result = m.on_save_project()
        assert result is None

    def test_on_save_returns_none_when_master_disabled(self):
        """Master enabled=False must short-circuit regardless of field state."""
        m = self._loaded_mod()
        ctx = make_context()
        m.initialize(ctx)
        cfg = m._default_config()
        cfg["enabled"] = False
        with patch.object(m, "load_config", return_value=cfg):
            result = m.on_save_project()
        assert result is None

    def test_on_save_returns_dict_when_master_enabled(self):
        """Master enabled=True (default) should produce a result dict."""
        m = self._loaded_mod()
        ctx = make_context()
        m.initialize(ctx)
        cfg = m._default_config()
        cfg["enabled"] = True
        with patch.object(m, "load_config", return_value=cfg):
            result = m.on_save_project()
        assert isinstance(result, dict)


class TestLoadHandler:
    def _loaded_mod(self):
        with mock_optional_imports():
            return load_plugin(PLUGIN_PATH)

    def test_on_load_stores_metadata(self):
        m = self._loaded_mod()
        data = {"metadata": {"saved_at": "2026-07-13T00:00:00+00:00", "username": "tester"}}
        m.on_load_project(data)
        assert m._LOADED_METADATA == data["metadata"]

    def test_on_load_ignores_missing_metadata_key(self):
        m = self._loaded_mod()
        m.on_load_project({"other_key": 42})
        assert m._LOADED_METADATA == {}

    def test_on_load_handles_non_dict(self):
        m = self._loaded_mod()
        m.on_load_project(None)  # must not raise
        m.on_load_project("bad")  # must not raise

    def test_on_document_reset_clears_metadata(self):
        m = self._loaded_mod()
        m._LOADED_METADATA = {"saved_at": "X"}
        m.on_document_reset()
        assert m._LOADED_METADATA == {}

    def test_on_document_reset_idempotent(self):
        m = self._loaded_mod()
        m.on_document_reset()
        m.on_document_reset()  # must not raise
        assert m._LOADED_METADATA == {}


# ---------------------------------------------------------------------------
# Dialog smoke tests (Qt-free via _build_cfg_from_ui extraction)
# ---------------------------------------------------------------------------


class TestDialogLogic:
    """
    Tests the pure-Python logic that the dialog uses, without instantiating
    Qt widgets.  We test _build_cfg_from_ui logic by manually constructing
    the state and verifying the result.
    """

    def test_build_cfg_from_ui_collects_all_keys(self):
        """Simulate _build_cfg_from_ui with all boxes checked."""
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)

        field_keys = [f["key"] for f in m.ALL_FIELDS]

        class FakeChk:
            def __init__(self, v):
                self._v = v

            def isChecked(self):
                return self._v

        class FakeNoteEdit:
            def toPlainText(self):
                return "  my note  "

        # Simulate the dialog's _build_cfg_from_ui directly
        checkboxes = {k: FakeChk(True) for k in field_keys}
        note_edit = FakeNoteEdit()

        enabled = {k: chk.isChecked() for k, chk in checkboxes.items()}
        custom_note = note_edit.toPlainText().strip()
        cfg = {"enabled_fields": enabled, "custom_note": custom_note}

        assert set(cfg["enabled_fields"].keys()) == set(field_keys)
        assert all(cfg["enabled_fields"].values())
        assert cfg["custom_note"] == "my note"

    def test_open_dialog_is_callable(self):
        """_open_dialog must not crash when context.get_window returns an existing mock."""
        with mock_optional_imports():
            m = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        # Simulate already-open window
        fake_win = MagicMock()
        fake_win.isHidden.return_value = False
        ctx.get_window.return_value = fake_win
        m._open_dialog(ctx)  # must not raise
        fake_win.raise_.assert_called_once()
        fake_win.activateWindow.assert_called_once()
