"""
Headless GUI tests for the Metadata Saver plugin.

Covers: MetadataSaverDialog + save/load handlers.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

META_PATH = PLUGINS_DIR / "Metadata_Saver" / "metadata_saver.py"

with mock_chemistry_imports():
    _meta = load_plugin_for_gui(META_PATH)


# ===========================================================================
# Metadata Saver — handlers and MetadataSaverDialog
# ===========================================================================


class TestMetadataHandlers:
    def test_default_config_shape(self):
        cfg = _meta._default_config()
        assert cfg["enabled"] is True
        assert cfg["silence_notice"] is False
        assert set(cfg["enabled_fields"]) == {f["key"] for f in _meta.ALL_FIELDS}

    def test_collect_metadata_skips_disabled_fields(self):
        cfg = _meta._default_config()
        cfg["enabled_fields"] = {k: False for k in cfg["enabled_fields"]}
        cfg["enabled_fields"]["os_name"] = True
        meta = _meta.collect_metadata(cfg, None)
        assert list(meta) == ["os_name"]

    def test_collect_metadata_skips_empty_note(self):
        cfg = _meta._default_config()
        cfg["enabled_fields"]["note"] = True
        cfg["custom_note"] = "   "
        assert "note" not in _meta.collect_metadata(cfg, None)

    def test_on_save_project_returns_versioned_payload(self):
        old_ctx = _meta._CONTEXT
        _meta._CONTEXT = MagicMock()
        try:
            data = _meta.on_save_project()
        finally:
            _meta._CONTEXT = old_ctx
        assert data["plugin_version"] == _meta.PLUGIN_VERSION
        assert "saved_at" in data["metadata"]

    def test_on_load_project_roundtrip_and_reset(self):
        _meta.on_load_project({"metadata": {"username": "alice"}})
        assert _meta._LOADED_METADATA == {"username": "alice"}
        _meta.on_document_reset()
        assert _meta._LOADED_METADATA == {}

    def test_on_load_project_ignores_malformed_data(self):
        _meta.on_load_project({"metadata": {"username": "alice"}})
        _meta.on_load_project("not a dict")
        assert _meta._LOADED_METADATA == {}

    def test_initialize_registers_menu_and_handlers(self, qapp):
        old_ctx = _meta._CONTEXT
        ctx = MagicMock()
        try:
            _meta.initialize(ctx)
        finally:
            _meta._CONTEXT = old_ctx
        assert (
            ctx.add_menu_action.call_args.args[0] == "Settings/Metadata Saver..."
        )
        ctx.register_save_handler.assert_called_once()
        ctx.register_load_handler.assert_called_once()
        ctx.register_document_reset_handler.assert_called_once()


class TestMetadataSaverDialog:
    @pytest.fixture
    def dlg(self, qapp):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        d = _meta.MetadataSaverDialog(ctx, parent=None)
        yield d
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Metadata Saver — Settings"

    def test_master_toggle_checked_by_default(self, dlg):
        assert dlg._chk_master.isChecked()

    def test_one_checkbox_per_field(self, dlg):
        assert set(dlg._checkboxes) == {f["key"] for f in _meta.ALL_FIELDS}

    def test_checkbox_defaults_match_registry(self, dlg):
        for field in _meta.ALL_FIELDS:
            assert dlg._checkboxes[field["key"]].isChecked() == field["default"]

    def test_master_toggle_disables_fields(self, dlg):
        dlg._chk_master.setChecked(False)
        assert not dlg._fields_container.isEnabled()
        dlg._chk_master.setChecked(True)
        assert dlg._fields_container.isEnabled()

    def test_build_cfg_reflects_ui_state(self, dlg):
        dlg._checkboxes["username"].setChecked(False)
        dlg._note_edit.setText("hello")
        cfg = dlg._build_cfg_from_ui()
        assert cfg["enabled_fields"]["username"] is False
        assert cfg["custom_note"] == "hello"

    def test_preview_fills_viewer(self, dlg):
        dlg._on_preview()
        text = dlg._metadata_view.toPlainText()
        assert "os_name:" in text

    def test_preview_with_nothing_selected(self, dlg):
        for chk in dlg._checkboxes.values():
            chk.setChecked(False)
        dlg._on_preview()
        assert dlg._metadata_view.toPlainText() == "(no fields selected)"
