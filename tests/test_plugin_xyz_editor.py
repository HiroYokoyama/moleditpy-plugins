"""
Tests for the XYZ Editor plugin (initialize -> add_menu_action + save/load/reset
handlers; save handler returns {} when no molecule).
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
XYZ_EDITOR_PATH = PLUGINS_DIR / "XYZ_Editor" / "xyz_editor.py"


class TestXYZEditor:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()

    def test_initialize_menu_path_contains_xyz(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            path = ctx.add_menu_action.call_args[0][0]
            assert "XYZ" in path

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_initialize_registers_save_handler(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_save_handler.assert_called_once()

    def test_initialize_registers_load_handler(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_load_handler.assert_called_once()

    def test_initialize_registers_document_reset_handler(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_document_reset_handler.assert_called_once()

    def test_save_handler_returns_empty_dict_when_no_molecule(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            ctx.current_molecule = None
            mod.initialize(ctx)
            save_fn = ctx.register_save_handler.call_args[0][0]
            result = save_fn()
            assert result == {} or result is None

    def test_load_handler_tolerates_empty_dict(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            ctx.current_molecule = None
            mod.initialize(ctx)
            load_fn = ctx.register_load_handler.call_args[0][0]
            load_fn({})
