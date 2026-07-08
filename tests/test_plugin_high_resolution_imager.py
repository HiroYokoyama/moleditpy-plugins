"""
Tests for the High Resolution Imager plugin (initialize -> add_export_action; PLUGIN_CONTEXT stored).
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
HI_RES_PATH = PLUGINS_DIR / "High_Resolution_Imager" / "high_res_imager.py"


class TestHighResolutionImager:
    def test_initialize_registers_export_action(self):
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_screenshot(self):
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            label = ctx.add_export_action.call_args[0][0]
            assert "Screenshot" in label

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_initialize_does_not_register_menu_action(self):
        """High Resolution Imager uses add_export_action, not add_menu_action."""
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_not_called()
