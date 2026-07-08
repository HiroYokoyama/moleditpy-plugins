"""
Tests for the Advanced Rendering plugin (get_icon; initialize -> 4 x register_3d_style + add_menu_action).
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
ADV_RENDER_PATH = PLUGINS_DIR / "Advanced_Rendering" / "advanced_rendering.py"


class TestAdvancedRendering:
    def test_get_icon_returns_none(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            assert mod.get_icon() is None

    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()

    def test_initialize_registers_four_3d_styles(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert ctx.register_3d_style.call_count == 4

    def test_initialize_menu_path_contains_advanced(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            path = ctx.add_menu_action.call_args[0][0]
            assert "Advanced" in path

    def test_initialize_style_names_contain_advanced_rendering(self):
        """All registered 3D style keys end with '(Advanced Rendering)'."""
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            style_keys = [
                call[0][0] for call in ctx.register_3d_style.call_args_list
            ]
            assert all("Advanced Rendering" in k for k in style_keys)
