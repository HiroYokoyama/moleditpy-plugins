"""
Tests for the Advanced Rendering plugin (get_icon; initialize -> 4 x register_3d_style + add_menu_action).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

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


class TestAdvancedRenderingStyleDrawers:
    def _drawers(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        return {
            call.args[0]: call.args[1]
            for call in ctx.register_3d_style.call_args_list
        }

    def test_four_styles_registered(self):
        drawers = self._drawers()
        assert set(drawers) == {
            "Ball & Stick (Advanced Rendering)",
            "CPK (Advanced Rendering)",
            "Wireframe (Advanced Rendering)",
            "Stick (Advanced Rendering)",
        }

    def test_drawer_uses_normalized_style_key(self):
        drawers = self._drawers()
        viewer = MagicMock()
        mw_obj = SimpleNamespace(
            view_3d_manager=MagicMock(), _adv_rendering_viewer=viewer
        )
        mw_obj.view_3d_manager.current_3d_style = "Ball & Stick (Advanced Rendering)"

        drawers["Ball & Stick (Advanced Rendering)"](mw_obj, "MOL")

        mw_obj.view_3d_manager.draw_standard_3d_style.assert_called_once_with(
            "MOL", style_override="ball_and_stick"
        )
        viewer.apply_pbr_forced.assert_called_once()
        viewer.update_lights.assert_called_once()
        viewer.sync_style_ui.assert_called_once_with(
            "Ball & Stick (Advanced Rendering)"
        )

    def test_cpk_drawer_key(self):
        drawers = self._drawers()
        mw_obj = SimpleNamespace(
            view_3d_manager=MagicMock(), _adv_rendering_viewer=MagicMock()
        )
        drawers["CPK (Advanced Rendering)"](mw_obj, "MOL")
        mw_obj.view_3d_manager.draw_standard_3d_style.assert_called_once_with(
            "MOL", style_override="cpk"
        )

    def test_drawer_without_viewer_does_not_crash(self):
        drawers = self._drawers()
        mw_obj = SimpleNamespace(view_3d_manager=MagicMock())  # no viewer attr
        drawers["Stick (Advanced Rendering)"](mw_obj, "MOL")
        mw_obj.view_3d_manager.draw_standard_3d_style.assert_called_once()

    def test_drawer_without_view_manager_does_not_crash(self):
        # Regression (2026.07.08): sync_style_ui read mw_obj.view_3d_manager
        # unguarded even though the draw call above is hasattr-guarded for it.
        drawers = self._drawers()
        viewer = MagicMock()
        mw_obj = SimpleNamespace(_adv_rendering_viewer=viewer)
        drawers["Wireframe (Advanced Rendering)"](mw_obj, "MOL")
        viewer.sync_style_ui.assert_called_once_with("")
