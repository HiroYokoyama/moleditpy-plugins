"""
Headless GUI tests for the Vector Viewer plugin.

Covers: VectorViewerPlugin.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

VECTOR_PATH = PLUGINS_DIR / "Vector_Viewer" / "vector_viewer.py"

with mock_chemistry_imports():
    _vector = load_plugin_for_gui(VECTOR_PATH)


def _vector_context():
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


# ===========================================================================
# VectorViewerPlugin  (Vector Viewer)
# ===========================================================================


class TestVectorViewerPlugin:
    @pytest.fixture
    def win(self, qapp):
        ctx = _vector_context()
        w = _vector.VectorViewerPlugin(ctx)
        yield w
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "Vector Viewer"

    def test_vector_input_placeholder(self, win):
        assert "vector" in win.vec_input.placeholderText().lower()

    def test_default_scale(self, win):
        assert win.scale_spin.value() == 1.0

    def test_default_resolution(self, win):
        assert win.res_spin.value() == 20

    def test_default_opacity(self, win):
        assert win.opacity_spin.value() == 0.5

    def test_reverse_unchecked_by_default(self, win):
        assert not win.reverse_chk.isChecked()

    def test_transparent_export_checked_by_default(self, win):
        assert win.trans_chk.isChecked()

    def test_color_button_shows_green_hex(self, win):
        assert win.color_btn.text() == "#008000"

    def test_registered_as_main_panel(self, win):
        win.context.register_window.assert_called_with("main_panel", win)

    def test_update_with_no_plotter_is_noop(self, win):
        win.context.plotter = None
        win.vec_input.setText("1 0 0")
        win.update_visualization()  # must not raise

    def test_update_with_empty_input_skips_plotter(self, win):
        plotter = MagicMock()
        win.context.plotter = plotter
        win.vec_input.setText("")
        win.update_visualization()
        plotter.add_mesh.assert_not_called()

    def test_update_with_invalid_input_skips_plotter(self, win):
        plotter = MagicMock()
        win.context.plotter = plotter
        win.vec_input.setText("1 2")  # fewer than 3 components
        win.update_visualization()
        plotter.add_mesh.assert_not_called()

    def test_close_removes_actor(self, win):
        plotter = MagicMock()
        win.context.plotter = plotter
        sentinel = object()
        win.vis_actor = sentinel
        win.close()
        plotter.remove_actor.assert_called_once_with(sentinel)
        assert win.vis_actor is None

    def test_initialize_sets_launch_fn(self, qapp):
        ctx = _vector_context()
        _vector._launch_fn = None
        _vector.initialize(ctx)
        assert callable(_vector._launch_fn)
        ctx.show_status_message.assert_called_once()
        _vector._launch_fn = None
