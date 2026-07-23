"""
Headless GUI tests for the Vector Viewer plugin.

Covers: VectorViewerPlugin.
"""

from __future__ import annotations

import contextlib
import sys
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


# Real numpy/pyvista pulled in up front so they're cached before the
# chemistry-blocking finder snapshots sys.modules — otherwise the first-ever
# import would be intercepted and replaced with a MagicMock. Guarded so CI's
# bare test-gui job (only pytest+PyQt6 installed) skips this instead of
# erroring at collection.
_np = pytest.importorskip("numpy")
_pv = pytest.importorskip("pyvista")


@contextlib.contextmanager
def _mock_chemistry_keep_real_numpy_pv():
    """Like mock_chemistry_imports(), but numpy/pyvista resolve to the real
    packages so the plugin's real vector math and pv.Arrow calls actually run.
    """
    keep_prefixes = ("numpy", "pyvista")
    real_mods = {
        k: v for k, v in sys.modules.items() if k.split(".")[0] in keep_prefixes
    }
    with mock_chemistry_imports():
        sys.modules.update(real_mods)
        yield


with _mock_chemistry_keep_real_numpy_pv():
    _vector_real = load_plugin_for_gui(VECTOR_PATH)


def _real_mol(coords):
    """MagicMock molecule exposing a real-numpy-friendly conformer."""
    from types import SimpleNamespace

    mol = MagicMock()
    mol.GetNumAtoms.return_value = len(coords)
    conf = mol.GetConformer.return_value
    conf.GetAtomPosition.side_effect = [
        SimpleNamespace(x=c[0], y=c[1], z=c[2]) for c in coords
    ]
    return mol


def _real_ctx(mol=None):
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = mol
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


# ===========================================================================
# Real numpy/pyvista paths — get_com, full arrow draw, color picker,
# closeEvent error handling, export_image
# ===========================================================================


class TestVectorViewerRealMath:
    @pytest.fixture
    def rwin(self, qapp):
        ctx = _real_ctx()
        w = _vector_real.VectorViewerPlugin(ctx)
        yield w
        w.destroy()

    def test_get_com_no_molecule_returns_origin(self, rwin):
        rwin.context.current_molecule = None
        com = rwin.get_com()
        assert list(com) == [0.0, 0.0, 0.0]

    def test_get_com_averages_positions(self, rwin):
        rwin.context.current_molecule = _real_mol([(0, 0, 0), (2, 4, 6)])
        com = rwin.get_com()
        assert list(com) == [1.0, 2.0, 3.0]

    def test_get_com_conformer_error_returns_origin(self, rwin):
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 1
        mol.GetConformer.side_effect = RuntimeError("no conformer")
        rwin.context.current_molecule = mol
        com = rwin.get_com()
        assert list(com) == [0.0, 0.0, 0.0]

    def test_full_draw_creates_arrow_and_adds_mesh(self, rwin):
        plotter = MagicMock()
        rwin.context.plotter = plotter
        rwin.context.current_molecule = None
        rwin.vec_input.setText("1 0 0")
        rwin.update_visualization()
        plotter.add_mesh.assert_called_once()
        assert plotter.add_mesh.call_args.kwargs["name"] == "vector_viewer_arrow"
        plotter.render.assert_called_once()
        assert rwin.vis_actor is plotter.add_mesh.return_value

    def test_reverse_checkbox_flips_vector(self, rwin):
        plotter = MagicMock()
        rwin.context.plotter = plotter
        rwin.context.current_molecule = None
        rwin.reverse_chk.setChecked(True)
        rwin.vec_input.setText("1 0 0")
        rwin.update_visualization()
        plotter.add_mesh.assert_called_once()

    def test_redraw_removes_previous_actor(self, rwin):
        plotter = MagicMock()
        rwin.context.plotter = plotter
        rwin.context.current_molecule = None
        rwin.vec_input.setText("1 0 0")
        rwin.update_visualization()
        first_actor = rwin.vis_actor
        rwin.vec_input.setText("2 0 0")
        rwin.update_visualization()
        plotter.remove_actor.assert_called_with(first_actor)

    def test_scale_and_opacity_applied(self, rwin):
        plotter = MagicMock()
        rwin.context.plotter = plotter
        rwin.context.current_molecule = None
        rwin.scale_spin.setValue(2.0)
        rwin.opacity_spin.setValue(0.75)
        rwin.vec_input.setText("1 0 0")
        rwin.update_visualization()
        kwargs = plotter.add_mesh.call_args.kwargs
        assert kwargs["opacity"] == 0.75

    def test_choose_color_valid_updates_button(self, rwin, monkeypatch):
        from PyQt6.QtGui import QColor

        new_color = QColor("red")
        monkeypatch.setattr(
            _vector_real.QColorDialog, "getColor", staticmethod(lambda *a, **k: new_color)
        )
        rwin.choose_color()
        assert rwin.arrow_color.name() == new_color.name()
        assert rwin.color_btn.text() == new_color.name().upper()

    def test_choose_color_invalid_leaves_unchanged(self, rwin, monkeypatch):
        from PyQt6.QtGui import QColor

        invalid_color = QColor()  # default-constructed QColor is invalid
        monkeypatch.setattr(
            _vector_real.QColorDialog,
            "getColor",
            staticmethod(lambda *a, **k: invalid_color),
        )
        original = rwin.arrow_color
        rwin.choose_color()
        assert rwin.arrow_color is original

    def test_close_event_swallows_remove_actor_exception(self, rwin):
        plotter = MagicMock()
        plotter.remove_actor.side_effect = RuntimeError("boom")
        rwin.context.plotter = plotter
        rwin.vis_actor = object()
        rwin.close()  # must not raise
        assert rwin.vis_actor is None

    def test_export_image_no_plotter_warns(self, rwin, monkeypatch):
        rwin.context.plotter = None
        warned = {}
        monkeypatch.setattr(
            _vector_real.QMessageBox,
            "warning",
            staticmethod(lambda *a, **k: warned.setdefault("called", True)),
        )
        rwin.export_image()
        assert warned.get("called") is True

    def test_export_image_no_filename_is_noop(self, rwin, monkeypatch):
        plotter = MagicMock()
        rwin.context.plotter = plotter
        monkeypatch.setattr(
            _vector_real.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        rwin.export_image()
        plotter.screenshot.assert_not_called()

    def test_export_image_success_reports_status(self, rwin, monkeypatch):
        plotter = MagicMock()
        rwin.context.plotter = plotter
        monkeypatch.setattr(
            _vector_real.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("C:/out.png", "")),
        )
        rwin.trans_chk.setChecked(True)
        rwin.export_image()
        plotter.screenshot.assert_called_once_with(
            "C:/out.png", transparent_background=True
        )
        rwin.context.show_status_message.assert_called_once()

    def test_export_image_screenshot_error_shows_critical(self, rwin, monkeypatch):
        plotter = MagicMock()
        plotter.screenshot.side_effect = RuntimeError("disk full")
        rwin.context.plotter = plotter
        monkeypatch.setattr(
            _vector_real.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("C:/out.png", "")),
        )
        critical = {}
        monkeypatch.setattr(
            _vector_real.QMessageBox,
            "critical",
            staticmethod(lambda *a, **k: critical.setdefault("called", True)),
        )
        rwin.export_image()
        assert critical.get("called") is True
