"""
Tests for the Vector Viewer plugin.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
VECTOR_VIEWER_PATH = PLUGINS_DIR / "Vector_Viewer" / "vector_viewer.py"


class TestVectorViewer:
    def test_initialize_sets_launch_fn(self):
        """initialize() stores launch in _launch_fn; run(mw) handles auto-registration."""
        with mock_optional_imports():
            mod = load_plugin(VECTOR_VIEWER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod._launch_fn is not None

    def test_initialize_shows_status_message(self):
        with mock_optional_imports():
            mod = load_plugin(VECTOR_VIEWER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.show_status_message.assert_called()

    def test_run_reuses_existing_window(self):
        """run() should show an existing window if one is registered."""
        with mock_optional_imports():
            mod = load_plugin(VECTOR_VIEWER_PATH)
            ctx = make_context()
            existing_win = MagicMock()
            ctx.get_window.return_value = existing_win
            mod.initialize(ctx)
            mod.run(ctx.get_main_window())




# ---------------------------------------------------------------------------
# get_com / update_visualization (extracted via AST)
# ---------------------------------------------------------------------------

import logging
from types import SimpleNamespace

import numpy as np

from conftest import extract_function

VECTOR_PATH = VECTOR_VIEWER_PATH


class TestVectorViewerGetCom:
    def _get_com(self):
        return extract_function(
            VECTOR_PATH, "VectorViewerPlugin", "get_com",
            extra_globals={"np": np, "logging": logging},
        )

    def _mol(self, coords):
        mol = MagicMock()
        mol.GetNumAtoms.return_value = len(coords)
        conf = mol.GetConformer.return_value
        conf.GetAtomPosition.side_effect = [
            SimpleNamespace(x=c[0], y=c[1], z=c[2]) for c in coords
        ]
        return mol

    def test_mean_of_positions(self):
        fn = self._get_com()
        s = MagicMock()
        s.context.current_molecule = self._mol([(0, 0, 0), (2, 4, 6)])
        com = fn(s)
        assert np.allclose(com, [1.0, 2.0, 3.0])

    def test_no_molecule_returns_origin(self):
        fn = self._get_com()
        s = MagicMock()
        s.context.current_molecule = None
        assert np.allclose(fn(s), [0.0, 0.0, 0.0])

    def test_conformer_error_returns_origin(self):
        fn = self._get_com()
        s = MagicMock()
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 1
        mol.GetConformer.side_effect = RuntimeError("no conformer")
        s.context.current_molecule = mol
        assert np.allclose(fn(s), [0.0, 0.0, 0.0])


class TestVectorViewerUpdateVisualization:
    def _fn(self, pv_mock):
        return extract_function(
            VECTOR_PATH, "VectorViewerPlugin", "update_visualization",
            extra_globals={"np": np, "pv": pv_mock},
        )

    def _self(self, text, scale=1.0, reverse=False, com=(0.0, 0.0, 0.0)):
        s = MagicMock()
        s.vec_input.text.return_value = text
        s.reverse_chk.isChecked.return_value = reverse
        s.scale_spin.value.return_value = scale
        s.res_spin.value.return_value = 20
        s.opacity_spin.value.return_value = 0.5
        s.arrow_color.name.return_value = "#00ff00"
        s.get_com.return_value = np.array(com)
        s.vis_actor = None
        return s

    def test_no_plotter_returns_early(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1 0 0")
        s.context.plotter = None
        fn(s)
        pv_mock.Arrow.assert_not_called()

    def test_empty_text_returns_early(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("   ")
        fn(s)
        pv_mock.Arrow.assert_not_called()

    def test_two_components_rejected(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        fn(self._self("1.0 2.0"))
        pv_mock.Arrow.assert_not_called()

    def test_non_numeric_rejected(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        fn(self._self("a b c"))
        pv_mock.Arrow.assert_not_called()

    def test_zero_vector_rejected(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        fn(self._self("0 0 0"))
        pv_mock.Arrow.assert_not_called()

    def test_comma_separated_vector_drawn(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1, 0, 0", scale=2.0)
        fn(s)
        kwargs = pv_mock.Arrow.call_args.kwargs
        assert np.allclose(kwargs["direction"], [2.0, 0.0, 0.0])
        # Arrow centred on COM: start = com - scaled/2
        assert np.allclose(kwargs["start"], [-1.0, 0.0, 0.0])
        s.context.plotter.add_mesh.assert_called_once()
        assert s.context.plotter.add_mesh.call_args.kwargs["name"] == \
            "vector_viewer_arrow"
        s.context.plotter.render.assert_called_once()

    def test_reverse_flips_direction(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1 2 3", reverse=True)
        fn(s)
        kwargs = pv_mock.Arrow.call_args.kwargs
        assert np.allclose(kwargs["direction"], [-1.0, -2.0, -3.0])

    def test_old_actor_removed_before_redraw(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1 0 0")
        old_actor = MagicMock()
        s.vis_actor = old_actor
        fn(s)
        s.context.plotter.remove_actor.assert_called_once_with(old_actor)

    def test_actor_stored_after_draw(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1 0 0")
        fn(s)
        assert s.vis_actor is s.context.plotter.add_mesh.return_value


class TestVectorViewerEntryPoints:
    def test_run_before_initialize_is_noop(self):
        with mock_optional_imports():
            mod = load_plugin(VECTOR_PATH)
        assert mod._launch_fn is None
        mod.run(MagicMock())  # must not raise

    def test_initialize_sets_launch_and_status(self):
        with mock_optional_imports():
            mod = load_plugin(VECTOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        assert callable(mod._launch_fn)
        ctx.show_status_message.assert_called_once()

    def test_launch_creates_new_window_when_none_registered(self):
        # Note: VectorViewerPlugin derives from the mocked QWidget, so a
        # direct instantiation doesn't run real __init__ logic (see the
        # extract_function-based tests below for real behaviour coverage).
        # This just guards that launch() doesn't crash on the "no existing
        # window" branch.
        with mock_optional_imports():
            mod = load_plugin(VECTOR_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            mod.initialize(ctx)
            mod._launch_fn()  # must not raise

    def test_run_dispatches_to_host_main_window(self):
        with mock_optional_imports():
            mod = load_plugin(VECTOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            called = {}
            mod._launch_fn = lambda: called.setdefault("ran", True)
            mw = SimpleNamespace(host=MagicMock())
            mod.run(mw)
        assert called.get("ran") is True


# ---------------------------------------------------------------------------
# export_image / closeEvent / choose_color / update_color_btn_style
# ---------------------------------------------------------------------------


_FAKE_SUPER_GLOBALS = {"super": lambda *a: SimpleNamespace(closeEvent=lambda e: None)}


class TestVectorViewerExportImage:
    def _fn(self, qfd=None, qmb=None):
        globs = {
            "QFileDialog": qfd or MagicMock(),
            "QMessageBox": qmb or MagicMock(),
        }
        return extract_function(
            VECTOR_PATH, "VectorViewerPlugin", "export_image", globs
        ), globs

    def test_no_plotter_shows_warning(self):
        fn, globs = self._fn()
        s = MagicMock()
        s.context.plotter = None
        fn(s)
        globs["QMessageBox"].warning.assert_called_once()
        globs["QFileDialog"].getSaveFileName.assert_not_called()

    def test_no_filename_chosen_is_noop(self):
        qfd = MagicMock()
        qfd.getSaveFileName.return_value = ("", "")
        fn, globs = self._fn(qfd=qfd)
        s = MagicMock()
        fn(s)
        s.context.plotter.screenshot.assert_not_called()

    def test_success_takes_screenshot_and_reports_status(self):
        qfd = MagicMock()
        qfd.getSaveFileName.return_value = ("C:/out.png", "")
        fn, globs = self._fn(qfd=qfd)
        s = MagicMock()
        s.trans_chk.isChecked.return_value = True
        fn(s)
        s.context.plotter.screenshot.assert_called_once_with(
            "C:/out.png", transparent_background=True
        )
        s.context.show_status_message.assert_called_once()

    def test_screenshot_exception_shows_critical(self):
        qfd = MagicMock()
        qfd.getSaveFileName.return_value = ("C:/out.png", "")
        fn, globs = self._fn(qfd=qfd)
        s = MagicMock()
        s.context.plotter.screenshot.side_effect = RuntimeError("disk full")
        fn(s)
        globs["QMessageBox"].critical.assert_called_once()
        s.context.show_status_message.assert_not_called()


class TestVectorViewerCloseEvent:
    def _fn(self):
        globs = dict(_FAKE_SUPER_GLOBALS)
        return extract_function(VECTOR_PATH, "VectorViewerPlugin", "closeEvent", globs)

    def test_removes_actor_and_resets(self):
        fn = self._fn()
        s = MagicMock()
        actor = MagicMock()
        s.vis_actor = actor
        fn(s, MagicMock())
        s.context.plotter.remove_actor.assert_called_once_with(actor)
        s.context.plotter.render.assert_called_once()
        assert s.vis_actor is None

    def test_no_plotter_does_not_crash(self):
        fn = self._fn()
        s = MagicMock()
        s.context.plotter = None
        fn(s, MagicMock())  # must not raise
        assert s.vis_actor is None

    def test_remove_actor_exception_is_swallowed(self):
        fn = self._fn()
        s = MagicMock()
        s.context.plotter.remove_actor.side_effect = RuntimeError("boom")
        fn(s, MagicMock())  # must not raise
        assert s.vis_actor is None


class TestVectorViewerChooseColor:
    def _fn(self, qcd):
        globs = {"QColorDialog": qcd}
        return extract_function(VECTOR_PATH, "VectorViewerPlugin", "choose_color", globs)

    def test_valid_color_updates_state(self):
        qcd = MagicMock()
        new_color = MagicMock()
        new_color.isValid.return_value = True
        qcd.getColor.return_value = new_color
        fn = self._fn(qcd)
        s = MagicMock()
        fn(s)
        assert s.arrow_color is new_color
        s.update_color_btn_style.assert_called_once()
        s.update_visualization.assert_called_once()

    def test_invalid_color_leaves_state_unchanged(self):
        qcd = MagicMock()
        new_color = MagicMock()
        new_color.isValid.return_value = False
        qcd.getColor.return_value = new_color
        fn = self._fn(qcd)
        s = MagicMock()
        original = s.arrow_color
        fn(s)
        assert s.arrow_color is original
        s.update_color_btn_style.assert_not_called()
        s.update_visualization.assert_not_called()


class TestVectorViewerUpdateColorBtnStyle:
    def _fn(self):
        return extract_function(
            VECTOR_PATH, "VectorViewerPlugin", "update_color_btn_style", {}
        )

    def test_light_color_uses_black_text(self):
        fn = self._fn()
        s = MagicMock()
        s.arrow_color.name.return_value = "#ffffff"
        s.arrow_color.lightness.return_value = 255
        fn(s)
        style = s.color_btn.setStyleSheet.call_args[0][0]
        assert "color: black" in style
        s.color_btn.setText.assert_called_once_with("#FFFFFF")

    def test_dark_color_uses_white_text(self):
        fn = self._fn()
        s = MagicMock()
        s.arrow_color.name.return_value = "#000000"
        s.arrow_color.lightness.return_value = 0
        fn(s)
        style = s.color_btn.setStyleSheet.call_args[0][0]
        assert "color: white" in style
