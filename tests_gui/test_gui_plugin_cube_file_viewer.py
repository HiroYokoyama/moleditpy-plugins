"""
Headless GUI tests for the Cube File Viewer plugin.

Covers: ChargeDialog (the "bond connectivity error" dialog).

Chemistry libs (pyvista, numpy, rdkit, …) are mocked; real PyQt6 is used.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CUBE_VIEWER_PATH = PLUGINS_DIR / "Cube_File_Viewer" / "cube_viewer.py"

with mock_chemistry_imports():
    _cube = load_plugin_for_gui(CUBE_VIEWER_PATH)


# ===========================================================================
# ChargeDialog  (Cube File Viewer)
# ===========================================================================


class TestCubeViewerChargeDialog:
    """ChargeDialog: the 'bond connectivity error' dialog in Cube File Viewer."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _cube.ChargeDialog(parent=None, current_charge=0)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Bond Connectivity Error"

    def test_default_result_action_is_cancel(self, dlg):
        assert dlg.result_action == "cancel"

    def test_charge_spinbox_range(self, dlg):
        assert dlg.spin.minimum() == -20
        assert dlg.spin.maximum() == 20

    def test_default_charge_value(self, dlg):
        assert dlg.spin.value() == 0

    def test_nonzero_charge_reflected_in_spin(self, qapp):
        d = _cube.ChargeDialog(parent=None, current_charge=2)
        assert d.spin.value() == 2
        d.destroy()

    def test_on_retry_sets_result_action(self, dlg):
        dlg.spin.setValue(1)
        dlg.on_retry()
        assert dlg.result_action == "retry"
        assert dlg.charge == 1

    def test_on_skip_sets_result_action(self, dlg):
        dlg.on_skip()
        assert dlg.result_action == "skip"


# ===========================================================================
# CubeViewerWidget — construction, controls, persistence, cleanup
# ===========================================================================

import os
import tempfile


def _make_main_window():
    from PyQt6.QtWidgets import QWidget

    mw = QWidget()
    mw.plotter = MagicMock()
    mw.removeDockWidget = MagicMock()
    mw.current_mol = None
    return mw


def _make_viewer(qapp, data_max=2.0):
    mw = _make_main_window()
    dock = MagicMock()
    grid = MagicMock()
    w = _cube.CubeViewerWidget(mw, dock, grid, data_max=data_max)
    # Heavy pyvista rendering; the control tests only care about wiring.
    w.update_iso = MagicMock()
    w.get_settings_path = lambda: os.path.join(
        tempfile.gettempdir(), "cube_viewer_test.json"
    )
    return w, mw, dock


def _teardown(w):
    from PyQt6.QtCore import QCoreApplication

    w._structure_watch_timer.stop()
    try:
        QCoreApplication.instance().aboutToQuit.disconnect(w.save_settings)
    except TypeError:
        pass
    w.destroy()


class TestCubeViewerWidgetConstruction:
    def test_data_max_floor(self, qapp):
        w, mw, dock = _make_viewer(qapp, data_max=0.0)
        assert w.data_max == pytest.approx(1e-6)
        _teardown(w)

    def test_isovalue_spin_range(self, qapp):
        w, mw, dock = _make_viewer(qapp, data_max=2.0)
        assert w.spin.minimum() == pytest.approx(0.00001)
        assert w.spin.maximum() == pytest.approx(2.0 * 1.2)
        _teardown(w)

    def test_default_isovalue(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        assert w.spin.value() == pytest.approx(0.05)
        _teardown(w)

    def test_slider_range(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        assert w.slider.minimum() == 0
        assert w.slider.maximum() == 1000
        _teardown(w)

    def test_default_colors(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        assert w.color_p == "#0000FF"
        assert w.color_n == "#FF0000"
        _teardown(w)

    def test_structure_watch_timer_running(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        assert w._structure_watch_timer.isActive()
        _teardown(w)


class TestCubeViewerControls:
    def test_spin_change_syncs_slider(self, qapp):
        w, mw, dock = _make_viewer(qapp, data_max=2.0)
        w.spin.setValue(1.2)  # 1.2 / (2*1.2) = 0.5 → 500
        assert w.slider.value() == 500
        w.update_iso.assert_called()
        _teardown(w)

    def test_slider_change_syncs_spin(self, qapp):
        w, mw, dock = _make_viewer(qapp, data_max=2.0)
        w.slider.setValue(250)  # 250/1000 * 2.4 = 0.6
        assert w.spin.value() == pytest.approx(0.6)
        _teardown(w)

    def test_opacity_slider_syncs_spin(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        w.opacity_slider.setValue(80)
        assert w.opacity_spin.value() == pytest.approx(0.8)
        _teardown(w)

    def test_opacity_spin_syncs_slider(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        w.opacity_spin.setValue(0.25)
        assert w.opacity_slider.value() == 25
        _teardown(w)

    def test_complementary_toggle_disables_neg_button(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        w.check_comp_color.setChecked(True)
        assert not w.btn_color_n.isEnabled()
        w.check_comp_color.setChecked(False)
        assert w.btn_color_n.isEnabled()
        _teardown(w)


class TestCubeViewerPersistence:
    def test_settings_round_trip(self, qapp, tmp_path):
        w, mw, dock = _make_viewer(qapp)
        path = str(tmp_path / "cube_viewer.json")
        w.get_settings_path = lambda: path
        w.spin.setValue(0.33)
        w.color_p = "#00ff00"
        w.opacity_spin.setValue(0.7)
        w.save_settings()

        w2, mw2, dock2 = _make_viewer(qapp)
        w2.get_settings_path = lambda: path
        w2.load_settings()
        assert w2.spin.value() == pytest.approx(0.33)
        assert w2.color_p == "#00ff00"
        assert w2.opacity_spin.value() == pytest.approx(0.7)
        _teardown(w)
        _teardown(w2)


class TestCubeViewerCleanup:
    def test_close_event_stops_timer(self, qapp):
        from PyQt6.QtGui import QCloseEvent

        w, mw, dock = _make_viewer(qapp)
        w.save_settings = MagicMock()
        w.closeEvent(QCloseEvent())
        assert not w._structure_watch_timer.isActive()
        w.save_settings.assert_called_once()
        w.destroy()

    def test_structure_change_detaches_stale_view(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        w.save_settings = MagicMock()
        sentinel = object()
        w.bind_structure(sentinel)
        mw.current_mol = object()  # a different molecule now shown
        w._check_structure_still_bound()
        assert not w._structure_watch_timer.isActive()
        mw.removeDockWidget.assert_called_once()
        w.destroy()

    def test_no_detach_while_structure_unchanged(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        sentinel = object()
        w.bind_structure(sentinel)
        mw.current_mol = sentinel
        w._check_structure_still_bound()
        assert w._structure_watch_timer.isActive()
        mw.removeDockWidget.assert_not_called()
        _teardown(w)
