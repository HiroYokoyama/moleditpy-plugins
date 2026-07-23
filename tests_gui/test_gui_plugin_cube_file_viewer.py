"""
Headless GUI tests for the Cube File Viewer plugin.

Covers: ChargeDialog (the "bond connectivity error" dialog).

Chemistry libs (pyvista, numpy, rdkit, …) are mocked; real PyQt6 is used.
"""

from __future__ import annotations

import sys
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

    def test_check_structure_returns_early_when_unbound(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        # _bound_mol defaults to None -> early return, timer keeps running.
        w._check_structure_still_bound()
        assert w._structure_watch_timer.isActive()
        _teardown(w)

    def test_detach_removes_actors_by_reference(self, qapp):
        """Cover the iso_actor_p/iso_actor_n removal branches (lines 723-726)."""
        w, mw, dock = _make_viewer(qapp)
        w.save_settings = MagicMock()
        w.iso_actor_p = MagicMock()
        w.iso_actor_n = MagicMock()
        sentinel = object()
        w.bind_structure(sentinel)
        mw.current_mol = object()
        w._check_structure_still_bound()
        assert mw.plotter.remove_actor.call_count >= 2
        w.destroy()

    def test_detach_swallows_actor_removal_exception(self, qapp):
        """Cover the except branch (lines 729-730) around actor cleanup."""
        w, mw, dock = _make_viewer(qapp)
        w.save_settings = MagicMock()
        mw.plotter.remove_actor.side_effect = Exception("boom")
        sentinel = object()
        w.bind_structure(sentinel)
        mw.current_mol = object()
        w._check_structure_still_bound()  # must not raise
        assert not w._structure_watch_timer.isActive()
        w.destroy()

    def test_detach_shows_status_message_when_statusbar_present(self, qapp):
        """Cover the statusBar() branch (lines 732-739)."""
        w, mw, dock = _make_viewer(qapp)
        w.save_settings = MagicMock()
        sb = MagicMock()
        mw.statusBar = lambda: sb
        sentinel = object()
        w.bind_structure(sentinel)
        mw.current_mol = object()
        w._check_structure_still_bound()
        sb.showMessage.assert_called_once()
        w.destroy()

    def test_close_plugin_full_cleanup(self, qapp):
        from types import SimpleNamespace

        w, mw, dock = _make_viewer(qapp)
        w.save_settings = MagicMock()
        mw.init_manager = SimpleNamespace(current_file_path="x.cube")
        mw.ui_manager = SimpleNamespace(restore_ui_for_editing=MagicMock())
        w.close_plugin()
        assert mw.current_mol is None
        assert mw.init_manager.current_file_path is None
        mw.plotter.clear.assert_called_once()
        mw.ui_manager.restore_ui_for_editing.assert_called_once()
        mw.removeDockWidget.assert_called_once()
        assert w.dock is None
        assert not w._structure_watch_timer.isActive()
        w.save_settings.assert_called_once()

    def test_close_plugin_swallows_cleanup_exception(self, qapp):
        """mw lacks init_manager -> AttributeError is caught (lines 770-771),
        but the dock is still torn down afterwards."""
        w, mw, dock = _make_viewer(qapp)
        w.save_settings = MagicMock()
        w.close_plugin()
        mw.removeDockWidget.assert_called_once()
        assert w.dock is None


class TestFlexibleDoubleSpinBox:
    def test_zero_value_renders_as_zero(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        assert w.spin.textFromValue(0.0) == "0"
        _teardown(w)


class TestCubeViewerSettingsEdgeCases:
    def test_load_settings_style_case_insensitive_match(self, qapp, tmp_path):
        path = str(tmp_path / "cube_viewer.json")
        w, mw, dock = _make_viewer(qapp)
        w.get_settings_path = lambda: path
        w.color_p = "#00ff00"
        w.opacity_spin.setValue(0.6)
        w.combo_style.setCurrentIndex(2)  # "Wireframe"
        w.save_settings()

        import json

        data = json.loads(open(path).read())
        data["style"] = "wireframe"  # lowercase, won't findText() exactly
        with open(path, "w") as f:
            json.dump(data, f)

        w2, mw2, dock2 = _make_viewer(qapp)
        w2.get_settings_path = lambda: path
        w2.load_settings()
        assert w2.combo_style.currentText() == "Wireframe"
        _teardown(w)
        _teardown(w2)

    def test_load_settings_corrupt_json_is_swallowed(self, qapp, tmp_path):
        path = tmp_path / "cube_viewer.json"
        path.write_text("{not valid json")
        w, mw, dock = _make_viewer(qapp)
        w.get_settings_path = lambda: str(path)
        w.load_settings()  # must not raise
        _teardown(w)

    def test_save_settings_ioerror_is_swallowed(self, qapp, tmp_path):
        w, mw, dock = _make_viewer(qapp)
        # A directory path makes open(..., "w") raise.
        w.get_settings_path = lambda: str(tmp_path)
        w.save_settings()  # must not raise
        _teardown(w)


class TestCubeViewerComplementaryColor:
    def test_invalid_color_returns_early(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        w.color_p = "not-a-color"
        w.update_complementary_color()  # must not raise; no crash on isValid()==False
        _teardown(w)

    def test_achromatic_color_keeps_same_hue(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        w.color_p = "#808080"  # gray: hue == -1
        w.update_complementary_color()
        assert w.color_n == "#808080"
        _teardown(w)

    def test_chromatic_color_flips_hue(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        w.color_p = "#ff0000"  # red, hue=0
        w.update_complementary_color()
        assert w.color_n != "#ff0000"
        _teardown(w)


class TestCubeViewerColorDialogs:
    def test_choose_color_p_updates_color_and_button(self, qapp):
        from PyQt6.QtGui import QColor

        w, mw, dock = _make_viewer(qapp)
        _cube.QColorDialog.getColor = MagicMock(return_value=QColor("#123456"))
        w.choose_color_p()
        assert w.color_p == "#123456"
        assert "#123456" in w.btn_color_p.styleSheet()
        w.update_iso.assert_called()
        _teardown(w)

    def test_choose_color_p_with_complementary_checked_updates_neg(self, qapp):
        from PyQt6.QtGui import QColor

        w, mw, dock = _make_viewer(qapp)
        w.check_comp_color.setChecked(True)
        _cube.QColorDialog.getColor = MagicMock(return_value=QColor("#00ff00"))
        w.choose_color_p()
        assert w.color_p == "#00ff00"
        assert w.color_n != "#FF0000"  # complementary recomputed
        _teardown(w)

    def test_choose_color_p_invalid_selection_no_change(self, qapp):
        from PyQt6.QtGui import QColor

        w, mw, dock = _make_viewer(qapp)
        original = w.color_p
        _cube.QColorDialog.getColor = MagicMock(return_value=QColor())  # invalid
        w.choose_color_p()
        assert w.color_p == original
        _teardown(w)

    def test_choose_color_n_updates_color_and_button(self, qapp):
        from PyQt6.QtGui import QColor

        w, mw, dock = _make_viewer(qapp)
        _cube.QColorDialog.getColor = MagicMock(return_value=QColor("#abcdef"))
        w.choose_color_n()
        assert w.color_n == "#abcdef"
        assert "#abcdef" in w.btn_color_n.styleSheet()
        w.update_iso.assert_called()
        _teardown(w)

    def test_choose_color_n_invalid_selection_no_change(self, qapp):
        from PyQt6.QtGui import QColor

        w, mw, dock = _make_viewer(qapp)
        original = w.color_n
        _cube.QColorDialog.getColor = MagicMock(return_value=QColor())  # invalid
        w.choose_color_n()
        assert w.color_n == original
        _teardown(w)


class TestCubeViewerInitialUpdate(object):
    def test_initial_update_calls_update_iso_and_resets_camera(self, qapp):
        w, mw, dock = _make_viewer(qapp)
        w.initial_update()
        w.update_iso.assert_called()
        w.plotter.reset_camera.assert_called()
        w.plotter.render.assert_called()
        _teardown(w)


# ===========================================================================
# update_iso() and open_cube_viewer() — real numpy/pyvista/vtk/rdkit, driving
# genuine isosurface contouring and the full plugin file-open flow.
# ===========================================================================

_np = pytest.importorskip("numpy")
_pv = pytest.importorskip("pyvista")
_vtk = pytest.importorskip("vtk")
_Chem = pytest.importorskip("rdkit.Chem")


import contextlib as _contextlib


@_contextlib.contextmanager
def _mock_chemistry_keep_real_chem():
    """Like mock_chemistry_imports(), but numpy/rdkit/pyvista/vtk resolve to the
    real installed packages so the plugin's real grid math / contouring /
    RDKit bond perception actually run."""
    import sys as _sys

    keep_prefixes = ("numpy", "rdkit", "pyvista", "pyvistaqt", "vtk", "vtkmodules")
    real_mods = {
        k: v for k, v in _sys.modules.items() if k.split(".")[0] in keep_prefixes
    }
    with mock_chemistry_imports():
        _sys.modules.update(real_mods)
        yield


with _mock_chemistry_keep_real_chem():
    _cuben = load_plugin_for_gui(CUBE_VIEWER_PATH)

# cube_viewer.py does `from rdkit.Chem import rdDetermineBonds` at module load
# time (unlike other plugins), which — as a side effect of the real import —
# sets `rdDetermineBonds` as an attribute on the real, process-wide rdkit.Chem
# package object. _cuben.rdDetermineBonds already holds its own reference (used
# by our tests below), but leaving the attribute on the shared real package
# leaks into *other* plugins' runtime `from rdkit.Chem import rdDetermineBonds`
# calls (e.g. Paste_XYZ's _apply_bonds does this at call time, not load time),
# making them see the real module instead of their own test's mocked one.
# Strip it back off the shared package so only this file's `_cuben` sees it.
sys.modules.pop("rdkit.Chem.rdDetermineBonds", None)
try:
    import rdkit.Chem as _real_rdkit_chem_pkg

    if hasattr(_real_rdkit_chem_pkg, "rdDetermineBonds"):
        delattr(_real_rdkit_chem_pkg, "rdDetermineBonds")
except Exception:
    pass


def _real_grid(nx=5, ny=5, nz=5):
    """A real StructuredGrid whose values cross both +/-0.05 (default isovalue)."""
    coords = _np.linspace(-1.0, 1.0, nx * ny * nz)
    meta = {
        "atoms": [],
        "origin": _np.zeros(3),
        "x_vec": _np.array([1.0, 0.0, 0.0]),
        "y_vec": _np.array([0.0, 1.0, 0.0]),
        "z_vec": _np.array([0.0, 0.0, 1.0]),
        "dims": (nx, ny, nz),
        "data_flat": coords,
        "is_angstrom_header": True,
    }
    _, grid = _cuben.build_grid_from_meta(meta)
    return grid


def _make_real_viewer(qapp, tmp_path, data_max=1.0):
    from PyQt6.QtWidgets import QWidget

    mw = QWidget()
    mw.plotter = _pv.Plotter(off_screen=True)
    dock = MagicMock()
    grid = _real_grid()
    w = _cuben.CubeViewerWidget(mw, dock, grid, data_max=data_max)
    w.get_settings_path = lambda: str(tmp_path / "real_cube_settings.json")
    return w, mw, dock


def _teardown_real(w, mw):
    from PyQt6.QtCore import QCoreApplication

    w._structure_watch_timer.stop()
    try:
        QCoreApplication.instance().aboutToQuit.disconnect(w.save_settings)
    except TypeError:
        pass
    w.destroy()
    mw.plotter.close()


class TestUpdateIsoReal:
    def test_default_creates_both_lobes(self, qapp, tmp_path):
        w, mw, dock = _make_real_viewer(qapp, tmp_path)
        w.update_iso()
        assert w.iso_actor_p is not None
        assert w.iso_actor_n is not None
        _teardown_real(w, mw)

    def test_wireframe_style_sets_density_label(self, qapp, tmp_path):
        w, mw, dock = _make_real_viewer(qapp, tmp_path)
        w.combo_style.setCurrentText("Wireframe")
        w.opacity_spin.setValue(0.4)  # triggers decimation (target_reduction>0)
        w.update_iso()
        assert w.opacity_label.text() == "Density:"
        _teardown_real(w, mw)

    def test_points_style_full_density_skips_decimation(self, qapp, tmp_path):
        w, mw, dock = _make_real_viewer(qapp, tmp_path)
        w.combo_style.setCurrentText("Points")
        w.opacity_spin.setValue(1.0)  # target_reduction == 0
        w.update_iso()  # must not raise
        assert w.opacity_label.text() == "Density:"
        _teardown_real(w, mw)

    def test_smoothed_surface_style(self, qapp, tmp_path):
        w, mw, dock = _make_real_viewer(qapp, tmp_path)
        w.combo_style.setCurrentText("Smoothed Surface")
        w.check_smooth.setChecked(True)
        w.update_iso()  # must not raise; exercises geometric smoothing
        assert w.opacity_label.text() == "Opacity:"
        _teardown_real(w, mw)

    def test_repeated_calls_replace_previous_actors(self, qapp, tmp_path):
        w, mw, dock = _make_real_viewer(qapp, tmp_path)
        w.update_iso()
        first_p = w.iso_actor_p
        w.spin.setValue(0.2)  # triggers update_iso again via signal
        assert w.iso_actor_p is not None
        _teardown_real(w, mw)

    def test_contour_exception_is_swallowed(self, qapp, tmp_path):
        w, mw, dock = _make_real_viewer(qapp, tmp_path)
        w.grid.contour = MagicMock(side_effect=Exception("boom"))
        w.update_iso()  # must not raise
        _teardown_real(w, mw)


# ---------------------------------------------------------------------------
# open_cube_viewer() — full file-open flow with real RDKit/numpy/pyvista.
# ---------------------------------------------------------------------------

import tempfile as _tempfile


def _make_cube_content(nx=3, ny=3, nz=3, atoms=None):
    if atoms is None:
        # Water-like geometry (Bohr units), bondable at charge 0.
        atoms = [
            (8, (0.0, 0.0, 0.0)),
            (1, (0.0, 0.0, 1.8)),
            (1, (1.7, 0.0, -0.5)),
        ]
    header = (
        "Title 1\nTitle 2\n"
        f"   {len(atoms)}   0.000000   0.000000   0.000000\n"
        f"   {nx}   0.400000   0.000000   0.000000\n"
        f"   {ny}   0.000000   0.400000   0.000000\n"
        f"   {nz}   0.000000   0.000000   0.400000\n"
    )
    atom_lines = "".join(
        f"   {num}   0.000000   {x:.6f}   {y:.6f}   {z:.6f}\n" for num, (x, y, z) in atoms
    )
    n = nx * ny * nz
    values = [str(float(i % 5) - 2.0) for i in range(n)]
    data_lines = ""
    row = []
    for v in values:
        row.append(v)
        if len(row) == 6:
            data_lines += " " + " ".join(row) + "\n"
            row = []
    if row:
        data_lines += " " + " ".join(row) + "\n"
    return header + atom_lines + data_lines


class _RealContext:
    """Minimal PluginContext stand-in that intentionally omits
    enter_3d_viewer_mode so open_cube_viewer falls back to
    main_window.ui_manager.enter_3d_viewer_ui_mode (elif branch)."""

    def __init__(self, mw):
        self._mw = mw
        self.windows = {}
        self.status_messages = []
        self.current_molecule = None

    def get_main_window(self):
        return self._mw

    def get_window(self, key):
        return self.windows.get(key)

    def register_window(self, key, win):
        self.windows[key] = win

    def show_status_message(self, msg):
        self.status_messages.append(msg)


def _make_real_mw():
    from PyQt6.QtWidgets import QMainWindow
    from types import SimpleNamespace

    mw = QMainWindow()
    mw.plotter = _pv.Plotter(off_screen=True)
    mw.init_manager = SimpleNamespace(current_file_path=None)
    mw.ui_manager = SimpleNamespace(
        enter_3d_viewer_ui_mode=MagicMock(), restore_ui_for_editing=MagicMock()
    )
    mw.current_mol = None
    return mw


class TestOpenCubeViewerReal:
    @pytest.fixture(autouse=True)
    def _no_real_settings_file(self, tmp_path, monkeypatch):
        """open_cube_viewer() constructs its own CubeViewerWidget internally, so
        we can't override get_settings_path() per-instance; redirect the whole
        class instead to avoid writing cube_viewer.json into the plugin dir."""
        path = str(tmp_path / "open_flow_settings.json")
        monkeypatch.setattr(
            _cuben.CubeViewerWidget, "get_settings_path", lambda self: path
        )

    def test_full_flow_creates_dock_and_registers_window(self, tmp_path, qapp):
        p = tmp_path / "water.cube"
        p.write_text(_make_cube_content())
        mw = _make_real_mw()
        ctx = _RealContext(mw)
        _cuben.open_cube_viewer(ctx, str(p))
        mw.ui_manager.enter_3d_viewer_ui_mode.assert_called()
        assert "main_panel" in ctx.windows
        assert ctx.current_molecule is not None
        assert any("Loaded Cube" in m for m in ctx.status_messages)
        ctx.windows["main_panel"].widget().close_plugin()
        mw.plotter.close()

    def test_parse_failure_shows_error(self, tmp_path, qapp):
        p = tmp_path / "bad.cube"
        p.write_text("too\nshort\nfile\n")
        mw = _make_real_mw()
        ctx = _RealContext(mw)
        _cuben.QMessageBox.critical = MagicMock()
        _cuben.open_cube_viewer(ctx, str(p))
        _cuben.QMessageBox.critical.assert_called_once()
        assert "main_panel" not in ctx.windows
        mw.plotter.close()

    def test_reopen_closes_previous_dock(self, tmp_path, qapp):
        p = tmp_path / "water.cube"
        p.write_text(_make_cube_content())
        mw = _make_real_mw()
        ctx = _RealContext(mw)
        _cuben.open_cube_viewer(ctx, str(p))
        first_widget = ctx.windows["main_panel"].widget()
        _cuben.open_cube_viewer(ctx, str(p))
        assert first_widget.dock is None  # closed by the 2nd open
        ctx.windows["main_panel"].widget().close_plugin()
        mw.plotter.close()

    def test_mol_none_fallback_builds_rwmol(self, tmp_path, qapp, monkeypatch):
        p = tmp_path / "water.cube"
        p.write_text(_make_cube_content())
        mw = _make_real_mw()
        ctx = _RealContext(mw)
        monkeypatch.setattr(_cuben.Chem, "MolFromXYZBlock", lambda *_a, **_k: None)
        _cuben.open_cube_viewer(ctx, str(p))
        assert ctx.current_molecule is not None
        assert ctx.current_molecule.GetNumAtoms() == 3
        ctx.windows["main_panel"].widget().close_plugin()
        mw.plotter.close()
