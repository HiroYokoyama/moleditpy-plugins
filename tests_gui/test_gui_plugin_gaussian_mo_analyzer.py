"""
Headless GUI tests for the Gaussian MO Analyzer plugin (package plugin).

Covers: initialize(), OrbitalWidget.
"""

from __future__ import annotations

import contextlib
import importlib
import os
import sys
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MO_PKG_DIR = PLUGINS_DIR / "Gaussian_MO_Analyzer"

with mock_chemistry_imports():
    sys.path.insert(0, str(MO_PKG_DIR))
    try:
        _mo = importlib.import_module("gaussian_fchk_mo_analyzer")
    finally:
        sys.path.remove(str(MO_PKG_DIR))


# ===========================================================================
# Gaussian MO Analyzer  (package plugin: OrbitalWidget + initialize)
# ===========================================================================


class TestMOAnalyzerInitialize:
    def test_registers_fchk_openers_and_drop_handler(self, qapp):
        ctx = MagicMock()
        _mo.initialize(ctx)
        exts = [c.args[0] for c in ctx.register_file_opener.call_args_list]
        assert exts == [".fchk", ".fck", ".fch"]
        ctx.register_drop_handler.assert_called_once()

    def test_drop_handler_accepts_only_fchk(self, qapp):
        ctx = MagicMock()
        _mo.initialize(ctx)
        handler = ctx.register_drop_handler.call_args.args[0]
        with patch.object(_mo, "OrbitalWidget") as ow:
            assert handler("mol.fchk") is True
            assert handler("mol.xyz") is False
            ow.assert_called_once()


class TestOrbitalWidget:
    @pytest.fixture
    def win(self, qapp, tmp_path):
        # Empty fchk → reader.data empty; basis prep patched out so the
        # widget builds with n_basis=0 (the "Error reading MOs" path).
        fchk = tmp_path / "water.fchk"
        fchk.write_text("dummy fchk\n", encoding="utf-8")
        ctx = MagicMock()
        with patch.object(
            _mo.analyzer.BasisSetEngine, "_prepare_basis_set", lambda self: None
        ):
            w = _mo.OrbitalWidget(None, ctx, str(fchk))
        yield w
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title_includes_filename(self, win):
        assert "Gaussian MO Analyzer" in win.windowTitle()
        assert "water.fchk" in win.windowTitle()

    def test_default_grid_points(self, win):
        assert win.spin_points.value() == 40

    def test_default_margin(self, win):
        assert win.spin_gap.value() == 4.0

    def test_default_isovalue(self, win):
        assert win.spin_iso.value() == 0.02

    def test_default_opacity(self, win):
        assert win.spin_opacity.value() == 0.4

    def test_smooth_shading_checked_by_default(self, win):
        assert win.check_smooth.isChecked()

    def test_progress_bar_hidden_initially(self, win):
        assert win.pbar.isHidden()

    def test_list_shows_error_without_mo_data(self, win):
        assert win.list_widget.count() == 1
        assert win.list_widget.item(0).text() == "Error reading MOs"

    def test_generate_with_no_selection_keeps_button_enabled(self, win):
        win.list_widget.clearSelection()
        win.generate_selected()
        assert win.btn_gen.isEnabled()

    def test_iso_change_without_cube_is_noop(self, win):
        win.last_cube_path = None
        win.on_iso_changed(0.05)  # must not raise

    def test_cube_path_layout(self, win):
        path = win.get_cube_path(7)
        assert path.endswith("water_MO7.cube")
        assert "water_cubes" in path

    def test_visualize_success_shows_status_message(self, win):
        with patch.object(_mo.gui, "CubeVisualizer") as VisCls:
            VisCls.return_value.load_file.return_value = True
            win.visualize("some.cube")
        VisCls.return_value.show_iso.assert_called_once()
        win.context.show_status_message.assert_called_once()

    def test_visualize_failure_shows_warning(self, win):
        with patch.object(_mo.gui, "CubeVisualizer") as VisCls, patch.object(
            _mo.gui, "QMessageBox"
        ) as MB:
            VisCls.return_value.load_file.return_value = False
            win.visualize("some.cube")
        MB.warning.assert_called_once()

    def test_on_iso_changed_with_existing_cube_revisualizes(self, win, tmp_path):
        cube = tmp_path / "x.cube"
        cube.write_text("data")
        win.last_cube_path = str(cube)
        win.visualize = MagicMock()
        win.on_iso_changed(0.05)
        win.visualize.assert_called_once_with(str(cube))

    def test_on_single_click_visualizes_when_cube_exists(self, win, tmp_path):
        cube_path = win.get_cube_path(1)
        os.makedirs(os.path.dirname(cube_path), exist_ok=True)
        with open(cube_path, "w") as f:
            f.write("data")
        item = MagicMock()
        item.data.return_value = 1
        win.visualize = MagicMock()
        win.on_single_click(item)
        assert win.last_cube_path == cube_path
        win.visualize.assert_called_once_with(cube_path)

    def test_on_single_click_noop_when_cube_missing(self, win):
        item = MagicMock()
        item.data.return_value = 999
        win.visualize = MagicMock()
        win.on_single_click(item)
        win.visualize.assert_not_called()

    def test_on_double_click_calls_generate_selected(self, win):
        win.generate_selected = MagicMock()
        win.on_double_click(MagicMock())
        win.generate_selected.assert_called_once()


# ===========================================================================
# OrbitalWidget — populate_list() / generate flow with real MO data.
#
# analyzer.py, gui.py and worker.py all do real numeric work (numpy); numpy
# is mocked package-wide in tests_gui, so real numpy is swapped into the
# already-imported submodules for the duration of these tests only (guarded
# with importorskip since tests_gui CI has no numpy installed).
# ===========================================================================


@contextlib.contextmanager
def _real_numpy_in_mo_modules():
    numpy = pytest.importorskip("numpy")
    saved = (_mo.analyzer.np, _mo.gui.np, _mo.worker.np)
    _mo.analyzer.np = _mo.gui.np = _mo.worker.np = numpy
    try:
        yield numpy
    finally:
        _mo.analyzer.np, _mo.gui.np, _mo.worker.np = saved


class _SyncCalcWorker(_mo.worker.CalcWorker):
    """CalcWorker whose start() runs synchronously (no real QThread spin-up)
    so tests can drive the full generate -> on_finished -> visualize flow
    deterministically; see project memory note on QThread.run() coverage."""

    def start(self):
        self.run()


_S_SHELL_MO_FCHK = textwrap.dedent(
    """\
    Atomic numbers                              I   N=           1
             1
    Current cartesian coordinates               R   N=           3
      0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
    Shell types                                 I   N=           1
             0
    Number of primitives per shell              I   N=           1
             1
    Shell to atom map                           I   N=           1
             1
    Primitive exponents                         R   N=           1
      1.000000000000E+00
    Contraction coefficients                    R   N=           1
      1.000000000000E+00
    Alpha MO coefficients                       R   N=           1
      1.000000000000E+00
    Alpha Orbital Energies                      R   N=           1
     -5.000000000000E-01
    Number of electrons                         I            2
    """
)


class TestOrbitalWidgetWithRealMOData:
    @pytest.fixture
    def win(self, qapp, tmp_path):
        with _real_numpy_in_mo_modules():
            fchk = tmp_path / "s.fchk"
            fchk.write_text(_S_SHELL_MO_FCHK)
            ctx = MagicMock()
            w = _mo.OrbitalWidget(None, ctx, str(fchk))
            w.context = ctx
            yield w
        w.destroy()

    def test_populate_list_uses_number_of_electrons_ceil_branch(self, win):
        # No "Number of alpha electrons" key -> n_occ from
        # ceil(Number of electrons / 2), exercising the np.ceil branch.
        assert win.list_widget.count() == 1
        assert "HOMO" in win.list_widget.item(0).text()

    def test_molecule_load_branch_reaches_reset_camera(self, win):
        # Constructor's rdkit/molecule-loading branch (real xyz block) runs
        # without raising and calls through to reset_3d_camera().
        win.context.reset_3d_camera.assert_called_once()

    def test_generate_selected_writes_real_cube_and_visualizes(self, win, tmp_path):
        win.spin_points.setValue(3)
        with patch.object(_mo.gui, "CalcWorker", _SyncCalcWorker), patch.object(
            _mo.gui, "CubeVisualizer"
        ) as VisCls:
            VisCls.return_value.load_file.return_value = True
            win.list_widget.setCurrentRow(0)
            win.list_widget.item(0).setSelected(True)
            win.generate_selected()
        assert os.path.exists(win.get_cube_path(1))
        VisCls.return_value.show_iso.assert_called()
        assert win.btn_gen.isEnabled()
        assert not win.pbar.isVisible()

    def test_generate_selected_worker_failure_shows_warning(self, win, tmp_path):
        # Corrupt the MO coefficients so evaluate_mo_on_grid raises inside
        # CalcWorker.run() -> finished_sig emits (False, msg).
        win.reader.data["Alpha MO coefficients"] = None
        win.spin_points.setValue(3)
        with patch.object(_mo.gui, "CalcWorker", _SyncCalcWorker), patch.object(
            _mo.gui, "QMessageBox"
        ) as MB:
            win.list_widget.setCurrentRow(0)
            win.list_widget.item(0).setSelected(True)
            win.generate_selected()
        MB.warning.assert_called_once()

    def test_generate_selected_skips_existing_cube(self, win, tmp_path):
        cube_path = win.get_cube_path(1)
        os.makedirs(os.path.dirname(cube_path), exist_ok=True)
        with open(cube_path, "w") as f:
            f.write("existing")
        with patch.object(_mo.gui, "CubeVisualizer") as VisCls:
            VisCls.return_value.load_file.return_value = True
            win.list_widget.setCurrentRow(0)
            win.list_widget.item(0).setSelected(True)
            win.generate_selected()
        VisCls.return_value.show_iso.assert_called()
        assert win.btn_gen.isEnabled()


class TestOrbitalWidgetAlphaElectronsBranch:
    @pytest.fixture
    def win(self, qapp, tmp_path):
        with _real_numpy_in_mo_modules():
            fchk_text = _S_SHELL_MO_FCHK + (
                "Number of alpha electrons                   I            1\n"
            )
            fchk = tmp_path / "s2.fchk"
            fchk.write_text(fchk_text)
            ctx = MagicMock()
            w = _mo.OrbitalWidget(None, ctx, str(fchk))
            yield w
        w.destroy()

    def test_populate_list_uses_alpha_electrons_directly(self, win):
        # "Number of alpha electrons" present -> n_occ = n_alpha[0] path
        # (no np.ceil call).
        assert win.list_widget.count() == 1
        assert "HOMO" in win.list_widget.item(0).text()


# ===========================================================================
# CalcWorker.run() — real QThread subclass (real PyQt6), run synchronously
# so coverage attributes to worker.py directly instead of a QThread.start().
# ===========================================================================


class TestCalcWorkerRealRun:
    def _engine_and_reader(self, tmp_path):
        fchk = tmp_path / "s.fchk"
        fchk.write_text(_S_SHELL_MO_FCHK)
        reader = _mo.analyzer.FCHKReader(str(fchk))
        engine = _mo.analyzer.BasisSetEngine(reader)
        return engine, reader

    def test_run_mo_mode_success(self, qapp, tmp_path):
        with _real_numpy_in_mo_modules() as numpy:
            engine, _ = self._engine_and_reader(tmp_path)
            out = tmp_path / "out.cube"
            worker = _mo.worker.CalcWorker(
                engine,
                0,
                3,
                2.0,
                numpy.array([[0.0, 0.0, 0.0]]),
                numpy.array([1.0]),
                str(out),
                mode="MO",
            )
            results = []
            worker.finished_sig.connect(lambda ok, msg: results.append((ok, msg)))
            worker.run()
        assert results == [(True, str(out))]
        assert out.exists()

    def test_run_basis_mode_success(self, qapp, tmp_path):
        with _real_numpy_in_mo_modules() as numpy:
            engine, _ = self._engine_and_reader(tmp_path)
            out = tmp_path / "basis.cube"
            worker = _mo.worker.CalcWorker(
                engine, 0, 3, 2.0, numpy.array([[0.0, 0.0, 0.0]]), None, str(out),
                mode="Basis",
            )
            results = []
            worker.finished_sig.connect(lambda ok, msg: results.append((ok, msg)))
            worker.run()
        assert results == [(True, str(out))]

    def test_run_cancelled_emits_nothing(self, qapp, tmp_path):
        with _real_numpy_in_mo_modules() as numpy:
            engine, _ = self._engine_and_reader(tmp_path)
            worker = _mo.worker.CalcWorker(
                engine,
                0,
                3,
                2.0,
                numpy.array([[0.0, 0.0, 0.0]]),
                numpy.array([1.0]),
                str(tmp_path / "cancelled.cube"),
            )
            worker._is_cancelled = True
            results = []
            worker.finished_sig.connect(lambda ok, msg: results.append((ok, msg)))
            worker.run()
        assert results == []

    def test_run_exception_emits_false(self, qapp, tmp_path):
        with _real_numpy_in_mo_modules() as numpy:
            engine, _ = self._engine_and_reader(tmp_path)
            worker = _mo.worker.CalcWorker(
                engine,
                0,
                3,
                2.0,
                numpy.array([[0.0, 0.0, 0.0]]),
                None,  # mode="MO" with coeffs=None -> TypeError inside engine
                str(tmp_path / "err.cube"),
            )
            results = []
            worker.finished_sig.connect(lambda ok, msg: results.append((ok, msg)))
            worker.run()
        assert len(results) == 1
        assert results[0][0] is False
