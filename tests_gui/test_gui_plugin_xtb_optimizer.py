"""
Headless GUI tests for the xTB Optimizer plugin.

Covers: XtbOptimizerDialog, XtbWorker.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_xtb_optimizer.py
"""

from __future__ import annotations

import math
import sys
import types
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

XTB_PATH = PLUGINS_DIR / "XTB_Optimizer" / "xtb_optimizer.py"

with mock_chemistry_imports():
    _xtb = load_plugin_for_gui(XTB_PATH)


def _ctx(setting="MMFF_RDKIT") -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.get_setting.return_value = setting
    return ctx


@pytest.fixture
def no_msgbox(monkeypatch):
    """Replace the module's QMessageBox so warnings/info don't block offscreen."""
    box = MagicMock()
    monkeypatch.setattr(_xtb, "QMessageBox", box)
    return {_xtb.__name__: box}


# ===========================================================================
# XtbOptimizerDialog  (visible plugin: "xTB Optimizer")
# ===========================================================================


class TestXtbOptimizerDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _xtb.XtbOptimizerDialog(context=_ctx(), parent=None)
        yield d
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "xTB Geometry Optimizer"

    def test_method_options(self, dlg):
        texts = [dlg.combo_method.itemText(i) for i in range(dlg.combo_method.count())]
        assert texts == ["GFN2-xTB", "GFN1-xTB"]

    def test_convergence_defaults(self, dlg):
        assert dlg.spin_fmax.value() == pytest.approx(0.05)
        assert dlg.spin_maxsteps.value() == 500

    def test_table_columns(self, dlg):
        assert dlg.table.columnCount() == 3
        headers = [dlg.table.horizontalHeaderItem(i).text() for i in range(3)]
        assert headers == ["Step", "Energy (eV)", "Fmax (eV/Å)"]

    def test_log_is_readonly(self, dlg):
        assert dlg.log_text.isReadOnly()

    def test_initial_state(self, dlg):
        assert dlg.lbl_status.text() == "Ready."
        assert dlg.btn_run.isEnabled()
        assert not dlg.btn_cancel.isEnabled()
        assert dlg.progress_bar.isHidden()

    def test_set_running_toggles_controls(self, dlg):
        dlg._set_running(True)
        assert not dlg.btn_run.isEnabled()
        assert dlg.btn_cancel.isEnabled()
        assert not dlg.combo_method.isEnabled()
        assert not dlg.progress_bar.isHidden()
        dlg._set_running(False)
        assert dlg.btn_run.isEnabled()
        assert not dlg.btn_cancel.isEnabled()

    def test_append_log(self, dlg):
        dlg._append_log("line one")
        dlg._append_log("line two")
        assert "line one" in dlg.log_text.toPlainText()
        assert "line two" in dlg.log_text.toPlainText()

    def test_step_update_fills_table(self, dlg):
        dlg._on_step_update(1, -123.456789, 0.0567)
        assert dlg.table.rowCount() == 1
        assert dlg.table.item(0, 0).text() == "1"
        assert dlg.table.item(0, 1).text() == "-123.456789"
        assert dlg.table.item(0, 2).text() == "0.0567"
        assert "Step 1" in dlg.lbl_status.text()

    def test_run_without_molecule_warns(self, dlg, no_msgbox):
        dlg._on_run()
        no_msgbox[_xtb.__name__].warning.assert_called_once()
        assert dlg.btn_run.isEnabled()

    def test_run_without_conformer_warns(self, dlg, no_msgbox):
        mol = MagicMock()
        mol.GetNumConformers.return_value = 0
        dlg.context.current_mol = mol
        dlg._on_run()
        msg = no_msgbox[_xtb.__name__].warning.call_args[0][2]
        assert "3D" in msg

    def test_finished_failure_reports_status(self, dlg):
        dlg._set_running(True)
        dlg._on_finished(False, "Cancelled by user.")
        assert dlg.lbl_status.text() == "Stopped: Cancelled by user."
        assert dlg.btn_run.isEnabled()

    def test_finished_success_applies_coordinates(self, dlg):
        mol = MagicMock()
        dlg.context.current_mol = mol
        dlg._step_count = 12
        dlg._on_finished(True, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        conf = mol.GetConformer.return_value
        assert conf.SetAtomPosition.call_count == 2
        dlg.context.refresh_3d_view.assert_called_once()
        dlg.context.push_undo_checkpoint.assert_called_once()
        assert "complete" in dlg.lbl_status.text()

    def test_finished_success_with_unloaded_molecule(self, dlg):
        dlg.context.current_mol = None
        dlg._on_finished(True, [[0.0, 0.0, 0.0]])
        assert "unloaded" in dlg.lbl_status.text()

    def test_accept_unregisters_window(self, dlg):
        dlg.accept()
        dlg.context.register_window.assert_called_with("main_panel", None)

    def test_worker_cancel_sets_flag(self, qapp):
        w = _xtb.XtbWorker(
            numbers=[1, 1], positions=[[0, 0, 0], [0.74, 0, 0]],
            method="GFN2-xTB", fmax=0.05, max_steps=10,
        )
        assert not w._cancelled
        w.cancel()
        assert w._cancelled

    def test_run_plugin_without_molecule_warns(self, no_msgbox):
        ctx = _ctx()
        _xtb.run_plugin(ctx)
        no_msgbox[_xtb.__name__].warning.assert_called_once()

    def test_initialize_registers_menu_action(self):
        ctx = _ctx()
        _xtb.initialize(ctx)
        assert ctx.add_menu_action.call_args[0][0] == "3D Edit/xTB Optimizer…"


# ===========================================================================
# XtbWorker.run() — fake ase/tblite/numpy stacks (real QThread.run() called
# synchronously, i.e. by direct call rather than .start(), for coverage
# attribution).
# ===========================================================================


class _FVec:
    """Minimal stand-in for a numpy array supporting the ops run() needs."""

    def __init__(self, val):
        self.val = val

    def __pow__(self, power):
        return _FVec(self.val**power)

    def sum(self, axis=None):
        return self

    def max(self):
        return self.val


class _FakePositions(list):
    def tolist(self):
        return [list(row) for row in self]


class _FakeAtoms:
    def __init__(self, numbers=None, positions=None):
        self.numbers = numbers
        self._positions = _FakePositions([list(p) for p in positions])
        self.calc = None
        self.forces_raise = False
        self.energy_raise = False

    def get_potential_energy(self):
        if self.energy_raise:
            raise RuntimeError("energy eval failed")
        return -42.0

    def get_forces(self):
        if self.forces_raise:
            raise RuntimeError("forces eval failed")
        return _FVec(2.0)

    def get_positions(self):
        return self._positions


class _FakeLBFGS:
    """Fake optimizer: calls the attached callback ``n_steps`` times then
    returns ``converged``. If ``raise_on_step`` is set, the callback error
    (e.g. cancellation) propagates out of run(), mirroring real ASE."""

    def __init__(self, atoms, logfile=None, *, n_steps=2, converged=True):
        self.atoms = atoms
        self._cb = None
        self.n_steps = n_steps
        self.converged = converged

    def attach(self, cb):
        self._cb = cb

    def run(self, fmax=None, steps=None):
        for _ in range(self.n_steps):
            self._cb()
        return self.converged


def _install_fake_ase_stack(monkeypatch, lbfgs_cls, atoms_cls=_FakeAtoms):
    """Patch sys.modules so XtbWorker.run()'s local imports resolve to fakes."""
    fake_numpy = types.ModuleType("numpy")
    fake_numpy.sqrt = math.sqrt
    monkeypatch.setitem(sys.modules, "numpy", fake_numpy)

    fake_ase = types.ModuleType("ase")
    fake_ase.Atoms = atoms_cls
    fake_ase_optimize = types.ModuleType("ase.optimize")
    fake_ase_optimize.LBFGS = lbfgs_cls
    fake_ase.optimize = fake_ase_optimize
    monkeypatch.setitem(sys.modules, "ase", fake_ase)
    monkeypatch.setitem(sys.modules, "ase.optimize", fake_ase_optimize)

    fake_tblite = types.ModuleType("tblite")
    fake_tblite_ase = types.ModuleType("tblite.ase")
    fake_tblite_ase.TBLite = lambda **kw: MagicMock()
    fake_tblite.ase = fake_tblite_ase
    monkeypatch.setitem(sys.modules, "tblite", fake_tblite)
    monkeypatch.setitem(sys.modules, "tblite.ase", fake_tblite_ase)


def _make_worker(**overrides):
    kwargs = dict(
        numbers=[6, 1],
        positions=[[0.0, 0.0, 0.0], [1.09, 0.0, 0.0]],
        method="GFN2-xTB",
        fmax=0.05,
        max_steps=10,
    )
    kwargs.update(overrides)
    return _xtb.XtbWorker(**kwargs)


class TestXtbWorkerRun:
    def test_run_converges(self, qapp, monkeypatch):
        _install_fake_ase_stack(
            monkeypatch, lambda *a, **k: _FakeLBFGS(*a, n_steps=3, converged=True, **k)
        )
        w = _make_worker()
        steps = []
        finished = {}
        w.step_update.connect(lambda s, e, f: steps.append((s, e, f)))
        w.finished.connect(lambda ok, payload: finished.update(ok=ok, payload=payload))
        w.run()
        assert len(steps) == 3
        assert finished["ok"] is True
        assert finished["payload"] == [[0.0, 0.0, 0.0], [1.09, 0.0, 0.0]]

    def test_run_reaches_max_steps_without_convergence(self, qapp, monkeypatch):
        _install_fake_ase_stack(
            monkeypatch, lambda *a, **k: _FakeLBFGS(*a, n_steps=2, converged=False, **k)
        )
        w = _make_worker()
        logs = []
        finished = {}
        w.log_message.connect(logs.append)
        w.finished.connect(lambda ok, payload: finished.update(ok=ok, payload=payload))
        w.run()
        assert finished["ok"] is True
        assert any("max steps" in m for m in logs)

    def test_run_step_callback_handles_energy_forces_exception(self, qapp, monkeypatch):
        atoms_cls = lambda **kw: _FakeAtoms(**kw)  # noqa: E731

        def _atoms_ctor(**kw):
            a = _FakeAtoms(**kw)
            a.forces_raise = True
            return a

        _install_fake_ase_stack(
            monkeypatch,
            lambda *a, **k: _FakeLBFGS(*a, n_steps=1, converged=True, **k),
            atoms_cls=_atoms_ctor,
        )
        w = _make_worker()
        steps = []
        w.step_update.connect(lambda s, e, f: steps.append((s, e, f)))
        w.run()
        # forces raised -> callback falls back to NaN energy/fmax
        assert len(steps) == 1
        assert math.isnan(steps[0][1]) and math.isnan(steps[0][2])

    def test_run_cancelled_mid_step(self, qapp, monkeypatch):
        _install_fake_ase_stack(
            monkeypatch, lambda *a, **k: _FakeLBFGS(*a, n_steps=1, converged=True, **k)
        )
        w = _make_worker()
        w._cancelled = True  # simulate cancel() called before the callback fires
        finished = {}
        logs = []
        w.log_message.connect(logs.append)
        w.finished.connect(lambda ok, payload: finished.update(ok=ok, payload=payload))
        w.run()
        assert finished == {"ok": False, "payload": "Cancelled by user."}
        assert any("cancelled" in m.lower() for m in logs)

    def test_run_runtime_error_not_cancellation(self, qapp, monkeypatch):
        class _BoomLBFGS(_FakeLBFGS):
            def run(self, fmax=None, steps=None):
                raise RuntimeError("SCF did not converge")

        _install_fake_ase_stack(monkeypatch, _BoomLBFGS)
        w = _make_worker()
        finished = {}
        w.finished.connect(lambda ok, payload: finished.update(ok=ok, payload=payload))
        w.run()
        assert finished["ok"] is False
        assert "SCF did not converge" in finished["payload"]

    def test_run_unexpected_exception(self, qapp, monkeypatch):
        def _boom_atoms(**kw):
            raise ValueError("bad atom data")

        _install_fake_ase_stack(monkeypatch, _FakeLBFGS, atoms_cls=_boom_atoms)
        w = _make_worker()
        finished = {}
        w.finished.connect(lambda ok, payload: finished.update(ok=ok, payload=payload))
        w.run()
        assert finished["ok"] is False
        assert "bad atom data" in finished["payload"]

    def test_run_missing_dependency_reports_import_error(self, qapp, monkeypatch):
        """tblite is never installed in this workspace; run() must hit the
        ImportError branch and report a helpful message."""
        monkeypatch.delitem(sys.modules, "tblite", raising=False)
        monkeypatch.delitem(sys.modules, "tblite.ase", raising=False)
        w = _make_worker()
        finished = {}
        w.finished.connect(lambda ok, payload: finished.update(ok=ok, payload=payload))
        w.run()
        assert finished["ok"] is False
        assert "Missing dependency" in finished["payload"]


# ===========================================================================
# XtbOptimizerDialog — _on_run / _on_cancel / closeEvent / accept / run_plugin
# ===========================================================================


class TestXtbOptimizerDialogFullLifecycle:
    @pytest.fixture
    def dlg(self, qapp):
        d = _xtb.XtbOptimizerDialog(context=_ctx(), parent=None)
        yield d
        d.destroy()

    def _mol_with_atoms(self, symbols):
        atoms = []
        for i, sym in enumerate(symbols):
            a = MagicMock()
            a.GetSymbol.return_value = sym
            a.GetIdx.return_value = i
            atoms.append(a)
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        mol.GetAtoms.return_value = atoms
        conf = MagicMock()
        conf.GetAtomPosition.side_effect = [
            MagicMock(x=float(i), y=0.0, z=0.0) for i in range(len(symbols))
        ]
        mol.GetConformer.return_value = conf
        return mol

    @staticmethod
    def _install_fake_rdkit_chem(monkeypatch):
        """rdkit is genuinely absent in CI; _on_run() imports it locally, so
        provide a minimal fake so the success/dummy-atom paths are reachable
        without relying on a real (or ImportError-catching-critical-dialog)
        rdkit install."""
        table = {"H": 1, "C": 6, "N": 7, "O": 8}

        class _FakePT:
            def GetAtomicNumber(self, symbol):
                return table.get(symbol, 0)

        fake_rdkit = types.ModuleType("rdkit")
        fake_chem = types.ModuleType("rdkit.Chem")
        fake_chem.GetPeriodicTable = lambda: _FakePT()
        fake_rdkit.Chem = fake_chem
        monkeypatch.setitem(sys.modules, "rdkit", fake_rdkit)
        monkeypatch.setitem(sys.modules, "rdkit.Chem", fake_chem)

    def test_on_run_gathers_atoms_and_starts_worker(self, dlg, monkeypatch):
        self._install_fake_rdkit_chem(monkeypatch)
        started = []
        monkeypatch.setattr(_xtb.XtbWorker, "start", lambda self: started.append(self))
        dlg.context.current_mol = self._mol_with_atoms(["C", "H"])
        dlg._on_run()
        assert started
        assert dlg._worker is not None
        assert not dlg.btn_run.isEnabled()
        assert dlg.btn_cancel.isEnabled()
        assert "Running" in dlg.lbl_status.text()

    def test_on_run_dummy_atom_warns_and_skips_worker(self, dlg, monkeypatch, no_msgbox):
        self._install_fake_rdkit_chem(monkeypatch)
        started = []
        monkeypatch.setattr(_xtb.XtbWorker, "start", lambda self: started.append(self))
        dlg.context.current_mol = self._mol_with_atoms(["C", "*"])
        dlg._on_run()
        no_msgbox[_xtb.__name__].warning.assert_called_once()
        assert not started
        assert dlg._worker is None

    def test_on_run_read_failure_shows_critical(self, dlg, no_msgbox):
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        mol.GetAtoms.side_effect = RuntimeError("corrupt molecule")
        dlg.context.current_mol = mol
        dlg._on_run()
        no_msgbox[_xtb.__name__].critical.assert_called_once()

    def test_on_cancel_requests_worker_cancellation(self, dlg):
        worker = MagicMock()
        worker.isRunning.return_value = True
        dlg._worker = worker
        dlg.btn_cancel.setEnabled(True)
        dlg._on_cancel()
        worker.cancel.assert_called_once()
        assert dlg.lbl_status.text() == "Cancelling…"
        assert not dlg.btn_cancel.isEnabled()

    def test_on_cancel_noop_without_worker(self, dlg):
        dlg._worker = None
        dlg._on_cancel()  # must not raise

    def test_finished_apply_exception_reports_status(self, dlg):
        mol = MagicMock()
        conf = mol.GetConformer.return_value
        conf.SetAtomPosition.side_effect = ValueError("boom")
        dlg.context.current_mol = mol
        dlg._on_finished(True, [[0.0, 0.0, 0.0]])
        assert "Error applying coordinates" in dlg.lbl_status.text()

    def test_close_event_no_worker_unregisters_and_accepts(self, dlg):
        event = MagicMock()
        dlg._worker = None
        dlg.closeEvent(event)
        dlg.context.register_window.assert_called_with("main_panel", None)
        event.accept.assert_called_once()

    def test_close_event_running_worker_confirm_yes_cancels_and_waits(self, dlg, no_msgbox):
        box = no_msgbox[_xtb.__name__]
        box.StandardButton.Yes = 1
        box.question.return_value = 1
        worker = MagicMock()
        worker.isRunning.return_value = True
        dlg._worker = worker
        event = MagicMock()
        dlg.closeEvent(event)
        worker.cancel.assert_called_once()
        worker.wait.assert_called_once_with(3000)
        assert event.accept.called  # accept() is called from both branches

    def test_close_event_running_worker_confirm_no_ignores(self, dlg, no_msgbox):
        box = no_msgbox[_xtb.__name__]
        box.StandardButton.Yes = 1
        box.StandardButton.No = 2
        box.question.return_value = 2
        worker = MagicMock()
        worker.isRunning.return_value = True
        dlg._worker = worker
        event = MagicMock()
        dlg.closeEvent(event)
        worker.cancel.assert_not_called()
        event.ignore.assert_called_once()
        event.accept.assert_not_called()

    def test_accept_defers_to_close_event_when_running(self, dlg, no_msgbox):
        box = no_msgbox[_xtb.__name__]
        box.StandardButton.Yes = 1
        box.StandardButton.No = 2
        box.question.return_value = 2
        worker = MagicMock()
        worker.isRunning.return_value = True
        dlg._worker = worker
        dlg.accept()  # must not raise; delegates to closeEvent's fake-event branch
        worker.cancel.assert_not_called()


class TestRunPluginCreatesDialog:
    def test_run_plugin_creates_and_registers_new_dialog(self, qapp):
        ctx = _ctx()
        ctx.current_mol = MagicMock()
        ctx.get_window.return_value = None
        _xtb.run_plugin(ctx)
        ctx.register_window.assert_called_once()
        key, win = ctx.register_window.call_args[0]
        assert key == "main_panel"
        assert isinstance(win, _xtb.XtbOptimizerDialog)
        win.destroy()
