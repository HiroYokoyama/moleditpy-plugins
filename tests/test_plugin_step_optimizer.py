"""
Tests for the Step Optimizer plugin.

Covers:
  1. module constants + version format
  2. initialize() registers exactly one menu action
  3. run_plugin no-molecule guard / singleton reuse
  4. _tick logic: converged, non-converged, molecule-changed guard
  5. _start logic: MMFF props missing, UFF path, ff.Initialize called
  6. _stop / reject behavior
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
STEP_OPT_PATH = PLUGINS_DIR / "Step_Optimizer" / "step_optimizer.py"


# ---------------------------------------------------------------------------
# Module constants
# ---------------------------------------------------------------------------


class TestStepOptimizerConstants:
    def test_constants_present(self):
        with mock_optional_imports():
            mod = load_plugin(STEP_OPT_PATH)
            assert mod.PLUGIN_NAME == "Step Optimizer"
            assert mod.PLUGIN_AUTHOR == "HiroYokoyama"
            assert mod.PLUGIN_SUPPORTED_MOLEDITPY_VERSION == ">=4.0.0, <5.0.0"
            assert mod.PLUGIN_DESCRIPTION

    def test_version_format(self):
        with mock_optional_imports():
            mod = load_plugin(STEP_OPT_PATH)
            assert re.match(r"^\d{4}\.\d{2}\.\d{2}$", mod.PLUGIN_VERSION)


# ---------------------------------------------------------------------------
# initialize() / menu registration
# ---------------------------------------------------------------------------


class TestStepOptimizerInitialize:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(STEP_OPT_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()
            call_args = ctx.add_menu_action.call_args
            assert call_args[0][0] == "3D Edit/Step Optimizer..."


# ---------------------------------------------------------------------------
# run_plugin
# ---------------------------------------------------------------------------


class TestStepOptimizerRunPlugin:
    def test_run_plugin_no_mol_warns(self):
        with mock_optional_imports():
            mod = load_plugin(STEP_OPT_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            mock_warn.assert_called_once()
            ctx.register_window.assert_not_called()

    def test_run_plugin_creates_and_registers_dialog(self):
        with mock_optional_imports():
            mod = load_plugin(STEP_OPT_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            ctx.current_mol = MagicMock()
            mod.run_plugin(ctx)
            ctx.register_window.assert_called_once()

    def test_run_plugin_reuses_existing_window(self):
        with mock_optional_imports():
            mod = load_plugin(STEP_OPT_PATH)
            ctx = make_context()
            ctx.current_mol = MagicMock()
            win = MagicMock()
            ctx.get_window.return_value = win
            mod.run_plugin(ctx)
            win.show.assert_called_once()
            win.raise_.assert_called_once()
            win.activateWindow.assert_called_once()
            ctx.register_window.assert_not_called()


# ---------------------------------------------------------------------------
# legacy run(mw)
# ---------------------------------------------------------------------------


class TestStepOptimizerLegacyRun:
    def test_run_invokes_launch_fn_when_set(self):
        with mock_optional_imports():
            mod = load_plugin(STEP_OPT_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            mod._launch_fn = MagicMock()
            mod.run(SimpleNamespace())
            mod._launch_fn.assert_called_once()

    def test_run_noop_when_launch_fn_none(self):
        with mock_optional_imports():
            mod = load_plugin(STEP_OPT_PATH)
            mod._launch_fn = None
            mod.run(SimpleNamespace())  # must not raise

    def test_run_unwraps_host_attribute(self):
        with mock_optional_imports():
            mod = load_plugin(STEP_OPT_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            mod._launch_fn = MagicMock()
            mw = SimpleNamespace(host=SimpleNamespace())
            mod.run(mw)
            mod._launch_fn.assert_called_once()


# ---------------------------------------------------------------------------
# _tick
# ---------------------------------------------------------------------------


class _FakeFF:
    def __init__(self, minimize_result=0, energy=1.2345):
        self.minimize_result = minimize_result
        self.energy = energy
        self.minimize_calls = []
        self.initialize_called = False

    def Initialize(self):
        self.initialize_called = True

    def Minimize(self, maxIts=None):
        self.minimize_calls.append(maxIts)
        return self.minimize_result

    def CalcEnergy(self):
        return self.energy


class _FakePoint:
    """Plain (x, y, z) holder mimicking rdkit.Geometry.Point3D's attributes."""

    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


def _fake_point3d(x, y, z):
    return _FakePoint(x, y, z)


def _tick_fn():
    return extract_function(
        STEP_OPT_PATH,
        "StepOptimizerDialog",
        "_tick",
        {
            "logging": MagicMock(),
            "math": math,
            "Geometry": SimpleNamespace(Point3D=_fake_point3d),
        },
    )


class _FakeConf:
    """Fake rdkit Conformer: positions list + SetAtomPosition call recording."""

    def __init__(self, positions):
        self.positions = list(positions)
        self.set_calls = []

    def GetAtomPosition(self, i):
        return self.positions[i]

    def SetAtomPosition(self, i, point):
        self.set_calls.append((i, point))
        self.positions[i] = point


class _FakeTickMol:
    def __init__(self, conf):
        self._conf = conf

    def GetConformer(self):
        return self._conf

    def GetNumAtoms(self):
        return len(self._conf.positions)


def _tick_self(mol, ff, steps=5, current_mol=None, max_move=0.0):
    if current_mol is None:
        current_mol = mol
    return SimpleNamespace(
        context=SimpleNamespace(
            current_mol=current_mol,
            refresh_3d_view=MagicMock(),
            push_undo_checkpoint=MagicMock(),
        ),
        target_mol=mol,
        ff=ff,
        spin_steps=MagicMock(value=lambda: steps),
        spin_max_move=MagicMock(value=lambda: max_move),
        step_count=0,
        steps_since_checkpoint=0,
        lbl_step=MagicMock(),
        lbl_energy=MagicMock(),
        lbl_state=MagicMock(),
        timer=MagicMock(),
        running=True,
        _pause_timer=None,
        _reenable_controls=None,
    )


def _wire_helpers(self_):
    """Wire the real (non-extracted) helper methods onto a SimpleNamespace self."""

    def pause_timer():
        self_.timer.stop()
        self_.running = False

    def reenable():
        self_.btn_start = getattr(self_, "btn_start", MagicMock())
        self_.btn_stop = getattr(self_, "btn_stop", MagicMock())
        self_.combo_ff = getattr(self_, "combo_ff", MagicMock())
        self_.btn_start.setEnabled(True)
        self_.btn_stop.setEnabled(False)
        self_.combo_ff.setEnabled(True)

    self_._pause_timer = pause_timer
    self_._reenable_controls = reenable
    return self_


class TestStepOptimizerTick:
    # Continuous-run redesign (2026-07-12): the run must keep going when
    # Minimize() reports convergence (res == 0), so the user can pull atoms
    # mid-run and watch the structure re-relax. The timer is only ever
    # stopped by Stop, Close, the molecule-changed guard, or an error - never
    # by convergence - and no undo checkpoint is pushed from _tick anymore.
    def test_converged_minimize_does_not_stop_the_run(self):
        fn = _tick_fn()
        mol = MagicMock()
        ff = _FakeFF(minimize_result=0, energy=3.14159)
        self_ = _wire_helpers(_tick_self(mol, ff, steps=5))
        fn(self_)
        self_.timer.stop.assert_not_called()
        self_.context.push_undo_checkpoint.assert_not_called()
        assert self_.step_count == 5
        assert self_.running is True

    def test_non_converged_keeps_running_and_refreshes_view(self):
        fn = _tick_fn()
        mol = MagicMock()
        ff = _FakeFF(minimize_result=1, energy=10.0)
        self_ = _wire_helpers(_tick_self(mol, ff, steps=3))
        fn(self_)
        self_.timer.stop.assert_not_called()
        self_.context.push_undo_checkpoint.assert_not_called()
        self_.context.refresh_3d_view.assert_called_once()
        assert self_.step_count == 3

    def test_state_label_untouched_regardless_of_convergence(self):
        # No "converged" indication anywhere: the state label is left as
        # whatever _start already set ("State: Running") for the whole run,
        # on both a converged (res == 0) and non-converged (res == 1) tick.
        for res in (0, 1):
            fn = _tick_fn()
            mol = MagicMock()
            ff = _FakeFF(minimize_result=res, energy=1.0)
            self_ = _wire_helpers(_tick_self(mol, ff, steps=1))
            fn(self_)
            self_.lbl_state.setText.assert_not_called()

    def test_molecule_changed_guard_stops_without_minimize(self):
        fn = _tick_fn()
        mol = MagicMock()
        other_mol = MagicMock()
        ff = _FakeFF()
        self_ = _wire_helpers(_tick_self(mol, ff, steps=3, current_mol=other_mol))
        fn(self_)
        assert ff.minimize_calls == []
        self_.timer.stop.assert_called_once()
        assert "Molecule changed" in self_.lbl_state.setText.call_args[0][0]

    def test_energy_label_formatted_to_4_decimals(self):
        fn = _tick_fn()
        mol = MagicMock()
        ff = _FakeFF(minimize_result=1, energy=1.23456789)
        self_ = _wire_helpers(_tick_self(mol, ff, steps=1))
        fn(self_)
        msg = self_.lbl_energy.setText.call_args[0][0]
        assert "1.2346" in msg

    def test_exception_during_tick_is_caught_and_logged(self):
        fn = _tick_fn()
        mol = MagicMock()

        class _RaisingFF(_FakeFF):
            def Minimize(self, maxIts=None):
                raise RuntimeError("boom")

        ff = _RaisingFF()
        self_ = _wire_helpers(_tick_self(mol, ff, steps=1))
        fn(self_)  # must not raise
        self_.timer.stop.assert_called_once()
        assert "Error" in self_.lbl_state.setText.call_args[0][0]


class _MovingFF(_FakeFF):
    """Fake force field whose Minimize() actually displaces conformer atoms,
    so the max-move clamp has real displacement vectors to rescale."""

    def __init__(self, conf, moves, **kwargs):
        super().__init__(**kwargs)
        self.conf = conf
        self.moves = moves  # list of (dx, dy, dz) applied to each atom

    def Minimize(self, maxIts=None):
        self.minimize_calls.append(maxIts)
        for i, (dx, dy, dz) in enumerate(self.moves):
            p = self.conf.positions[i]
            self.conf.positions[i] = _FakePoint(p.x + dx, p.y + dy, p.z + dz)
        return self.minimize_result


class TestStepOptimizerMaxMoveClamp:
    def test_large_displacement_rescaled_to_exactly_max_move(self):
        conf = _FakeConf([_FakePoint(0.0, 0.0, 0.0)])
        mol = _FakeTickMol(conf)
        # atom 0 moves 3.0 A along x; clamp is 0.5 A
        ff = _MovingFF(conf, moves=[(3.0, 0.0, 0.0)])
        self_ = _wire_helpers(_tick_self(mol, ff, steps=1, max_move=0.5))
        fn = _tick_fn()
        fn(self_)

        assert len(conf.set_calls) == 1
        i, point = conf.set_calls[0]
        assert i == 0
        dist = math.sqrt(point.x**2 + point.y**2 + point.z**2)
        assert math.isclose(dist, 0.5, rel_tol=1e-9, abs_tol=1e-9)

    def test_small_displacement_left_untouched(self):
        conf = _FakeConf([_FakePoint(0.0, 0.0, 0.0)])
        mol = _FakeTickMol(conf)
        # atom 0 moves only 0.1 A; clamp is 0.5 A -> no rescale needed
        ff = _MovingFF(conf, moves=[(0.1, 0.0, 0.0)])
        self_ = _wire_helpers(_tick_self(mol, ff, steps=1, max_move=0.5))
        fn = _tick_fn()
        fn(self_)

        assert conf.set_calls == []
        assert conf.positions[0].x == pytest.approx(0.1)

    def test_clamp_off_skips_snapshot_and_set_position(self):
        conf = _FakeConf([_FakePoint(0.0, 0.0, 0.0)])
        mol = _FakeTickMol(conf)
        ff = _MovingFF(conf, moves=[(5.0, 0.0, 0.0)])
        self_ = _wire_helpers(_tick_self(mol, ff, steps=1, max_move=0.0))
        fn = _tick_fn()
        fn(self_)

        assert ff.minimize_calls == [1]
        assert conf.set_calls == []
        assert conf.positions[0].x == pytest.approx(5.0)


# ---------------------------------------------------------------------------
# _start
# ---------------------------------------------------------------------------


class _FakeStartConf:
    def __init__(self, n):
        self.n = n

    def GetAtomPosition(self, i):
        return f"p{i}"


class _FakeStartMol:
    def __init__(self, n=3, n_confs=1):
        self.n = n
        self.n_confs = n_confs

    def GetNumConformers(self):
        return self.n_confs

    def GetNumAtoms(self):
        return self.n

    def GetConformer(self):
        return _FakeStartConf(self.n)


_UNSET = object()


class _FakeStartAllChem:
    def __init__(self, mmff_props=_UNSET, ff=None):
        self.mmff_props = mmff_props if mmff_props is not _UNSET else object()
        self.ff = ff if ff is not None else _FakeFF()
        self.mmff_calls = []
        self.uff_calls = []

    def MMFFGetMoleculeProperties(self, mol, mmffVariant=None):
        self.mmff_calls.append(mmffVariant)
        return self.mmff_props

    def MMFFGetMoleculeForceField(self, mol, props):
        return self.ff

    def UFFGetMoleculeForceField(self, mol):
        self.uff_calls.append(mol)
        return self.ff


def _start_fn(allchem, qmessage=None):
    return extract_function(
        STEP_OPT_PATH,
        "StepOptimizerDialog",
        "_start",
        {
            "AllChem": allchem,
            "QMessageBox": qmessage if qmessage is not None else MagicMock(),
            "PLUGIN_NAME": "Step Optimizer",
            "logging": MagicMock(),
        },
    )


def _start_self(mol, ff_text="MMFF94"):
    return SimpleNamespace(
        context=SimpleNamespace(current_mol=mol),
        combo_ff=MagicMock(currentText=lambda: ff_text, setEnabled=MagicMock()),
        btn_start=MagicMock(),
        btn_stop=MagicMock(),
        lbl_state=MagicMock(),
        timer=MagicMock(),
        target_mol=None,
        ff=None,
        step_count=0,
        steps_since_checkpoint=0,
        original_coords=[],
        running=False,
    )


class TestStepOptimizerStart:
    def test_no_molecule_warns(self):
        qmsg = MagicMock()
        allchem = _FakeStartAllChem()
        fn = _start_fn(allchem, qmessage=qmsg)
        self_ = _start_self(None)
        fn(self_)
        qmsg.warning.assert_called_once()
        self_.timer.start.assert_not_called()

    def test_no_conformer_warns(self):
        qmsg = MagicMock()
        allchem = _FakeStartAllChem()
        fn = _start_fn(allchem, qmessage=qmsg)
        mol = _FakeStartMol(n_confs=0)
        self_ = _start_self(mol)
        fn(self_)
        qmsg.warning.assert_called_once()
        self_.timer.start.assert_not_called()

    def test_mmff_props_none_warns_and_no_timer_start(self):
        qmsg = MagicMock()
        allchem = _FakeStartAllChem(mmff_props=None)
        fn = _start_fn(allchem, qmessage=qmsg)
        mol = _FakeStartMol()
        self_ = _start_self(mol, ff_text="MMFF94")
        fn(self_)
        qmsg.warning.assert_called_once()
        self_.timer.start.assert_not_called()
        assert self_.running is False

    def test_uff_path_builds_uff_force_field_and_initializes(self):
        ff = _FakeFF()
        allchem = _FakeStartAllChem(ff=ff)
        fn = _start_fn(allchem)
        mol = _FakeStartMol()
        self_ = _start_self(mol, ff_text="UFF")
        fn(self_)
        assert allchem.uff_calls == [mol]
        assert allchem.mmff_calls == []
        assert ff.initialize_called is True
        self_.timer.start.assert_called_once()
        assert self_.running is True
        assert self_.target_mol is mol

    def test_mmff_success_path_starts_timer(self):
        ff = _FakeFF()
        allchem = _FakeStartAllChem(ff=ff)
        fn = _start_fn(allchem)
        mol = _FakeStartMol()
        self_ = _start_self(mol, ff_text="MMFF94")
        fn(self_)
        assert allchem.mmff_calls == ["MMFF94"]
        assert ff.initialize_called is True
        self_.timer.start.assert_called_once()

    def test_original_coords_captured(self):
        ff = _FakeFF()
        allchem = _FakeStartAllChem(ff=ff)
        fn = _start_fn(allchem)
        mol = _FakeStartMol(n=2)
        self_ = _start_self(mol)
        fn(self_)
        assert self_.original_coords == ["p0", "p1"]


# ---------------------------------------------------------------------------
# _stop
# ---------------------------------------------------------------------------


def _stop_fn():
    return extract_function(STEP_OPT_PATH, "StepOptimizerDialog", "_stop", {})


class TestStepOptimizerStop:
    def test_stop_pushes_undo_checkpoint(self):
        fn = _stop_fn()
        ctx = SimpleNamespace(push_undo_checkpoint=MagicMock())
        self_ = SimpleNamespace(
            timer=MagicMock(),
            running=True,
            lbl_state=MagicMock(),
            context=ctx,
            target_mol=MagicMock(),
            steps_since_checkpoint=5,
            btn_start=MagicMock(),
            btn_stop=MagicMock(),
            combo_ff=MagicMock(),
        )
        # wire helper methods as bound-like closures
        def pause_timer():
            self_.timer.stop()
            self_.running = False

        def reenable():
            self_.btn_start.setEnabled(True)
            self_.btn_stop.setEnabled(False)
            self_.combo_ff.setEnabled(True)

        self_._pause_timer = pause_timer
        self_._reenable_controls = reenable
        fn(self_)
        self_.timer.stop.assert_called_once()
        ctx.push_undo_checkpoint.assert_called_once()
        assert self_.steps_since_checkpoint == 0
        self_.btn_start.setEnabled.assert_called_with(True)
        self_.btn_stop.setEnabled.assert_called_with(False)

    def test_stop_without_target_mol_does_not_push_checkpoint(self):
        fn = _stop_fn()
        ctx = SimpleNamespace(push_undo_checkpoint=MagicMock())
        self_ = SimpleNamespace(
            timer=MagicMock(),
            running=True,
            lbl_state=MagicMock(),
            context=ctx,
            target_mol=None,
            steps_since_checkpoint=0,
            btn_start=MagicMock(),
            btn_stop=MagicMock(),
            combo_ff=MagicMock(),
        )

        def pause_timer():
            self_.timer.stop()
            self_.running = False

        def reenable():
            pass

        self_._pause_timer = pause_timer
        self_._reenable_controls = reenable
        fn(self_)
        ctx.push_undo_checkpoint.assert_not_called()


# ---------------------------------------------------------------------------
# accept / reject / closeEvent
# ---------------------------------------------------------------------------


_FAKE_SUPER_GLOBALS = {
    "super": lambda *a: SimpleNamespace(accept=lambda: None, reject=lambda: None)
}


class _RSMol:
    def __init__(self, n, positions):
        self.n = n
        self.positions = list(positions)
        self.set_calls = []

    def GetNumAtoms(self):
        return self.n

    def GetConformer(self):
        return self

    def GetAtomPosition(self, i):
        return self.positions[i]

    def SetAtomPosition(self, i, pos):
        self.set_calls.append((i, pos))


class TestStepOptimizerAcceptReject:
    def test_accept_pushes_checkpoint_when_steps_ran(self):
        fn = extract_function(
            STEP_OPT_PATH, "StepOptimizerDialog", "accept", _FAKE_SUPER_GLOBALS
        )
        ctx = make_context()
        self_ = SimpleNamespace(
            context=ctx,
            target_mol=MagicMock(),
            steps_since_checkpoint=3,
            timer=MagicMock(),
            running=True,
        )

        def pause_timer():
            self_.timer.stop()
            self_.running = False

        self_._pause_timer = pause_timer
        fn(self_)
        ctx.push_undo_checkpoint.assert_called_once()
        ctx.register_window.assert_called_with("main_panel", None)
        self_.timer.stop.assert_called_once()

    def test_accept_no_checkpoint_when_no_steps_ran(self):
        fn = extract_function(
            STEP_OPT_PATH, "StepOptimizerDialog", "accept", _FAKE_SUPER_GLOBALS
        )
        ctx = make_context()
        self_ = SimpleNamespace(
            context=ctx,
            target_mol=MagicMock(),
            steps_since_checkpoint=0,
            timer=MagicMock(),
            running=False,
        )
        self_._pause_timer = lambda: None
        fn(self_)
        ctx.push_undo_checkpoint.assert_not_called()
        ctx.register_window.assert_called_with("main_panel", None)

    def test_reject_restores_original_coordinates(self):
        fn = extract_function(
            STEP_OPT_PATH, "StepOptimizerDialog", "reject", _FAKE_SUPER_GLOBALS
        )
        ctx = make_context()
        mol = _RSMol(2, ["p0_new", "p1_new"])
        self_ = SimpleNamespace(
            context=ctx,
            target_mol=mol,
            original_coords=["p0", "p1"],
            timer=MagicMock(),
            running=True,
        )
        self_._pause_timer = lambda: (self_.timer.stop(), setattr(self_, "running", False))
        fn(self_)
        assert mol.set_calls == [(0, "p0"), (1, "p1")]
        ctx.refresh_3d_view.assert_called_once()
        ctx.register_window.assert_called_with("main_panel", None)

    def test_reject_without_run_does_not_raise(self):
        fn = extract_function(
            STEP_OPT_PATH, "StepOptimizerDialog", "reject", _FAKE_SUPER_GLOBALS
        )
        ctx = make_context()
        self_ = SimpleNamespace(
            context=ctx,
            target_mol=None,
            original_coords=[],
            timer=MagicMock(),
            running=False,
        )
        self_._pause_timer = lambda: None
        fn(self_)  # must not raise
        ctx.register_window.assert_called_with("main_panel", None)

    def test_close_event_delegates_to_accept_and_ignores_event(self):
        fn = extract_function(STEP_OPT_PATH, "StepOptimizerDialog", "closeEvent")
        calls = []
        self_ = SimpleNamespace(accept=lambda: calls.append("accept"))
        event = MagicMock()
        fn(self_, event)
        assert calls == ["accept"]
        event.ignore.assert_called_once()
