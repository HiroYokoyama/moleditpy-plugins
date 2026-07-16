"""
Headless GUI tests for the Step Optimizer plugin.

Covers: StepOptimizerDialog.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_step_optimizer.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

STEP_OPT_PATH = PLUGINS_DIR / "Step_Optimizer" / "step_optimizer.py"

with mock_chemistry_imports():
    _step = load_plugin_for_gui(STEP_OPT_PATH)


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
    monkeypatch.setattr(_step, "QMessageBox", box)
    return {_step.__name__: box}


# ===========================================================================
# StepOptimizerDialog  (visible plugin: "Step Optimizer")
# ===========================================================================


class TestStepOptimizerDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _step.StepOptimizerDialog(context=_ctx(), parent=None)
        yield d
        d.timer.stop()
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Step Optimizer"

    def test_force_field_options(self, dlg):
        texts = [dlg.combo_ff.itemText(i) for i in range(dlg.combo_ff.count())]
        assert texts == ["MMFF94", "MMFF94s", "UFF"]

    def test_default_ff_from_mmff_setting(self, dlg):
        assert dlg.combo_ff.currentText() == "MMFF94"

    def test_default_ff_from_uff_setting(self, qapp):
        d = _step.StepOptimizerDialog(context=_ctx(setting="UFF_RDKIT"), parent=None)
        assert d.combo_ff.currentText() == "UFF"
        d.timer.stop()
        d.destroy()

    def test_steps_per_frame_defaults(self, dlg):
        assert dlg.spin_steps.value() == 1
        assert (dlg.spin_steps.minimum(), dlg.spin_steps.maximum()) == (1, 200)

    def test_max_move_defaults(self, dlg):
        assert dlg.spin_max_move.value() == pytest.approx(0.5)
        assert dlg.spin_max_move.specialValueText() == "Off"

    def test_initial_labels(self, dlg):
        assert dlg.lbl_step.text() == "Step: 0"
        assert dlg.lbl_energy.text() == "Energy: --"
        assert dlg.lbl_state.text() == "State: Idle"

    def test_initial_button_states(self, dlg):
        assert dlg.btn_start.isEnabled()
        assert not dlg.btn_stop.isEnabled()

    def test_timer_configured_but_idle(self, dlg):
        assert dlg.timer.interval() == 50
        assert not dlg.timer.isActive()
        assert not dlg.running

    def test_start_without_molecule_warns(self, dlg, no_msgbox):
        dlg._start()
        no_msgbox[_step.__name__].warning.assert_called_once()
        assert not dlg.running
        assert not dlg.timer.isActive()

    def test_start_without_conformer_warns(self, dlg, no_msgbox):
        mol = MagicMock()
        mol.GetNumConformers.return_value = 0
        dlg.context.current_mol = mol
        dlg._start()
        msg = no_msgbox[_step.__name__].warning.call_args[0][2]
        assert "3D" in msg
        assert not dlg.running

    def test_stop_pushes_undo_checkpoint(self, dlg):
        dlg.target_mol = MagicMock()
        dlg.steps_since_checkpoint = 7
        dlg._stop()
        dlg.context.push_undo_checkpoint.assert_called_once()
        assert dlg.steps_since_checkpoint == 0
        assert dlg.lbl_state.text() == "State: Stopped"
        assert dlg.btn_start.isEnabled()
        assert not dlg.btn_stop.isEnabled()

    def test_stop_without_run_skips_checkpoint(self, dlg):
        dlg._stop()
        dlg.context.push_undo_checkpoint.assert_not_called()

    def test_accept_unregisters_window(self, dlg):
        dlg.accept()
        dlg.context.register_window.assert_called_with("main_panel", None)

    def test_accept_checkpoints_pending_steps(self, dlg):
        dlg.target_mol = MagicMock()
        dlg.steps_since_checkpoint = 3
        dlg.accept()
        dlg.context.push_undo_checkpoint.assert_called_once()

    def test_reject_restores_original_coordinates(self, dlg):
        dlg.target_mol = MagicMock()
        dlg.original_coords = [MagicMock(), MagicMock()]
        dlg.reject()
        conf = dlg.target_mol.GetConformer.return_value
        assert conf.SetAtomPosition.call_count == 2
        dlg.context.refresh_3d_view.assert_called_once()
        dlg.context.register_window.assert_called_with("main_panel", None)

    def test_run_plugin_without_molecule_warns(self, no_msgbox):
        ctx = _ctx()
        _step.run_plugin(ctx)
        no_msgbox[_step.__name__].warning.assert_called_once()
        ctx.register_window.assert_not_called()

    def test_run_plugin_raises_existing_window(self, no_msgbox):
        ctx = _ctx()
        ctx.current_mol = MagicMock()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _step.run_plugin(ctx)
        existing.show.assert_called_once()
        existing.raise_.assert_called_once()

    def test_initialize_registers_menu_action(self):
        ctx = _ctx()
        _step.initialize(ctx)
        assert ctx.add_menu_action.call_args[0][0] == "3D Edit/Step Optimizer..."
