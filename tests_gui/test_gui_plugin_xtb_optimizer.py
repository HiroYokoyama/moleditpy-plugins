"""
Headless GUI tests for the xTB Optimizer plugin.

Covers: XtbOptimizerDialog, XtbWorker.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_xtb_optimizer.py
"""

from __future__ import annotations

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
