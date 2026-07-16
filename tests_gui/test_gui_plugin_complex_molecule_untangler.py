"""
Headless GUI tests for the Complex Molecule Untangler plugin.

Covers: UntanglerDialog, run_plugin.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_complex_molecule_untangler.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

UNTANGLER_PATH = (
    PLUGINS_DIR / "Complex_Molecule_Untangler" / "complex_molecule_untangler.py"
)

with mock_chemistry_imports():
    _untangler = load_plugin_for_gui(UNTANGLER_PATH)


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
    monkeypatch.setattr(_untangler, "QMessageBox", box)
    return {_untangler.__name__: box}


# ===========================================================================
# UntanglerDialog  (visible plugin: "Complex Molecule Untangler")
# ===========================================================================


class TestUntanglerDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _untangler.UntanglerDialog(context=_ctx(), parent=None)
        yield d
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Complex Molecule Untangler"

    def test_force_field_options(self, dlg):
        texts = [dlg.combo_ff.itemText(i) for i in range(dlg.combo_ff.count())]
        assert texts == ["MMFF94", "UFF"]

    def test_default_ff_from_uff_setting(self, qapp):
        d = _untangler.UntanglerDialog(context=_ctx(setting="UFF"), parent=None)
        assert d.combo_ff.currentText() == "UFF"
        d.destroy()

    def test_iterations_defaults(self, dlg):
        assert dlg.spin_iter.value() == 500
        assert (dlg.spin_iter.minimum(), dlg.spin_iter.maximum()) == (100, 10000)
        assert dlg.spin_iter.singleStep() == 100

    def test_progress_bar_initial(self, dlg):
        assert dlg.pbar.value() == 0
        assert not dlg.pbar.isTextVisible()

    def test_run_button(self, dlg):
        assert dlg.btn_run.text() == "Untangle Molecule"
        assert dlg.btn_run.isEnabled()

    def test_run_without_molecule_warns(self, dlg, no_msgbox):
        dlg.run_untangle()
        no_msgbox[_untangler.__name__].warning.assert_called_once()
        assert dlg.btn_run.isEnabled()
        assert dlg.worker is None

    def test_on_finished_error_restores_button(self, dlg, no_msgbox):
        dlg.btn_run.setEnabled(False)
        dlg.btn_run.setText("Processing...")
        dlg.on_finished(None, "No rotatable bonds found.")
        assert dlg.btn_run.isEnabled()
        assert dlg.btn_run.text() == "Untangle Molecule"
        no_msgbox[_untangler.__name__].warning.assert_called_once()

    def test_on_finished_success_updates_molecule(self, dlg, no_msgbox):
        new_mol = MagicMock()
        dlg.on_finished(new_mol, "Processed 5 bonds.")
        assert dlg.context.current_mol is new_mol
        dlg.context.refresh_3d_view.assert_called_once()
        dlg.context.push_undo_checkpoint.assert_called_once()
        no_msgbox[_untangler.__name__].information.assert_called_once()

    def test_worker_stores_parameters(self, qapp):
        mol = MagicMock()
        w = _untangler.UntangleWorker(mol, max_iter=250, force_field="UFF")
        assert w.mol is mol
        assert w.max_iter == 250
        assert w.force_field == "UFF"

    def test_run_plugin_registers_new_dialog(self):
        ctx = _ctx()
        ctx.get_window.return_value = None
        _untangler.run_plugin(ctx)
        args = ctx.register_window.call_args[0]
        assert args[0] == "main_panel"
        assert isinstance(args[1], _untangler.UntanglerDialog)
        args[1].destroy()

    def test_run_plugin_raises_existing_window(self):
        ctx = _ctx()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _untangler.run_plugin(ctx)
        existing.show.assert_called_once()
        ctx.register_window.assert_not_called()

    def test_initialize_then_run_launches(self):
        ctx = _ctx()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _untangler.initialize(ctx)
        _untangler.run(MagicMock())
        existing.show.assert_called_once()
