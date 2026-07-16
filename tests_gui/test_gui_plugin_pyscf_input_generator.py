"""
Headless GUI tests for the PySCF Input Generator plugin.

Covers: PyscfSetupDialog.

Follows the same pattern as MOPAC/GAMESS/Psi4:
- parent=None, mol=None is safe (calc_initial_charge_mult guards with `if not self.mol`)
- setup_ui() is pure Qt widget construction
- load_presets_from_file() reads a settings file only if it exists

Run locally:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/test_gui_plugin_pyscf_input_generator.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
PYSCF_PATH = PLUGINS_DIR / "PySCF_Input_Generator" / "pyscf_input_generator.py"

with mock_chemistry_imports():
    _pyscf = load_plugin_for_gui(PYSCF_PATH)


# ===========================================================================
# PyscfSetupDialog  (visible plugin: "PySCF Input Generator")
# ===========================================================================


class TestPyscfSetupDialog:
    """PyscfSetupDialog with mol=None — exercises the no-molecule UI path."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _pyscf.PyscfSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "PySCF Input Generator"

    def test_category_combo_default(self, dlg):
        assert dlg.category_combo.currentText() == "Hartree-Fock"

    def test_basis_combo_default(self, dlg):
        assert dlg.basis_combo.currentText() == "def2-svp"

    def test_charge_spinbox_range(self, dlg):
        assert dlg.charge_spin.minimum() == -10
        assert dlg.charge_spin.maximum() == 10

    def test_mult_spinbox_minimum_is_one(self, dlg):
        assert dlg.mult_spin.minimum() == 1

    def test_symmetry_checked_by_default(self, dlg):
        assert dlg.symmetry_check.isChecked()

    def test_functional_combo_initially_disabled(self, dlg):
        assert not dlg.functional_combo.isEnabled()

    def test_post_hf_combo_initially_disabled(self, dlg):
        assert not dlg.post_hf_combo.isEnabled()

    def test_post_hf_ref_combo_initially_disabled(self, dlg):
        assert not dlg.post_hf_ref_combo.isEnabled()

    def test_preset_combo_initially_empty(self, dlg):
        assert dlg.preset_combo.count() == 0

    def test_preview_area_is_editable(self, dlg):
        assert not dlg.preview_text.isReadOnly()

    def test_basis_combo_includes_sto_3g(self, dlg):
        items = [dlg.basis_combo.itemText(i) for i in range(dlg.basis_combo.count())]
        assert "sto-3g" in items

    def test_category_combo_includes_dft(self, dlg):
        items = [dlg.category_combo.itemText(i) for i in range(dlg.category_combo.count())]
        assert "DFT" in items

    def test_save_button_label(self, dlg):
        assert dlg.btn_save.text() == "Save Python Script..."
