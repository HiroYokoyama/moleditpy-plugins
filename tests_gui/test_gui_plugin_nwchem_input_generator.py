"""
Headless GUI tests for the NWChem Input Generator plugin.

Covers: NwchemSetupDialog.

Follows the same pattern as MOPAC/GAMESS/Psi4:
- parent=None, mol=None is safe (calc_initial_charge_mult guards with `if not self.mol`)
- setup_ui() is pure Qt widget construction
- load_presets_from_file() reads a settings file only if it exists

Run locally:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/test_gui_plugin_nwchem_input_generator.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
NWCHEM_PATH = PLUGINS_DIR / "NWChem_Input_Generator" / "nwchem_input_generator.py"

with mock_chemistry_imports():
    _nwchem = load_plugin_for_gui(NWCHEM_PATH)


# ===========================================================================
# NwchemSetupDialog  (visible plugin: "NWChem Input Generator")
# ===========================================================================


class TestNwchemSetupDialog:
    """NwchemSetupDialog with mol=None — exercises the no-molecule UI path."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _nwchem.NwchemSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "NWChem Input Generator"

    def test_default_title_field(self, dlg):
        assert dlg.title_edit.text() == "NWChem Job"

    def test_module_combo_default(self, dlg):
        assert dlg.module_combo.currentText() == "dft"

    def test_functional_combo_default(self, dlg):
        assert dlg.functional_combo.currentText() == "b3lyp"

    def test_task_combo_default(self, dlg):
        assert dlg.task_combo.currentText() == "optimize"

    def test_basis_combo_default(self, dlg):
        assert dlg.basis_combo.currentText() == "6-31G*"

    def test_charge_spinbox_range(self, dlg):
        assert dlg.charge_spin.minimum() == -10
        assert dlg.charge_spin.maximum() == 10

    def test_mult_spinbox_minimum_is_one(self, dlg):
        assert dlg.mult_spin.minimum() == 1

    def test_preset_combo_initially_empty(self, dlg):
        assert dlg.preset_combo.count() == 0

    def test_preview_area_is_editable(self, dlg):
        assert not dlg.preview_text.isReadOnly()

    def test_module_combo_items_include_scf(self, dlg):
        items = [dlg.module_combo.itemText(i) for i in range(dlg.module_combo.count())]
        assert "scf" in items

    def test_basis_combo_includes_cc_pvdz(self, dlg):
        items = [dlg.basis_combo.itemText(i) for i in range(dlg.basis_combo.count())]
        assert "cc-pvdz" in items

    def test_task_combo_includes_freq(self, dlg):
        items = [dlg.task_combo.itemText(i) for i in range(dlg.task_combo.count())]
        assert "freq" in items
