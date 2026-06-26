"""
Headless GUI tests for analysis/editor plugin dialogs.

Plugins covered (all registry-visible):
  - XYZ Editor          → XYZEditorWindow
  - Symmetry Analyzer   → SymmetryAnalysisPlugin
  - Molecule Comparator → MoleculeComparator
  - Conformational Search → ConformerSearchDialog

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_analysis.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

XYZ_EDITOR_PATH = PLUGINS_DIR / "XYZ_Editor" / "xyz_editor.py"
SYMMETRY_PATH = PLUGINS_DIR / "Symmetry_Analyzer" / "symmetry_analyzer.py"
COMPARATOR_PATH = PLUGINS_DIR / "Molecule_Comparator" / "molecule_comparator.py"
CONF_SEARCH_PATH = PLUGINS_DIR / "Conformational_Search" / "conf_search.py"

# ---------------------------------------------------------------------------
# Load plugins once at collection time — real Qt, chemistry mocked.
# ---------------------------------------------------------------------------

with mock_chemistry_imports():
    _xyz_editor = load_plugin_for_gui(XYZ_EDITOR_PATH)
    _symmetry = load_plugin_for_gui(SYMMETRY_PATH)
    _comparator = load_plugin_for_gui(COMPARATOR_PATH)
    _conf_search = load_plugin_for_gui(CONF_SEARCH_PATH)


# ---------------------------------------------------------------------------
# Context helpers
# ---------------------------------------------------------------------------


def _ctx_no_mol() -> MagicMock:
    """Context with no main window and no active molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# XYZEditorWindow  (visible plugin: "XYZ Editor")
# ===========================================================================


class TestXYZEditorWindow:
    """XYZEditorWindow with no main window and no molecule."""

    @pytest.fixture
    def win(self, qapp):
        ctx = _ctx_no_mol()
        w = _xyz_editor.XYZEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "XYZ Editor"

    def test_table_has_five_columns(self, win):
        assert win.table.columnCount() == 5

    def test_table_column_headers(self, win):
        headers = [
            win.table.horizontalHeaderItem(i).text() for i in range(5)
        ]
        assert headers == ["Index", "Symbol", "X", "Y", "Z"]

    def test_table_empty_with_no_molecule(self, win):
        assert win.table.rowCount() == 0

    def test_apply_button_exists(self, win):
        assert win.apply_btn.text() == "Apply to View"

    def test_save_button_exists(self, win):
        assert win.save_btn.text() == "Save as XYZ..."

    def test_add_atom_button_exists(self, win):
        assert win.add_btn.text() == "Add Atom"

    def test_timer_is_active(self, win):
        assert win.update_timer.isActive()

    def test_resize_to_expected_size(self, win):
        assert win.width() >= 400
        assert win.height() >= 300


# ===========================================================================
# SymmetryAnalysisPlugin  (visible plugin: "Symmetry Analyzer")
# ===========================================================================


class TestSymmetryAnalysisPlugin:
    """SymmetryAnalysisPlugin — pymatgen is mocked so HAS_PYMATGEN=True."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        d = _symmetry.SymmetryAnalysisPlugin(context=ctx)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert "Symmetry" in dlg.windowTitle() or dlg.windowTitle() == ""

    def test_groups_list_is_empty_initially(self, dlg):
        assert dlg.groups_list.count() == 0

    def test_ops_list_is_empty_initially(self, dlg):
        assert dlg.ops_list.count() == 0

    def test_point_group_label_default(self, dlg):
        assert dlg.selected_group_label.text() == "Point Group: -"

    def test_op_details_is_readonly(self, dlg):
        assert dlg.op_details.isReadOnly()

    def test_max_tol_spin_default(self, dlg):
        assert dlg.max_tol_spin.value() == pytest.approx(1.0)

    def test_symmetrize_button_initially_disabled(self, dlg):
        assert not dlg.sym_btn.isEnabled()

    def test_analyze_button_exists(self, dlg):
        assert dlg.calc_btn is not None


# ===========================================================================
# MoleculeComparator  (visible plugin: "Molecule Comparator")
# ===========================================================================


class TestMoleculeComparator:
    """MoleculeComparator with no main window."""

    @pytest.fixture
    def comp(self, qapp):
        ctx = _ctx_no_mol()
        w = _comparator.MoleculeComparator(context=ctx)
        yield w
        w.destroy()

    def test_creates_without_error(self, comp):
        assert comp is not None

    def test_window_title(self, comp):
        assert comp.windowTitle() == "Molecule Comparator"

    def test_mol_list_initially_empty(self, comp):
        assert comp.mol_list.count() == 0

    def test_molecules_list_initially_empty(self, comp):
        assert comp.molecules == []

    def test_style_combo_default_sticks(self, comp):
        assert comp.combo_style.currentText() == "Sticks"

    def test_style_combo_has_four_options(self, comp):
        assert comp.combo_style.count() == 4

    def test_align_method_default(self, comp):
        assert comp.combo_align_method.currentText() == "Substructure (MCS)"

    def test_ignore_hydrogens_unchecked_by_default(self, comp):
        assert not comp.check_ignore_hs.isChecked()

    def test_wireframe_lighting_unchecked_by_default(self, comp):
        assert not comp.check_wireframe_lighting.isChecked()

    def test_add_current_button_exists(self, comp):
        assert comp.btn_add_current.text() == "Add Current"

    def test_align_button_exists(self, comp):
        assert comp.btn_align.text() == "Align & Calculate RMSD"


# ===========================================================================
# ConformerSearchDialog  (visible plugin: "Conformational Search")
# ===========================================================================


class TestConformerSearchDialog:
    """ConformerSearchDialog with no molecule loaded."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        d = _conf_search.ConformerSearchDialog(context=ctx, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Conformational Search & Preview"

    def test_table_has_two_columns(self, dlg):
        assert dlg.table.columnCount() == 2

    def test_table_column_headers(self, dlg):
        assert dlg.table.horizontalHeaderItem(0).text() == "Rank"
        assert dlg.table.horizontalHeaderItem(1).text() == "Energy (kcal/mol)"

    def test_table_initially_empty(self, dlg):
        assert dlg.table.rowCount() == 0

    def test_force_field_combo_has_options(self, dlg):
        texts = [dlg.combo_ff.itemText(i) for i in range(dlg.combo_ff.count())]
        assert "MMFF94" in texts
        assert "UFF" in texts

    def test_show_all_checkbox_unchecked_by_default(self, dlg):
        assert not dlg.cb_show_all.isChecked()

    def test_run_button_exists(self, dlg):
        assert dlg.btn_run.text() == "Run Search"

    def test_close_button_exists(self, dlg):
        assert dlg.btn_close.text() == "Close"

    def test_conformer_data_initially_empty(self, dlg):
        assert dlg.conformer_data == []
