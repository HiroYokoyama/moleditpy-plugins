"""
Headless GUI tests for the Molecule Comparator plugin.

Covers: MoleculeComparator.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_molecule_comparator.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

COMPARATOR_PATH = PLUGINS_DIR / "Molecule_Comparator" / "molecule_comparator.py"

with mock_chemistry_imports():
    _comparator = load_plugin_for_gui(COMPARATOR_PATH)


def _ctx_no_mol() -> MagicMock:
    """Context with no main window and no active molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


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
# Helpers: comparator with visualization side effects stubbed out
# ===========================================================================


def _make_comp(qapp_unused=None):
    ctx = _ctx_no_mol()
    w = _comparator.MoleculeComparator(context=ctx)
    w.update_visualization = MagicMock()
    w.update_wireframe_lighting = MagicMock()
    w.reset_view = MagicMock()
    return w, ctx


def _fake_entry_mol(name=None):
    mol = MagicMock()
    mol.HasProp.return_value = name is not None
    mol.GetProp.return_value = name
    return mol


def _add_mols(w, n):
    for i in range(n):
        w.molecules.append(
            {
                "name": f"M{i}",
                "mol": MagicMock(),
                "color": "#ff0000",
                "scope": "Carbon Only",
                "rms": None,
            }
        )
    w.update_list()


# ===========================================================================
# Molecule list management
# ===========================================================================


class TestComparatorListManagement:
    def test_add_current_without_molecule_warns(self, qapp):
        w, ctx = _make_comp()
        with pytest.MonkeyPatch.context() as mp:
            warned = MagicMock()
            mp.setattr(_comparator.QMessageBox, "warning", warned)
            w.add_current_molecule()
            warned.assert_called_once()
        assert w.molecules == []
        w.destroy()

    def test_add_current_uses_mol_name_prop(self, qapp):
        w, ctx = _make_comp()
        ctx.current_molecule = _fake_entry_mol("benzene")
        w.add_current_molecule()
        assert w.molecules[0]["name"] == "benzene"
        w.destroy()

    def test_add_current_falls_back_to_numbered_name(self, qapp):
        w, ctx = _make_comp()
        ctx.current_molecule = _fake_entry_mol(None)
        w.add_current_molecule()
        assert w.molecules[0]["name"] == "Mol 1"
        w.destroy()

    def test_add_current_cycles_default_colors(self, qapp):
        w, ctx = _make_comp()
        ctx.current_molecule = _fake_entry_mol(None)
        w.add_current_molecule()
        w.add_current_molecule()
        colors = _comparator.DEFAULT_COLORS
        assert w.molecules[0]["color"] == colors[0]
        assert w.molecules[1]["color"] == colors[1]
        w.destroy()

    def test_update_list_marks_reference_and_selects_last(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 3)
        assert w.mol_list.item(0).text() == "M0 (Ref)"
        assert w.mol_list.item(1).text() == "M1"
        assert w.mol_list.currentRow() == 2
        w.destroy()

    def test_set_as_reference_moves_selected_to_top(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 3)
        w.mol_list.setCurrentRow(2)
        w.set_as_reference()
        assert [m["name"] for m in w.molecules] == ["M2", "M0", "M1"]
        assert w.mol_list.item(0).text() == "M2 (Ref)"
        w.destroy()

    def test_set_as_reference_noop_for_reference_row(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 2)
        w.mol_list.setCurrentRow(0)
        w.set_as_reference()
        assert [m["name"] for m in w.molecules] == ["M0", "M1"]
        w.destroy()

    def test_remove_molecule_removes_current_row(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 2)
        w.mol_list.setCurrentRow(0)
        w.remove_molecule()
        assert [m["name"] for m in w.molecules] == ["M1"]
        w.destroy()

    def test_remove_molecule_noop_without_selection(self, qapp):
        w, ctx = _make_comp()
        w.remove_molecule()
        assert w.molecules == []
        w.destroy()


# ===========================================================================
# Selection-driven UI state
# ===========================================================================


class TestComparatorSelectionState:
    def test_invalid_row_disables_color_and_scope(self, qapp):
        w, ctx = _make_comp()
        w.on_selection_changed(-1)
        assert not w.btn_color.isEnabled()
        assert not w.combo_scope.isEnabled()
        w.destroy()

    def test_valid_row_enables_and_reflects_scope(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 1)
        w.molecules[0]["scope"] = "All Atoms"
        w.on_selection_changed(0)
        assert w.btn_color.isEnabled()
        assert w.combo_scope.isEnabled()
        assert w.combo_scope.currentText() == "All Atoms"
        w.destroy()

    def test_change_scope_updates_entry(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 1)
        w.mol_list.setCurrentRow(0)
        w.combo_scope.setCurrentText("All Atoms")
        w.change_scope()
        assert w.molecules[0]["scope"] == "All Atoms"
        w.destroy()

    def test_text_color_black_on_light_white_on_dark(self, qapp):
        w, ctx = _make_comp()
        assert w._get_text_color_for_background("#ffffff") == "black"
        assert w._get_text_color_for_background("#000000") == "white"
        assert w._get_text_color_for_background("garbage") == "black"
        w.destroy()


# ===========================================================================
# Results table / clipboard
# ===========================================================================


class TestComparatorResults:
    def test_results_table_formats_rms_states(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 3)
        w.molecules[0]["rms"] = 0.0
        w.molecules[1]["rms"] = -1.0
        w.molecules[2]["rms"] = None
        w.update_results_table()
        assert w.table_results.item(0, 1).text() == "0.0000"
        assert w.table_results.item(1, 1).text() == "N/A"
        assert w.table_results.item(2, 1).text() == "-"
        w.destroy()

    def test_copy_results_noop_when_empty(self, qapp):
        from PyQt6.QtWidgets import QApplication

        w, ctx = _make_comp()
        QApplication.clipboard().setText("sentinel")
        w.copy_results_to_clipboard()
        assert QApplication.clipboard().text() == "sentinel"
        w.destroy()

    def test_copy_results_writes_tsv(self, qapp):
        from PyQt6.QtWidgets import QApplication

        w, ctx = _make_comp()
        _add_mols(w, 2)
        w.molecules[0]["rms"] = 0.0
        w.molecules[1]["rms"] = 1.23456
        w.copy_results_to_clipboard()
        text = QApplication.clipboard().text()
        assert text.splitlines() == [
            "Molecule\tRMSD (Å)",
            "M0\t0.0000",
            "M1\t1.2346",
        ]
        ctx.show_status_message.assert_called_with(
            "Results copied to clipboard", 3000
        )
        w.destroy()

    def test_run_alignment_needs_two_molecules(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 1)
        w.run_alignment()
        assert not hasattr(w, "progress_dialog")
        w.destroy()


# ===========================================================================
# Close / cleanup behavior
# ===========================================================================


class TestComparatorCleanup:
    def _closable(self):
        w, ctx = _make_comp()
        w.mw = MagicMock()  # cleanup path touches mw.view_3d_manager
        w.exit_3d_only_mode = MagicMock()
        ctx.draw_molecule_3d.reset_mock()  # construction already draws once
        return w, ctx

    def test_cleanup_redraws_current_molecule(self, qapp):
        w, ctx = self._closable()
        ctx.current_molecule = MagicMock()
        w.cleanup_and_close()
        ctx.draw_molecule_3d.assert_called_once_with(ctx.current_molecule)
        w.exit_3d_only_mode.assert_called_once()
        w.destroy()

    def test_cleanup_clears_plotter_without_molecule(self, qapp):
        w, ctx = self._closable()
        ctx.current_molecule = None
        w.cleanup_and_close()
        ctx.plotter.clear.assert_called_once()
        w.destroy()

    def test_cleanup_resets_color_overrides(self, qapp):
        w, ctx = self._closable()
        w.mw.view_3d_manager._plugin_color_overrides = {"stale": 1}
        w.cleanup_and_close()
        assert w.mw.view_3d_manager._plugin_color_overrides == {}
        w.destroy()

    def test_close_plugin_hides_window(self, qapp):
        w, ctx = self._closable()
        w.show()
        w.close_plugin()
        assert w.isHidden()
        w.destroy()


# ===========================================================================
# initialize() — document reset closes the window
# ===========================================================================


class TestComparatorInitialize:
    def test_registers_document_reset_handler(self, qapp):
        ctx = _ctx_no_mol()
        _comparator.initialize(ctx)
        ctx.register_document_reset_handler.assert_called_once()

    def test_reset_closes_open_comparator_window(self, qapp):
        ctx = _ctx_no_mol()
        handlers = []
        ctx.register_document_reset_handler.side_effect = handlers.append
        _comparator.initialize(ctx)
        mw = MagicMock()
        ctx.get_main_window.return_value = mw
        handlers[0]()
        mw.molecule_comparator_window.close.assert_called_once()
