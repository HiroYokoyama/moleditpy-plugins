"""
GUI tests for the OpenBabel Conversion Tool plugin.

Covers: MoleculeSelectionDialog — openbabel/pybel are mocked by
mock_chemistry_imports(), so OBABEL_AVAILABLE=True and pybel is a MagicMock.
The dialog only needs a list of mock "molecule" objects with .title / .formula.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

OPENBABEL_PATH = PLUGINS_DIR / "OpenBabel_Conversion_Tool" / "openbabel_conversion_tool.py"

with mock_chemistry_imports():
    _openbabel = load_plugin_for_gui(OPENBABEL_PATH)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_mock_mols(n: int = 3) -> list:
    mols = []
    for i in range(n):
        m = MagicMock()
        m.title = f"Molecule {i + 1}"
        m.formula = f"C{i + 1}H{2 * (i + 1)}"
        mols.append(m)
    return mols


# ===========================================================================
# MoleculeSelectionDialog  (OpenBabel Conversion Tool)
# ===========================================================================


class TestMoleculeSelectionDialog:
    """MoleculeSelectionDialog with mock pybel molecule objects."""

    @pytest.fixture
    def dlg(self, qapp):
        mols = _make_mock_mols(3)
        d = _openbabel.MoleculeSelectionDialog(mols, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Select Molecule"

    def test_list_widget_has_correct_item_count(self, dlg):
        assert dlg.list_widget.count() == 3

    def test_first_row_selected_by_default(self, dlg):
        assert dlg.list_widget.currentRow() == 0

    def test_get_selected_index_returns_zero_initially(self, dlg):
        assert dlg.get_selected_index() == 0

    def test_get_selected_index_after_selection_change(self, dlg):
        dlg.list_widget.setCurrentRow(2)
        assert dlg.get_selected_index() == 2

    def test_single_selection_mode(self, dlg):
        from PyQt6.QtWidgets import QAbstractItemView
        assert dlg.list_widget.selectionMode() == QAbstractItemView.SelectionMode.SingleSelection

    def test_empty_mol_list_creates_empty_widget(self, qapp):
        d = _openbabel.MoleculeSelectionDialog([], parent=None)
        assert d.list_widget.count() == 0
        d.destroy()

    def test_get_selected_index_returns_none_when_nothing_selected(self, qapp):
        d = _openbabel.MoleculeSelectionDialog([], parent=None)
        d.list_widget.clearSelection()
        assert d.get_selected_index() is None
        d.destroy()

    def test_single_mol_list(self, qapp):
        d = _openbabel.MoleculeSelectionDialog(_make_mock_mols(1), parent=None)
        assert d.list_widget.count() == 1
        d.destroy()

    def test_obabel_available_flag(self):
        # Under mock_chemistry_imports, openbabel is mocked → OBABEL_AVAILABLE=True
        assert _openbabel.OBABEL_AVAILABLE is True
