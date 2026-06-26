"""
GUI tests for OpenBabel Conversion Tool and PubChem Structure Identifier.

OpenBabel: MoleculeSelectionDialog — openbabel/pybel are mocked by
mock_chemistry_imports(), so OBABEL_AVAILABLE=True and pybel is a MagicMock.
The dialog only needs a list of mock "molecule" objects with .title / .formula.

PubChem: MoleculeDetailsDialog — pure Qt, builds HTML from a details dict.
No network calls in __init__. PubChemResolver static methods are tested
separately (pure Python, no Qt).
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

OPENBABEL_PATH = PLUGINS_DIR / "OpenBabel_Conversion_Tool" / "openbabel_conversion_tool.py"
PUBCHEM_PATH = PLUGINS_DIR / "PubChem_Structure_Identifier" / "pubchem_structure_identifier.py"

with mock_chemistry_imports():
    _openbabel = load_plugin_for_gui(OPENBABEL_PATH)
    _pubchem = load_plugin_for_gui(PUBCHEM_PATH)


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


_FULL_DETAILS = {
    "Common Name": "Aspirin",
    "IUPAC Name": "2-acetoxybenzoic acid",
    "Formula": "C9H8O4",
    "Mol. Weight": "180.16",
    "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
}


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


# ===========================================================================
# MoleculeDetailsDialog  (PubChem Structure Identifier)
# ===========================================================================


class TestMoleculeDetailsDialog:
    """MoleculeDetailsDialog — displays PubChem properties, no network calls."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _pubchem.MoleculeDetailsDialog(_FULL_DETAILS, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Molecule Details - PubChem"

    def test_details_stored(self, dlg):
        assert dlg.details is _FULL_DETAILS

    def test_copy_button_exists(self, dlg):
        assert dlg.btn_copy.text() == "Copy Info"

    def test_close_button_exists(self, dlg):
        assert dlg.btn_close.text() == "Close"

    def test_info_browser_is_readonly(self, dlg):
        assert dlg.info_browser.isReadOnly()

    def test_formula_appears_in_rendered_html(self, dlg):
        html = dlg.info_browser.toHtml()
        assert "C9H8O4" in html

    def test_iupac_name_appears_in_html(self, dlg):
        html = dlg.info_browser.toHtml()
        assert "acetoxybenzoic" in html

    def test_empty_details_creates_without_error(self, qapp):
        d = _pubchem.MoleculeDetailsDialog({}, parent=None)
        assert d is not None
        d.destroy()

    def test_partial_details_no_error(self, qapp):
        d = _pubchem.MoleculeDetailsDialog({"Formula": "C6H6"}, parent=None)
        html = d.info_browser.toHtml()
        assert "C6H6" in html
        d.destroy()

    def test_dialog_size(self, dlg):
        assert dlg.width() >= 400
        assert dlg.height() >= 200


# ===========================================================================
# PubChemResolver  (pure Python — no Qt, no network for guard paths)
# ===========================================================================


class TestPubChemResolver:
    """Static methods on PubChemResolver that return early without network calls."""

    def test_resolve_name_empty_string_returns_error(self):
        smiles, err = _pubchem.PubChemResolver.resolve_name_to_smiles("")
        assert smiles is None
        assert "Empty" in err

    def test_resolve_name_none_returns_error(self):
        smiles, err = _pubchem.PubChemResolver.resolve_name_to_smiles(None)
        assert smiles is None

    def test_get_compound_details_empty_inchikey_returns_error(self):
        details, err = _pubchem.PubChemResolver.get_compound_details("")
        assert details is None
        assert "Empty" in err

    def test_get_compound_details_none_inchikey_returns_error(self):
        details, err = _pubchem.PubChemResolver.get_compound_details(None)
        assert details is None

    def test_base_url_is_pubchem(self):
        assert "pubchem" in _pubchem.PubChemResolver.BASE_URL
