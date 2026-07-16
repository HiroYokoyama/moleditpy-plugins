"""
GUI tests for the PubChem Structure Identifier plugin.

Covers: MoleculeDetailsDialog — pure Qt, builds HTML from a details dict.
No network calls in __init__. PubChemResolver static methods are tested
separately (pure Python, no Qt).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PUBCHEM_PATH = PLUGINS_DIR / "PubChem_Structure_Identifier" / "pubchem_structure_identifier.py"

with mock_chemistry_imports():
    _pubchem = load_plugin_for_gui(PUBCHEM_PATH)


_FULL_DETAILS = {
    "Common Name": "Aspirin",
    "IUPAC Name": "2-acetoxybenzoic acid",
    "Formula": "C9H8O4",
    "Mol. Weight": "180.16",
    "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
}


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
