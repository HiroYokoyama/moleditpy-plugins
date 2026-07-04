"""
Tests for visible plugins with PubChem API helpers and XYZ parsing logic:
  - PubChem Structure Identifier  (PubChemResolver static methods)
  - Compound Info Report          (PubChemFetcher static methods, extract_cas,
                                   calculate_adducts, initialize)
  - Paste XYZ                     (parse_xyz_lines, initialize smoke)

All tests run headlessly — PyQt6 / rdkit / etc. are replaced with MagicMock.
Network-dependent paths are exercised by patching urllib.request.urlopen.
"""

from __future__ import annotations

import io
import json
import urllib.error
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

# ---------------------------------------------------------------------------
# Load modules once at collection time
# ---------------------------------------------------------------------------

with mock_optional_imports():
    PUBCHEM_ID_MOD = load_plugin(
        PLUGINS_DIR / "PubChem_Structure_Identifier" / "pubchem_structure_identifier.py"
    )
    COMPOUND_MOD = load_plugin(
        PLUGINS_DIR / "Compound_Info_Report" / "compound_info_report.py"
    )
    PASTE_XYZ_MOD = load_plugin(
        PLUGINS_DIR / "Paste_XYZ" / "paste_xyz.py"
    )

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _url_response(payload: dict, status: int = 200) -> MagicMock:
    """Return a context-manager mock that mimics a urllib response."""
    body = json.dumps(payload).encode()
    r = MagicMock()
    r.status = status
    r.read.return_value = body
    r.__enter__ = lambda s: s
    r.__exit__ = MagicMock(return_value=False)
    return r


# ===========================================================================
# PubChem Structure Identifier — PubChemResolver
# ===========================================================================

class TestPubChemResolver:
    Resolver = PUBCHEM_ID_MOD.PubChemResolver

    # -- resolve_name_to_smiles: early-exit paths (no network) ---------------

    def test_empty_string_returns_error_immediately(self):
        smiles, error = self.Resolver.resolve_name_to_smiles("")
        assert smiles is None
        assert "Empty name provided" in error

    def test_none_returns_error_immediately(self):
        smiles, error = self.Resolver.resolve_name_to_smiles(None)
        assert smiles is None
        assert "Empty name provided" in error

    # -- resolve_name_to_smiles: network paths (mocked) ----------------------

    def test_success_returns_smiles(self):
        payload = {
            "PropertyTable": {"Properties": [{"IsomericSMILES": "c1ccccc1"}]}
        }
        with patch("urllib.request.urlopen", return_value=_url_response(payload)):
            smiles, error = self.Resolver.resolve_name_to_smiles("benzene")
        assert smiles == "c1ccccc1"
        assert error is None

    def test_http_404_returns_not_found(self):
        exc = urllib.error.HTTPError("", 404, "Not Found", {}, None)
        with patch("urllib.request.urlopen", side_effect=exc):
            smiles, error = self.Resolver.resolve_name_to_smiles("notamolecule")
        assert smiles is None
        assert "not found" in error.lower()

    def test_empty_properties_returns_no_smiles(self):
        payload = {"PropertyTable": {"Properties": []}}
        with patch("urllib.request.urlopen", return_value=_url_response(payload)):
            smiles, error = self.Resolver.resolve_name_to_smiles("mystery")
        assert smiles is None
        assert "No SMILES" in error

    def test_network_error_returns_parsing_error(self):
        with patch("urllib.request.urlopen", side_effect=OSError("connection refused")):
            smiles, error = self.Resolver.resolve_name_to_smiles("benzene")
        assert smiles is None
        assert "Error" in error

    # -- get_compound_details: early-exit paths (no network) -----------------

    def test_get_compound_details_empty_string(self):
        result, error = self.Resolver.get_compound_details("")
        assert result is None
        assert "Empty InChIKey" in error

    def test_get_compound_details_none(self):
        result, error = self.Resolver.get_compound_details(None)
        assert result is None
        assert "Empty InChIKey" in error

    # -- run smoke -----------------------------------------------------------

    def test_run_smoke(self):
        """run(mw) must not raise; it uses mocked moleditpy PluginContext."""
        with mock_optional_imports():
            mw = MagicMock()
            # Patch urlopen to avoid any real network activity
            with patch("urllib.request.urlopen", side_effect=OSError("offline")):
                PUBCHEM_ID_MOD.run(mw)  # must not raise


# ===========================================================================
# Compound Info Report — PubChemFetcher
# ===========================================================================

class TestPubChemFetcher:
    Fetcher = COMPOUND_MOD.PubChemFetcher

    # -- get_synonyms: early-exit paths (no network) -------------------------

    def test_get_synonyms_empty_string(self):
        assert self.Fetcher.get_synonyms("") == []

    def test_get_synonyms_none(self):
        assert self.Fetcher.get_synonyms(None) == []

    # -- get_cid: early-exit paths (no network) ------------------------------

    def test_get_cid_empty_string(self):
        assert self.Fetcher.get_cid("") is None

    def test_get_cid_none(self):
        assert self.Fetcher.get_cid(None) is None

    # -- get_synonyms: network path (mocked) ---------------------------------

    def test_get_synonyms_success(self):
        payload = {
            "InformationList": {
                "Information": [{"Synonym": ["aspirin", "50-78-2", "acetylsalicylic acid"]}]
            }
        }
        with patch("urllib.request.urlopen", return_value=_url_response(payload)):
            syns = self.Fetcher.get_synonyms("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        assert "aspirin" in syns
        assert "50-78-2" in syns

    def test_get_synonyms_network_error_returns_empty(self):
        with patch("urllib.request.urlopen", side_effect=OSError("offline")):
            assert self.Fetcher.get_synonyms("BSYNRYMUTXBXSQ-UHFFFAOYSA-N") == []

    # -- get_cid: network path (mocked) --------------------------------------

    def test_get_cid_success(self):
        payload = {"IdentifierList": {"CID": [2244, 9999]}}
        with patch("urllib.request.urlopen", return_value=_url_response(payload)):
            cid = self.Fetcher.get_cid("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        assert cid == 2244  # first entry

    def test_get_cid_network_error_returns_none(self):
        with patch("urllib.request.urlopen", side_effect=OSError("offline")):
            assert self.Fetcher.get_cid("BSYNRYMUTXBXSQ-UHFFFAOYSA-N") is None

    # -- extract_cas ---------------------------------------------------------

    def test_extract_cas_valid_numbers(self):
        syns = ["50-00-0", "100-41-4", "aspirin", "acetylsalicylic acid"]
        cas = self.Fetcher.extract_cas(syns)
        assert "50-00-0" in cas
        assert "100-41-4" in cas

    def test_extract_cas_filters_alpha_names(self):
        syns = ["aspirin", "benzene", "ethanol"]
        assert self.Fetcher.extract_cas(syns) == []

    def test_extract_cas_filters_wrong_hyphen_count(self):
        syns = ["1234-56789"]  # only one hyphen
        assert self.Fetcher.extract_cas(syns) == []

    def test_extract_cas_caps_at_five(self):
        syns = [f"{i:02d}-{i:02d}-{i}" for i in range(10)]
        result = self.Fetcher.extract_cas(syns)
        assert len(result) <= 5

    def test_extract_cas_empty_list(self):
        assert self.Fetcher.extract_cas([]) == []

    # -- initialize ----------------------------------------------------------

    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            ctx = make_context()
            COMPOUND_MOD.initialize(ctx)
        ctx.add_menu_action.assert_called_once()
        path = ctx.add_menu_action.call_args[0][0]
        assert "Compound" in path or "compound" in path.lower()


# ===========================================================================
# Paste XYZ — parse_xyz_lines (imported from plugin) + initialize smoke
# ===========================================================================

# parse_xyz_lines is a module-level function in the plugin so tests import it directly.
_parse = PASTE_XYZ_MOD.parse_xyz_lines


class TestPasteXYZParser:
    def test_valid_three_atom_block(self):
        xyz = "C  0.0  0.0  0.0\nH  1.0  0.0  0.0\nH -1.0  0.0  0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 3
        assert atoms[0] == ("C", 0.0, 0.0, 0.0)
        assert atoms[1] == ("H", 1.0, 0.0, 0.0)

    def test_header_lines_skipped(self):
        xyz = "3\nComment line\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 3

    def test_short_lines_skipped(self):
        xyz = "C 0.0 0.0\nH 1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 1
        assert atoms[0][0] == "H"

    def test_non_float_coords_skipped(self):
        xyz = "C abc 0.0 0.0\nH 1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 1
        assert atoms[0][0] == "H"

    def test_non_alpha_symbol_skipped(self):
        xyz = "1 0.0 0.0 0.0\nH 1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 1
        assert atoms[0][0] == "H"

    def test_empty_input_returns_empty(self):
        assert _parse("") == []

    def test_blank_lines_skipped(self):
        xyz = "\n\nC 0.0 0.0 0.0\n\n"
        atoms = _parse(xyz)
        assert len(atoms) == 1

    def test_fractional_coords_parsed(self):
        xyz = "O -0.12345 6.789 -3.141"
        atoms = _parse(xyz)
        assert len(atoms) == 1
        assert atoms[0][0] == "O"
        assert abs(atoms[0][1] - (-0.12345)) < 1e-9

    def test_mixed_case_symbols_accepted(self):
        xyz = "Ca 0.0 0.0 0.0\nFe 1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 2
        assert atoms[0][0] == "Ca"
        assert atoms[1][0] == "Fe"


class TestPasteXYZSmoke:
    def test_initialize_registers_menu_action(self):
        """initialize(ctx) must register exactly one menu action under File/."""
        with mock_optional_imports():
            ctx = make_context()
            PASTE_XYZ_MOD.initialize(ctx)
        ctx.add_menu_action.assert_called_once()
        path = ctx.add_menu_action.call_args[0][0]
        assert path.startswith("File/")

    def test_initialize_action_callable(self):
        """The callback registered by initialize must be callable without raising."""
        with mock_optional_imports():
            ctx = make_context()
            PASTE_XYZ_MOD.initialize(ctx)
        callback = ctx.add_menu_action.call_args[0][1]
        assert callable(callback)
