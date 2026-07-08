"""
Tests for the PubChem Structure Identifier plugin (PubChemResolver static
methods).

All tests run headlessly — PyQt6 / rdkit / etc. are replaced with MagicMock.
Network-dependent paths are exercised by patching urllib.request.urlopen.
"""

from __future__ import annotations

import json
import urllib.error
from pathlib import Path
from unittest.mock import MagicMock, patch

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

with mock_optional_imports():
    PUBCHEM_ID_MOD = load_plugin(
        PLUGINS_DIR / "PubChem_Structure_Identifier" / "pubchem_structure_identifier.py"
    )


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
