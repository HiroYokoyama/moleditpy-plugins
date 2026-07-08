"""
Tests for the Compound Info Report plugin (PubChemFetcher static methods,
extract_cas, initialize).

All tests run headlessly — PyQt6 / rdkit / etc. are replaced with MagicMock.
Network-dependent paths are exercised by patching urllib.request.urlopen.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

with mock_optional_imports():
    COMPOUND_MOD = load_plugin(
        PLUGINS_DIR / "Compound_Info_Report" / "compound_info_report.py"
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

    # -- extract_cas -----------------------------------------------------------

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
