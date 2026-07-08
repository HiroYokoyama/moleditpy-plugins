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




# ---------------------------------------------------------------------------
# fetch_experimental_properties (mocked urllib)
# ---------------------------------------------------------------------------


def _pugview_payload(props):
    return {
        "Record": {
            "Section": [
                {
                    "TOCHeading": "Chemical and Physical Properties",
                    "Section": [
                        {
                            "TOCHeading": "Experimental Properties",
                            "Section": props,
                        }
                    ],
                }
            ]
        }
    }


class _FakeResponse:
    def __init__(self, payload, status=200):
        self.status = status
        self._body = json.dumps(payload).encode("utf-8")

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class TestFetchExperimentalProperties:
    def _patch(self, monkeypatch, payload, status=200):
        monkeypatch.setattr(
            COMPOUND_MOD.urllib.request,
            "urlopen",
            lambda url: _FakeResponse(payload, status=status),
        )

    def test_density_string_extracted(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Density",
                    "Information": [
                        {"Value": {"StringWithMarkup": [{"String": "1.234 g/cm3"}]}}
                    ],
                }
            ]
        )
        self._patch(monkeypatch, payload)
        density, desc = COMPOUND_MOD.PubChemFetcher.fetch_experimental_properties(241)
        assert density == "1.234 g/cm3"
        assert desc is None

    def test_density_number_value(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Density",
                    "Information": [{"Value": {"Number": [0.879]}}],
                }
            ]
        )
        self._patch(monkeypatch, payload)
        density, _ = COMPOUND_MOD.PubChemFetcher.fetch_experimental_properties(241)
        assert density == "0.879"

    def test_physical_description_matching_heuristic(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Physical Description",
                    "Information": [
                        {
                            "Value": {
                                "StringWithMarkup": [
                                    {"String": "White crystalline powder"}
                                ]
                            }
                        }
                    ],
                }
            ]
        )
        self._patch(monkeypatch, payload)
        _, desc = COMPOUND_MOD.PubChemFetcher.fetch_experimental_properties(241)
        assert desc == "White crystalline powder"

    def test_description_without_keywords_rejected(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Physical Description",
                    "Information": [
                        {"Value": {"StringWithMarkup": [{"String": "Pungent smell"}]}}
                    ],
                }
            ]
        )
        self._patch(monkeypatch, payload)
        _, desc = COMPOUND_MOD.PubChemFetcher.fetch_experimental_properties(241)
        assert desc is None

    def test_non_200_status_returns_none_pair(self, monkeypatch):
        self._patch(monkeypatch, _pugview_payload([]), status=404)
        result = COMPOUND_MOD.PubChemFetcher.fetch_experimental_properties(241)
        assert result == (None, None)

    def test_no_cid_short_circuits_without_network(self, monkeypatch):
        def _explode(url):
            raise AssertionError("network must not be touched")

        monkeypatch.setattr(COMPOUND_MOD.urllib.request, "urlopen", _explode)
        assert COMPOUND_MOD.PubChemFetcher.fetch_experimental_properties(None) == (
            None,
            None,
        )

    def test_first_matching_density_wins(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Density",
                    "Information": [
                        {"Value": {"StringWithMarkup": [{"String": "1.0 g/cm3"}]}}
                    ],
                },
                {
                    "TOCHeading": "Density",
                    "Information": [
                        {"Value": {"StringWithMarkup": [{"String": "2.0 g/cm3"}]}}
                    ],
                },
            ]
        )
        self._patch(monkeypatch, payload)
        density, _ = COMPOUND_MOD.PubChemFetcher.fetch_experimental_properties(241)
        assert density == "1.0 g/cm3"
