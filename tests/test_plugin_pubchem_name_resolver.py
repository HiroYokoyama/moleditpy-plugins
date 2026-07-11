"""
Tests for the PubChem Name Resolver plugin.
"""

from __future__ import annotations

import types
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
PUBCHEM_PATH = PLUGINS_DIR / "PubChem_Name_Ressolver" / "pubchem_ressolver.py"


class TestPubChemNameResolver:
    def test_initialize_does_not_raise(self):
        """initialize() only stores the context; run(mw) handles auto-registration."""
        with mock_optional_imports():
            mod = load_plugin(PUBCHEM_PATH)
            ctx = make_context()
            mod.initialize(ctx)  # must not raise

    def test_run_search_empty_input_does_nothing(self):
        """run_search() with empty text must return without setting lbl_info."""
        with mock_optional_imports():
            mod = load_plugin(PUBCHEM_PATH)
            ctx = make_context()
            dialog = mod.MoleculeResolverDialog(ctx)
            dialog.line_input = MagicMock()
            dialog.line_input.text.return_value = "   "
            dialog.lbl_info = MagicMock()
            dialog.run_search()
            dialog.lbl_info.setText.assert_not_called()

    def test_run_creates_dialog_on_first_call(self):
        """run(mw) creates a MoleculeResolverDialog when no window is registered."""
        with mock_optional_imports():
            mod = load_plugin(PUBCHEM_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            mod.initialize(ctx)
            # run(mw) is the auto-registered entry point
            mod.run(ctx.get_main_window())


def _make_requests_stub():
    """Stub requests module with a real RequestException class."""
    m = types.ModuleType("requests")

    class RequestException(Exception):
        pass

    m.exceptions = SimpleNamespace(RequestException=RequestException)
    m.get = MagicMock()
    return m


def _pubchem_search_env(requests_stub, chem=None):
    globs = {
        "requests": requests_stub,
        "Chem": chem or MagicMock(),
        "QMessageBox": MagicMock(),
        "QApplication": MagicMock(),
        "Qt": MagicMock(),
        "PLUGIN_NAME": "PubChem Name Resolver",
    }
    fn = extract_function(PUBCHEM_PATH, "MoleculeResolverDialog", "run_search", globs)
    return fn, globs


def _pubchem_self(query, search_type="Auto (Name/CAS)"):
    s = SimpleNamespace()
    s.line_input = MagicMock()
    s.line_input.text.return_value = query
    s.combo_type = MagicMock()
    s.combo_type.currentText.return_value = search_type
    s.lbl_info = MagicMock()
    s.btn_search = MagicMock()
    s.table = MagicMock()
    s.context = MagicMock()
    s.candidates_data = []
    s.update_table = MagicMock()
    return s


def _http_response(status_code, payload=None):
    return SimpleNamespace(status_code=status_code, json=lambda: payload or {})


class TestPubChemRunSearch:
    def test_empty_query_returns_immediately(self):
        req = _make_requests_stub()
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("   ")
        fn(s)
        s.btn_search.setEnabled.assert_not_called()
        req.get.assert_not_called()

    def test_http_200_populates_candidates(self):
        req = _make_requests_stub()
        req.get.return_value = _http_response(200, {
            "PropertyTable": {"Properties": [
                {"IsomericSMILES": "CCO", "Title": "Ethanol",
                 "MolecularFormula": "C2H6O"},
            ]}
        })
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("ethanol")
        fn(s)
        assert s.candidates_data == [
            {"name": "Ethanol", "smiles": "CCO", "formula": "C2H6O"}
        ]
        s.update_table.assert_called_once()
        assert "Found 1 candidates" in s.lbl_info.setText.call_args[0][0]

    def test_canonical_smiles_fallback(self):
        req = _make_requests_stub()
        req.get.return_value = _http_response(200, {
            "PropertyTable": {"Properties": [
                {"CanonicalSMILES": "CC", "Title": "Ethane",
                 "MolecularFormula": "C2H6"},
            ]}
        })
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("ethane")
        fn(s)
        assert s.candidates_data[0]["smiles"] == "CC"

    def test_any_smiles_key_fallback(self):
        req = _make_requests_stub()
        req.get.return_value = _http_response(200, {
            "PropertyTable": {"Properties": [
                {"ConnectivitySMILES": "CN", "Title": "X",
                 "MolecularFormula": "CH5N"},
            ]}
        })
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("methylamine")
        fn(s)
        assert s.candidates_data[0]["smiles"] == "CN"

    def test_http_404_reports_not_found(self):
        req = _make_requests_stub()
        req.get.return_value = _http_response(404)
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("nosuchcompound")
        fn(s)
        assert s.lbl_info.setText.call_args[0][0] == "Not found in PubChem search."
        s.update_table.assert_not_called()

    def test_network_error_shows_dialog(self):
        req = _make_requests_stub()
        req.get.side_effect = req.exceptions.RequestException("offline")
        fn, globs = _pubchem_search_env(req)
        s = _pubchem_self("water")
        fn(s)
        globs["QMessageBox"].critical.assert_called_once()
        assert s.lbl_info.setText.call_args[0][0] == "Network error."

    def test_search_button_reenabled_after_error(self):
        req = _make_requests_stub()
        req.get.side_effect = req.exceptions.RequestException("offline")
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("water")
        fn(s)
        s.btn_search.setEnabled.assert_called_with(True)

    def test_smiles_mode_valid_input(self):
        req = _make_requests_stub()
        chem = MagicMock()
        chem.rdMolDescriptors.CalcMolFormula.return_value = "C6H6"
        fn, _ = _pubchem_search_env(req, chem=chem)
        s = _pubchem_self("c1ccccc1", search_type="SMILES")
        fn(s)
        assert s.candidates_data == [
            {"name": "User Input SMILES", "smiles": "c1ccccc1", "formula": "C6H6"}
        ]
        req.get.assert_not_called()  # no network in SMILES mode

    def test_smiles_mode_invalid_input_reports_error(self):
        req = _make_requests_stub()
        chem = MagicMock()
        chem.MolFromSmiles.return_value = None
        fn, globs = _pubchem_search_env(req, chem=chem)
        s = _pubchem_self("not_smiles", search_type="SMILES")
        fn(s)
        globs["QMessageBox"].critical.assert_called_once()
        assert "Invalid SMILES" in globs["QMessageBox"].critical.call_args[0][2]
        assert s.lbl_info.setText.call_args[0][0] == "Error occurred."


class TestPubChemLoadMolecule:
    def _env(self):
        globs = {
            "QMessageBox": MagicMock(),
            "QApplication": MagicMock(),
            "Qt": MagicMock(),
            "PLUGIN_NAME": "PubChem Name Resolver",
        }
        fn = extract_function(
            PUBCHEM_PATH, "MoleculeResolverDialog", "load_molecule", globs
        )
        return fn, globs

    def _self(self, selected_row=0, smiles="CCO"):
        s = SimpleNamespace()
        item = MagicMock()
        item.row.return_value = selected_row
        s.table = MagicMock()
        s.table.selectedItems.return_value = [item] if selected_row is not None else []
        s.candidates_data = [{"name": "Ethanol", "smiles": smiles}]
        s.lbl_info = MagicMock()
        s.context = MagicMock()
        s.accept = MagicMock()
        return s

    def test_no_selection_warns(self):
        fn, globs = self._env()
        s = self._self(selected_row=None)
        fn(s)
        globs["QMessageBox"].warning.assert_called_once()
        s.accept.assert_not_called()

    def test_empty_smiles_warns(self):
        fn, globs = self._env()
        s = self._self(smiles="")
        fn(s)
        globs["QMessageBox"].warning.assert_called_once()
        s.accept.assert_not_called()

    def test_success_loads_via_string_importer(self):
        fn, globs = self._env()
        s = self._self()
        mw = MagicMock()
        s.context.get_main_window.return_value = mw
        fn(s)
        mw.string_importer_manager.load_from_smiles.assert_called_once_with("CCO")
        s.accept.assert_called_once()
        globs["QMessageBox"].information.assert_called_once()

    def test_missing_importer_reports_error(self):
        fn, globs = self._env()
        s = self._self()
        mw = SimpleNamespace()  # no string_importer_manager attribute
        s.context.get_main_window.return_value = mw
        fn(s)
        globs["QMessageBox"].critical.assert_called_once()
        s.accept.assert_not_called()
