"""
Headless GUI tests for the PubChem Name Resolver plugin.

Covers: MoleculeResolverDialog.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PUBCHEM_PATH = PLUGINS_DIR / "PubChem_Name_Ressolver" / "pubchem_ressolver.py"

with mock_chemistry_imports():
    _pubchem = load_plugin_for_gui(PUBCHEM_PATH)


# ===========================================================================
# MoleculeResolverDialog  (PubChem Name Resolver)
# ===========================================================================


class TestMoleculeResolverDialog:
    @pytest.fixture
    def dlg(self, qapp):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        d = _pubchem.MoleculeResolverDialog(ctx)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "PubChem Name Resolver"

    def test_search_type_options(self, dlg):
        items = [dlg.combo_type.itemText(i) for i in range(dlg.combo_type.count())]
        assert items == ["Auto (Name/CAS)", "SMILES"]

    def test_table_columns(self, dlg):
        assert dlg.table.columnCount() == 3
        headers = [
            dlg.table.horizontalHeaderItem(i).text() for i in range(3)
        ]
        assert headers == ["Name/Synonym", "Formula", "SMILES"]

    def test_registered_as_main_panel(self, dlg):
        dlg.context.register_window.assert_called_with("main_panel", dlg)

    def test_empty_search_is_noop(self, dlg):
        dlg.line_input.setText("   ")
        dlg.run_search()
        assert dlg.btn_search.isEnabled()
        assert dlg.lbl_info.text() == "Enter a chemical identifier and click Search."

    def test_smiles_search_populates_table(self, dlg):
        dlg.combo_type.setCurrentText("SMILES")
        dlg.line_input.setText("CCO")
        dlg.run_search()
        assert dlg.table.rowCount() == 1
        assert dlg.table.item(0, 0).text() == "User Input SMILES"
        assert dlg.table.item(0, 2).text() == "CCO"
        assert "Found 1 candidates" in dlg.lbl_info.text()
        assert dlg.btn_search.isEnabled()

    def test_invalid_smiles_shows_error(self, dlg):
        dlg.combo_type.setCurrentText("SMILES")
        dlg.line_input.setText("not-a-smiles")

        # The mocked requests module's RequestException is not a real exception
        # class; the except clause needs one to let the ValueError pass through.
        class _FakeRequestException(Exception):
            pass

        with patch.object(_pubchem.Chem, "MolFromSmiles", return_value=None), \
                patch.object(
                    _pubchem.requests.exceptions,
                    "RequestException",
                    _FakeRequestException,
                ), \
                patch.object(_pubchem, "QMessageBox") as mb:
            dlg.run_search()
        mb.critical.assert_called_once()
        assert "Invalid SMILES" in mb.critical.call_args.args[2]
        assert dlg.lbl_info.text() == "Error occurred."

    def test_load_without_selection_warns(self, dlg):
        with patch.object(_pubchem, "QMessageBox") as mb:
            dlg.load_molecule()
        mb.warning.assert_called_once()

    def test_load_selected_row_uses_string_importer(self, dlg):
        dlg.candidates_data = [
            {"name": "Ethanol", "smiles": "CCO", "formula": "C2H6O"}
        ]
        dlg.update_table()
        dlg.table.selectRow(0)
        mw = MagicMock()
        dlg.context.get_main_window.return_value = mw
        with patch.object(_pubchem, "QMessageBox"):
            dlg.load_molecule()
        mw.string_importer_manager.load_from_smiles.assert_called_once_with("CCO")

    def test_load_selected_row_empty_smiles_warns(self, dlg):
        dlg.candidates_data = [{"name": "Mystery", "smiles": "", "formula": ""}]
        dlg.update_table()
        dlg.table.selectRow(0)
        with patch.object(_pubchem, "QMessageBox") as mb:
            dlg.load_molecule()
        mb.warning.assert_called_once()

    def test_load_selected_row_missing_importer_reports_error(self, dlg):
        dlg.candidates_data = [
            {"name": "Ethanol", "smiles": "CCO", "formula": "C2H6O"}
        ]
        dlg.update_table()
        dlg.table.selectRow(0)

        class _NoManager:
            pass

        dlg.context.get_main_window.return_value = _NoManager()
        with patch.object(_pubchem, "QMessageBox") as mb:
            dlg.load_molecule()
        mb.critical.assert_called_once()
        assert "does not support" in mb.critical.call_args.args[2]

    def test_load_selected_row_importer_exception_reports_error(self, dlg):
        dlg.candidates_data = [
            {"name": "Ethanol", "smiles": "CCO", "formula": "C2H6O"}
        ]
        dlg.update_table()
        dlg.table.selectRow(0)
        mw = MagicMock()
        mw.string_importer_manager.load_from_smiles.side_effect = RuntimeError("boom")
        dlg.context.get_main_window.return_value = mw
        with patch.object(_pubchem, "QMessageBox") as mb:
            dlg.load_molecule()
        mb.critical.assert_called_once()
        assert "boom" in mb.critical.call_args.args[2]
        assert dlg.lbl_info.text() == "Load failed."

    def test_auto_search_http_200_populates_table(self, dlg):
        dlg.combo_type.setCurrentText("Auto (Name/CAS)")
        dlg.line_input.setText("ethanol")
        response = MagicMock()
        response.status_code = 200
        response.json.return_value = {
            "PropertyTable": {
                "Properties": [
                    {
                        "SMILES": "CCO",
                        "Title": "Ethanol",
                        "MolecularFormula": "C2H6O",
                    }
                ]
            }
        }
        with patch.object(_pubchem.requests, "get", return_value=response):
            dlg.run_search()
        assert dlg.table.rowCount() == 1
        assert dlg.table.item(0, 0).text() == "Ethanol"
        assert dlg.table.item(0, 2).text() == "CCO"
        assert "Found 1 candidates" in dlg.lbl_info.text()

    def test_auto_search_http_404_reports_not_found(self, dlg):
        dlg.combo_type.setCurrentText("Auto (Name/CAS)")
        dlg.line_input.setText("nosuchcompound")
        response = MagicMock()
        response.status_code = 404
        with patch.object(_pubchem.requests, "get", return_value=response):
            dlg.run_search()
        assert dlg.lbl_info.text() == "Not found in PubChem search."
        assert dlg.table.rowCount() == 0

    def test_auto_search_network_error_shows_dialog(self, dlg):
        dlg.combo_type.setCurrentText("Auto (Name/CAS)")
        dlg.line_input.setText("water")

        class _FakeRequestException(Exception):
            pass

        with patch.object(
            _pubchem.requests, "get", side_effect=_FakeRequestException("offline")
        ), patch.object(
            _pubchem.requests.exceptions, "RequestException", _FakeRequestException
        ), patch.object(_pubchem, "QMessageBox") as mb:
            dlg.run_search()
        mb.critical.assert_called_once()
        assert dlg.lbl_info.text() == "Network error."

    def test_auto_search_generic_exception_reports_error(self, dlg):
        dlg.combo_type.setCurrentText("Auto (Name/CAS)")
        dlg.line_input.setText("water")

        class _FakeRequestException(Exception):
            pass

        with patch.object(
            _pubchem.requests, "get", side_effect=ValueError("kaboom")
        ), patch.object(
            _pubchem.requests.exceptions, "RequestException", _FakeRequestException
        ), patch.object(_pubchem, "QMessageBox") as mb:
            dlg.run_search()
        mb.critical.assert_called_once()
        assert "kaboom" in mb.critical.call_args.args[2]
        assert dlg.lbl_info.text() == "Error occurred."


class TestPubChemRunEntryPoint:
    def test_run_with_no_context_does_nothing(self):
        original_context = _pubchem.PLUGIN_CONTEXT
        _pubchem.PLUGIN_CONTEXT = None
        try:
            mw = MagicMock()
            _pubchem.run(mw)  # must not raise, must return early
        finally:
            _pubchem.PLUGIN_CONTEXT = original_context
