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
