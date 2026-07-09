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

from conftest import extract_function, load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
PUBCHEM_ID_PATH = PLUGINS_DIR / "PubChem_Structure_Identifier" / "pubchem_structure_identifier.py"

with mock_optional_imports():
    PUBCHEM_ID_MOD = load_plugin(PUBCHEM_ID_PATH)


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

    # -- get_compound_details: network paths ----------------------------------

    def test_get_compound_details_success(self):
        desc_payload = {"InformationList": {"Information": [{"Title": "Aspirin"}]}}
        props_payload = {
            "PropertyTable": {
                "Properties": [
                    {
                        "MolecularFormula": "C9H8O4",
                        "MolecularWeight": "180.16",
                        "IUPACName": "2-acetyloxybenzoic acid",
                    }
                ]
            }
        }
        with patch(
            "urllib.request.urlopen",
            side_effect=[_url_response(desc_payload), _url_response(props_payload)],
        ):
            details, error = self.Resolver.get_compound_details("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        assert error is None
        assert details["Common Name"] == "Aspirin"
        assert details["Formula"] == "C9H8O4"
        assert details["Mol. Weight"] == "180.16"
        assert details["IUPAC Name"] == "2-acetyloxybenzoic acid"
        assert details["InChIKey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    def test_get_compound_details_description_failure_falls_back_to_unknown(self):
        props_payload = {
            "PropertyTable": {"Properties": [{"MolecularFormula": "C9H8O4"}]}
        }
        with patch(
            "urllib.request.urlopen",
            side_effect=[OSError("desc unavailable"), _url_response(props_payload)],
        ):
            details, error = self.Resolver.get_compound_details("XYZ")
        assert error is None
        assert details["Common Name"] == "Unknown"
        assert details["Formula"] == "C9H8O4"

    def test_get_compound_details_no_properties_returns_error(self):
        desc_payload = {"InformationList": {"Information": []}}
        props_payload = {"PropertyTable": {"Properties": []}}
        with patch(
            "urllib.request.urlopen",
            side_effect=[_url_response(desc_payload), _url_response(props_payload)],
        ):
            details, error = self.Resolver.get_compound_details("XYZ")
        assert details is None
        assert "No properties" in error

    def test_get_compound_details_http_404(self):
        exc = urllib.error.HTTPError("", 404, "Not Found", {}, None)
        desc_payload = {"InformationList": {"Information": []}}
        with patch(
            "urllib.request.urlopen",
            side_effect=[_url_response(desc_payload), exc],
        ):
            details, error = self.Resolver.get_compound_details("XYZ")
        assert details is None
        assert "not found" in error.lower()


# ===========================================================================
# resolve_and_load() — callback for "Import from PubChem"
# ===========================================================================


class TestResolveAndLoad:
    def _ctx(self):
        ctx = MagicMock()
        mw = MagicMock()
        ctx.get_main_window.return_value = mw
        return ctx, mw

    def test_user_cancels_dialog_does_nothing(self):
        ctx, mw = self._ctx()
        with patch.object(PUBCHEM_ID_MOD.QInputDialog, "getText", return_value=("benzene", False)):
            PUBCHEM_ID_MOD.resolve_and_load(ctx)
        assert not mw.string_importer_manager.load_from_smiles.called

    def test_empty_name_does_nothing(self):
        ctx, mw = self._ctx()
        with patch.object(PUBCHEM_ID_MOD.QInputDialog, "getText", return_value=("   ", True)):
            PUBCHEM_ID_MOD.resolve_and_load(ctx)
        assert not mw.string_importer_manager.load_from_smiles.called

    def test_resolver_error_shows_warning(self):
        ctx, mw = self._ctx()
        with patch.object(PUBCHEM_ID_MOD.QInputDialog, "getText", return_value=("nonsense", True)):
            with patch("urllib.request.urlopen", side_effect=OSError("offline")):
                with patch.object(PUBCHEM_ID_MOD.QMessageBox, "warning") as mock_warn:
                    PUBCHEM_ID_MOD.resolve_and_load(ctx)
        mock_warn.assert_called_once()
        assert not mw.string_importer_manager.load_from_smiles.called

    def test_success_loads_smiles_and_shows_status(self):
        ctx, mw = self._ctx()
        payload = {"PropertyTable": {"Properties": [{"SMILES": "c1ccccc1"}]}}
        with patch.object(PUBCHEM_ID_MOD.QInputDialog, "getText", return_value=("benzene", True)):
            with patch("urllib.request.urlopen", return_value=_url_response(payload)):
                PUBCHEM_ID_MOD.resolve_and_load(ctx)
        mw.string_importer_manager.load_from_smiles.assert_called_once_with("c1ccccc1")
        ctx.show_status_message.assert_called_once()


# ===========================================================================
# identify_current_molecule() — callback for "Identify Molecule"
# ===========================================================================


class TestIdentifyCurrentMolecule:
    def _ctx(self, mol=MagicMock()):
        ctx = MagicMock()
        mw = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = mol
        return ctx, mw

    def test_no_rdkit_warns(self):
        ctx, mw = self._ctx()
        with patch.object(PUBCHEM_ID_MOD, "Chem", None):
            with patch.object(PUBCHEM_ID_MOD.QMessageBox, "warning") as mock_warn:
                PUBCHEM_ID_MOD.identify_current_molecule(ctx)
        mock_warn.assert_called_once()
        args = mock_warn.call_args[0]
        assert "RDKit" in args[2]

    def test_no_molecule_warns(self):
        ctx, mw = self._ctx(mol=None)
        with patch.object(PUBCHEM_ID_MOD.QMessageBox, "warning") as mock_warn:
            PUBCHEM_ID_MOD.identify_current_molecule(ctx)
        mock_warn.assert_called_once()
        args = mock_warn.call_args[0]
        assert "No valid molecule" in args[2]

    def test_inchikey_computation_failure_warns(self):
        ctx, mw = self._ctx()
        with patch.object(PUBCHEM_ID_MOD.Chem, "MolToInchiKey", side_effect=RuntimeError("bad mol")):
            with patch.object(PUBCHEM_ID_MOD.QMessageBox, "warning") as mock_warn:
                PUBCHEM_ID_MOD.identify_current_molecule(ctx)
        mock_warn.assert_called_once()
        args = mock_warn.call_args[0]
        assert "Failed to calculate InChIKey" in args[2]

    def test_success_shows_details_dialog(self):
        ctx, mw = self._ctx()
        payload = {"PropertyTable": {"Properties": [{"MolecularFormula": "C6H6"}]}}
        with patch.object(PUBCHEM_ID_MOD.Chem, "MolToInchiKey", return_value="XYZ-UHFFFAOYSA-N"):
            with patch(
                "urllib.request.urlopen",
                side_effect=[OSError("no desc"), _url_response(payload)],
            ):
                with patch.object(PUBCHEM_ID_MOD, "MoleculeDetailsDialog") as dlg_cls:
                    PUBCHEM_ID_MOD.identify_current_molecule(ctx)
        dlg_cls.assert_called_once()
        dlg_cls.return_value.exec.assert_called_once()

    def test_failure_shows_information_dialog(self):
        ctx, mw = self._ctx()
        with patch.object(PUBCHEM_ID_MOD.Chem, "MolToInchiKey", return_value="XYZ"):
            with patch("urllib.request.urlopen", side_effect=OSError("offline")):
                with patch.object(PUBCHEM_ID_MOD.QMessageBox, "information") as mock_info:
                    PUBCHEM_ID_MOD.identify_current_molecule(ctx)
        mock_info.assert_called_once()


# ===========================================================================
# MoleculeDetailsDialog.copy_to_clipboard() — ordered text formatting
# ===========================================================================


class TestCopyToClipboard:
    def _fn(self, qapp):
        return extract_function(
            PUBCHEM_ID_PATH,
            "MoleculeDetailsDialog",
            "copy_to_clipboard",
            extra_globals={"QApplication": qapp, "QTimer": MagicMock()},
        )

    def test_fields_ordered_and_formatted(self):
        qapp = MagicMock()
        fn = self._fn(qapp)
        self_stub = MagicMock()
        self_stub.details = {
            "InChIKey": "ABC-XYZ",
            "Formula": "C6H6",
            "Common Name": "Benzene",
        }
        with mock_optional_imports():
            fn(self_stub)
        text = qapp.clipboard.return_value.setText.call_args[0][0]
        lines = text.splitlines()
        assert lines[0] == "Common Name: Benzene"
        assert lines[1] == "Formula: C6H6"
        assert lines[2] == "InChIKey: ABC-XYZ"

    def test_missing_fields_skipped(self):
        qapp = MagicMock()
        fn = self._fn(qapp)
        self_stub = MagicMock()
        self_stub.details = {"Formula": "H2O"}
        with mock_optional_imports():
            fn(self_stub)
        text = qapp.clipboard.return_value.setText.call_args[0][0]
        assert text == "Formula: H2O"
