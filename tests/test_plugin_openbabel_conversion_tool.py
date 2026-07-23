"""
Tests for the OpenBabel Conversion Tool plugin (initialize -> export + drop
handler; OBABEL_AVAILABLE guard; open_file_with_openbabel no-extension early
exit; export_with_openbabel success/error paths).
"""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
OBABEL_PATH = PLUGINS_DIR / "OpenBabel_Conversion_Tool" / "openbabel_conversion_tool.py"


class TestOpenBabelConversionTool:
    def test_initialize_obabel_available_registers_export_action(self):
        """When openbabel is mocked (available), the export action is registered."""
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_openbabel(self):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            label = ctx.add_export_action.call_args[0][0]
            assert "OpenBabel" in label

    def test_initialize_obabel_available_registers_drop_handler(self):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_drop_handler.assert_called_once()

    def test_initialize_obabel_unavailable_returns_early(self):
        """When OBABEL_AVAILABLE=False initialize() prints a warning and returns
        without calling any context registration method."""
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            mod.OBABEL_AVAILABLE = False
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_export_action.assert_not_called()
            ctx.register_file_opener.assert_not_called()
            ctx.register_drop_handler.assert_not_called()

    def test_open_file_no_extension_warns_and_returns(self):
        """open_file_with_openbabel() warns when the path has no file extension."""
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.open_file_with_openbabel("noextension", ctx)
            mock_warn.assert_called_once()

    def test_open_file_extensionless_does_not_reach_pybel(self):
        """Extensionless path exits before calling pybel.readfile."""
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            with patch.object(mod.QMessageBox, "warning"):
                mod.open_file_with_openbabel("noextension", ctx)


class TestOpenFileWithOpenBabelFullFlow:
    """Exercise the full read -> convert -> load pipeline in
    open_file_with_openbabel(); pybel/RDKit are mocked so all branches can be
    driven directly by controlling return values."""

    def _mod(self):
        with mock_optional_imports():
            return load_plugin(OBABEL_PATH)

    def test_no_molecules_in_file_warns(self):
        mod = self._mod()
        ctx = make_context()
        mod.pybel.readfile.return_value = []
        with patch.object(mod.QMessageBox, "warning") as warn:
            mod.open_file_with_openbabel("mol.xyz", ctx)
        warn.assert_called_once()
        args, _ = warn.call_args
        assert "No molecules" in args[2]

    def test_single_molecule_success_full_happy_path(self):
        mod = self._mod()
        ctx = make_context()
        mw = ctx.get_main_window.return_value

        fake_mol = MagicMock()
        fake_mol.write.return_value = "MOLBLOCK"
        mod.pybel.readfile.return_value = [fake_mol]

        rd_mol = MagicMock()
        # first check: no conformer -> Compute2DCoords; second check: has one -> WedgeMolBonds
        rd_mol.GetNumConformers.side_effect = [0, 1]
        mod.Chem.MolFromMolBlock.return_value = rd_mol

        mod.open_file_with_openbabel("mol.xyz", ctx)

        assert ctx.current_molecule is rd_mol
        ctx.push_undo_checkpoint.assert_called_once()
        ctx.reset_3d_camera.assert_called_once()
        ctx.enter_3d_viewer_mode.assert_called_once()
        ctx.show_status_message.assert_called_once()
        assert mw.init_manager.current_file_path == "mol.xyz"
        assert mw.state_manager.has_unsaved_changes is False
        ctx.refresh_ui.assert_called_once()
        mod.AllChem.Compute2DCoords.assert_called_once_with(rd_mol)
        mod.AllChem.WedgeMolBonds.assert_called_once()

    def test_ui_manager_fallback_used_when_no_enter_3d_viewer_mode(self):
        mod = self._mod()
        ctx = make_context()
        del ctx.enter_3d_viewer_mode
        mw = ctx.get_main_window.return_value

        fake_mol = MagicMock()
        fake_mol.write.return_value = "MOLBLOCK"
        mod.pybel.readfile.return_value = [fake_mol]

        rd_mol = MagicMock()
        rd_mol.GetNumConformers.return_value = 1
        mod.Chem.MolFromMolBlock.return_value = rd_mol

        mod.open_file_with_openbabel("mol.xyz", ctx)

        mw.ui_manager.enter_3d_viewer_ui_mode.assert_called_once()

    def test_3d_mode_switch_exception_is_swallowed(self):
        mod = self._mod()
        ctx = make_context()
        ctx.enter_3d_viewer_mode.side_effect = RuntimeError("boom")

        fake_mol = MagicMock()
        fake_mol.write.return_value = "MOLBLOCK"
        mod.pybel.readfile.return_value = [fake_mol]

        rd_mol = MagicMock()
        rd_mol.GetNumConformers.return_value = 1
        mod.Chem.MolFromMolBlock.return_value = rd_mol

        # Should not raise; the outer try/except logs and continues.
        mod.open_file_with_openbabel("mol.xyz", ctx)
        ctx.show_status_message.assert_called_once()

    def test_molblock_conversion_failure_shows_critical(self):
        mod = self._mod()
        ctx = make_context()
        fake_mol = MagicMock()
        fake_mol.write.side_effect = Exception("obabel write failed")
        mod.pybel.readfile.return_value = [fake_mol]
        with patch.object(mod.QMessageBox, "critical") as crit:
            mod.open_file_with_openbabel("mol.xyz", ctx)
        crit.assert_called_once()
        args, _ = crit.call_args
        assert "MolBlock" in args[2]

    def test_rdkit_mol_creation_failure_shows_critical(self):
        mod = self._mod()
        ctx = make_context()
        fake_mol = MagicMock()
        fake_mol.write.return_value = "MOLBLOCK"
        mod.pybel.readfile.return_value = [fake_mol]
        mod.Chem.MolFromMolBlock.return_value = None
        with patch.object(mod.QMessageBox, "critical") as crit:
            mod.open_file_with_openbabel("mol.xyz", ctx)
        crit.assert_called_once()
        args, _ = crit.call_args
        assert "Failed to create RDKit" in args[2]

    def test_multi_molecule_dialog_accepted_valid_index(self):
        mod = self._mod()
        ctx = make_context()
        mol0, mol1 = MagicMock(), MagicMock()
        mol1.write.return_value = "MOLBLOCK"
        mod.pybel.readfile.return_value = [mol0, mol1]

        rd_mol = MagicMock()
        rd_mol.GetNumConformers.return_value = 1
        mod.Chem.MolFromMolBlock.return_value = rd_mol

        fake_dialog = MagicMock()
        fake_dialog.exec.return_value = mod.QDialog.DialogCode.Accepted
        fake_dialog.get_selected_index.return_value = 1
        with patch.object(mod, "MoleculeSelectionDialog", return_value=fake_dialog):
            mod.open_file_with_openbabel("multi.sdf", ctx)

        mol1.write.assert_called_once_with("mol")
        assert ctx.current_molecule is rd_mol

    def test_multi_molecule_dialog_accepted_invalid_index_returns(self):
        mod = self._mod()
        ctx = make_context()
        mol0, mol1 = MagicMock(), MagicMock()
        mod.pybel.readfile.return_value = [mol0, mol1]

        fake_dialog = MagicMock()
        fake_dialog.exec.return_value = mod.QDialog.DialogCode.Accepted
        fake_dialog.get_selected_index.return_value = None
        with patch.object(mod, "MoleculeSelectionDialog", return_value=fake_dialog):
            mod.open_file_with_openbabel("multi.sdf", ctx)

        mol0.write.assert_not_called()
        mol1.write.assert_not_called()
        assert ctx.push_undo_checkpoint.call_count == 0

    def test_multi_molecule_dialog_cancelled_returns(self):
        mod = self._mod()
        ctx = make_context()
        mol0, mol1 = MagicMock(), MagicMock()
        mod.pybel.readfile.return_value = [mol0, mol1]

        fake_dialog = MagicMock()
        fake_dialog.exec.return_value = mod.QDialog.DialogCode.Rejected
        with patch.object(mod, "MoleculeSelectionDialog", return_value=fake_dialog):
            mod.open_file_with_openbabel("multi.sdf", ctx)

        mol0.write.assert_not_called()
        mol1.write.assert_not_called()

    def test_non_ascii_path_temp_file_created_and_removed(self, tmp_path):
        """Non-ASCII path triggers the temp-file workaround; the source file
        exists so shutil.copy2 succeeds and the temp file is cleaned up."""
        mod = self._mod()
        ctx = make_context()

        src_dir = tmp_path / "café"
        src_dir.mkdir()
        src = src_dir / "café.xyz"
        src.write_text("dummy")

        fake_mol = MagicMock()
        fake_mol.write.return_value = "MOLBLOCK"
        mod.pybel.readfile.return_value = [fake_mol]

        rd_mol = MagicMock()
        rd_mol.GetNumConformers.return_value = 1
        mod.Chem.MolFromMolBlock.return_value = rd_mol

        mod.open_file_with_openbabel(str(src), ctx)

        # readfile was called with a temp path, not the original non-ascii one
        called_path = mod.pybel.readfile.call_args[0][1]
        assert called_path != str(src)
        assert not os.path.exists(called_path)  # cleaned up after read
        assert ctx.current_molecule is rd_mol

    def test_non_ascii_path_copy_failure_falls_back_to_original_path(self):
        """If the temp-file copy fails (e.g. source doesn't exist), the code
        falls back to opening the original path directly."""
        mod = self._mod()
        ctx = make_context()

        fake_mol = MagicMock()
        fake_mol.write.return_value = "MOLBLOCK"
        mod.pybel.readfile.return_value = [fake_mol]

        rd_mol = MagicMock()
        rd_mol.GetNumConformers.return_value = 1
        mod.Chem.MolFromMolBlock.return_value = rd_mol

        missing_path = "does_not_exist_café.xyz"
        mod.open_file_with_openbabel(missing_path, ctx)

        called_path = mod.pybel.readfile.call_args[0][1]
        assert called_path == missing_path
        assert ctx.current_molecule is rd_mol

    def test_unexpected_exception_shows_import_error(self):
        mod = self._mod()
        ctx = make_context()
        mod.pybel.readfile.side_effect = RuntimeError("kaboom")
        with patch.object(mod.QMessageBox, "critical") as crit:
            mod.open_file_with_openbabel("mol.xyz", ctx)
        crit.assert_called_once()
        args, _ = crit.call_args
        assert "Import Error" in args[1]

    def test_open_file_wrapper_forwards_to_open_file_with_openbabel(self):
        """Covers the closure at line 61 registered as the file opener."""
        mod = self._mod()
        mod.pybel.informats = {"pdb": "PDB"}
        ctx = make_context()
        mod.initialize(ctx)
        wrapper = ctx.register_file_opener.call_args[0][1]
        with patch.object(mod, "open_file_with_openbabel") as inner:
            wrapper("thing.pdb")
        inner.assert_called_once_with("thing.pdb", ctx)


class TestExportWithOpenBabel:
    """
    Regression coverage for a real bug: the success branch referenced an
    undefined global `PLUGIN_CONTEXT` instead of the `context` parameter,
    which raised NameError, was swallowed by the enclosing `except
    Exception`, and surfaced as a spurious "Export Error" dialog even
    though the file had already been written successfully.
    """

    def _mod(self):
        with mock_optional_imports():
            return load_plugin(OBABEL_PATH)

    def test_obabel_unavailable_warns_and_returns(self):
        mod = self._mod()
        mod.OBABEL_AVAILABLE = False
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        with patch.object(mod.QMessageBox, "warning") as mock_warn:
            mod.export_with_openbabel(ctx)
        mock_warn.assert_called_once()

    def test_no_molecule_warns_and_returns(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = None
        with patch.object(mod.QMessageBox, "warning") as mock_warn:
            mod.export_with_openbabel(ctx)
        mock_warn.assert_called_once()
        args, _ = mock_warn.call_args
        assert "Export Error" in args[1]

    def test_user_cancels_dialog_no_further_action(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        mod.pybel.outformats = {"xyz": "XYZ format"}
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("", "")):
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.export_with_openbabel(ctx)
        mock_warn.assert_not_called()

    def test_no_extension_warns(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        mod.pybel.outformats = {"xyz": "XYZ format"}
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("no_ext_file", "")):
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.export_with_openbabel(ctx)
        mock_warn.assert_called_once()
        args, _ = mock_warn.call_args
        assert "extension" in args[2].lower()

    def test_unsupported_format_warns(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        mod.pybel.outformats = {"xyz": "XYZ format"}
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("out.pdb", "")):
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.export_with_openbabel(ctx)
        mock_warn.assert_called_once()
        args, _ = mock_warn.call_args
        assert "not supported" in args[2].lower()

    def test_successful_export_reports_status_not_error(self):
        """Regression test for the PLUGIN_CONTEXT NameError bug (2026.07.10 fix)."""
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        mod.pybel.outformats = {"xyz": "XYZ format"}
        pybel_mol = MagicMock()
        mod.pybel.readstring.return_value = pybel_mol
        mod.Chem.MolToMolBlock.return_value = "MOLBLOCK"
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("out.xyz", "")):
            with patch.object(mod.QMessageBox, "critical") as mock_critical:
                mod.export_with_openbabel(ctx)
        pybel_mol.write.assert_called_once_with("xyz", "out.xyz", overwrite=True)
        ctx.show_status_message.assert_called_once_with("Exported to out.xyz", 3000)
        mock_critical.assert_not_called()

    def test_rdkit_to_molblock_failure_shows_critical(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        mod.pybel.outformats = {"xyz": "XYZ format"}
        mod.Chem.MolToMolBlock.side_effect = Exception("rdkit failure")
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("out.xyz", "")):
            with patch.object(mod.QMessageBox, "critical") as mock_critical:
                mod.export_with_openbabel(ctx)
        mock_critical.assert_called_once()
        args, _ = mock_critical.call_args
        assert "MolBlock" in args[2]

    def test_pybel_readstring_failure_shows_critical(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        mod.pybel.outformats = {"xyz": "XYZ format"}
        mod.Chem.MolToMolBlock.return_value = "MOLBLOCK"
        mod.pybel.readstring.side_effect = Exception("obabel read failure")
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("out.xyz", "")):
            with patch.object(mod.QMessageBox, "critical") as mock_critical:
                mod.export_with_openbabel(ctx)
        mock_critical.assert_called_once()
        args, _ = mock_critical.call_args
        assert "read MolBlock" in args[2]

    def test_unexpected_export_exception_shows_export_error(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        mod.pybel.outformats = {"xyz": "XYZ format"}
        mod.Chem.MolToMolBlock.return_value = "MOLBLOCK"
        pybel_mol = MagicMock()
        pybel_mol.write.side_effect = RuntimeError("disk full")
        mod.pybel.readstring.return_value = pybel_mol
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("out.xyz", "")):
            with patch.object(mod.QMessageBox, "critical") as mock_critical:
                mod.export_with_openbabel(ctx)
        mock_critical.assert_called_once()
        args, _ = mock_critical.call_args
        assert "Export Error" in args[1]

    def test_export_wrapper_forwards_to_export_with_openbabel(self):
        """Covers the closure at line 106 registered as the export action."""
        mod = self._mod()
        ctx = make_context()
        mod.initialize(ctx)
        wrapper = ctx.add_export_action.call_args[0][1]
        with patch.object(mod, "export_with_openbabel") as inner:
            wrapper()
        inner.assert_called_once_with(ctx)


class TestFileOpenerRegistration:
    def _init_with_formats(self, informats):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            mod.pybel.informats = informats
            ctx = make_context()
            mod.initialize(ctx)
        return mod, ctx

    def test_registers_opener_per_format_except_native(self):
        mod, ctx = self._init_with_formats(
            {"sdf": "MDL SDF", "pdb": "PDB", "mol": "MDL MOL", "xyz": "XYZ"}
        )
        registered = {call[0][0] for call in ctx.register_file_opener.call_args_list}
        assert registered == {".sdf", ".pdb"}  # .mol/.xyz left to the native loader

    def test_openers_registered_with_low_priority(self):
        mod, ctx = self._init_with_formats({"sdf": "MDL SDF"})
        for call in ctx.register_file_opener.call_args_list:
            assert call[1].get("priority") == -1


class TestDropHandler:
    def _get_handler(self, informats):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            mod.pybel.informats = informats
            ctx = make_context()
            mod.initialize(ctx)
        mod.open_file_with_openbabel = MagicMock()
        handler = ctx.register_drop_handler.call_args[0][0]
        return mod, handler

    def test_supported_extension_handled(self):
        mod, handler = self._get_handler({"pdb": "PDB"})
        assert handler("protein.pdb") is True
        mod.open_file_with_openbabel.assert_called_once()

    def test_extension_case_insensitive(self):
        mod, handler = self._get_handler({"pdb": "PDB"})
        assert handler("PROTEIN.PDB") is True

    def test_native_mol_extension_not_handled(self):
        mod, handler = self._get_handler({"mol": "MDL MOL"})
        assert handler("thing.mol") is False
        mod.open_file_with_openbabel.assert_not_called()

    def test_unknown_extension_not_handled(self):
        mod, handler = self._get_handler({"pdb": "PDB"})
        assert handler("notes.foo") is False

    def test_sdf_always_claimed_even_if_missing_from_informats(self):
        # INTENTIONAL special case, do not "clean up": this plugin must always
        # claim dropped .sdf files (multi-structure SDF selection dialog) —
        # the main app cannot handle multi-structure SDF itself.
        mod, handler = self._get_handler({"pdb": "PDB"})
        assert handler("multi.sdf") is True
        mod.open_file_with_openbabel.assert_called_once()


class TestMoleculeSelectionDialog:
    def test_get_selected_index_none_when_no_selection(self):
        from conftest import extract_function

        fn = extract_function(
            OBABEL_PATH, "MoleculeSelectionDialog", "get_selected_index"
        )
        fake = MagicMock()
        fake.list_widget.selectedIndexes.return_value = []
        assert fn(fake) is None

    def test_get_selected_index_returns_row(self):
        from conftest import extract_function

        fn = extract_function(
            OBABEL_PATH, "MoleculeSelectionDialog", "get_selected_index"
        )
        fake = MagicMock()
        row = MagicMock()
        row.row.return_value = 3
        fake.list_widget.selectedIndexes.return_value = [row]
        assert fn(fake) == 3


class TestExportExtensionHandling:
    def test_uppercase_extension_lowercased_and_accepted(self):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            mod.pybel.outformats = {"xyz": "XYZ format"}
            ctx = make_context()
            ctx.current_molecule = MagicMock()
            with patch.object(
                mod.QFileDialog, "getSaveFileName", return_value=("out.XYZ", "")
            ), patch.object(mod.QMessageBox, "warning") as warn, patch.object(
                mod.QMessageBox, "critical"
            ) as crit:
                mod.export_with_openbabel(ctx)
            warn.assert_not_called()
            crit.assert_not_called()
            ctx.show_status_message.assert_called_once()
