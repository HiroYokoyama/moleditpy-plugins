"""
Tests for the OpenBabel Conversion Tool plugin (initialize -> export + drop
handler; OBABEL_AVAILABLE guard; open_file_with_openbabel no-extension early
exit; export_with_openbabel success/error paths).
"""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

from conftest import load_plugin as _load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
OBABEL_PATH = PLUGINS_DIR / "OpenBabel_Conversion_Tool" / "openbabel_conversion_tool.py"


def load_plugin(path):
    """Supply plain discovered formats for tests of the existing UI flows."""
    mod = _load_plugin(path)
    mod.pybel = MagicMock()
    mod.OBABEL_AVAILABLE = True
    return mod


class TestNativeEngineIsolation:
    def test_import_never_imports_native_engine(self):
        import builtins
        original_import = builtins.__import__

        def guarded_import(name, *args, **kwargs):
            assert name.split('.')[0] not in {'openbabel', 'pybel'}
            return original_import(name, *args, **kwargs)

        with mock_optional_imports(), patch('builtins.__import__', guarded_import):
            mod = _load_plugin(OBABEL_PATH)
        assert mod.OBABEL_AVAILABLE is None
        assert mod.pybel.informats == {}

    def test_format_discovery_uses_worker_once(self, tmp_path):
        import json
        with mock_optional_imports():
            mod = _load_plugin(OBABEL_PATH)

        def discover(script, args, parent):
            assert script == mod._FORMAT_SCRIPT
            Path(args[0]).write_text(json.dumps({
                'informats': {'pdb': 'PDB'}, 'outformats': {'xyz': 'XYZ'},
            }), encoding='utf-8')

        with patch.object(mod, '_run_conversion', side_effect=discover) as worker:
            assert mod._ensure_formats()
            assert mod._ensure_formats()
        worker.assert_called_once()
        assert mod.pybel.informats == {'pdb': 'PDB'}
        assert mod.pybel.outformats == {'xyz': 'XYZ'}

    def test_worker_failure_disables_handlers(self):
        with mock_optional_imports():
            mod = _load_plugin(OBABEL_PATH)
        ctx = make_context()
        with patch.object(mod, '_run_conversion', side_effect=ValueError('native crash')):
            mod.initialize(ctx)
        assert mod.OBABEL_AVAILABLE is False
        ctx.register_file_opener.assert_not_called()
        ctx.register_drop_handler.assert_not_called()

    def test_frozen_executable_never_recursively_launches_editor(self):
        import pytest
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
        with patch.object(mod.sys, 'frozen', True, create=True), patch.object(mod, 'QProcess') as process:
            with pytest.raises(ValueError, match='Python-installed MoleditPy'):
                mod._run_conversion('print(1)', [], None)
        process.assert_not_called()


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
            mod = load_plugin(OBABEL_PATH)
        mod._read_molecules = lambda path, fmt, parent: list(mod.pybel.readfile(fmt, path))
        return mod

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
            mod = load_plugin(OBABEL_PATH)
        def convert(script, args, parent):
            mol = mod.pybel.readstring('mol', Path(args[2]).read_text(encoding='utf-8'))
            mol.write(args[1], args[3], overwrite=True)
            Path(args[3]).write_text('converted', encoding='utf-8')
        mod._run_conversion = MagicMock(side_effect=convert)
        return mod

    def test_obabel_unavailable_warns_and_returns(self):
        mod = self._mod()
        mod.OBABEL_AVAILABLE = False
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        ctx.current_molecule.GetNumAtoms.return_value = 3
        ctx.current_molecule.GetNumAtoms.return_value = 3
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
        ctx.current_molecule.GetNumAtoms.return_value = 3
        ctx.current_molecule.GetNumAtoms.return_value = 3
        mod.pybel.outformats = {"xyz": "XYZ format"}
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("", "")):
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.export_with_openbabel(ctx)
        mock_warn.assert_not_called()

    def test_no_extension_warns(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        ctx.current_molecule.GetNumAtoms.return_value = 3
        ctx.current_molecule.GetNumAtoms.return_value = 3
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
        ctx.current_molecule.GetNumAtoms.return_value = 3
        ctx.current_molecule.GetNumAtoms.return_value = 3
        mod.pybel.outformats = {"xyz": "XYZ format"}
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("out.pdb", "")):
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.export_with_openbabel(ctx)
        mock_warn.assert_called_once()
        args, _ = mock_warn.call_args
        assert "not supported" in args[2].lower()

    def test_successful_export_reports_status_not_error(self, tmp_path):
        """Regression test for the PLUGIN_CONTEXT NameError bug (2026.07.10 fix)."""
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        ctx.current_molecule.GetNumAtoms.return_value = 3
        ctx.current_molecule.GetNumAtoms.return_value = 3
        mod.pybel.outformats = {"xyz": "XYZ format"}
        pybel_mol = MagicMock()
        mod.pybel.readstring.return_value = pybel_mol
        mod.Chem.MolToMolBlock.return_value = "MOLBLOCK"
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=(str(tmp_path / "out.xyz"), "")):
            with patch.object(mod.QMessageBox, "critical") as mock_critical:
                mod.export_with_openbabel(ctx)
        assert pybel_mol.write.call_args.args[0] == "xyz"
        assert (tmp_path / "out.xyz").read_text() == "converted"
        ctx.show_status_message.assert_called_once_with(f"Exported to {tmp_path / 'out.xyz'}", 3000)
        mock_critical.assert_not_called()

    def test_rdkit_to_molblock_failure_shows_critical(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        ctx.current_molecule.GetNumAtoms.return_value = 3
        ctx.current_molecule.GetNumAtoms.return_value = 3
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
        ctx.current_molecule.GetNumAtoms.return_value = 3
        ctx.current_molecule.GetNumAtoms.return_value = 3
        mod.pybel.outformats = {"xyz": "XYZ format"}
        mod.Chem.MolToMolBlock.return_value = "MOLBLOCK"
        mod.pybel.readstring.side_effect = Exception("obabel read failure")
        with patch.object(mod.QFileDialog, "getSaveFileName", return_value=("out.xyz", "")):
            with patch.object(mod.QMessageBox, "critical") as mock_critical:
                mod.export_with_openbabel(ctx)
        mock_critical.assert_called_once()
        args, _ = mock_critical.call_args
        assert "obabel read failure" in args[2]

    def test_unexpected_export_exception_shows_export_error(self):
        mod = self._mod()
        ctx = make_context()
        ctx.current_molecule = MagicMock()
        ctx.current_molecule.GetNumAtoms.return_value = 3
        ctx.current_molecule.GetNumAtoms.return_value = 3
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
    def _init_with_formats(self, informats, host_sdf=False):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            mod.pybel.informats = informats
            ctx = make_context()
            with patch.object(mod, "host_handles_sdf", return_value=host_sdf):
                mod.initialize(ctx)
        return mod, ctx

    def test_registers_opener_per_format_except_native(self):
        mod, ctx = self._init_with_formats(
            {"sdf": "MDL SDF", "pdb": "PDB", "mol": "MDL MOL", "xyz": "XYZ"}
        )
        registered = {call[0][0] for call in ctx.register_file_opener.call_args_list}
        assert registered == {".sdf", ".pdb"}  # .mol/.xyz left to the native loader

    def test_sdf_left_to_host_with_its_own_selector(self):
        # MoleditPy 4.11.0+ has a multi-molecule SDF selector of its own.
        mod, ctx = self._init_with_formats(
            {"sdf": "MDL SDF", "pdb": "PDB"}, host_sdf=True
        )
        registered = {call[0][0] for call in ctx.register_file_opener.call_args_list}
        assert registered == {".pdb"}

    def test_openers_registered_with_low_priority(self):
        mod, ctx = self._init_with_formats({"sdf": "MDL SDF"})
        for call in ctx.register_file_opener.call_args_list:
            assert call[1].get("priority") == -1


class TestDropHandler:
    def _get_handler(self, informats, host_sdf=False):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            mod.pybel.informats = informats
            ctx = make_context()
            with patch.object(mod, "host_handles_sdf", return_value=host_sdf):
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

    def test_sdf_claimed_on_old_host_even_if_missing_from_informats(self):
        # INTENTIONAL special case, do not "clean up": a host before 4.11.0
        # cannot pick a structure from a multi-molecule SDF, so the plugin
        # claims dropped .sdf files there (its own selection dialog).
        mod, handler = self._get_handler({"pdb": "PDB"})
        assert handler("multi.sdf") is True
        mod.open_file_with_openbabel.assert_called_once()

    def test_sdf_left_to_host_with_its_own_selector(self):
        mod, handler = self._get_handler({"sdf": "MDL SDF"}, host_sdf=True)
        assert handler("multi.sdf") is False
        mod.open_file_with_openbabel.assert_not_called()


class TestHostVersion:
    def _mod(self):
        with mock_optional_imports():
            return load_plugin(OBABEL_PATH)

    def test_version_parsed_from_host_constants(self):
        mod = self._mod()
        fake = MagicMock(VERSION="4.11.0")
        with patch.dict("sys.modules", {"moleditpy.utils.constants": fake}):
            assert mod._host_version() == (4, 11, 0)
            assert mod.host_handles_sdf() is True

    def test_older_host_keeps_plugin_sdf_handling(self):
        mod = self._mod()
        fake = MagicMock(VERSION="4.10.3")
        with patch.dict("sys.modules", {"moleditpy.utils.constants": fake}):
            assert mod.host_handles_sdf() is False

    def test_unparseable_or_missing_version_is_unknown(self):
        mod = self._mod()
        fake = MagicMock(VERSION="dev")
        with patch.dict("sys.modules", {"moleditpy.utils.constants": fake}):
            assert mod._host_version() is None
            assert mod.host_handles_sdf() is False
        with patch.dict(
            "sys.modules",
            {"moleditpy.utils.constants": None, "moleditpy_linux.utils.constants": None},
        ):
            assert mod._host_version() is None


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
    def test_uppercase_extension_lowercased_and_accepted(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            mod.pybel.outformats = {"xyz": "XYZ format"}
            ctx = make_context()
            ctx.current_molecule = MagicMock()
            ctx.current_molecule.GetNumAtoms.return_value = 3
            mod.Chem.MolToMolBlock.return_value = "MOLBLOCK"
            def convert(script, args, parent):
                Path(args[3]).write_text("converted")
            mod._run_conversion = convert
            with patch.object(
                mod.QFileDialog, "getSaveFileName", return_value=(str(tmp_path / "out.XYZ"), "")
            ), patch.object(mod.QMessageBox, "warning") as warn, patch.object(
                mod.QMessageBox, "critical"
            ) as crit:
                mod.export_with_openbabel(ctx)
            warn.assert_not_called()
            crit.assert_not_called()
            ctx.show_status_message.assert_called_once()
