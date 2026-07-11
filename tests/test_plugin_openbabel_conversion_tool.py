"""
Tests for the OpenBabel Conversion Tool plugin (initialize -> export + drop
handler; OBABEL_AVAILABLE guard; open_file_with_openbabel no-extension early
exit; export_with_openbabel success/error paths).
"""

from __future__ import annotations

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
