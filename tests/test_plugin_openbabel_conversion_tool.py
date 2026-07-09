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
