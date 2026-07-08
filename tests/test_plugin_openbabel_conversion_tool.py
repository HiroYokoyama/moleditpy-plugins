"""
Tests for the OpenBabel Conversion Tool plugin (initialize -> export + drop
handler; OBABEL_AVAILABLE guard; open_file_with_openbabel no-extension early exit).
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

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
