"""
Tests for the Conformational Search plugin.

Covers:
  1. initialize() must register at least one menu/export/plugin action
  2. No-molecule guard paths (run_plugin with mol=None)
  3. Dialog accept/reject round-trips
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CONF_SEARCH_PATH = PLUGINS_DIR / "Conformational_Search" / "conf_search.py"


def _menu_registered(ctx: MagicMock) -> bool:
    """Return True if initialize() called any recognised registration method."""
    return (
        ctx.add_menu_action.called
        or ctx.add_export_action.called
        or ctx.add_plugin_menu.called
        or ctx.add_analysis_tool.called
        or ctx.add_toolbar_action.called
    )


class TestConformationalSearch:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert _menu_registered(ctx), (
                "Conformational Search initialize() must call add_menu_action()"
            )

    def test_initialize_menu_path_contains_conformational(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            call_args = ctx.add_menu_action.call_args
            assert call_args is not None
            path = call_args[0][0]
            assert "Conformational" in path or "conformational" in path.lower()

    def test_run_plugin_no_mol_warns(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            mock_warn.assert_called_once()

    def test_dialog_accept_does_not_raise(self):
        """ConformerSearchDialog.accept() must not raise when target_mol is None."""
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.current_mol = None
            # QDialog base is mocked — create a real instance normally
            dialog = mod.ConformerSearchDialog(ctx)
            dialog.accept()  # super().accept() calls mocked QDialog.accept → OK

    def test_run_plugin_registers_dialog_window(self):
        """run_plugin() when a mol is present creates and registers the conformer dialog."""
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            ctx.current_mol = MagicMock()
            mod.run_plugin(ctx)
            ctx.register_window.assert_called_once()
