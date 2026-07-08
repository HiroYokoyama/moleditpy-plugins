"""
Tests for the PubChem Name Resolver plugin.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
PUBCHEM_PATH = PLUGINS_DIR / "PubChem_Name_Ressolver" / "pubchem_ressolver.py"


class TestPubChemNameResolver:
    def test_initialize_does_not_raise(self):
        """initialize() defines run_resolver locally; run(mw) handles auto-registration."""
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
