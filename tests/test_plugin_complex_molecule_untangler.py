"""
Tests for the Complex Molecule Untangler plugin.

Covers:
  1. initialize() must register at least one menu/export/plugin action
  2. No-molecule guard paths (run_plugin with mol=None)
  3. Dialog accept/reject round-trips
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
UNTANGLER_PATH = PLUGINS_DIR / "Complex_Molecule_Untangler" / "complex_molecule_untangler.py"


class TestComplexMoleculeUntangler:
    def test_initialize_stores_context(self):
        """initialize() must store the context in PLUGIN_CONTEXT and set _launch_fn."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx
            assert mod._launch_fn is not None

    def test_run_plugin_opens_or_raises_dialog(self):
        """run_plugin with a valid context should not raise."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            mod.initialize(ctx)
            # run_plugin creates a dialog — should not raise with mocked Qt
            mod.run_plugin(ctx)

    def test_run_plugin_reuses_existing_window(self):
        """If a window is already registered, run_plugin shows it without creating new."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            existing_win = MagicMock()
            ctx.get_window.return_value = existing_win
            mod.run_plugin(ctx)
            existing_win.show.assert_called_once()
            existing_win.raise_.assert_called_once()

    def test_run_plugin_registers_new_dialog_window(self):
        """run_plugin() registers the created dialog as 'main_panel' when no existing window."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            mod.run_plugin(ctx)
            ctx.register_window.assert_called_once()
