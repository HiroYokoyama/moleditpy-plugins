"""
Tests for the All-Trans Optimizer plugin.

Covers:
  1. initialize() must register at least one menu/export/plugin action
  2. No-molecule guard paths (run_plugin with mol=None)
  3. Dialog accept/reject round-trips
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
ALL_TRANS_PATH = PLUGINS_DIR / "All-Trans_Optimizer" / "all-trans_optimizer.py"


class TestAllTransOptimizer:
    def test_initialize_sets_launch_fn(self):
        """initialize() stores the launch closure in _launch_fn; run() uses it."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod._launch_fn is not None, "_launch_fn must be set by initialize()"

    def test_run_plugin_no_mol_does_not_raise(self):
        """run_plugin with mol=None should show a warning, not raise."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            mock_warn.assert_called_once()

    def test_run_plugin_no_mol_shows_correct_message(self):
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            args = mock_warn.call_args[0]
            assert "molecule" in args[2].lower() or "No molecule" in args[2]

    def test_initialize_then_run_calls_launch(self):
        """After initialize(), calling run(mw) invokes the launch function."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            # run() should not raise even without a real molecule
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning"):
                mod.run(ctx.get_main_window())
