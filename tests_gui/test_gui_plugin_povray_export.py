"""
Headless GUI tests for the POV-Ray Export plugin.

Covers module constants, initialize(), and run().

All tests use real PyQt6 (QT_QPA_PLATFORM=offscreen).
Chemistry libraries are mocked via mock_chemistry_imports().
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

POVRAY_PATH = PLUGINS_DIR / "POV-Ray_Export" / "povray_export.py"

with mock_chemistry_imports():
    _povray = load_plugin_for_gui(POVRAY_PATH)


def _make_ctx() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# POV-Ray Export — module-level and initialize()/run() tests
# ===========================================================================


class TestPovRayExportModule:
    """Module constants, initialize(), and run() for POV-Ray Export."""

    def test_plugin_name(self):
        assert _povray.PLUGIN_NAME == "POV-Ray Export"

    def test_plugin_version_is_string(self):
        assert isinstance(_povray.PLUGIN_VERSION, str)

    def test_plugin_context_initially_none(self):
        assert _povray.PLUGIN_CONTEXT is None

    def test_initialize_registers_export_action(self):
        ctx = _make_ctx()
        _povray.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_povray(self):
        ctx = _make_ctx()
        _povray.initialize(ctx)
        label = ctx.add_export_action.call_args.args[0]
        assert "POV" in label

    def test_initialize_callback_is_callable(self):
        ctx = _make_ctx()
        _povray.initialize(ctx)
        cb = ctx.add_export_action.call_args.args[1]
        assert callable(cb)

    def test_run_noop_without_plugin_manager(self):
        mw = MagicMock(spec=[])
        _povray.run(mw)  # must not raise

    def test_run_noop_when_context_is_none(self):
        _povray.PLUGIN_CONTEXT = None
        mw = MagicMock()
        _povray.run(mw)  # PLUGIN_CONTEXT is None → early return
