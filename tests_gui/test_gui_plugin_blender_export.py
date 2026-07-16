"""
Headless GUI tests for the Blender Export plugin.

Covers module constants, initialize(), and run() (no dialog class; export
logic is callable-only).

All tests use real PyQt6 (QT_QPA_PLATFORM=offscreen).
Chemistry libraries are mocked via mock_chemistry_imports().
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

BLENDER_PATH = PLUGINS_DIR / "Blender_Export" / "blender_export.py"

with mock_chemistry_imports():
    _blender = load_plugin_for_gui(BLENDER_PATH)


def _make_ctx() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# Blender Export — module-level and initialize()/run() tests
# (no dialog class; export logic is callable-only)
# ===========================================================================


class TestBlenderExportModule:
    """Module constants, initialize(), and run() for Blender Export."""

    def test_plugin_name(self):
        assert _blender.PLUGIN_NAME == "Blender Export"

    def test_plugin_version_is_string(self):
        assert isinstance(_blender.PLUGIN_VERSION, str)

    def test_plugin_context_initially_none(self):
        assert _blender.PLUGIN_CONTEXT is None

    def test_initialize_registers_export_action(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_blender(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        label = ctx.add_export_action.call_args.args[0]
        assert "Blender" in label

    def test_initialize_callback_is_callable(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        cb = ctx.add_export_action.call_args.args[1]
        assert callable(cb)

    def test_run_noop_without_plugin_manager(self):
        mw = MagicMock(spec=[])  # no plugin_manager attribute
        _blender.run(mw)  # must not raise

    def test_run_noop_when_context_is_none(self):
        _blender.PLUGIN_CONTEXT = None
        mw = MagicMock()  # has plugin_manager (MagicMock auto-attribute)
        _blender.run(mw)  # PLUGIN_CONTEXT is None → early return, no error
