"""
Tests for the Blender Export plugin (plugins/Blender_Export/blender_export.py).

All heavy deps (PyQt6, rdkit, numpy) are stubbed via mock_optional_imports().
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
BLENDER_PATH = PLUGINS_DIR / "Blender_Export" / "blender_export.py"

with mock_optional_imports():
    _blender = load_plugin(BLENDER_PATH)


class TestBlenderExportInitialize:
    def test_initialize_calls_add_export_action(self):
        with mock_optional_imports():
            mod = load_plugin(BLENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_export_action_label_contains_blender(self):
        with mock_optional_imports():
            mod = load_plugin(BLENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        label = ctx.add_export_action.call_args[0][0]
        assert "Blender" in label

    def test_export_action_callback_is_callable(self):
        with mock_optional_imports():
            mod = load_plugin(BLENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        callback = ctx.add_export_action.call_args[0][1]
        assert callable(callback)

    def test_plugin_version_constant_present(self):
        assert hasattr(_blender, "PLUGIN_VERSION")
        assert _blender.PLUGIN_VERSION
