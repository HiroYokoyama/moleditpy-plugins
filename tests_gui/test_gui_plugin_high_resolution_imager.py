"""
Headless GUI tests for the High Resolution Imager plugin.

Covers module constants and initialize() (no dialog class — the dialog is
built inline in take_screenshot()).

All tests use real PyQt6 (QT_QPA_PLATFORM=offscreen).
Chemistry libraries are mocked via mock_chemistry_imports().
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

HI_RES_PATH = PLUGINS_DIR / "High_Resolution_Imager" / "high_res_imager.py"

with mock_chemistry_imports():
    _hi_res = load_plugin_for_gui(HI_RES_PATH)


def _make_ctx() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# High Resolution Imager — module-level and initialize() tests
# (dialog is built inline in take_screenshot(); not a separate class)
# ===========================================================================


class TestHighResImagerModule:
    """Module constants and initialize() for High Resolution Imager."""

    def test_plugin_name(self):
        assert _hi_res.PLUGIN_NAME == "High Resolution Imager"

    def test_plugin_version_is_date_string(self):
        assert isinstance(_hi_res.PLUGIN_VERSION, str)
        assert len(_hi_res.PLUGIN_VERSION) == 10  # YYYY.MM.DD

    def test_plugin_context_initially_none(self):
        assert _hi_res.PLUGIN_CONTEXT is None

    def test_initialize_registers_export_action(self):
        ctx = _make_ctx()
        _hi_res.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_screenshot(self):
        ctx = _make_ctx()
        _hi_res.initialize(ctx)
        label = ctx.add_export_action.call_args.args[0]
        assert "Screenshot" in label

    def test_initialize_callback_is_callable(self):
        ctx = _make_ctx()
        _hi_res.initialize(ctx)
        cb = ctx.add_export_action.call_args.args[1]
        assert callable(cb)

    def test_initialize_stores_context(self):
        ctx = _make_ctx()
        _hi_res.PLUGIN_CONTEXT = None
        _hi_res.initialize(ctx)
        assert _hi_res.PLUGIN_CONTEXT is ctx
        _hi_res.PLUGIN_CONTEXT = None  # reset
