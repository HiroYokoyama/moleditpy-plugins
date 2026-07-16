"""
Headless GUI tests for the Atom Colorizer plugin.

Covers: AtomColorizerWindow.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ATOM_COLORIZER_PATH = PLUGINS_DIR / "Atom_Colorizer" / "atom_colorizer.py"

with mock_chemistry_imports():
    _atom_colorizer = load_plugin_for_gui(ATOM_COLORIZER_PATH)


# ===========================================================================
# AtomColorizerWindow  (visible plugin: "Atom Colorizer")
# ===========================================================================


def _colorizer_context() -> MagicMock:
    """Minimal stub: get_main_window() returns None so no QObject parent issues."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


class TestAtomColorizerWindow:
    """AtomColorizerWindow with a no-main-window context."""

    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "Atom Colorizer"

    def test_indices_placeholder(self, win):
        assert win.le_indices.placeholderText() == "e.g. 0, 1, 5"

    def test_choose_color_button_exists(self, win):
        assert win.btn_color.text() == "Choose Color"

    def test_timer_is_active(self, win):
        assert win.sel_timer.isActive()

    def test_is_non_modal(self, win):
        assert not win.isModal()
