"""
Headless GUI tests for the XYZ Editor plugin.

Covers: XYZEditorWindow.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_xyz_editor.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

XYZ_EDITOR_PATH = PLUGINS_DIR / "XYZ_Editor" / "xyz_editor.py"

with mock_chemistry_imports():
    _xyz_editor = load_plugin_for_gui(XYZ_EDITOR_PATH)


def _ctx_no_mol() -> MagicMock:
    """Context with no main window and no active molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# XYZEditorWindow  (visible plugin: "XYZ Editor")
# ===========================================================================


class TestXYZEditorWindow:
    """XYZEditorWindow with no main window and no molecule."""

    @pytest.fixture
    def win(self, qapp):
        ctx = _ctx_no_mol()
        w = _xyz_editor.XYZEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "XYZ Editor"

    def test_table_has_five_columns(self, win):
        assert win.table.columnCount() == 5

    def test_table_column_headers(self, win):
        headers = [
            win.table.horizontalHeaderItem(i).text() for i in range(5)
        ]
        assert headers == ["Index", "Symbol", "X", "Y", "Z"]

    def test_table_empty_with_no_molecule(self, win):
        assert win.table.rowCount() == 0

    def test_apply_button_exists(self, win):
        assert win.apply_btn.text() == "Apply to View"

    def test_save_button_exists(self, win):
        assert win.save_btn.text() == "Save as XYZ..."

    def test_add_atom_button_exists(self, win):
        assert win.add_btn.text() == "Add Atom"

    def test_timer_is_active(self, win):
        assert win.update_timer.isActive()

    def test_resize_to_expected_size(self, win):
        assert win.width() >= 400
        assert win.height() >= 300
