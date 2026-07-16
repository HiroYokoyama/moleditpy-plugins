"""
Headless GUI tests for the Paste XYZ plugin.

Covers: PasteXYZDialog.

These tests use real PyQt6 (QT_QPA_PLATFORM=offscreen) rather than mocking
all Qt classes.  Chemistry/scientific libraries (rdkit, numpy, …) are still
replaced with MagicMock via mock_chemistry_imports() so no installed chemistry
stack is required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PASTE_XYZ_PATH = PLUGINS_DIR / "Paste_XYZ" / "paste_xyz.py"

with mock_chemistry_imports():
    _paste_xyz = load_plugin_for_gui(PASTE_XYZ_PATH)


# ===========================================================================
# PasteXYZDialog  (visible plugin: "Paste XYZ")
# ===========================================================================


class TestPasteXYZDialog:
    """PasteXYZDialog — pure-Qt dialog; rdkit guarded by try/except."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _paste_xyz.PasteXYZDialog(parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Paste XYZ"

    def test_default_size_wide_enough(self, dlg):
        assert dlg.width() >= 400

    def test_get_data_initially_empty(self, dlg):
        assert dlg.get_data() == ""

    def test_get_data_returns_typed_text(self, dlg):
        dlg.text_edit.setPlainText("C 0.0 0.0 0.0\nH 1.0 0.0 0.0")
        result = dlg.get_data()
        assert "C 0.0 0.0 0.0" in result
        assert "H 1.0 0.0 0.0" in result

    def test_multiline_xyz_round_trip(self, dlg):
        xyz = "C 0.0 0.0 0.0\nH 1.089 0.0 0.0\nH -0.363 1.027 0.0"
        dlg.text_edit.setPlainText(xyz)
        assert dlg.get_data().strip() == xyz

    def test_load_button_label(self, dlg):
        assert dlg.load_btn.text() == "Load"

    def test_cancel_button_label(self, dlg):
        assert dlg.cancel_btn.text() == "Cancel"

    def test_text_edit_placeholder_mentions_xyz(self, dlg):
        ph = dlg.text_edit.placeholderText().lower()
        assert "xyz" in ph

    def test_clear_resets_get_data(self, dlg):
        dlg.text_edit.setPlainText("something")
        dlg.text_edit.clear()
        assert dlg.get_data() == ""
