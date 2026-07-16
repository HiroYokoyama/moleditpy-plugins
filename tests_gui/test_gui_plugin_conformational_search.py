"""
Headless GUI tests for the Conformational Search plugin.

Covers: ConformerSearchDialog.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_conformational_search.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CONF_SEARCH_PATH = PLUGINS_DIR / "Conformational_Search" / "conf_search.py"

with mock_chemistry_imports():
    _conf_search = load_plugin_for_gui(CONF_SEARCH_PATH)


def _ctx_no_mol() -> MagicMock:
    """Context with no main window and no active molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# ConformerSearchDialog  (visible plugin: "Conformational Search")
# ===========================================================================


class TestConformerSearchDialog:
    """ConformerSearchDialog with no molecule loaded."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        d = _conf_search.ConformerSearchDialog(context=ctx, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Conformational Search & Preview"

    def test_table_has_two_columns(self, dlg):
        assert dlg.table.columnCount() == 2

    def test_table_column_headers(self, dlg):
        assert dlg.table.horizontalHeaderItem(0).text() == "Rank"
        assert dlg.table.horizontalHeaderItem(1).text() == "Energy (kcal/mol)"

    def test_table_initially_empty(self, dlg):
        assert dlg.table.rowCount() == 0

    def test_force_field_combo_has_options(self, dlg):
        texts = [dlg.combo_ff.itemText(i) for i in range(dlg.combo_ff.count())]
        assert "MMFF94" in texts
        assert "UFF" in texts

    def test_show_all_checkbox_unchecked_by_default(self, dlg):
        assert not dlg.cb_show_all.isChecked()

    def test_run_button_exists(self, dlg):
        assert dlg.btn_run.text() == "Run Search"

    def test_close_button_exists(self, dlg):
        assert dlg.btn_close.text() == "Close"

    def test_conformer_data_initially_empty(self, dlg):
        assert dlg.conformer_data == []
