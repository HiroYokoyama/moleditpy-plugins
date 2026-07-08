"""
Tests for the Conformational Search plugin.

Covers:
  1. initialize() must register at least one menu/export/plugin action
  2. No-molecule guard paths (run_plugin with mol=None)
  3. Dialog accept/reject round-trips
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CONF_SEARCH_PATH = PLUGINS_DIR / "Conformational_Search" / "conf_search.py"


def _menu_registered(ctx: MagicMock) -> bool:
    """Return True if initialize() called any recognised registration method."""
    return (
        ctx.add_menu_action.called
        or ctx.add_export_action.called
        or ctx.add_plugin_menu.called
        or ctx.add_analysis_tool.called
        or ctx.add_toolbar_action.called
    )


class TestConformationalSearch:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert _menu_registered(ctx), (
                "Conformational Search initialize() must call add_menu_action()"
            )

    def test_initialize_menu_path_contains_conformational(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            call_args = ctx.add_menu_action.call_args
            assert call_args is not None
            path = call_args[0][0]
            assert "Conformational" in path or "conformational" in path.lower()

    def test_run_plugin_no_mol_warns(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            mock_warn.assert_called_once()

    def test_dialog_accept_does_not_raise(self):
        """ConformerSearchDialog.accept() must not raise when target_mol is None."""
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.current_mol = None
            # QDialog base is mocked — create a real instance normally
            dialog = mod.ConformerSearchDialog(ctx)
            dialog.accept()  # super().accept() calls mocked QDialog.accept → OK

    def test_run_plugin_registers_dialog_window(self):
        """run_plugin() when a mol is present creates and registers the conformer dialog."""
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            ctx.current_mol = MagicMock()
            mod.run_plugin(ctx)
            ctx.register_window.assert_called_once()


def _conf_filter_fn():
    return extract_function(
        CONF_SEARCH_PATH, "ConformerSearchDialog", "apply_filter_and_update", {}
    )


def _conf_self(results_raw, show_all=False):
    s = SimpleNamespace()
    s.results_raw = results_raw
    s.conformer_data = []
    s.cb_show_all = MagicMock()
    s.cb_show_all.isChecked.return_value = show_all
    s.update_table = MagicMock()
    s.lbl_info = MagicMock()
    return s


class TestConfSearchFilter:
    def test_empty_results_returns_early(self):
        fn = _conf_filter_fn()
        s = _conf_self([])
        fn(s)
        s.update_table.assert_not_called()
        assert s.conformer_data == []

    def test_show_all_keeps_everything(self):
        fn = _conf_filter_fn()
        raw = [(1.0, 0), (1.0, 1), (2.0, 2)]
        s = _conf_self(raw, show_all=True)
        fn(s)
        assert s.conformer_data == raw
        s.update_table.assert_called_once()

    def test_duplicate_energies_deduplicated(self):
        fn = _conf_filter_fn()
        # 1.0 and 1.00005 are within the 1e-4 window -> one survivor
        raw = [(1.0, 0), (1.00005, 1), (2.0, 2)]
        s = _conf_self(raw)
        fn(s)
        assert s.conformer_data == [(1.0, 0), (2.0, 2)]

    def test_distinct_energies_all_kept(self):
        fn = _conf_filter_fn()
        raw = [(1.0, 0), (1.5, 1), (2.0, 2)]
        s = _conf_self(raw)
        fn(s)
        assert s.conformer_data == raw

    def test_info_label_shows_counts(self):
        fn = _conf_filter_fn()
        raw = [(1.0, 0), (1.00001, 1), (3.0, 2)]
        s = _conf_self(raw)
        fn(s)
        msg = s.lbl_info.setText.call_args[0][0]
        assert "Showing 2 conformers" in msg
        assert "Total found: 3" in msg


class _FakeTable:
    """Recorder stand-in for the installer's QTableWidget."""

    def __init__(self):
        self.rows = 0
        self.items = {}  # (row, col) -> _Item
        self.cell_widgets = {}  # (row, col) -> widget
        self.hidden = {}  # row -> bool

    def setRowCount(self, n):
        self.rows = n
        if n == 0:
            self.items.clear()
            self.cell_widgets.clear()

    def rowCount(self):
        return self.rows

    def insertRow(self, row):
        self.rows += 1

    def setItem(self, row, col, item):
        self.items[(row, col)] = item

    def item(self, row, col):
        return self.items.get((row, col))

    def setCellWidget(self, row, col, widget):
        self.cell_widgets[(row, col)] = widget

    def setUpdatesEnabled(self, flag):
        pass

    def setRowHidden(self, row, hidden):
        self.hidden[row] = hidden


class TestConfSearchUpdateTable:
    def _run(self, conformer_data):
        items = []

        class _RecItem:
            def __init__(self, text):
                self.text = text
                self.user_data = None
                items.append(self)

            def setData(self, role, value):
                self.user_data = value

        globs = {"QTableWidgetItem": _RecItem, "Qt": MagicMock()}
        fn = extract_function(
            CONF_SEARCH_PATH, "ConformerSearchDialog", "update_table", globs
        )
        s = SimpleNamespace()
        s.conformer_data = conformer_data
        s.table = _FakeTable()
        fn(s)
        return s, items

    def test_rows_ranked_and_energy_formatted(self):
        s, _ = self._run([(1.23456789, 7), (2.5, 3)])
        assert s.table.rows == 2
        assert s.table.items[(0, 0)].text == "1"
        assert s.table.items[(0, 1)].text == "1.2346"  # 4 decimal places
        assert s.table.items[(1, 0)].text == "2"
        assert s.table.items[(1, 1)].text == "2.5000"

    def test_conformer_id_stored_as_user_data(self):
        s, _ = self._run([(1.0, 42)])
        assert s.table.items[(0, 0)].user_data == 42
