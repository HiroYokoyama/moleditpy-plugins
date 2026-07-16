"""
Headless GUI tests for the Charge Editor plugin.

Covers: ChargeEditorWindow, initialize().

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_charge_editor.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports
from gui_test_helpers import FakeAtom as _FakeAtom
from gui_test_helpers import FakeMol as _FakeMol
from gui_test_helpers import ctx_no_mol as _ctx_no_mol

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CHARGE_EDITOR_PATH = PLUGINS_DIR / "Charge_Editor" / "charge_editor.py"

with mock_chemistry_imports():
    _charge = load_plugin_for_gui(CHARGE_EDITOR_PATH)


# ===========================================================================
# ChargeEditorWindow — empty context
# ===========================================================================


class TestChargeEditorWindowEmpty:
    @pytest.fixture
    def win(self, qapp):
        ctx = _ctx_no_mol()
        w = _charge.ChargeEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_window_title(self, win):
        assert win.windowTitle() == "Charge Editor"

    def test_table_columns(self, win):
        assert win.table.columnCount() == 5
        headers = [win.table.horizontalHeaderItem(i).text() for i in range(5)]
        assert headers == ["#", "Atom", "Formal Charge", "Radical e⁻", "H count"]

    def test_summary_default_text(self, win):
        assert win.summary_label.text() == "Total charge: 0    Multiplicity: 1"

    def test_buttons_exist(self, win):
        assert win.clear_charges_btn.text() == "Clear All Charges"
        assert win.clear_radicals_btn.text() == "Clear All Radicals"
        assert win.unselect_btn.text() == "Unselect"

    def test_update_timer_running(self, win):
        assert win.update_timer.isActive()

    def test_table_empty_without_molecule(self, win):
        assert win.table.rowCount() == 0

    def test_row_atom_index_empty_table(self, win):
        assert win._row_atom_index(0) is None

    def test_clear_charges_without_molecule_warns(self, win):
        win.clear_all_charges()
        msg = win.context.show_status_message.call_args[0][0]
        assert "No molecule" in msg

    def test_clear_radicals_without_molecule_warns(self, win):
        win.clear_all_radicals()
        msg = win.context.show_status_message.call_args[0][0]
        assert "No molecule" in msg

    def test_charge_change_without_molecule_is_noop(self, win):
        win.on_charge_changed(0, 1)  # must not raise
        win.context.push_undo_checkpoint.assert_not_called()

    def test_close_stops_timer_and_unregisters(self, qapp):
        ctx = _ctx_no_mol()
        w = _charge.ChargeEditorWindow(context=ctx)
        w.close()
        assert not w.update_timer.isActive()
        ctx.register_window.assert_called_with("main_panel", None)
        w.destroy()

    def test_charge_and_radical_range_constants(self):
        assert (_charge.CHARGE_MIN, _charge.CHARGE_MAX) == (-4, 4)
        assert (_charge.RADICAL_MIN, _charge.RADICAL_MAX) == (0, 4)


# ===========================================================================
# ChargeEditorWindow — populated from a fake molecule
# ===========================================================================


class TestChargeEditorWindowWithMol:
    @pytest.fixture
    def win(self, qapp):
        ctx = _ctx_no_mol()
        mol = _FakeMol(
            atoms=[
                _FakeAtom(0, "N", charge=1, radicals=0, hs=4),
                _FakeAtom(1, "O", charge=-1, radicals=0, hs=0),
                _FakeAtom(2, "C", charge=0, radicals=1, hs=3),
            ]
        )
        ctx.current_mol = mol
        ctx.current_molecule = mol
        w = _charge.ChargeEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_one_row_per_atom(self, win):
        assert win.table.rowCount() == 3

    def test_symbols_shown(self, win):
        symbols = [win.table.item(r, win.COL_SYMBOL).text() for r in range(3)]
        assert symbols == ["N", "O", "C"]

    def test_charge_spinboxes_match_atoms(self, win):
        values = [win.table.cellWidget(r, win.COL_CHARGE).value() for r in range(3)]
        assert values == [1, -1, 0]

    def test_radical_spinboxes_match_atoms(self, win):
        values = [win.table.cellWidget(r, win.COL_RADICAL).value() for r in range(3)]
        assert values == [0, 0, 1]

    def test_h_count_column(self, win):
        counts = [win.table.item(r, win.COL_HS).text() for r in range(3)]
        assert counts == ["4", "0", "3"]

    def test_summary_reports_net_charge_and_multiplicity(self, win):
        text = win.summary_label.text()
        assert "Total charge: +0" in text
        assert "Multiplicity: 2" in text
        assert "1 unpaired" in text

    def test_spinbox_ranges(self, win):
        spin = win.table.cellWidget(0, win.COL_CHARGE)
        assert (spin.minimum(), spin.maximum()) == (-4, 4)
        rad = win.table.cellWidget(0, win.COL_RADICAL)
        assert (rad.minimum(), rad.maximum()) == (0, 4)

    def test_row_atom_index_parses(self, win):
        assert win._row_atom_index(2) == 2

    def test_signature_tracks_charges(self, win):
        sig = win.get_mol_signature(win.context.current_molecule)
        assert sig == win.last_seen_signature
        other = _FakeMol(atoms=[_FakeAtom(0, "N", charge=0)])
        assert win.get_mol_signature(other) != sig


# ===========================================================================
# Charge Editor — initialize()
# ===========================================================================


class TestChargeEditorInitialize:
    def test_registers_menu_and_reset_handler(self, qapp):
        ctx = MagicMock()
        _charge.initialize(ctx)
        path = ctx.add_menu_action.call_args[0][0]
        assert path == "3D Edit/Charge Editor..."
        ctx.register_document_reset_handler.assert_called_once()

    def test_reset_handler_ignores_missing_window(self, qapp):
        ctx = MagicMock()
        _charge.initialize(ctx)
        handler = ctx.register_document_reset_handler.call_args[0][0]
        ctx.get_window.return_value = None
        handler()  # must not raise
