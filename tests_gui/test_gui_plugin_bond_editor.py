"""
Headless GUI tests for the Bond Editor plugin.

Covers: BondEditorWindow, _ClickFilter, bond-type label maps, initialize().

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_bond_editor.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports
from gui_test_helpers import FakeAtom as _FakeAtom
from gui_test_helpers import FakeBond as _FakeBond
from gui_test_helpers import FakeMol as _FakeMol
from gui_test_helpers import ctx_no_mol as _ctx_no_mol

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

BOND_EDITOR_PATH = PLUGINS_DIR / "Bond_Editor" / "bond_editor.py"

with mock_chemistry_imports():
    _bond = load_plugin_for_gui(BOND_EDITOR_PATH)


def _ethane_like():
    return _FakeMol(
        atoms=[_FakeAtom(0), _FakeAtom(1), _FakeAtom(2, symbol="O")],
        bonds=[_FakeBond(0, 0, 1, "SINGLE"), _FakeBond(1, 1, 2, "DOUBLE")],
    )


# ===========================================================================
# Bond Editor — bond-type label maps (module functions)
# ===========================================================================


class TestBondTypeMaps:
    def test_each_label_maps_to_distinct_bond_type(self):
        types = {lbl: _bond.bond_type_from_label(lbl) for lbl in _bond.BOND_TYPE_LABELS}
        assert types["Single"] is _bond.Chem.BondType.SINGLE
        assert types["Double"] is _bond.Chem.BondType.DOUBLE
        assert types["Triple"] is _bond.Chem.BondType.TRIPLE
        assert types["Aromatic"] is _bond.Chem.BondType.AROMATIC

    def test_unknown_label_defaults_to_single(self):
        assert _bond.bond_type_from_label("Quadruple") is _bond.Chem.BondType.SINGLE

    @pytest.mark.parametrize(
        "name,label",
        [("SINGLE", "Single"), ("DOUBLE", "Double"),
         ("TRIPLE", "Triple"), ("AROMATIC", "Aromatic")],
    )
    def test_label_from_bond_type_names(self, name, label):
        assert _bond.label_from_bond_type(name) == label

    def test_label_from_qualified_rdkit_repr(self):
        assert _bond.label_from_bond_type("Chem.BondType.DOUBLE") == "Double"

    def test_label_from_unknown_type_is_single(self):
        assert _bond.label_from_bond_type("DATIVE") == "Single"

    def test_bond_type_labels_constant(self):
        assert _bond.BOND_TYPE_LABELS == ["Single", "Double", "Triple", "Aromatic"]


# ===========================================================================
# Bond Editor — _ClickFilter with real QMouseEvents
# ===========================================================================


class TestClickFilter:
    def _events(self, press_xy, release_xy):
        from PyQt6.QtCore import QEvent, QPointF, Qt
        from PyQt6.QtGui import QMouseEvent

        press = QMouseEvent(
            QEvent.Type.MouseButtonPress,
            QPointF(*press_xy),
            QPointF(*press_xy),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
        release = QMouseEvent(
            QEvent.Type.MouseButtonRelease,
            QPointF(*release_xy),
            QPointF(*release_xy),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
        )
        return press, release

    def test_click_within_5px_invokes_callback(self, qapp):
        calls = []
        f = _bond._ClickFilter(lambda x, y, obj, mods: calls.append((x, y)))
        press, release = self._events((10, 10), (12, 13))
        obj = MagicMock()
        f.eventFilter(obj, press)
        f.eventFilter(obj, release)
        assert calls == [(12, 13)]

    def test_drag_beyond_5px_is_ignored(self, qapp):
        calls = []
        f = _bond._ClickFilter(lambda x, y, obj, mods: calls.append((x, y)))
        press, release = self._events((10, 10), (30, 30))
        obj = MagicMock()
        f.eventFilter(obj, press)
        f.eventFilter(obj, release)
        assert calls == []

    def test_never_consumes_events(self, qapp):
        f = _bond._ClickFilter(lambda *a: None)
        press, release = self._events((0, 0), (0, 0))
        obj = MagicMock()
        assert f.eventFilter(obj, press) is False
        assert f.eventFilter(obj, release) is False

    def test_release_without_press_is_ignored(self, qapp):
        calls = []
        f = _bond._ClickFilter(lambda x, y, obj, mods: calls.append((x, y)))
        _, release = self._events((0, 0), (5, 5))
        f.eventFilter(MagicMock(), release)
        assert calls == []


# ===========================================================================
# BondEditorWindow — empty context
# ===========================================================================


class TestBondEditorWindowEmpty:
    @pytest.fixture
    def win(self, qapp):
        ctx = _ctx_no_mol()
        w = _bond.BondEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_window_title(self, win):
        assert win.windowTitle() == "Bond Editor"

    def test_table_columns(self, win):
        assert win.table.columnCount() == 5
        headers = [win.table.horizontalHeaderItem(i).text() for i in range(5)]
        assert headers == ["#", "Atom 1", "Atom 2", "Type", "Length (Å)"]

    def test_table_empty_without_molecule(self, win):
        assert win.table.rowCount() == 0

    def test_click_mode_combo_options(self, win):
        texts = [win.click_mode_combo.itemText(i) for i in range(win.click_mode_combo.count())]
        assert texts == ["Select bond", "Create bond"]

    def test_add_type_combo_options(self, win):
        texts = [win.add_type_combo.itemText(i) for i in range(win.add_type_combo.count())]
        assert texts == _bond.BOND_TYPE_LABELS

    def test_buttons_exist(self, win):
        assert win.delete_btn.text() == "Delete Selected Bonds"
        assert win.unselect_btn.text() == "Unselect"

    def test_update_timer_running(self, win):
        assert win.update_timer.isActive()

    def test_registers_main_panel_window(self, win):
        win.context.register_window.assert_any_call("main_panel", win)

    def test_no_click_filter_without_plotter(self, win):
        assert win._click_filter is None

    def test_row_bond_atoms_empty_table(self, win):
        assert win._row_bond_atoms(0) is None

    def test_add_bond_without_molecule_warns(self, win):
        win.add_bond(0, 1)
        msg = win.context.show_status_message.call_args[0][0]
        assert "No molecule" in msg

    def test_mode_change_resets_pick_and_reports(self, win):
        win._first_pick_idx = 3
        win.click_mode_combo.setCurrentText("Create bond")
        assert win._first_pick_idx is None
        msg = win.context.show_status_message.call_args[0][0]
        assert "two atoms" in msg

    def test_create_bond_pick_state_machine(self, win):
        win.click_mode_combo.blockSignals(True)
        win.click_mode_combo.setCurrentText("Create bond")
        win.click_mode_combo.blockSignals(False)
        win._create_bond_pick(3)
        assert win._first_pick_idx == 3
        assert win._picked_atoms == {"Atom 1": 3}
        # Picking the same atom cancels
        win._create_bond_pick(3)
        assert win._first_pick_idx is None
        assert win._picked_atoms == {}

    def test_second_pick_attempts_add_bond(self, win):
        win._create_bond_pick(0)
        win._create_bond_pick(1)  # no molecule → status message, pick cleared
        assert win._first_pick_idx is None
        msg = win.context.show_status_message.call_args[0][0]
        assert "No molecule" in msg

    def test_unselect_all_without_pick_is_noop(self, win):
        win.unselect_all()
        assert win._first_pick_idx is None

    def test_close_stops_timer_and_unregisters(self, qapp):
        ctx = _ctx_no_mol()
        w = _bond.BondEditorWindow(context=ctx)
        w.close()
        assert not w.update_timer.isActive()
        ctx.register_window.assert_called_with("main_panel", None)
        w.destroy()


# ===========================================================================
# BondEditorWindow — populated from a fake molecule
# ===========================================================================


class TestBondEditorWindowWithMol:
    @pytest.fixture
    def win(self, qapp):
        ctx = _ctx_no_mol()
        mol = _ethane_like()
        ctx.current_mol = mol
        ctx.current_molecule = mol
        w = _bond.BondEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_one_row_per_bond(self, win):
        assert win.table.rowCount() == 2

    def test_atom_labels_in_rows(self, win):
        assert win.table.item(0, win.COL_A1).text() == "0 (C)"
        assert win.table.item(1, win.COL_A2).text() == "2 (O)"

    def test_type_combo_reflects_bond_type(self, win):
        assert win.table.cellWidget(0, win.COL_TYPE).currentText() == "Single"
        assert win.table.cellWidget(1, win.COL_TYPE).currentText() == "Double"

    def test_length_na_without_conformer(self, win):
        assert win.table.item(0, win.COL_LEN).text() == "n/a"

    def test_row_bond_atoms_parses_labels(self, win):
        assert win._row_bond_atoms(0) == (0, 1)
        assert win._row_bond_atoms(1) == (1, 2)

    def test_signature_captured(self, win):
        assert win.last_seen_signature is not None
        assert win.last_seen_signature[1] == 3  # atoms
        assert win.last_seen_signature[2] == 2  # bonds

    def test_delete_with_no_selection_warns(self, win):
        win.context.show_status_message.reset_mock()
        win.delete_selected_bonds()
        msg = win.context.show_status_message.call_args[0][0]
        assert "No bonds selected" in msg

    def test_signature_changes_when_molecule_changes(self, win):
        sig_before = win.last_seen_signature
        bigger = _FakeMol(
            atoms=[_FakeAtom(i) for i in range(4)],
            bonds=[_FakeBond(0, 0, 1), _FakeBond(1, 1, 2), _FakeBond(2, 2, 3)],
        )
        assert win.get_mol_signature(bigger) != sig_before


# ===========================================================================
# Bond Editor — initialize()
# ===========================================================================


class TestBondEditorInitialize:
    def test_registers_menu_and_reset_handler(self, qapp):
        ctx = MagicMock()
        _bond.initialize(ctx)
        path = ctx.add_menu_action.call_args[0][0]
        assert path == "3D Edit/Bond Editor..."
        ctx.register_document_reset_handler.assert_called_once()

    def test_reset_handler_reloads_open_window(self, qapp):
        ctx = MagicMock()
        _bond.initialize(ctx)
        handler = ctx.register_document_reset_handler.call_args[0][0]
        fake_win = MagicMock()
        ctx.get_window.return_value = fake_win
        handler()
        fake_win.load_molecule.assert_called_once()
