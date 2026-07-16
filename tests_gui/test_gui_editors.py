"""
Headless GUI tests for the table-based editor plugins.

Plugins covered (all registry-visible):
  - Bond Editor   → BondEditorWindow, _ClickFilter, bond-type label maps
  - Charge Editor → ChargeEditorWindow

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_editors.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

BOND_EDITOR_PATH = PLUGINS_DIR / "Bond_Editor" / "bond_editor.py"
CHARGE_EDITOR_PATH = PLUGINS_DIR / "Charge_Editor" / "charge_editor.py"

with mock_chemistry_imports():
    _bond = load_plugin_for_gui(BOND_EDITOR_PATH)
    _charge = load_plugin_for_gui(CHARGE_EDITOR_PATH)


# ---------------------------------------------------------------------------
# Context / fake-molecule helpers
# ---------------------------------------------------------------------------


def _ctx_no_mol() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.plotter = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


class _FakeAtom:
    def __init__(self, idx, symbol="C", charge=0, radicals=0, hs=0):
        self._idx = idx
        self._symbol = symbol
        self._charge = charge
        self._radicals = radicals
        self._hs = hs

    def GetIdx(self):
        return self._idx

    def GetSymbol(self):
        return self._symbol

    def HasProp(self, name):
        return False

    def GetFormalCharge(self):
        return self._charge

    def GetNumRadicalElectrons(self):
        return self._radicals

    def GetTotalNumHs(self):
        return self._hs


class _FakeBond:
    def __init__(self, idx, begin, end, bond_type="SINGLE"):
        self._idx = idx
        self._begin = begin
        self._end = end
        self._type = bond_type

    def GetIdx(self):
        return self._idx

    def GetBeginAtomIdx(self):
        return self._begin

    def GetEndAtomIdx(self):
        return self._end

    def GetBondType(self):
        return self._type


class _FakeMol:
    """No conformer — keeps the plugins away from mocked numpy."""

    def __init__(self, atoms, bonds=()):
        self._atoms = list(atoms)
        self._bonds = list(bonds)

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetNumBonds(self):
        return len(self._bonds)

    def GetNumConformers(self):
        return 0

    def GetAtoms(self):
        return list(self._atoms)

    def GetBonds(self):
        return list(self._bonds)

    def GetAtomWithIdx(self, idx):
        return self._atoms[idx]


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
