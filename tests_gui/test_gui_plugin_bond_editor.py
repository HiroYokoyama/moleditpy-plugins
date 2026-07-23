"""
Headless GUI tests for the Bond Editor plugin.

Covers: BondEditorWindow, _ClickFilter, bond-type label maps, initialize().

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_bond_editor.py
"""

from __future__ import annotations

import contextlib
import sys
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


# Import the real packages up front so they are already cached in
# sys.modules before _mock_chemistry_keep_real_chem() snapshots it below —
# otherwise the chemistry-blocking meta path finder would intercept the
# first-ever `import rdkit`/`import vtk` and hand back a MagicMock instead.
import numpy as _np
import pyvista as _pv
import vtk as _vtk
from rdkit import Chem as _Chem
from rdkit.Geometry import Point3D as _Point3D


@contextlib.contextmanager
def _mock_chemistry_keep_real_chem():
    """Like mock_chemistry_imports(), but numpy/rdkit/pyvista/vtk resolve to the
    real packages (all installed in this environment).

    Needed so the plugin's real math (np.linalg.norm, RWMol sanitize, PolyData
    tube generation, vtkCellPicker) actually runs instead of chasing MagicMock
    attribute chains.
    """
    keep_prefixes = ("numpy", "rdkit", "pyvista", "pyvistaqt", "vtk", "vtkmodules")
    real_mods = {
        k: v
        for k, v in sys.modules.items()
        if k.split(".")[0] in keep_prefixes
    }
    with mock_chemistry_imports():
        sys.modules.update(real_mods)
        yield


# Second module instance with real chemistry libs, used for tests that drive
# real conformer/vector math and RDKit sanitize/undo behavior through the
# plugin's bound methods (separate from `_bond` above, which keeps everything
# mocked for the plain widget-construction tests).
with _mock_chemistry_keep_real_chem():
    _bondrn = load_plugin_for_gui(BOND_EDITOR_PATH)


def _real_mol(ring=False):
    """Real RDKit RWMol with a conformer: C0-C1(single)-C2=O(double), or a
    4-membered all-carbon ring when ring=True."""
    rw = _Chem.RWMol()
    if ring:
        for _ in range(4):
            rw.AddAtom(_Chem.Atom(6))
        rw.AddBond(0, 1, _Chem.BondType.SINGLE)
        rw.AddBond(1, 2, _Chem.BondType.SINGLE)
        rw.AddBond(2, 3, _Chem.BondType.SINGLE)
        rw.AddBond(3, 0, _Chem.BondType.SINGLE)
        n = 4
        coords = [(0, 0, 0), (1.5, 0, 0), (1.5, 1.5, 0), (0, 1.5, 0)]
    else:
        rw.AddAtom(_Chem.Atom(6))
        rw.AddAtom(_Chem.Atom(6))
        rw.AddAtom(_Chem.Atom(8))
        rw.AddBond(0, 1, _Chem.BondType.SINGLE)
        rw.AddBond(1, 2, _Chem.BondType.DOUBLE)
        n = 3
        coords = [(0, 0, 0), (1.5, 0, 0), (1.5, 1.2, 0)]
    conf = _Chem.Conformer(n)
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, _Point3D(x, y, z))
    rw.AddConformer(conf, assignId=True)
    mol = rw.GetMol()
    _Chem.SanitizeMol(mol)
    return mol


def _real_ctx(mol=None):
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.plotter = None
    ctx.current_mol = mol
    ctx.current_molecule = mol
    return ctx


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

    def test_show_editor_creates_new_window(self, qapp):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        ctx.get_window.return_value = None
        _bond.initialize(ctx)
        show_editor = ctx.add_menu_action.call_args[0][1]
        show_editor()
        win = ctx.register_window.call_args[0][1]
        assert isinstance(win, _bond.BondEditorWindow)
        win.close()
        win.destroy()

    def test_show_editor_reuses_existing_window(self, qapp):
        ctx = MagicMock()
        _bond.initialize(ctx)
        fake_win = MagicMock()
        ctx.get_window.return_value = fake_win
        show_editor = ctx.add_menu_action.call_args[0][1]
        show_editor()
        fake_win.show.assert_called_once()
        fake_win.raise_.assert_called_once()
        fake_win.activateWindow.assert_called_once()
        fake_win.load_molecule.assert_called_once()


# ===========================================================================
# Real RDKit / numpy / pyvista / vtk — conformer-dependent code paths
# ===========================================================================


class TestRealMolLoadAndSignature:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _bondrn.BondEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_bond_lengths_computed_from_conformer(self, win):
        assert win.table.item(0, win.COL_LEN).text() == "1.5000"
        # sqrt(0^2 + 1.2^2 + 0^2)
        assert win.table.item(1, win.COL_LEN).text() == "1.2000"

    def test_signature_includes_position_hash(self, win):
        sig = win.get_mol_signature(win.context.current_molecule)
        assert len(sig) == 5  # id, natoms, nbonds, bond_sig_hash, pos_hash

    def test_get_mol_signature_returns_none_on_error(self, win):
        broken = object()
        assert win.get_mol_signature(broken) is None

    def test_atom_label_uses_custom_symbol(self, win):
        mol = win.context.current_molecule
        mol.GetAtomWithIdx(0).SetProp("custom_symbol", "Cx")
        assert win._atom_label(mol, 0) == "0 (Cx)"

    def test_check_molecule_update_reloads_on_change(self, win):
        new_mol = _real_mol(ring=True)
        win.context.current_molecule = new_mol
        win.check_molecule_update()
        assert win.table.rowCount() == 4
        assert win.last_seen_signature == win.get_mol_signature(new_mol)

    def test_check_molecule_update_noop_when_unchanged(self, win):
        win.check_molecule_update()
        assert win.table.rowCount() == 2


class TestPlotterPickingRealInteractor:
    def test_enable_installs_event_filter_on_real_interactor(self, qapp):
        interactor = MagicMock()
        plotter = MagicMock()
        plotter.interactor = interactor
        ctx = _real_ctx(mol=None)
        ctx.plotter = plotter
        w = _bondrn.BondEditorWindow(context=ctx)
        assert w._click_filter is not None
        filt = w._click_filter
        interactor.installEventFilter.assert_called_once_with(filt)
        w.close()
        interactor.removeEventFilter.assert_called_once_with(filt)
        assert w._click_filter is None
        w.destroy()

    def test_enable_noop_when_interactor_is_none(self, qapp):
        plotter = MagicMock()
        plotter.interactor = None
        ctx = _real_ctx(mol=None)
        ctx.plotter = plotter
        w = _bondrn.BondEditorWindow(context=ctx)
        assert w._click_filter is None
        w.close()
        w.destroy()


class _FakePicker:
    """Stand-in for vtk.vtkCellPicker with a scripted result."""

    def __init__(self, actor=None, pos=(0.0, 0.0, 0.0)):
        self._actor = actor
        self._pos = pos

    def SetTolerance(self, tol):
        pass

    def Pick(self, x, y, z, renderer):
        pass

    def GetActor(self):
        return self._actor

    def GetPickPosition(self):
        return self._pos


class TestOnPlotterClick:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        ctx.plotter = MagicMock()
        w = _bondrn.BondEditorWindow(context=ctx)  # parent must be None, not a Mock
        mw = MagicMock()
        mw.view_3d_manager = MagicMock()
        ctx.get_main_window.return_value = mw
        yield w
        w.update_timer.stop()
        w.destroy()

    def _widget(self, qapp):
        from PyQt6.QtWidgets import QWidget

        widget = QWidget()
        widget.resize(400, 300)
        return widget

    def test_empty_space_click_select_mode_clears_selection(self, win, qapp, monkeypatch):
        widget = self._widget(qapp)
        win.table.selectRow(0)
        monkeypatch.setattr(_vtk, "vtkCellPicker", lambda: _FakePicker(actor=None))
        win._on_plotter_click(10, 10, widget, None)
        assert len(win.table.selectedIndexes()) == 0

    def test_empty_space_click_create_mode_cancels_pick(self, win, qapp, monkeypatch):
        widget = self._widget(qapp)
        win.click_mode_combo.setCurrentText("Create bond")
        win._first_pick_idx = 0
        win._picked_atoms = {"Atom 1": 0}
        monkeypatch.setattr(_vtk, "vtkCellPicker", lambda: _FakePicker(actor=None))
        win._on_plotter_click(10, 10, widget, None)
        assert win._first_pick_idx is None

    def test_select_bond_mode_selects_nearest_bond(self, win, qapp, monkeypatch):
        widget = self._widget(qapp)
        # pick position on the C0-C1 bond axis (midpoint)
        monkeypatch.setattr(
            _vtk,
            "vtkCellPicker",
            lambda: _FakePicker(actor="something", pos=(0.75, 0.0, 0.0)),
        )
        win._on_plotter_click(10, 10, widget, None)
        selected_rows = sorted({i.row() for i in win.table.selectedIndexes()})
        assert selected_rows == [0]

    def test_create_bond_mode_picks_nearest_atom(self, win, qapp, monkeypatch):
        widget = self._widget(qapp)
        win.click_mode_combo.setCurrentText("Create bond")
        win.context.current_mol.GetAtomWithIdx  # sanity: real mol present
        monkeypatch.setattr(
            _vtk,
            "vtkCellPicker",
            lambda: _FakePicker(
                actor=win.context.get_main_window().view_3d_manager.atom_actor,
                pos=(0.0, 0.0, 0.0),
            ),
        )
        win._on_plotter_click(10, 10, widget, None)
        assert win._first_pick_idx == 0

    def test_create_bond_mode_wrong_actor_cancels_pick(self, win, qapp, monkeypatch):
        widget = self._widget(qapp)
        win.click_mode_combo.setCurrentText("Create bond")
        win._first_pick_idx = 1
        win._picked_atoms = {"Atom 1": 1}
        win.context.get_main_window().view_3d_manager.atom_actor = "atom-actor"
        monkeypatch.setattr(
            _vtk, "vtkCellPicker", lambda: _FakePicker(actor="other-actor")
        )
        win._on_plotter_click(10, 10, widget, None)
        assert win._first_pick_idx is None

    def test_no_molecule_returns_early(self, win, qapp, monkeypatch):
        widget = self._widget(qapp)
        win.context.current_mol = None
        monkeypatch.setattr(_vtk, "vtkCellPicker", lambda: _FakePicker(actor="x"))
        win._on_plotter_click(10, 10, widget, None)  # should not raise

    def test_no_view3d_manager_returns_early(self, win, qapp, monkeypatch):
        widget = self._widget(qapp)
        win.context.get_main_window.return_value.view_3d_manager = None
        monkeypatch.setattr(_vtk, "vtkCellPicker", lambda: _FakePicker(actor="x"))
        win._on_plotter_click(10, 10, widget, None)  # should not raise


class TestModeOverlayAndPickedLabels:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        ctx.plotter = MagicMock()
        w = _bondrn.BondEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_overlay_select_mode_removes_label(self, win):
        win.click_mode_combo.setCurrentText("Select bond")
        win._update_mode_overlay()
        win.context.plotter.remove_actor.assert_any_call("bond_editor_mode_label")

    def test_overlay_create_mode_first_pick_pending(self, win):
        win.click_mode_combo.setCurrentText("Create bond")
        win._first_pick_idx = None
        win._update_mode_overlay()
        text = win.context.plotter.add_text.call_args[0][0]
        assert "click the first atom" in text

    def test_overlay_create_mode_second_pick_pending(self, win):
        win.click_mode_combo.setCurrentText("Create bond")
        win._first_pick_idx = 2
        win._update_mode_overlay()
        text = win.context.plotter.add_text.call_args[0][0]
        assert "atom 2 picked" in text

    def test_picked_atom_labels_added_for_picked_atoms(self, win):
        win._picked_atoms = {"Atom 1": 0, "Atom 2": 2}
        win._update_picked_atom_labels()
        args, kwargs = win.context.plotter.add_point_labels.call_args
        assert len(args[0]) == 2
        assert "Atom 1: 0 (C)" in args[1]

    def test_picked_atom_labels_removed_when_empty(self, win):
        win._picked_atoms = {}
        win._update_picked_atom_labels()
        win.context.plotter.remove_actor.assert_any_call("bond_editor_atom_labels")

    def test_cancel_bond_pick_with_pending_pick(self, win):
        win._first_pick_idx = 1
        win._picked_atoms = {"Atom 1": 1}
        win._cancel_bond_pick()
        assert win._first_pick_idx is None
        msg = win.context.show_status_message.call_args[0][0]
        assert "cleared" in msg


class TestNearestHelpersDirect:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _bondrn.BondEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_nearest_atom_to_point(self, win):
        mol = win.context.current_molecule
        idx = win._nearest_atom_to_point(mol, (1.4, 0.0, 0.0))
        assert idx == 1

    def test_nearest_bond_to_point_within_tolerance(self, win):
        mol = win.context.current_molecule
        pair = win._nearest_bond_to_point(mol, (0.75, 0.0, 0.0))
        assert pair == (0, 1)

    def test_nearest_bond_to_point_out_of_range(self, win):
        mol = win.context.current_molecule
        pair = win._nearest_bond_to_point(mol, (50.0, 50.0, 50.0))
        assert pair is None

    def test_select_bond_row_by_pair_scrolls_and_reports(self, win):
        win._select_bond_row_by_pair((1, 2))
        selected_rows = sorted({i.row() for i in win.table.selectedIndexes()})
        assert selected_rows == [1]
        msg = win.context.show_status_message.call_args[0][0]
        assert "Selected bond 1-2" in msg

    def test_select_bond_row_by_pair_no_match(self, win):
        win.context.show_status_message.reset_mock()
        win._select_bond_row_by_pair((5, 6))
        assert len(win.table.selectedIndexes()) == 0


class TestEditOperationsRealChem:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _bondrn.BondEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_add_bond_success(self, win):
        win.add_bond(0, 2)
        new_mol = win.context.current_molecule
        assert new_mol.GetBondBetweenAtoms(0, 2) is not None
        win.context.push_undo_checkpoint.assert_called()

    def test_add_bond_self_bond_warns(self, win):
        win.context.show_status_message.reset_mock()
        win.add_bond(0, 0)
        msg = win.context.show_status_message.call_args[0][0]
        assert "itself" in msg

    def test_add_bond_existing_bond_warns(self, win):
        win.context.show_status_message.reset_mock()
        win.add_bond(0, 1)
        msg = win.context.show_status_message.call_args[0][0]
        assert "already exists" in msg

    def test_add_bond_no_refresh_3d_view_falls_back_to_reset_camera(self, win):
        del win.context.refresh_3d_view
        win.add_bond(0, 2)
        win.context.reset_3d_camera.assert_called_once()

    def test_add_bond_exception_shows_critical(self, win, qapp, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _bondrn.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.add_bond(0, 2)
        assert calls
        assert "Failed to add bond" in calls[0][2]

    def test_delete_selected_bonds(self, win):
        win.table.selectRow(0)
        win.delete_selected_bonds()
        new_mol = win.context.current_molecule
        assert new_mol.GetNumBonds() == 1

    def test_delete_selected_bonds_exception_shows_critical(self, win, qapp, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _bondrn.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.table.selectRow(0)
        win.delete_selected_bonds()
        assert calls
        assert "Failed to delete bonds" in calls[0][2]

    def test_on_type_changed_updates_bond(self, win):
        win.on_type_changed(0, "Double")
        new_mol = win.context.current_molecule
        bond = new_mol.GetBondBetweenAtoms(0, 1)
        assert str(bond.GetBondType()).rsplit(".", 1)[-1] == "DOUBLE"

    def test_on_type_changed_aromatic_flags_atoms(self, win):
        win.on_type_changed(0, "Aromatic")
        new_mol = win.context.current_molecule
        assert new_mol.GetAtomWithIdx(0).GetIsAromatic()
        assert new_mol.GetAtomWithIdx(1).GetIsAromatic()

    def test_on_type_changed_no_pair_is_noop(self, win):
        win.table.setRowCount(0)
        win.on_type_changed(0, "Double")  # no exception

    def test_on_type_changed_exception_shows_critical(self, win, qapp, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _bondrn.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.on_type_changed(0, "Double")
        assert calls
        assert "Failed to change bond type" in calls[0][2]

    def test_commit_sanitize_fallback_on_bad_valence(self, win):
        rw = _Chem.RWMol()
        c = rw.AddAtom(_Chem.Atom(6))
        for _ in range(5):
            h = rw.AddAtom(_Chem.Atom(1))
            rw.AddBond(c, h, _Chem.BondType.SINGLE)
        win._commit(rw, "overvalent test")
        assert win.context.current_molecule.GetNumAtoms() == 6


class TestMovingSideAndBondLength:
    def test_moving_side_returns_none_for_ring_bond(self, qapp):
        ctx = _real_ctx(mol=_real_mol(ring=True))
        w = _bondrn.BondEditorWindow(context=ctx)
        mol = ctx.current_molecule
        assert w._moving_side(mol, 0, 1) is None
        w.update_timer.stop()
        w.destroy()

    def test_moving_side_non_ring_returns_far_side(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _bondrn.BondEditorWindow(context=ctx)
        mol = ctx.current_molecule
        side = w._moving_side(mol, 0, 1)
        assert side == {1, 2}
        w.update_timer.stop()
        w.destroy()

    def test_on_item_changed_non_numeric_reloads(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _bondrn.BondEditorWindow(context=ctx)
        item = w.table.item(0, w.COL_LEN)
        item.setText("abc")
        assert w.table.rowCount() == 2  # reload restores rows
        w.update_timer.stop()
        w.destroy()

    def test_on_item_changed_non_length_column_ignored(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _bondrn.BondEditorWindow(context=ctx)
        item = w.table.item(0, w.COL_A1)
        w.on_item_changed(item)  # no exception, no-op
        w.update_timer.stop()
        w.destroy()

    def test_on_item_changed_negative_length_warns(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _bondrn.BondEditorWindow(context=ctx)
        w.context.show_status_message.reset_mock()
        item = w.table.item(0, w.COL_LEN)
        item.setText("-1.0")
        msg = w.context.show_status_message.call_args[0][0]
        assert "positive" in msg
        w.update_timer.stop()
        w.destroy()

    def test_set_bond_length_moves_far_side(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _bondrn.BondEditorWindow(context=ctx)
        w.set_bond_length(0, 2.0)
        new_mol = w.context.current_molecule
        conf = new_mol.GetConformer()
        p0 = conf.GetAtomPosition(0)
        p1 = conf.GetAtomPosition(1)
        dist = ((p0.x - p1.x) ** 2 + (p0.y - p1.y) ** 2 + (p0.z - p1.z) ** 2) ** 0.5
        assert dist == pytest.approx(2.0, abs=1e-6)
        w.update_timer.stop()
        w.destroy()

    def test_set_bond_length_ring_bond_refused(self, qapp):
        ctx = _real_ctx(mol=_real_mol(ring=True))
        w = _bondrn.BondEditorWindow(context=ctx)
        w.context.show_status_message.reset_mock()
        w.set_bond_length(0, 2.0)
        msg = w.context.show_status_message.call_args[0][0]
        assert "ring" in msg
        w.update_timer.stop()
        w.destroy()

    def test_set_bond_length_coincident_atoms_refused(self, qapp):
        mol = _real_mol()
        conf = mol.GetConformer()
        conf.SetAtomPosition(1, _Point3D(0.0, 0.0, 0.0))  # same as atom 0
        ctx = _real_ctx(mol=mol)
        w = _bondrn.BondEditorWindow(context=ctx)
        w.context.show_status_message.reset_mock()
        w.set_bond_length(0, 2.0)
        msg = w.context.show_status_message.call_args[0][0]
        assert "coincident" in msg
        w.update_timer.stop()
        w.destroy()

    def test_set_bond_length_no_conformer_returns_early(self, qapp):
        ctx = _real_ctx(mol=_FakeMol(atoms=[_FakeAtom(0), _FakeAtom(1)], bonds=[_FakeBond(0, 0, 1)]))
        w = _bondrn.BondEditorWindow(context=ctx)
        w.set_bond_length(0, 2.0)  # should not raise
        w.update_timer.stop()
        w.destroy()

    def test_set_bond_length_exception_shows_critical(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        ctx = _real_ctx(mol=_real_mol())
        w = _bondrn.BondEditorWindow(context=ctx)
        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _bondrn.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        w.set_bond_length(0, 2.0)
        assert calls
        assert "Failed to set bond length" in calls[0][2]
        w.update_timer.stop()
        w.destroy()


class TestHighlightSelectedBonds:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        ctx.plotter = MagicMock()
        w = _bondrn.BondEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_no_selection_removes_actor(self, win):
        win.table.clearSelection()
        win.highlight_selected_bonds()
        win.context.plotter.remove_actor.assert_any_call("bond_editor_selection")

    def test_selection_adds_tube_mesh(self, win):
        win.table.selectRow(0)
        win.highlight_selected_bonds()
        assert win.context.plotter.add_mesh.called
        name = win.context.plotter.add_mesh.call_args[1]["name"]
        assert name == "bond_editor_selection"

    def test_camera_position_restored(self, win):
        win.context.plotter.camera_position = "cam-state"
        win.table.selectRow(0)
        win.highlight_selected_bonds()
        assert win.context.plotter.camera_position == "cam-state"

    def test_no_plotter_is_noop(self, win):
        win.context.plotter = None
        win.highlight_selected_bonds()  # should not raise


class TestCloseEventWithPlotter:
    def test_close_removes_all_actors(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        ctx.plotter = MagicMock()
        w = _bondrn.BondEditorWindow(context=ctx)
        w.close()
        ctx.plotter.remove_actor.assert_any_call("bond_editor_selection")
        ctx.plotter.remove_actor.assert_any_call("bond_editor_mode_label")
        ctx.plotter.remove_actor.assert_any_call("bond_editor_atom_labels")
        ctx.plotter.render.assert_called()
        w.destroy()
