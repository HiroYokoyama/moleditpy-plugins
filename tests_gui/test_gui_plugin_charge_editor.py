"""
Headless GUI tests for the Charge Editor plugin.

Covers: ChargeEditorWindow, initialize().

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_charge_editor.py
"""

from __future__ import annotations

import contextlib
import sys
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


# Real numpy/rdkit/pyvista/vtk pulled in up front so they're cached before the
# chemistry-blocking finder snapshots sys.modules — otherwise the first-ever
# import would be intercepted and replaced with a MagicMock. Guarded so CI's
# bare test-gui job (only pytest+PyQt6 installed) skips this instead of
# erroring at collection.
_np = pytest.importorskip("numpy")
_pv = pytest.importorskip("pyvista")
_vtk = pytest.importorskip("vtk")
_Chem = pytest.importorskip("rdkit.Chem")
_Point3D = pytest.importorskip("rdkit.Geometry").Point3D


@contextlib.contextmanager
def _mock_chemistry_keep_real_chem():
    """Like mock_chemistry_imports(), but numpy/rdkit/pyvista/vtk resolve to the
    real packages (all installed in this environment), so the plugin's real
    conformer math, RWMol sanitize, and vtkCellPicker calls actually run.
    """
    keep_prefixes = ("numpy", "rdkit", "pyvista", "pyvistaqt", "vtk", "vtkmodules")
    real_mods = {
        k: v for k, v in sys.modules.items() if k.split(".")[0] in keep_prefixes
    }
    with mock_chemistry_imports():
        sys.modules.update(real_mods)
        yield


with _mock_chemistry_keep_real_chem():
    _chargern = load_plugin_for_gui(CHARGE_EDITOR_PATH)


def _real_mol():
    """Real RDKit RWMol with a conformer: C0-C1-O2."""
    rw = _Chem.RWMol()
    rw.AddAtom(_Chem.Atom(6))
    rw.AddAtom(_Chem.Atom(6))
    rw.AddAtom(_Chem.Atom(8))
    rw.AddBond(0, 1, _Chem.BondType.SINGLE)
    rw.AddBond(1, 2, _Chem.BondType.SINGLE)
    conf = _Chem.Conformer(3)
    for i, (x, y, z) in enumerate([(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (3.0, 0.0, 0.0)]):
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

    def test_show_editor_creates_new_window(self, qapp):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        ctx.get_window.return_value = None
        _charge.initialize(ctx)
        show_editor = ctx.add_menu_action.call_args[0][1]
        show_editor()
        win = ctx.register_window.call_args[0][1]
        assert isinstance(win, _charge.ChargeEditorWindow)
        win.update_timer.stop()
        win.close()
        win.destroy()

    def test_show_editor_reuses_existing_window(self, qapp):
        ctx = MagicMock()
        _charge.initialize(ctx)
        fake_win = MagicMock()
        ctx.get_window.return_value = fake_win
        show_editor = ctx.add_menu_action.call_args[0][1]
        show_editor()
        fake_win.show.assert_called_once()
        fake_win.raise_.assert_called_once()
        fake_win.activateWindow.assert_called_once()
        fake_win.load_molecule.assert_called_once()


# ===========================================================================
# _ClickFilter — real QMouseEvents
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
        f = _charge._ClickFilter(lambda x, y, obj, mods: calls.append((x, y)))
        press, release = self._events((10, 10), (12, 13))
        obj = MagicMock()
        f.eventFilter(obj, press)
        f.eventFilter(obj, release)
        assert calls == [(12, 13)]

    def test_drag_beyond_5px_is_ignored(self, qapp):
        calls = []
        f = _charge._ClickFilter(lambda x, y, obj, mods: calls.append((x, y)))
        press, release = self._events((10, 10), (30, 30))
        obj = MagicMock()
        f.eventFilter(obj, press)
        f.eventFilter(obj, release)
        assert calls == []

    def test_never_consumes_events(self, qapp):
        f = _charge._ClickFilter(lambda *a: None)
        press, release = self._events((0, 0), (0, 0))
        obj = MagicMock()
        assert f.eventFilter(obj, press) is False
        assert f.eventFilter(obj, release) is False

    def test_release_without_press_is_ignored(self, qapp):
        calls = []
        f = _charge._ClickFilter(lambda x, y, obj, mods: calls.append((x, y)))
        _, release = self._events((0, 0), (5, 5))
        f.eventFilter(MagicMock(), release)
        assert calls == []


# ===========================================================================
# Real RDKit / numpy / pyvista / vtk — conformer-dependent code paths
# ===========================================================================


class TestPlotterPickingRealInteractor:
    def test_enable_installs_event_filter_on_real_interactor(self, qapp):
        interactor = MagicMock()
        plotter = MagicMock()
        plotter.interactor = interactor
        ctx = _real_ctx(mol=None)
        ctx.plotter = plotter
        w = _chargern.ChargeEditorWindow(context=ctx)
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
        w = _chargern.ChargeEditorWindow(context=ctx)
        assert w._click_filter is None
        w.close()
        w.destroy()

    def test_enable_noop_when_plotter_is_none(self, qapp):
        ctx = _real_ctx(mol=None)
        w = _chargern.ChargeEditorWindow(context=ctx)
        assert w._click_filter is None
        w.update_timer.stop()
        w.close()
        w.destroy()

    def test_disable_swallows_missing_plotter(self, qapp):
        interactor = MagicMock()
        plotter = MagicMock()
        plotter.interactor = interactor
        ctx = _real_ctx(mol=None)
        ctx.plotter = plotter
        w = _chargern.ChargeEditorWindow(context=ctx)
        ctx.plotter = None
        w._disable_plotter_picking()  # must not raise
        assert w._click_filter is None
        w.update_timer.stop()
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
        w = _chargern.ChargeEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def _widget(self):
        from PyQt6.QtWidgets import QWidget

        widget = QWidget()
        widget.resize(400, 300)
        return widget

    def test_no_plotter_returns_early(self, win, qapp):
        win.context.plotter = None
        win._on_plotter_click(10, 10, self._widget(), None)  # must not raise

    def test_no_actor_returns_early(self, win, qapp, monkeypatch):
        monkeypatch.setattr(_vtk, "vtkCellPicker", lambda: _FakePicker(actor=None))
        win._on_plotter_click(10, 10, self._widget(), None)
        assert len(win.table.selectedIndexes()) == 0

    def test_no_molecule_returns_early(self, win, qapp, monkeypatch):
        win.context.current_mol = None
        monkeypatch.setattr(_vtk, "vtkCellPicker", lambda: _FakePicker(actor="x"))
        win._on_plotter_click(10, 10, self._widget(), None)  # must not raise

    def test_no_conformer_returns_early(self, win, qapp, monkeypatch):
        mock_mol = MagicMock()
        mock_mol.GetNumConformers.return_value = 0
        win.context.current_mol = mock_mol
        monkeypatch.setattr(_vtk, "vtkCellPicker", lambda: _FakePicker(actor="x"))
        win._on_plotter_click(10, 10, self._widget(), None)  # must not raise

    def test_picks_nearest_atom_and_selects_row(self, win, qapp, monkeypatch):
        from PyQt6.QtCore import Qt

        monkeypatch.setattr(
            _vtk, "vtkCellPicker", lambda: _FakePicker(actor="x", pos=(1.4, 0.0, 0.0))
        )
        win._on_plotter_click(10, 10, self._widget(), Qt.KeyboardModifier.NoModifier)
        rows = sorted({i.row() for i in win.table.selectedIndexes()})
        assert rows == [1]

    def test_ctrl_click_adds_to_selection(self, win, qapp, monkeypatch):
        from PyQt6.QtCore import Qt

        monkeypatch.setattr(
            _vtk, "vtkCellPicker", lambda: _FakePicker(actor="x", pos=(0.0, 0.0, 0.0))
        )
        win._on_plotter_click(10, 10, self._widget(), Qt.KeyboardModifier.NoModifier)
        monkeypatch.setattr(
            _vtk, "vtkCellPicker", lambda: _FakePicker(actor="x", pos=(3.0, 0.0, 0.0))
        )
        win._on_plotter_click(10, 10, self._widget(), Qt.KeyboardModifier.ControlModifier)
        rows = sorted({i.row() for i in win.table.selectedIndexes()})
        assert rows == [0, 2]

    def test_exception_in_pick_is_silenced(self, win, qapp, monkeypatch):
        def boom():
            raise RuntimeError("picker boom")

        monkeypatch.setattr(_vtk, "vtkCellPicker", boom)
        win._on_plotter_click(10, 10, self._widget(), None)  # must not raise


class TestNearestAtomAndSelectRow:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _chargern.ChargeEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_nearest_atom_to_point(self, win):
        mol = win.context.current_molecule
        assert win._nearest_atom_to_point(mol, (1.4, 0.0, 0.0)) == 1
        assert win._nearest_atom_to_point(mol, (0.1, 0.0, 0.0)) == 0

    def test_select_atom_row_selects_matching_row(self, win):
        win._select_atom_row(2, False)
        rows = sorted({i.row() for i in win.table.selectedIndexes()})
        assert rows == [2]

    def test_select_atom_row_no_match_leaves_selection_empty(self, win):
        win._select_atom_row(99, False)
        assert len(win.table.selectedIndexes()) == 0

    def test_ctrl_reselect_same_atom_toggles_off(self, win):
        win._select_atom_row(0, False)
        win._select_atom_row(0, True)
        rows = sorted({i.row() for i in win.table.selectedIndexes()})
        assert rows == []


class TestSignatureAndUpdateReal:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _chargern.ChargeEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_signature_exception_returns_none(self, win):
        bad_mol = MagicMock()
        bad_mol.GetNumAtoms.side_effect = RuntimeError("boom")
        assert win.get_mol_signature(bad_mol) is None

    def test_check_molecule_update_reloads_on_change(self, win):
        new_mol = _real_mol()
        win.context.current_molecule = new_mol
        win.check_molecule_update()
        assert win.last_seen_signature == win.get_mol_signature(new_mol)

    def test_check_molecule_update_noop_when_unchanged(self, win):
        rows_before = win.table.rowCount()
        win.check_molecule_update()
        assert win.table.rowCount() == rows_before

    def test_check_molecule_update_exception_silenced(self, win, monkeypatch):
        monkeypatch.setattr(
            win, "get_mol_signature", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.check_molecule_update()  # must not raise


class TestLoadMoleculeCustomSymbolAndHCountFallback:
    def test_custom_symbol_shown_in_table(self, qapp):
        mol = _real_mol()
        mol.GetAtomWithIdx(0).SetProp("custom_symbol", "Xx")
        ctx = _real_ctx(mol=mol)
        w = _chargern.ChargeEditorWindow(context=ctx)
        assert w.table.item(0, w.COL_SYMBOL).text() == "Xx"
        w.update_timer.stop()
        w.destroy()

    def test_h_count_falls_back_to_zero_on_error(self, qapp):
        rw = _Chem.RWMol()
        rw.AddAtom(_Chem.Atom(6))  # not sanitized -> GetTotalNumHs() raises
        mol = rw.GetMol()
        ctx = _real_ctx(mol=mol)
        w = _chargern.ChargeEditorWindow(context=ctx)
        assert w.table.item(0, w.COL_HS).text() == "0"
        w.update_timer.stop()
        w.destroy()


class TestEditOperationsRealChem:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _chargern.ChargeEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_on_charge_changed_success(self, win):
        win.on_charge_changed(1, -1)
        new_mol = win.context.current_molecule
        assert new_mol.GetAtomWithIdx(1).GetFormalCharge() == -1
        win.context.push_undo_checkpoint.assert_called()

    def test_on_charge_changed_exception_shows_critical(self, win, qapp, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _chargern.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.on_charge_changed(0, 1)
        assert calls
        assert "Failed to set charge" in calls[0][2]

    def test_on_radical_changed_success(self, win):
        # Atom 0 is a terminal, open-valence carbon: RDKit's sanitize step
        # recomputes radicals to satisfy valence, so we assert the commit
        # succeeded rather than a specific post-sanitize radical count.
        win.on_radical_changed(0, 2)
        new_mol = win.context.current_molecule
        assert new_mol.GetNumAtoms() == 3
        win.context.push_undo_checkpoint.assert_called()

    def test_on_radical_changed_exception_shows_critical(self, win, qapp, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _chargern.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.on_radical_changed(0, 1)
        assert calls
        assert "Failed to set radicals" in calls[0][2]

    def test_clear_all_charges_success(self, win):
        win.on_charge_changed(0, 2)
        win.clear_all_charges()
        new_mol = win.context.current_molecule
        assert all(a.GetFormalCharge() == 0 for a in new_mol.GetAtoms())

    def test_clear_all_charges_exception_shows_critical(self, win, qapp, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _chargern.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.clear_all_charges()
        assert calls
        assert "Failed to clear charges" in calls[0][2]

    def test_clear_all_radicals_success(self, win):
        win.clear_all_radicals()
        new_mol = win.context.current_molecule
        assert new_mol.GetNumAtoms() == 3
        win.context.push_undo_checkpoint.assert_called()

    def test_clear_all_radicals_exception_shows_critical(self, win, qapp, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _chargern.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.clear_all_radicals()
        assert calls
        assert "Failed to clear radicals" in calls[0][2]

    def test_commit_sanitize_fallback_on_bad_valence(self, win):
        rw = _Chem.RWMol()
        c = rw.AddAtom(_Chem.Atom(6))
        for _ in range(5):
            h = rw.AddAtom(_Chem.Atom(1))
            rw.AddBond(c, h, _Chem.BondType.SINGLE)
        win._commit(rw, "overvalent test")
        assert win.context.current_molecule.GetNumAtoms() == 6

    def test_commit_no_refresh_3d_view_attr_is_fine(self, win):
        del win.context.refresh_3d_view
        win.on_charge_changed(0, 1)  # must not raise
        assert win.context.current_molecule.GetAtomWithIdx(0).GetFormalCharge() == 1

    def test_commit_secondary_refresh_exception_silenced(self, win):
        win.context.refresh_2d_scene.side_effect = RuntimeError("boom")
        win.on_charge_changed(0, 1)  # must not raise
        assert win.context.current_molecule.GetAtomWithIdx(0).GetFormalCharge() == 1


class TestHighlightRadiusReal:
    def test_scales_with_real_periodic_table(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _chargern.ChargeEditorWindow(context=ctx)
        mol = ctx.current_molecule
        r_c = w._highlight_radius(mol, 0)
        r_o = w._highlight_radius(mol, 2)
        assert r_c > 0
        assert r_o < r_c  # RDKit's table gives O a smaller VDW radius than C
        w.update_timer.stop()
        w.destroy()

    def test_ghost_atom_uses_fixed_fallback(self, qapp):
        rw = _Chem.RWMol()
        rw.AddAtom(_Chem.Atom(0))  # dummy/wildcard atom
        mol = rw.GetMol()
        ctx = _real_ctx(mol=mol)
        w = _chargern.ChargeEditorWindow(context=ctx)
        assert abs(w._highlight_radius(mol, 0) - 0.45) < 1e-9
        w.update_timer.stop()
        w.destroy()


class TestHighlightSelectedAtoms:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        ctx.plotter = MagicMock()
        w = _chargern.ChargeEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_no_selection_removes_actor(self, win):
        win.table.clearSelection()
        win.highlight_selected_atoms()
        win.context.plotter.remove_actor.assert_any_call("charge_editor_selection")

    def test_selection_adds_glyph_mesh(self, win):
        win.table.selectRow(0)
        win.highlight_selected_atoms()
        assert win.context.plotter.add_mesh.called
        name = win.context.plotter.add_mesh.call_args[1]["name"]
        assert name == "charge_editor_selection"

    def test_camera_position_restored(self, win):
        win.context.plotter.camera_position = "cam-state"
        win.table.selectRow(0)
        win.highlight_selected_atoms()
        assert win.context.plotter.camera_position == "cam-state"

    def test_no_plotter_is_noop(self, win):
        win.context.plotter = None
        win.highlight_selected_atoms()  # must not raise


class TestCloseEventWithPlotter:
    def test_close_removes_selection_actor(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        ctx.plotter = MagicMock()
        w = _chargern.ChargeEditorWindow(context=ctx)
        w.close()
        ctx.plotter.remove_actor.assert_any_call("charge_editor_selection")
        ctx.plotter.render.assert_called()
