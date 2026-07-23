"""
Headless GUI tests for the XYZ Editor plugin.

Covers: XYZEditorWindow.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_xyz_editor.py
"""

from __future__ import annotations

import contextlib
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

XYZ_EDITOR_PATH = PLUGINS_DIR / "XYZ_Editor" / "xyz_editor.py"

with mock_chemistry_imports():
    _xyz_editor = load_plugin_for_gui(XYZ_EDITOR_PATH)


# ---------------------------------------------------------------------------
# Second module instance with real numpy/rdkit/pyvista/vtk kept, so tests that
# need real molecule reconstruction (apply_changes, duplicate_atoms,
# adjust_hydrogens, _on_plotter_click, highlight_selected_atoms) exercise the
# actual code paths instead of chasing MagicMock attribute chains. Mirrors
# tests_gui/test_gui_plugin_bond_editor.py.
# ---------------------------------------------------------------------------

_np = pytest.importorskip("numpy")
_pv = pytest.importorskip("pyvista")
_vtk = pytest.importorskip("vtk")
_Chem = pytest.importorskip("rdkit.Chem")
_Point3D = pytest.importorskip("rdkit.Geometry").Point3D


@contextlib.contextmanager
def _mock_chemistry_keep_real_chem():
    keep_prefixes = ("numpy", "rdkit", "pyvista", "pyvistaqt", "vtk", "vtkmodules")
    real_mods = {
        k: v for k, v in sys.modules.items() if k.split(".")[0] in keep_prefixes
    }
    with mock_chemistry_imports():
        sys.modules.update(real_mods)
        yield


with _mock_chemistry_keep_real_chem():
    _xyzrn = load_plugin_for_gui(XYZ_EDITOR_PATH)


def _real_mol():
    """Real RDKit RWMol with a conformer: C0-C1(single)-O2(single)."""
    rw = _Chem.RWMol()
    rw.AddAtom(_Chem.Atom(6))
    rw.AddAtom(_Chem.Atom(6))
    rw.AddAtom(_Chem.Atom(8))
    rw.AddBond(0, 1, _Chem.BondType.SINGLE)
    rw.AddBond(1, 2, _Chem.BondType.SINGLE)
    conf = _Chem.Conformer(3)
    for i, (x, y, z) in enumerate([(0, 0, 0), (1.5, 0, 0), (2.5, 1.0, 0)]):
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


# ===========================================================================
# Fakes with a conformer (load_molecule calls mol.GetConformer())
# ===========================================================================

from gui_test_helpers import FakeAtom, FakeMol


class _Pos:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class _FakeConf:
    def __init__(self, coords):
        self._coords = coords

    def GetAtomPosition(self, i):
        return _Pos(*self._coords[i])

    def GetPositions(self):
        return self._coords


class _MolWithConf(FakeMol):
    def __init__(self, atoms, bonds=(), coords=None):
        super().__init__(atoms, bonds)
        self._conf = _FakeConf(
            coords or [(float(i), 0.5, -0.5) for i in range(len(atoms))]
        )

    def GetNumConformers(self):
        return 1

    def GetConformer(self):
        return self._conf


def _win_with_mol(symbols=("C", "O", "H")):
    ctx = _ctx_no_mol()
    ctx.plotter = None
    mol = _MolWithConf([FakeAtom(i, s) for i, s in enumerate(symbols)])
    ctx.current_molecule = mol
    w = _xyz_editor.XYZEditorWindow(context=ctx)
    return w, ctx, mol


def _teardown(w):
    w.update_timer.stop()
    w.destroy()


# ===========================================================================
# Table population from a molecule
# ===========================================================================


class TestXYZEditorLoadMolecule:
    def test_rows_match_atom_count(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O", "H"))
        assert w.table.rowCount() == 3
        _teardown(w)

    def test_symbols_and_indices_in_table(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O"))
        assert w.table.item(0, 0).text() == "0"
        assert w.table.item(0, 1).text() == "C"
        assert w.table.item(1, 1).text() == "O"
        _teardown(w)

    def test_coordinates_formatted_5_decimals(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        assert w.table.item(0, 2).text() == "0.00000"
        assert w.table.item(0, 3).text() == "0.50000"
        assert w.table.item(0, 4).text() == "-0.50000"
        _teardown(w)

    def test_index_column_not_editable(self, qapp):
        from PyQt6.QtCore import Qt

        w, ctx, mol = _win_with_mol(("C",))
        flags = w.table.item(0, 0).flags()
        assert not flags & Qt.ItemFlag.ItemIsEditable
        _teardown(w)

    def test_reload_replaces_rows(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O", "H"))
        ctx.current_molecule = _MolWithConf([FakeAtom(0, "N")])
        w.load_molecule()
        assert w.table.rowCount() == 1
        assert w.table.item(0, 1).text() == "N"
        _teardown(w)

    def test_check_molecule_update_detects_new_mol(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        ctx.current_molecule = _MolWithConf([FakeAtom(0, "S"), FakeAtom(1, "S")])
        w.check_molecule_update()
        assert w.table.rowCount() == 2
        _teardown(w)

    def test_check_molecule_update_noop_when_unchanged(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        w.load_molecule = MagicMock()
        w.check_molecule_update()
        w.load_molecule.assert_not_called()
        _teardown(w)


# ===========================================================================
# Row editing: add / remove / select
# ===========================================================================


class TestXYZEditorRowOps:
    def test_add_atom_row_appends_placeholder(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        w.add_atom_row()
        assert w.table.rowCount() == 2
        assert w.table.item(1, 0).text() == "+"
        assert w.table.item(1, 1).text() == "C"
        _teardown(w)

    def test_selected_atom_indices_skips_placeholder_rows(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O"))
        w.add_atom_row()
        w.table.selectAll()
        assert w._selected_atom_indices() == {0, 1}
        _teardown(w)

    def test_selected_atom_indices_empty_without_selection(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O"))
        assert w._selected_atom_indices() == set()
        _teardown(w)

    def test_unselect_all_clears_selection(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O"))
        w.table.selectAll()
        assert w.table.selectedIndexes()
        w.unselect_all()
        assert not w.table.selectedIndexes()
        _teardown(w)

    def test_delete_without_selection_warns(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O"))
        w.apply_changes = MagicMock()
        w.delete_selected_atoms()
        ctx.show_status_message.assert_called_with("No atoms selected to delete.")
        w.apply_changes.assert_not_called()
        _teardown(w)

    def test_delete_with_selection_removes_and_applies(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O"))
        w.apply_changes = MagicMock()
        w.table.selectRow(0)
        w.delete_selected_atoms()
        assert w.table.rowCount() == 1
        w.apply_changes.assert_called_once()
        _teardown(w)


# ===========================================================================
# Duplicate / Adjust-H guard paths (recent features)
# ===========================================================================


class TestXYZEditorDuplicateGuards:
    def test_duplicate_without_molecule_warns(self, qapp):
        ctx = _ctx_no_mol()
        ctx.plotter = None
        w = _xyz_editor.XYZEditorWindow(context=ctx)
        w.duplicate_atoms()
        ctx.show_status_message.assert_called_with("No molecule to duplicate.")
        _teardown(w)

    def test_duplicate_with_unapplied_selection_warns(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        w.add_atom_row()
        w.table.clearSelection()
        w.table.selectRow(1)  # only the unapplied '+' row
        w.duplicate_atoms()
        ctx.show_status_message.assert_called_with(
            "Selected rows are not applied yet — press Apply first."
        )
        _teardown(w)

    def test_adjust_h_without_molecule_warns(self, qapp):
        ctx = _ctx_no_mol()
        ctx.plotter = None
        w = _xyz_editor.XYZEditorWindow(context=ctx)
        w.adjust_hydrogens()
        ctx.show_status_message.assert_called_with(
            "No molecule to adjust hydrogens on."
        )
        _teardown(w)

    def test_adjust_h_with_unapplied_selection_warns(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        w.add_atom_row()
        w.table.clearSelection()
        w.table.selectRow(1)
        w.adjust_hydrogens()
        ctx.show_status_message.assert_called_with(
            "Selected rows are not applied yet — press Apply first."
        )
        _teardown(w)


# ===========================================================================
# XYZ generation / clipboard / save
# ===========================================================================


class TestXYZEditorExport:
    def test_generate_xyz_content_layout(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O"))
        lines = w._generate_xyz_content()
        assert lines[0] == "2"
        assert lines[1] == "Generated by MoleditPy XYZ Editor"
        assert lines[2].split() == ["C", "0.00000", "0.50000", "-0.50000"]
        assert len(lines) == 4
        _teardown(w)

    def test_copy_to_clipboard_sets_text(self, qapp):
        from PyQt6.QtGui import QGuiApplication

        w, ctx, mol = _win_with_mol(("C",))
        w.copy_to_clipboard()
        text = QGuiApplication.clipboard().text()
        assert text.startswith("1\nGenerated by MoleditPy XYZ Editor")
        ctx.show_status_message.assert_called_with("XYZ data copied to clipboard.")
        _teardown(w)

    def test_save_as_xyz_cancel_writes_nothing(self, qapp, tmp_path, monkeypatch):
        w, ctx, mol = _win_with_mol(("C",))
        monkeypatch.setattr(
            _xyz_editor.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        w.save_as_xyz()
        assert list(tmp_path.iterdir()) == []
        _teardown(w)

    def test_save_as_xyz_appends_extension_and_writes(self, qapp, tmp_path, monkeypatch):
        w, ctx, mol = _win_with_mol(("C",))
        target = tmp_path / "out"
        monkeypatch.setattr(
            _xyz_editor.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(target), "")),
        )
        w.save_as_xyz()
        written = tmp_path / "out.xyz"
        assert written.exists()
        assert written.read_text(encoding="utf-8").startswith("1\n")
        _teardown(w)


# ===========================================================================
# Close behavior — timer stop + window-registry unregister (reopen fix)
# ===========================================================================


class TestXYZEditorClose:
    def test_close_stops_timer(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        assert w.update_timer.isActive()
        w.close()
        assert not w.update_timer.isActive()
        w.destroy()

    def test_close_unregisters_main_panel(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        w.close()
        ctx.register_window.assert_called_with("main_panel", None)
        w.destroy()

    def test_close_removes_selection_actor(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        plotter = MagicMock()
        ctx.plotter = plotter
        w.close()
        plotter.remove_actor.assert_called_with("xyz_selection")
        plotter.render.assert_called()
        w.destroy()


# ===========================================================================
# initialize() — menu, persistence, document reset
# ===========================================================================


class TestXYZEditorInitialize:
    def _init_ctx(self):
        ctx = _ctx_no_mol()
        ctx.get_window.return_value = None
        handlers = {}
        ctx.register_save_handler.side_effect = lambda fn: handlers.setdefault(
            "save", fn
        )
        ctx.register_load_handler.side_effect = lambda fn: handlers.setdefault(
            "load", fn
        )
        ctx.register_document_reset_handler.side_effect = (
            lambda fn: handlers.setdefault("reset", fn)
        )
        _xyz_editor.initialize(ctx)
        return ctx, handlers

    def test_registers_menu_action(self, qapp):
        ctx, handlers = self._init_ctx()
        path = ctx.add_menu_action.call_args[0][0]
        assert path == "3D Edit/XYZ Editor..."

    def test_registers_all_handlers(self, qapp):
        ctx, handlers = self._init_ctx()
        assert set(handlers) == {"save", "load", "reset"}

    def test_save_state_empty_without_molecule(self, qapp):
        ctx, handlers = self._init_ctx()
        ctx.current_molecule = None
        assert handlers["save"]() == {}

    def test_save_state_collects_custom_labels(self, qapp):
        class _Labelled(FakeAtom):
            def HasProp(self, name):
                return name == "custom_symbol"

            def GetProp(self, name):
                return "Xx"

        ctx, handlers = self._init_ctx()
        ctx.current_molecule = _MolWithConf([_Labelled(0, "C"), FakeAtom(1, "O")])
        assert handlers["save"]() == {"custom_labels": {0: "Xx"}}

    def test_load_state_ignores_non_dict_labels(self, qapp):
        ctx, handlers = self._init_ctx()
        ctx.current_molecule = _MolWithConf([FakeAtom(0, "C")])
        handlers["load"]({"custom_labels": ["not", "a", "dict"]})  # must not raise

    def test_document_reset_reloads_open_window(self, qapp):
        ctx, handlers = self._init_ctx()
        win = MagicMock()
        ctx.get_window.return_value = win
        handlers["reset"]()
        win.load_molecule.assert_called_once()

    def test_show_editor_reuses_existing_window(self, qapp):
        ctx, handlers = self._init_ctx()
        win = MagicMock()
        ctx.get_window.return_value = win
        show_editor = ctx.add_menu_action.call_args[0][1]
        show_editor()
        win.show.assert_called_once()
        win.raise_.assert_called_once()
        win.activateWindow.assert_called_once()
        win.load_molecule.assert_called_once()

    def test_show_editor_creates_new_window_when_none_open(self, qapp):
        ctx, handlers = self._init_ctx()
        ctx.get_window.return_value = None
        show_editor = ctx.add_menu_action.call_args[0][1]
        show_editor()
        new_win = ctx.register_window.call_args[0][1]
        assert new_win is not None
        new_win.update_timer.stop()
        new_win.destroy()

    def test_load_state_defer_singleshot_exception_silenced(self, qapp, monkeypatch):
        ctx, handlers = self._init_ctx()
        ctx.current_molecule = _MolWithConf([FakeAtom(0, "C")])
        monkeypatch.setattr(
            _xyz_editor.QTimer, "singleShot", MagicMock(side_effect=RuntimeError("boom"))
        )
        handlers["load"]({"custom_labels": {"0": "C13"}})  # must not raise


# ===========================================================================
# run(mw) legacy entry point
# ===========================================================================


class TestXYZEditorRun:
    def test_run_without_context_is_noop(self, qapp):
        _xyz_editor.PLUGIN_CONTEXT = None
        _xyz_editor.run(MagicMock())  # must not raise

    def test_run_unwraps_host_and_creates_window(self, qapp):
        ctx = _ctx_no_mol()
        ctx.get_window.return_value = None
        _xyz_editor.PLUGIN_CONTEXT = ctx
        mw = MagicMock()
        mw.host = MagicMock()
        _xyz_editor.run(mw)
        new_win = ctx.register_window.call_args_list[-1][0][1]
        assert new_win is not None
        assert new_win.isVisible()
        new_win.update_timer.stop()
        new_win.destroy()

    def test_run_reuses_existing_window(self, qapp):
        ctx = _ctx_no_mol()
        win = MagicMock()
        ctx.get_window.return_value = win
        _xyz_editor.PLUGIN_CONTEXT = ctx
        _xyz_editor.run(MagicMock())
        win.show.assert_called_once()
        win.raise_.assert_called_once()
        win.activateWindow.assert_called_once()


# ===========================================================================
# _ClickFilter — real QMouseEvent objects
# ===========================================================================


class TestXYZClickFilterRealEvents:
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

    def test_click_within_threshold_invokes_callback(self, qapp):
        calls = []
        f = _xyz_editor._ClickFilter(lambda x, y, obj, mods: calls.append((x, y)))
        press, release = self._events((10, 10), (12, 13))
        obj = MagicMock()
        f.eventFilter(obj, press)
        f.eventFilter(obj, release)
        assert calls == [(12, 13)]

    def test_drag_beyond_threshold_ignored(self, qapp):
        calls = []
        f = _xyz_editor._ClickFilter(lambda x, y, obj, mods: calls.append((x, y)))
        press, release = self._events((10, 10), (30, 30))
        f.eventFilter(MagicMock(), press)
        f.eventFilter(MagicMock(), release)
        assert calls == []

    def test_never_consumes_events(self, qapp):
        f = _xyz_editor._ClickFilter(lambda *a: None)
        press, release = self._events((0, 0), (0, 0))
        assert f.eventFilter(MagicMock(), press) is False
        assert f.eventFilter(MagicMock(), release) is False

    def test_release_without_press_ignored(self, qapp):
        calls = []
        f = _xyz_editor._ClickFilter(lambda x, y, obj, mods: calls.append((x, y)))
        _, release = self._events((0, 0), (5, 5))
        f.eventFilter(MagicMock(), release)
        assert calls == []


# ===========================================================================
# Plotter-picking install/remove — real interactor, exception branches
# ===========================================================================


class TestXYZPlotterPickingRealInteractor:
    def test_enable_installs_and_disable_removes(self, qapp):
        interactor = MagicMock()
        plotter = MagicMock()
        plotter.interactor = interactor
        ctx = _real_ctx(mol=None)
        ctx.plotter = plotter
        w = _xyzrn.XYZEditorWindow(context=ctx)
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
        w = _xyzrn.XYZEditorWindow(context=ctx)
        assert w._click_filter is None
        w.close()
        w.destroy()

    def test_enable_exception_silenced(self, qapp):
        plotter = MagicMock()
        type(plotter).interactor = property(lambda self: (_ for _ in ()).throw(RuntimeError("boom")))
        ctx = _real_ctx(mol=None)
        ctx.plotter = plotter
        w = _xyzrn.XYZEditorWindow(context=ctx)  # must not raise
        assert w._click_filter is None
        w.close()
        w.destroy()

    def test_disable_exception_silenced(self, qapp):
        ctx = _real_ctx(mol=None)
        ctx.plotter = None
        w = _xyzrn.XYZEditorWindow(context=ctx)
        w._click_filter = MagicMock()
        # context.plotter access itself won't raise (MagicMock ctx), but
        # removeEventFilter on a bogus interactor will.
        bad_plotter = MagicMock()
        bad_plotter.interactor.removeEventFilter.side_effect = RuntimeError("boom")
        w.context.plotter = bad_plotter
        w._disable_plotter_picking()  # must not raise
        assert w._click_filter is None
        w.destroy()


# ===========================================================================
# _on_plotter_click — real molecule + fragment/row-map/select_rows helpers
# ===========================================================================


class _FakePicker:
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


class TestXYZOnPlotterClick:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        ctx.plotter = MagicMock()
        w = _xyzrn.XYZEditorWindow(context=ctx)
        mw = MagicMock()
        mw.view_3d_manager = MagicMock()
        ctx.get_main_window.return_value = mw
        yield w
        w.update_timer.stop()
        w.destroy()

    def _widget(self):
        from PyQt6.QtWidgets import QWidget

        widget = QWidget()
        widget.resize(400, 300)
        return widget

    def test_no_molecule_returns_early(self, win, monkeypatch):
        widget = self._widget()
        win.context.current_mol = None
        monkeypatch.setattr(_vtk, "vtkCellPicker", lambda: _FakePicker(actor="x"))
        win._on_plotter_click(10, 10, widget, None)  # must not raise

    def test_no_view3d_manager_returns_early(self, win, monkeypatch):
        widget = self._widget()
        win.context.get_main_window.return_value.view_3d_manager = None
        monkeypatch.setattr(_vtk, "vtkCellPicker", lambda: _FakePicker(actor="x"))
        win._on_plotter_click(10, 10, widget, None)  # must not raise

    def test_wrong_actor_returns_early(self, win, monkeypatch):
        widget = self._widget()
        win.context.get_main_window.return_value.view_3d_manager.atom_actor = "atom-actor"
        monkeypatch.setattr(
            _vtk, "vtkCellPicker", lambda: _FakePicker(actor="other-actor")
        )
        win._on_plotter_click(10, 10, widget, None)
        assert len(win.table.selectedIndexes()) == 0

    def test_click_selects_nearest_atom(self, win, monkeypatch):
        from PyQt6.QtCore import Qt

        widget = self._widget()
        atom_actor = win.context.get_main_window().view_3d_manager.atom_actor
        monkeypatch.setattr(
            _vtk,
            "vtkCellPicker",
            lambda: _FakePicker(actor=atom_actor, pos=(1.5, 0.0, 0.0)),
        )
        win._on_plotter_click(10, 10, widget, Qt.KeyboardModifier.NoModifier)
        rows = {i.row() for i in win.table.selectedIndexes()}
        assert rows == {1}

    def test_ctrl_click_adds_fragment_selection_whole_mol_mode(self, win, monkeypatch):
        from PyQt6.QtCore import Qt

        widget = self._widget()
        win.whole_mol_cb.setChecked(True)
        atom_actor = win.context.get_main_window().view_3d_manager.atom_actor
        monkeypatch.setattr(
            _vtk,
            "vtkCellPicker",
            lambda: _FakePicker(actor=atom_actor, pos=(0.0, 0.0, 0.0)),
        )
        win._on_plotter_click(10, 10, widget, Qt.KeyboardModifier.NoModifier)
        rows = {i.row() for i in win.table.selectedIndexes()}
        assert rows == {0, 1, 2}  # whole connected fragment


class TestXYZFragmentHelpersRealChem:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _xyzrn.XYZEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_fragment_atom_indices_whole_molecule(self, win):
        mol = win.context.current_molecule
        assert win._fragment_atom_indices(mol, 0) == {0, 1, 2}

    def test_atom_index_to_row_map(self, win):
        assert win._atom_index_to_row_map() == {0: 0, 1: 1, 2: 2}

    def test_select_rows_ctrl_toggle(self, win):
        win._select_rows({0, 1}, False, 0)
        rows = {i.row() for i in win.table.selectedIndexes()}
        assert rows == {0, 1}
        win._select_rows({0}, True, 0)  # ctrl-click already-selected row -> deselect
        rows = {i.row() for i in win.table.selectedIndexes()}
        assert rows == {1}


# ===========================================================================
# closeEvent / get_mol_signature / check_molecule_update exception branches
# ===========================================================================


class TestXYZExceptionBranches:
    def test_close_event_timer_stop_exception_silenced(self, qapp):
        ctx = _ctx_no_mol()
        w = _xyz_editor.XYZEditorWindow(context=ctx)
        bad_timer = MagicMock()
        bad_timer.isActive.side_effect = RuntimeError("boom")
        w.update_timer = bad_timer
        w.close()  # must not raise
        w.destroy()

    def test_get_mol_signature_exception_returns_none(self, qapp):
        ctx = _ctx_no_mol()
        w = _xyz_editor.XYZEditorWindow(context=ctx)
        bad_mol = MagicMock()
        bad_mol.GetNumAtoms.side_effect = RuntimeError("boom")
        assert w.get_mol_signature(bad_mol) is None
        w.update_timer.stop()
        w.destroy()

    def test_check_molecule_update_exception_silenced(self, qapp):
        ctx = _ctx_no_mol()
        w = _xyz_editor.XYZEditorWindow(context=ctx)
        type(ctx).current_molecule = property(
            lambda self: (_ for _ in ()).throw(RuntimeError("boom"))
        )
        w.check_molecule_update()  # must not raise
        w.update_timer.stop()
        w.destroy()
        del type(ctx).current_molecule


# ===========================================================================
# load_molecule: dummyLabel fallback
# ===========================================================================


class TestXYZLoadMoleculeDummyLabel:
    def test_dummy_label_used_when_no_custom_symbol(self, qapp):
        class _DummyAtom(FakeAtom):
            def HasProp(self, name):
                return name == "dummyLabel"

            def GetProp(self, name):
                return "Xx"

        w, ctx, mol = _win_with_mol(("C",))
        ctx.current_molecule = _MolWithConf([_DummyAtom(0, "C")])
        w.load_molecule()
        assert w.table.item(0, 1).text() == "Xx"
        _teardown(w)


# ===========================================================================
# remove_selected_rows: no-selection guard + noncontiguous batches
# ===========================================================================


class TestXYZRemoveSelectedRowsReal:
    def test_no_selection_is_noop(self, qapp):
        w, ctx, mol = _win_with_mol(("C", "O"))
        w.remove_selected_rows()
        assert w.table.rowCount() == 2
        _teardown(w)

    def test_noncontiguous_rows_removed(self, qapp):
        from PyQt6.QtWidgets import QTableWidgetSelectionRange

        w, ctx, mol = _win_with_mol(("C", "O", "N", "H"))
        last_col = w.table.columnCount() - 1
        w.table.setRangeSelected(
            QTableWidgetSelectionRange(0, 0, 0, last_col), True
        )
        w.table.setRangeSelected(
            QTableWidgetSelectionRange(2, 0, 2, last_col), True
        )
        w.remove_selected_rows()
        assert w.table.rowCount() == 2
        remaining = {w.table.item(r, 1).text() for r in range(w.table.rowCount())}
        assert remaining == {"O", "H"}
        _teardown(w)


class TestXYZSelectedAtomIndicesInvalidText:
    def test_non_integer_index_skipped(self, qapp):
        w, ctx, mol = _win_with_mol(("C",))
        w.table.item(0, 0).setText("not_an_int")
        w.table.selectRow(0)
        assert w._selected_atom_indices() == set()
        _teardown(w)


# ===========================================================================
# duplicate_atoms / adjust_hydrogens / apply_changes — real RDKit flows
# ===========================================================================


class TestXYZDuplicateAtomsRealChem:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _xyzrn.XYZEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_duplicate_whole_molecule(self, win):
        win.duplicate_atoms()
        new_mol = win.context.current_molecule
        assert new_mol.GetNumAtoms() == 6
        win.context.push_undo_checkpoint.assert_called_once()
        msg = win.context.show_status_message.call_args[0][0]
        assert "Duplicated 3 atom" in msg

    def test_duplicate_selected_subset(self, win):
        from PyQt6.QtWidgets import QTableWidgetSelectionRange

        last_col = win.table.columnCount() - 1
        win.table.setRangeSelected(
            QTableWidgetSelectionRange(0, 0, 1, last_col), True
        )
        win.duplicate_atoms()
        new_mol = win.context.current_molecule
        assert new_mol.GetNumAtoms() == 5

    def test_duplicate_custom_offset(self, win):
        for spin, val in zip(win.dup_offset, (2.0, 0.0, 0.0)):
            spin.setValue(val)
        win.duplicate_atoms()
        new_mol = win.context.current_molecule
        conf = new_mol.GetConformer()
        p = conf.GetAtomPosition(3)
        assert p.x == pytest.approx(2.0)

    def test_duplicate_exception_shows_critical(self, win, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _xyzrn.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.duplicate_atoms()
        assert calls
        assert "Failed to duplicate" in calls[0][2]


class TestXYZAdjustHydrogensRealChem:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _xyzrn.XYZEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_add_missing_hydrogens_whole_molecule(self, win):
        win.adjust_hydrogens()
        new_mol = win.context.current_molecule
        assert new_mol.GetNumAtoms() > 3
        win.context.push_undo_checkpoint.assert_called_once()

    def test_selection_passes_only_on_atoms(self, win):
        win.table.selectRow(0)
        win.adjust_hydrogens()
        new_mol = win.context.current_molecule
        assert new_mol.GetNumAtoms() > 3

    def test_exception_shows_critical(self, win, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        monkeypatch.setattr(
            _xyzrn.Chem, "RWMol", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.adjust_hydrogens()
        assert calls
        assert "Failed to adjust hydrogens" in calls[0][2]

    def test_excess_hydrogen_removed(self, win):
        # Build a methane-like overvalent carbon with 5 H's, then adjust down to 4.
        rw = _Chem.RWMol()
        c = rw.AddAtom(_Chem.Atom(6))
        h_idxs = [rw.AddAtom(_Chem.Atom(1)) for _ in range(5)]
        for h in h_idxs:
            rw.AddBond(c, h, _Chem.BondType.SINGLE)
        conf = _Chem.Conformer(6)
        coords = [(0, 0, 0)] + [(1.0, i * 0.5, 0.0) for i in range(5)]
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, _Point3D(x, y, z))
        rw.AddConformer(conf, assignId=True)
        rw.UpdatePropertyCache(strict=False)
        mol = rw.GetMol()
        win.context.current_molecule = mol
        win.adjust_hydrogens()
        new_mol = win.context.current_molecule
        assert new_mol.GetNumAtoms() == 5


class TestXYZApplyChangesRealChem:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        w = _xyzrn.XYZEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_edit_coordinate_and_apply(self, win):
        win.table.blockSignals(True)
        win.table.item(0, 2).setText("5.0")
        win.table.blockSignals(False)
        win.apply_changes()
        new_mol = win.context.current_molecule
        conf = new_mol.GetConformer()
        assert conf.GetAtomPosition(0).x == pytest.approx(5.0)
        win.context.push_undo_checkpoint.assert_called_once()

    def test_unrecognized_symbol_becomes_dummy(self, win):
        win.table.blockSignals(True)
        win.table.item(0, 1).setText("Xx")
        win.table.blockSignals(False)
        win.apply_changes()
        new_mol = win.context.current_molecule
        assert new_mol.GetAtomWithIdx(0).HasProp("custom_symbol")

    def test_prefix_match_ghost_suffix(self, win):
        win.table.blockSignals(True)
        win.table.item(0, 1).setText("Ag*")
        win.table.blockSignals(False)
        win.apply_changes()
        new_mol = win.context.current_molecule
        assert new_mol.GetAtomWithIdx(0).GetAtomicNum() == 47
        assert new_mol.GetAtomWithIdx(0).GetProp("custom_symbol") == "Ag*"

    def test_malformed_coordinate_aborts(self, win):
        win.table.blockSignals(True)
        win.table.item(0, 2).setText("abc")
        win.table.blockSignals(False)
        before = win.context.current_molecule
        win.apply_changes()
        assert win.context.current_molecule is before
        win.context.push_undo_checkpoint.assert_not_called()

    def test_bonds_preserved_for_surviving_atoms(self, win):
        win.apply_changes()
        new_mol = win.context.current_molecule
        assert new_mol.GetBondBetweenAtoms(0, 1) is not None
        assert new_mol.GetBondBetweenAtoms(1, 2) is not None

    def test_apply_after_row_removed_drops_bond(self, win):
        win.table.selectRow(2)
        win.remove_selected_rows()
        win.apply_changes()
        new_mol = win.context.current_molecule
        assert new_mol.GetNumAtoms() == 2

    def test_no_molecule_still_builds_from_table(self, qapp):
        ctx = _real_ctx(mol=None)
        w = _xyzrn.XYZEditorWindow(context=ctx)
        w.add_atom_row()
        w.apply_changes()
        new_mol = w.context.current_molecule
        assert new_mol.GetNumAtoms() == 1
        w.update_timer.stop()
        w.destroy()

    def test_exception_shows_critical(self, win, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        calls = []
        monkeypatch.setattr(
            QMessageBox, "critical", lambda *a, **k: calls.append(a) or None
        )
        # Conformer() is inside apply_changes' try block (RWMol() itself is
        # constructed before the try, so patching it would raise uncaught).
        monkeypatch.setattr(
            _xyzrn.Chem, "Conformer", MagicMock(side_effect=RuntimeError("boom"))
        )
        win.apply_changes()
        assert calls
        assert "Failed to apply changes" in calls[0][2]


# ===========================================================================
# highlight_selected_atoms — real pyvista PolyData/glyph path
# ===========================================================================


class TestXYZHighlightSelectedAtomsRealChem:
    @pytest.fixture
    def win(self, qapp):
        ctx = _real_ctx(mol=_real_mol())
        ctx.plotter = MagicMock()
        w = _xyzrn.XYZEditorWindow(context=ctx)
        yield w
        w.update_timer.stop()
        w.destroy()

    def test_no_selection_removes_actor(self, win):
        win.table.clearSelection()
        win.highlight_selected_atoms()
        win.context.plotter.remove_actor.assert_called_with("xyz_selection")
        win.context.plotter.render.assert_called()

    def test_selection_adds_sphere_mesh(self, win):
        win.table.selectRow(0)  # itemSelectionChanged already ran highlight once
        win.context.plotter.reset_mock()
        win.highlight_selected_atoms()
        win.context.plotter.add_mesh.assert_called_once()
        kwargs = win.context.plotter.add_mesh.call_args[1]
        assert kwargs["name"] == "xyz_selection"

    def test_unknown_symbol_uses_default_radius(self, win):
        win.table.item(0, 1).setText("Xx")
        win.table.selectRow(0)
        win.context.plotter.reset_mock()
        win.highlight_selected_atoms()  # must not raise
        win.context.plotter.add_mesh.assert_called_once()

    def test_malformed_row_skipped_falls_to_no_points(self, win):
        win.table.blockSignals(True)
        win.table.item(0, 2).setText("not_a_float")
        win.table.blockSignals(False)
        win.table.selectRow(0)
        win.highlight_selected_atoms()  # must not raise
        win.context.plotter.remove_actor.assert_called_with("xyz_selection")

    def test_camera_position_restored_after_highlight(self, win):
        win.context.plotter.camera_position = "cam1"
        win.table.selectRow(0)
        win.highlight_selected_atoms()
        assert win.context.plotter.camera_position == "cam1"

    def test_on_item_changed_updates_highlight_for_selected_row(self, win):
        win.table.selectRow(0)
        win.context.plotter.reset_mock()
        item = win.table.item(0, 2)
        item.setText("9.0")
        win.context.plotter.add_mesh.assert_called()

    def test_on_item_changed_ignores_unselected_row(self, win):
        win.context.plotter.reset_mock()
        item = win.table.item(0, 2)
        item.setText("9.0")
        win.context.plotter.add_mesh.assert_not_called()
