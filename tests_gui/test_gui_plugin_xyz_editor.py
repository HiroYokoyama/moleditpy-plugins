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
