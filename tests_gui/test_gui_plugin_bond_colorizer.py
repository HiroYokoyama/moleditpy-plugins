"""
Headless GUI tests for the Bond Colorizer plugin.

Covers: BondColorizerWindow, initialize().
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import ANY, MagicMock, patch

import pytest
from PyQt6.QtGui import QCloseEvent, QColor
from PyQt6.QtWidgets import QWidget

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

BOND_PATH = PLUGINS_DIR / "Bond_Colorizer" / "bond_colorizer.py"

with mock_chemistry_imports():
    _bond = load_plugin_for_gui(BOND_PATH)


class FakeBond:
    def __init__(self, idx, a1, a2):
        self._idx, self._a1, self._a2 = idx, a1, a2

    def GetIdx(self):
        return self._idx

    def GetBeginAtomIdx(self):
        return self._a1

    def GetEndAtomIdx(self):
        return self._a2


def _bond_context(mol=None):
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = mol
    return ctx


def _mw_with_edit3d(measurement_mode=False, edit_mode=False):
    """A real QWidget (valid as a QDialog parent) carrying mocked managers,
    for the enter/restore select-mode paths."""
    mw = QWidget()
    mw.edit_3d_manager = MagicMock()
    mw.edit_3d_manager.measurement_mode = measurement_mode
    mw.edit_3d_manager.is_3d_edit_mode = edit_mode
    mw.init_manager = MagicMock()
    mw.ui_manager = MagicMock()
    return mw


# ===========================================================================
# BondColorizerWindow  (Bond Colorizer)
# ===========================================================================


@pytest.fixture
def win(qapp):
    ctx = _bond_context()
    w = _bond.BondColorizerWindow(ctx)
    yield w
    w.sel_timer.stop()
    w.destroy()


class TestBondColorizerWindow:
    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "Bond Colorizer"

    def test_default_color_is_red(self, win):
        assert win.current_color.name() == "#ff0000"

    def test_bond_ids_placeholder(self, win):
        assert "Bond IDs" in win.le_bond_ids.placeholderText()

    def test_atom_pairs_placeholder(self, win):
        assert "Atom pairs" in win.le_atom_pairs.placeholderText()

    def test_selection_timer_running(self, win):
        assert win.sel_timer.isActive()
        assert win.sel_timer.interval() == 200

    def test_registered_as_main_panel(self, win):
        win.context.register_window.assert_called_with("main_panel", win)

    def test_apply_with_no_selection_warns(self, win):
        with patch.object(_bond, "QMessageBox") as mb:
            win.apply_color()
        mb.warning.assert_called_once()
        assert "No bonds selected" in mb.warning.call_args.args[2]

    def test_apply_with_no_molecule_warns(self, win):
        win.le_bond_ids.setText("0")
        win.context.current_molecule = None
        with patch.object(_bond, "QMessageBox") as mb:
            win.apply_color()
        assert "No molecule" in mb.warning.call_args.args[2]

    def test_apply_invalid_bond_id_format_warns(self, win):
        win.le_bond_ids.setText("abc")
        win.context.current_molecule = MagicMock()
        with patch.object(_bond, "QMessageBox") as mb:
            win.apply_color()
        assert "Invalid bond ID" in mb.warning.call_args.args[2]

    def test_apply_valid_bond_ids_colors_bonds(self, win):
        mol = MagicMock()
        mol.GetNumBonds.return_value = 5
        win.context.current_molecule = mol
        win.le_bond_ids.setText("0, 2")
        controller = win.context.get_3d_controller.return_value
        win.apply_color()
        assert controller.set_bond_color.call_args_list == [
            ((0, "#ff0000"),),
            ((2, "#ff0000"),),
        ]
        win.context.refresh_3d_view.assert_called_once()

    def test_apply_atom_pair_resolves_bond_index(self, win):
        mol = MagicMock()
        mol.GetNumBonds.return_value = 5
        bond = MagicMock()
        bond.GetIdx.return_value = 3
        mol.GetBondBetweenAtoms.return_value = bond
        win.context.current_molecule = mol
        win.le_atom_pairs.setText("0-1")
        controller = win.context.get_3d_controller.return_value
        win.apply_color()
        controller.set_bond_color.assert_called_once_with(3, "#ff0000")

    def test_apply_nonexistent_pair_reports_skipped(self, win):
        mol = MagicMock()
        mol.GetNumBonds.return_value = 5
        mol.GetBondBetweenAtoms.return_value = None
        win.context.current_molecule = mol
        win.le_atom_pairs.setText("7-8")
        with patch.object(_bond, "QMessageBox") as mb:
            win.apply_color()
        assert "Skipped invalid" in mb.warning.call_args.args[2]

    def test_reset_colors_clears_fields_and_overrides(self, win):
        mol = MagicMock()
        mol.GetNumBonds.return_value = 2
        win.context.current_molecule = mol
        win.le_bond_ids.setText("0")
        controller = win.context.get_3d_controller.return_value
        win.reset_colors()
        assert controller.set_bond_color.call_args_list == [
            ((0, None),),
            ((1, None),),
        ]
        assert win.le_bond_ids.text() == ""

    def test_selection_sync_fills_atom_pairs(self, win):
        mol = MagicMock()
        mol.GetBonds.return_value = [FakeBond(0, 0, 1), FakeBond(1, 1, 2)]
        win.context.current_molecule = mol
        win.context.get_selected_atom_indices.return_value = {0, 1}
        win.get_selection_from_viewer()
        assert win.le_atom_pairs.text() == "0-1"

    def test_close_stops_timer_and_unregisters(self, win):
        win.close()
        assert not win.sel_timer.isActive()
        win.context.register_window.assert_called_with("main_panel", None)


class TestBondColorizerInitialize:
    def test_registers_all_handlers(self, qapp):
        ctx = MagicMock()
        _bond.initialize(ctx)
        ctx.register_save_handler.assert_called_once()
        ctx.register_load_handler.assert_called_once()
        ctx.register_document_reset_handler.assert_called_once()

    def test_load_handler_restores_colors(self, qapp):
        ctx = MagicMock()
        _bond.initialize(ctx)
        load_handler = ctx.register_load_handler.call_args.args[0]
        controller = ctx.get_3d_controller.return_value
        load_handler({"bond_colors": {"2": "#00ff00"}})
        controller.set_bond_color.assert_called_once_with(2, "#00ff00")
        ctx.refresh_3d_view.assert_called_once()

    def test_save_handler_returns_empty_without_overrides(self, qapp):
        ctx = MagicMock()
        mw = MagicMock(spec=[])  # no view_3d_manager attribute
        ctx.get_main_window.return_value = mw
        _bond.initialize(ctx)
        save_handler = ctx.register_save_handler.call_args.args[0]
        assert save_handler() == {}


class TestEnterSelectMode:
    def test_forces_measurement_mode_when_not_active(self, qapp):
        ctx = _bond_context()
        mw = _mw_with_edit3d(measurement_mode=False)
        ctx.get_main_window.return_value = mw
        win = _bond.BondColorizerWindow(ctx)
        try:
            assert win._forced_measurement_mode is True
            mw.init_manager.measurement_action.setChecked.assert_called_with(True)
            mw.edit_3d_manager.toggle_measurement_mode.assert_called_with(True)
        finally:
            win.sel_timer.stop()
            win.destroy()

    def test_skips_toggle_when_already_active(self, qapp):
        ctx = _bond_context()
        mw = _mw_with_edit3d(measurement_mode=True)
        ctx.get_main_window.return_value = mw
        win = _bond.BondColorizerWindow(ctx)
        try:
            assert win._forced_measurement_mode is False
            mw.edit_3d_manager.toggle_measurement_mode.assert_not_called()
        finally:
            win.sel_timer.stop()
            win.destroy()

    def test_exception_is_silenced(self, qapp):
        ctx = _bond_context()
        mw = _mw_with_edit3d(measurement_mode=False)
        mw.init_manager.measurement_action.setChecked.side_effect = RuntimeError("boom")
        ctx.get_main_window.return_value = mw
        win = _bond.BondColorizerWindow(ctx)  # must not raise
        try:
            assert win is not None
        finally:
            win.sel_timer.stop()
            win.destroy()


class TestRestoreSelectMode:
    def test_close_restores_full_path(self, qapp):
        ctx = _bond_context()
        mw = _mw_with_edit3d()
        mw.edit_3d_manager.selected_atoms_3d = {1, 2, 3}
        ctx.get_main_window.return_value = mw
        win = _bond.BondColorizerWindow(ctx)
        win.close()
        mw.edit_3d_manager.clear_measurement_selection.assert_called_once()
        mw.edit_3d_manager.update_3d_selection_display.assert_called_once()
        assert mw.edit_3d_manager.selected_atoms_3d == set()
        mw.edit_3d_manager.toggle_measurement_mode.assert_called_with(
            win._restore_measurement_mode
        )
        mw.ui_manager.toggle_3d_edit_mode.assert_called_with(win._restore_edit_mode)

    def test_clear_selection_exception_is_silenced(self, qapp):
        ctx = _bond_context()
        mw = _mw_with_edit3d()
        mw.edit_3d_manager.selected_atoms_3d = {1}
        mw.edit_3d_manager.clear_measurement_selection.side_effect = RuntimeError("boom")
        ctx.get_main_window.return_value = mw
        win = _bond.BondColorizerWindow(ctx)
        win.close()  # must not raise

    def test_toggle_exception_is_silenced(self, qapp):
        ctx = _bond_context()
        mw = _mw_with_edit3d()
        mw.edit_3d_manager.toggle_measurement_mode.side_effect = RuntimeError("boom")
        ctx.get_main_window.return_value = mw
        win = _bond.BondColorizerWindow(ctx)
        win.close()  # must not raise


class TestSelectionSyncRealMainWindow:
    def test_merges_3d_and_measurement_picks(self, qapp):
        ctx = _bond_context()
        mw = _mw_with_edit3d()
        mw.edit_3d_manager.selected_atoms_3d = {2}
        mw.edit_3d_manager.selected_atoms_for_measurement = [3]
        ctx.get_main_window.return_value = mw
        ctx.get_selected_atom_indices.return_value = set()
        win = _bond.BondColorizerWindow(ctx)
        try:
            mol = MagicMock()
            mol.GetBonds.return_value = [FakeBond(0, 2, 3)]
            win.context.current_molecule = mol
            win.get_selection_from_viewer()
            assert win.le_atom_pairs.text() == "2-3"
        finally:
            win.sel_timer.stop()
            win.destroy()

    def test_no_molecule_with_real_mw_returns_early(self, qapp):
        ctx = _bond_context()
        mw = _mw_with_edit3d()
        ctx.get_main_window.return_value = mw
        win = _bond.BondColorizerWindow(ctx)
        try:
            win.context.current_molecule = None
            win.get_selection_from_viewer()  # must not raise
            assert win.le_atom_pairs.text() == ""
        finally:
            win.sel_timer.stop()
            win.destroy()

    def test_exception_is_silenced(self, win):
        win.context.current_molecule = MagicMock()
        win.context.current_molecule.GetBonds.side_effect = RuntimeError("boom")
        win.get_selection_from_viewer()  # must not raise


class TestAutoUpdateSelectionRealWidget:
    def test_runs_when_neither_field_focused(self, win):
        win._auto_update_selection()  # must not raise


class TestChooseColor:
    def test_valid_color_updates_button_and_current_color(self, win):
        with patch.object(_bond, "QColorDialog") as qcd:
            qcd.getColor.return_value = QColor(10, 20, 30)
            win.choose_color()
        assert win.current_color.name() == "#0a141e"

    def test_invalid_color_leaves_current_color_unchanged(self, win):
        before = win.current_color.name()
        with patch.object(_bond, "QColorDialog") as qcd:
            qcd.getColor.return_value = QColor()  # invalid (default-constructed)
            win.choose_color()
        assert win.current_color.name() == before


class TestApplyColorExtraPaths:
    def test_out_of_range_bond_id_skipped_and_warned(self, win):
        mol = MagicMock()
        mol.GetNumBonds.return_value = 2
        win.context.current_molecule = mol
        win.le_bond_ids.setText("0, 50")
        controller = win.context.get_3d_controller.return_value
        with patch.object(_bond, "QMessageBox") as mb:
            win.apply_color()
        controller.set_bond_color.assert_called_once_with(0, "#ff0000")
        assert "Skipped invalid" in mb.warning.call_args.args[2]

    def test_invalid_atom_pair_token_shape_warns(self, win):
        win.context.current_molecule = MagicMock()
        win.le_atom_pairs.setText("1-2-3")
        with patch.object(_bond, "QMessageBox") as mb:
            win.apply_color()
        assert "Invalid atom-pair format" in mb.warning.call_args.args[2]

    def test_invalid_atom_pair_non_int_warns(self, win):
        win.context.current_molecule = MagicMock()
        win.le_atom_pairs.setText("a-b")
        with patch.object(_bond, "QMessageBox") as mb:
            win.apply_color()
        assert "Invalid atom-pair format" in mb.warning.call_args.args[2]

    def test_controller_exception_reports_critical(self, win):
        mol = MagicMock()
        mol.GetNumBonds.return_value = 5
        win.context.current_molecule = mol
        win.le_bond_ids.setText("0")
        controller = win.context.get_3d_controller.return_value
        controller.set_bond_color.side_effect = RuntimeError("boom")
        with patch.object(_bond, "QMessageBox") as mb:
            win.apply_color()
        mb.critical.assert_called_once()
        assert "Failed to apply color" in mb.critical.call_args.args[2]


class TestResetColorsExtraPaths:
    def test_no_molecule_is_noop(self, win):
        win.context.current_molecule = None
        win.reset_colors()  # must not raise
        win.context.get_3d_controller.assert_not_called()

    def test_controller_exception_reports_critical(self, win):
        mol = MagicMock()
        mol.GetNumBonds.return_value = 2
        win.context.current_molecule = mol
        controller = win.context.get_3d_controller.return_value
        controller.set_bond_color.side_effect = RuntimeError("boom")
        with patch.object(_bond, "QMessageBox") as mb:
            win.reset_colors()
        mb.critical.assert_called_once()
        assert "Failed to reset colors" in mb.critical.call_args.args[2]


class TestCloseEventTimerException:
    def test_timer_isactive_exception_is_silenced(self, win):
        win.sel_timer = MagicMock()
        win.sel_timer.isActive.side_effect = RuntimeError("boom")
        win.closeEvent(QCloseEvent())  # must not raise


class TestLaunch:
    def test_creates_new_window_when_none_registered(self, qapp):
        ctx = _bond_context()
        ctx.get_window.return_value = None
        _bond.launch(ctx)
        ctx.register_window.assert_any_call("main_panel", ANY)
        win = ctx.register_window.call_args.args[1]
        win.sel_timer.stop()
        win.destroy()

    def test_reuses_existing_window(self, qapp):
        ctx = _bond_context()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _bond.launch(ctx)
        existing.show.assert_called_once()
        existing.raise_.assert_called_once()
        existing.activateWindow.assert_called_once()


class TestRunLegacyEntryPoint:
    def teardown_method(self):
        _bond.PLUGIN_CONTEXT = None

    def test_without_plugin_manager_is_noop(self, qapp):
        mw = MagicMock(spec=[])  # no host, no plugin_manager
        _bond.run(mw)  # returns before touching PLUGIN_CONTEXT

    def test_launches_via_plugin_context(self, qapp):
        ctx = _bond_context()
        ctx.get_window.return_value = None
        _bond.PLUGIN_CONTEXT = ctx
        mw = MagicMock()
        mw.plugin_manager = MagicMock()
        _bond.run(mw)
        ctx.register_window.assert_any_call("main_panel", ANY)
        win = ctx.register_window.call_args.args[1]
        win.sel_timer.stop()
        win.destroy()
