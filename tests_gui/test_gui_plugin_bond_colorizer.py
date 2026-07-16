"""
Headless GUI tests for the Bond Colorizer plugin.

Covers: BondColorizerWindow, initialize().
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

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


# ===========================================================================
# BondColorizerWindow  (Bond Colorizer)
# ===========================================================================


class TestBondColorizerWindow:
    @pytest.fixture
    def win(self, qapp):
        ctx = _bond_context()
        w = _bond.BondColorizerWindow(ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

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
