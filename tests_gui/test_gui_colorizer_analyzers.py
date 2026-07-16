"""
Headless GUI tests for:
  - Bond Colorizer        → BondColorizerWindow
  - Vector Viewer         → VectorViewerPlugin
  - Gaussian MO Analyzer  → OrbitalWidget (package plugin)
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

BOND_PATH = PLUGINS_DIR / "Bond_Colorizer" / "bond_colorizer.py"
VECTOR_PATH = PLUGINS_DIR / "Vector_Viewer" / "vector_viewer.py"
MO_PKG_DIR = PLUGINS_DIR / "Gaussian_MO_Analyzer"

with mock_chemistry_imports():
    _bond = load_plugin_for_gui(BOND_PATH)
    _vector = load_plugin_for_gui(VECTOR_PATH)
    sys.path.insert(0, str(MO_PKG_DIR))
    try:
        _mo = importlib.import_module("gaussian_fchk_mo_analyzer")
    finally:
        sys.path.remove(str(MO_PKG_DIR))


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


# ===========================================================================
# VectorViewerPlugin  (Vector Viewer)
# ===========================================================================


def _vector_context():
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


class TestVectorViewerPlugin:
    @pytest.fixture
    def win(self, qapp):
        ctx = _vector_context()
        w = _vector.VectorViewerPlugin(ctx)
        yield w
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "Vector Viewer"

    def test_vector_input_placeholder(self, win):
        assert "vector" in win.vec_input.placeholderText().lower()

    def test_default_scale(self, win):
        assert win.scale_spin.value() == 1.0

    def test_default_resolution(self, win):
        assert win.res_spin.value() == 20

    def test_default_opacity(self, win):
        assert win.opacity_spin.value() == 0.5

    def test_reverse_unchecked_by_default(self, win):
        assert not win.reverse_chk.isChecked()

    def test_transparent_export_checked_by_default(self, win):
        assert win.trans_chk.isChecked()

    def test_color_button_shows_green_hex(self, win):
        assert win.color_btn.text() == "#008000"

    def test_registered_as_main_panel(self, win):
        win.context.register_window.assert_called_with("main_panel", win)

    def test_update_with_no_plotter_is_noop(self, win):
        win.context.plotter = None
        win.vec_input.setText("1 0 0")
        win.update_visualization()  # must not raise

    def test_update_with_empty_input_skips_plotter(self, win):
        plotter = MagicMock()
        win.context.plotter = plotter
        win.vec_input.setText("")
        win.update_visualization()
        plotter.add_mesh.assert_not_called()

    def test_update_with_invalid_input_skips_plotter(self, win):
        plotter = MagicMock()
        win.context.plotter = plotter
        win.vec_input.setText("1 2")  # fewer than 3 components
        win.update_visualization()
        plotter.add_mesh.assert_not_called()

    def test_close_removes_actor(self, win):
        plotter = MagicMock()
        win.context.plotter = plotter
        sentinel = object()
        win.vis_actor = sentinel
        win.close()
        plotter.remove_actor.assert_called_once_with(sentinel)
        assert win.vis_actor is None

    def test_initialize_sets_launch_fn(self, qapp):
        ctx = _vector_context()
        _vector._launch_fn = None
        _vector.initialize(ctx)
        assert callable(_vector._launch_fn)
        ctx.show_status_message.assert_called_once()
        _vector._launch_fn = None


# ===========================================================================
# Gaussian MO Analyzer  (package plugin: OrbitalWidget + initialize)
# ===========================================================================


class TestMOAnalyzerInitialize:
    def test_registers_fchk_openers_and_drop_handler(self, qapp):
        ctx = MagicMock()
        _mo.initialize(ctx)
        exts = [c.args[0] for c in ctx.register_file_opener.call_args_list]
        assert exts == [".fchk", ".fck", ".fch"]
        ctx.register_drop_handler.assert_called_once()

    def test_drop_handler_accepts_only_fchk(self, qapp):
        ctx = MagicMock()
        _mo.initialize(ctx)
        handler = ctx.register_drop_handler.call_args.args[0]
        with patch.object(_mo, "OrbitalWidget") as ow:
            assert handler("mol.fchk") is True
            assert handler("mol.xyz") is False
            ow.assert_called_once()


class TestOrbitalWidget:
    @pytest.fixture
    def win(self, qapp, tmp_path):
        # Empty fchk → reader.data empty; basis prep patched out so the
        # widget builds with n_basis=0 (the "Error reading MOs" path).
        fchk = tmp_path / "water.fchk"
        fchk.write_text("dummy fchk\n", encoding="utf-8")
        ctx = MagicMock()
        with patch.object(
            _mo.analyzer.BasisSetEngine, "_prepare_basis_set", lambda self: None
        ):
            w = _mo.OrbitalWidget(None, ctx, str(fchk))
        yield w
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title_includes_filename(self, win):
        assert "Gaussian MO Analyzer" in win.windowTitle()
        assert "water.fchk" in win.windowTitle()

    def test_default_grid_points(self, win):
        assert win.spin_points.value() == 40

    def test_default_margin(self, win):
        assert win.spin_gap.value() == 4.0

    def test_default_isovalue(self, win):
        assert win.spin_iso.value() == 0.02

    def test_default_opacity(self, win):
        assert win.spin_opacity.value() == 0.4

    def test_smooth_shading_checked_by_default(self, win):
        assert win.check_smooth.isChecked()

    def test_progress_bar_hidden_initially(self, win):
        assert win.pbar.isHidden()

    def test_list_shows_error_without_mo_data(self, win):
        assert win.list_widget.count() == 1
        assert win.list_widget.item(0).text() == "Error reading MOs"

    def test_generate_with_no_selection_keeps_button_enabled(self, win):
        win.list_widget.clearSelection()
        win.generate_selected()
        assert win.btn_gen.isEnabled()

    def test_iso_change_without_cube_is_noop(self, win):
        win.last_cube_path = None
        win.on_iso_changed(0.05)  # must not raise

    def test_cube_path_layout(self, win):
        path = win.get_cube_path(7)
        assert path.endswith("water_MO7.cube")
        assert "water_cubes" in path
