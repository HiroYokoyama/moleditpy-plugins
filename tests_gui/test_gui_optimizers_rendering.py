"""
Headless GUI tests for optimizer plugins and Advanced Rendering.

Plugins covered (all registry-visible):
  - Step Optimizer              → StepOptimizerDialog
  - All-Trans Optimizer         → _select_torsions, run_plugin guards
  - xTB Optimizer               → XtbOptimizerDialog, XtbWorker
  - Complex Molecule Untangler  → UntanglerDialog, run_plugin
  - Advanced Rendering          → AdvancedGraphicsWidget, HideOnCloseDialog

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_optimizers_rendering.py
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

STEP_OPT_PATH = PLUGINS_DIR / "Step_Optimizer" / "step_optimizer.py"
ALL_TRANS_PATH = PLUGINS_DIR / "All-Trans_Optimizer" / "all-trans_optimizer.py"
XTB_PATH = PLUGINS_DIR / "XTB_Optimizer" / "xtb_optimizer.py"
UNTANGLER_PATH = (
    PLUGINS_DIR / "Complex_Molecule_Untangler" / "complex_molecule_untangler.py"
)
ADV_RENDER_PATH = PLUGINS_DIR / "Advanced_Rendering" / "advanced_rendering.py"

with mock_chemistry_imports():
    _step = load_plugin_for_gui(STEP_OPT_PATH)
    _all_trans = load_plugin_for_gui(ALL_TRANS_PATH)
    _xtb = load_plugin_for_gui(XTB_PATH)
    _untangler = load_plugin_for_gui(UNTANGLER_PATH)
    _adv = load_plugin_for_gui(ADV_RENDER_PATH)

# save_settings' sanitize() does isinstance(x, np.generic); the MagicMock numpy
# would make that a TypeError, so give the module a shim with valid class args.
_adv.np = SimpleNamespace(generic=(), ndarray=())


def _ctx(setting="MMFF_RDKIT") -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.get_setting.return_value = setting
    return ctx


@pytest.fixture
def no_msgbox(monkeypatch):
    """Replace each module's QMessageBox so warnings/info don't block offscreen."""
    boxes = {}
    for mod in (_step, _all_trans, _xtb, _untangler):
        box = MagicMock()
        monkeypatch.setattr(mod, "QMessageBox", box)
        boxes[mod.__name__] = box
    return boxes


# ===========================================================================
# StepOptimizerDialog  (visible plugin: "Step Optimizer")
# ===========================================================================


class TestStepOptimizerDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _step.StepOptimizerDialog(context=_ctx(), parent=None)
        yield d
        d.timer.stop()
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Step Optimizer"

    def test_force_field_options(self, dlg):
        texts = [dlg.combo_ff.itemText(i) for i in range(dlg.combo_ff.count())]
        assert texts == ["MMFF94", "MMFF94s", "UFF"]

    def test_default_ff_from_mmff_setting(self, dlg):
        assert dlg.combo_ff.currentText() == "MMFF94"

    def test_default_ff_from_uff_setting(self, qapp):
        d = _step.StepOptimizerDialog(context=_ctx(setting="UFF_RDKIT"), parent=None)
        assert d.combo_ff.currentText() == "UFF"
        d.timer.stop()
        d.destroy()

    def test_steps_per_frame_defaults(self, dlg):
        assert dlg.spin_steps.value() == 1
        assert (dlg.spin_steps.minimum(), dlg.spin_steps.maximum()) == (1, 200)

    def test_max_move_defaults(self, dlg):
        assert dlg.spin_max_move.value() == pytest.approx(0.5)
        assert dlg.spin_max_move.specialValueText() == "Off"

    def test_initial_labels(self, dlg):
        assert dlg.lbl_step.text() == "Step: 0"
        assert dlg.lbl_energy.text() == "Energy: --"
        assert dlg.lbl_state.text() == "State: Idle"

    def test_initial_button_states(self, dlg):
        assert dlg.btn_start.isEnabled()
        assert not dlg.btn_stop.isEnabled()

    def test_timer_configured_but_idle(self, dlg):
        assert dlg.timer.interval() == 50
        assert not dlg.timer.isActive()
        assert not dlg.running

    def test_start_without_molecule_warns(self, dlg, no_msgbox):
        dlg._start()
        no_msgbox[_step.__name__].warning.assert_called_once()
        assert not dlg.running
        assert not dlg.timer.isActive()

    def test_start_without_conformer_warns(self, dlg, no_msgbox):
        mol = MagicMock()
        mol.GetNumConformers.return_value = 0
        dlg.context.current_mol = mol
        dlg._start()
        msg = no_msgbox[_step.__name__].warning.call_args[0][2]
        assert "3D" in msg
        assert not dlg.running

    def test_stop_pushes_undo_checkpoint(self, dlg):
        dlg.target_mol = MagicMock()
        dlg.steps_since_checkpoint = 7
        dlg._stop()
        dlg.context.push_undo_checkpoint.assert_called_once()
        assert dlg.steps_since_checkpoint == 0
        assert dlg.lbl_state.text() == "State: Stopped"
        assert dlg.btn_start.isEnabled()
        assert not dlg.btn_stop.isEnabled()

    def test_stop_without_run_skips_checkpoint(self, dlg):
        dlg._stop()
        dlg.context.push_undo_checkpoint.assert_not_called()

    def test_accept_unregisters_window(self, dlg):
        dlg.accept()
        dlg.context.register_window.assert_called_with("main_panel", None)

    def test_accept_checkpoints_pending_steps(self, dlg):
        dlg.target_mol = MagicMock()
        dlg.steps_since_checkpoint = 3
        dlg.accept()
        dlg.context.push_undo_checkpoint.assert_called_once()

    def test_reject_restores_original_coordinates(self, dlg):
        dlg.target_mol = MagicMock()
        dlg.original_coords = [MagicMock(), MagicMock()]
        dlg.reject()
        conf = dlg.target_mol.GetConformer.return_value
        assert conf.SetAtomPosition.call_count == 2
        dlg.context.refresh_3d_view.assert_called_once()
        dlg.context.register_window.assert_called_with("main_panel", None)

    def test_run_plugin_without_molecule_warns(self, no_msgbox):
        ctx = _ctx()
        _step.run_plugin(ctx)
        no_msgbox[_step.__name__].warning.assert_called_once()
        ctx.register_window.assert_not_called()

    def test_run_plugin_raises_existing_window(self, no_msgbox):
        ctx = _ctx()
        ctx.current_mol = MagicMock()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _step.run_plugin(ctx)
        existing.show.assert_called_once()
        existing.raise_.assert_called_once()

    def test_initialize_registers_menu_action(self):
        ctx = _ctx()
        _step.initialize(ctx)
        assert ctx.add_menu_action.call_args[0][0] == "3D Edit/Step Optimizer..."


# ===========================================================================
# All-Trans Optimizer  (no dialog — torsion selection + guards)
# ===========================================================================


class TestAllTransSelectTorsions:
    def test_keeps_one_match_per_central_bond(self):
        matches = [(0, 1, 2, 3), (4, 1, 2, 5), (6, 1, 2, 7)]
        assert _all_trans._select_torsions(matches) == [(0, 1, 2, 3)]

    def test_reversed_central_bond_is_same_bond(self):
        matches = [(0, 1, 2, 3), (5, 2, 1, 6)]
        assert _all_trans._select_torsions(matches) == [(0, 1, 2, 3)]

    def test_distinct_bonds_all_kept_in_order(self):
        matches = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
        assert _all_trans._select_torsions(matches) == matches

    def test_empty_matches(self):
        assert _all_trans._select_torsions([]) == []


class TestAllTransRunPlugin:
    def test_no_molecule_warns(self, no_msgbox):
        ctx = _ctx()
        _all_trans.run_plugin(ctx)
        no_msgbox[_all_trans.__name__].warning.assert_called_once()

    def test_no_conformer_warns(self, no_msgbox):
        ctx = _ctx()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 0
        ctx.current_mol = mol
        _all_trans.run_plugin(ctx)
        msg = no_msgbox[_all_trans.__name__].warning.call_args[0][2]
        assert "3D" in msg

    def test_no_torsions_informs(self, no_msgbox):
        ctx = _ctx()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        mol.GetSubstructMatches.return_value = ()
        ctx.current_mol = mol
        _all_trans.run_plugin(ctx)
        no_msgbox[_all_trans.__name__].information.assert_called_once()
        ctx.push_undo_checkpoint.assert_not_called()

    def test_applies_one_dihedral_per_bond(self, no_msgbox, monkeypatch):
        set_dihedral = MagicMock()
        monkeypatch.setattr(
            _all_trans.rdMolTransforms, "SetDihedralDeg", set_dihedral
        )
        ctx = _ctx()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        mol.GetSubstructMatches.return_value = ((0, 1, 2, 3), (4, 1, 2, 5))
        ctx.current_mol = mol
        _all_trans.run_plugin(ctx)
        assert set_dihedral.call_count == 1
        assert set_dihedral.call_args[0][-1] == 180.0
        ctx.refresh_3d_view.assert_called_once()
        ctx.push_undo_checkpoint.assert_called_once()
        assert "1 torsions" in ctx.show_status_message.call_args[0][0]

    def test_initialize_registers_menu_action(self):
        ctx = _ctx()
        _all_trans.initialize(ctx)
        assert ctx.add_menu_action.call_args[0][0] == "3D Edit/All-Trans Optimizer"

    def test_plugin_name(self):
        assert _all_trans.PLUGIN_NAME == "All-Trans Optimizer"


# ===========================================================================
# XtbOptimizerDialog  (visible plugin: "xTB Optimizer")
# ===========================================================================


class TestXtbOptimizerDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _xtb.XtbOptimizerDialog(context=_ctx(), parent=None)
        yield d
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "xTB Geometry Optimizer"

    def test_method_options(self, dlg):
        texts = [dlg.combo_method.itemText(i) for i in range(dlg.combo_method.count())]
        assert texts == ["GFN2-xTB", "GFN1-xTB"]

    def test_convergence_defaults(self, dlg):
        assert dlg.spin_fmax.value() == pytest.approx(0.05)
        assert dlg.spin_maxsteps.value() == 500

    def test_table_columns(self, dlg):
        assert dlg.table.columnCount() == 3
        headers = [dlg.table.horizontalHeaderItem(i).text() for i in range(3)]
        assert headers == ["Step", "Energy (eV)", "Fmax (eV/Å)"]

    def test_log_is_readonly(self, dlg):
        assert dlg.log_text.isReadOnly()

    def test_initial_state(self, dlg):
        assert dlg.lbl_status.text() == "Ready."
        assert dlg.btn_run.isEnabled()
        assert not dlg.btn_cancel.isEnabled()
        assert dlg.progress_bar.isHidden()

    def test_set_running_toggles_controls(self, dlg):
        dlg._set_running(True)
        assert not dlg.btn_run.isEnabled()
        assert dlg.btn_cancel.isEnabled()
        assert not dlg.combo_method.isEnabled()
        assert not dlg.progress_bar.isHidden()
        dlg._set_running(False)
        assert dlg.btn_run.isEnabled()
        assert not dlg.btn_cancel.isEnabled()

    def test_append_log(self, dlg):
        dlg._append_log("line one")
        dlg._append_log("line two")
        assert "line one" in dlg.log_text.toPlainText()
        assert "line two" in dlg.log_text.toPlainText()

    def test_step_update_fills_table(self, dlg):
        dlg._on_step_update(1, -123.456789, 0.0567)
        assert dlg.table.rowCount() == 1
        assert dlg.table.item(0, 0).text() == "1"
        assert dlg.table.item(0, 1).text() == "-123.456789"
        assert dlg.table.item(0, 2).text() == "0.0567"
        assert "Step 1" in dlg.lbl_status.text()

    def test_run_without_molecule_warns(self, dlg, no_msgbox):
        dlg._on_run()
        no_msgbox[_xtb.__name__].warning.assert_called_once()
        assert dlg.btn_run.isEnabled()

    def test_run_without_conformer_warns(self, dlg, no_msgbox):
        mol = MagicMock()
        mol.GetNumConformers.return_value = 0
        dlg.context.current_mol = mol
        dlg._on_run()
        msg = no_msgbox[_xtb.__name__].warning.call_args[0][2]
        assert "3D" in msg

    def test_finished_failure_reports_status(self, dlg):
        dlg._set_running(True)
        dlg._on_finished(False, "Cancelled by user.")
        assert dlg.lbl_status.text() == "Stopped: Cancelled by user."
        assert dlg.btn_run.isEnabled()

    def test_finished_success_applies_coordinates(self, dlg):
        mol = MagicMock()
        dlg.context.current_mol = mol
        dlg._step_count = 12
        dlg._on_finished(True, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        conf = mol.GetConformer.return_value
        assert conf.SetAtomPosition.call_count == 2
        dlg.context.refresh_3d_view.assert_called_once()
        dlg.context.push_undo_checkpoint.assert_called_once()
        assert "complete" in dlg.lbl_status.text()

    def test_finished_success_with_unloaded_molecule(self, dlg):
        dlg.context.current_mol = None
        dlg._on_finished(True, [[0.0, 0.0, 0.0]])
        assert "unloaded" in dlg.lbl_status.text()

    def test_accept_unregisters_window(self, dlg):
        dlg.accept()
        dlg.context.register_window.assert_called_with("main_panel", None)

    def test_worker_cancel_sets_flag(self, qapp):
        w = _xtb.XtbWorker(
            numbers=[1, 1], positions=[[0, 0, 0], [0.74, 0, 0]],
            method="GFN2-xTB", fmax=0.05, max_steps=10,
        )
        assert not w._cancelled
        w.cancel()
        assert w._cancelled

    def test_run_plugin_without_molecule_warns(self, no_msgbox):
        ctx = _ctx()
        _xtb.run_plugin(ctx)
        no_msgbox[_xtb.__name__].warning.assert_called_once()

    def test_initialize_registers_menu_action(self):
        ctx = _ctx()
        _xtb.initialize(ctx)
        assert ctx.add_menu_action.call_args[0][0] == "3D Edit/xTB Optimizer…"


# ===========================================================================
# UntanglerDialog  (visible plugin: "Complex Molecule Untangler")
# ===========================================================================


class TestUntanglerDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _untangler.UntanglerDialog(context=_ctx(), parent=None)
        yield d
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Complex Molecule Untangler"

    def test_force_field_options(self, dlg):
        texts = [dlg.combo_ff.itemText(i) for i in range(dlg.combo_ff.count())]
        assert texts == ["MMFF94", "UFF"]

    def test_default_ff_from_uff_setting(self, qapp):
        d = _untangler.UntanglerDialog(context=_ctx(setting="UFF"), parent=None)
        assert d.combo_ff.currentText() == "UFF"
        d.destroy()

    def test_iterations_defaults(self, dlg):
        assert dlg.spin_iter.value() == 500
        assert (dlg.spin_iter.minimum(), dlg.spin_iter.maximum()) == (100, 10000)
        assert dlg.spin_iter.singleStep() == 100

    def test_progress_bar_initial(self, dlg):
        assert dlg.pbar.value() == 0
        assert not dlg.pbar.isTextVisible()

    def test_run_button(self, dlg):
        assert dlg.btn_run.text() == "Untangle Molecule"
        assert dlg.btn_run.isEnabled()

    def test_run_without_molecule_warns(self, dlg, no_msgbox):
        dlg.run_untangle()
        no_msgbox[_untangler.__name__].warning.assert_called_once()
        assert dlg.btn_run.isEnabled()
        assert dlg.worker is None

    def test_on_finished_error_restores_button(self, dlg, no_msgbox):
        dlg.btn_run.setEnabled(False)
        dlg.btn_run.setText("Processing...")
        dlg.on_finished(None, "No rotatable bonds found.")
        assert dlg.btn_run.isEnabled()
        assert dlg.btn_run.text() == "Untangle Molecule"
        no_msgbox[_untangler.__name__].warning.assert_called_once()

    def test_on_finished_success_updates_molecule(self, dlg, no_msgbox):
        new_mol = MagicMock()
        dlg.on_finished(new_mol, "Processed 5 bonds.")
        assert dlg.context.current_mol is new_mol
        dlg.context.refresh_3d_view.assert_called_once()
        dlg.context.push_undo_checkpoint.assert_called_once()
        no_msgbox[_untangler.__name__].information.assert_called_once()

    def test_worker_stores_parameters(self, qapp):
        mol = MagicMock()
        w = _untangler.UntangleWorker(mol, max_iter=250, force_field="UFF")
        assert w.mol is mol
        assert w.max_iter == 250
        assert w.force_field == "UFF"

    def test_run_plugin_registers_new_dialog(self):
        ctx = _ctx()
        ctx.get_window.return_value = None
        _untangler.run_plugin(ctx)
        args = ctx.register_window.call_args[0]
        assert args[0] == "main_panel"
        assert isinstance(args[1], _untangler.UntanglerDialog)
        args[1].destroy()

    def test_run_plugin_raises_existing_window(self):
        ctx = _ctx()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _untangler.run_plugin(ctx)
        existing.show.assert_called_once()
        ctx.register_window.assert_not_called()

    def test_initialize_then_run_launches(self):
        ctx = _ctx()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _untangler.initialize(ctx)
        _untangler.run(MagicMock())
        existing.show.assert_called_once()


# ===========================================================================
# Advanced Rendering — HideOnCloseDialog + AdvancedGraphicsWidget
# ===========================================================================


class TestHideOnCloseDialog:
    def test_close_hides_instead_of_destroying(self, qapp):
        dlg = _adv.HideOnCloseDialog()
        dlg.show()
        dlg.close()
        assert dlg.isHidden()

    def test_close_saves_viewer_settings(self, qapp):
        dlg = _adv.HideOnCloseDialog()
        viewer = MagicMock()
        dlg.set_viewer(viewer)
        dlg.show()
        dlg.close()
        viewer.save_settings.assert_called_once()
        dlg.destroy()

    def test_get_icon_returns_none(self):
        assert _adv.get_icon() is None


def _make_adv_widget(tmp_path):
    """AdvancedGraphicsWidget on a stub main window, settings redirected to tmp."""
    from PyQt6.QtWidgets import QWidget

    mw = QWidget()
    mw.plotter = None
    mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
    w = _adv.AdvancedGraphicsWidget(mw)
    w.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
    return mw, w


def _teardown_adv_widget(mw, w):
    from PyQt6.QtCore import QCoreApplication

    w._sync_timer.stop()
    try:
        QCoreApplication.instance().aboutToQuit.disconnect(w.save_settings)
    except TypeError:
        pass
    w.destroy()
    mw.destroy()


class TestAdvancedGraphicsWidget:
    @pytest.fixture
    def widget(self, qapp, tmp_path):
        mw, w = _make_adv_widget(tmp_path)
        yield w
        _teardown_adv_widget(mw, w)

    def test_default_state(self, widget):
        assert widget.light_intensity == pytest.approx(1.0)
        assert widget.light_elevation == 45
        assert not widget.use_shadows
        assert not widget.use_ssao
        assert not widget.use_atom_pbr
        assert widget.atom_roughness == pytest.approx(0.5)

    def test_default_preset_exists(self, widget):
        assert "Default" in widget.presets
        assert "Default" in widget.default_preset_names

    def test_preset_combo_lists_default(self, widget):
        texts = [
            widget.combo_presets.itemText(i)
            for i in range(widget.combo_presets.count())
        ]
        assert "Default" in texts

    def test_slider_defaults(self, widget):
        assert widget.slider_light.value() == 100
        assert widget.slider_ele.value() == 45
        assert widget.slider_azi.value() == 0
        assert not widget.slider_edl.isEnabled()
        assert not widget.slider_atom_metallic.isEnabled()
        assert not widget.slider_atom_roughness.isEnabled()

    def test_gather_settings_matches_default_preset_keys(self, widget):
        assert set(widget.gather_settings_dict()) == set(widget.presets["Default"])

    def test_gather_settings_reflects_state(self, widget):
        widget.light_intensity = 2.5
        widget.use_aa = True
        data = widget.gather_settings_dict()
        assert data["light_intensity"] == pytest.approx(2.5)
        assert data["use_aa"] is True

    def test_sync_style_ui_advanced_enables_pbr(self, widget):
        widget.sync_style_ui("Ball & Stick (Advanced Rendering)")
        assert widget.use_atom_pbr
        assert widget.check_atom_pbr.isChecked()
        assert widget.slider_atom_metallic.isEnabled()
        assert widget.slider_atom_roughness.isEnabled()

    def test_sync_style_ui_standard_disables_pbr(self, widget):
        widget.sync_style_ui("Ball & Stick (Advanced Rendering)")
        widget.sync_style_ui("ball_and_stick")
        assert not widget.use_atom_pbr
        assert not widget.check_atom_pbr.isChecked()
        assert not widget.slider_atom_metallic.isEnabled()

    def test_ssao_refused_on_standard_style(self, widget):
        widget.check_ssao.setChecked(True)
        assert not widget.check_ssao.isChecked()
        assert not widget.use_ssao

    def test_edl_refused_on_standard_style(self, widget):
        widget.check_edl.setChecked(True)
        assert not widget.check_edl.isChecked()
        assert not widget.use_edl

    def test_depth_peeling_refused_on_standard_style(self, widget):
        widget.check_depth.setChecked(True)
        assert not widget.check_depth.isChecked()
        assert not widget.use_depth_peeling

    def test_metallic_slider_updates_state(self, widget):
        widget.slider_atom_metallic.setValue(40)
        assert widget.atom_metallic == pytest.approx(0.4)

    def test_light_slider_updates_intensity(self, widget):
        widget.slider_light.setValue(250)
        assert widget.light_intensity == pytest.approx(2.5)

    def test_save_settings_writes_json(self, widget, tmp_path):
        import json

        widget.use_aa = True
        widget.save_settings()
        data = json.loads((tmp_path / "adv_settings.json").read_text())
        assert data["use_aa"] is True
        assert "presets" in data

    def test_delete_default_preset_is_locked(self, widget, monkeypatch):
        box = MagicMock()
        monkeypatch.setattr(_adv, "QMessageBox", box)
        widget.combo_presets.setCurrentText("Default")
        widget.delete_preset()
        assert "Default" in widget.presets
        box.warning.assert_called_once()

    def test_update_preset_combo_shows_new_presets(self, widget):
        widget.presets["MyLook"] = widget.gather_settings_dict()
        widget.update_preset_combo()
        texts = [
            widget.combo_presets.itemText(i)
            for i in range(widget.combo_presets.count())
        ]
        assert "MyLook" in texts


class TestAdvancedRenderingInitialize:
    def test_registers_styles_and_menu(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw

        _adv.initialize(ctx)
        viewer = mw._adv_rendering_viewer
        viewer.get_settings_path = lambda: str(tmp_path / "adv_settings.json")

        assert ctx.register_3d_style.call_count == 4
        style_names = [c[0][0] for c in ctx.register_3d_style.call_args_list]
        assert "Ball & Stick (Advanced Rendering)" in style_names
        assert ctx.add_menu_action.call_args[0][0] == "Settings/Advanced Graphics Settings"
        assert isinstance(viewer, _adv.AdvancedGraphicsWidget)

        _teardown_adv_widget(mw, viewer)

    def test_style_drawer_survives_bare_main_window(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw

        _adv.initialize(ctx)
        viewer = mw._adv_rendering_viewer
        viewer.get_settings_path = lambda: str(tmp_path / "adv_settings.json")

        drawer = ctx.register_3d_style.call_args_list[0][0][1]
        drawer(SimpleNamespace(), None)  # no view_3d_manager, no viewer → no crash

        _teardown_adv_widget(mw, viewer)
