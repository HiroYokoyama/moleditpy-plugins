"""
Headless GUI tests for the Cube File Viewer Advanced plugin.

Covers: ChargeDialog, FlexibleDoubleSpinBox.

Chemistry libs (pyvista, numpy, rdkit, …) are mocked; real PyQt6 is used.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CUBE_ADV_PATH = PLUGINS_DIR / "Cube_File_Viewer_Advanced" / "cube_viewer_advanced.py"

with mock_chemistry_imports():
    _cube_adv = load_plugin_for_gui(CUBE_ADV_PATH)


# ===========================================================================
# ChargeDialog  (Cube File Viewer Advanced)
# ===========================================================================


class TestCubeAdvancedChargeDialog:
    """ChargeDialog in Cube File Viewer Advanced — same contract as the basic version."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _cube_adv.ChargeDialog(parent=None, current_charge=-1)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Bond Connectivity Error"

    def test_initial_charge_value(self, dlg):
        assert dlg.spin.value() == -1

    def test_spin_range(self, dlg):
        assert dlg.spin.minimum() <= -10
        assert dlg.spin.maximum() >= 10

    def test_default_result_action(self, dlg):
        assert dlg.result_action == "cancel"

    def test_on_skip_sets_action(self, dlg):
        dlg.on_skip()
        assert dlg.result_action == "skip"


# ===========================================================================
# FlexibleDoubleSpinBox  (Cube File Viewer Advanced)
# ===========================================================================


class TestFlexibleDoubleSpinBox:
    """FlexibleDoubleSpinBox: QDoubleSpinBox subclass with custom text formatting."""

    @pytest.fixture
    def spin(self, qapp):
        s = _cube_adv.FlexibleDoubleSpinBox()
        yield s
        s.destroy()

    def test_creates_without_error(self, spin):
        assert spin is not None

    def test_high_decimal_places(self, spin):
        assert spin.decimals() >= 5

    def test_text_from_value_strips_trailing_zeros(self, spin):
        assert spin.textFromValue(1.5) == "1.5"

    def test_text_from_value_integer(self, spin):
        assert spin.textFromValue(2.0) == "2"

    def test_text_from_value_zero(self, spin):
        assert spin.textFromValue(0.0) == "0"

    def test_text_from_value_negative(self, spin):
        assert spin.textFromValue(-1.25) == "-1.25"

    def test_decimals_setting(self, spin):
        assert spin.decimals() == 10


# ===========================================================================
# ChargeDialog — actions
# ===========================================================================


class TestCubeAdvancedChargeDialogActions:
    @pytest.fixture
    def dlg(self, qapp):
        d = _cube_adv.ChargeDialog(parent=None, current_charge=2)
        yield d
        d.destroy()

    def test_on_retry_sets_action_and_charge_from_spin(self, dlg):
        dlg.spin.setValue(5)
        dlg.on_retry()
        assert dlg.result_action == "retry"
        assert dlg.charge == 5

    def test_buttons_exist(self, dlg):
        texts = []
        for btns in dlg.findChildren(_cube_adv.QPushButton):
            texts.append(btns.text())
        assert set(texts) == {"Retry", "Skip Chemistry", "Cancel"}


# ===========================================================================
# CubeViewerWidget — construction, no molecule
# ===========================================================================

from unittest.mock import MagicMock, patch
from types import SimpleNamespace

_cube_adv.np = SimpleNamespace(generic=(), ndarray=())


def _make_main_window():
    # Real QWidget so it can be a Qt parent; mocked methods for assertions.
    from PyQt6.QtWidgets import QWidget

    mw = QWidget()
    mw.removeDockWidget = MagicMock()
    mw.current_mol = None
    mw.plotter = MagicMock()
    return mw


def _make_viewer(data_max=2.0, mol=None):
    ctx = MagicMock()
    mw = _make_main_window()
    ctx.get_main_window.return_value = mw
    ctx.plotter = MagicMock()
    ctx.current_molecule = mol
    dock = MagicMock()
    grid = MagicMock()
    w = _cube_adv.CubeViewerWidget(ctx, dock, grid, data_max=data_max)
    return w, ctx, mw, dock, grid


def _teardown_viewer(w):
    from PyQt6.QtCore import QCoreApplication

    w._structure_watch_timer.stop()
    try:
        QCoreApplication.instance().aboutToQuit.disconnect(w.save_settings)
    except TypeError:
        pass
    w.destroy()


class TestCubeViewerWidgetConstruction:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w
        _teardown_viewer(w)

    def test_creates_without_error(self, viewer):
        assert viewer is not None

    def test_data_max_floor(self, qapp):
        w, *_ = _make_viewer(data_max=0.0)
        assert w.data_max == pytest.approx(1e-6)
        _teardown_viewer(w)

    def test_isovalue_spin_range(self, viewer):
        assert viewer.spin.minimum() == pytest.approx(0.00001)
        assert viewer.spin.maximum() == pytest.approx(2.0 * 1.2)

    def test_isovalue_default_value(self, viewer):
        assert viewer.spin.value() == pytest.approx(0.05)

    def test_slider_range(self, viewer):
        assert viewer.slider.minimum() == 0
        assert viewer.slider.maximum() == 1000

    def test_opacity_spin_default(self, viewer):
        assert viewer.opacity_spin.value() == pytest.approx(0.4)

    def test_opacity_slider_default(self, viewer):
        assert viewer.opacity_slider.value() == 40

    def test_style_combo_options(self, viewer):
        texts = [viewer.combo_style.itemText(i) for i in range(viewer.combo_style.count())]
        assert texts == ["Surface", "Smoothed Surface", "Wireframe", "Points"]

    def test_smooth_checkbox_checked_by_default(self, viewer):
        assert viewer.check_smooth.isChecked()

    def test_pbr_unchecked_and_sliders_disabled(self, viewer):
        assert not viewer.check_pbr.isChecked()
        assert not viewer.slider_metallic.isEnabled()
        assert not viewer.slider_roughness.isEnabled()

    def test_edl_slider_disabled_initially(self, viewer):
        assert not viewer.slider_edl.isEnabled()

    def test_preset_combo_has_default(self, viewer):
        texts = [viewer.combo_presets.itemText(i) for i in range(viewer.combo_presets.count())]
        assert "Default" in texts

    def test_close_button_text(self, viewer):
        found = [b for b in viewer.findChildren(_cube_adv.QPushButton) if b.text() == "Close Plugin"]
        assert len(found) == 1

    def test_reset_button_text(self, viewer):
        found = [b for b in viewer.findChildren(_cube_adv.QPushButton) if b.text() == "Reset All Settings"]
        assert len(found) == 1

    def test_env_path_placeholder(self, viewer):
        assert "texture" in viewer.line_env_path.placeholderText().lower()

    def test_structure_watch_timer_running(self, viewer):
        assert viewer._structure_watch_timer.isActive()

    def test_bound_mol_initially_none(self, viewer):
        assert viewer._bound_mol is None

    def test_default_colors(self, viewer):
        assert viewer.color_p == (0, 0, 255)
        assert viewer.color_n == (255, 0, 0)


# ===========================================================================
# CubeViewerWidget — isovalue / opacity signal wiring
# ===========================================================================


class TestCubeViewerWidgetSignals:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w
        _teardown_viewer(w)

    def test_slider_change_updates_spin(self, viewer):
        viewer.slider.setValue(500)
        assert viewer.spin.value() == pytest.approx((500 / 1000) * viewer.max_val, rel=1e-3)

    def test_spin_change_updates_slider(self, viewer):
        viewer.spin.setValue(viewer.max_val)
        assert viewer.slider.value() == 1000

    def test_opacity_slider_change_updates_spin(self, viewer):
        viewer.opacity_slider.setValue(75)
        assert viewer.opacity_spin.value() == pytest.approx(0.75)

    def test_opacity_spin_change_updates_slider(self, viewer):
        viewer.opacity_spin.setValue(0.9)
        assert viewer.opacity_slider.value() == 90

    def test_pbr_toggle_enables_sliders(self, viewer):
        viewer.check_pbr.setChecked(True)
        assert viewer.use_pbr
        assert viewer.slider_metallic.isEnabled()
        assert viewer.slider_roughness.isEnabled()

    def test_pbr_untoggle_disables_sliders(self, viewer):
        viewer.check_pbr.setChecked(True)
        viewer.check_pbr.setChecked(False)
        assert not viewer.use_pbr
        assert not viewer.slider_metallic.isEnabled()

    def test_metallic_slider_updates_state(self, viewer):
        viewer.slider_metallic.setValue(70)
        assert viewer.metallic == pytest.approx(0.7)

    def test_roughness_slider_updates_state(self, viewer):
        viewer.slider_roughness.setValue(30)
        assert viewer.roughness == pytest.approx(0.3)

    def test_style_wireframe_switches_label_to_density(self, viewer):
        viewer.combo_style.setCurrentText("Wireframe")
        assert viewer.opacity_label.text() == "Density:"

    def test_style_surface_keeps_opacity_label(self, viewer):
        viewer.combo_style.setCurrentText("Wireframe")
        viewer.combo_style.setCurrentText("Surface")
        assert viewer.opacity_label.text() == "Opacity:"

    def test_comp_color_toggle_disables_neg_button(self, viewer):
        viewer.check_comp_color.setChecked(True)
        assert not viewer.btn_color_n.isEnabled()

    def test_comp_color_untoggle_enables_neg_button(self, viewer):
        viewer.check_comp_color.setChecked(True)
        viewer.check_comp_color.setChecked(False)
        assert viewer.btn_color_n.isEnabled()

    def test_silhouette_checkbox_updates_state(self, viewer):
        viewer.check_silhouette.setChecked(True)
        assert viewer.use_silhouette

    def test_aa_checkbox_updates_state(self, viewer):
        viewer.check_aa.setChecked(True)
        assert viewer.use_aa

    def test_edl_checkbox_enables_slider(self, viewer):
        viewer.check_edl.setChecked(True)
        assert viewer.slider_edl.isEnabled()
        assert viewer.use_edl

    def test_edl_strength_slider_updates_state(self, viewer):
        viewer.check_edl.setChecked(True)
        viewer.slider_edl.setValue(60)
        assert viewer.edl_strength == pytest.approx(0.6)

    def test_light_intensity_slider_updates_state(self, viewer):
        viewer.slider_light.setValue(200)
        assert viewer.light_intensity == pytest.approx(2.0)

    def test_shadows_checkbox_updates_state(self, viewer):
        viewer.check_shadows.setChecked(True)
        assert viewer.use_shadows


# ===========================================================================
# CubeViewerWidget — mutually-exclusive rendering effects
# ===========================================================================


class TestCubeViewerMutualExclusivity:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w
        _teardown_viewer(w)

    def test_depth_peeling_disables_edl_and_shadows(self, viewer):
        viewer.check_edl.setChecked(True)
        viewer.check_shadows.setChecked(True)
        viewer.check_depth.setChecked(True)
        assert not viewer.check_edl.isChecked()
        assert not viewer.check_shadows.isChecked()
        assert not viewer.use_edl
        assert not viewer.use_shadows

    def test_edl_disables_depth_peeling(self, viewer):
        viewer.check_depth.setChecked(True)
        viewer.check_edl.setChecked(True)
        assert not viewer.check_depth.isChecked()
        assert not viewer.use_depth_peeling


# ===========================================================================
# CubeViewerWidget — presets (QInputDialog / QMessageBox mocked to avoid
# blocking modal exec() under offscreen)
# ===========================================================================


class TestCubeViewerPresets:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        w.get_settings_path = lambda: str(_tmp_settings_path)
        yield w
        _teardown_viewer(w)

    @pytest.fixture(autouse=True)
    def _tmp_path_setup(self, tmp_path):
        global _tmp_settings_path
        _tmp_settings_path = tmp_path / "cube_viewer_advanced.json"

    def test_save_preset_adds_to_presets_and_combo(self, viewer, monkeypatch):
        monkeypatch.setattr(_cube_adv.QInputDialog, "getText", staticmethod(lambda *a, **kw: ("MyPreset", True)))
        viewer.save_preset()
        assert "MyPreset" in viewer.presets
        texts = [viewer.combo_presets.itemText(i) for i in range(viewer.combo_presets.count())]
        assert "MyPreset" in texts

    def test_save_preset_cancelled_does_nothing(self, viewer, monkeypatch):
        monkeypatch.setattr(_cube_adv.QInputDialog, "getText", staticmethod(lambda *a, **kw: ("Ignored", False)))
        viewer.save_preset()
        assert "Ignored" not in viewer.presets

    def test_save_preset_default_name_blocked(self, viewer, monkeypatch):
        monkeypatch.setattr(_cube_adv.QInputDialog, "getText", staticmethod(lambda *a, **kw: ("Default", True)))
        box = MagicMock()
        monkeypatch.setattr(_cube_adv, "QMessageBox", box)
        viewer.spin.setValue(0.9)
        viewer.save_preset()
        box.warning.assert_called_once()
        assert viewer.presets["Default"]["isovalue"] != pytest.approx(0.9)

    def test_delete_preset_removes_custom(self, viewer, monkeypatch):
        monkeypatch.setattr(_cube_adv.QInputDialog, "getText", staticmethod(lambda *a, **kw: ("Temp", True)))
        viewer.save_preset()
        viewer.combo_presets.setCurrentText("Temp")
        viewer.delete_preset()
        assert "Temp" not in viewer.presets

    def test_delete_default_preset_blocked(self, viewer, monkeypatch):
        box = MagicMock()
        monkeypatch.setattr(_cube_adv, "QMessageBox", box)
        viewer.combo_presets.setCurrentText("Default")
        viewer.delete_preset()
        assert "Default" in viewer.presets
        box.warning.assert_called_once()

    def test_load_preset_settings_applies_isovalue_and_style(self, viewer):
        viewer.load_preset_settings({"isovalue": 0.33, "style": "Points"})
        assert viewer.spin.value() == pytest.approx(0.33)
        assert viewer.combo_style.currentText() == "Points"

    def test_on_preset_activated_loads_settings(self, viewer, monkeypatch):
        monkeypatch.setattr(_cube_adv.QInputDialog, "getText", staticmethod(lambda *a, **kw: ("Temp2", True)))
        viewer.spin.setValue(0.77)
        viewer.save_preset()
        viewer.spin.setValue(0.01)
        idx = viewer.combo_presets.findText("Temp2")
        viewer.combo_presets.setCurrentIndex(idx)
        viewer.on_preset_activated(idx)
        assert viewer.spin.value() == pytest.approx(0.77)


# ===========================================================================
# CubeViewerWidget — settings persistence round trip
# ===========================================================================


class TestCubeViewerSettingsPersistence:
    @pytest.fixture
    def viewer(self, qapp, tmp_path):
        w, ctx, mw, dock, grid = _make_viewer()
        w.get_settings_path = lambda: str(tmp_path / "cube_viewer_advanced.json")
        yield w
        _teardown_viewer(w)

    def test_save_settings_writes_json(self, viewer, tmp_path):
        import json

        viewer.spin.setValue(0.42)
        viewer.save_settings()
        data = json.loads((tmp_path / "cube_viewer_advanced.json").read_text())
        assert data["isovalue"] == pytest.approx(0.42)
        assert "presets" in data

    def test_load_settings_restores_isovalue(self, qapp, tmp_path):
        w, ctx, mw, dock, grid = _make_viewer()
        w.get_settings_path = lambda: str(tmp_path / "cube_viewer_advanced.json")
        w.spin.setValue(0.66)
        w.save_settings()
        w.spin.setValue(0.01)
        w.load_settings()
        assert w.spin.value() == pytest.approx(0.66)
        _teardown_viewer(w)


# ===========================================================================
# CubeViewerWidget — color pickers (QColorDialog mocked; static exec avoided)
# ===========================================================================


class TestCubeViewerColorPickers:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w
        _teardown_viewer(w)

    def test_choose_color_p_valid_updates_color(self, viewer, monkeypatch):
        from PyQt6.QtGui import QColor

        monkeypatch.setattr(_cube_adv.QColorDialog, "getColor", staticmethod(lambda **kw: QColor(10, 20, 30)))
        viewer.choose_color_p()
        assert viewer.color_p == (10, 20, 30)

    def test_choose_color_p_invalid_leaves_color(self, viewer, monkeypatch):
        from PyQt6.QtGui import QColor

        monkeypatch.setattr(_cube_adv.QColorDialog, "getColor", staticmethod(lambda **kw: QColor()))
        before = viewer.color_p
        viewer.choose_color_p()
        assert viewer.color_p == before

    def test_choose_color_n_valid_updates_color(self, viewer, monkeypatch):
        from PyQt6.QtGui import QColor

        monkeypatch.setattr(_cube_adv.QColorDialog, "getColor", staticmethod(lambda **kw: QColor(1, 2, 3)))
        viewer.choose_color_n()
        assert viewer.color_n == (1, 2, 3)

    def test_choose_color_p_with_complementary_updates_neg(self, viewer, monkeypatch):
        from PyQt6.QtGui import QColor

        viewer.check_comp_color.setChecked(True)
        monkeypatch.setattr(_cube_adv.QColorDialog, "getColor", staticmethod(lambda **kw: QColor(0, 255, 0)))
        viewer.choose_color_p()
        assert viewer.color_n != (255, 0, 0)  # complementary recomputed


# ===========================================================================
# CubeViewerWidget — reset-all
# ===========================================================================


class TestCubeViewerResetAll:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w
        _teardown_viewer(w)

    def test_reset_restores_isovalue_and_opacity(self, viewer):
        viewer.spin.setValue(0.9)
        viewer.opacity_spin.setValue(0.1)
        viewer.on_reset_all()
        assert viewer.spin.value() == pytest.approx(0.02)
        assert viewer.opacity_spin.value() == pytest.approx(0.4)

    def test_reset_clears_advanced_effects(self, viewer):
        viewer.check_pbr.setChecked(True)
        viewer.check_ssao.setChecked(True)
        viewer.check_aa.setChecked(True)
        viewer.on_reset_all()
        assert not viewer.check_pbr.isChecked()
        assert not viewer.check_ssao.isChecked()
        assert not viewer.check_aa.isChecked()
        assert not viewer.use_ssao
        assert not viewer.use_aa

    def test_reset_restores_style_and_smoothing(self, viewer):
        viewer.combo_style.setCurrentText("Points")
        viewer.check_smooth.setChecked(False)
        viewer.on_reset_all()
        assert viewer.combo_style.currentText() == "Surface"
        assert viewer.check_smooth.isChecked()


# ===========================================================================
# CubeViewerWidget — structure-change watchdog
# ===========================================================================


class TestCubeViewerStructureWatch:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w, mw
        _teardown_viewer(w)

    def test_bind_structure_records_mol(self, viewer):
        w, mw = viewer
        mol = object()
        w.bind_structure(mol)
        assert w._bound_mol is mol

    def test_check_still_bound_noop_without_binding(self, viewer):
        w, mw = viewer
        w._check_structure_still_bound()
        assert w._structure_watch_timer.isActive()

    def test_check_detects_mismatch_and_detaches(self, viewer):
        w, mw = viewer
        mol = object()
        w.bind_structure(mol)
        mw.current_mol = object()  # a different molecule now shown
        w._check_structure_still_bound()
        assert not w._structure_watch_timer.isActive()
        mw.removeDockWidget.assert_called_once()

    def test_check_no_mismatch_keeps_timer_running(self, viewer):
        w, mw = viewer
        mol = object()
        w.bind_structure(mol)
        mw.current_mol = mol
        w._check_structure_still_bound()
        assert w._structure_watch_timer.isActive()


# ===========================================================================
# CubeViewerWidget — close / cleanup
# ===========================================================================


class TestCubeViewerClose:
    def test_close_plugin_removes_dock_and_clears_state(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        w._structure_watch_timer.stop()
        w.close_plugin()
        mw.removeDockWidget.assert_called_once_with(dock)
        assert w.dock is None
        assert ctx.current_molecule is None

    def test_close_event_stops_timer_and_saves(self, qapp, tmp_path):
        w, ctx, mw, dock, grid = _make_viewer()
        w.get_settings_path = lambda: str(tmp_path / "cube_viewer_advanced.json")
        from PyQt6.QtGui import QCloseEvent

        w.closeEvent(QCloseEvent())
        assert not w._structure_watch_timer.isActive()
        assert (tmp_path / "cube_viewer_advanced.json").exists()
        w.destroy()


# ===========================================================================
# CubeViewerWidget — environment texture path field
# ===========================================================================


class TestCubeViewerTexturePath:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w
        _teardown_viewer(w)

    def test_texture_path_entered_missing_file_does_not_crash(self, viewer, capsys):
        viewer.line_env_path.setText("Z:/does/not/exist.png")
        viewer.on_texture_path_entered()
        assert "not found" in capsys.readouterr().out

    def test_texture_path_entered_blank_is_noop(self, viewer):
        viewer.line_env_path.setText("")
        viewer.on_texture_path_entered()  # no crash

    def test_remove_env_texture_clears_field(self, viewer):
        viewer.line_env_path.setText("some/path.png")
        viewer.env_texture_path = "some/path.png"
        viewer.on_remove_env_texture()
        assert viewer.line_env_path.text() == ""
        assert viewer.env_texture_path == ""


# ===========================================================================
# initialize() / run_plugin() / open_cube_viewer entry points
# ===========================================================================


class TestCubeAdvancedInitialize:
    def test_registers_cube_and_cub_openers(self, qapp):
        ctx = MagicMock()
        _cube_adv.initialize(ctx)
        exts = [c.args[0] for c in ctx.register_file_opener.call_args_list]
        assert exts == [".cube", ".cub"]

    def test_registers_drop_handler_priority_10(self, qapp):
        ctx = MagicMock()
        _cube_adv.initialize(ctx)
        ctx.register_drop_handler.assert_called_once()
        assert ctx.register_drop_handler.call_args.kwargs.get("priority") == 10

    def test_file_opener_wrapper_calls_open_cube_viewer(self, qapp):
        ctx = MagicMock()
        _cube_adv.initialize(ctx)
        wrapper = ctx.register_file_opener.call_args_list[0].args[1]
        with patch.object(_cube_adv, "open_cube_viewer") as ocv:
            wrapper("foo.cube")
            ocv.assert_called_once_with(ctx, "foo.cube")

    def test_drop_handler_accepts_only_cube_extensions(self, qapp):
        ctx = MagicMock()
        _cube_adv.initialize(ctx)
        handler = ctx.register_drop_handler.call_args.args[0]
        with patch.object(_cube_adv, "open_cube_viewer") as ocv:
            assert handler("orbital.cube") is True
            assert handler("orbital.CUB") is True
            assert handler("structure.xyz") is False
            assert ocv.call_count == 2


class TestRunPlugin:
    def test_run_plugin_cancelled_dialog_does_not_open_viewer(self, qapp, monkeypatch):
        ctx = MagicMock()
        monkeypatch.setattr(_cube_adv.QFileDialog, "getOpenFileName", staticmethod(lambda *a, **kw: ("", "")))
        with patch.object(_cube_adv, "open_cube_viewer") as ocv:
            _cube_adv.run_plugin(ctx)
            ocv.assert_not_called()

    def test_run_plugin_opens_selected_file(self, qapp, monkeypatch):
        ctx = MagicMock()
        monkeypatch.setattr(
            _cube_adv.QFileDialog, "getOpenFileName", staticmethod(lambda *a, **kw: ("orbital.cube", ""))
        )
        with patch.object(_cube_adv, "open_cube_viewer") as ocv:
            _cube_adv.run_plugin(ctx)
            ocv.assert_called_once_with(ctx, "orbital.cube")

    def test_run_returns_early_without_plugin_manager(self, qapp):
        mw = MagicMock(spec=[])
        _cube_adv.run(mw)  # no plugin_manager attr → early return, no crash
