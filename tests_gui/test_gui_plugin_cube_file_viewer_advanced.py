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


@pytest.fixture(autouse=True)
def _no_qt_event_pump(monkeypatch):
    """Prevent the offscreen Qt fatal-abort seen intermittently on CI.

    CubeViewerWidget.__init__ schedules QTimer.singleShot(100, initial_update);
    load_env_texture() then calls QCoreApplication.processEvents(), which drains
    that leaked timer against a plotter a later test has already swapped/destroyed
    -> `QTimerInfoList::activateTimers -> qt_assert -> abort` (Fatal Python error:
    Aborted, core dumped). No test relies on the deferred callback or the pump, so
    neutralize both at the source module.
    """
    monkeypatch.setattr(_cube_adv.QTimer, "singleShot",
                        staticmethod(lambda *a, **k: None), raising=False)
    monkeypatch.setattr(_cube_adv.QCoreApplication, "processEvents",
                        staticmethod(lambda *a, **k: None), raising=False)


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
    # contour().n_points must be a real int (not a bare MagicMock) so the
    # `> 0` comparison in update_iso() doesn't raise and short-circuit into
    # the outer except before reaching add_mesh/process_mesh.
    grid.contour.return_value.n_points = 5
    w = _cube_adv.CubeViewerWidget(ctx, dock, grid, data_max=data_max)
    # Keep test-triggered saves (closeEvent, aboutToQuit) out of plugins/
    import os
    import tempfile

    w.get_settings_path = lambda: os.path.join(
        tempfile.gettempdir(), "cube_viewer_advanced_test.json"
    )
    return w, ctx, mw, dock, grid


def _teardown_viewer(w):
    from PyQt6.QtCore import QCoreApplication

    w._structure_watch_timer.stop()
    try:
        QCoreApplication.instance().aboutToQuit.disconnect(w.save_settings)
    except TypeError:
        pass
    w.destroy()


def _fallback_plotter_getattr(exclude):
    """__getattr__ builder for hand-rolled fake plotters used mid-test.

    CubeViewerWidget.__init__ schedules ``QTimer.singleShot(100,
    self.initial_update)``; that callback fires whenever the Qt event queue
    is next drained (e.g. by another test's ``QCoreApplication.processEvents()``
    call inside load_env_texture), which can land well after this test has
    already swapped ``self.plotter`` for a deliberately incomplete fake.
    Any attribute not in ``exclude`` resolves to a no-op so that late-firing
    callback never crashes on an unrelated missing method, while attributes
    in ``exclude`` still raise AttributeError so hasattr() checks under test
    stay accurate.
    """

    def _getattr(self, name):
        if name in exclude:
            raise AttributeError(name)
        return lambda *a, **kw: None

    return _getattr


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


# ===========================================================================
# CubeViewerWidget — plotter-capability gating at construction
# ===========================================================================


class TestCubeViewerPlotterCapabilityGating:
    def test_ssao_and_depth_checkboxes_disabled_without_plotter_methods(self, qapp):
        class BarePlotter:
            renderer = None
            __getattr__ = _fallback_plotter_getattr({"enable_ssao", "enable_depth_peeling"})

        ctx = MagicMock()
        mw = _make_main_window()
        ctx.get_main_window.return_value = mw
        ctx.plotter = BarePlotter()
        ctx.current_molecule = None
        dock = MagicMock()
        grid = MagicMock()
        grid.contour.return_value.n_points = 5
        w = _cube_adv.CubeViewerWidget(ctx, dock, grid, data_max=1.0)
        assert not w.check_ssao.isEnabled()
        assert not w.check_depth.isEnabled()
        _teardown_viewer(w)


# ===========================================================================
# CubeViewerWidget — load_settings extra branches (presets merge, style
# case-insensitive match, env-texture sync paths, malformed JSON)
# ===========================================================================


class TestCubeViewerLoadSettingsBranches:
    @pytest.fixture
    def viewer(self, qapp, tmp_path):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w, mw, tmp_path
        _teardown_viewer(w)

    def test_load_settings_merges_user_presets(self, viewer):
        import json

        w, mw, tmp_path = viewer
        settings_path = tmp_path / "s1.json"
        settings_path.write_text(json.dumps({"presets": {"Custom": {"isovalue": 0.5}}}))
        w.get_settings_path = lambda: str(settings_path)
        w.load_settings()
        assert "Custom" in w.presets
        texts = [w.combo_presets.itemText(i) for i in range(w.combo_presets.count())]
        assert "Custom" in texts

    def test_load_settings_style_case_insensitive_match(self, viewer):
        import json

        w, mw, tmp_path = viewer
        settings_path = tmp_path / "s2.json"
        settings_path.write_text(json.dumps({"style": "points"}))
        w.get_settings_path = lambda: str(settings_path)
        w.load_settings()
        assert w.combo_style.currentText() == "Points"

    def test_load_settings_style_no_match_leaves_style_unset(self, viewer):
        import json

        w, mw, tmp_path = viewer
        settings_path = tmp_path / "s3.json"
        settings_path.write_text(json.dumps({"style": "not-a-real-style"}))
        w.get_settings_path = lambda: str(settings_path)
        before = w.combo_style.currentText()
        w.load_settings()
        assert w.combo_style.currentText() == before

    def test_load_settings_malformed_json_does_not_crash(self, viewer):
        w, mw, tmp_path = viewer
        settings_path = tmp_path / "bad.json"
        settings_path.write_text("{not valid json")
        w.get_settings_path = lambda: str(settings_path)
        w.load_settings()  # caught by the outer except, no crash


class TestCubeViewerLoadSettingsTexture:
    @pytest.fixture
    def viewer(self, qapp, tmp_path):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w, mw, tmp_path
        _teardown_viewer(w)

    def test_load_settings_sets_text_field_for_any_path(self, viewer):
        import json

        w, mw, tmp_path = viewer
        settings_path = tmp_path / "s.json"
        settings_path.write_text(json.dumps({"env_texture_path": "Z:/nowhere.png"}))
        w.get_settings_path = lambda: str(settings_path)
        w.load_settings()
        assert w.line_env_path.text() == "Z:/nowhere.png"

    def test_load_settings_path_with_space_is_skipped(self, viewer):
        import json

        w, mw, tmp_path = viewer
        settings_path = tmp_path / "s.json"
        settings_path.write_text(json.dumps({"env_texture_path": "path with space.png"}))
        w.get_settings_path = lambda: str(settings_path)
        w.load_settings()  # no crash; texture not loaded (VTK can't handle spaces)
        assert w.line_env_path.text() == "path with space.png"

    def test_load_settings_existing_path_triggers_load(self, viewer):
        import json

        w, mw, tmp_path = viewer
        tex = tmp_path / "tex.png"
        tex.write_bytes(b"fake")
        settings_path = tmp_path / "s.json"
        settings_path.write_text(json.dumps({"env_texture_path": str(tex)}))
        w.get_settings_path = lambda: str(settings_path)
        w.load_settings()
        assert w.env_texture_path == str(tex)

    def test_load_settings_no_path_syncs_from_adv_viewer_instance(self, viewer):
        import json

        w, mw, tmp_path = viewer
        tex = tmp_path / "shared.png"
        tex.write_bytes(b"fake")
        mw._adv_rendering_viewer = SimpleNamespace(env_texture_path=str(tex))
        settings_path = tmp_path / "s.json"
        settings_path.write_text(json.dumps({"env_texture_path": ""}))
        w.get_settings_path = lambda: str(settings_path)
        w.load_settings()
        assert w.env_texture_path == str(tex)
        assert w.line_env_path.text() == str(tex)

    def test_load_settings_no_path_syncs_from_global_shared_path(self, viewer):
        import json

        w, mw, tmp_path = viewer
        tex = tmp_path / "global.png"
        tex.write_bytes(b"fake")
        mw.current_env_texture_path = str(tex)
        settings_path = tmp_path / "s.json"
        settings_path.write_text(json.dumps({"env_texture_path": ""}))
        w.get_settings_path = lambda: str(settings_path)
        w.load_settings()
        assert w.env_texture_path == str(tex)


# ===========================================================================
# CubeViewerWidget — load_preset_settings extra branches
# ===========================================================================


class TestCubeViewerLoadPresetSettingsBranches:
    @pytest.fixture
    def viewer(self, qapp, tmp_path):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w, tmp_path
        _teardown_viewer(w)

    def test_load_preset_settings_style_case_insensitive_match(self, viewer):
        w, _tmp_path = viewer
        w.load_preset_settings({"style": "wireframe"})
        assert w.combo_style.currentText() == "Wireframe"

    def test_load_preset_settings_loads_existing_env_texture(self, viewer):
        w, tmp_path = viewer
        tex = tmp_path / "preset_tex.png"
        tex.write_bytes(b"fake")
        w.load_preset_settings({"env_texture_path": str(tex)})
        assert w.env_texture_path == str(tex)
        assert w.line_env_path.text() == str(tex)


# ===========================================================================
# CubeViewerWidget — update_iso remaining branches (real widget)
# ===========================================================================


class TestCubeViewerUpdateIsoRealWidget:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w, grid
        _teardown_viewer(w)

    def test_smoothed_surface_style_smooths_mesh(self, viewer):
        w, grid = viewer
        w.combo_style.setCurrentText("Smoothed Surface")
        mesh = grid.contour.return_value
        mesh.smooth.assert_called()

    def test_remove_actor_exception_is_swallowed(self, viewer):
        w, grid = viewer
        w.update_iso()  # first call populates iso_actor_p/n
        w.plotter.remove_actor.side_effect = RuntimeError("stale actor")
        w.update_iso()  # second call hits the guarded removal try/except

    def test_contour_exception_is_caught(self, viewer):
        w, grid = viewer
        grid.contour.side_effect = RuntimeError("boom")
        w.update_iso()  # must not raise


# ===========================================================================
# CubeViewerWidget — on_remove_env_texture branches
# ===========================================================================


class TestCubeViewerRemoveEnvTextureBranches:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w
        _teardown_viewer(w)

    def test_fallback_to_set_environment_texture_none(self, viewer):
        class SetOnlyPlotter:
            __getattr__ = _fallback_plotter_getattr({"remove_environment_texture"})

            def __init__(self):
                self.calls = []

            def set_environment_texture(self, tex):
                self.calls.append(tex)

        plotter = SetOnlyPlotter()
        viewer.plotter = plotter
        viewer.on_remove_env_texture()
        assert plotter.calls == [None]
        assert viewer.env_texture_path == ""

    def test_fallback_set_environment_texture_exception_is_swallowed(self, viewer):
        class RaisingSetOnlyPlotter:
            __getattr__ = _fallback_plotter_getattr({"remove_environment_texture"})

            def set_environment_texture(self, tex):
                raise RuntimeError("boom")

        viewer.plotter = RaisingSetOnlyPlotter()
        viewer.on_remove_env_texture()  # inner except swallows, still finishes
        assert viewer.env_texture_path == ""

    def test_outer_exception_is_caught(self, viewer):
        class RenderFailsPlotter:
            __getattr__ = _fallback_plotter_getattr(set())

            def __init__(self):
                self.should_fail = True

            def remove_environment_texture(self):
                pass

            def render(self):
                if self.should_fail:
                    # Fail only once: a still-pending initial_update() QTimer
                    # from __init__ may call render() again much later (once
                    # some other test's processEvents() flushes the queue) —
                    # that stray call must not blow up an unrelated test.
                    self.should_fail = False
                    raise RuntimeError("render fail")

        viewer.plotter = RenderFailsPlotter()
        viewer.on_remove_env_texture()  # must not raise


# ===========================================================================
# CubeViewerWidget — load_env_texture (full branch coverage)
# ===========================================================================


class TestCubeViewerLoadEnvTexture:
    @pytest.fixture
    def viewer(self, qapp, tmp_path):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w, mw, tmp_path
        _teardown_viewer(w)

    def test_missing_file_raises_and_shows_critical(self, viewer, monkeypatch):
        w, mw, tmp_path = viewer
        box = MagicMock()
        monkeypatch.setattr(_cube_adv, "QMessageBox", box)
        with pytest.raises(FileNotFoundError):
            w.load_env_texture(str(tmp_path / "nope.png"))
        box.critical.assert_called_once()
        assert w.env_texture_path == ""

    def test_space_in_path_raises_valueerror(self, viewer, monkeypatch):
        w, mw, tmp_path = viewer
        p = tmp_path / "has space.png"
        p.write_bytes(b"x")
        box = MagicMock()
        monkeypatch.setattr(_cube_adv, "QMessageBox", box)
        with pytest.raises(ValueError):
            w.load_env_texture(str(p))
        box.critical.assert_called_once()

    def test_wrong_extension_raises_valueerror(self, viewer, monkeypatch):
        w, mw, tmp_path = viewer
        p = tmp_path / "notatexture.txt"
        p.write_text("hi")
        box = MagicMock()
        monkeypatch.setattr(_cube_adv, "QMessageBox", box)
        with pytest.raises(ValueError):
            w.load_env_texture(str(p))
        box.critical.assert_called_once()

    def test_large_file_user_declines(self, viewer, monkeypatch):
        w, mw, tmp_path = viewer
        p = tmp_path / "big.png"
        p.write_bytes(b"x")
        monkeypatch.setattr(_cube_adv.os.path, "getsize", lambda _p: 20 * 1024 * 1024)
        from PyQt6.QtWidgets import QMessageBox

        monkeypatch.setattr(
            QMessageBox, "question",
            staticmethod(lambda *a, **kw: QMessageBox.StandardButton.No),
        )
        w.load_env_texture(str(p))  # returns early, no crash
        assert w.env_texture_path != str(p)

    def test_large_file_user_continues(self, viewer, monkeypatch):
        w, mw, tmp_path = viewer
        p = tmp_path / "big2.png"
        p.write_bytes(b"x")
        monkeypatch.setattr(_cube_adv.os.path, "getsize", lambda _p: 20 * 1024 * 1024)
        from PyQt6.QtWidgets import QMessageBox

        monkeypatch.setattr(
            QMessageBox, "question",
            staticmethod(lambda *a, **kw: QMessageBox.StandardButton.Yes),
        )
        w.load_env_texture(str(p))
        assert w.env_texture_path == str(p)

    def test_progress_cancelled_before_read_texture(self, viewer, monkeypatch, capsys):
        w, mw, tmp_path = viewer
        p = tmp_path / "cancel.png"
        p.write_bytes(b"x")
        from PyQt6.QtWidgets import QProgressDialog

        monkeypatch.setattr(QProgressDialog, "wasCanceled", lambda self: True)
        w.load_env_texture(str(p))
        assert "cancelled" in capsys.readouterr().out

    def test_progress_cancelled_after_read_texture(self, viewer, monkeypatch, capsys):
        w, mw, tmp_path = viewer
        p = tmp_path / "cancel2.png"
        p.write_bytes(b"x")
        from PyQt6.QtWidgets import QProgressDialog

        calls = {"n": 0}

        def fake_cancel(_self):
            calls["n"] += 1
            return calls["n"] > 1  # not cancelled first check, cancelled afterwards

        monkeypatch.setattr(QProgressDialog, "wasCanceled", fake_cancel)
        w.load_env_texture(str(p))
        assert "cancelled" in capsys.readouterr().out

    def test_read_texture_raises_shows_critical(self, viewer, monkeypatch):
        w, mw, tmp_path = viewer
        p = tmp_path / "readerr.png"
        p.write_bytes(b"x")
        box = MagicMock()
        monkeypatch.setattr(_cube_adv, "QMessageBox", box)
        monkeypatch.setattr(_cube_adv.pv, "read_texture", MagicMock(side_effect=RuntimeError("bad tex")))
        with pytest.raises(ValueError):
            w.load_env_texture(str(p))
        box.critical.assert_called_once()

    def test_set_environment_texture_raises_shows_critical(self, viewer, monkeypatch):
        w, mw, tmp_path = viewer
        p = tmp_path / "seterr.png"
        p.write_bytes(b"x")
        box = MagicMock()
        monkeypatch.setattr(_cube_adv, "QMessageBox", box)
        w.plotter.set_environment_texture = MagicMock(side_effect=RuntimeError("vtk fail"))
        with pytest.raises(ValueError):
            w.load_env_texture(str(p))
        box.critical.assert_called_once()

    def test_successful_load_updates_state_and_syncs_mw(self, viewer):
        w, mw, tmp_path = viewer
        p = tmp_path / "ok.png"
        p.write_bytes(b"x")
        w.load_env_texture(str(p))
        assert w.env_texture_path == str(p)
        assert w.line_env_path.text() == str(p)
        assert mw.current_env_texture_path == str(p)

    def test_successful_load_syncs_adv_viewer_instance(self, viewer):
        w, mw, tmp_path = viewer
        mw._adv_rendering_viewer = SimpleNamespace(env_texture_path="")
        p = tmp_path / "ok2.png"
        p.write_bytes(b"x")
        w.load_env_texture(str(p))
        assert mw._adv_rendering_viewer.env_texture_path == str(p)


# ===========================================================================
# CubeViewerWidget — on_texture_path_entered (found-file branch)
# ===========================================================================


class TestCubeViewerTexturePathEnteredFound:
    @pytest.fixture
    def viewer(self, qapp, tmp_path):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w, tmp_path
        _teardown_viewer(w)

    def test_existing_file_triggers_load(self, viewer):
        w, tmp_path = viewer
        p = tmp_path / "entered.png"
        p.write_bytes(b"x")
        w.line_env_path.setText(str(p))
        w.on_texture_path_entered()
        assert w.env_texture_path == str(p)

    def test_load_error_is_caught(self, viewer, monkeypatch):
        w, tmp_path = viewer
        p = tmp_path / "entered_bad.png"
        p.write_bytes(b"x")
        w.line_env_path.setText(str(p))
        monkeypatch.setattr(w, "load_env_texture", MagicMock(side_effect=RuntimeError("boom")))
        w.on_texture_path_entered()  # must not raise


# ===========================================================================
# CubeViewerWidget — on_reset_all stale-flag disable branches / exception
# ===========================================================================


class TestCubeViewerResetAllExtra:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w
        _teardown_viewer(w)

    def test_reset_disables_stale_effect_flags(self, viewer):
        # Checkboxes already unchecked, so on_reset_all's setChecked(False)
        # calls are no-ops (no signal emitted since state doesn't change) —
        # this simulates stale use_X flags and exercises the manual disable
        # branches inside on_reset_all directly.
        viewer.use_ssao = True
        viewer.use_depth_peeling = True
        viewer.use_edl = True
        viewer.use_aa = True
        viewer.use_shadows = True
        viewer.on_reset_all()
        assert not viewer.use_ssao
        assert not viewer.use_depth_peeling
        assert not viewer.use_edl
        assert not viewer.use_aa
        assert not viewer.use_shadows

    def test_reset_all_exception_is_caught_and_signals_unblocked(self, viewer):
        viewer.btn_color_p = None  # AttributeError inside the try block
        viewer.on_reset_all()  # must not raise
        assert not viewer.signalsBlocked()


# ===========================================================================
# CubeViewerWidget — close_plugin / _detach_stale_view effect branches
# ===========================================================================


class TestCubeViewerCloseAndDetachEffects:
    def test_close_plugin_disables_effects_and_restores_ui(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        w._structure_watch_timer.stop()
        w.use_edl = True
        w.use_shadows = True
        w.use_ssao = True
        mw.init_manager = SimpleNamespace(current_file_path="x.cube")
        mw.ui_manager = SimpleNamespace(restore_ui_for_editing=MagicMock())
        w.close_plugin()
        mw.ui_manager.restore_ui_for_editing.assert_called_once()
        assert mw.init_manager.current_file_path is None

    def test_detach_stale_view_disables_effects(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        w.use_edl = True
        w.use_shadows = True
        w.use_ssao = True
        w.iso_actor_p = MagicMock()
        mw.statusBar = MagicMock(return_value=MagicMock())
        w._detach_stale_view()
        mw.statusBar.return_value.showMessage.assert_called_once()
        assert w.dock is None


# ===========================================================================
# CubeViewerWidget — effect-toggle exception / gating branches
# ===========================================================================


class TestCubeViewerEffectToggleExceptions:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w
        _teardown_viewer(w)

    def test_ssao_toggle_exception_is_caught(self, viewer):
        viewer.plotter.enable_ssao = MagicMock(side_effect=RuntimeError("boom"))
        viewer.check_ssao.setChecked(True)  # must not raise

    def test_depth_peeling_toggle_exception_is_caught(self, viewer):
        viewer.plotter.enable_depth_peeling = MagicMock(side_effect=RuntimeError("boom"))
        viewer.check_depth.setChecked(True)  # must not raise

    def test_shadows_toggle_exception_is_caught(self, viewer):
        viewer.plotter.enable_shadows = MagicMock(side_effect=RuntimeError("boom"))
        viewer.check_shadows.setChecked(True)  # must not raise

    def test_aa_toggle_exception_is_caught(self, viewer):
        viewer.plotter.enable_anti_aliasing = MagicMock(side_effect=RuntimeError("boom"))
        viewer.check_aa.setChecked(True)  # must not raise

    def test_clean_render_pipeline_returns_early_without_renderer(self, viewer):
        class NoRendererPlotter:
            __getattr__ = _fallback_plotter_getattr({"renderer"})

        viewer.plotter = NoRendererPlotter()
        viewer._clean_render_pipeline()  # returns early, no crash

    def test_light_intensity_updates_real_lights_list(self, viewer):
        light = MagicMock()
        viewer.plotter.renderer.GetLights.return_value = [light]
        viewer.slider_light.setValue(150)
        light.SetIntensity.assert_called_once_with(1.5)

    def test_light_intensity_exception_is_caught(self, viewer):
        viewer.plotter.render.side_effect = RuntimeError("boom")
        viewer.slider_light.setValue(120)  # must not raise


# ===========================================================================
# CubeViewerWidget — _disable_conflicting_effects: mol-redraw / restore /
# final-render branches (reached through the effect-toggle slots)
# ===========================================================================


class TestCubeViewerDisableConflictingEffectsExtra:
    @pytest.fixture
    def viewer(self, qapp):
        w, ctx, mw, dock, grid = _make_viewer()
        yield w, ctx
        _teardown_viewer(w)

    def test_mol_redraw_happens_when_current_molecule_set(self, viewer):
        w, ctx = viewer
        ctx.current_molecule = MagicMock()
        w.check_ssao.setChecked(True)
        ctx.draw_molecule_3d.assert_called()

    def test_mol_redraw_exception_is_caught(self, viewer):
        w, ctx = viewer
        ctx.current_molecule = MagicMock()
        ctx.draw_molecule_3d.side_effect = RuntimeError("boom")
        w.check_ssao.setChecked(True)  # must not raise

    def test_update_iso_exception_inside_disable_conflicting_is_caught(self, viewer, monkeypatch):
        w, ctx = viewer
        monkeypatch.setattr(w, "update_iso", MagicMock(side_effect=RuntimeError("boom")))
        w.check_ssao.setChecked(True)  # must not raise

    def test_edl_and_shadows_disable_restore_exceptions_are_caught(self, viewer):
        w, ctx = viewer
        w.use_edl = True
        w.use_shadows = True
        w.plotter.disable_eye_dome_lighting = MagicMock(side_effect=RuntimeError("boom"))
        w.plotter.disable_shadows = MagicMock(side_effect=RuntimeError("boom"))
        w.plotter.enable_eye_dome_lighting = MagicMock(side_effect=RuntimeError("boom"))
        w.plotter.enable_shadows = MagicMock(side_effect=RuntimeError("boom"))
        w.check_ssao.setChecked(True)  # must not raise despite all the failures

    def test_final_render_exception_is_caught(self, viewer):
        w, ctx = viewer
        w.plotter.render.side_effect = RuntimeError("boom")
        w.check_ssao.setChecked(True)  # must not raise


# ===========================================================================
# CubeViewerWidget — on_preset_activated with a current molecule set
# ===========================================================================


class TestCubeViewerPresetActivatedWithMolecule:
    @pytest.fixture
    def viewer(self, qapp, tmp_path):
        w, ctx, mw, dock, grid = _make_viewer()
        w.get_settings_path = lambda: str(tmp_path / "cube_viewer_advanced.json")
        yield w, ctx
        _teardown_viewer(w)

    def test_activation_redraws_molecule_when_set(self, viewer, monkeypatch):
        w, ctx = viewer
        monkeypatch.setattr(_cube_adv.QInputDialog, "getText", staticmethod(lambda *a, **kw: ("Preset1", True)))
        w.save_preset()
        ctx.current_molecule = MagicMock()
        idx = w.combo_presets.findText("Preset1")
        w.combo_presets.setCurrentIndex(idx)
        w.on_preset_activated(idx)
        ctx.draw_molecule_3d.assert_called()

    def test_activation_mol_redraw_exception_is_caught(self, viewer, monkeypatch):
        w, ctx = viewer
        monkeypatch.setattr(_cube_adv.QInputDialog, "getText", staticmethod(lambda *a, **kw: ("Preset2", True)))
        w.save_preset()
        ctx.current_molecule = MagicMock()
        ctx.draw_molecule_3d.side_effect = RuntimeError("boom")
        idx = w.combo_presets.findText("Preset2")
        w.combo_presets.setCurrentIndex(idx)
        w.on_preset_activated(idx)  # must not raise
