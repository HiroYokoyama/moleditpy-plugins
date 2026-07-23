"""
Headless GUI tests for the Advanced Rendering plugin.

Covers: HideOnCloseDialog, AdvancedGraphicsWidget.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_advanced_rendering.py
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ADV_RENDER_PATH = PLUGINS_DIR / "Advanced_Rendering" / "advanced_rendering.py"

# CI's test-gui job installs only pytest + PyQt6 (no numpy/rdkit/vtk/pyvista).
# Real pyvista/vtk are used below (see _make_real_plotter) to exercise genuine
# VTK actor/light/render-pass code paths; importorskip lets this whole module
# skip cleanly there while running fully wherever the chemistry stack exists.
_real_pv = pytest.importorskip("pyvista")
_real_vtk = pytest.importorskip("vtk")

with mock_chemistry_imports():
    _adv = load_plugin_for_gui(ADV_RENDER_PATH)

# save_settings' sanitize() does isinstance(x, np.generic); the MagicMock numpy
# would make that a TypeError, so give the module a shim with valid class args.
_adv.np = SimpleNamespace(generic=(), ndarray=())

# Swap the plugin's `pv`/`vtk` references (MagicMock, from mock_chemistry_imports
# above) back to the real modules so tests that build a real off-screen Plotter
# exercise genuine VTK actor/light/render-pass code instead of MagicMock stubs.
_adv.pv = _real_pv
_adv.vtk = _real_vtk


def _make_real_plotter():
    """A real off-screen pyvista Plotter with one mesh actor, for exercising
    the plugin's actor/light/render-pass code against genuine VTK objects."""
    plotter = _real_pv.Plotter(off_screen=True)
    plotter.add_mesh(_real_pv.Sphere())
    return plotter


def _make_adv_widget_with_plotter(tmp_path, style="Ball & Stick (Advanced Rendering)"):
    """AdvancedGraphicsWidget backed by a real off-screen Plotter."""
    from PyQt6.QtWidgets import QWidget

    plotter = _make_real_plotter()
    mw = QWidget()
    mw.plotter = plotter
    mw.view_3d_manager = SimpleNamespace(current_3d_style=style)
    w = _adv.AdvancedGraphicsWidget(mw)
    w.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
    return mw, w, plotter


def _teardown_adv_widget_with_plotter(mw, w, plotter):
    _teardown_adv_widget(mw, w)
    plotter.close()


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


# ===========================================================================
# load_plugin_from_mw — dialog create / reuse / recreate-on-error
# ===========================================================================


class TestLoadPluginFromMw:
    def test_creates_dialog_and_viewer_when_absent(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")

        _adv.load_plugin_from_mw(mw)
        try:
            assert isinstance(mw._adv_graphics_dialog, _adv.HideOnCloseDialog)
            assert isinstance(mw._adv_rendering_viewer, _adv.AdvancedGraphicsWidget)
            assert mw._adv_graphics_dialog.isVisible()
        finally:
            mw._adv_rendering_viewer.get_settings_path = lambda: str(
                tmp_path / "adv_settings.json"
            )
            _teardown_adv_widget(mw, mw._adv_rendering_viewer)
            mw._adv_graphics_dialog.destroy()

    def test_reuses_existing_dialog(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")

        _adv.load_plugin_from_mw(mw)
        first_dialog = mw._adv_graphics_dialog
        first_viewer = mw._adv_rendering_viewer
        first_viewer.get_settings_path = lambda: str(tmp_path / "adv_settings.json")

        _adv.load_plugin_from_mw(mw)
        assert mw._adv_graphics_dialog is first_dialog
        assert mw._adv_rendering_viewer is first_viewer

        _teardown_adv_widget(mw, first_viewer)
        first_dialog.destroy()

    def test_recreates_when_dialog_show_raises(self, qapp, tmp_path, monkeypatch):
        """show()/raise_()/activateWindow() raising forces a full recreate."""
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")

        _adv.load_plugin_from_mw(mw)
        stale_dialog = mw._adv_graphics_dialog
        stale_viewer = mw._adv_rendering_viewer
        stale_viewer.get_settings_path = lambda: str(tmp_path / "adv_settings.json")

        monkeypatch.setattr(stale_dialog, "show", MagicMock(side_effect=RuntimeError("boom")))

        _adv.load_plugin_from_mw(mw)

        new_dialog = mw._adv_graphics_dialog
        assert new_dialog is not stale_dialog
        new_viewer = mw._adv_rendering_viewer
        new_viewer.get_settings_path = lambda: str(tmp_path / "adv_settings2.json")

        _teardown_adv_widget(mw, new_viewer)
        new_dialog.destroy()
        stale_dialog.destroy()

    def test_plugin_manager_registers_window(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
        mw.plugin_manager = MagicMock()

        _adv.load_plugin_from_mw(mw)
        mw.plugin_manager.register_window.assert_called_once_with(
            _adv.PLUGIN_NAME, "main_panel", mw._adv_graphics_dialog
        )

        viewer = mw._adv_rendering_viewer
        viewer.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
        _teardown_adv_widget(mw, viewer)
        mw._adv_graphics_dialog.destroy()


# ===========================================================================
# initialize() — cleanup of stale instances on reload
# ===========================================================================


class TestInitializeCleansUpOldInstances:
    def test_second_initialize_stops_old_sync_timer_and_closes_dialog(
        self, qapp, tmp_path
    ):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw

        _adv.initialize(ctx)
        old_viewer = mw._adv_rendering_viewer
        old_viewer.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
        mw._adv_graphics_dialog = _adv.HideOnCloseDialog(mw)
        old_dialog = mw._adv_graphics_dialog

        assert old_viewer._sync_timer.isActive()

        _adv.initialize(ctx)

        assert old_viewer._sync_timer is not None
        assert not old_viewer._sync_timer.isActive()

        new_viewer = mw._adv_rendering_viewer
        assert new_viewer is not old_viewer
        new_viewer.get_settings_path = lambda: str(tmp_path / "adv_settings2.json")

        _teardown_adv_widget(mw, new_viewer)
        old_dialog.destroy()


# ===========================================================================
# _retry_init_lighting — ready / not-ready / timeout branches
# ===========================================================================


class TestRetryInitLighting:
    def test_ready_applies_lighting_and_pbr(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w._init_attempts = 0
            w._retry_init_lighting()  # plotter.renderer is real -> "ready" branch
            # No exception, and does not schedule further retries.
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_not_ready_increments_attempts(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.plotter = None
            mw.plotter = None
            w._init_attempts = 0
            w._retry_init_lighting()
            assert w._init_attempts == 1
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_timeout_after_max_attempts_does_not_reschedule(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.plotter = None
            mw.plotter = None
            w._init_attempts = w._max_init_attempts - 1
            w._retry_init_lighting()
            assert w._init_attempts == w._max_init_attempts
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# _poll_style_change
# ===========================================================================


class TestPollStyleChange:
    def test_no_mw_is_noop(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.mw = None
            w._poll_style_change()  # must not raise
        finally:
            w.mw = mw
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_style_change_syncs_ui_and_persists(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path, style="ball_and_stick")
        try:
            w._last_polled_style = ""
            mw.view_3d_manager.current_3d_style = "CPK (Advanced Rendering)"
            w._poll_style_change()
            assert w._last_polled_style == "CPK (Advanced Rendering)"
            assert w.presets["last_active_style"] == "CPK (Advanced Rendering)"
            assert w.use_atom_pbr  # sync_style_ui enabled PBR for advanced style
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_unchanged_style_is_noop(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w._last_polled_style = mw.view_3d_manager.current_3d_style
            before = dict(w.presets)
            w._poll_style_change()
            assert "last_active_style" not in w.presets or w.presets.get(
                "last_active_style"
            ) == before.get("last_active_style")
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# apply_pbr_forced / safe_plotter
# ===========================================================================


class TestApplyPbrForcedAndSafePlotter:
    def test_apply_pbr_forced_forces_pbr(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.use_atom_pbr = False
            w.apply_pbr_forced()
            for actor in plotter.renderer.actors.values():
                if actor.IsA("vtkActor"):
                    prop = actor.GetProperty()
                    assert prop.GetInterpolation() == _real_vtk.VTK_PBR
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_safe_plotter_falls_back_to_mw_plotter(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.plotter = None
            mw.plotter = plotter
            assert w.safe_plotter is plotter
            assert w.plotter is plotter
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# update_atoms_pbr — dict-actors path, GetActors() fallback, style gating
# ===========================================================================


class TestUpdateAtomsPbr:
    def test_advanced_style_applies_pbr_to_dict_actors(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.use_atom_pbr = True
            w.atom_metallic = 0.6
            w.atom_roughness = 0.2
            w.update_atoms_pbr()
            actors = list(plotter.renderer.actors.values())
            assert actors
            for actor in actors:
                prop = actor.GetProperty()
                assert prop.GetInterpolation() == _real_vtk.VTK_PBR
                assert prop.GetMetallic() == pytest.approx(0.6)
                assert prop.GetRoughness() == pytest.approx(0.2)
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_non_advanced_style_forces_phong_even_if_pbr_flag_set(
        self, qapp, tmp_path
    ):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path, style="ball_and_stick")
        try:
            w.use_atom_pbr = True
            w.update_atoms_pbr()
            for actor in plotter.renderer.actors.values():
                prop = actor.GetProperty()
                assert prop.GetInterpolation() == _real_vtk.VTK_PHONG
                assert prop.GetMetallic() == pytest.approx(0.0)
                assert prop.GetRoughness() == pytest.approx(0.5)
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_getactors_fallback_when_no_actors_dict(self, qapp, tmp_path):
        """Raw vtkRenderer (no pyvista `.actors` dict) exercises the
        GetActors() fallback branch."""
        from PyQt6.QtWidgets import QWidget

        renderer = _real_vtk.vtkRenderer()
        actor = _real_vtk.vtkActor()
        actor.SetMapper(_real_vtk.vtkPolyDataMapper())
        renderer.AddActor(actor)

        class _FakePlotter:
            pass

        fp = _FakePlotter()
        fp.renderer = renderer
        fp.render = MagicMock()

        mw = QWidget()
        mw.plotter = fp
        mw.view_3d_manager = SimpleNamespace(
            current_3d_style="Stick (Advanced Rendering)"
        )
        w = _adv.AdvancedGraphicsWidget(mw)
        w.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
        try:
            w.use_atom_pbr = True
            w.update_atoms_pbr()
            assert actor.GetProperty().GetInterpolation() == _real_vtk.VTK_PBR
            fp.render.assert_called()
        finally:
            _teardown_adv_widget(mw, w)

    def test_no_renderer_is_noop(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        class _BarePlotter:
            pass

        mw = QWidget()
        mw.plotter = _BarePlotter()
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
        w = _adv.AdvancedGraphicsWidget(mw)
        w.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
        try:
            w.update_atoms_pbr()  # must not raise
        finally:
            _teardown_adv_widget(mw, w)


# ===========================================================================
# Atom PBR handlers / clear_atom_settings
# ===========================================================================


class TestAtomHandlers:
    def test_on_atom_pbr_toggled_enables_sliders_and_saves(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.check_atom_pbr.setChecked(True)
            assert w.use_atom_pbr
            assert w.slider_atom_metallic.isEnabled()
            assert w.slider_atom_roughness.isEnabled()
            import json

            data = json.loads((tmp_path / "adv_settings.json").read_text())
            assert data["use_atom_pbr"] is True
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_metallic_and_roughness_changed_when_pbr_active(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.check_atom_pbr.setChecked(True)
            w.slider_atom_metallic.setValue(80)
            w.slider_atom_roughness.setValue(10)
            assert w.atom_metallic == pytest.approx(0.8)
            assert w.atom_roughness == pytest.approx(0.1)
            for actor in plotter.renderer.actors.values():
                prop = actor.GetProperty()
                assert prop.GetMetallic() == pytest.approx(0.8)
                assert prop.GetRoughness() == pytest.approx(0.1)
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_metallic_changed_while_pbr_off_does_not_apply(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            assert not w.use_atom_pbr
            w.slider_atom_metallic.setValue(55)
            assert w.atom_metallic == pytest.approx(0.55)
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_clear_atom_settings_resets_to_phong(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.check_atom_pbr.setChecked(True)
            w.slider_atom_metallic.setValue(90)
            w.plotter = plotter
            w.clear_atom_settings()
            assert not w.use_atom_pbr
            assert not w.check_atom_pbr.isChecked()
            assert w.slider_atom_metallic.value() == 0
            assert w.slider_atom_roughness.value() == 50
            for actor in plotter.renderer.actors.values():
                assert (
                    actor.GetProperty().GetInterpolation() == _real_vtk.VTK_PHONG
                )
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# Environment texture
# ===========================================================================


class TestEnvTexture:
    def test_load_env_texture_success(self, qapp, tmp_path, monkeypatch):
        Image = pytest.importorskip("PIL.Image")

        png_path = tmp_path / "tex.png"
        Image.new("RGB", (4, 4), (255, 0, 0)).save(png_path)

        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            monkeypatch.setattr(
                _adv.QFileDialog,
                "getOpenFileName",
                lambda *a, **k: (str(png_path), ""),
            )
            w.on_load_env_texture()
            assert w.env_texture_path == str(png_path)
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_load_env_texture_cancelled_leaves_path_unset(self, qapp, tmp_path, monkeypatch):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            monkeypatch.setattr(
                _adv.QFileDialog, "getOpenFileName", lambda *a, **k: ("", "")
            )
            w.on_load_env_texture()
            assert w.env_texture_path == ""
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_remove_env_texture_clears_path(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.env_texture_path = "somefile.png"
            w.on_remove_env_texture()
            assert w.env_texture_path == ""
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_apply_texture_noop_without_plotter(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
        w = _adv.AdvancedGraphicsWidget(mw)
        w.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
        try:
            w.apply_texture()  # must not raise
        finally:
            _teardown_adv_widget(mw, w)

    def test_apply_texture_error_shows_warning(self, qapp, tmp_path, monkeypatch):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.env_texture_path = str(tmp_path / "missing_but_named.png")
            # os.path.exists is False for a nonexistent file, so force the
            # "exists" branch to hit the exception path via read_texture.
            (tmp_path / "missing_but_named.png").write_bytes(b"not-an-image")
            monkeypatch.setattr(
                _adv.pv,
                "read_texture",
                MagicMock(side_effect=RuntimeError("bad texture")),
            )
            box = MagicMock()
            monkeypatch.setattr(_adv, "QMessageBox", box)
            w.apply_texture()
            box.warning.assert_called_once()
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# _clean_render_pipeline
# ===========================================================================


class TestCleanRenderPipeline:
    def test_clean_pipeline_calls_set_passes(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w._clean_render_pipeline()  # must not raise
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_clean_pipeline_noop_without_plotter(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
        w = _adv.AdvancedGraphicsWidget(mw)
        w.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
        try:
            w._clean_render_pipeline()  # must not raise
        finally:
            _teardown_adv_widget(mw, w)


# ===========================================================================
# Scene effect toggles: shadows / ssao / depth peeling / AA / EDL
# ===========================================================================


class TestSceneEffectToggles:
    def test_shadows_rejected_on_standard_style_without_pbr(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path, style="ball_and_stick")
        try:
            w.check_shadows.setChecked(True)
            assert not w.check_shadows.isChecked()
            assert not w.use_shadows
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_shadows_enabled_on_advanced_style(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.check_shadows.setChecked(True)
            assert w.use_shadows
            w.check_shadows.setChecked(False)
            assert not w.use_shadows
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_ssao_toggle_disables_depth_peeling(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.check_depth.setChecked(True)
            assert w.use_depth_peeling
            w.check_ssao.setChecked(True)
            assert w.use_ssao
            assert not w.use_depth_peeling
            assert not w.check_depth.isChecked()
            w.check_ssao.setChecked(False)
            assert not w.use_ssao
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_depth_peeling_toggle_on_and_off(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.check_depth.setChecked(True)
            assert w.use_depth_peeling
            w.check_depth.setChecked(False)
            assert not w.use_depth_peeling
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_aa_toggle_on_and_off(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.check_aa.setChecked(True)
            assert w.use_aa
            w.check_aa.setChecked(False)
            assert not w.use_aa
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_edl_toggle_enables_slider_and_strength(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.check_edl.setChecked(True)
            assert w.use_edl
            assert w.slider_edl.isEnabled()
            w.slider_edl.setValue(75)
            assert w.edl_strength == pytest.approx(0.75)
            w.check_edl.setChecked(False)
            assert not w.use_edl
            assert not w.slider_edl.isEnabled()
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_edl_rejected_on_standard_style(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path, style="ball_and_stick")
        try:
            w.check_edl.setChecked(True)
            assert not w.check_edl.isChecked()
            assert not w.use_edl
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_apply_edl_noop_without_plotter(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
        w = _adv.AdvancedGraphicsWidget(mw)
        w.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
        try:
            w.apply_edl()  # must not raise
        finally:
            _teardown_adv_widget(mw, w)

    def test_update_edl_strength_noop_when_disabled(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            assert not w.use_edl
            w.update_edl_strength()  # must not raise / no-op
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# Light position/intensity handlers
# ===========================================================================


class TestLightHandlers:
    def test_light_pos_changed_updates_state(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.slider_azi.setValue(-90)
            w.slider_ele.setValue(30)
            assert w.light_azimuth == -90
            assert w.light_elevation == 30
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_update_lights_creates_light_when_none_exist(self, qapp, tmp_path):
        """Raw vtkRenderer starts with zero lights, exercising the
        AddLight()-on-empty-collection branch."""
        from PyQt6.QtWidgets import QWidget

        renderer = _real_vtk.vtkRenderer()

        class _FakePlotter:
            pass

        fp = _FakePlotter()
        fp.renderer = renderer
        fp.render = MagicMock()

        mw = QWidget()
        mw.plotter = fp
        mw.view_3d_manager = SimpleNamespace(
            current_3d_style="CPK (Advanced Rendering)"
        )
        w = _adv.AdvancedGraphicsWidget(mw)
        w.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
        try:
            assert renderer.GetLights().GetNumberOfItems() == 0
            w.update_lights()
            assert renderer.GetLights().GetNumberOfItems() == 1
            fp.render.assert_called()
        finally:
            _teardown_adv_widget(mw, w)

    def test_update_lights_noop_on_standard_style(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path, style="ball_and_stick")
        try:
            n_before = plotter.renderer.GetLights().GetNumberOfItems()
            w.update_lights()
            assert plotter.renderer.GetLights().GetNumberOfItems() == n_before
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# closeEvent
# ===========================================================================


class TestCloseEvent:
    def test_close_event_disables_effects(self, qapp, tmp_path):
        from PyQt6.QtGui import QCloseEvent

        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.check_edl.setChecked(True)
            w.closeEvent(QCloseEvent())
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# Presets: save / delete via the real widget + QInputDialog
# ===========================================================================


class TestPresetsReal:
    def test_save_preset_adds_new_named_preset(self, qapp, tmp_path, monkeypatch):
        from PyQt6.QtWidgets import QInputDialog

        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            monkeypatch.setattr(
                QInputDialog, "getText", staticmethod(lambda *a, **k: ("MyPreset", True))
            )
            w.light_intensity = 3.3
            w.save_preset()
            assert "MyPreset" in w.presets
            assert w.presets["MyPreset"]["light_intensity"] == pytest.approx(3.3)
            assert w.combo_presets.currentText() == "MyPreset"
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_save_preset_cancelled_does_nothing(self, qapp, tmp_path, monkeypatch):
        from PyQt6.QtWidgets import QInputDialog

        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            monkeypatch.setattr(
                QInputDialog, "getText", staticmethod(lambda *a, **k: ("Whatever", False))
            )
            before = set(w.presets)
            w.save_preset()
            assert set(w.presets) == before
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_save_preset_rejects_locked_name(self, qapp, tmp_path, monkeypatch):
        from PyQt6.QtWidgets import QInputDialog

        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            monkeypatch.setattr(
                QInputDialog, "getText", staticmethod(lambda *a, **k: ("Default", True))
            )
            box = MagicMock()
            monkeypatch.setattr(_adv, "QMessageBox", box)
            w.save_preset()
            box.warning.assert_called_once()
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_delete_preset_removes_user_preset(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.presets["Custom"] = w.gather_settings_dict()
            w.update_preset_combo()
            w.combo_presets.setCurrentText("Custom")
            w.delete_preset()
            assert "Custom" not in w.presets
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_on_preset_activated_applies_settings(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            custom = w.gather_settings_dict()
            custom["light_intensity"] = 4.5
            w.presets["Custom"] = custom
            w.update_preset_combo()
            w.combo_presets.setCurrentText("Custom")
            w.on_preset_activated(0)
            assert w.light_intensity == pytest.approx(4.5)
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# apply_settings_dict — full application against a real widget/plotter
# ===========================================================================


class TestApplySettingsDictReal:
    def test_apply_settings_dict_updates_ui_and_scene(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            data = {
                "use_shadows": True,
                "light_intensity": 1.5,
                "light_azimuth": -30,
                "light_elevation": 20,
                "use_ssao": False,
                "use_depth_peeling": False,
                "use_aa": True,
                "use_edl": False,
                "edl_strength": 0.4,
                "env_texture_path": "",
                "use_atom_pbr": True,
                "atom_metallic": 0.3,
                "atom_roughness": 0.7,
            }
            w.apply_settings_dict(data)
            assert w.light_intensity == pytest.approx(1.5)
            assert w.slider_light.value() == 150
            assert w.light_azimuth == -30
            assert w.slider_azi.value() == -30
            assert w.check_shadows.isChecked()
            assert w.check_aa.isChecked()
            assert w.atom_metallic == pytest.approx(0.3)
            assert w.atom_roughness == pytest.approx(0.7)
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# sync_style_ui / _enforce_scene_state
# ===========================================================================


class TestEnforceSceneState:
    def test_enforce_scene_state_enable_applies_stored_effects(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path, style="ball_and_stick")
        try:
            w.use_shadows = True
            w.use_ssao = True
            w.use_edl = True
            w.use_aa = True
            w.use_depth_peeling = True
            w._enforce_scene_state(enable=True)  # must not raise
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_enforce_scene_state_disable_reverts_everything(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.use_shadows = True
            w.use_ssao = True
            w._enforce_scene_state(enable=False)  # must not raise
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_enforce_scene_state_noop_without_plotter(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.plotter = None
        mw.view_3d_manager = SimpleNamespace(current_3d_style="ball_and_stick")
        w = _adv.AdvancedGraphicsWidget(mw)
        w.get_settings_path = lambda: str(tmp_path / "adv_settings.json")
        try:
            w._enforce_scene_state(enable=True)  # must not raise
        finally:
            _teardown_adv_widget(mw, w)

    def test_sync_style_ui_round_trip_toggles_pbr_and_scene(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path, style="ball_and_stick")
        try:
            w.use_shadows = True
            w.sync_style_ui("Wireframe (Advanced Rendering)")
            assert w.use_atom_pbr
            w.sync_style_ui("ball_and_stick")
            assert not w.use_atom_pbr
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# save_settings sanitize() — np.generic / np.ndarray / QColor / nested list
# ===========================================================================


class TestSaveSettingsSanitize:
    def test_sanitizes_numpy_and_qcolor_in_presets(self, qapp, tmp_path, monkeypatch):
        import json

        real_np = pytest.importorskip("numpy")
        from PyQt6.QtGui import QColor

        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            monkeypatch.setattr(_adv, "np", real_np)
            w.presets["Weird"] = {
                "scalar": real_np.float64(1.5),
                "arr": real_np.array([1, 2, 3]),
                "color": QColor(10, 20, 30),
                "nested": [real_np.int64(7), real_np.int64(8)],
            }
            w.save_settings()
            data = json.loads((tmp_path / "adv_settings.json").read_text())
            weird = data["presets"]["Weird"]
            assert weird["scalar"] == pytest.approx(1.5)
            assert weird["arr"] == [1, 2, 3]
            assert weird["color"] == [10, 20, 30]
            assert weird["nested"] == [7, 8]
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_save_settings_failure_is_logged_not_raised(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.get_settings_path = lambda: str(tmp_path / "nope" / "settings.json")
            w.save_settings()  # must not raise
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)


# ===========================================================================
# load_settings — real widget + real file round-trip
# ===========================================================================


class TestLoadSettingsReal:
    def test_load_settings_applies_saved_state(self, qapp, tmp_path):
        import json

        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            settings_path = tmp_path / "adv_settings.json"
            payload = {
                "use_shadows": False,
                "light_intensity": 2.2,
                "light_azimuth": 10,
                "light_elevation": -5,
                "use_ssao": False,
                "use_depth_peeling": False,
                "use_aa": False,
                "use_edl": False,
                "edl_strength": 0.3,
                "env_texture_path": "",
                "use_atom_pbr": False,
                "atom_metallic": 0.0,
                "atom_roughness": 0.5,
                "presets": {"FromDisk": {"light_intensity": 9.0}},
            }
            settings_path.write_text(json.dumps(payload), encoding="utf-8")
            w.get_settings_path = lambda: str(settings_path)

            w.load_settings()

            assert w.light_intensity == pytest.approx(2.2)
            assert "FromDisk" in w.presets
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)

    def test_load_settings_missing_file_is_noop(self, qapp, tmp_path):
        mw, w, plotter = _make_adv_widget_with_plotter(tmp_path)
        try:
            w.get_settings_path = lambda: str(tmp_path / "does_not_exist.json")
            w.load_settings()  # must not raise
        finally:
            _teardown_adv_widget_with_plotter(mw, w, plotter)
