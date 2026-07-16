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

with mock_chemistry_imports():
    _adv = load_plugin_for_gui(ADV_RENDER_PATH)

# save_settings' sanitize() does isinstance(x, np.generic); the MagicMock numpy
# would make that a TypeError, so give the module a shim with valid class args.
_adv.np = SimpleNamespace(generic=(), ndarray=())


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
