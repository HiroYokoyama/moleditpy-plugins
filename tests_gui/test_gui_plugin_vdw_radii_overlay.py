"""
Headless GUI tests for the VDW Radii Overlay plugin.

Covers: VDWConfigWindow.

Chemistry libs (pyvista, numpy, rdkit, …) are mocked; real PyQt6 is used.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

VDW_PATH = PLUGINS_DIR / "VDW_Radii_Overlay" / "vdw_radii_overlay.py"

with mock_chemistry_imports():
    _vdw = load_plugin_for_gui(VDW_PATH)


def _vdw_context() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


# ===========================================================================
# VDWConfigWindow  (VDW Radii Overlay)
# ===========================================================================


class TestVDWConfigWindow:
    """VDWConfigWindow: slider/spinbox settings dialog for VDW overlay."""

    @pytest.fixture
    def win(self, qapp, tmp_path, monkeypatch):
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "vdw.json"))
        # Isolate the module-global settings dict per test so mutations from
        # interactive tests (below) don't leak across tests.
        monkeypatch.setattr(
            _vdw,
            "_vdw_settings",
            {"occupancy": 0.3, "resolution": 0.125, "base_style": "default"},
        )
        ctx = _vdw_context()
        w = _vdw.VDWConfigWindow(context=ctx)
        yield w
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "VDW Overlay Settings"

    def test_is_non_modal(self, win):
        assert not win.isModal()

    def test_occupancy_slider_exists(self, win):
        assert win.slider_occ is not None

    def test_occupancy_spinbox_range(self, win):
        assert win.spin_occ.minimum() == 0.0
        assert win.spin_occ.maximum() == 1.0

    def test_resolution_slider_exists(self, win):
        assert win.slider_res is not None

    def test_resolution_spinbox_range(self, win):
        assert win.spin_res.minimum() == pytest.approx(0.05)
        assert win.spin_res.maximum() == pytest.approx(0.50)

    def test_default_occupancy_matches_settings(self, win):
        expected = _vdw._vdw_settings.get("occupancy", 0.3)
        assert win.spin_occ.value() == pytest.approx(expected)

    def test_default_resolution_matches_settings(self, win):
        expected = _vdw._vdw_settings.get("resolution", 0.125)
        assert win.spin_res.value() == pytest.approx(expected)

    def test_size_is_reasonable(self, win):
        assert win.width() >= 200
        assert win.height() >= 100

    def test_occupancy_slider_updates_spin_and_settings(self, win):
        win.slider_occ.setValue(80)
        assert win.spin_occ.value() == pytest.approx(0.8)
        assert _vdw._vdw_settings["occupancy"] == pytest.approx(0.8)
        win.context.refresh_3d_view.assert_called()

    def test_occupancy_spin_updates_slider_and_settings(self, win):
        win.spin_occ.setValue(0.65)
        assert win.slider_occ.value() == 65
        assert _vdw._vdw_settings["occupancy"] == pytest.approx(0.65)

    def test_resolution_slider_updates_spin_and_settings(self, win):
        win.slider_res.setValue(30)
        assert win.spin_res.value() == pytest.approx(0.3)
        assert _vdw._vdw_settings["resolution"] == pytest.approx(0.3)

    def test_resolution_spin_updates_slider_and_settings(self, win):
        win.spin_res.setValue(0.2)
        assert win.slider_res.value() == 20
        assert _vdw._vdw_settings["resolution"] == pytest.approx(0.2)

    def test_base_style_changed_to_stick(self, win):
        win.combo_base_style.setCurrentIndex(1)
        assert _vdw._vdw_settings["base_style"] == "stick"

    def test_base_style_changed_back_to_default(self, win):
        win.combo_base_style.setCurrentIndex(1)
        win.combo_base_style.setCurrentIndex(0)
        assert _vdw._vdw_settings["base_style"] == "default"

    def test_settings_persist_to_disk_on_change(self, win):
        win.slider_occ.setValue(50)
        on_disk = json.loads(Path(_vdw.SETTINGS_FILE).read_text())
        assert on_disk["occupancy"] == pytest.approx(0.5)

    def test_reset_defaults_restores_all_widgets(self, win):
        win.slider_occ.setValue(90)
        win.slider_res.setValue(45)
        win.combo_base_style.setCurrentIndex(1)

        win.reset_defaults()

        assert win.spin_occ.value() == pytest.approx(0.3)
        assert win.slider_occ.value() == 30
        assert win.spin_res.value() == pytest.approx(0.125)
        assert win.slider_res.value() == 12
        assert win.combo_base_style.currentIndex() == 0
        assert _vdw._vdw_settings["occupancy"] == pytest.approx(0.3)
        assert _vdw._vdw_settings["resolution"] == pytest.approx(0.125)
        assert _vdw._vdw_settings["base_style"] == "default"

    def test_reset_defaults_calls_refresh(self, win):
        win.context.refresh_3d_view.reset_mock()
        win.reset_defaults()
        win.context.refresh_3d_view.assert_called()

    def test_refresh_ui_values_reads_global_settings(self, win):
        _vdw._vdw_settings["occupancy"] = 0.6
        _vdw._vdw_settings["resolution"] = 0.3
        _vdw._vdw_settings["base_style"] = "stick"

        win.refresh_ui_values()

        assert win.spin_occ.value() == pytest.approx(0.6)
        assert win.slider_occ.value() == 60
        assert win.spin_res.value() == pytest.approx(0.3)
        assert win.slider_res.value() == 30
        assert win.combo_base_style.currentIndex() == 1

    def test_refresh_ui_values_default_base_style(self, win):
        _vdw._vdw_settings["base_style"] = "default"
        win.refresh_ui_values()
        assert win.combo_base_style.currentIndex() == 0

    def test_update_view_calls_context_refresh(self, win):
        win.context.refresh_3d_view.reset_mock()
        win.update_view()
        win.context.refresh_3d_view.assert_called_once()

    def test_registers_window_with_context(self, qapp, tmp_path, monkeypatch):
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "vdw.json"))
        ctx = _vdw_context()
        w = _vdw.VDWConfigWindow(context=ctx)
        try:
            ctx.register_window.assert_called_once_with("main_panel", w)
        finally:
            w.destroy()


class TestOpenSettingsGui:
    def test_creates_and_shows_new_window(self, qapp, tmp_path, monkeypatch):
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "vdw.json"))
        ctx = _vdw_context()
        ctx.get_window.return_value = None
        _vdw.open_settings(ctx)
        win = ctx.register_window.call_args[0][1]
        try:
            assert win.isVisible()
        finally:
            win.destroy()

    def test_reuses_existing_registered_window(self, qapp, tmp_path, monkeypatch):
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "vdw.json"))
        ctx = _vdw_context()
        w = _vdw.VDWConfigWindow(context=ctx)
        ctx.get_window.return_value = w
        try:
            _vdw.open_settings(ctx)
            assert w.isVisible()
        finally:
            w.destroy()
