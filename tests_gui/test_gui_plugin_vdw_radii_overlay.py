"""
Headless GUI tests for the VDW Radii Overlay plugin.

Covers: VDWConfigWindow.

Chemistry libs (pyvista, numpy, rdkit, …) are mocked; real PyQt6 is used.
"""

from __future__ import annotations

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
    def win(self, qapp):
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
