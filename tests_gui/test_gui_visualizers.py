"""
Headless GUI tests for visualizer plugin dialogs.

Covers:
  - ChargeDialog          (Cube File Viewer)
  - VDWConfigWindow       (VDW Radii Overlay)
  - MappedCubeSetupDialog (Mapped Cube Viewer)
  - ChargeDialog          (Cube File Viewer Advanced)
  - FlexibleDoubleSpinBox (Cube File Viewer Advanced)

Chemistry libs (pyvista, numpy, rdkit, …) are mocked; real PyQt6 is used.
Run locally:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/test_gui_visualizers.py -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CUBE_VIEWER_PATH = PLUGINS_DIR / "Cube_File_Viewer" / "cube_viewer.py"
VDW_PATH = PLUGINS_DIR / "VDW_Radii_Overlay" / "vdw_radii_overlay.py"
MAPPED_PATH = PLUGINS_DIR / "Mapped_Cube_Viewer" / "mapped_cube_viewer.py"
CUBE_ADV_PATH = PLUGINS_DIR / "Cube_File_Viewer_Advanced" / "cube_viewer_advanced.py"

with mock_chemistry_imports():
    _cube = load_plugin_for_gui(CUBE_VIEWER_PATH)
    _vdw = load_plugin_for_gui(VDW_PATH)
    _mapped = load_plugin_for_gui(MAPPED_PATH)
    _cube_adv = load_plugin_for_gui(CUBE_ADV_PATH)


def _vdw_context() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


# ===========================================================================
# ChargeDialog  (Cube File Viewer)
# ===========================================================================


class TestCubeViewerChargeDialog:
    """ChargeDialog: the 'bond connectivity error' dialog in Cube File Viewer."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _cube.ChargeDialog(parent=None, current_charge=0)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Bond Connectivity Error"

    def test_default_result_action_is_cancel(self, dlg):
        assert dlg.result_action == "cancel"

    def test_charge_spinbox_range(self, dlg):
        assert dlg.spin.minimum() == -20
        assert dlg.spin.maximum() == 20

    def test_default_charge_value(self, dlg):
        assert dlg.spin.value() == 0

    def test_nonzero_charge_reflected_in_spin(self, qapp):
        d = _cube.ChargeDialog(parent=None, current_charge=2)
        assert d.spin.value() == 2
        d.destroy()

    def test_on_retry_sets_result_action(self, dlg):
        dlg.spin.setValue(1)
        dlg.on_retry()
        assert dlg.result_action == "retry"
        assert dlg.charge == 1

    def test_on_skip_sets_result_action(self, dlg):
        dlg.on_skip()
        assert dlg.result_action == "skip"


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


# ===========================================================================
# MappedCubeSetupDialog  (Mapped Cube Viewer)
# ===========================================================================


class TestMappedCubeSetupDialog:
    """MappedCubeSetupDialog: file-picker dialog for surface + property cubes."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _mapped.MappedCubeSetupDialog(parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Mapped Cube Setup"

    def test_surface_file_initially_none(self, dlg):
        assert dlg.surface_file is None

    def test_property_file_initially_none(self, dlg):
        assert dlg.property_file is None

    def test_surface_line_edit_empty(self, dlg):
        assert dlg.le_surf.text() == ""

    def test_property_line_edit_empty(self, dlg):
        assert dlg.le_prop.text() == ""

    def test_line_edits_accept_text(self, dlg):
        dlg.le_surf.setText("/some/surf.cube")
        dlg.le_prop.setText("/some/prop.cube")
        assert dlg.le_surf.text() == "/some/surf.cube"
        assert dlg.le_prop.text() == "/some/prop.cube"


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
