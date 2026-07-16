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
