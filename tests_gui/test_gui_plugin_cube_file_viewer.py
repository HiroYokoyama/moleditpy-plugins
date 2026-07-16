"""
Headless GUI tests for the Cube File Viewer plugin.

Covers: ChargeDialog (the "bond connectivity error" dialog).

Chemistry libs (pyvista, numpy, rdkit, …) are mocked; real PyQt6 is used.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CUBE_VIEWER_PATH = PLUGINS_DIR / "Cube_File_Viewer" / "cube_viewer.py"

with mock_chemistry_imports():
    _cube = load_plugin_for_gui(CUBE_VIEWER_PATH)


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
