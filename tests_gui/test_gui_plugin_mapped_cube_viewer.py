"""
Headless GUI tests for the Mapped Cube Viewer plugin.

Covers: MappedCubeSetupDialog.

Chemistry libs (pyvista, numpy, rdkit, …) are mocked; real PyQt6 is used.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MAPPED_PATH = PLUGINS_DIR / "Mapped_Cube_Viewer" / "mapped_cube_viewer.py"

with mock_chemistry_imports():
    _mapped = load_plugin_for_gui(MAPPED_PATH)


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
