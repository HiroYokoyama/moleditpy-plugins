"""
Headless GUI tests for the Symmetry Analyzer plugin.

Covers: SymmetryAnalysisPlugin.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_symmetry_analyzer.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

SYMMETRY_PATH = PLUGINS_DIR / "Symmetry_Analyzer" / "symmetry_analyzer.py"

with mock_chemistry_imports():
    _symmetry = load_plugin_for_gui(SYMMETRY_PATH)


def _ctx_no_mol() -> MagicMock:
    """Context with no main window and no active molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# SymmetryAnalysisPlugin  (visible plugin: "Symmetry Analyzer")
# ===========================================================================


class TestSymmetryAnalysisPlugin:
    """SymmetryAnalysisPlugin — pymatgen is mocked so HAS_PYMATGEN=True."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        d = _symmetry.SymmetryAnalysisPlugin(context=ctx)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert "Symmetry" in dlg.windowTitle() or dlg.windowTitle() == ""

    def test_groups_list_is_empty_initially(self, dlg):
        assert dlg.groups_list.count() == 0

    def test_ops_list_is_empty_initially(self, dlg):
        assert dlg.ops_list.count() == 0

    def test_point_group_label_default(self, dlg):
        assert dlg.selected_group_label.text() == "Point Group: -"

    def test_op_details_is_readonly(self, dlg):
        assert dlg.op_details.isReadOnly()

    def test_max_tol_spin_default(self, dlg):
        assert dlg.max_tol_spin.value() == pytest.approx(1.0)

    def test_symmetrize_button_initially_disabled(self, dlg):
        assert not dlg.sym_btn.isEnabled()

    def test_analyze_button_exists(self, dlg):
        assert dlg.calc_btn is not None
