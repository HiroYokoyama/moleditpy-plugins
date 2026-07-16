"""
Headless GUI tests for the Molecule Comparator plugin.

Covers: MoleculeComparator.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_molecule_comparator.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

COMPARATOR_PATH = PLUGINS_DIR / "Molecule_Comparator" / "molecule_comparator.py"

with mock_chemistry_imports():
    _comparator = load_plugin_for_gui(COMPARATOR_PATH)


def _ctx_no_mol() -> MagicMock:
    """Context with no main window and no active molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# MoleculeComparator  (visible plugin: "Molecule Comparator")
# ===========================================================================


class TestMoleculeComparator:
    """MoleculeComparator with no main window."""

    @pytest.fixture
    def comp(self, qapp):
        ctx = _ctx_no_mol()
        w = _comparator.MoleculeComparator(context=ctx)
        yield w
        w.destroy()

    def test_creates_without_error(self, comp):
        assert comp is not None

    def test_window_title(self, comp):
        assert comp.windowTitle() == "Molecule Comparator"

    def test_mol_list_initially_empty(self, comp):
        assert comp.mol_list.count() == 0

    def test_molecules_list_initially_empty(self, comp):
        assert comp.molecules == []

    def test_style_combo_default_sticks(self, comp):
        assert comp.combo_style.currentText() == "Sticks"

    def test_style_combo_has_four_options(self, comp):
        assert comp.combo_style.count() == 4

    def test_align_method_default(self, comp):
        assert comp.combo_align_method.currentText() == "Substructure (MCS)"

    def test_ignore_hydrogens_unchecked_by_default(self, comp):
        assert not comp.check_ignore_hs.isChecked()

    def test_wireframe_lighting_unchecked_by_default(self, comp):
        assert not comp.check_wireframe_lighting.isChecked()

    def test_add_current_button_exists(self, comp):
        assert comp.btn_add_current.text() == "Add Current"

    def test_align_button_exists(self, comp):
        assert comp.btn_align.text() == "Align & Calculate RMSD"
