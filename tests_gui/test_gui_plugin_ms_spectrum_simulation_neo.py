"""
Headless GUI tests for the MS Spectrum Simulation Neo plugin.

Covers: MSSpectrumDialog.

These tests use real PyQt6 (QT_QPA_PLATFORM=offscreen) rather than mocking
all Qt classes.  Chemistry/scientific libraries (rdkit, numpy, …) are still
replaced with MagicMock via mock_chemistry_imports() so no installed chemistry
stack is required.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MS_NEO_PATH = PLUGINS_DIR / "MS_Spectrum_Simulation_Neo" / "ms_spectrum_neo.py"

with mock_chemistry_imports():
    _ms_neo = load_plugin_for_gui(MS_NEO_PATH)


def _ms_context() -> MagicMock:
    """Minimal stub context for MSSpectrumDialog (Neo).

    get_main_window() returns None → dialog has no parent, sync disabled.
    current_molecule = None  → formula_input stays blank, no Chem calls.
    """
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# MSSpectrumDialog  (visible plugin: "MS Spectrum Simulation Neo")
# ===========================================================================


class TestMSSpectrumDialogNeo:
    """MSSpectrumDialog (Neo) with mol=None context — exercises no-molecule UI."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ms_context()
        d = _ms_neo.MSSpectrumDialog(context=ctx)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "MS Spectrum Simulation Neo"

    def test_default_size_tall_enough(self, dlg):
        assert dlg.height() >= 400

    def test_formula_input_initially_empty(self, dlg):
        assert dlg.formula_input.text() == ""

    def test_export_button_exists(self, dlg):
        assert dlg.export_btn.text() == "Export to Image"

    def test_csv_button_exists(self, dlg):
        assert dlg.btn_export_csv.text() == "Export CSV"

    def test_sync_check_disabled_when_no_main_window(self, dlg):
        # With get_main_window()=None, sync_check must be disabled
        assert not dlg.sync_check.isEnabled()

    def test_charge_spinbox_default_positive(self, dlg):
        assert dlg.charge_spin.value() == 1

    def test_adduct_combo_populated(self, dlg):
        assert dlg.adduct_combo.count() > 0

    def test_formula_input_accepts_text(self, dlg):
        dlg.formula_input.setText("C6H12O6")
        assert dlg.formula_input.text() == "C6H12O6"
