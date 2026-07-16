"""
Headless GUI tests for the Compound Info Report plugin.

Covers: ReportDialog.

generate_initial_report() calls build_html() which calls re.sub() on a
MagicMock rdkit formula — that raises TypeError.  Patch it to a no-op so
setup_ui() runs (creating all widgets) without triggering rdkit calls.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

COMPOUND_PATH = PLUGINS_DIR / "Compound_Info_Report" / "compound_info_report.py"

with mock_chemistry_imports():
    _compound = load_plugin_for_gui(COMPOUND_PATH)


# ===========================================================================
# ReportDialog  (Compound Info Report)
# ===========================================================================


class TestReportDialog:
    """ReportDialog with mol=None and generate_initial_report patched."""

    @pytest.fixture
    def dlg(self, qapp):
        with patch.object(
            _compound.ReportDialog, "generate_initial_report", lambda self: None
        ):
            d = _compound.ReportDialog(mol=None, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Compound Info Report"

    def test_default_size_tall_enough(self, dlg):
        assert dlg.height() >= 400

    def test_pubchem_checkbox_initially_unchecked(self, dlg):
        assert not dlg.chk_pubchem.isChecked()

    def test_preview_is_text_browser(self, dlg):
        from PyQt6.QtWidgets import QTextBrowser
        assert isinstance(dlg.preview, QTextBrowser)

    def test_preview_has_external_links_enabled(self, dlg):
        assert dlg.preview.openExternalLinks()

    def test_print_button_exists(self, dlg):
        assert dlg.btn_print.text() == "Print..."

    def test_pdf_button_exists(self, dlg):
        assert dlg.btn_pdf.text() == "Save PDF..."

    def test_close_button_exists(self, dlg):
        assert dlg.btn_close.text() == "Close"

    def test_calculate_adducts_returns_list(self, dlg):
        adducts = dlg.calculate_adducts(180.0634)
        assert isinstance(adducts, list)
        assert len(adducts) > 0

    def test_adduct_mh_plus_correct(self, dlg):
        adducts = dlg.calculate_adducts(100.0)
        names = [a[0] for a in adducts]
        assert any("M+H" in n for n in names)
