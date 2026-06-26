"""
Headless GUI tests for:
  - Compound Info Report  → ReportDialog
  - Settings Saver        → SettingsSaverDialog
  - Gaussian FCHK Loader  → FCHKLoaderDialog
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

COMPOUND_PATH = PLUGINS_DIR / "Compound_Info_Report" / "compound_info_report.py"
SETTINGS_PATH = PLUGINS_DIR / "Settings_Saver" / "settings_saver.py"
FCHK_PATH = PLUGINS_DIR / "Gaussian_FCHK_Loader" / "gaussian_fchk_loader.py"

with mock_chemistry_imports():
    _compound = load_plugin_for_gui(COMPOUND_PATH)
    _settings = load_plugin_for_gui(SETTINGS_PATH)
    _fchk = load_plugin_for_gui(FCHK_PATH)


# ===========================================================================
# ReportDialog  (Compound Info Report)
# ===========================================================================
# generate_initial_report() calls build_html() which calls re.sub() on a
# MagicMock rdkit formula — that raises TypeError.  Patch it to a no-op so
# setup_ui() runs (creating all widgets) without triggering rdkit calls.


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


# ===========================================================================
# SettingsSaverDialog  (Settings Saver)
# ===========================================================================


def _settings_context() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


class TestSettingsSaverDialog:
    """SettingsSaverDialog with MagicMock context and parent=None."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _settings_context()
        d = _settings.SettingsSaverDialog(context=ctx, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Settings Saver Manager"

    def test_preset_list_initially_empty(self, dlg):
        assert dlg.preset_list.count() == 0

    def test_load_button_exists(self, dlg):
        assert dlg.btn_load.text() == "Load Preset"

    def test_save_button_exists(self, dlg):
        assert dlg.btn_save.text() == "Save New..."

    def test_delete_button_exists(self, dlg):
        assert dlg.btn_delete.text() == "Delete"

    def test_embed_checkbox_unchecked_by_default(self, dlg):
        assert not dlg.chk_embed.isChecked()

    def test_global_default_button_disabled_when_embed_off(self, dlg):
        assert not dlg.btn_set_global.isEnabled()

    def test_export_button_has_menu(self, dlg):
        assert dlg.btn_export.menu() is not None

    def test_close_button_exists(self, dlg):
        assert dlg.btn_close.text() == "Close"


# ===========================================================================
# FCHKLoaderDialog  (Gaussian FCHK Loader)
# ===========================================================================


def _fchk_context() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = MagicMock()
    return ctx


class TestFCHKLoaderDialog:
    """FCHKLoaderDialog with a fake fchk path and MagicMock context."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _fchk_context()
        d = _fchk.FCHKLoaderDialog(
            parent=None, context=ctx, fchk_path="/tmp/test_molecule.fchk"
        )
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Select FCHK Analysis Mode"

    def test_freq_button_exists(self, dlg):
        assert "Frequency" in dlg.btn_freq.text()

    def test_mo_button_exists(self, dlg):
        assert "MO" in dlg.btn_mo.text() or "Orbital" in dlg.btn_mo.text()

    def test_freq_button_disabled_when_plugin_absent(self, dlg):
        # In test environment the freq analyzer plugin is not found via
        # find_file_recursive, so the button should be disabled.
        if dlg.freq_analyzer_path is None:
            assert not dlg.btn_freq.isEnabled()

    def test_fchk_path_stored(self, dlg):
        assert dlg.fchk_path == "/tmp/test_molecule.fchk"

    def test_basename_shown_in_label(self, dlg):
        # The label text is set during init_ui from os.path.basename(fchk_path)
        # Verify we can find text containing the filename in the dialog's children
        from PyQt6.QtWidgets import QLabel
        labels = dlg.findChildren(QLabel)
        label_texts = " ".join(lbl.text() for lbl in labels)
        assert "test_molecule.fchk" in label_texts
