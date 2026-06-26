"""
Headless GUI tests for plugin dialog widgets.

These tests use real PyQt6 (QT_QPA_PLATFORM=offscreen) rather than mocking
all Qt classes.  Chemistry/scientific libraries (rdkit, numpy, …) are still
replaced with MagicMock via mock_chemistry_imports() so no installed chemistry
stack is required.

Only visible plugins (registry visible=true) are tested here.

Run locally:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/ -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PASTE_XYZ_PATH = PLUGINS_DIR / "Paste_XYZ" / "paste_xyz.py"
GAUSSIAN_NEO_PATH = (
    PLUGINS_DIR / "Gaussian_Input_Generator_Neo" / "gaussian_input_generator_neo.py"
)
MS_NEO_PATH = PLUGINS_DIR / "MS_Spectrum_Simulation_Neo" / "ms_spectrum_neo.py"


# ---------------------------------------------------------------------------
# Load each plugin once at collection time.
# Real Qt classes are used; chemistry deps are MagicMock.
# ---------------------------------------------------------------------------

with mock_chemistry_imports():
    _paste_xyz = load_plugin_for_gui(PASTE_XYZ_PATH)
    _gaussian = load_plugin_for_gui(GAUSSIAN_NEO_PATH)
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
# PasteXYZDialog  (visible plugin: "Paste XYZ")
# ===========================================================================


class TestPasteXYZDialog:
    """PasteXYZDialog — pure-Qt dialog; rdkit guarded by try/except."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _paste_xyz.PasteXYZDialog(parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Paste XYZ"

    def test_default_size_wide_enough(self, dlg):
        assert dlg.width() >= 400

    def test_get_data_initially_empty(self, dlg):
        assert dlg.get_data() == ""

    def test_get_data_returns_typed_text(self, dlg):
        dlg.text_edit.setPlainText("C 0.0 0.0 0.0\nH 1.0 0.0 0.0")
        result = dlg.get_data()
        assert "C 0.0 0.0 0.0" in result
        assert "H 1.0 0.0 0.0" in result

    def test_multiline_xyz_round_trip(self, dlg):
        xyz = "C 0.0 0.0 0.0\nH 1.089 0.0 0.0\nH -0.363 1.027 0.0"
        dlg.text_edit.setPlainText(xyz)
        assert dlg.get_data().strip() == xyz

    def test_load_button_label(self, dlg):
        assert dlg.load_btn.text() == "Load"

    def test_cancel_button_label(self, dlg):
        assert dlg.cancel_btn.text() == "Cancel"

    def test_text_edit_placeholder_mentions_xyz(self, dlg):
        ph = dlg.text_edit.placeholderText().lower()
        assert "xyz" in ph

    def test_clear_resets_get_data(self, dlg):
        dlg.text_edit.setPlainText("something")
        dlg.text_edit.clear()
        assert dlg.get_data() == ""


# ===========================================================================
# RouteBuilderDialog  (visible plugin: "Gaussian Input Generator Neo")
# ===========================================================================


class TestRouteBuilderDialog:
    """RouteBuilderDialog — tabbed QDialog with combo-driven route preview."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _gaussian.RouteBuilderDialog(parent=None, current_route="")
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Route Builder"

    def test_has_four_tabs(self, dlg):
        assert dlg.tabs.count() == 4

    def test_tab_method_basis_present(self, dlg):
        names = [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())]
        assert any("Method" in n or "Basis" in n for n in names)

    def test_tab_job_type_present(self, dlg):
        names = [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())]
        assert any("Job" in n for n in names)

    def test_tab_properties_present(self, dlg):
        names = [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())]
        assert any("Propert" in n for n in names)

    def test_preview_label_exists(self, dlg):
        assert hasattr(dlg, "preview_label")

    def test_preview_non_empty_after_init(self, dlg):
        assert dlg.preview_label.text() != ""

    def test_preview_contains_default_basis_set(self, dlg):
        assert "6-31G" in dlg.preview_label.text()

    def test_ok_button_label(self, dlg):
        assert dlg.btn_ok.text() == "Apply to Job"

    def test_cancel_button_label(self, dlg):
        assert dlg.btn_cancel.text() == "Cancel"

    def test_method_type_combo_has_dft(self, dlg):
        items = [dlg.method_type.itemText(i) for i in range(dlg.method_type.count())]
        assert any("DFT" in it for it in items)

    def test_dft_method_list_contains_b3lyp(self, dlg):
        dlg.method_type.setCurrentText("DFT")
        items = [dlg.method_name.itemText(i) for i in range(dlg.method_name.count())]
        assert "B3LYP" in items

    def test_mp2_method_list_contains_ccsd(self, dlg):
        dlg.method_type.setCurrentText("MP2")
        items = [dlg.method_name.itemText(i) for i in range(dlg.method_name.count())]
        assert "CCSD" in items

    def test_hf_method_list_contains_hf(self, dlg):
        dlg.method_type.setCurrentText("Hartree-Fock")
        items = [dlg.method_name.itemText(i) for i in range(dlg.method_name.count())]
        assert "HF" in items

    def test_switching_method_updates_preview(self, dlg):
        dlg.method_type.setCurrentText("Hartree-Fock")
        dlg.method_name.setCurrentText("HF")
        assert "HF" in dlg.preview_label.text()

    def test_basis_set_change_updates_preview(self, dlg):
        dlg.method_type.setCurrentText("DFT")
        dlg.basis_set.setCurrentText("def2TZVP")
        assert "def2TZVP" in dlg.preview_label.text()

    def test_job_opt_only_shows_opt_group_hides_freq(self, dlg):
        # "Optimization Only (Opt)" is index 1
        dlg.job_type.setCurrentIndex(1)
        # Use isHidden() — isVisible() requires the window to be .show()n first
        assert not dlg.opt_group.isHidden()
        assert dlg.freq_group.isHidden()

    def test_job_freq_only_shows_freq_group_hides_opt(self, dlg):
        # "Frequency Only (Freq)" is index 2
        dlg.job_type.setCurrentIndex(2)
        assert dlg.opt_group.isHidden()
        assert not dlg.freq_group.isHidden()

    def test_job_sp_hides_both_groups(self, dlg):
        # "Single Point Energy (SP)" is index 3
        dlg.job_type.setCurrentIndex(3)
        assert dlg.opt_group.isHidden()
        assert dlg.freq_group.isHidden()

    def test_opt_freq_job_shows_both_groups(self, dlg):
        # "Optimization + Freq (Opt Freq)" is index 0
        dlg.job_type.setCurrentIndex(0)
        assert not dlg.opt_group.isHidden()
        assert not dlg.freq_group.isHidden()


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
