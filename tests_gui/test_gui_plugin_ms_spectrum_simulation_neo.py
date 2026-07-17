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


def _make_dlg(qapp_unused=None):
    ctx = _ms_context()
    return _ms_neo.MSSpectrumDialog(context=ctx), ctx


# ===========================================================================
# Formula parsing (pure logic)
# ===========================================================================


class TestMSFormulaParsing:
    def test_simple_formula(self, qapp):
        d, ctx = _make_dlg()
        assert d.parse_formula_str("C6H6") == {"C": 6, "H": 6}
        d.destroy()

    def test_parenthesized_group_with_multiplier(self, qapp):
        d, ctx = _make_dlg()
        assert d.parse_formula_str("Ca(OH)2") == {"Ca": 1, "O": 2, "H": 2}
        d.destroy()

    def test_nested_parentheses(self, qapp):
        d, ctx = _make_dlg()
        assert d.parse_formula_str("((CH3)2N)2") == {"C": 4, "H": 12, "N": 2}
        d.destroy()

    def test_charge_symbols_ignored(self, qapp):
        d, ctx = _make_dlg()
        assert d.parse_formula_str("C6H5+") == {"C": 6, "H": 5}
        d.destroy()

    def test_whitespace_allowed(self, qapp):
        d, ctx = _make_dlg()
        assert d.parse_formula_str("C6 H6") == {"C": 6, "H": 6}
        d.destroy()

    def test_invalid_characters_return_none(self, qapp):
        d, ctx = _make_dlg()
        assert d.parse_formula_str("C6H6$") is None
        d.destroy()

    def test_two_letter_elements(self, qapp):
        d, ctx = _make_dlg()
        assert d.parse_formula_str("NaCl") == {"Na": 1, "Cl": 1}
        d.destroy()


# ===========================================================================
# Superscript / subscript / adduct deltas (pure logic)
# ===========================================================================


class TestMSAdductLogic:
    def test_superscript_mapping(self, qapp):
        d, ctx = _make_dlg()
        assert d.to_superscript("2+") == "²⁺"
        assert d.to_superscript("10-") == "¹⁰⁻"
        d.destroy()

    def test_subscript_mapping(self, qapp):
        d, ctx = _make_dlg()
        assert d.to_subscript("43") == "₄₃"
        d.destroy()

    def test_positive_deltas(self, qapp):
        d, ctx = _make_dlg()
        assert d.get_adduct_delta(0, "Positive", 1) == {}
        assert d.get_adduct_delta(1, "Positive", 2) == {"H": 2}
        assert d.get_adduct_delta(2, "Positive", 1) == {"Na": 1}
        assert d.get_adduct_delta(4, "Positive", 1) == {"N": 1, "H": 4}
        d.destroy()

    def test_negative_deltas(self, qapp):
        d, ctx = _make_dlg()
        assert d.get_adduct_delta(1, "Negative", 1) == {"H": -1}
        assert d.get_adduct_delta(2, "Negative", 1) == {"Cl": 1}
        assert d.get_adduct_delta(3, "Negative", 1) == {"C": 1, "H": 1, "O": 2}
        d.destroy()


# ===========================================================================
# Adduct combo reacts to charge
# ===========================================================================


class TestMSAdductOptions:
    def test_neutral_charge_offers_neutral_items(self, qapp):
        d, ctx = _make_dlg()
        d.charge_spin.setValue(0)
        assert [d.adduct_combo.itemText(i) for i in range(d.adduct_combo.count())] == [
            "M (Neutral) [M]",
            "None",
        ]
        d.destroy()

    def test_positive_charge_offers_six_adducts(self, qapp):
        d, ctx = _make_dlg()
        d.charge_spin.setValue(1)
        assert d.adduct_combo.count() == 6
        assert d.adduct_combo.itemText(0) == "M [M]⁺"
        assert d.adduct_combo.itemText(1) == "H (Proton) [M+H]⁺"
        d.destroy()

    def test_negative_charge_offers_five_adducts(self, qapp):
        d, ctx = _make_dlg()
        d.charge_spin.setValue(-1)
        assert d.adduct_combo.count() == 5
        assert d.adduct_combo.itemText(1) == "H (Deprotonation) [M-H]⁻"
        d.destroy()

    def test_multiple_charge_shows_count_and_superscript(self, qapp):
        d, ctx = _make_dlg()
        d.charge_spin.setValue(2)
        assert d.adduct_combo.itemText(0) == "M [M]²⁺"
        assert d.adduct_combo.itemText(1) == "H (Proton) [M+2H]²⁺"
        d.destroy()

    def test_selection_preserved_across_charge_change(self, qapp):
        d, ctx = _make_dlg()
        d.charge_spin.setValue(1)
        d.adduct_combo.setCurrentIndex(2)
        d.charge_spin.setValue(2)
        assert d.adduct_combo.currentIndex() == 2
        d.destroy()


# ===========================================================================
# Sync timer lifecycle
# ===========================================================================


class TestMSSyncTimer:
    def test_toggle_sync_checked_starts_timer(self, qapp):
        d, ctx = _make_dlg()
        d.check_update = MagicMock()
        d.toggle_sync(2)
        assert d.timer.isActive()
        d.check_update.assert_called_once()
        d.timer.stop()
        d.destroy()

    def test_toggle_sync_unchecked_stops_timer(self, qapp):
        d, ctx = _make_dlg()
        d.check_update = MagicMock()
        d.toggle_sync(2)
        d.toggle_sync(0)
        assert not d.timer.isActive()
        d.destroy()

    def test_close_event_stops_running_timer(self, qapp):
        d, ctx = _make_dlg()
        d.timer.start(500)
        d.close()
        assert not d.timer.isActive()
        d.destroy()


# ===========================================================================
# HistogramWidget view management
# ===========================================================================


class TestMSHistogramWidget:
    def test_reset_view_defaults_without_peaks(self, qapp):
        h = _ms_neo.HistogramWidget([])
        h.reset_view()
        assert (h.view_min, h.view_max) == (0.0, 100.0)
        h.destroy()

    def test_reset_view_pads_peak_range(self, qapp):
        h = _ms_neo.HistogramWidget([(100.0, 50.0), (200.0, 100.0)])
        h.reset_view()
        assert (h.view_min, h.view_max) == (95.0, 205.0)
        h.destroy()

    def test_reset_view_single_mass_gets_symmetric_window(self, qapp):
        h = _ms_neo.HistogramWidget([(120.0, 100.0)])
        h.reset_view()
        assert (h.view_min, h.view_max) == (115.0, 125.0)
        h.destroy()

    def test_reset_view_prefers_stick_peaks(self, qapp):
        h = _ms_neo.HistogramWidget([(500.0, 1.0)])
        h.stick_peaks = [(100.0, 50.0)]
        h.reset_view()
        assert (h.view_min, h.view_max) == (95.0, 105.0)
        h.destroy()

    def test_default_draw_mode_is_stick(self, qapp):
        h = _ms_neo.HistogramWidget([])
        assert h.draw_mode == "stick"
        h.destroy()


# ===========================================================================
# CSV export cancel path
# ===========================================================================


class TestMSBroadeningZoomPreserved:
    def test_gauss_toggle_keeps_zoom(self, qapp):
        """Regression (v2026.07.17): gauss_check was double-connected and the
        direct connect passed the Qt state as reset= -> zoom wiped."""
        d, ctx = _make_dlg()
        d.plot_widget.view_min, d.plot_widget.view_max = 10.0, 20.0
        d.gauss_check.setChecked(True)
        assert (d.plot_widget.view_min, d.plot_widget.view_max) == (10.0, 20.0)
        d.destroy()

    def test_width_change_keeps_zoom(self, qapp):
        d, ctx = _make_dlg()
        d.plot_widget.view_min, d.plot_widget.view_max = 10.0, 20.0
        d.width_spin.setValue(0.5)
        assert (d.plot_widget.view_min, d.plot_widget.view_max) == (10.0, 20.0)
        d.destroy()


class TestMSProfileLocalMaxima:
    def test_single_apex_found(self, qapp):
        h = _ms_neo.HistogramWidget(
            [(100.0, 0.0), (100.1, 50.0), (100.2, 100.0), (100.3, 50.0), (100.4, 0.0)]
        )
        assert h._profile_local_maxima() == [(100.2, 100.0)]
        h.destroy()

    def test_two_separate_apexes_found(self, qapp):
        h = _ms_neo.HistogramWidget(
            [
                (100.0, 0.0),
                (100.1, 100.0),
                (100.2, 5.0),
                (100.3, 60.0),
                (100.4, 0.0),
            ]
        )
        assert h._profile_local_maxima() == [(100.1, 100.0), (100.3, 60.0)]
        h.destroy()

    def test_merged_curve_gives_single_apex(self, qapp):
        # Two overlapping components already summed into one hump
        h = _ms_neo.HistogramWidget(
            [(100.0, 10.0), (100.1, 80.0), (100.15, 100.0), (100.2, 80.0), (100.3, 10.0)]
        )
        assert h._profile_local_maxima() == [(100.15, 100.0)]
        h.destroy()

    def test_threshold_filters_noise(self, qapp):
        h = _ms_neo.HistogramWidget(
            [(100.0, 0.0), (100.1, 0.01), (100.2, 0.0), (100.3, 50.0), (100.4, 0.0)]
        )
        assert h._profile_local_maxima() == [(100.3, 50.0)]
        h.destroy()

    def test_flat_plateau_labelled_once(self, qapp):
        h = _ms_neo.HistogramWidget(
            [(100.0, 0.0), (100.1, 100.0), (100.2, 100.0), (100.3, 0.0)]
        )
        assert h._profile_local_maxima() == [(100.2, 100.0)]
        h.destroy()

    def test_empty_profile_has_no_maxima(self, qapp):
        h = _ms_neo.HistogramWidget([])
        assert h._profile_local_maxima() == []
        h.destroy()


class TestMSProfileRendering:
    def test_profile_mode_with_stick_labels_renders(self, qapp):
        # Executes the broadened-peak-top label path end to end.
        h = _ms_neo.HistogramWidget(
            [(100.0 + 0.01 * i, min(100.0, 10.0 * i)) for i in range(30)]
        )
        h.draw_mode = "profile"
        h.stick_peaks = [(100.1, 80.0), (100.2, 100.0), (100.25, 0.01)]
        h.resize(600, 400)
        h.reset_view()
        pixmap = h.grab()
        assert not pixmap.isNull()
        h.destroy()


class TestMSScaledMargins:
    def test_margins_scale_with_width(self, qapp):
        """Regression: wheel/pan used fixed 60/40 margins while the painter
        scales them, so the zoom anchor drifted on resized windows."""
        h = _ms_neo.HistogramWidget([])
        h.resize(1200, 400)
        assert h._scaled_x_margins() == (120, 80)
        h.destroy()


class TestMSExport:
    def test_export_csv_cancel_is_noop(self, qapp, tmp_path, monkeypatch):
        d, ctx = _make_dlg()
        monkeypatch.setattr(
            _ms_neo.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        d.export_csv()
        assert list(tmp_path.iterdir()) == []
        d.destroy()

class TestMSGaussianBroadening:
    def test_zero_sigma_does_not_crash(self, qapp):
        """Regression for zero division error when sigma is extremely small or 0"""
        d, ctx = _make_dlg()
        peaks = [(100.0, 50.0), (101.0, 100.0)]
        # This should not raise ZeroDivisionError
        result = d.apply_gaussian_broadening(peaks, 0.0)
        assert result is not None
        d.destroy()
        
    def test_spinbox_decimals(self, qapp):
        d, ctx = _make_dlg()
        # The spinbox should allow 3 decimals so 0.001 doesn't get rounded to 0.00
        assert d.width_spin.decimals() == 3
        d.destroy()
