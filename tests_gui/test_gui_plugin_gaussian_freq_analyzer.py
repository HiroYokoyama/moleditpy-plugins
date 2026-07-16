"""
Headless GUI tests for the Gaussian Freq Analyzer plugin.

Covers:
  - GaussianFCHKFreqAnalyzer — main analyser widget (QWidget subclass), safe
    to construct when context.get_main_window()=None.
  - SpectrumDialog (QDialog) — on_range_changed() is called at end of
    __init__; it calls recalc_curve() which guards on len(freqs)==0
    (MagicMock.__len__ returns 0), so no crash.
  - FCHKParser — pure-Python parser object that needs no Qt.

Run:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/test_gui_plugin_gaussian_freq_analyzer.py -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

GAUSSIAN_PATH = PLUGINS_DIR / "Gaussian_Freq_Analyzer" / "gaussian_fchk_freq_analyzer.py"

with mock_chemistry_imports():
    _gaussian = load_plugin_for_gui(GAUSSIAN_PATH)


def _gaussian_context() -> MagicMock:
    """Minimal stub context — get_main_window() returns None."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


# ===========================================================================
# GaussianFCHKFreqAnalyzer  (visible plugin: "Gaussian Freq Analyzer")
# ===========================================================================


class TestGaussianFCHKFreqAnalyzer:
    """Main analyser widget — safe to construct when context.get_main_window()=None."""

    @pytest.fixture
    def widget(self, qapp):
        ctx = _gaussian_context()
        w = _gaussian.GaussianFCHKFreqAnalyzer(ctx)
        yield w
        w.timer.stop()
        w.destroy()

    def test_creates_without_error(self, widget):
        assert widget is not None

    def test_info_label_text(self, widget):
        assert widget.lbl_info.text() == "Drop .fchk file here or click Open"

    def test_freq_list_column_count(self, widget):
        assert widget.list_freq.columnCount() == 3

    def test_freq_list_header_labels(self, widget):
        header = widget.list_freq.headerItem()
        assert header.text(0) == "No."
        assert "Frequency" in header.text(1)
        assert "Intensity" in header.text(2)

    def test_freq_list_initially_empty(self, widget):
        assert widget.list_freq.topLevelItemCount() == 0

    def test_scaling_factor_default(self, widget):
        assert widget.spin_sf.value() == 1.0

    def test_scaling_factor_range(self, widget):
        assert widget.spin_sf.minimum() == 0.0
        assert widget.spin_sf.maximum() == 2.0

    def test_show_vectors_unchecked(self, widget):
        assert not widget.chk_vectors.isChecked()

    def test_vector_scale_default(self, widget):
        assert widget.spin_vec_scale.value() == 2.0

    def test_amplitude_slider_default(self, widget):
        assert widget.slider_amp.value() == 5

    def test_speed_slider_default(self, widget):
        assert widget.slider_speed.value() == 20

    def test_play_button_text(self, widget):
        assert widget.btn_play.text() == "Play"

    def test_stop_button_text(self, widget):
        assert widget.btn_stop.text() == "Stop"

    def test_timer_not_active_initially(self, widget):
        assert not widget.timer.isActive()

    def test_parser_initially_none(self, widget):
        assert widget.parser is None

    def test_meta_label_initially_empty(self, widget):
        assert widget.lbl_meta.text() == ""

    def test_is_not_playing_initially(self, widget):
        assert not widget.is_playing

    def test_accepts_drops(self, widget):
        assert widget.acceptDrops()


# ===========================================================================
# SpectrumDialog (Gaussian, module-level)
# ===========================================================================


class TestGaussianSpectrumDialog:
    """Top-level SpectrumDialog from the Gaussian plugin.

    on_range_changed() is called at end of __init__; it calls recalc_curve()
    which guards on len(freqs)==0 (MagicMock.__len__ returns 0), so no crash.
    """

    @pytest.fixture
    def dlg(self, qapp):
        d = _gaussian.SpectrumDialog([], [], parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_default_window_title(self, dlg):
        assert dlg.windowTitle() == "IR Spectrum"

    def test_custom_title(self, qapp):
        d = _gaussian.SpectrumDialog([], [], title="FCHK Raman")
        assert d.windowTitle() == "FCHK Raman"
        d.destroy()

    def test_fwhm_spinbox_default(self, dlg):
        assert dlg.spin_fwhm.value() == 50

    def test_min_wavenumber_default(self, dlg):
        assert dlg.spin_min.value() == 0

    def test_max_wavenumber_default(self, dlg):
        assert dlg.spin_max.value() == 4000

    def test_size_at_least_600_wide(self, dlg):
        assert dlg.width() >= 600


# ===========================================================================
# FCHKParser  — pure Python, no Qt required
# ===========================================================================


class TestFCHKParser:
    """FCHKParser parses Gaussian FCHK files; tested in isolation."""

    def test_creates_without_error(self):
        assert _gaussian.FCHKParser() is not None

    def test_initial_charge(self):
        assert _gaussian.FCHKParser().charge == 0

    def test_initial_multiplicity(self):
        assert _gaussian.FCHKParser().multiplicity == 1

    def test_initial_frequencies_empty(self):
        assert _gaussian.FCHKParser().frequencies == []

    def test_initial_atoms_empty(self):
        assert _gaussian.FCHKParser().atoms == []

    def test_initial_coords_empty(self):
        assert _gaussian.FCHKParser().coords == []

    def test_initial_vib_modes_empty(self):
        assert _gaussian.FCHKParser().vib_modes == []
