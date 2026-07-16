"""
Headless GUI tests for the ORCA Freq Analyzer plugin.

Covers:
  - OrcaOutFreqAnalyzer — main analyser widget (QWidget subclass), safe to
    construct with no main-window / no loaded file.
  - SpectrumDialog (QDialog) — accepts empty freq/intensity lists; numpy is
    mocked so np.array([]) returns a MagicMock whose __len__ returns 0,
    letting recalc_curve() exit early without crashing.
  - OrcaParser — pure-Python parser object that needs no Qt.

Run:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/test_gui_plugin_orca_freq_analyzer.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ORCA_PATH = PLUGINS_DIR / "ORCA_Freq_Analyzer" / "orca_out_freq_analyzer.py"

with mock_chemistry_imports():
    _orca = load_plugin_for_gui(ORCA_PATH)


# ===========================================================================
# OrcaOutFreqAnalyzer  (visible plugin: "ORCA Freq Analyzer")
# ===========================================================================


class TestOrcaOutFreqAnalyzer:
    """Main analyser widget — safe to construct with main_window=None."""

    @pytest.fixture
    def widget(self, qapp):
        w = _orca.OrcaOutFreqAnalyzer(None)
        yield w
        w.timer.stop()
        w.destroy()

    def test_creates_without_error(self, widget):
        assert widget is not None

    def test_info_label_text(self, widget):
        assert widget.lbl_info.text() == "Drop .out file here or click Open"

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

    def test_amplitude_slider_range(self, widget):
        assert widget.slider_amp.minimum() == 1
        assert widget.slider_amp.maximum() == 20

    def test_speed_slider_default(self, widget):
        assert widget.slider_speed.value() == 20

    def test_speed_slider_range(self, widget):
        assert widget.slider_speed.minimum() == 1
        assert widget.slider_speed.maximum() == 60

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
# SpectrumDialog (ORCA, module-level)  — numpy mocked; empty lists are safe
# ===========================================================================


class TestOrcaSpectrumDialog:
    """Top-level SpectrumDialog from the ORCA plugin.

    np.array([]) returns a MagicMock; SpectrumPlotWidget.__init__ only stores
    the arrays without operating on them, so no crash occurs.
    """

    @pytest.fixture
    def dlg(self, qapp):
        d = _orca.SpectrumDialog([], [], parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_default_window_title(self, dlg):
        assert dlg.windowTitle() == "IR Spectrum"

    def test_custom_title(self, qapp):
        d = _orca.SpectrumDialog([], [], title="Raman Spectrum", parent=None)
        assert d.windowTitle() == "Raman Spectrum"
        d.destroy()

    def test_fwhm_spinbox_default(self, dlg):
        assert dlg.spin_fwhm.value() == 50

    def test_fwhm_spinbox_range(self, dlg):
        assert dlg.spin_fwhm.minimum() == 1
        assert dlg.spin_fwhm.maximum() == 500

    def test_min_wavenumber_default(self, dlg):
        assert dlg.spin_min.value() == 0

    def test_max_wavenumber_default(self, dlg):
        assert dlg.spin_max.value() == 4000

    def test_size_at_least_600_wide(self, dlg):
        assert dlg.width() >= 600


# ===========================================================================
# OrcaParser  — pure Python, no Qt required
# ===========================================================================


class TestOrcaParser:
    """OrcaParser parses .out files; tested in isolation without opening any file."""

    def test_creates_without_error(self):
        assert _orca.OrcaParser() is not None

    def test_initial_charge(self):
        assert _orca.OrcaParser().charge == 0

    def test_initial_multiplicity(self):
        assert _orca.OrcaParser().multiplicity == 1

    def test_initial_frequencies_empty(self):
        assert _orca.OrcaParser().frequencies == []

    def test_initial_atoms_empty(self):
        assert _orca.OrcaParser().atoms == []

    def test_initial_coords_empty(self):
        assert _orca.OrcaParser().coords == []

    def test_initial_vib_modes_empty(self):
        assert _orca.OrcaParser().vib_modes == []

    def test_initial_intensities_empty(self):
        assert _orca.OrcaParser().intensities == []
