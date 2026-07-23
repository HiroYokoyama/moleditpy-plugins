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
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from PyQt6.QtCore import QEvent, QPoint, QPointF, Qt
from PyQt6.QtGui import QMouseEvent, QWheelEvent

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MS_NEO_PATH = PLUGINS_DIR / "MS_Spectrum_Simulation_Neo" / "ms_spectrum_neo.py"

with mock_chemistry_imports():
    _ms_neo = load_plugin_for_gui(MS_NEO_PATH)


# ===========================================================================
# Fake rdkit.Chem — real isotope math (a handful of published masses/
# abundances) so the real (un-extracted) _calculate_peaks / check_update
# bodies run to completion instead of chasing MagicMock attribute chains.
# ===========================================================================

_ISOTOPE_DATA = {
    "H": (1, [(1, 1.0078250319, 0.999885), (2, 2.0141017781, 0.000115)], 1.008),
    "C": (6, [(12, 12.0000000000, 0.9893), (13, 13.0033548378, 0.0107)], 12.011),
    "N": (7, [(14, 14.0030740052, 0.99636), (15, 15.0001088984, 0.00364)], 14.007),
    "O": (
        8,
        [
            (16, 15.9949146221, 0.99757),
            (17, 16.9991315, 0.00038),
            (18, 17.9991604, 0.00205),
        ],
        15.999,
    ),
    "Cl": (17, [(35, 34.96885268, 0.7576), (37, 36.96590259, 0.2424)], 35.453),
    "Na": (11, [(23, 22.98976928, 1.0)], 22.98977),
}
_ATOMIC_NUM_TO_SYM = {v[0]: k for k, v in _ISOTOPE_DATA.items()}


class _FakePT:
    def GetAtomicNumber(self, sym):
        if sym not in _ISOTOPE_DATA:
            return 0
        return _ISOTOPE_DATA[sym][0]

    def GetAtomicWeight(self, sym):
        return _ISOTOPE_DATA[sym][2]

    def GetMostCommonIsotope(self, atomic_num):
        sym = _ATOMIC_NUM_TO_SYM[atomic_num]
        return max(_ISOTOPE_DATA[sym][1], key=lambda t: t[2])[0]

    def GetAbundanceForIsotope(self, atomic_num, mass_number):
        sym = _ATOMIC_NUM_TO_SYM[atomic_num]
        for m, _mass, abundance in _ISOTOPE_DATA[sym][1]:
            if m == mass_number:
                return abundance
        return 0.0

    def GetMassForIsotope(self, atomic_num, mass_number):
        sym = _ATOMIC_NUM_TO_SYM[atomic_num]
        for m, mass, _abundance in _ISOTOPE_DATA[sym][1]:
            if m == mass_number:
                return mass
        raise RuntimeError(f"unknown isotope {mass_number} for {sym}")


class _FakeRDMolDescriptors:
    @staticmethod
    def CalcMolFormula(mol):
        return getattr(mol, "formula", "CH4")


class _FakeChem:
    rdMolDescriptors = _FakeRDMolDescriptors()

    @staticmethod
    def GetPeriodicTable():
        return _FakePT()

    @staticmethod
    def AddHs(mol):
        return mol


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
        pytest.importorskip("numpy")
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


# ===========================================================================
# check_update — real molecule sync paths (state_manager, fallback, silenced
# exceptions)
# ===========================================================================


class TestMSCheckUpdateReal:
    def test_uses_state_manager_to_rdkit_mol_when_2d_checked(self, qapp, monkeypatch):
        # get_main_window() is used directly as the QDialog parent -> must be
        # a real QWidget, not a SimpleNamespace.
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        fakemol = SimpleNamespace(formula="C6H6")
        mw = _ms_neo.QWidget()
        mw.state_manager = SimpleNamespace(
            data=SimpleNamespace(to_rdkit_mol=lambda: fakemol)
        )
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = None
        d = _ms_neo.MSSpectrumDialog(context=ctx)
        assert d.formula_input.text() == "C6H6"
        d.destroy()
        mw.destroy()

    def test_falls_back_to_current_molecule_without_state_manager(self, qapp, monkeypatch):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        mw = _ms_neo.QWidget()  # no state_manager attr
        fakemol = SimpleNamespace(formula="H2O")
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = fakemol
        d = _ms_neo.MSSpectrumDialog(context=ctx)
        assert d.formula_input.text() == "H2O"
        d.destroy()
        mw.destroy()

    def test_no_new_mol_returns_without_update(self, qapp, monkeypatch):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        mw = _ms_neo.QWidget()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = None
        d = _ms_neo.MSSpectrumDialog(context=ctx)
        assert d.formula_input.text() == ""
        d.destroy()
        mw.destroy()

    def test_exception_in_formula_calc_is_silenced(self, qapp, monkeypatch):
        class _BadChem:
            class rdMolDescriptors:
                @staticmethod
                def CalcMolFormula(mol):
                    raise ValueError("boom")

            @staticmethod
            def AddHs(mol):
                return mol

        monkeypatch.setattr(_ms_neo, "Chem", _BadChem)
        fakemol = SimpleNamespace()
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        ctx.current_molecule = fakemol
        d = _ms_neo.MSSpectrumDialog(context=ctx)  # must not raise
        assert d is not None
        d.destroy()


# ===========================================================================
# _calculate_peaks — real isotope math via _FakeChem, driven through the
# bound dialog widgets (formula_input, charge_spin, adduct_combo)
# ===========================================================================


class TestMSCalculatePeaksReal:
    def test_full_recalc_water_positive(self, qapp, monkeypatch):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        d, ctx = _make_dlg()
        d.formula_input.setText("H2O")
        assert d.peaks
        assert "Monoisotopic" in d.lbl_ion.text()
        d.destroy()

    def test_deuterium_formula(self, qapp, monkeypatch):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        d, ctx = _make_dlg()
        d.formula_input.setText("D2O")
        assert d.peaks
        d.destroy()

    def test_unknown_element_hides_spectrum(self, qapp, monkeypatch):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        d, ctx = _make_dlg()
        d.formula_input.setText("Xx2")
        assert d.peaks == []
        d.destroy()

    def test_invalid_formula_syntax_hides_spectrum(self, qapp, monkeypatch):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        d, ctx = _make_dlg()
        d.formula_input.setText("C6?H6")
        assert d.peaks == []
        d.destroy()

    def test_charge_zero_hides_spectrum_but_keeps_mass(self, qapp, monkeypatch):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        d, ctx = _make_dlg()
        d.formula_input.setText("H2O")
        d.charge_spin.setValue(0)
        assert d.peaks == []
        assert "Neutral Avg Mass" in d.lbl_mw.text()
        d.destroy()

    def test_deprotonation_removing_only_atom_yields_empty(self, qapp, monkeypatch):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        d, ctx = _make_dlg()
        d.charge_spin.setValue(-1)
        d.adduct_combo.setCurrentIndex(1)  # H (Deprotonation)
        d.formula_input.setText("H")
        assert d.peaks == []
        d.destroy()

    def test_gaussian_broadening_via_ui_checkbox(self, qapp, monkeypatch):
        pytest.importorskip("numpy")
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        d, ctx = _make_dlg()
        d.formula_input.setText("Cl2")
        d.gauss_check.setChecked(True)
        assert d.plot_widget.draw_mode == "profile"
        d.destroy()


# ===========================================================================
# get_adduct_delta / parse_formula_str — remaining bound branches
# ===========================================================================


class TestMSAdductLogicExtra:
    def test_positive_potassium_and_acetonitrile(self, qapp):
        d, ctx = _make_dlg()
        assert d.get_adduct_delta(3, "Positive", 1) == {"K": 1}
        assert d.get_adduct_delta(5, "Positive", 1) == {"C": 2, "H": 4, "N": 1}
        d.destroy()

    def test_negative_m_and_acetate(self, qapp):
        d, ctx = _make_dlg()
        assert d.get_adduct_delta(0, "Negative", 1) == {}
        assert d.get_adduct_delta(4, "Negative", 1) == {"C": 2, "H": 3, "O": 2}
        d.destroy()

    def test_leading_digit_token(self, qapp):
        d, ctx = _make_dlg()
        assert d.parse_formula_str("2H") == {"H": 1}
        d.destroy()

    def test_unmatched_open_paren_still_merges(self, qapp):
        d, ctx = _make_dlg()
        assert d.parse_formula_str("(CH2") == {"C": 1, "H": 2}
        d.destroy()

    def test_adduct_index_out_of_range_resets_to_zero(self, qapp):
        # 6 items at charge +2 -> pick index 5, then drop to charge -1 (5
        # items, no index 5) -> update_adduct_options must clamp to 0.
        d, ctx = _make_dlg()
        d.charge_spin.setValue(2)
        d.adduct_combo.setCurrentIndex(5)
        d.charge_spin.setValue(-1)
        assert d.adduct_combo.currentIndex() == 0
        d.destroy()


# ===========================================================================
# Export buttons — real file I/O
# ===========================================================================


class TestMSExportImage:
    def test_export_image_success(self, qapp, tmp_path, monkeypatch):
        d, ctx = _make_dlg()
        out = tmp_path / "spectrum.png"
        monkeypatch.setattr(
            _ms_neo.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(out), "")),
        )
        info_calls = []
        monkeypatch.setattr(
            _ms_neo.QMessageBox,
            "information",
            staticmethod(lambda *a, **k: info_calls.append(a)),
        )
        d.export_image()
        assert out.exists()
        assert info_calls
        d.destroy()

    def test_export_image_cancel_is_noop(self, qapp, tmp_path, monkeypatch):
        d, ctx = _make_dlg()
        monkeypatch.setattr(
            _ms_neo.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        d.export_image()
        assert list(tmp_path.iterdir()) == []
        d.destroy()

    def test_export_image_save_failure_shows_critical(self, qapp, tmp_path, monkeypatch):
        d, ctx = _make_dlg()
        monkeypatch.setattr(
            _ms_neo.QFileDialog,
            "getSaveFileName",
            staticmethod(
                lambda *a, **k: (str(tmp_path / "no_such_dir" / "x.png"), "")
            ),
        )
        crit_calls = []
        monkeypatch.setattr(
            _ms_neo.QMessageBox,
            "critical",
            staticmethod(lambda *a, **k: crit_calls.append(a)),
        )
        d.export_image()
        assert crit_calls
        d.destroy()


class TestMSExportCsv:
    def test_export_csv_writes_file(self, qapp, tmp_path, monkeypatch):
        d, ctx = _make_dlg()
        d.plot_widget.peaks = [(100.0, 50.0), (101.0, 100.0)]
        out = tmp_path / "spectrum.csv"
        monkeypatch.setattr(
            _ms_neo.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(out), "")),
        )
        info_calls = []
        monkeypatch.setattr(
            _ms_neo.QMessageBox,
            "information",
            staticmethod(lambda *a, **k: info_calls.append(a)),
        )
        d.export_csv()
        assert out.exists()
        assert "m/z,Intensity" in out.read_text()
        assert info_calls
        d.destroy()


# ===========================================================================
# PDF report generation — real QPrinter/QTextDocument
# ===========================================================================


class TestMSCreateReport:
    def test_create_report_writes_pdf(self, qapp, tmp_path, monkeypatch):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        d, ctx = _make_dlg()
        d.formula_input.setText("H2O")
        out = tmp_path / "report.pdf"
        monkeypatch.setattr(
            _ms_neo.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(out), "")),
        )
        info_calls = []
        monkeypatch.setattr(
            _ms_neo.QMessageBox,
            "information",
            staticmethod(lambda *a, **k: info_calls.append(a)),
        )
        d.create_report()
        assert out.exists()
        assert info_calls
        d.destroy()

    def test_create_report_with_molecule_hits_scene_and_draw_fallback(
        self, qapp, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(_ms_neo, "Chem", _FakeChem)
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        ctx.current_molecule = SimpleNamespace(formula="C6H6")
        # ctx.scene is a MagicMock: scene.items() isn't iterable -> the
        # scene-capture try/except is exercised, then falls back to Draw.
        d = _ms_neo.MSSpectrumDialog(context=ctx)
        out = tmp_path / "report_mol.pdf"
        monkeypatch.setattr(
            _ms_neo.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(out), "")),
        )
        monkeypatch.setattr(
            _ms_neo.QMessageBox, "information", staticmethod(lambda *a, **k: None)
        )
        d.create_report()
        assert out.exists()
        d.destroy()

    def test_create_report_cancel_is_noop(self, qapp, tmp_path, monkeypatch):
        d, ctx = _make_dlg()
        monkeypatch.setattr(
            _ms_neo.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        d.create_report()
        assert list(tmp_path.iterdir()) == []
        d.destroy()

    def test_create_report_pdf_error_shows_critical(self, qapp, tmp_path, monkeypatch):
        d, ctx = _make_dlg()
        out = tmp_path / "report.pdf"
        monkeypatch.setattr(
            _ms_neo.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(out), "")),
        )
        monkeypatch.setattr(
            _ms_neo, "QPrinter", MagicMock(side_effect=RuntimeError("boom"))
        )
        crit_calls = []
        monkeypatch.setattr(
            _ms_neo.QMessageBox,
            "critical",
            staticmethod(lambda *a, **k: crit_calls.append(a)),
        )
        d.create_report()
        assert crit_calls
        d.destroy()


# ===========================================================================
# HistogramWidget — real QWheelEvent / QMouseEvent driving
# ===========================================================================


class TestMSHistogramRealEvents:
    def test_wheel_zoom_in_updates_view(self, qapp):
        h = _ms_neo.HistogramWidget([(100.0, 50.0), (110.0, 80.0)])
        h.resize(600, 400)
        h.reset_view()
        orig_range = h.view_max - h.view_min
        event = QWheelEvent(
            QPointF(300, 200),
            QPointF(300, 200),
            QPoint(0, 0),
            QPoint(0, 120),
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
            Qt.ScrollPhase.NoScrollPhase,
            False,
        )
        h.wheelEvent(event)
        assert (h.view_max - h.view_min) == pytest.approx(orig_range * 0.9, rel=1e-6)
        h.destroy()

    def test_wheel_noop_before_reset_view(self, qapp):
        h = _ms_neo.HistogramWidget([])
        event = QWheelEvent(
            QPointF(10, 10),
            QPointF(10, 10),
            QPoint(0, 0),
            QPoint(0, -120),
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
            Qt.ScrollPhase.NoScrollPhase,
            False,
        )
        h.wheelEvent(event)  # must not raise
        assert h.view_min is None
        h.destroy()

    def test_mouse_drag_pans_and_releases(self, qapp):
        h = _ms_neo.HistogramWidget([(100.0, 50.0)])
        h.resize(600, 400)
        h.reset_view()
        press = QMouseEvent(
            QEvent.Type.MouseButtonPress,
            QPointF(300, 200),
            QPointF(300, 200),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
        h.mousePressEvent(press)
        assert h.last_mouse_x == 300.0
        move = QMouseEvent(
            QEvent.Type.MouseMove,
            QPointF(350, 200),
            QPointF(350, 200),
            Qt.MouseButton.NoButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
        h.mouseMoveEvent(move)
        assert h.view_min < 95.0  # dragging right pans the view left
        release = QMouseEvent(
            QEvent.Type.MouseButtonRelease,
            QPointF(350, 200),
            QPointF(350, 200),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
        )
        h.mouseReleaseEvent(release)
        assert h.last_mouse_x is None
        h.destroy()

    def test_double_click_resets_view(self, qapp):
        h = _ms_neo.HistogramWidget([(100.0, 50.0)])
        h.resize(600, 400)
        h.reset_view()
        h.view_min, h.view_max = -999.0, 999.0
        dbl = QMouseEvent(
            QEvent.Type.MouseButtonDblClick,
            QPointF(300, 200),
            QPointF(300, 200),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
        h.mouseDoubleClickEvent(dbl)
        assert h.view_min != -999.0
        h.destroy()


# ===========================================================================
# draw_spectrum — stick mode, no-data message, info text
# ===========================================================================


class TestMSPaintEventModes:
    def test_stick_mode_with_labels_renders(self, qapp):
        h = _ms_neo.HistogramWidget([(100.1234, 100.0), (101.5, 50.0)])
        h.stick_peaks = h.peaks
        h.info_text = "C6H6\n[M]+"
        h.resize(600, 400)
        # view_min/view_max left None -> draw_spectrum must self-initialize.
        pixmap = h.grab()
        assert not pixmap.isNull()
        assert h.view_min is not None
        h.destroy()

    def test_no_data_message_renders(self, qapp):
        h = _ms_neo.HistogramWidget([])
        h.resize(600, 400)
        pixmap = h.grab()
        assert not pixmap.isNull()
        h.destroy()


# ===========================================================================
# initialize() / toggle_window — window lifecycle
# ===========================================================================


class TestMSInitializeToggleWindow:
    def test_toggle_window_creates_then_reuses(self, qapp):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        ctx.current_molecule = None
        created = {}
        ctx.get_window.side_effect = lambda key: created.get(key)
        ctx.register_window.side_effect = created.__setitem__

        _ms_neo.initialize(ctx)
        callback = ctx.add_analysis_tool.call_args[0][1]

        callback()  # no existing window -> creates + registers one
        assert created.get("main_panel") is not None
        win = created["main_panel"]

        # Enable sync without firing the stateChanged signal so the timer
        # stays inactive, forcing toggle_window's reuse branch to start it.
        win.sync_check.blockSignals(True)
        win.sync_check.setChecked(True)
        win.sync_check.blockSignals(False)
        assert not win.timer.isActive()

        callback()  # existing window -> show/raise/activate + start timer
        assert win.timer.isActive()
        win.timer.stop()
        win.destroy()
