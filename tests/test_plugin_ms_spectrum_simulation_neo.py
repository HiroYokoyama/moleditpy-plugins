"""
Tests for the MS Spectrum Simulation Neo plugin: parse_formula_str,
to_subscript, to_superscript, get_adduct_delta.

MSSpectrumDialog inherits a Qt base, so its methods are extracted from the
plugin source via AST and compiled standalone (self is unused in all of them).
"""

from __future__ import annotations

import ast
import textwrap
from pathlib import Path

import pytest

from conftest import extract_function, load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
MS_PATH = PLUGINS_DIR / "MS_Spectrum_Simulation_Neo" / "ms_spectrum_neo.py"

with mock_optional_imports():
    _ms = load_plugin(MS_PATH)

# MSSpectrumDialog methods — self is unused in all of them
_parse_formula_str = extract_function(MS_PATH, "MSSpectrumDialog", "parse_formula_str")
_to_subscript = extract_function(MS_PATH, "MSSpectrumDialog", "to_subscript")
_to_superscript = extract_function(MS_PATH, "MSSpectrumDialog", "to_superscript")
_get_adduct_delta = extract_function(MS_PATH, "MSSpectrumDialog", "get_adduct_delta")


class TestParseFormulaStr:
    def test_water(self):
        assert _parse_formula_str(None, "H2O") == {"H": 2, "O": 1}

    def test_glucose(self):
        assert _parse_formula_str(None, "C6H12O6") == {"C": 6, "H": 12, "O": 6}

    def test_single_atom_no_count(self):
        assert _parse_formula_str(None, "C") == {"C": 1}

    def test_parentheses_multiplier(self):
        assert _parse_formula_str(None, "(CH2)3") == {"C": 3, "H": 6}

    def test_multichar_element(self):
        assert _parse_formula_str(None, "NaCl") == {"Na": 1, "Cl": 1}

    def test_caffeine(self):
        assert _parse_formula_str(None, "C8H10N4O2") == {"C": 8, "H": 10, "N": 4, "O": 2}

    def test_empty_returns_falsy(self):
        result = _parse_formula_str(None, "")
        assert not result

    def test_invalid_char_returns_none(self):
        assert _parse_formula_str(None, "C6?H6") is None

    def test_charge_symbol_stripped(self):
        # '+' is in the tokenizer — should not produce None
        result = _parse_formula_str(None, "C6H5+")
        assert result is not None
        assert result.get("C") == 6


class TestSuperSubscript:
    def test_subscript_digits(self):
        assert _to_subscript(None, "123") == "₁₂₃"

    def test_superscript_positive(self):
        assert _to_superscript(None, "2+") == "²⁺"

    def test_superscript_minus(self):
        assert _to_superscript(None, "1-") == "¹⁻"

    def test_unknown_chars_pass_through(self):
        assert _to_subscript(None, "abc") == "abc"

    def test_empty_string(self):
        assert _to_subscript(None, "") == ""


class TestGetAdductDelta:
    # Positive mode
    def test_positive_m_no_delta(self):
        assert _get_adduct_delta(None, 0, "Positive", 1) == {}

    def test_positive_proton(self):
        assert _get_adduct_delta(None, 1, "Positive", 1) == {"H": 1}

    def test_positive_sodium(self):
        assert _get_adduct_delta(None, 2, "Positive", 1) == {"Na": 1}

    def test_positive_potassium(self):
        assert _get_adduct_delta(None, 3, "Positive", 1) == {"K": 1}

    def test_positive_double_charge_proton(self):
        assert _get_adduct_delta(None, 1, "Positive", 2) == {"H": 2}

    # Negative mode
    def test_negative_m_no_delta(self):
        assert _get_adduct_delta(None, 0, "Negative", 1) == {}

    def test_negative_deprotonation(self):
        assert _get_adduct_delta(None, 1, "Negative", 1) == {"H": -1}

    def test_negative_chloride(self):
        assert _get_adduct_delta(None, 2, "Negative", 1) == {"Cl": 1}

    def test_negative_formate(self):
        assert _get_adduct_delta(None, 3, "Negative", 1) == {"C": 1, "H": 1, "O": 2}




# ---------------------------------------------------------------------------
# Gaussian broadening (requires real numpy)
# ---------------------------------------------------------------------------


def _extract_method_source(path: Path, class_name: str, method_name: str) -> str:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == method_name:
                    return textwrap.dedent(ast.get_source_segment(source, item))
    raise AssertionError(f"{class_name}.{method_name} not found in {path.name}")


def _make_function(src: str, namespace: dict):
    exec(src, namespace)  # noqa: S102 - test-only extraction
    name = src.split("def ", 1)[1].split("(", 1)[0]
    return namespace[name]


class TestGaussianBroadening:
    def _fn(self):
        pytest.importorskip("numpy")
        src = _extract_method_source(MS_PATH, "MSSpectrumDialog", "apply_gaussian_broadening")
        return _make_function(src, {})

    def test_empty_peaks_returns_empty(self):
        assert self._fn()(None, [], 0.05) == []

    def test_profile_normalized_to_100(self):
        result = self._fn()(None, [(100.0, 55.0)], 0.05)
        max_y = max(y for _, y in result)
        assert max_y == pytest.approx(100.0, abs=1e-6)

    def test_peak_maximum_near_center(self):
        result = self._fn()(None, [(200.0, 10.0)], 0.05)
        x_at_max = max(result, key=lambda p: p[1])[0]
        assert x_at_max == pytest.approx(200.0, abs=0.01)

    def test_two_isotope_peaks_keep_relative_height(self):
        result = self._fn()(None, [(100.0, 100.0), (101.0, 50.0)], 0.02)
        near_100 = max(y for x, y in result if abs(x - 100.0) < 0.1)
        near_101 = max(y for x, y in result if abs(x - 101.0) < 0.1)
        assert near_100 == pytest.approx(100.0, abs=1e-6)
        assert near_101 == pytest.approx(50.0, rel=0.02)

    def test_apex_exactly_at_stick_mass(self):
        """v2026.07.17: the grid contains the exact stick masses, so an
        isolated peak's apex is exactly the theoretical m/z at 100%
        (previously the apex was the nearest arange point)."""
        result = self._fn()(None, [(123.4567, 100.0)], 0.03)
        apex_x, apex_y = max(result, key=lambda p: p[1])
        assert apex_x == 123.4567
        assert apex_y == pytest.approx(100.0, abs=1e-9)

    def test_every_stick_mass_is_a_grid_point(self):
        sticks = [(100.0123, 100.0), (101.0156, 32.0), (102.0189, 5.0)]
        result = self._fn()(None, sticks, 0.02)
        xs = {x for x, _ in result}
        for s_mass, _ in sticks:
            assert s_mass in xs

    def test_overlapping_sticks_merge_to_single_apex_between(self):
        # 0.5 sigma apart: broadened components merge into one hump whose
        # apex lies between the sticks (the "real" observed peak).
        sigma = 0.05
        result = self._fn()(None, [(100.00, 100.0), (100.025, 100.0)], sigma)
        apex_x, _ = max(result, key=lambda p: p[1])
        assert 100.00 < apex_x < 100.025


# ---------------------------------------------------------------------------
# _calculate_peaks — isotope-pattern convolution, formula/adduct integration
#
# Exercised with a lightweight fake rdkit periodic table carrying real
# published isotope masses/abundances for a handful of common elements, so
# the assertions check real chemistry (e.g. the classic Cl2 9:6:1-ish
# isotope envelope) rather than mocked numbers.
# ---------------------------------------------------------------------------

import logging as _logging

_ISOTOPE_DATA = {
    # symbol: (atomic_number, [(mass_number, exact_mass, abundance), ...], avg_weight)
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
    "S": (
        16,
        [
            (32, 31.9720710, 0.9499),
            (33, 32.97145876, 0.0075),
            (34, 33.9678669, 0.0425),
            (36, 35.96708076, 0.0001),
        ],
        32.065,
    ),
}
_ATOMIC_NUM_TO_SYM = {v[0]: k for k, v in _ISOTOPE_DATA.items()}
_ELECTRON_MASS = 0.00054858


class _FakePeriodicTable:
    def GetAtomicNumber(self, sym):
        if sym not in _ISOTOPE_DATA:
            return 0
        return _ISOTOPE_DATA[sym][0]

    def GetAtomicWeight(self, sym):
        return _ISOTOPE_DATA[sym][2]  # KeyError -> caller's except Exception

    def GetMostCommonIsotope(self, atomic_num):
        sym = _ATOMIC_NUM_TO_SYM[atomic_num]
        isotopes = _ISOTOPE_DATA[sym][1]
        return max(isotopes, key=lambda t: t[2])[0]

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


class _FakeChemModule:
    @staticmethod
    def GetPeriodicTable():
        return _FakePeriodicTable()


class _Val:
    """Stand-in Qt widget exposing whichever accessor _calculate_peaks calls."""

    def __init__(self, v):
        self.v = v

    def text(self):
        return self.v

    def value(self):
        return self.v

    def currentIndex(self):
        return self.v


class _FakeDialogForPeaks:
    def __init__(self, formula, charge, species_idx=0):
        self.formula_input = _Val(formula)
        self.charge_spin = _Val(charge)
        self.adduct_combo = _Val(species_idx)

    def parse_formula_str(self, formula):
        return _parse_formula_str(None, formula)

    def get_adduct_delta(self, species_idx, mode, charge):
        return _get_adduct_delta(None, species_idx, mode, charge)


_calculate_peaks = extract_function(
    MS_PATH,
    "MSSpectrumDialog",
    "_calculate_peaks",
    extra_globals={"Chem": _FakeChemModule, "logging": _logging},
)


class TestCalculatePeaksNeutralMass:
    def test_water_average_and_exact_mass_neutral(self):
        fake = _FakeDialogForPeaks("H2O", charge=0)
        peaks, avg_mw, exact_mw, ion_mz = _calculate_peaks(fake)
        assert peaks == []
        assert ion_mz == 0.0
        assert avg_mw == pytest.approx(1.008 * 2 + 15.999, abs=1e-6)
        assert exact_mw == pytest.approx(
            2 * 1.0078250319 + 15.9949146221, abs=1e-6
        )

    def test_dimethyl_ether_parentheses_average_mass(self):
        # (CH3)2O -> C2H6O; exercises parenthesis-multiplier parsing feeding
        # straight into the isotope/mass calculation.
        fake = _FakeDialogForPeaks("(CH3)2O", charge=0)
        peaks, avg_mw, exact_mw, ion_mz = _calculate_peaks(fake)
        assert peaks == []
        expected_avg = 2 * 12.011 + 6 * 1.008 + 15.999
        assert avg_mw == pytest.approx(expected_avg, abs=1e-6)

    def test_deuterium_special_cased_mass(self):
        fake = _FakeDialogForPeaks("D2O", charge=0)
        _peaks, avg_mw, exact_mw, _ion_mz = _calculate_peaks(fake)
        assert avg_mw == pytest.approx(2 * 2.0141017781 + 15.999, abs=1e-6)
        assert exact_mw == pytest.approx(
            2 * 2.0141017781 + 15.9949146221, abs=1e-6
        )

    def test_unrecognized_element_returns_empty(self):
        fake = _FakeDialogForPeaks("Xx2", charge=0)
        result = _calculate_peaks(fake)
        assert result == ([], 0.0, 0.0, 0.0)

    def test_empty_formula_returns_empty(self):
        fake = _FakeDialogForPeaks("", charge=1)
        result = _calculate_peaks(fake)
        assert result == ([], 0.0, 0.0, 0.0)

    def test_deprotonation_removing_only_atom_yields_empty(self):
        # "H" with [M-H]- leaves zero atoms -> final_counts empty.
        fake = _FakeDialogForPeaks("H", charge=-1, species_idx=1)
        result = _calculate_peaks(fake)
        assert result == ([], 0.0, 0.0, 0.0)


class TestCalculatePeaksIsotopePattern:
    def test_cl2_isotope_envelope_matches_known_ratios(self):
        # Classic Cl2+ isotope envelope: M : M+2 : M+4 ~= 100 : 64 : 10.2,
        # driven by real 35Cl/37Cl abundances (75.76% / 24.24%).
        fake = _FakeDialogForPeaks("Cl2", charge=1, species_idx=0)
        peaks, _avg_mw, _exact_mw, _ion_mz = _calculate_peaks(fake)

        assert len(peaks) == 3
        masses = [m for m, _p in peaks]
        intensities = [p for _m, p in peaks]

        # Sorted ascending by mass (m/z).
        assert masses == sorted(masses)
        # Isotope spacing ~= mass difference between 37Cl and 35Cl (~1.997).
        assert masses[1] - masses[0] == pytest.approx(1.99705, abs=0.01)
        assert masses[2] - masses[1] == pytest.approx(1.99705, abs=0.01)

        # Base peak normalized to 100.
        assert max(intensities) == pytest.approx(100.0, abs=1e-6)
        assert intensities[0] == pytest.approx(100.0, abs=1e-6)
        assert intensities[1] == pytest.approx(64.0, rel=0.01)
        assert intensities[2] == pytest.approx(10.24, rel=0.02)

    def test_positive_charge_one_mz_shifted_by_electron_mass(self):
        fake = _FakeDialogForPeaks("Cl2", charge=1, species_idx=0)
        peaks, _avg_mw, _exact_mw, ion_mz = _calculate_peaks(fake)
        base_peak_mz = max(peaks, key=lambda p: p[1])[0]
        # Neutral monoisotopic Cl2 mass minus one electron (charge +1, mode M).
        expected_neutral = 2 * 34.96885268
        assert base_peak_mz == pytest.approx(
            expected_neutral - _ELECTRON_MASS, abs=1e-4
        )

    def test_doubly_charged_species_divides_mz_by_charge(self):
        # C2H4, [M]2+ -> m/z should be roughly half the neutral mass.
        fake = _FakeDialogForPeaks("C2H4", charge=2, species_idx=0)
        peaks, _avg_mw, _exact_mw, ion_mz = _calculate_peaks(fake)
        base_peak_mz = max(peaks, key=lambda p: p[1])[0]

        neutral_monoisotopic = 2 * 12.0 + 4 * 1.0078250319
        expected_mz = (neutral_monoisotopic - 2 * _ELECTRON_MASS) / 2
        assert base_peak_mz == pytest.approx(expected_mz, abs=1e-4)
        assert ion_mz == pytest.approx(expected_mz, abs=1e-4)

    def test_sodium_adduct_shifts_mass_by_sodium_minus_electron(self):
        # C2H4 + Na adduct, [M+Na]+ (species_idx=2 in positive mode).
        fake = _FakeDialogForPeaks("C2H4", charge=1, species_idx=2)
        peaks, _avg_mw, _exact_mw, _ion_mz = _calculate_peaks(fake)
        base_peak_mz = max(peaks, key=lambda p: p[1])[0]

        neutral_monoisotopic = 2 * 12.0 + 4 * 1.0078250319
        expected_mz = neutral_monoisotopic + 22.98976928 - _ELECTRON_MASS
        assert base_peak_mz == pytest.approx(expected_mz, abs=1e-4)

    def test_intensities_all_normalized_le_100(self):
        fake = _FakeDialogForPeaks("Cl2", charge=1, species_idx=0)
        peaks, *_ = _calculate_peaks(fake)
        assert all(0.0 < p <= 100.0 for _m, p in peaks)


# ---------------------------------------------------------------------------
# HistogramWidget — pure view-state math (reset_view, wheelEvent, panning)
# ---------------------------------------------------------------------------


class _FakeEventPos:
    def __init__(self, x):
        self._x = x

    def x(self):
        return self._x


class _FakeAngleDelta:
    def __init__(self, y):
        self._y = y

    def y(self):
        return self._y


class _FakeWheelEvent:
    def __init__(self, delta_y, mouse_x):
        self._delta_y = delta_y
        self._mouse_x = mouse_x

    def angleDelta(self):
        return _FakeAngleDelta(self._delta_y)

    def position(self):
        return _FakeEventPos(self._mouse_x)


class _FakeMouseEvent:
    def __init__(self, x, button=1):
        self._x = x
        self._button = button

    def position(self):
        return _FakeEventPos(self._x)

    def button(self):
        return self._button


class _FakeHistogramWidget:
    """Minimal stand-in exposing only what reset_view/wheelEvent touch."""

    def __init__(self, peaks=None, stick_peaks=None, width=600):
        self.peaks = peaks or []
        self.stick_peaks = stick_peaks or []
        self.view_min = None
        self.view_max = None
        self.last_mouse_x = None
        self._width = width
        self.update_calls = 0

    def width(self):
        return self._width

    def _scaled_x_margins(self):
        # mirrors HistogramWidget._scaled_x_margins (v2026.07.17)
        scale = self.width() / 600.0
        return int(60 * scale), int(40 * scale)

    def update(self):
        self.update_calls += 1


class _FakeQt:
    class MouseButton:
        LeftButton = 1


_reset_view = extract_function(MS_PATH, "HistogramWidget", "reset_view")
_wheel_event = extract_function(MS_PATH, "HistogramWidget", "wheelEvent")
_mouse_press_event = extract_function(
    MS_PATH,
    "HistogramWidget",
    "mousePressEvent",
    extra_globals={"Qt": _FakeQt},
)
_mouse_move_event = extract_function(MS_PATH, "HistogramWidget", "mouseMoveEvent")


class TestHistogramResetView:
    def test_uses_stick_peaks_with_padding(self):
        w = _FakeHistogramWidget(stick_peaks=[(100.0, 100.0), (102.0, 50.0)])
        _reset_view(w)
        assert w.view_min == pytest.approx(100.0 - 5.0)
        assert w.view_max == pytest.approx(102.0 + 5.0)
        assert w.update_calls == 1

    def test_single_peak_gets_symmetric_padding(self):
        w = _FakeHistogramWidget(stick_peaks=[(50.0, 100.0)])
        _reset_view(w)
        assert w.view_min == pytest.approx(45.0)
        assert w.view_max == pytest.approx(55.0)

    def test_no_peaks_defaults_to_0_100(self):
        w = _FakeHistogramWidget()
        _reset_view(w)
        assert w.view_min == 0.0
        assert w.view_max == 100.0

    def test_falls_back_to_peaks_when_no_stick_peaks(self):
        w = _FakeHistogramWidget(peaks=[(10.0, 20.0), (30.0, 40.0)])
        _reset_view(w)
        assert w.view_min == pytest.approx(5.0)
        assert w.view_max == pytest.approx(35.0)


class TestHistogramWheelZoom:
    def test_scroll_up_zooms_in(self):
        w = _FakeHistogramWidget()
        w.view_min, w.view_max = 0.0, 100.0
        original_range = w.view_max - w.view_min
        event = _FakeWheelEvent(delta_y=120, mouse_x=300)  # center of 600px widget
        _wheel_event(w, event)
        new_range = w.view_max - w.view_min
        assert new_range == pytest.approx(original_range * 0.9, rel=1e-6)

    def test_scroll_down_zooms_out(self):
        w = _FakeHistogramWidget()
        w.view_min, w.view_max = 0.0, 100.0
        original_range = w.view_max - w.view_min
        event = _FakeWheelEvent(delta_y=-120, mouse_x=300)
        _wheel_event(w, event)
        new_range = w.view_max - w.view_min
        assert new_range == pytest.approx(original_range * 1.1, rel=1e-6)

    def test_zoom_clamped_to_minimum_range(self):
        w = _FakeHistogramWidget()
        w.view_min, w.view_max = 49.5, 50.5  # range = 1.0, already at floor
        event = _FakeWheelEvent(delta_y=120, mouse_x=300)
        _wheel_event(w, event)
        assert (w.view_max - w.view_min) == pytest.approx(1.0)

    def test_no_op_when_view_not_initialized(self):
        w = _FakeHistogramWidget()
        event = _FakeWheelEvent(delta_y=120, mouse_x=300)
        _wheel_event(w, event)  # must not raise
        assert w.view_min is None and w.view_max is None


class TestHistogramPan:
    def test_drag_right_shifts_view_left(self):
        w = _FakeHistogramWidget()
        w.view_min, w.view_max = 0.0, 100.0
        _mouse_press_event(w, _FakeMouseEvent(x=300, button=_FakeQt.MouseButton.LeftButton))
        assert w.last_mouse_x == 300
        _mouse_move_event(w, _FakeMouseEvent(x=350, button=_FakeQt.MouseButton.LeftButton))
        # Dragging right (dx > 0) pans the view left (both bounds decrease).
        assert w.view_min < 0.0
        assert w.view_max < 100.0
        assert (w.view_max - w.view_min) == pytest.approx(100.0)
