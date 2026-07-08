"""
Tests for the MS Spectrum Simulation Neo plugin: parse_formula_str,
to_subscript, to_superscript, get_adduct_delta.

MSSpectrumDialog inherits a Qt base, so its methods are extracted from the
plugin source via AST and compiled standalone (self is unused in all of them).
"""

from __future__ import annotations

from pathlib import Path

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
