"""
Tests for the Elemental Analysis plugin's module-level logic: formula parsing
(including solvate notation), Hill ordering/formatting, composition math, the
'Anal. Calcd for' line and the TSV copy text, plus initialize() registration.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGIN_PATH = (
    Path(__file__).resolve().parents[1]
    / "plugins"
    / "Elemental_Analysis"
    / "elemental_analysis.py"
)

with mock_optional_imports():
    ea = load_plugin(PLUGIN_PATH)


# IUPAC-style standard atomic weights (RDKit's values for these elements).
_WEIGHTS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.067, "Cl": 35.453}
_ISOTOPES = {"H": (1, 1.00782503207), "C": (6, 12.0), "N": (7, 14.0030740048), "O": (8, 15.99491461956)}


class _FakePT:
    def GetAtomicWeight(self, sym):
        if sym not in _WEIGHTS:
            raise RuntimeError(f"Element '{sym}' not found")
        return _WEIGHTS[sym]

    def GetAtomicNumber(self, sym):
        return _ISOTOPES[sym][0]

    def GetMostCommonIsotope(self, anum):
        return anum

    def GetMassForIsotope(self, anum, _mass_number):
        return next(m for n, m in _ISOTOPES.values() if n == anum)


PT = _FakePT()


class TestParseFormula:
    def test_simple(self):
        assert ea.parse_formula("C8H10N4O2") == {"C": 8, "H": 10, "N": 4, "O": 2}

    def test_parentheses(self):
        assert ea.parse_formula("(CH3)2SO") == {"C": 2, "H": 6, "S": 1, "O": 1}

    def test_nested_parentheses(self):
        assert ea.parse_formula("C(C(CH3)3)2") == {"C": 9, "H": 18}

    def test_charge_ignored(self):
        assert ea.parse_formula("C6H5+") == {"C": 6, "H": 5}

    def test_whitespace(self):
        assert ea.parse_formula(" C6 H6 ") == {"C": 6, "H": 6}

    def test_empty(self):
        assert ea.parse_formula("") == {}

    @pytest.mark.parametrize("bad", ["C6?H6", "c6h6", "C6H6·", "·H2O"])
    def test_invalid(self, bad):
        assert ea.parse_formula(bad) is None

    def test_hemihydrate(self):
        assert ea.parse_formula("C6H6·0.5H2O") == {"C": 6, "H": 7, "O": 0.5}

    def test_integer_solvate_and_star_separator(self):
        assert ea.parse_formula("CuSO4*5H2O") == {"Cu": 1, "S": 1, "O": 9, "H": 10}

    def test_integral_float_counts_become_int(self):
        counts = ea.parse_formula("C2H6O·2H2O")
        assert counts == {"C": 2, "H": 10, "O": 3}
        assert all(isinstance(v, int) for v in counts.values())


class TestFormatting:
    def test_hill_order_with_carbon(self):
        assert ea.hill_order(["O", "N", "H", "C", "Br"]) == ["C", "H", "Br", "N", "O"]

    def test_hill_order_without_carbon(self):
        assert ea.hill_order(["O", "S", "H", "Cu"]) == ["Cu", "H", "O", "S"]

    def test_format_formula_omits_ones(self):
        assert ea.format_formula({"O": 1, "C": 1, "H": 4}) == "CH4O"

    def test_format_count(self):
        assert ea.format_count(3) == "3"
        assert ea.format_count(2.0) == "2"
        assert ea.format_count(0.5) == "0.5"
        assert ea.format_count(1.25) == "1.25"


class TestComposition:
    def test_caffeine(self):
        rows, mw = ea.calc_composition(PT, {"C": 8, "H": 10, "N": 4, "O": 2})
        assert mw == pytest.approx(194.194, abs=1e-3)
        pct = {r[0]: r[3] for r in rows}
        assert pct["C"] == pytest.approx(49.48, abs=0.005)
        assert pct["H"] == pytest.approx(5.19, abs=0.005)
        assert pct["N"] == pytest.approx(28.85, abs=0.005)
        assert sum(pct.values()) == pytest.approx(100.0)
        assert [r[0] for r in rows] == ["C", "H", "N", "O"]

    def test_fractional_solvate(self):
        rows, mw = ea.calc_composition(PT, ea.parse_formula("C8H10N4O2·0.5H2O"))
        assert mw == pytest.approx(194.194 + 0.5 * (2 * 1.008 + 15.999), abs=1e-3)

    def test_deuterium_uses_isotope_mass(self):
        rows, mw = ea.calc_composition(PT, {"C": 1, "D": 4})
        assert mw == pytest.approx(12.011 + 4 * 2.0141017781)

    def test_unknown_element_raises(self):
        with pytest.raises(RuntimeError):
            ea.calc_composition(PT, {"Xx": 1})

    def test_empty(self):
        assert ea.calc_composition(PT, {}) == ([], 0.0)

    def test_monoisotopic_mass(self):
        assert ea.monoisotopic_mass(PT, {"C": 8, "H": 10, "N": 4, "O": 2}) == pytest.approx(
            194.0804, abs=1e-4
        )


class TestCopyText:
    def _rows(self):
        return ea.calc_composition(PT, {"C": 8, "H": 10, "N": 4, "O": 2})[0]

    def test_analysis_line(self):
        assert (
            ea.analysis_line("C8H10N4O2", self._rows())
            == "Anal. Calcd for C8H10N4O2: C, 49.48; H, 5.19; N, 28.85."
        )

    def test_analysis_line_includes_sulfur(self):
        rows = ea.calc_composition(PT, {"C": 2, "H": 6, "O": 1, "S": 1})[0]
        line = ea.analysis_line("C2H6OS", rows)
        assert "; S, " in line and "O," not in line

    def test_analysis_line_empty_without_chns(self):
        rows = ea.calc_composition(PT, {"Cl": 2})[0]
        assert ea.analysis_line("Cl2", rows) == ""

    def test_rows_to_tsv(self):
        tsv = ea.rows_to_tsv(self._rows())
        lines = tsv.split("\n")
        assert lines[0] == "Element\tCount\tMass %"
        assert lines[1] == "C\t8\t49.48"
        assert len(lines) == 5

    def test_rows_to_tsv_no_header(self):
        assert ea.rows_to_tsv(self._rows(), header=False).startswith("C\t8\t")


class TestMetadataAndInitialize:
    def test_metadata(self):
        assert ea.PLUGIN_NAME == "Elemental Analysis"
        assert ea.PLUGIN_TAGS == ["Analysis"]
        assert ea.PLUGIN_SUPPORTED_MOLEDITPY_VERSION == ">=4.0.0, <5.0.0"

    def test_solvate_formulas_parse(self):
        for _label, formula in ea.SOLVATES:
            assert ea.parse_formula(formula), formula

    def test_initialize_registers_analysis_tool(self):
        ctx = make_context()
        ea.initialize(ctx)
        ctx.add_analysis_tool.assert_called_once()
        assert ctx.add_analysis_tool.call_args[0][0] == "Elemental Analysis"

    def test_toggle_reuses_existing_window(self):
        ctx = make_context()
        ea.initialize(ctx)
        toggle = ctx.add_analysis_tool.call_args[0][1]
        win = ctx.get_window.return_value
        toggle()
        win.show.assert_called_once()
        win.check_update.assert_called_once()

    def test_standalone_context(self):
        sc = ea.StandaloneContext()
        assert sc.get_main_window() is None
        sc.register_window("main_panel", "w")
        assert sc.get_window("main_panel") == "w"
