"""
Extended tests for the analysis plugins:
  - ORCA_Freq_Analyzer     / orca_out_freq_analyzer.py
  - Gaussian_Freq_Analyzer / gaussian_fchk_freq_analyzer.py
  - Gaussian_FCHK_Loader   / gaussian_fchk_loader.py
  - Gaussian_MO_Analyzer   / gaussian_fchk_mo_analyzer/ (folder plugin)
  - Symmetry_Analyzer      / symmetry_analyzer.py
  - Molecule_Comparator    / molecule_comparator.py
  - MS_Spectrum_Simulation_Neo / ms_spectrum_neo.py
  - Compound_Info_Report   / compound_info_report.py

Complements test_plugin_parsers_exports.py (basic OrcaParser / FCHKParser),
test_plugin_ui_misc.py (MS formula helpers), test_plugin_pubchem_and_paste.py
(Compound Info Report synonym/CID fetchers) and test_plugin_advanced_and_misc.py
(Symmetry worker smoke). Focus here: normal-mode / IR-intensity alignment,
imaginary frequencies, unit conversions, MO-analyzer FCHK reading and cube
writing, symmetry-operation classification, and PubChem experimental-property
extraction.

Methods that live on Qt-derived classes (mocked bases) are extracted via
``ast.get_source_segment`` + ``exec`` into an isolated namespace with a small
pure-Python numpy stub, following the technique used in
test_plugin_advanced_and_misc.py.
"""

from __future__ import annotations

import ast
import importlib.util
import io
import json
import math
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ORCA_PATH = PLUGINS_DIR / "ORCA_Freq_Analyzer" / "orca_out_freq_analyzer.py"
GAUSS_FREQ_PATH = PLUGINS_DIR / "Gaussian_Freq_Analyzer" / "gaussian_fchk_freq_analyzer.py"
FCHK_LOADER_PATH = PLUGINS_DIR / "Gaussian_FCHK_Loader" / "gaussian_fchk_loader.py"
MO_PKG_DIR = PLUGINS_DIR / "Gaussian_MO_Analyzer" / "gaussian_fchk_mo_analyzer"
MO_ANALYZER_PATH = MO_PKG_DIR / "analyzer.py"
SYMMETRY_PATH = PLUGINS_DIR / "Symmetry_Analyzer" / "symmetry_analyzer.py"
COMPARATOR_PATH = PLUGINS_DIR / "Molecule_Comparator" / "molecule_comparator.py"
MS_NEO_PATH = PLUGINS_DIR / "MS_Spectrum_Simulation_Neo" / "ms_spectrum_neo.py"
COMPOUND_PATH = PLUGINS_DIR / "Compound_Info_Report" / "compound_info_report.py"

with mock_optional_imports():
    _orca = load_plugin(ORCA_PATH)
    _gfreq = load_plugin(GAUSS_FREQ_PATH)
    _fchk_loader = load_plugin(FCHK_LOADER_PATH)
    _compound = load_plugin(COMPOUND_PATH)

# Avoid the lazy `from rdkit.Chem import ...` inside FCHKParser.parse firing
# outside the mock context (same trick as test_plugin_parsers_exports.py).
_gfreq.Chem = None


# ---------------------------------------------------------------------------
# Source-extraction helpers (for methods on Qt-derived classes)
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


def _extract_class_source(path: Path, class_name: str) -> str:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            return ast.get_source_segment(source, node)
    raise AssertionError(f"class {class_name} not found in {path.name}")


def _make_function(src: str, namespace: dict):
    exec(src, namespace)  # noqa: S102 - test-only extraction
    name = src.split("def ", 1)[1].split("(", 1)[0]
    return namespace[name]


# ===========================================================================
# ORCA Freq Analyzer — normal modes, IR intensities, imaginary frequencies
# ===========================================================================

_ORCA_FULL = """\
O   R   C   A

Total Charge           Charge          ....   -1
Multiplicity           Mult            ....    3

CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.000000    0.000000    0.000000
  H      0.000000    0.759337    0.596043
  H      0.000000   -0.759337    0.596043

-----------------------
VIBRATIONAL FREQUENCIES
-----------------------
Scaling factor for frequencies =  1.000000000 (already applied!)

   0:       0.00 cm**-1
   1:       0.00 cm**-1
   2:       0.00 cm**-1
   3:       0.00 cm**-1
   4:       0.00 cm**-1
   5:       0.00 cm**-1
   6:   -1609.85 cm**-1
   7:    3681.19 cm**-1
   8:    3811.12 cm**-1

------------
NORMAL MODES
------------
These modes are printed in Cartesian displacements
mass weighted comment line
filler line three
filler line four
filler line five

                  0          1          2          3          4          5
      0       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      1       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      2       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      3       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      4       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      5       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      6       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      7       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      8       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
                  6          7          8
      0       0.010000   0.020000   0.030000
      1       0.040000   0.050000   0.060000
      2       0.070000   0.080000   0.090000
      3       0.100000   0.110000   0.120000
      4       0.130000   0.140000   0.150000
      5       0.160000   0.170000   0.180000
      6       0.190000   0.200000   0.210000
      7       0.220000   0.230000   0.240000
      8       0.250000   0.260000   0.270000

-----------
IR SPECTRUM
-----------

 Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
   6:  -1609.85   0.015725   79.47  0.002871  (-0.001018 -0.053574 -0.000000)
   7:   3681.19   0.005000   10.00  0.001000  ( 0.001000  0.002000  0.000000)
   8:   3811.12   0.001000    5.00  0.000500  ( 0.000000  0.001000  0.000000)

The first frequency considered to be a vibration is 6
"""


class TestOrcaFinalModes:
    def _parsed(self, tmp_path, content=_ORCA_FULL):
        f = tmp_path / "full.out"
        f.write_text(content)
        p = _orca.OrcaParser()
        p.parse(str(f))
        return p

    def test_charge_parsed(self, tmp_path):
        assert self._parsed(tmp_path).charge == -1

    def test_multiplicity_parsed(self, tmp_path):
        assert self._parsed(tmp_path).multiplicity == 3

    def test_three_vibrational_modes_kept(self, tmp_path):
        p = self._parsed(tmp_path)
        assert len(p.final_modes) == 3

    def test_translational_modes_filtered(self, tmp_path):
        p = self._parsed(tmp_path)
        assert all(abs(m["freq"]) >= 10.0 for m in p.final_modes)

    def test_imaginary_frequency_preserved_negative(self, tmp_path):
        p = self._parsed(tmp_path)
        assert p.final_modes[0]["freq"] == pytest.approx(-1609.85)

    def test_intensities_aligned_with_modes(self, tmp_path):
        p = self._parsed(tmp_path)
        intens = [m["intensity"] for m in p.final_modes]
        assert intens == pytest.approx([79.47, 10.00, 5.00])

    def test_mode_vector_has_one_triple_per_atom(self, tmp_path):
        p = self._parsed(tmp_path)
        for mode in p.final_modes:
            assert len(mode["vector"]) == 3
            assert all(len(v) == 3 for v in mode["vector"])

    def test_mode_vector_values_column_aligned(self, tmp_path):
        p = self._parsed(tmp_path)
        # Mode 6 = first column of the second NORMAL MODES block
        vec = p.final_modes[0]["vector"]
        assert vec[0] == pytest.approx((0.01, 0.04, 0.07))
        assert vec[1] == pytest.approx((0.10, 0.13, 0.16))
        assert vec[2] == pytest.approx((0.19, 0.22, 0.25))

    def test_second_mode_vector_values(self, tmp_path):
        p = self._parsed(tmp_path)
        vec = p.final_modes[1]["vector"]
        assert vec[0] == pytest.approx((0.02, 0.05, 0.08))

    def test_last_ir_spectrum_block_wins(self, tmp_path):
        extra = textwrap.dedent("""\

            -----------
            IR SPECTRUM
            -----------

             Mode   freq       eps      Int      T**2         TX        TY        TZ
                   cm**-1   L/(mol*cm) km/mol    a.u.
            ----------------------------------------------------------------------------
               6:  -1609.85   0.015725   99.00  0.002871  (-0.001018 -0.053574 -0.000000)
               7:   3681.19   0.005000   11.00  0.001000  ( 0.001000  0.002000  0.000000)
               8:   3811.12   0.001000    6.00  0.000500  ( 0.000000  0.001000  0.000000)

            The first frequency considered to be a vibration is 6
        """)
        p = self._parsed(tmp_path, _ORCA_FULL + extra)
        intens = [m["intensity"] for m in p.final_modes]
        assert intens == pytest.approx([99.00, 11.00, 6.00])

    def test_mode_without_ir_entry_gets_none_intensity(self, tmp_path):
        # Drop mode 8 from the IR block: its intensity must be None, not shifted
        content = _ORCA_FULL.replace(
            "   8:   3811.12   0.001000    5.00  0.000500  ( 0.000000  0.001000  0.000000)\n",
            "",
        )
        p = self._parsed(tmp_path, content)
        assert p.final_modes[2]["intensity"] is None
        assert p.final_modes[0]["intensity"] == pytest.approx(79.47)


class TestOrcaContextRegression:
    """Regression: OrcaOutFreqAnalyzer used self.context (never assigned) for
    plotter access; vectors / GIF export / geometry reset raised or silently
    failed. Fixed to use the module-level PLUGIN_CONTEXT (2026.07.07)."""

    def test_no_self_context_attribute_in_analyzer_class(self):
        source = ORCA_PATH.read_text(encoding="utf-8")
        tree = ast.parse(source)
        offenders = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == "OrcaOutFreqAnalyzer":
                for sub in ast.walk(node):
                    if (
                        isinstance(sub, ast.Attribute)
                        and sub.attr == "context"
                        and isinstance(sub.value, ast.Name)
                        and sub.value.id == "self"
                    ):
                        offenders.append(sub.lineno)
        assert offenders == [], f"self.context used at lines {offenders}"

    def test_remove_vectors_uses_plugin_context_plotter(self):
        src = _extract_method_source(ORCA_PATH, "OrcaOutFreqAnalyzer", "remove_vectors")
        ctx = MagicMock()
        ns = {"PLUGIN_CONTEXT": ctx, "logging": MagicMock()}
        remove_vectors = _make_function(src, ns)
        fake_self = SimpleNamespace(vector_actor="actor-1")
        remove_vectors(fake_self)
        ctx.plotter.remove_actor.assert_called_once_with("actor-1")
        assert fake_self.vector_actor is None

    def test_remove_vectors_noop_without_actor(self):
        src = _extract_method_source(ORCA_PATH, "OrcaOutFreqAnalyzer", "remove_vectors")
        ctx = MagicMock()
        ns = {"PLUGIN_CONTEXT": ctx, "logging": MagicMock()}
        remove_vectors = _make_function(src, ns)
        fake_self = SimpleNamespace(vector_actor=None)
        remove_vectors(fake_self)
        ctx.plotter.remove_actor.assert_not_called()

    def test_plugin_version_bumped(self):
        # Must be newer than the last release with the self.context bug.
        assert _orca.PLUGIN_VERSION > "2026.06.27"


class TestOrcaInitializeExtras:
    def test_log_opener_registered(self):
        with mock_optional_imports():
            mod = load_plugin(ORCA_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        exts = [c[0][0] for c in ctx.register_file_opener.call_args_list]
        assert ".log" in exts

    def test_drop_handler_registered_with_priority(self):
        with mock_optional_imports():
            mod = load_plugin(ORCA_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        assert ctx.register_drop_handler.call_args.kwargs.get("priority") == 10

    def test_drop_handler_rejects_other_extensions(self):
        with mock_optional_imports():
            mod = load_plugin(ORCA_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        handler = ctx.register_drop_handler.call_args[0][0]
        assert handler("molecule.xyz") is False

    def test_drop_handler_rejects_non_orca_out_file(self, tmp_path):
        f = tmp_path / "gauss.out"
        f.write_text("Entering Gaussian System\n")
        with mock_optional_imports():
            mod = load_plugin(ORCA_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        handler = ctx.register_drop_handler.call_args[0][0]
        assert handler(str(f)) is False


# ===========================================================================
# Gaussian Freq Analyzer — FCHKParser extended
# ===========================================================================

_FCHK_TWO_MODES = textwrap.dedent("""\
    Charge                 I                -2
    Multiplicity           I                3
    Atomic numbers                             I   N=           3
         8     1     1
    Real atomic weights                        R   N=           3
       1.5999491E+01  1.0078250E+00  1.0078250E+00
    Current cartesian coordinates              R   N=           9
       0.000000E+00  0.000000E+00  0.000000E+00
       0.000000E+00  1.434600E+00  1.126430E+00
       0.000000E+00 -1.434600E+00  1.126430E+00
    Vib-Modes                                  R   N=          18
       0.01  0.02  0.03  0.04  0.05  0.06  0.07  0.08  0.09
       0.11  0.12  0.13  0.14  0.15  0.16  0.17  0.18  0.19
    Vib-E2                                     R   N=           8
       1609.85  3681.19
       1.000000  2.000000
       10.000000  20.000000
       90.0  80.0
""")


class TestFCHKParserExtended:
    def _parse(self, tmp_path, content):
        f = tmp_path / "mol.fchk"
        f.write_text(content)
        p = _gfreq.FCHKParser()
        p.parse(str(f))
        return p

    def test_n_modes_derived_from_vib_modes_size(self, tmp_path):
        p = self._parse(tmp_path, _FCHK_TWO_MODES)
        assert len(p.frequencies) == 2
        assert p.frequencies == pytest.approx([1609.85, 3681.19])

    def test_intensities_fourth_block_of_two_modes(self, tmp_path):
        p = self._parse(tmp_path, _FCHK_TWO_MODES)
        assert p.intensities == pytest.approx([90.0, 80.0])

    def test_vib_modes_structured_per_atom(self, tmp_path):
        p = self._parse(tmp_path, _FCHK_TWO_MODES)
        assert len(p.vib_modes) == 2
        assert p.vib_modes[0] == [
            pytest.approx((0.01, 0.02, 0.03)),
            pytest.approx((0.04, 0.05, 0.06)),
            pytest.approx((0.07, 0.08, 0.09)),
        ]
        assert p.vib_modes[1][0] == pytest.approx((0.11, 0.12, 0.13))

    def test_negative_charge_parsed(self, tmp_path):
        p = self._parse(tmp_path, _FCHK_TWO_MODES)
        assert p.charge == -2

    def test_multiplicity_parsed(self, tmp_path):
        p = self._parse(tmp_path, _FCHK_TWO_MODES)
        assert p.multiplicity == 3

    def test_masses_from_real_atomic_weights(self, tmp_path):
        p = self._parse(tmp_path, _FCHK_TWO_MODES)
        assert len(p.masses) == 3
        assert p.masses[0] == pytest.approx(15.999491, abs=1e-5)

    def test_vibe2_fallback_without_vib_modes_uses_3n_minus_6(self, tmp_path):
        content = textwrap.dedent("""\
            Atomic numbers                             I   N=           3
                 8     1     1
            Current cartesian coordinates              R   N=           9
               0.0  0.0  0.0
               0.0  1.4  1.1
               0.0 -1.4  1.1
            Vib-E2                                     R   N=          12
               1609.85  3681.19  3811.12
               1.0  2.0  3.0
               10.0  20.0  30.0
               90.0  80.0  70.0
        """)
        p = self._parse(tmp_path, content)
        # 3 atoms -> 3N-6 = 3 modes
        assert p.frequencies == pytest.approx([1609.85, 3681.19, 3811.12])
        assert p.intensities == pytest.approx([90.0, 80.0, 70.0])

    def test_separate_ir_inten_section_converted_from_au(self, tmp_path):
        content = textwrap.dedent("""\
            Atomic numbers                             I   N=           2
                 6     8
            Current cartesian coordinates              R   N=           6
               0.0  0.0  0.0
               0.0  0.0  2.1
            IR Inten                                   R   N=           1
               2.0
        """)
        p = self._parse(tmp_path, content)
        assert p.intensities == pytest.approx([2.0 * 974.868])

    def test_incomplete_coordinate_triple_dropped(self, tmp_path):
        content = textwrap.dedent("""\
            Atomic numbers                             I   N=           3
                 8     1     1
            Current cartesian coordinates              R   N=           8
               0.0  0.0  0.0
               0.0  1.4  1.1
               0.0 -1.4
        """)
        p = self._parse(tmp_path, content)
        assert len(p.coords) == 2

    def test_bohr_conversion_exact_value(self, tmp_path):
        p = self._parse(tmp_path, _FCHK_TWO_MODES)
        assert p.coords[1][1] == pytest.approx(1.434600 * 0.529177210903, abs=1e-9)


# ===========================================================================
# Gaussian FCHK Loader — load_module_from_path + initialize contract
# ===========================================================================

class TestLoadModuleFromPath:
    def test_loads_valid_module(self, tmp_path):
        mod_file = tmp_path / "mini_plugin_mod.py"
        mod_file.write_text("ANSWER = 42\n")
        mod = _fchk_loader.load_module_from_path("mini_plugin_mod_t1", str(mod_file))
        assert mod is not None
        assert mod.ANSWER == 42

    def test_module_registered_in_sys_modules(self, tmp_path):
        mod_file = tmp_path / "mini_plugin_mod2.py"
        mod_file.write_text("X = 1\n")
        _fchk_loader.load_module_from_path("mini_plugin_mod_t2", str(mod_file))
        assert "mini_plugin_mod_t2" in sys.modules
        del sys.modules["mini_plugin_mod_t2"]

    def test_syntax_error_returns_none(self, tmp_path):
        mod_file = tmp_path / "broken_mod.py"
        mod_file.write_text("def broken(:\n")
        assert _fchk_loader.load_module_from_path("broken_mod_t", str(mod_file)) is None

    def test_missing_file_returns_none(self, tmp_path):
        missing = tmp_path / "ghost.py"
        assert _fchk_loader.load_module_from_path("ghost_mod_t", str(missing)) is None


class TestFCHKLoaderInitialize:
    def _init(self):
        with mock_optional_imports():
            mod = load_plugin(FCHK_LOADER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        return ctx

    def test_registers_three_extensions_priority_100(self):
        ctx = self._init()
        calls = ctx.register_file_opener.call_args_list
        exts = {c[0][0] for c in calls}
        assert exts == {".fchk", ".fck", ".fch"}
        assert all(c.kwargs.get("priority") == 100 for c in calls)

    def test_drop_handler_accepts_fchk_case_insensitive(self):
        ctx = self._init()
        handler = ctx.register_drop_handler.call_args[0][0]
        assert handler("C:/data/RESULT.FCHK") is True

    def test_drop_handler_rejects_other_files(self):
        ctx = self._init()
        handler = ctx.register_drop_handler.call_args[0][0]
        assert handler("C:/data/result.xyz") is False

    def test_drop_handler_priority_100(self):
        ctx = self._init()
        assert ctx.register_drop_handler.call_args.kwargs.get("priority") == 100


# ===========================================================================
# Gaussian MO Analyzer — FCHKReader (extracted with numpy stub)
# ===========================================================================


class _FakeNP:
    """Minimal numpy stand-in for FCHKReader / normalization prefactor."""

    pi = math.pi

    @staticmethod
    def array(values, dtype=None):
        if dtype is int:
            return [int(v) for v in values]
        if dtype is float:
            return [float(v) for v in values]
        return list(values)

    @staticmethod
    def sqrt(x):
        return math.sqrt(x)


def _make_fchk_reader_class():
    src = _extract_class_source(MO_ANALYZER_PATH, "FCHKReader")
    ns = {"np": _FakeNP(), "re": __import__("re")}
    exec(src, ns)  # noqa: S102
    return ns["FCHKReader"]


class TestMOFCHKReader:
    def _read(self, tmp_path, content):
        f = tmp_path / "mo.fchk"
        f.write_text(content)
        return _make_fchk_reader_class()(str(f))

    def test_integer_array_parsed(self, tmp_path):
        r = self._read(
            tmp_path,
            "Atomic numbers                             I   N=           3\n"
            "           8           1           1\n",
        )
        assert r.get("Atomic numbers") == [8, 1, 1]

    def test_float_array_parsed(self, tmp_path):
        r = self._read(
            tmp_path,
            "Alpha Orbital Energies                     R   N=           2\n"
            "  -1.02500000E+00   3.75000000E-01\n",
        )
        assert r.get("Alpha Orbital Energies") == pytest.approx([-1.025, 0.375])

    def test_scalar_integer_parsed(self, tmp_path):
        r = self._read(tmp_path, "Number of electrons                        I           10\n")
        assert r.get("Number of electrons") == [10]

    def test_scalar_negative_integer(self, tmp_path):
        r = self._read(tmp_path, "Charge                                     I           -2\n")
        assert r.get("Charge") == [-2]

    def test_fortran_d_exponent_converted(self, tmp_path):
        r = self._read(
            tmp_path,
            "Primitive exponents                        R   N=           2\n"
            "  1.30000000D+01  2.50000000D-01\n",
        )
        assert r.get("Primitive exponents") == pytest.approx([13.0, 0.25])

    def test_packed_negative_numbers_split(self, tmp_path):
        # Fortran occasionally packs values: "0.5-0.25" must split into two
        r = self._read(
            tmp_path,
            "MO coefficients                            R   N=           2\n"
            "  0.5-0.25\n",
        )
        assert r.get("MO coefficients") == pytest.approx([0.5, -0.25])

    def test_values_truncated_to_declared_count(self, tmp_path):
        r = self._read(
            tmp_path,
            "Shell types                                I   N=           2\n"
            "           0           1           2\n",
        )
        assert r.get("Shell types") == [0, 1]

    def test_get_default_for_missing_key(self, tmp_path):
        r = self._read(tmp_path, "Number of electrons                        I           10\n")
        assert r.get("Nope", default="fallback") == "fallback"

    def test_multiline_array_accumulated(self, tmp_path):
        r = self._read(
            tmp_path,
            "Current cartesian coordinates              R   N=           6\n"
            "  1.0  2.0  3.0\n"
            "  4.0  5.0  6.0\n",
        )
        assert r.get("Current cartesian coordinates") == pytest.approx(
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        )


class TestMONormalizationPrefactor:
    def _prefactor(self):
        src = _extract_method_source(
            MO_ANALYZER_PATH, "BasisSetEngine", "_normalization_prefactor"
        )
        return _make_function(src, {"np": _FakeNP()})

    def test_s_function_value(self):
        fn = self._prefactor()
        alpha = 0.5
        expected = (2 * alpha / math.pi) ** 0.75
        assert fn(None, alpha, 0, 0, 0) == pytest.approx(expected)

    def test_p_function_value(self):
        fn = self._prefactor()
        alpha = 1.2
        expected = (2 * alpha / math.pi) ** 0.75 * math.sqrt(8 * alpha / 2)
        assert fn(None, alpha, 1, 0, 0) == pytest.approx(expected)

    def test_d_xx_function_value(self):
        fn = self._prefactor()
        alpha = 0.8
        expected = (2 * alpha / math.pi) ** 0.75 * math.sqrt((8 * alpha) ** 2 * 2 / 24)
        assert fn(None, alpha, 2, 0, 0) == pytest.approx(expected)

    def test_symmetric_in_lmn_permutation(self):
        fn = self._prefactor()
        assert fn(None, 0.9, 1, 1, 0) == pytest.approx(fn(None, 0.9, 0, 1, 1))


class _FakeGrid:
    def __init__(self, shape, values):
        self.shape = shape
        self._values = values

    def flatten(self):
        return list(self._values)


class _FakeMat:
    def __init__(self, rows):
        self._rows = rows

    def __getitem__(self, idx):
        i, j = idx
        return self._rows[i][j]


class TestMOCubeWriter:
    def _write(self, tmp_path, n_vals):
        src = _extract_class_source(MO_ANALYZER_PATH, "CubeWriter")
        ns = {}
        exec(src, ns)  # noqa: S102
        writer = ns["CubeWriter"]
        out = tmp_path / "orbital.cube"
        data = _FakeGrid((1, 2, 3) if n_vals == 6 else (1, 1, n_vals),
                         [float(i) for i in range(n_vals)])
        vectors = _FakeMat([[0.2, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.2]])
        writer.write(
            str(out),
            atoms=[(0.0, 0.0, 0.0)],
            atom_nos=[8],
            origin=(-1.0, -1.0, -1.0),
            vectors=vectors,
            data=data,
            comment="HOMO",
        )
        return out.read_text().splitlines()

    def test_header_and_atom_lines(self, tmp_path):
        lines = self._write(tmp_path, 6)
        assert "HOMO" in lines[0]
        # natoms + origin line
        parts = lines[2].split()
        assert int(parts[0]) == 1
        assert float(parts[1]) == pytest.approx(-1.0)
        # atom line: Z, charge, x, y, z
        atom_parts = lines[6].split()
        assert int(atom_parts[0]) == 8
        assert float(atom_parts[1]) == pytest.approx(8.0)

    def test_grid_dimension_lines(self, tmp_path):
        lines = self._write(tmp_path, 6)
        assert int(lines[3].split()[0]) == 1
        assert int(lines[4].split()[0]) == 2
        assert int(lines[5].split()[0]) == 3

    def test_six_values_per_data_line(self, tmp_path):
        lines = self._write(tmp_path, 6)
        data_lines = lines[7:]
        assert len(data_lines) == 1
        assert len(data_lines[0].split()) == 6

    def test_partial_last_line_when_not_multiple_of_six(self, tmp_path):
        lines = self._write(tmp_path, 7)
        data_lines = lines[7:]
        assert len(data_lines) == 2
        assert len(data_lines[0].split()) == 6
        assert len(data_lines[1].split()) == 1

    def test_values_scientific_format_roundtrip(self, tmp_path):
        lines = self._write(tmp_path, 6)
        vals = [float(v) for v in lines[7].split()]
        assert vals == pytest.approx([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])


# ===========================================================================
# Gaussian MO Analyzer — package initialize() contract
# ===========================================================================


def _load_mo_package():
    pkg_name = "gaussian_fchk_mo_analyzer_testpkg"
    spec = importlib.util.spec_from_file_location(
        pkg_name,
        MO_PKG_DIR / "__init__.py",
        submodule_search_locations=[str(MO_PKG_DIR)],
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[pkg_name] = mod
    try:
        spec.loader.exec_module(mod)
    except Exception:
        sys.modules.pop(pkg_name, None)
        raise
    return pkg_name, mod


class TestMOPackageInitialize:
    def _init(self):
        with mock_optional_imports():
            pkg_name, mod = _load_mo_package()
            try:
                ctx = make_context()
                mod.initialize(ctx)
            finally:
                for k in list(sys.modules):
                    if k.startswith(pkg_name):
                        del sys.modules[k]
        return mod, ctx

    def test_registers_three_extensions_priority_10(self):
        _, ctx = self._init()
        calls = ctx.register_file_opener.call_args_list
        exts = {c[0][0] for c in calls}
        assert exts == {".fchk", ".fck", ".fch"}
        assert all(c.kwargs.get("priority") == 10 for c in calls)

    def test_drop_handler_rejects_non_fchk(self):
        _, ctx = self._init()
        handler = ctx.register_drop_handler.call_args[0][0]
        assert handler("thing.cube") is False

    def test_drop_handler_accepts_fchk_and_opens_widget(self):
        mod, ctx = self._init()
        mod.OrbitalWidget = MagicMock()
        handler = ctx.register_drop_handler.call_args[0][0]
        assert handler("job.fchk") is True
        mod.OrbitalWidget.assert_called_once()

    def test_version_constant_present(self):
        with mock_optional_imports():
            pkg_name, mod = _load_mo_package()
            for k in list(sys.modules):
                if k.startswith(pkg_name):
                    del sys.modules[k]
        assert mod.PLUGIN_VERSION
        assert mod.PLUGIN_NAME == "Gaussian MO Analyzer"


# ===========================================================================
# Symmetry Analyzer — operation classification (numpy stub)
# ===========================================================================


class _SymNP:
    @staticmethod
    def trace(m):
        return m[0][0] + m[1][1] + m[2][2]

    @staticmethod
    def eye(n):
        return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    @staticmethod
    def allclose(a, b, atol=1e-8):
        return all(
            abs(a[i][j] - b[i][j]) <= atol for i in range(3) for j in range(3)
        )

    @staticmethod
    def isclose(a, b, atol=1e-8):
        return abs(a - b) <= atol

    @staticmethod
    def clip(v, lo, hi):
        return max(lo, min(hi, v))

    @staticmethod
    def degrees(x):
        return math.degrees(x)

    @staticmethod
    def arccos(x):
        return math.acos(x)

    class linalg:
        @staticmethod
        def det(m):
            return (
                m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
            )


_S3 = math.sqrt(3.0) / 2.0

_MAT_E = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
_MAT_C2 = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]
_MAT_C3 = [[-0.5, -_S3, 0.0], [_S3, -0.5, 0.0], [0.0, 0.0, 1.0]]
_MAT_SIGMA = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]
_MAT_INV = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]
_MAT_S4 = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]


def _op(matrix):
    return SimpleNamespace(rotation_matrix=matrix)


class TestSymmetrySortKey:
    def _key_fn(self):
        src = _extract_method_source(
            SYMMETRY_PATH, "SymmetryAnalysisPlugin", "_get_op_sort_key"
        )
        return _make_function(src, {"np": _SymNP()})

    def test_identity_first(self):
        assert self._key_fn()(None, _op(_MAT_E)) == (0, 0)

    def test_c2_rotation(self):
        assert self._key_fn()(None, _op(_MAT_C2)) == (1, -2)

    def test_c3_rotation(self):
        assert self._key_fn()(None, _op(_MAT_C3)) == (1, -3)

    def test_reflection(self):
        assert self._key_fn()(None, _op(_MAT_SIGMA)) == (2, 0)

    def test_inversion(self):
        assert self._key_fn()(None, _op(_MAT_INV)) == (3, 0)

    def test_improper_s4(self):
        assert self._key_fn()(None, _op(_MAT_S4)) == (4, -4)

    def test_overall_ordering(self):
        fn = self._key_fn()
        ops = [_MAT_S4, _MAT_INV, _MAT_C2, _MAT_SIGMA, _MAT_E, _MAT_C3]
        keys = sorted(fn(None, _op(m)) for m in ops)
        # E -> C3 -> C2 -> sigma -> i -> S4 (higher rotation order first)
        assert keys == [(0, 0), (1, -3), (1, -2), (2, 0), (3, 0), (4, -4)]


class TestSymmetryOpLabel:
    def _label_fn(self):
        src = _extract_method_source(
            SYMMETRY_PATH, "SymmetryAnalysisPlugin", "_get_op_label"
        )
        return _make_function(src, {"np": _SymNP()})

    def test_identity_label(self):
        assert self._label_fn()(None, _op(_MAT_E), 0) == "#1: E (Identity)"

    def test_c2_label_with_subscript(self):
        assert self._label_fn()(None, _op(_MAT_C2), 1) == "#2: C₂ (Rotation)"

    def test_c3_label(self):
        assert self._label_fn()(None, _op(_MAT_C3), 0) == "#1: C₃ (Rotation)"

    def test_reflection_label(self):
        assert self._label_fn()(None, _op(_MAT_SIGMA), 2) == "#3: σ (Reflection)"

    def test_inversion_label(self):
        assert self._label_fn()(None, _op(_MAT_INV), 0) == "#1: i (Inversion)"

    def test_s4_label(self):
        assert self._label_fn()(None, _op(_MAT_S4), 0) == "#1: S₄ (Improper Rotation)"


class TestFormatSymmetrySymbol:
    def _fmt(self):
        src = _extract_method_source(
            SYMMETRY_PATH, "SymmetryAnalysisPlugin", "_format_symmetry_symbol"
        )
        return _make_function(src, {})

    def test_c2v(self):
        assert self._fmt()(None, "C2v") == "<html><i>C</i><sub>2v</sub></html>"

    def test_td(self):
        assert self._fmt()(None, "Td") == "<html><i>T</i><sub>d</sub></html>"

    def test_infinity_axis(self):
        assert self._fmt()(None, "C*v") == "<html><i>C</i><sub>∞v</sub></html>"

    def test_non_matching_symbol_returned_as_is(self):
        assert self._fmt()(None, "1?") == "1?"


# ===========================================================================
# Molecule Comparator — document-reset handler behaviour
# ===========================================================================

class TestComparatorResetHandler:
    def _handler(self, mw):
        with mock_optional_imports():
            mod = load_plugin(COMPARATOR_PATH)
            ctx = make_context()
            ctx.get_main_window.return_value = mw
            mod.initialize(ctx)
        return ctx.register_document_reset_handler.call_args[0][0]

    def test_reset_closes_open_window(self):
        mw = SimpleNamespace(molecule_comparator_window=MagicMock())
        handler = self._handler(mw)
        handler()
        mw.molecule_comparator_window.close.assert_called_once()

    def test_reset_noop_without_window(self):
        mw = SimpleNamespace()  # no molecule_comparator_window attribute
        handler = self._handler(mw)
        handler()  # must not raise


# ===========================================================================
# MS Spectrum Simulation Neo — Gaussian broadening (requires real numpy)
# ===========================================================================

class TestGaussianBroadening:
    def _fn(self):
        pytest.importorskip("numpy")
        src = _extract_method_source(MS_NEO_PATH, "MSSpectrumDialog", "apply_gaussian_broadening")
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


# ===========================================================================
# Compound Info Report — fetch_experimental_properties (mocked urllib)
# ===========================================================================


def _pugview_payload(props):
    return {
        "Record": {
            "Section": [
                {
                    "TOCHeading": "Chemical and Physical Properties",
                    "Section": [
                        {
                            "TOCHeading": "Experimental Properties",
                            "Section": props,
                        }
                    ],
                }
            ]
        }
    }


class _FakeResponse:
    def __init__(self, payload, status=200):
        self.status = status
        self._body = json.dumps(payload).encode("utf-8")

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class TestFetchExperimentalProperties:
    def _patch(self, monkeypatch, payload, status=200):
        monkeypatch.setattr(
            _compound.urllib.request,
            "urlopen",
            lambda url: _FakeResponse(payload, status=status),
        )

    def test_density_string_extracted(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Density",
                    "Information": [
                        {"Value": {"StringWithMarkup": [{"String": "1.234 g/cm3"}]}}
                    ],
                }
            ]
        )
        self._patch(monkeypatch, payload)
        density, desc = _compound.PubChemFetcher.fetch_experimental_properties(241)
        assert density == "1.234 g/cm3"
        assert desc is None

    def test_density_number_value(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Density",
                    "Information": [{"Value": {"Number": [0.879]}}],
                }
            ]
        )
        self._patch(monkeypatch, payload)
        density, _ = _compound.PubChemFetcher.fetch_experimental_properties(241)
        assert density == "0.879"

    def test_physical_description_matching_heuristic(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Physical Description",
                    "Information": [
                        {
                            "Value": {
                                "StringWithMarkup": [
                                    {"String": "White crystalline powder"}
                                ]
                            }
                        }
                    ],
                }
            ]
        )
        self._patch(monkeypatch, payload)
        _, desc = _compound.PubChemFetcher.fetch_experimental_properties(241)
        assert desc == "White crystalline powder"

    def test_description_without_keywords_rejected(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Physical Description",
                    "Information": [
                        {"Value": {"StringWithMarkup": [{"String": "Pungent smell"}]}}
                    ],
                }
            ]
        )
        self._patch(monkeypatch, payload)
        _, desc = _compound.PubChemFetcher.fetch_experimental_properties(241)
        assert desc is None

    def test_non_200_status_returns_none_pair(self, monkeypatch):
        self._patch(monkeypatch, _pugview_payload([]), status=404)
        result = _compound.PubChemFetcher.fetch_experimental_properties(241)
        assert result == (None, None)

    def test_no_cid_short_circuits_without_network(self, monkeypatch):
        def _explode(url):
            raise AssertionError("network must not be touched")

        monkeypatch.setattr(_compound.urllib.request, "urlopen", _explode)
        assert _compound.PubChemFetcher.fetch_experimental_properties(None) == (
            None,
            None,
        )

    def test_first_matching_density_wins(self, monkeypatch):
        payload = _pugview_payload(
            [
                {
                    "TOCHeading": "Density",
                    "Information": [
                        {"Value": {"StringWithMarkup": [{"String": "1.0 g/cm3"}]}}
                    ],
                },
                {
                    "TOCHeading": "Density",
                    "Information": [
                        {"Value": {"StringWithMarkup": [{"String": "2.0 g/cm3"}]}}
                    ],
                },
            ]
        )
        self._patch(monkeypatch, payload)
        density, _ = _compound.PubChemFetcher.fetch_experimental_properties(241)
        assert density == "1.0 g/cm3"
