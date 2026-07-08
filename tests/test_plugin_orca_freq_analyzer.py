"""
Tests for the ORCA Freq Analyzer plugin (plugins/ORCA_Freq_Analyzer/orca_out_freq_analyzer.py).

All heavy deps (PyQt6, rdkit, numpy) are stubbed via mock_optional_imports().
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
ORCA_PATH = PLUGINS_DIR / "ORCA_Freq_Analyzer" / "orca_out_freq_analyzer.py"

with mock_optional_imports():
    _orca = load_plugin(ORCA_PATH)


_WATER_ORCA = textwrap.dedent("""\
    O   R   C   A

    CARTESIAN COORDINATES (ANGSTROEM)
    ---------------------------------
      O      0.000000    0.000000    0.000000
      H      0.000000    0.759337    0.596043
      H      0.000000   -0.759337    0.596043

    VIBRATIONAL FREQUENCIES
    -----------------------
    Scaling factor for frequencies =  1.000000000 (already applied!)

       0:       0.00 cm**-1
       1:       0.00 cm**-1
       2:       0.00 cm**-1
       3:    1609.85 cm**-1
       4:    3681.19 cm**-1
       5:    3811.12 cm**-1

    NORMAL MODES
    end
""")


class TestIsValidOrcaFile:
    def test_orca_keyword_makes_file_valid(self, tmp_path):
        f = tmp_path / "calc.out"
        f.write_text("Some header line\nORCA version 5.0.3\nmore text\n")
        assert _orca.is_valid_orca_file(str(f)) is True

    def test_spaced_orca_header_makes_file_valid(self, tmp_path):
        # Alternative header format: "O   R   C   A"
        f = tmp_path / "calc.out"
        f.write_text("Some preamble\nO   R   C   A\n\n")
        assert _orca.is_valid_orca_file(str(f)) is True

    def test_no_orca_keyword_returns_false(self, tmp_path):
        f = tmp_path / "gaussian.log"
        f.write_text("Entering Gaussian System\n Link0: mem=4gb\n")
        assert _orca.is_valid_orca_file(str(f)) is False

    def test_nonexistent_file_returns_false(self, tmp_path):
        assert _orca.is_valid_orca_file(str(tmp_path / "ghost.out")) is False

    def test_empty_file_returns_false(self, tmp_path):
        f = tmp_path / "empty.out"
        f.write_text("")
        assert _orca.is_valid_orca_file(str(f)) is False

    def test_orca_within_first_500_lines(self, tmp_path):
        lines = ["line\n"] * 498 + ["ORCA job\n", "end\n"]
        f = tmp_path / "late.out"
        f.write_text("".join(lines))
        assert _orca.is_valid_orca_file(str(f)) is True

    def test_orca_after_500_lines_not_detected(self, tmp_path):
        # 501 filler lines, then ORCA — should NOT be detected
        lines = ["filler\n"] * 501 + ["ORCA\n"]
        f = tmp_path / "toolate.out"
        f.write_text("".join(lines))
        assert _orca.is_valid_orca_file(str(f)) is False

    def test_case_sensitive_orca_lowercase_not_detected(self, tmp_path):
        f = tmp_path / "lower.out"
        f.write_text("orca calculation\n")
        assert _orca.is_valid_orca_file(str(f)) is False


class TestOrcaParserParse:
    def _parsed_water(self, tmp_path):
        f = tmp_path / "water.out"
        f.write_text(_WATER_ORCA)
        p = _orca.OrcaParser()
        p.parse(str(f))
        return p

    def test_atom_symbols_parsed(self, tmp_path):
        p = self._parsed_water(tmp_path)
        assert p.atoms == ["O", "H", "H"]

    def test_coord_count_matches_atoms(self, tmp_path):
        p = self._parsed_water(tmp_path)
        assert len(p.coords) == len(p.atoms) == 3

    def test_oxygen_at_origin(self, tmp_path):
        p = self._parsed_water(tmp_path)
        x, y, z = p.coords[0]
        assert abs(x) < 1e-9 and abs(y) < 1e-9 and abs(z) < 1e-9

    def test_hydrogen_coords_correct(self, tmp_path):
        p = self._parsed_water(tmp_path)
        _, y1, _ = p.coords[1]
        _, y2, _ = p.coords[2]
        assert abs(y1 - 0.759337) < 1e-4
        assert abs(y2 + 0.759337) < 1e-4

    def test_all_six_frequencies_parsed(self, tmp_path):
        p = self._parsed_water(tmp_path)
        assert len(p.frequencies) == 6

    def test_zero_frequencies_included(self, tmp_path):
        p = self._parsed_water(tmp_path)
        assert p.frequencies[0] == 0.0
        assert p.frequencies[2] == 0.0

    def test_nonzero_frequencies_parsed(self, tmp_path):
        p = self._parsed_water(tmp_path)
        assert abs(p.frequencies[3] - 1609.85) < 0.01
        assert abs(p.frequencies[4] - 3681.19) < 0.01
        assert abs(p.frequencies[5] - 3811.12) < 0.01

    def test_multiple_coord_blocks_uses_last(self, tmp_path):
        """Parser must use the LAST CARTESIAN COORDINATES block (geometry optimisation)."""
        content = textwrap.dedent("""\
            ORCA

            CARTESIAN COORDINATES (ANGSTROEM)
            ---------------------------------
              C      0.000000    0.000000    0.000000

            CARTESIAN COORDINATES (ANGSTROEM)
            ---------------------------------
              C      1.000000    2.000000    3.000000
        """)
        f = tmp_path / "opt.out"
        f.write_text(content)
        p = _orca.OrcaParser()
        p.parse(str(f))
        assert len(p.atoms) == 1
        assert p.atoms[0] == "C"
        assert abs(p.coords[0][0] - 1.0) < 1e-9

    def test_no_coords_no_atoms(self, tmp_path):
        f = tmp_path / "nocoord.out"
        f.write_text("ORCA\nsome output without coordinates\n")
        p = _orca.OrcaParser()
        p.parse(str(f))
        assert p.atoms == []
        assert p.coords == []

    def test_no_freq_block_gives_empty_frequencies(self, tmp_path):
        content = textwrap.dedent("""\
            ORCA

            CARTESIAN COORDINATES (ANGSTROEM)
            ---------------------------------
              O      0.0  0.0  0.0
        """)
        f = tmp_path / "nofreq.out"
        f.write_text(content)
        p = _orca.OrcaParser()
        p.parse(str(f))
        assert p.frequencies == []


class TestOrcaFreqInitialize:
    def test_initialize_registers_out_file_opener(self):
        with mock_optional_imports():
            mod = load_plugin(ORCA_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        # Should register at least one file opener (.out)
        ctx.register_file_opener.assert_called()
        exts = [call[0][0] for call in ctx.register_file_opener.call_args_list]
        assert ".out" in exts

    def test_plugin_version_constant_present(self):
        assert hasattr(_orca, "PLUGIN_VERSION")
        assert _orca.PLUGIN_VERSION




# ---------------------------------------------------------------------------
# normal modes, IR intensities, imaginary frequencies, context regression
# ---------------------------------------------------------------------------

import ast
from types import SimpleNamespace
from unittest.mock import MagicMock


def _extract_method_source(path, class_name, method_name):
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == method_name:
                    return textwrap.dedent(ast.get_source_segment(source, item))
    raise AssertionError(f"{class_name}.{method_name} not found in {path.name}")


def _make_function(src, namespace):
    exec(src, namespace)  # noqa: S102 - test-only extraction
    name = src.split("def ", 1)[1].split("(", 1)[0]
    return namespace[name]


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
