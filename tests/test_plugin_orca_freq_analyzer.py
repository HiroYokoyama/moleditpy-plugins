"""
Tests for the ORCA Freq Analyzer plugin (plugins/ORCA_Freq_Analyzer/orca_out_freq_analyzer.py).

All heavy deps (PyQt6, rdkit, numpy) are stubbed via mock_optional_imports().
"""

from __future__ import annotations

import textwrap
from pathlib import Path

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
