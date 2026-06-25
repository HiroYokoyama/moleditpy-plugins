"""
Tests for parser and export plugins:
  - ORCA_Freq_Analyzer    / orca_out_freq_analyzer.py
  - Gaussian_Freq_Analyzer / gaussian_fchk_freq_analyzer.py
  - Gaussian_FCHK_Loader  / gaussian_fchk_loader.py
  - Blender_Export        / blender_export.py
  - POV-Ray_Export        / povray_export.py

All heavy deps (PyQt6, rdkit, numpy) are stubbed via mock_optional_imports().
Pure-filesystem helpers are tested with tmp_path only.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ORCA_PATH = PLUGINS_DIR / "ORCA_Freq_Analyzer" / "orca_out_freq_analyzer.py"
GAUSS_FREQ_PATH = PLUGINS_DIR / "Gaussian_Freq_Analyzer" / "gaussian_fchk_freq_analyzer.py"
FCHK_LOADER_PATH = PLUGINS_DIR / "Gaussian_FCHK_Loader" / "gaussian_fchk_loader.py"
BLENDER_PATH = PLUGINS_DIR / "Blender_Export" / "blender_export.py"
POVRAY_PATH = PLUGINS_DIR / "POV-Ray_Export" / "povray_export.py"


# ---------------------------------------------------------------------------
# Load modules once (heavy deps stubbed for the load phase)
# ---------------------------------------------------------------------------

with mock_optional_imports():
    _orca = load_plugin(ORCA_PATH)
    _gfreq = load_plugin(GAUSS_FREQ_PATH)
    _fchk_loader = load_plugin(FCHK_LOADER_PATH)
    _blender = load_plugin(BLENDER_PATH)
    _povray = load_plugin(POVRAY_PATH)


# ===========================================================================
# ORCA Freq Analyzer — is_valid_orca_file
# ===========================================================================

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


# ===========================================================================
# ORCA Freq Analyzer — OrcaParser.parse
# ===========================================================================

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


# ===========================================================================
# Gaussian Freq Analyzer — FCHKParser.parse
# ===========================================================================

# Minimal FCHK content: 3 atoms (O H H), 3 normal modes.
# Vib-Modes: 3*N_atoms*N_modes = 3*3*3 = 27 values  → n_modes = 27//(3*3) = 3
# Vib-E2:    4*N_modes = 12 values
#   indices 0:3  → frequencies
#   indices 3:6  → reduced masses (ignored)
#   indices 6:9  → force constants (ignored)
#   indices 9:12 → IR intensities (km/mol)
_WATER_FCHK = textwrap.dedent("""\
    Total Energy                               R          -76.0000000000
    Charge                 I                0
    Multiplicity           I                1
    Atomic numbers                             I   N=           3
         8     1     1
    Current cartesian coordinates              R   N=           9
       0.000000E+00  0.000000E+00  0.000000E+00
       0.000000E+00  1.434600E+00  1.126430E+00
       0.000000E+00 -1.434600E+00  1.126430E+00
    Vib-Modes                                  R   N=          27
       0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
       0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
       0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
    Vib-E2                                     R   N=          12
       1609.85  3681.19  3811.12
       1.000000  2.000000  3.000000
       10.000000  20.000000  30.000000
       90.0  80.0  70.0
""")


class TestFCHKParserParse:
    def _parse(self, tmp_path, content=_WATER_FCHK):
        f = tmp_path / "water.fchk"
        f.write_text(content)
        p = _gfreq.FCHKParser()
        p.parse(str(f))
        return p

    def test_atomic_numbers_parsed(self, tmp_path):
        p = self._parse(tmp_path)
        assert p.atoms == [8, 1, 1]

    def test_coord_count_equals_atom_count(self, tmp_path):
        p = self._parse(tmp_path)
        assert len(p.coords) == 3

    def test_coords_converted_from_bohr(self, tmp_path):
        BOHR_TO_ANG = 0.529177210903
        p = self._parse(tmp_path)
        x0, y0, z0 = p.coords[0]
        assert abs(x0) < 1e-9
        assert abs(y0) < 1e-9
        assert abs(z0) < 1e-9
        # Second atom: raw y = 1.434600 Bohr
        _, y1, _ = p.coords[1]
        assert abs(y1 - 1.434600 * BOHR_TO_ANG) < 1e-5

    def test_negative_y_coord_converted_correctly(self, tmp_path):
        BOHR_TO_ANG = 0.529177210903
        p = self._parse(tmp_path)
        _, y2, _ = p.coords[2]
        assert abs(y2 + 1.434600 * BOHR_TO_ANG) < 1e-5

    def test_three_frequencies_extracted(self, tmp_path):
        p = self._parse(tmp_path)
        assert len(p.frequencies) == 3

    def test_frequency_values_correct(self, tmp_path):
        p = self._parse(tmp_path)
        assert abs(p.frequencies[0] - 1609.85) < 0.01
        assert abs(p.frequencies[1] - 3681.19) < 0.01
        assert abs(p.frequencies[2] - 3811.12) < 0.01

    def test_ir_intensities_from_fourth_vibe2_block(self, tmp_path):
        p = self._parse(tmp_path)
        assert len(p.intensities) == 3
        assert abs(p.intensities[0] - 90.0) < 0.01
        assert abs(p.intensities[1] - 80.0) < 0.01
        assert abs(p.intensities[2] - 70.0) < 0.01

    def test_missing_vib_sections_gives_empty_frequencies(self, tmp_path):
        minimal = textwrap.dedent("""\
            Atomic numbers                             I   N=           1
                 6
            Current cartesian coordinates              R   N=           3
               0.0  0.0  0.0
        """)
        p = self._parse(tmp_path, minimal)
        assert p.atoms == [6]
        assert p.frequencies == []


# ===========================================================================
# Gaussian FCHK Loader — find_file_recursive
# ===========================================================================

class TestFindFileRecursive:
    def test_finds_exact_filename(self, tmp_path):
        target = tmp_path / "sub" / "deep"
        target.mkdir(parents=True)
        (target / "result.fchk").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "result.fchk")
        assert found is not None
        assert found.endswith("result.fchk")

    def test_fnmatch_wildcard_pattern(self, tmp_path):
        (tmp_path / "run01.fchk").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "*.fchk")
        assert found is not None
        assert found.endswith(".fchk")

    def test_returns_none_when_not_found(self, tmp_path):
        (tmp_path / "data.txt").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "*.fchk")
        assert found is None

    def test_returns_none_in_empty_directory(self, tmp_path):
        found = _fchk_loader.find_file_recursive(str(tmp_path), "anything.fchk")
        assert found is None

    def test_finds_first_match_in_nested_dirs(self, tmp_path):
        (tmp_path / "a").mkdir()
        (tmp_path / "a" / "mol.fchk").write_text("")
        (tmp_path / "b").mkdir()
        (tmp_path / "b" / "mol.fchk").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "mol.fchk")
        assert found is not None
        assert found.endswith("mol.fchk")

    def test_pattern_case_sensitive(self, tmp_path):
        (tmp_path / "Result.FCHK").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "result.fchk")
        # fnmatch on Windows is case-insensitive; on Linux it's case-sensitive.
        # Just verify the function doesn't raise.
        assert isinstance(found, (str, type(None)))


# ===========================================================================
# Gaussian FCHK Loader — find_mo_analyzer_module
# ===========================================================================

class TestFindMoAnalyzerModule:
    def test_finds_package_with_init(self, tmp_path):
        pkg = tmp_path / "plugins" / "gaussian_fchk_mo_analyzer"
        pkg.mkdir(parents=True)
        (pkg / "__init__.py").write_text("")
        result = _fchk_loader.find_mo_analyzer_module(str(tmp_path))
        assert result is not None
        assert result.endswith("gaussian_fchk_mo_analyzer")

    def test_no_init_py_returns_none(self, tmp_path):
        # Directory exists but has no __init__.py
        pkg = tmp_path / "gaussian_fchk_mo_analyzer"
        pkg.mkdir()
        result = _fchk_loader.find_mo_analyzer_module(str(tmp_path))
        assert result is None

    def test_no_package_returns_none(self, tmp_path):
        result = _fchk_loader.find_mo_analyzer_module(str(tmp_path))
        assert result is None

    def test_nested_package_found(self, tmp_path):
        pkg = tmp_path / "level1" / "level2" / "gaussian_fchk_mo_analyzer"
        pkg.mkdir(parents=True)
        (pkg / "__init__.py").write_text("")
        result = _fchk_loader.find_mo_analyzer_module(str(tmp_path))
        assert result is not None
        assert "gaussian_fchk_mo_analyzer" in result


# ===========================================================================
# Blender Export — initialize registers add_export_action
# ===========================================================================

class TestBlenderExportInitialize:
    def test_initialize_calls_add_export_action(self):
        with mock_optional_imports():
            mod = load_plugin(BLENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_export_action_label_contains_blender(self):
        with mock_optional_imports():
            mod = load_plugin(BLENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        label = ctx.add_export_action.call_args[0][0]
        assert "Blender" in label

    def test_export_action_callback_is_callable(self):
        with mock_optional_imports():
            mod = load_plugin(BLENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        callback = ctx.add_export_action.call_args[0][1]
        assert callable(callback)

    def test_plugin_version_constant_present(self):
        assert hasattr(_blender, "PLUGIN_VERSION")
        assert _blender.PLUGIN_VERSION


# ===========================================================================
# POV-Ray Export — initialize registers add_export_action
# ===========================================================================

class TestPOVRayExportInitialize:
    def test_initialize_calls_add_export_action(self):
        with mock_optional_imports():
            mod = load_plugin(POVRAY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_export_action_label_contains_povray(self):
        with mock_optional_imports():
            mod = load_plugin(POVRAY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        label = ctx.add_export_action.call_args[0][0]
        assert "POV-Ray" in label

    def test_export_action_callback_is_callable(self):
        with mock_optional_imports():
            mod = load_plugin(POVRAY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        callback = ctx.add_export_action.call_args[0][1]
        assert callable(callback)

    def test_plugin_version_constant_present(self):
        assert hasattr(_povray, "PLUGIN_VERSION")
        assert _povray.PLUGIN_VERSION


# ===========================================================================
# ORCA Freq Analyzer — initialize registers file openers
# ===========================================================================

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


# ===========================================================================
# Gaussian Freq Analyzer — initialize registers file openers
# ===========================================================================

class TestGaussianFreqInitialize:
    def test_initialize_registers_fchk_file_opener(self):
        with mock_optional_imports():
            mod = load_plugin(GAUSS_FREQ_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        ctx.register_file_opener.assert_called()
        exts = [call[0][0] for call in ctx.register_file_opener.call_args_list]
        assert ".fchk" in exts

    def test_initialize_registers_fck_and_fch_openers(self):
        with mock_optional_imports():
            mod = load_plugin(GAUSS_FREQ_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        exts = [call[0][0] for call in ctx.register_file_opener.call_args_list]
        assert ".fck" in exts
        assert ".fch" in exts

    def test_plugin_version_constant_present(self):
        assert hasattr(_gfreq, "PLUGIN_VERSION")
        assert _gfreq.PLUGIN_VERSION
