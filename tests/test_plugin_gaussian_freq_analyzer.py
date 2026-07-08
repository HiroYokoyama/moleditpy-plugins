"""
Tests for the Gaussian Freq Analyzer plugin (plugins/Gaussian_Freq_Analyzer/gaussian_fchk_freq_analyzer.py).

All heavy deps (PyQt6, rdkit, numpy) are stubbed via mock_optional_imports().
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
GAUSS_FREQ_PATH = PLUGINS_DIR / "Gaussian_Freq_Analyzer" / "gaussian_fchk_freq_analyzer.py"

with mock_optional_imports():
    _gfreq = load_plugin(GAUSS_FREQ_PATH)


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
        # The module was loaded with rdkit mocked (Chem = MagicMock, truthy).
        # The lazy `from rdkit.Chem import GetPeriodicTable` inside parse()
        # fires after the mock context exits, hitting the real absent rdkit.
        # Setting Chem = None makes the `if Chem:` guard skip that branch.
        _gfreq.Chem = None
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




# ---------------------------------------------------------------------------
# FCHKParser extended edge cases
# ---------------------------------------------------------------------------


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
