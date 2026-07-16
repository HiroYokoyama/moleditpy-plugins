"""
Tests for the Gaussian MO Analyzer plugin (plugins/Gaussian_MO_Analyzer/
gaussian_fchk_mo_analyzer/, a folder plugin): FCHKReader, normalization
prefactor, cube writer, and the package initialize() contract.

Methods that live on Qt-derived classes (mocked bases) are extracted via
``ast.get_source_segment`` + ``exec`` into an isolated namespace with a
small pure-Python numpy stub.
"""

from __future__ import annotations

import ast
import importlib.util
import math
import sys
import textwrap
import types
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from conftest import load_plugin, make_context, mock_optional_imports, mocks_with_real_numpy

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
MO_PKG_DIR = PLUGINS_DIR / "Gaussian_MO_Analyzer" / "gaussian_fchk_mo_analyzer"
MO_ANALYZER_PATH = MO_PKG_DIR / "analyzer.py"

# analyzer.py has no PyQt6/rdkit imports at module scope (rdkit is imported
# lazily inside get_xyz_block/get_basis_labels), so it can be loaded as a
# real module with real numpy for genuine coverage of BasisSetEngine, unlike
# the AST-extracted classes above which stand in for methods on Qt-derived
# classes elsewhere in this repo.
with mocks_with_real_numpy():
    _mo_real = load_plugin(MO_ANALYZER_PATH)
RealFCHKReader = _mo_real.FCHKReader
RealBasisSetEngine = _mo_real.BasisSetEngine


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


VIS_PATH = MO_PKG_DIR / "vis.py"


class TestMOCubeParserNegativeAtoms:
    """CubeVisualizer._parse_cube (vis.py) must skip the DSET_IDS block that
    a Gaussian MO cube (negative atom count) inserts between the atom lines
    and the volumetric data; otherwise the data is shifted / misreshaped."""

    _HEADER = (
        "MO cube\n"
        "test\n"
        "   -1    0.0 0.0 0.0\n"
        "    2    1.0 0.0 0.0\n"
        "    2    0.0 1.0 0.0\n"
        "    2    0.0 0.0 1.0\n"
        "    1    1.0    0.0 0.0 0.0\n"
    )
    _DATA = "1.0 2.0 3.0 4.0\n5.0 6.0 7.0 8.0\n"

    def _parse(self, tmp_path, text):
        src = _extract_method_source(VIS_PATH, "CubeVisualizer", "_parse_cube")
        fn = _make_function(src, {"np": np})
        path = tmp_path / "mo.cube"
        path.write_text(text, encoding="utf-8")
        return fn(types.SimpleNamespace(), str(path))

    def test_positive_atom_count_reads_all_data(self, tmp_path):
        header = (
            "cube\ntest\n"
            "    1    0.0 0.0 0.0\n"
            "    2    1.0 0.0 0.0\n"
            "    2    0.0 1.0 0.0\n"
            "    2    0.0 0.0 1.0\n"
            "    8    8.0    0.0 0.0 0.0\n"
        )
        meta = self._parse(tmp_path, header + self._DATA)
        assert meta["dims"] == (2, 2, 2)
        assert list(meta["data"]) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

    def test_dset_ids_line_is_skipped(self, tmp_path):
        meta = self._parse(tmp_path, self._HEADER + "    1    5\n" + self._DATA)
        assert list(meta["data"]) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

    def test_wrapped_dset_ids_block_is_skipped(self, tmp_path):
        ids = "    3    5 6\n    7\n"  # 3 ids wrapped over two lines
        meta = self._parse(tmp_path, self._HEADER + ids + self._DATA)
        assert list(meta["data"]) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]


def _s_shell_fchk(exponent=1.0, coeff=1.0):
    return textwrap.dedent(f"""\
        Atomic numbers                              I   N=           1
                 1
        Current cartesian coordinates               R   N=           3
          0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
        Shell types                                 I   N=           1
                 0
        Number of primitives per shell              I   N=           1
                 1
        Shell to atom map                           I   N=           1
                 1
        Primitive exponents                         R   N=           1
          {exponent:.12E}
        Contraction coefficients                    R   N=           1
          {coeff:.12E}
        """)


_P_SHELL_FCHK = textwrap.dedent("""\
    Atomic numbers                              I   N=           1
             1
    Current cartesian coordinates               R   N=           3
      0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
    Shell types                                 I   N=           1
             1
    Number of primitives per shell              I   N=           1
             1
    Shell to atom map                           I   N=           1
             1
    Primitive exponents                         R   N=           1
      5.000000000000E-01
    Contraction coefficients                    R   N=           1
      1.000000000000E+00
    """)

_SP_SHELL_FALLBACK_FCHK = textwrap.dedent("""\
    Atomic numbers                              I   N=           1
             1
    Current cartesian coordinates               R   N=           3
      0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
    Shell types                                 I   N=           1
            -1
    Number of primitives per shell              I   N=           1
             1
    Shell to atom map                           I   N=           1
             1
    Primitive exponents                         R   N=           1
      5.000000000000E-01
    Contraction coefficients                    R   N=           2
      1.000000000000E+00  5.000000000000E-01
    """)

_SP_SHELL_PMOD_FCHK = textwrap.dedent("""\
    Atomic numbers                              I   N=           1
             1
    Current cartesian coordinates               R   N=           3
      0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
    Shell types                                 I   N=           1
            -1
    Number of primitives per shell              I   N=           1
             1
    Shell to atom map                           I   N=           1
             1
    Primitive exponents                         R   N=           1
      5.000000000000E-01
    Contraction coefficients                    R   N=           1
      1.000000000000E+00
    P(S=P) Contraction coefficients             R   N=           1
      3.000000000000E-01
    """)

_UNSUPPORTED_THEN_S_FCHK = textwrap.dedent("""\
    Atomic numbers                              I   N=           1
             1
    Current cartesian coordinates               R   N=           3
      0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
    Shell types                                 I   N=           2
             9           0
    Number of primitives per shell              I   N=           2
             1           1
    Shell to atom map                           I   N=           2
             1           1
    Primitive exponents                         R   N=           2
      1.000000000000E+00  2.000000000000E+00
    Contraction coefficients                    R   N=           2
      1.000000000000E+00  1.000000000000E+00
    """)


def _d_shell_fchk(pure_flag=None):
    extra = ""
    if pure_flag is not None:
        extra = f"Pure/Cartesian d shells                     I                {pure_flag}\n"
    return textwrap.dedent("""\
        Atomic numbers                              I   N=           1
                 1
        Current cartesian coordinates               R   N=           3
          0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
        Shell types                                 I   N=           1
                 2
        Number of primitives per shell              I   N=           1
                 1
        Shell to atom map                           I   N=           1
                 1
        Primitive exponents                         R   N=           1
          1.000000000000E+00
        Contraction coefficients                    R   N=           1
          1.000000000000E+00
        """) + extra


def _neg_d_shell_fchk():
    return textwrap.dedent("""\
        Atomic numbers                              I   N=           1
                 1
        Current cartesian coordinates               R   N=           3
          0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
        Shell types                                 I   N=           1
                -2
        Number of primitives per shell              I   N=           1
                 1
        Shell to atom map                           I   N=           1
                 1
        Primitive exponents                         R   N=           1
          1.000000000000E+00
        Contraction coefficients                    R   N=           1
          1.000000000000E+00
        Pure/Cartesian d shells                     I                0
        """)


def _make_engine(tmp_path, content, name="mo.fchk"):
    f = tmp_path / name
    f.write_text(content)
    reader = RealFCHKReader(str(f))
    return RealBasisSetEngine(reader), reader


class TestMOBasisSetEnginePrepareBasisSet:
    def test_s_shell_single_basis_function(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _s_shell_fchk())
        assert engine.n_basis == 1
        assert len(engine.shells) == 1

    def test_p_shell_three_basis_functions(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _P_SHELL_FCHK)
        assert engine.n_basis == 3
        assert len(engine.shells) == 1

    def test_sp_shell_fallback_coeffs_four_basis_functions(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _SP_SHELL_FALLBACK_FCHK)
        # SP shell expands into an S sub-shell (1 func) + P sub-shell (3 funcs)
        assert engine.n_basis == 4
        assert len(engine.shells) == 2
        assert engine.shells[0]["type"] == 0
        assert engine.shells[1]["type"] == 1

    def test_sp_shell_p_modifiers_path(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _SP_SHELL_PMOD_FCHK)
        assert engine.n_basis == 4
        assert len(engine.shells) == 2

    def test_unsupported_shell_type_skipped_pointers_advance(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _UNSUPPORTED_THEN_S_FCHK)
        # Only the trailing S shell (type 0) should survive; the unsupported
        # type-9 shell is skipped but its primitive/coeff slots are consumed.
        assert engine.n_basis == 1
        assert len(engine.shells) == 1
        assert engine.shells[0]["type"] == 0

    def test_d_shell_defaults_to_cartesian_six_functions(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _d_shell_fchk(pure_flag=None))
        assert engine.n_basis == 6

    def test_d_shell_cartesian_flag_one_six_functions(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _d_shell_fchk(pure_flag=1))
        assert engine.n_basis == 6

    def test_d_shell_pure_flag_zero_five_functions(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _d_shell_fchk(pure_flag=0))
        assert engine.n_basis == 5

    def test_negative_shell_type_uses_abs_for_definitions(self, tmp_path):
        # Shell type -2 with the pure/spherical D flag set: effective_type
        # abs(-2) == 2, same 5-function spherical definition as type 2 pure.
        engine, _ = _make_engine(tmp_path, _neg_d_shell_fchk())
        assert engine.n_basis == 5


class TestMOBasisSetEngineEvaluateMoOnGrid:
    def test_s_function_value_at_center(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _s_shell_fchk(exponent=1.0, coeff=1.0))
        grid = np.array([[0.0, 0.0, 0.0]])
        coeffs = np.array([1.0])
        phi = engine.evaluate_mo_on_grid(0, grid, coeffs)
        expected = (2 * 1.0 / np.pi) ** 0.75
        assert phi[0] == pytest.approx(expected)

    def test_s_function_decays_with_distance(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _s_shell_fchk(exponent=1.0, coeff=1.0))
        grid = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        coeffs = np.array([1.0])
        phi = engine.evaluate_mo_on_grid(0, grid, coeffs)
        assert phi[1] < phi[0]

    def test_small_coefficient_below_threshold_contributes_zero(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _s_shell_fchk())
        grid = np.array([[0.0, 0.0, 0.0]])
        coeffs = np.array([1e-10])
        phi = engine.evaluate_mo_on_grid(0, grid, coeffs)
        assert phi[0] == pytest.approx(0.0)

    def test_out_of_range_mo_index_raises_value_error(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _s_shell_fchk())
        grid = np.array([[0.0, 0.0, 0.0]])
        coeffs = np.array([1.0])
        with pytest.raises(ValueError):
            engine.evaluate_mo_on_grid(1, grid, coeffs)  # only 1 basis fn -> mo 1 out of range


class TestMOBasisSetEngineEvaluateBasisOnGrid:
    def test_single_basis_function_matches_manual_gaussian(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _s_shell_fchk(exponent=1.0, coeff=1.0))
        grid = np.array([[1.0, 0.0, 0.0]])
        phi = engine.evaluate_basis_on_grid(0, grid)
        prefactor = (2 * 1.0 / np.pi) ** 0.75
        expected = prefactor * math.exp(-1.0 * 1.0**2)
        assert phi[0] == pytest.approx(expected)

    def test_out_of_range_basis_index_returns_zeros(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _s_shell_fchk())
        grid = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        phi = engine.evaluate_basis_on_grid(5, grid)
        assert list(phi) == [0.0, 0.0]


class TestMOBasisSetEngineGetBasisLabels:
    def test_p_shell_labels(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _P_SHELL_FCHK)
        with mock_optional_imports():
            labels = engine.get_basis_labels()
        assert len(labels) == 3
        assert labels[0].endswith("Px") or "L1_0" in labels[0]

    def test_negative_shell_type_label_count_matches_n_basis(self, tmp_path):
        # Regression test: get_basis_labels() previously produced zero
        # labels for shell types below -1 (e.g. -2 for pure/spherical D),
        # even though _prepare_basis_set correctly builds 5 basis functions
        # for such a shell via abs(stype). This mismatch (0 labels vs
        # n_basis == 5) breaks any caller that zips labels with MO
        # coefficients (e.g. a GUI table).
        engine, _ = _make_engine(tmp_path, _neg_d_shell_fchk())
        with mock_optional_imports():
            labels = engine.get_basis_labels()
        assert len(labels) == engine.n_basis == 5

    def test_sp_shell_labels_four_components(self, tmp_path):
        engine, _ = _make_engine(tmp_path, _SP_SHELL_FALLBACK_FCHK)
        with mock_optional_imports():
            labels = engine.get_basis_labels()
        assert len(labels) == 4


class _FakePeriodicTable:
    _MAP = {8: "O", 1: "H"}

    def GetElementSymbol(self, z):
        return self._MAP.get(int(z), f"Z{z}")


def _install_fake_rdkit():
    fake_chem = types.ModuleType("rdkit.Chem")
    fake_chem.GetPeriodicTable = lambda: _FakePeriodicTable()
    fake_rdkit = types.ModuleType("rdkit")
    fake_rdkit.Chem = fake_chem
    saved = {k: sys.modules.get(k) for k in ("rdkit", "rdkit.Chem")}
    sys.modules["rdkit"] = fake_rdkit
    sys.modules["rdkit.Chem"] = fake_chem
    return saved


def _restore_rdkit(saved):
    for k, v in saved.items():
        if v is None:
            sys.modules.pop(k, None)
        else:
            sys.modules[k] = v


class TestMOFCHKReaderGetXyzBlock:
    def test_returns_none_when_atoms_missing(self, tmp_path):
        f = tmp_path / "empty.fchk"
        f.write_text("Charge                                      I                0\n")
        reader = RealFCHKReader(str(f))
        assert reader.get_xyz_block() is None

    def test_returns_none_when_coords_missing(self, tmp_path):
        f = tmp_path / "atoms_only.fchk"
        f.write_text(
            "Atomic numbers                              I   N=           1\n"
            "         1\n"
        )
        reader = RealFCHKReader(str(f))
        assert reader.get_xyz_block() is None

    def test_header_and_bohr_to_angstrom_conversion(self, tmp_path):
        f = tmp_path / "water.fchk"
        f.write_text(
            "Atomic numbers                              I   N=           2\n"
            "         8           1\n"
            "Current cartesian coordinates               R   N=           6\n"
            "  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00\n"
            "  1.000000000000E+00  0.000000000000E+00  0.000000000000E+00\n"
        )
        reader = RealFCHKReader(str(f))
        saved = _install_fake_rdkit()
        try:
            block = reader.get_xyz_block()
        finally:
            _restore_rdkit(saved)
        lines = block.splitlines()
        assert lines[0] == "2"
        assert lines[1] == "Generated from FCHK"
        assert lines[2].split()[0] == "O"
        h_parts = lines[3].split()
        assert h_parts[0] == "H"
        assert float(h_parts[1]) == pytest.approx(0.529177249)


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
