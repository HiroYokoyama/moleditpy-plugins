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
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
MO_PKG_DIR = PLUGINS_DIR / "Gaussian_MO_Analyzer" / "gaussian_fchk_mo_analyzer"
MO_ANALYZER_PATH = MO_PKG_DIR / "analyzer.py"


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
