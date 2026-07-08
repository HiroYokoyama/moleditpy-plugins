"""
Tests for the Symmetry Analyzer plugin (initialize -> add_menu_action;
PLUGIN_CONTEXT stored; SymmetryAnalysisWorker.run emits finished with (dict, bool)).
"""

from __future__ import annotations

import ast
import logging
import textwrap
from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
SYMMETRY_PATH = PLUGINS_DIR / "Symmetry_Analyzer" / "symmetry_analyzer.py"


def _extract_method_as_fn(path: Path, class_name: str, method_name: str, extra_globals: dict | None = None):
    """
    Use AST to extract a class method as a standalone callable.

    Needed because Qt base classes are MagicMock instances whose metaclass call
    returns a MagicMock instead of a real type, so the class definition doesn't
    produce a usable type and object.__new__ fails.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == method_name:
                    func_src = ast.get_source_segment(source, item)
                    if func_src:
                        local_ns: dict = {}
                        globs = {"logging": logging, **(extra_globals or {})}
                        exec(textwrap.dedent(func_src), globs, local_ns)
                        return local_ns[method_name]
    return None


# Extract SymmetryAnalysisWorker.run() as a standalone function.
# SymmetryAnalysisWorker inherits QThread (mocked -> MagicMock metaclass),
# so the class itself becomes a MagicMock; object.__new__ cannot be used.
# The run() body uses np.arange (returns a MagicMock, iterable as empty)
# and PointGroupAnalyzer (inside a try/except, never reached when loop is empty).
_sym_run_raw = _extract_method_as_fn(
    SYMMETRY_PATH, "SymmetryAnalysisWorker", "run",
    extra_globals={"np": MagicMock()},  # np.arange returns MagicMock -> iter([])
)


class _FakeWorker:
    """Minimal self-object for the extracted run() — carries only what run() reads."""
    def __init__(self, min_tol, max_tol):
        self.mol_pmg = MagicMock()
        self.min_tol = min_tol
        self.max_tol = max_tol
        self.finished = MagicMock()


class TestSymmetryAnalyzer:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()

    def test_initialize_menu_path_contains_symmetr(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            path = ctx.add_menu_action.call_args[0][0]
            assert "Symmetr" in path

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_worker_run_emits_finished(self):
        """run() must always emit finished — with mocked numpy the tolerance
        loop iterates 0 times (MagicMock __iter__ returns empty iterator)."""
        worker = _FakeWorker(0.1, 0.1)
        with mock_optional_imports():
            _sym_run_raw(worker)
        worker.finished.emit.assert_called_once()

    def test_worker_run_emits_dict_and_bool(self):
        """finished.emit(group_data, found_any): group_data is dict, found_any is bool."""
        worker = _FakeWorker(0.05, 0.2)
        with mock_optional_imports():
            _sym_run_raw(worker)
        args = worker.finished.emit.call_args[0]
        assert len(args) == 2
        group_data, found_any = args
        assert isinstance(group_data, dict)
        assert isinstance(found_any, bool)

    def test_worker_found_any_false_when_no_real_analysis(self):
        """With mocked pymatgen the loop body is never entered → found_any=False."""
        worker = _FakeWorker(0.1, 0.5)
        with mock_optional_imports():
            _sym_run_raw(worker)
        _, found_any = worker.finished.emit.call_args[0]
        assert found_any is False




# ---------------------------------------------------------------------------
# operation classification (numpy stub)
# ---------------------------------------------------------------------------

import math
from types import SimpleNamespace


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
