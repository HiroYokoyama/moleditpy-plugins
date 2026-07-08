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
