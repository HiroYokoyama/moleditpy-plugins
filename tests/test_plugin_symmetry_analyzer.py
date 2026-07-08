"""
Tests for the Symmetry Analyzer plugin (initialize -> add_menu_action;
PLUGIN_CONTEXT stored; SymmetryAnalysisWorker.run emits finished with (dict, bool)).
"""

from __future__ import annotations

import ast
import logging
import textwrap
from pathlib import Path
from unittest import mock
from unittest.mock import MagicMock

import numpy as np

from conftest import (
    FakeMol,
    load_plugin,
    make_context,
    mock_optional_imports,
    mocks_with_real_numpy,
)

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


# ---------------------------------------------------------------------------
# Shared helpers for the extended coverage below
# ---------------------------------------------------------------------------

from types import SimpleNamespace  # noqa: E402  (kept near first usage below)


def _extract_raw(name, extra_globals=None):
    """Extract a SymmetryAnalysisPlugin method with no bound self-behavior."""
    src = _extract_method_source(SYMMETRY_PATH, "SymmetryAnalysisPlugin", name)
    return _make_function(src, dict(extra_globals or {}))


class _FakeListWidget:
    """Minimal QListWidget stand-in: tracks addItem/clear/setCurrentRow calls."""

    def __init__(self):
        self.items: list[str] = []
        self.current_row = None

    def clear(self):
        self.items = []

    def addItem(self, text):
        self.items.append(text)

    def count(self):
        return len(self.items)

    def setCurrentRow(self, row):
        self.current_row = row


class _FakeItem:
    def __init__(self, text):
        self._text = text

    def text(self):
        return self._text


# ---------------------------------------------------------------------------
# get_pymatgen_molecule
# ---------------------------------------------------------------------------


class _RdAtomZ:
    def __init__(self, z):
        self._z = z

    def GetAtomicNum(self):
        return self._z


class _RdMolStub:
    """rdkit-like stand-in exposing only what get_pymatgen_molecule needs."""

    def __init__(self, atomic_nums, coords, no_conformer=False):
        self.atoms = [_RdAtomZ(z) for z in atomic_nums]
        self._conf = FakeMol(atomic_nums, coords).conf
        self.no_conformer = no_conformer

    def GetConformer(self):
        if self.no_conformer:
            raise ValueError("no conformer")
        return self._conf

    def GetAtoms(self):
        return list(self.atoms)

    def GetNumAtoms(self):
        return len(self.atoms)


class TestGetPymatgenMolecule:
    def _fn(self, molecule_cls):
        return _extract_raw("get_pymatgen_molecule", {"Molecule": molecule_cls})

    def test_no_current_molecule_returns_none(self):
        fn = self._fn(MagicMock())
        s = SimpleNamespace(context=SimpleNamespace(current_molecule=None))
        assert fn(s) is None

    def test_no_conformer_returns_none(self):
        fn = self._fn(MagicMock())
        rd_mol = _RdMolStub([8, 1, 1], [(0, 0, 0), (1, 0, 0), (0, 1, 0)], no_conformer=True)
        s = SimpleNamespace(context=SimpleNamespace(current_molecule=rd_mol))
        assert fn(s) is None

    def test_species_and_coords_extracted_and_passed_to_molecule(self):
        captured = {}

        def fake_molecule(species, coords):
            captured["species"] = species
            captured["coords"] = coords
            return "MOL"

        rd_mol = _RdMolStub(
            [8, 1, 1],
            [(0.0, 0.0, 0.3), (0.8, 0.0, -0.3), (-0.8, 0.0, -0.3)],
        )
        s = SimpleNamespace(context=SimpleNamespace(current_molecule=rd_mol))
        fn = self._fn(fake_molecule)
        result = fn(s)
        assert result == "MOL"
        assert captured["species"] == [8, 1, 1]
        assert [list(c) for c in captured["coords"]] == [
            [0.0, 0.0, 0.3],
            [0.8, 0.0, -0.3],
            [-0.8, 0.0, -0.3],
        ]


# ---------------------------------------------------------------------------
# on_analysis_finished — sort/select logic
# ---------------------------------------------------------------------------


class TestOnAnalysisFinished:
    def _self(self):
        return SimpleNamespace(
            calc_btn=MagicMock(),
            groups_list=_FakeListWidget(),
            group_data={},
        )

    def _fn(self):
        return _extract_raw("on_analysis_finished")

    def test_not_found_any_adds_placeholder(self):
        fn = self._fn()
        s = self._self()
        fn(s, {}, False)
        assert s.groups_list.items == ["No point groups found."]
        s.calc_btn.setEnabled.assert_any_call(True)

    def test_sorted_by_min_tol_then_max_tol_desc(self):
        fn = self._fn()
        s = self._self()
        group_data = {
            "C2v": {"tols": [0.2, 0.3]},
            "C1": {"tols": [0.05, 0.9]},
            "D6h": {"tols": [0.05, 0.4]},
        }
        fn(s, group_data, True)
        texts = s.groups_list.items
        # C1 and D6h tie on min(tol)=0.05; sort key is (min_tol, -max_tol), so
        # the entry with the *larger* max tolerance (C1: 0.9 > D6h: 0.4) has
        # the more negative key and sorts first.
        assert texts[0].startswith("C1")
        assert texts[1].startswith("D6h")
        assert texts[2].startswith("C2v")

    def test_selects_first_non_c1_row(self):
        fn = self._fn()
        s = self._self()
        group_data = {
            "C1": {"tols": [0.05, 0.05]},
            "C2v": {"tols": [0.10, 0.20]},
        }
        fn(s, group_data, True)
        # sorted: C1 (min .05) then C2v (min .10) -> row 1 is the first non-C1.
        assert s.groups_list.current_row == 1

    def test_all_c1_selects_row_zero(self):
        fn = self._fn()
        s = self._self()
        fn(s, {"C1": {"tols": [0.1, 0.2]}}, True)
        assert s.groups_list.current_row == 0

    def test_item_text_format(self):
        fn = self._fn()
        s = self._self()
        fn(s, {"Td": {"tols": [0.1, 0.35]}}, True)
        assert s.groups_list.items == ["Td  (Tol: 0.10 - 0.35 Å)"]


# ---------------------------------------------------------------------------
# on_group_selected
# ---------------------------------------------------------------------------


class TestOnGroupSelected:
    def _make_self(self, item_text, group_data):
        fmt_fn = _extract_raw("_format_symmetry_symbol")
        sort_fn = _extract_raw("_get_op_sort_key", {"np": _SymNP()})
        return SimpleNamespace(
            groups_list=SimpleNamespace(currentItem=lambda: _FakeItem(item_text) if item_text else None),
            _format_symmetry_symbol=lambda sym: fmt_fn(None, sym),
            _get_op_sort_key=lambda op: sort_fn(None, op),
            selected_group_label=MagicMock(),
            group_data=group_data,
            analyzer=None,
            symmetry_ops=None,
            update_ops_list=MagicMock(),
            sym_btn=MagicMock(),
            op_details=MagicMock(),
        )

    def _fn(self):
        return _extract_raw("on_group_selected")

    def test_no_current_item_returns(self):
        fn = self._fn()
        s = self._make_self(None, {})
        fn(s)  # should not raise
        s.update_ops_list.assert_not_called()
        s.selected_group_label.setText.assert_not_called()

    def test_symbol_not_in_group_data_sets_label_but_no_ops(self):
        fn = self._fn()
        s = self._make_self("No point groups found.", {})
        fn(s)
        s.selected_group_label.setText.assert_called_once()
        s.update_ops_list.assert_not_called()
        s.sym_btn.setEnabled.assert_not_called()

    def test_group_found_sorts_ops_and_enables_button(self):
        fn = self._fn()
        analyzer = MagicMock()
        analyzer.get_symmetry_operations.return_value = [
            _op(_MAT_S4),
            _op(_MAT_E),
            _op(_MAT_C2),
        ]
        s = self._make_self("Td  (Tol: 0.10 - 1.00 Å)", {"Td": {"analyzer": analyzer}})
        fn(s)
        assert s.analyzer is analyzer
        mats = [op.rotation_matrix for op in s.symmetry_ops]
        assert mats == [_MAT_E, _MAT_C2, _MAT_S4]
        s.update_ops_list.assert_called_once()
        s.sym_btn.setEnabled.assert_called_once_with(True)
        s.op_details.clear.assert_called_once()

    def test_label_uses_html_format_for_html_symbols(self):
        fn = self._fn()
        s = self._make_self("C2v  (Tol: 0.10 - 1.00 Å)", {})
        fn(s)
        text = s.selected_group_label.setText.call_args[0][0]
        assert "Point Group:" in text
        assert "<i>C</i>" in text


# ---------------------------------------------------------------------------
# update_ops_list
# ---------------------------------------------------------------------------


class TestUpdateOpsList:
    def test_populates_list_in_order(self):
        fn = _extract_raw("update_ops_list")
        label_fn = _extract_raw("_get_op_label", {"np": _SymNP()})
        ops = [_op(_MAT_E), _op(_MAT_C2)]
        s = SimpleNamespace(
            symmetry_ops=ops,
            _get_op_label=lambda op, i: label_fn(None, op, i),
            ops_list=_FakeListWidget(),
            op_details=MagicMock(),
        )
        fn(s)
        assert s.ops_list.items == ["#1: E (Identity)", "#2: C₂ (Rotation)"]
        s.op_details.clear.assert_called_once()


# ---------------------------------------------------------------------------
# on_op_selection_changed
# ---------------------------------------------------------------------------


class _FakeSelectedItem:
    def __init__(self, row):
        self.row_val = row


class TestOnOpSelectionChanged:
    def _make_self(self, selected_rows, n_ops):
        items = [_FakeSelectedItem(r) for r in selected_rows]
        ops_list = SimpleNamespace(
            selectedItems=lambda: items,
            row=lambda item: item.row_val,
        )
        return SimpleNamespace(
            ops_list=ops_list,
            symmetry_ops=[_op(_MAT_E) for _ in range(n_ops)],
            visualize_ops=MagicMock(),
            op_details=MagicMock(),
            _display_single_op_details=MagicMock(),
        )

    def _fn(self):
        return _extract_raw("on_op_selection_changed")

    def test_zero_selected_clears_details(self):
        fn = self._fn()
        s = self._make_self([], 3)
        fn(s)
        s.visualize_ops.assert_called_once_with([])
        s.op_details.clear.assert_called_once()
        s._display_single_op_details.assert_not_called()

    def test_one_selected_shows_details(self):
        fn = self._fn()
        s = self._make_self([1], 3)
        fn(s)
        s._display_single_op_details.assert_called_once_with(s.symmetry_ops[1])

    def test_multiple_selected_shows_count(self):
        fn = self._fn()
        s = self._make_self([0, 2], 3)
        fn(s)
        text = s.op_details.setText.call_args[0][0]
        assert "2 operations selected" in text

    def test_out_of_range_row_is_skipped(self):
        fn = self._fn()
        s = self._make_self([5], 3)  # row 5 invalid for n_ops=3
        fn(s)
        s.visualize_ops.assert_called_once_with([])
        s.op_details.clear.assert_called_once()


# ---------------------------------------------------------------------------
# _display_single_op_details (real numpy formatting)
# ---------------------------------------------------------------------------


class TestDisplaySingleOpDetails:
    def _fn(self):
        return _extract_raw("_display_single_op_details", {"np": np})

    def test_sets_text_with_matrix_and_translation(self):
        fn = self._fn()
        op = SimpleNamespace(
            rotation_matrix=np.eye(3),
            translation_vector=np.zeros(3),
            as_xyz_str=lambda: "x, y, z",
        )
        s = SimpleNamespace(op_details=MagicMock())
        fn(s, op)
        text = s.op_details.setText.call_args[0][0]
        assert "Rotation Matrix" in text
        assert "Translation Vector" in text
        assert "Type: x, y, z" in text


# ---------------------------------------------------------------------------
# _add_op_visualization — real geometric math (eigen-decomposition of
# rotation/reflection matrices to find axes/normals)
# ---------------------------------------------------------------------------


class _FakePlotter:
    def __init__(self):
        self.added = []
        self.removed = []

    def add_mesh(self, mesh, **kwargs):
        actor = SimpleNamespace(mesh=mesh, kwargs=kwargs)
        self.added.append(actor)
        return actor

    def remove_actor(self, actor):
        self.removed.append(actor)

    def render(self):
        pass


class _FakePV:
    def Sphere(self, radius, center):
        return SimpleNamespace(kind="sphere", radius=radius, center=center)

    def Line(self, start, end):
        return SimpleNamespace(kind="line", start=start, end=end)

    def Disc(self, center, inner, outer, normal, c_res):
        return SimpleNamespace(kind="disc", center=center, normal=normal)


def _rot_z(deg):
    th = np.radians(deg)
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


class TestAddOpVisualization:
    def _fn(self):
        return _extract_raw("_add_op_visualization", {"np": np})

    def test_identity_adds_nothing(self):
        fn = self._fn()
        plotter = _FakePlotter()
        s = SimpleNamespace(vis_actors=[])
        fn(s, plotter, _op(np.eye(3)), np.zeros(3), _FakePV(), 1.0)
        assert plotter.added == []
        assert s.vis_actors == []

    def test_inversion_adds_orange_sphere(self):
        fn = self._fn()
        plotter = _FakePlotter()
        s = SimpleNamespace(vis_actors=[])
        com = np.array([0.1, 0.2, 0.3])
        fn(s, plotter, _op(-np.eye(3)), com, _FakePV(), 2.0)
        assert len(s.vis_actors) == 1
        actor = plotter.added[0]
        assert actor.kwargs["color"] == "orange"
        assert actor.mesh.kind == "sphere"
        assert np.allclose(actor.mesh.center, com)

    def test_c2_rotation_adds_line_along_rotation_axis(self):
        fn = self._fn()
        plotter = _FakePlotter()
        s = SimpleNamespace(vis_actors=[])
        com = np.zeros(3)
        fn(s, plotter, _op(_rot_z(180)), com, _FakePV(), 1.0)
        assert len(s.vis_actors) == 1
        line = plotter.added[0].mesh
        assert line.kind == "line"
        # C2 about z: the drawn axis line must run parallel to z (x=y=0).
        assert abs(line.start[0]) < 1e-8 and abs(line.start[1]) < 1e-8
        assert abs(line.end[0]) < 1e-8 and abs(line.end[1]) < 1e-8

    def test_mirror_reflection_adds_disc_with_correct_normal(self):
        fn = self._fn()
        plotter = _FakePlotter()
        s = SimpleNamespace(vis_actors=[])
        com = np.zeros(3)
        # sigma_xy mirror plane: reflect z -> -z, normal is +/-z.
        m = np.diag([1.0, 1.0, -1.0])
        fn(s, plotter, _op(m), com, _FakePV(), 1.0)
        assert len(s.vis_actors) == 1
        disc = plotter.added[0].mesh
        assert disc.kind == "disc"
        normal = np.array(disc.normal)
        assert abs(abs(normal[2]) - 1.0) < 1e-6
        assert abs(normal[0]) < 1e-6 and abs(normal[1]) < 1e-6


# ---------------------------------------------------------------------------
# visualize_ops — COM/radius computation + dispatch to _add_op_visualization
# ---------------------------------------------------------------------------


class TestVisualizeOps:
    def _fn(self):
        return _extract_raw("visualize_ops")

    def test_no_plotter_returns_immediately(self):
        fn = self._fn()
        s = SimpleNamespace(context=SimpleNamespace(plotter=None))
        with mocks_with_real_numpy():
            fn(s, [_op(np.eye(3))])  # should not raise

    def test_empty_ops_clears_actors_and_renders(self):
        fn = self._fn()
        plotter = _FakePlotter()
        actor_stub = object()
        s = SimpleNamespace(
            context=SimpleNamespace(plotter=plotter, current_molecule=None),
            vis_actors=[actor_stub],
        )
        with mocks_with_real_numpy():
            fn(s, [])
        assert s.vis_actors == []
        assert plotter.removed == [actor_stub]

    def test_dispatches_to_add_op_visualization_per_op_with_com(self):
        fn = self._fn()
        plotter = _FakePlotter()
        rd_mol = FakeMol(["O", "H", "H"], [(0, 0, 1), (1, 0, -1), (-1, 0, -1)])
        calls = []

        def fake_add(plotter_, op, com, pv, scale):
            calls.append((op, tuple(com), scale))

        s = SimpleNamespace(
            context=SimpleNamespace(plotter=plotter, current_molecule=rd_mol),
            vis_actors=[],
            _add_op_visualization=fake_add,
        )
        ops = [_op(np.eye(3)), _op(-np.eye(3))]
        with mocks_with_real_numpy():
            fn(s, ops)
        assert len(calls) == 2
        # COM of the three points (mean) -> (0, 0, -1/3)
        com = calls[0][1]
        assert abs(com[0]) < 1e-8 and abs(com[1]) < 1e-8
        assert abs(com[2] - (-1.0 / 3.0)) < 1e-8


# ---------------------------------------------------------------------------
# symmetrize_structure — real vector math over synthetic water / chiral
# geometries, exact and perturbed.
# ---------------------------------------------------------------------------


class _SymOp:
    """pymatgen SymmOp stand-in: pure rotation about an already-centered origin."""

    def __init__(self, matrix):
        self.matrix = np.array(matrix, dtype=float)

    def operate(self, p):
        return self.matrix @ np.asarray(p, dtype=float)


class _Species:
    def __init__(self, symbol):
        self.symbol = symbol


def _refl(axis):
    d = [1.0, 1.0, 1.0]
    d[axis] = -1.0
    return np.diag(d)


# Water-like C2v geometry: C2 about z, mirror planes xz (y->-y) and yz (x->-x).
_WATER_COORDS = [(0.0, 0.0, 0.3), (0.8, 0.0, -0.3), (-0.8, 0.0, -0.3)]
_WATER_SPECIES = ["O", "H", "H"]
_WATER_OPS = [
    _SymOp(np.eye(3)),
    _SymOp(_rot_z(180)),
    _SymOp(_refl(1)),
    _SymOp(_refl(0)),
]


class TestSymmetrizeStructure:
    def _fn(self):
        return _extract_raw("symmetrize_structure", {"np": np, "QMessageBox": MagicMock()})

    def _make_self(self, coords, species, ops, sch_symbol="C2v"):
        mol_pmg = SimpleNamespace(
            cart_coords=np.array(coords, dtype=float),
            species=[_Species(s) for s in species],
        )
        analyzer = SimpleNamespace(
            sch_symbol=sch_symbol,
            get_symmetry_operations=lambda: ops,
        )
        return SimpleNamespace(
            analyzer=analyzer,
            get_pymatgen_molecule=lambda: mol_pmg,
            update_rdkit_coords=MagicMock(),
        )

    def test_no_analyzer_returns_immediately(self):
        fn = self._fn()
        s = SimpleNamespace(analyzer=None, get_pymatgen_molecule=MagicMock())
        fn(s)
        s.get_pymatgen_molecule.assert_not_called()

    def test_no_molecule_returns(self):
        fn = self._fn()
        s = SimpleNamespace(
            analyzer=SimpleNamespace(get_symmetry_operations=MagicMock()),
            get_pymatgen_molecule=lambda: None,
            update_rdkit_coords=MagicMock(),
        )
        fn(s)
        s.update_rdkit_coords.assert_not_called()

    def test_no_symmetry_ops_returns(self):
        fn = self._fn()
        s = self._make_self(_WATER_COORDS, _WATER_SPECIES, [])
        fn(s)
        s.update_rdkit_coords.assert_not_called()

    def test_exact_symmetric_geometry_is_reproduced(self):
        fn = self._fn()
        s = self._make_self(_WATER_COORDS, _WATER_SPECIES, _WATER_OPS)
        with mock.patch.dict("sys.modules", {"scipy": None}):
            fn(s)
        s.update_rdkit_coords.assert_called_once()
        final_coords = s.update_rdkit_coords.call_args[0][0]
        assert np.allclose(final_coords, np.array(_WATER_COORDS), atol=1e-8)

    def test_slightly_perturbed_water_still_symmetrizes_cleanly(self):
        fn = self._fn()
        perturbed = [
            (0.0, 0.0, 0.3),
            (0.82, 0.01, -0.29),  # H1 nudged slightly off-symmetric
            (-0.79, -0.01, -0.31),
        ]
        s = self._make_self(perturbed, _WATER_SPECIES, _WATER_OPS)
        with mock.patch.dict("sys.modules", {"scipy": None}):
            fn(s)
        s.update_rdkit_coords.assert_called_once()
        final_coords = s.update_rdkit_coords.call_args[0][0]
        # Averaging over the (still-applicable) C2v ops should pull the
        # structure back close to the exact symmetric geometry.
        assert np.allclose(final_coords, np.array(_WATER_COORDS), atol=0.05)
        # No mapping-error warning should have fired for a small perturbation.
        s.update_rdkit_coords.assert_called_once()

    def test_mismatched_symmetry_op_on_chiral_geometry_triggers_warning(self):
        fn = self._fn()
        chiral_coords = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.2),
            (0.0, 1.0, 0.5),
            (-1.0, -0.3, 0.7),
        ]
        chiral_species = ["C", "H", "H", "H"]
        # A C2(z) op is not a real symmetry of this chiral geometry -> the
        # rotated positions won't line up with any real atom.
        bad_ops = [_SymOp(np.eye(3)), _SymOp(_rot_z(180))]
        s = self._make_self(chiral_coords, chiral_species, bad_ops, sch_symbol="C1")
        qmb = MagicMock()
        fn2 = _extract_raw("symmetrize_structure", {"np": np, "QMessageBox": qmb})
        with mock.patch.dict("sys.modules", {"scipy": None}):
            fn2(s)
        qmb.warning.assert_called_once()


# ---------------------------------------------------------------------------
# update_rdkit_coords
# ---------------------------------------------------------------------------


class _Point3DStub:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class _FakeConfSettable:
    def __init__(self, n):
        self.positions = [None] * n

    def SetAtomPosition(self, i, p):
        self.positions[i] = p


class _RdMolSettable:
    def __init__(self, n):
        self.conf = _FakeConfSettable(n)

    def GetNumAtoms(self):
        return len(self.conf.positions)

    def GetConformer(self):
        return self.conf


class TestUpdateRdkitCoords:
    def _fn(self):
        return _extract_raw("update_rdkit_coords", {"Point3D": _Point3DStub})

    def test_no_molecule_returns(self):
        fn = self._fn()
        ctx = SimpleNamespace(current_molecule=None, push_undo_checkpoint=MagicMock())
        s = SimpleNamespace(context=ctx)
        fn(s, np.zeros((0, 3)))
        ctx.push_undo_checkpoint.assert_not_called()

    def test_updates_conformer_positions_and_refreshes_view(self):
        fn = self._fn()
        rd_mol = _RdMolSettable(2)
        ctx = SimpleNamespace(
            current_molecule=rd_mol,
            push_undo_checkpoint=MagicMock(),
            refresh_3d_view=MagicMock(),
        )
        s = SimpleNamespace(context=ctx)
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        fn(s, coords)
        ctx.push_undo_checkpoint.assert_called_once()
        assert rd_mol.conf.positions[0].x == 1.0
        assert rd_mol.conf.positions[1].z == 6.0
        assert ctx.current_mol is rd_mol
        ctx.refresh_3d_view.assert_called_once()


# ---------------------------------------------------------------------------
# closeEvent
# ---------------------------------------------------------------------------


class TestCloseEvent:
    def _fn(self):
        src = _extract_method_source(SYMMETRY_PATH, "SymmetryAnalysisPlugin", "closeEvent")
        # super().closeEvent(event) requires a real __class__ closure cell that
        # only exists when compiled inside an actual class body; strip it for
        # standalone extraction (QDialog.closeEvent is a mocked no-op anyway).
        src = src.replace("super().closeEvent(event)", "pass  # super() stripped for standalone extraction")
        return _make_function(src, {})

    def _make_self(self, worker=None, plotter=None):
        return SimpleNamespace(
            worker=worker,
            context=SimpleNamespace(plotter=plotter),
            vis_actors=["a1", "a2"] if plotter else None,
            groups_list=MagicMock(),
            ops_list=MagicMock(),
            selected_group_label=MagicMock(),
            op_details=MagicMock(),
            sym_btn=MagicMock(),
            group_data={"x": 1},
            symmetry_ops=[1, 2],
            analyzer=object(),
        )

    def test_no_worker_no_plotter_resets_state(self):
        fn = self._fn()
        s = self._make_self()
        fn(s, MagicMock())
        assert s.group_data == {}
        assert s.symmetry_ops == []
        assert s.analyzer is None
        s.sym_btn.setEnabled.assert_called_once_with(False)

    def test_running_worker_is_stopped_cleanly(self):
        fn = self._fn()
        worker = MagicMock()
        worker.isRunning.side_effect = [True, False]
        s = self._make_self(worker=worker)
        fn(s, MagicMock())
        worker.quit.assert_called_once()
        worker.wait.assert_called_once_with(1000)
        worker.terminate.assert_not_called()

    def test_worker_still_running_after_wait_is_terminated(self):
        fn = self._fn()
        worker = MagicMock()
        worker.isRunning.return_value = True
        s = self._make_self(worker=worker)
        fn(s, MagicMock())
        worker.terminate.assert_called_once()

    def test_plotter_actors_are_removed(self):
        fn = self._fn()
        plotter = _FakePlotter()
        s = self._make_self(plotter=plotter)
        fn(s, MagicMock())
        assert plotter.removed == ["a1", "a2"]
        assert s.vis_actors == []


# ---------------------------------------------------------------------------
# initialize()'s toggle_window closure and the legacy run(mw) entry point
# ---------------------------------------------------------------------------


class TestInitializeToggleWindow:
    def test_toggle_shows_existing_window(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            callback = ctx.add_menu_action.call_args[0][1]
            existing = MagicMock()
            ctx.get_window.return_value = existing
            callback()
            existing.show.assert_called_once()
            existing.raise_.assert_called_once()
            existing.activateWindow.assert_called_once()

    def test_toggle_creates_new_window_when_absent(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            callback = ctx.add_menu_action.call_args[0][1]
            ctx.get_window.return_value = None
            callback()  # should not raise


class TestRunEntryPoint:
    def test_run_without_plugin_manager_returns(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            mw = SimpleNamespace()
            mod.run(mw)  # no plugin_manager attr -> early return, no crash

    def test_run_without_context_returns(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            mod.PLUGIN_CONTEXT = None
            mw = SimpleNamespace(plugin_manager=MagicMock())
            mod.run(mw)  # PLUGIN_CONTEXT falsy -> early return, no crash

    def test_run_shows_existing_window(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            existing = MagicMock()
            ctx.get_window.return_value = existing
            mw = SimpleNamespace(plugin_manager=MagicMock())
            mod.run(mw)
            existing.show.assert_called_once()
            existing.raise_.assert_called_once()
            existing.activateWindow.assert_called_once()

    def test_run_creates_new_window(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.get_window.return_value = None
            mw = SimpleNamespace(plugin_manager=MagicMock())
            mod.run(mw)  # should not raise
