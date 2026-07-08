"""
Tests for the XYZ Editor plugin (initialize -> add_menu_action + save/load/reset
handlers; save handler returns {} when no molecule).
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
XYZ_EDITOR_PATH = PLUGINS_DIR / "XYZ_Editor" / "xyz_editor.py"


class TestXYZEditor:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()

    def test_initialize_menu_path_contains_xyz(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            path = ctx.add_menu_action.call_args[0][0]
            assert "XYZ" in path

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_initialize_registers_save_handler(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_save_handler.assert_called_once()

    def test_initialize_registers_load_handler(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_load_handler.assert_called_once()

    def test_initialize_registers_document_reset_handler(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_document_reset_handler.assert_called_once()

    def test_save_handler_returns_empty_dict_when_no_molecule(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            ctx.current_molecule = None
            mod.initialize(ctx)
            save_fn = ctx.register_save_handler.call_args[0][0]
            result = save_fn()
            assert result == {} or result is None

    def test_load_handler_tolerates_empty_dict(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            ctx.current_molecule = None
            mod.initialize(ctx)
            load_fn = ctx.register_load_handler.call_args[0][0]
            load_fn({})




# ---------------------------------------------------------------------------
# serialization, mol signature, persistence handlers
# ---------------------------------------------------------------------------

import ast
import logging
import os
import json
import textwrap
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as real_numpy


def _extract_method_as_fn(
    path: Path,
    class_name: str | None,
    method_name: str,
    extra_globals: dict | None = None,
):
    """
    Use AST to extract a method (or module-level function) as a standalone
    callable, exec'd into a controlled namespace.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    scope = tree
    if class_name is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == class_name:
                scope = node
                break
        else:
            raise AssertionError(f"class {class_name} not found in {path}")
    for node in scope.body:
        if (
            isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name == method_name
        ):
            func_src = ast.get_source_segment(source, node)
            assert func_src is not None
            globs = {"logging": logging, "os": os, "json": json}
            globs.update(extra_globals or {})
            local_ns: dict = {}
            exec(textwrap.dedent(func_src), globs, local_ns)
            return local_ns[method_name]
    raise AssertionError(f"{class_name}.{method_name} not found in {path}")


class FakePos:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class FakeAtom:
    def __init__(self, symbol, num=6, props=None, radicals=0, idx=0):
        self._symbol = symbol
        self._num = num
        self._props = props or {}
        self._radicals = radicals
        self._idx = idx

    def GetSymbol(self):
        return self._symbol

    def GetAtomicNum(self):
        return self._num

    def GetNumRadicalElectrons(self):
        return self._radicals

    def GetIdx(self):
        return self._idx

    def HasProp(self, key):
        return key in self._props

    def GetProp(self, key):
        return self._props[key]

    def SetProp(self, key, val):
        self._props[key] = val


class FakeConf:
    def __init__(self, coords):
        self._coords = coords

    def GetAtomPosition(self, i):
        return FakePos(*self._coords[i])

    def GetPositions(self):
        return real_numpy.array(self._coords, dtype=float)


class FakeMol:
    def __init__(self, atoms, coords, bonds=0):
        self._atoms = atoms
        self._coords = coords
        self._bonds = bonds
        for i, a in enumerate(self._atoms):
            a._idx = i

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetNumBonds(self):
        return self._bonds

    def GetAtoms(self):
        return list(self._atoms)

    def GetAtomWithIdx(self, i):
        return self._atoms[i]

    def GetConformer(self):
        return FakeConf(self._coords)


class _FakeItem:
    def __init__(self, text):
        self._text = text

    def text(self):
        return self._text


class _FakeTable:
    def __init__(self, rows):
        # rows: list of (idx, symbol, x, y, z) as strings
        self._rows = rows

    def rowCount(self):
        return len(self._rows)

    def item(self, row, col):
        return _FakeItem(str(self._rows[row][col]))


_xyz_generate = _extract_method_as_fn(
    XYZ_EDITOR_PATH, "XYZEditorWindow", "_generate_xyz_content"
)


class TestXYZGenerateContent:
    def test_header_lines(self):
        fake = SimpleNamespace(
            table=_FakeTable([("0", "C", "0.00000", "0.00000", "0.00000")])
        )
        lines = _xyz_generate(fake)
        assert lines[0] == "1"
        assert "MoleditPy" in lines[1]

    def test_empty_table(self):
        lines = _xyz_generate(SimpleNamespace(table=_FakeTable([])))
        assert lines == ["0", "Generated by MoleditPy XYZ Editor"]

    def test_atom_line_formatting(self):
        fake = SimpleNamespace(
            table=_FakeTable([("0", "Cl", "1.50000", "-2.00000", "0.12500")])
        )
        line = _xyz_generate(fake)[2]
        assert line.startswith("Cl  ")
        assert "1.50000" in line and "-2.00000" in line

    def test_symbol_whitespace_stripped(self):
        fake = SimpleNamespace(
            table=_FakeTable([("0", "  N ", "0.0", "0.0", "0.0")])
        )
        assert _xyz_generate(fake)[2].startswith("N ")

    def test_row_count_matches(self):
        rows = [(str(i), "C", "0.0", "0.0", "0.0") for i in range(5)]
        lines = _xyz_generate(SimpleNamespace(table=_FakeTable(rows)))
        assert lines[0] == "5"
        assert len(lines) == 7


_xyz_signature = _extract_method_as_fn(
    XYZ_EDITOR_PATH,
    "XYZEditorWindow",
    "get_mol_signature",
    extra_globals={"np": real_numpy},
)


class TestXYZMolSignature:
    def test_none_mol(self):
        assert _xyz_signature(SimpleNamespace(), None) is None

    def test_same_mol_same_signature(self):
        mol = FakeMol([FakeAtom("O", 8)], [(0.0, 0.0, 0.117)], bonds=0)
        s1 = _xyz_signature(SimpleNamespace(), mol)
        s2 = _xyz_signature(SimpleNamespace(), mol)
        assert s1 == s2 and s1 is not None

    def test_coordinate_change_changes_signature(self):
        atoms = [FakeAtom("O", 8)]
        m1 = FakeMol(atoms, [(0.0, 0.0, 0.0)])
        s1 = _xyz_signature(SimpleNamespace(), m1)
        m1._coords = [(0.5, 0.0, 0.0)]
        s2 = _xyz_signature(SimpleNamespace(), m1)
        assert s1 != s2

    def test_sub_rounding_change_ignored(self):
        atoms = [FakeAtom("O", 8)]
        m1 = FakeMol(atoms, [(0.00001, 0.0, 0.0)])
        s1 = _xyz_signature(SimpleNamespace(), m1)
        m1._coords = [(0.00002, 0.0, 0.0)]
        s2 = _xyz_signature(SimpleNamespace(), m1)
        assert s1 == s2  # rounded to 4 decimals

    def test_exception_returns_none(self):
        mol = MagicMock()
        mol.GetNumAtoms.side_effect = RuntimeError("boom")
        assert _xyz_signature(SimpleNamespace(), mol) is None


class TestXYZPersistenceHandlers:
    def _handlers(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            save = ctx.register_save_handler.call_args[0][0]
            load = ctx.register_load_handler.call_args[0][0]
            reset = ctx.register_document_reset_handler.call_args[0][0]
            return ctx, save, load, reset

    def test_save_no_molecule_returns_empty(self):
        ctx, save, _, _ = self._handlers()
        ctx.current_molecule = None
        assert save() == {}

    def test_save_collects_custom_symbols(self):
        ctx, save, _, _ = self._handlers()
        ctx.current_molecule = FakeMol(
            [
                FakeAtom("C", props={"custom_symbol": "C13"}),
                FakeAtom("H", 1),
            ],
            [(0, 0, 0), (1, 0, 0)],
        )
        assert save() == {"custom_labels": {0: "C13"}}

    def test_load_applies_labels(self):
        ctx, _, load, _ = self._handlers()
        mol = FakeMol([FakeAtom("C"), FakeAtom("N", 7)], [(0, 0, 0), (1, 0, 0)])
        ctx.current_molecule = mol
        load({"custom_labels": {"1": "N15"}})
        assert mol.GetAtomWithIdx(1).GetProp("custom_symbol") == "N15"

    def test_load_empty_data_no_raise(self):
        ctx, _, load, _ = self._handlers()
        ctx.current_molecule = FakeMol([FakeAtom("C")], [(0, 0, 0)])
        load({})

    def test_reset_reloads_open_window(self):
        ctx, _, _, reset = self._handlers()
        win = MagicMock()
        ctx.get_window.return_value = win
        reset()
        win.load_molecule.assert_called_once()


class TestXYZClickFilter:
    def _filter(self):
        qt = SimpleNamespace(MouseButton=SimpleNamespace(LeftButton=1))
        qevent = SimpleNamespace(
            Type=SimpleNamespace(
                MouseButtonPress="press", MouseButtonRelease="release"
            )
        )
        return _extract_method_as_fn(
            XYZ_EDITOR_PATH,
            "_ClickFilter",
            "eventFilter",
            extra_globals={"Qt": qt, "QEvent": qevent},
        )

    def _event(self, kind, x, y):
        pos = SimpleNamespace(
            toPoint=lambda: SimpleNamespace(x=lambda: x, y=lambda: y)
        )
        return SimpleNamespace(
            type=lambda: kind,
            button=lambda: 1,
            position=lambda: pos,
            modifiers=lambda: 0,
        )

    def test_click_within_threshold_fires_callback(self):
        fn = self._filter()
        calls = []
        self_ = SimpleNamespace(_press_pos=None, _callback=lambda *a: calls.append(a))
        assert fn(self_, "widget", self._event("press", 10, 10)) is False
        assert fn(self_, "widget", self._event("release", 13, 13)) is False
        assert len(calls) == 1
        assert calls[0][0] == 13 and calls[0][1] == 13

    def test_drag_beyond_threshold_ignored(self):
        fn = self._filter()
        calls = []
        self_ = SimpleNamespace(_press_pos=None, _callback=lambda *a: calls.append(a))
        fn(self_, "w", self._event("press", 10, 10))
        fn(self_, "w", self._event("release", 30, 30))
        assert calls == []

    def test_release_without_press_ignored(self):
        fn = self._filter()
        calls = []
        self_ = SimpleNamespace(_press_pos=None, _callback=lambda *a: calls.append(a))
        fn(self_, "w", self._event("release", 5, 5))
        assert calls == []
