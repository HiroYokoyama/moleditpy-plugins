"""
Tests for the Psi4 Input Generator plugin: generate_content (input
assembly, auto reference), output-file save logic.

Psi4SetupDialog inherits QDialog (mocked -> MagicMock), so methods are
extracted standalone via AST and invoked with fake self objects that mimic
only the widgets each method touches.
"""

from __future__ import annotations

import ast
import json
import logging
import os
import textwrap
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
PSI4_PATH = PLUGINS_DIR / "Psi4_Input_Generator" / "psi4_input_generator.py"


def _extract_method_as_fn(
    path: Path,
    class_name: str,
    method_name: str,
    extra_globals: dict | None = None,
):
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
                if (
                    isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and item.name == method_name
                ):
                    func_src = ast.get_source_segment(source, item)
                    if func_src:
                        local_ns: dict = {}
                        globs = {"logging": logging, "os": os, "json": json}
                        globs.update(extra_globals or {})
                        exec(textwrap.dedent(func_src), globs, local_ns)
                        return local_ns[method_name]
    return None


class FakeSpin:
    def __init__(self, val):
        self._val = val

    def value(self):
        return self._val


class FakeCombo:
    def __init__(self, text, index=0):
        self._text = text
        self._index = index

    def currentText(self):
        return self._text

    def currentIndex(self):
        return self._index


class FakeCheck:
    def __init__(self, checked=False):
        self._checked = checked

    def isChecked(self):
        return self._checked


class FakeLine:
    def __init__(self, text=""):
        self._text = text

    def text(self):
        return self._text


class FakeText:
    def __init__(self, text=""):
        self._text = text

    def toPlainText(self):
        return self._text


class FakePos:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class FakeAtom:
    def __init__(self, symbol, num=6, props=None):
        self._symbol = symbol
        self._num = num
        self._props = props or {}

    def GetSymbol(self):
        return self._symbol

    def GetAtomicNum(self):
        return self._num

    def HasProp(self, key):
        return key in self._props

    def GetProp(self, key):
        return self._props[key]


class FakeConf:
    def __init__(self, coords):
        self._coords = coords

    def GetAtomPosition(self, i):
        return FakePos(*self._coords[i])


class FakeMol:
    def __init__(self, atoms, coords):
        assert len(atoms) == len(coords)
        self._atoms = atoms
        self._coords = coords

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetAtomWithIdx(self, i):
        return self._atoms[i]

    def GetConformer(self):
        return FakeConf(self._coords)


def water_mol():
    return FakeMol(
        [FakeAtom("O", 8), FakeAtom("H", 1), FakeAtom("H", 1)],
        [(0.0, 0.0, 0.117), (0.0, 0.757, -0.469), (0.0, -0.757, -0.469)],
    )

_psi4_generate = _extract_method_as_fn(
    PSI4_PATH, "Psi4SetupDialog", "generate_content"
)
_psi4_auto_ref = _extract_method_as_fn(
    PSI4_PATH, "Psi4SetupDialog", "update_auto_reference"
)


def _psi4_self(
    mem=2,
    threads=4,
    task="energy",
    method="b3lyp",
    basis="def2-svp",
    reference="rks",
    charge=0,
    mult=1,
    mol=None,
):
    self = MagicMock()
    self.mem_spin = FakeSpin(mem)
    self.thread_spin = FakeSpin(threads)
    self.task_combo = FakeCombo(task)
    self.method_combo = FakeCombo(method)
    self.basis_combo = FakeCombo(basis)
    self.ref_combo = FakeCombo(reference)
    self.charge_spin = FakeSpin(charge)
    self.mult_spin = FakeSpin(mult)
    self.mol = mol
    return self


class TestPsi4GenerateContent:
    def test_extraction_succeeded(self):
        assert _psi4_generate is not None

    def test_resources_lines(self):
        content = _psi4_generate(_psi4_self(mem=8, threads=16))
        assert "memory 8 GB" in content
        assert "set_num_threads(16)" in content

    def test_output_file_uses_hint(self):
        content = _psi4_generate(_psi4_self(), filename_hint="run01")
        assert "psi4.core.set_output_file('run01.out', False)" in content

    def test_molecule_block_charge_mult(self):
        content = _psi4_generate(_psi4_self(charge=1, mult=2, mol=water_mol()))
        lines = content.split("\n")
        mol_idx = lines.index("molecule {")
        assert lines[mol_idx + 1] == "1 2"
        assert lines[mol_idx + 2].split()[0] == "O"
        assert "}" in lines[mol_idx + 5]

    def test_set_block(self):
        content = _psi4_generate(_psi4_self(basis="cc-pvtz", reference="uhf"))
        assert "  basis cc-pvtz" in content
        assert "  reference uhf" in content

    def test_task_line_last(self):
        content = _psi4_generate(_psi4_self(task="optimize", method="mp2"))
        assert content.split("\n")[-1] == "optimize('mp2')"


class TestPsi4AutoReference:
    def _run(self, method, mult):
        self = MagicMock()
        self.method_combo = FakeCombo(method)
        self.mult_spin = FakeSpin(mult)
        _psi4_auto_ref(self)
        calls = [c.args[0] for c in self.ref_combo.setCurrentText.call_args_list]
        assert len(calls) == 1
        return calls[0]

    @pytest.mark.parametrize(
        "method,mult,expected",
        [
            ("b3lyp", 1, "rks"),
            ("b3lyp", 2, "uks"),
            ("wB97X-D", 3, "uks"),
            ("scf", 1, "rhf"),
            ("scf", 2, "uhf"),
            ("mp2", 1, "rhf"),
            ("mp2", 2, "uhf"),
            ("ccsd(t)", 1, "rhf"),
        ],
    )
    def test_reference_selection(self, method, mult, expected):
        assert self._run(method, mult) == expected

    def test_signals_blocked_and_restored(self):
        self = MagicMock()
        self.method_combo = FakeCombo("scf")
        self.mult_spin = FakeSpin(1)
        _psi4_auto_ref(self)
        assert self.ref_combo.blockSignals.call_args_list[0].args == (True,)
        assert self.ref_combo.blockSignals.call_args_list[-1].args == (False,)
        self.update_preview.assert_called_once()


class TestPsi4SaveFile:
    def _save(self, tmp_path, content, target_name="job.dat"):
        target = tmp_path / target_name
        qfd = MagicMock()
        qfd.getSaveFileName.return_value = (str(target), "")
        fn = _extract_method_as_fn(
            PSI4_PATH,
            "Psi4SetupDialog",
            "save_file",
            extra_globals={"QFileDialog": qfd, "QMessageBox": MagicMock()},
        )
        self = MagicMock()
        self.filename = None
        self.preview_text = FakeText(content)
        fn(self)
        return target.read_text(encoding="utf-8")

    def test_output_placeholder_replaced(self, tmp_path):
        content = _psi4_generate(_psi4_self(), filename_hint="[filename]")
        saved = self._save(tmp_path, content, "benzene_sp.dat")
        assert "psi4.core.set_output_file('benzene_sp.out', False)" in saved
        assert "[filename]" not in saved

    def test_missing_output_line_prepended(self, tmp_path):
        saved = self._save(tmp_path, "molecule { 0 1 }", "abc.dat")
        assert saved.startswith("psi4.core.set_output_file('abc.out', False)")
        assert "molecule { 0 1 }" in saved
