"""
Tests for the GAMESS Input Generator plugin: generate_content
($CONTRL/$SYSTEM/$BASIS/$DATA assembly), preset save/load round-trip.

GamessSetupDialog inherits QDialog (mocked -> MagicMock), so methods are
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
GAMESS_PATH = PLUGINS_DIR / "Gamess_Input_Generator" / "gamess_input_generator.py"


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

_gamess_generate = _extract_method_as_fn(
    GAMESS_PATH, "GamessSetupDialog", "generate_content"
)


def _gamess_self(
    runtyp="OPTIMIZE",
    scftyp="RHF",
    dfttyp="None (Hartree-Fock)",
    charge=0,
    mult=1,
    nosym=True,
    gbasis="N31",
    ngauss=6,
    ndfunc=1,
    npfunc=0,
    mwords=100,
    mol=None,
):
    self = MagicMock()
    self.run_type = FakeCombo(runtyp)
    self.scf_type = FakeCombo(scftyp)
    self.dft_type = FakeCombo(dfttyp)
    self.charge_spin = FakeSpin(charge)
    self.mult_spin = FakeSpin(mult)
    self.chk_nosym = FakeCheck(nosym)
    self.basis_gbasis = FakeCombo(gbasis)
    self.basis_ngauss = FakeSpin(ngauss)
    self.basis_ndfunc = FakeSpin(ndfunc)
    self.basis_npfunc = FakeSpin(npfunc)
    self.mem_spin = FakeSpin(mwords)
    self.mol = mol
    return self


class TestGamessGenerateContent:
    def test_extraction_succeeded(self):
        assert _gamess_generate is not None

    def test_contrl_line_core_keywords(self):
        content = _gamess_generate(_gamess_self(charge=-1, mult=2, scftyp="UHF"))
        contrl = content.split("\n")[0]
        assert "$CONTRL" in contrl
        assert "SCFTYP=UHF" in contrl
        assert "RUNTYP=OPTIMIZE" in contrl
        assert "ICHARG=-1" in contrl
        assert "MULT=2" in contrl
        assert "COORD=CART" in contrl
        assert contrl.rstrip().endswith("$END")

    def test_hf_omits_dfttyp(self):
        contrl = _gamess_generate(_gamess_self()).split("\n")[0]
        assert "DFTTYP" not in contrl

    def test_dft_adds_dfttyp(self):
        contrl = _gamess_generate(_gamess_self(dfttyp="B3LYP")).split("\n")[0]
        assert "DFTTYP=B3LYP" in contrl

    def test_nosym_flag(self):
        assert "NOSYM=1" in _gamess_generate(_gamess_self(nosym=True)).split("\n")[0]
        assert "NOSYM" not in _gamess_generate(_gamess_self(nosym=False)).split("\n")[0]

    def test_system_mwords(self):
        content = _gamess_generate(_gamess_self(mwords=250))
        assert " $SYSTEM MWORDS=250 $END" in content

    def test_basis_line_with_polarization(self):
        content = _gamess_generate(_gamess_self(ndfunc=1, npfunc=2))
        basis = [ln for ln in content.split("\n") if "$BASIS" in ln][0]
        assert "GBASIS=N31" in basis
        assert "NGAUSS=6" in basis
        assert "NDFUNC=1" in basis
        assert "NPFUNC=2" in basis

    def test_basis_line_omits_zero_polarization(self):
        content = _gamess_generate(_gamess_self(ndfunc=0, npfunc=0))
        basis = [ln for ln in content.split("\n") if "$BASIS" in ln][0]
        assert "NDFUNC" not in basis
        assert "NPFUNC" not in basis

    def test_data_block_structure(self):
        content = _gamess_generate(_gamess_self(mol=water_mol()))
        lines = content.split("\n")
        idx = lines.index(" $DATA")
        assert lines[idx + 1] == "GAMESS Job generated by MoleditPy"
        assert lines[idx + 2] == "C1"
        assert lines[-1] == " $END"

    def test_atom_lines_have_atomic_number(self):
        content = _gamess_generate(_gamess_self(mol=water_mol()))
        atom_lines = [
            ln for ln in content.split("\n") if ln.startswith(("O", "H"))
        ]
        assert len(atom_lines) == 3
        o_parts = atom_lines[0].split()
        assert o_parts[0] == "O"
        assert o_parts[1] == "8"
        assert len(o_parts) == 5

    def test_no_molecule_still_has_data_block(self):
        content = _gamess_generate(_gamess_self(mol=None))
        lines = content.split("\n")
        idx = lines.index(" $DATA")
        # header, title, C1, then immediately $END
        assert lines[idx + 3] == " $END"


class TestGamessPresets:
    def test_round_trip(self, tmp_path):
        settings = tmp_path / "gamess.json"
        save = _extract_method_as_fn(
            GAMESS_PATH,
            "GamessSetupDialog",
            "save_presets_to_file",
            extra_globals={"SETTINGS_FILE": str(settings), "QMessageBox": MagicMock()},
        )
        load = _extract_method_as_fn(
            GAMESS_PATH,
            "GamessSetupDialog",
            "load_presets_from_file",
            extra_globals={"SETTINGS_FILE": str(settings)},
        )
        writer = MagicMock()
        writer.presets_data = {"p1": {"run_type": "ENERGY", "memory": 50}}
        save(writer)

        reader = MagicMock()
        load(reader)
        assert reader.presets_data == writer.presets_data

    def test_corrupt_json_gives_empty(self, tmp_path):
        settings = tmp_path / "gamess.json"
        settings.write_text("[[[", encoding="utf-8")
        load = _extract_method_as_fn(
            GAMESS_PATH,
            "GamessSetupDialog",
            "load_presets_from_file",
            extra_globals={"SETTINGS_FILE": str(settings)},
        )
        reader = MagicMock()
        load(reader)
        assert reader.presets_data == {}
