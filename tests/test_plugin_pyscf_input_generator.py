"""
Tests for the PySCF Input Generator plugin: generate_content (script
generation, method selection, chkfile save logic), save_file placeholder
replacement.

PyscfSetupDialog inherits QDialog (mocked -> MagicMock), so methods are
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
PYSCF_PATH = PLUGINS_DIR / "PySCF_Input_Generator" / "pyscf_input_generator.py"


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

_pyscf_generate = _extract_method_as_fn(
    PYSCF_PATH, "PyscfSetupDialog", "generate_content"
)


def _pyscf_self(
    category="Hartree-Fock",
    functional="b3lyp",
    post_hf="MP2",
    post_hf_ref="RHF",
    basis="def2-svp",
    charge=0,
    mult=1,
    symmetry=True,
    mol=None,
):
    self = MagicMock()
    self.category_combo = FakeCombo(category)
    self.functional_combo = FakeCombo(functional)
    self.post_hf_combo = FakeCombo(post_hf)
    self.post_hf_ref_combo = FakeCombo(post_hf_ref)
    self.basis_combo = FakeCombo(basis)
    self.charge_spin = FakeSpin(charge)
    self.mult_spin = FakeSpin(mult)
    self.symmetry_check = FakeCheck(symmetry)
    self.mol = mol
    return self


class TestPyscfGenerateContent:
    def test_extraction_succeeded(self):
        assert _pyscf_generate is not None

    def test_imports_line_first(self):
        content = _pyscf_generate(_pyscf_self())
        assert content.split("\n")[0] == "from pyscf import gto, scf, dft, mp, cc"

    def test_molecule_block_params(self):
        content = _pyscf_generate(_pyscf_self(charge=-1, mult=3, basis="cc-pvtz"))
        assert "basis='cc-pvtz'," in content
        assert "charge=-1," in content
        assert "spin=2," in content  # 2S = mult - 1
        assert "symmetry=True" in content
        assert "mol.build()" in content

    def test_symmetry_false(self):
        content = _pyscf_generate(_pyscf_self(symmetry=False))
        assert "symmetry=False" in content

    def test_atom_coordinates_embedded(self):
        content = _pyscf_generate(_pyscf_self(mol=water_mol()))
        assert "O 0.00000000 0.00000000 0.11700000" in content
        assert "H 0.00000000 0.75700000 -0.46900000" in content

    def test_hf_singlet_uses_rhf(self):
        content = _pyscf_generate(_pyscf_self(category="Hartree-Fock", mult=1))
        assert "mf = scf.RHF(mol)" in content

    def test_hf_open_shell_uses_rohf(self):
        content = _pyscf_generate(_pyscf_self(category="Hartree-Fock", mult=2))
        assert "mf = scf.ROHF(mol)" in content

    def test_dft_singlet_rks_with_xc(self):
        content = _pyscf_generate(_pyscf_self(category="DFT", functional="m06"))
        assert "mf = dft.RKS(mol)" in content
        assert "mf.xc = 'm06'" in content

    def test_dft_open_shell_roks(self):
        content = _pyscf_generate(_pyscf_self(category="DFT", mult=2))
        assert "mf = dft.ROKS(mol)" in content

    @pytest.mark.parametrize(
        "ref,expected",
        [("RHF", "scf.RHF"), ("UHF", "scf.UHF"), ("ROHF", "scf.ROHF")],
    )
    def test_post_hf_reference_selection(self, ref, expected):
        content = _pyscf_generate(
            _pyscf_self(category="Post-HF (MP2, CCSD)", post_hf_ref=ref)
        )
        assert f"mf = {expected}(mol)" in content

    def test_mp2_block(self):
        content = _pyscf_generate(
            _pyscf_self(category="Post-HF (MP2, CCSD)", post_hf="MP2")
        )
        assert "pt = mp.MP2(mf)" in content
        assert "pt.kernel()" in content

    def test_ccsd_t_block(self):
        content = _pyscf_generate(
            _pyscf_self(category="Post-HF (MP2, CCSD)", post_hf="CCSD(T)")
        )
        assert "from pyscf.cc import ccsd_t" in content
        assert "e_t = ccsd_t.kernel(mycc, mycc.ao2mo())" in content
        assert "mycc.e_tot + e_t" in content

    def test_chkfile_uses_filename_hint(self):
        content = _pyscf_generate(_pyscf_self(), filename_hint="myjob")
        assert "mf.chkfile = 'myjob.chk'" in content


class TestPyscfSaveFile:
    def _save(self, tmp_path, content, target_name="final_job.py"):
        target = tmp_path / target_name
        qfd = MagicMock()
        qfd.getSaveFileName.return_value = (str(target), "")
        fn = _extract_method_as_fn(
            PYSCF_PATH,
            "PyscfSetupDialog",
            "save_file",
            extra_globals={"QFileDialog": qfd, "QMessageBox": MagicMock()},
        )
        self = MagicMock()
        self.filename = None
        self.preview_text = FakeText(content)
        fn(self)
        return target.read_text(encoding="utf-8")

    def test_chkfile_placeholder_replaced_with_basename(self, tmp_path):
        content = _pyscf_generate(_pyscf_self())  # contains job.chk hint default
        saved = self._save(tmp_path, content, "water_opt.py")
        assert "mf.chkfile = 'water_opt.chk'" in saved
        assert "job.chk" not in saved

    def test_preview_placeholder_replaced(self, tmp_path):
        content = _pyscf_generate(_pyscf_self(), filename_hint="[filename]")
        saved = self._save(tmp_path, content, "water_opt.py")
        assert "mf.chkfile = 'water_opt.chk'" in saved
        assert "[filename]" not in saved

    def test_missing_chkfile_inserted_before_kernel(self, tmp_path):
        content = "mf = scf.RHF(mol)\nmf.kernel()\n"
        saved = self._save(tmp_path, content, "abc.py")
        assert "mf.chkfile = 'abc.chk'\nmf.kernel()" in saved

    def test_no_kernel_appends_chk_line(self, tmp_path):
        content = "print('hello')"
        saved = self._save(tmp_path, content, "abc.py")
        assert "mf.chkfile = 'abc.chk'" in saved
