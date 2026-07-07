"""
Extended tests for the QM input generator plugins:
  - MOPAC Input Generator      (keyword line, multiplicity keywords, 1SCF opt flags)
  - GAMESS Input Generator     ($CONTRL/$SYSTEM/$BASIS/$DATA assembly)
  - PySCF Input Generator      (script generation, method selection, chkfile save logic)
  - Psi4 Input Generator       (input assembly, auto reference, output-file save logic)
  - NWChem Input Generator     (block assembly, scf multiplicity, start-line save logic)
  - Gaussian Input Generator Neo (Link0/route/tail assembly, Route Builder, save_file)

All generator dialogs inherit from QDialog (mocked -> MagicMock), so methods are
extracted standalone via AST (`_extract_method_as_fn`) and invoked with fake
`self` objects that mimic only the widgets each method touches.
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

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MOPAC_PATH = PLUGINS_DIR / "Mopac_Input_Generator" / "mopac_input_generator.py"
GAMESS_PATH = PLUGINS_DIR / "Gamess_Input_Generator" / "gamess_input_generator.py"
PYSCF_PATH = PLUGINS_DIR / "PySCF_Input_Generator" / "pyscf_input_generator.py"
PSI4_PATH = PLUGINS_DIR / "Psi4_Input_Generator" / "psi4_input_generator.py"
NWCHEM_PATH = PLUGINS_DIR / "NWChem_Input_Generator" / "nwchem_input_generator.py"
GAUSSIAN_PATH = (
    PLUGINS_DIR / "Gaussian_Input_Generator_Neo" / "gaussian_input_generator_neo.py"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# MOPAC Input Generator — generate_content
# ---------------------------------------------------------------------------

_mopac_generate = _extract_method_as_fn(
    MOPAC_PATH, "MopacSetupDialog", "generate_content"
)


def _mopac_self(
    keywords="PM7 PRECISE",
    charge=0,
    mult=1,
    nosym=False,
    title="MOPAC Calculation",
    comment="Generated by MoleditPy",
    mol=None,
):
    self = MagicMock()
    self.keywords_edit = FakeLine(keywords)
    self.charge_spin = FakeSpin(charge)
    self.mult_spin = FakeSpin(mult)
    self.chk_nosym = FakeCheck(nosym)
    self.title_edit = FakeLine(title)
    self.comment_edit = FakeLine(comment)
    self.mol = mol
    return self


class TestMopacGenerateContent:
    def test_extraction_succeeded(self):
        assert _mopac_generate is not None

    def test_header_lines(self):
        content = _mopac_generate(_mopac_self())
        lines = content.split("\n")
        assert lines[0].startswith("PM7 PRECISE CHARGE=0")
        assert lines[1] == "MOPAC Calculation"
        assert lines[2] == "Generated by MoleditPy"

    def test_negative_charge_in_keywords(self):
        content = _mopac_generate(_mopac_self(charge=-2))
        assert "CHARGE=-2" in content.split("\n")[0]

    @pytest.mark.parametrize(
        "mult,expected",
        [
            (1, "SINGLET"),
            (2, "DOUBLET"),
            (3, "TRIPLET"),
            (4, "QUARTET"),
            (5, "QUINTET"),
            (6, "SEXTET"),
        ],
    )
    def test_named_multiplicity_keywords(self, mult, expected):
        content = _mopac_generate(_mopac_self(mult=mult))
        assert expected in content.split("\n")[0]

    @pytest.mark.parametrize("mult,ms", [(7, "MS=3.0"), (8, "MS=3.5"), (10, "MS=4.5")])
    def test_high_multiplicity_uses_ms(self, mult, ms):
        content = _mopac_generate(_mopac_self(mult=mult))
        assert ms in content.split("\n")[0]

    def test_nosym_appended_when_checked(self):
        assert "NOSYM" in _mopac_generate(_mopac_self(nosym=True)).split("\n")[0]
        assert "NOSYM" not in _mopac_generate(_mopac_self(nosym=False)).split("\n")[0]

    def test_no_molecule_gives_three_lines(self):
        content = _mopac_generate(_mopac_self(mol=None))
        assert len(content.split("\n")) == 3

    def test_geometry_lines_with_opt_flag_1(self):
        content = _mopac_generate(_mopac_self(mol=water_mol()))
        lines = content.split("\n")
        assert len(lines) == 6
        assert lines[3].startswith("O")
        # optimization flag 1 between coordinates
        assert lines[3].split()[2] == "1"
        assert lines[3].split()[4] == "1"
        assert lines[3].split()[6] == "1"

    def test_1scf_gives_fixed_geometry_flag_0(self):
        content = _mopac_generate(_mopac_self(keywords="PM7 1SCF", mol=water_mol()))
        atom_line = content.split("\n")[3]
        assert atom_line.split()[2] == "0"
        assert atom_line.split()[4] == "0"

    def test_1scf_detection_case_insensitive(self):
        content = _mopac_generate(_mopac_self(keywords="pm7 1scf", mol=water_mol()))
        assert content.split("\n")[3].split()[2] == "0"

    def test_coordinate_formatting(self):
        content = _mopac_generate(_mopac_self(mol=water_mol()))
        line = content.split("\n")[4]
        parts = line.split()
        assert parts[0] == "H"
        assert float(parts[1]) == pytest.approx(0.0)
        assert float(parts[3]) == pytest.approx(0.757)
        assert float(parts[5]) == pytest.approx(-0.469)


# ---------------------------------------------------------------------------
# MOPAC — presets save/load round-trip
# ---------------------------------------------------------------------------


class TestMopacPresets:
    def _load_fn(self, settings_file):
        return _extract_method_as_fn(
            MOPAC_PATH,
            "MopacSetupDialog",
            "load_presets_from_file",
            extra_globals={"SETTINGS_FILE": str(settings_file)},
        )

    def _save_fn(self, settings_file):
        return _extract_method_as_fn(
            MOPAC_PATH,
            "MopacSetupDialog",
            "save_presets_to_file",
            extra_globals={"SETTINGS_FILE": str(settings_file), "QMessageBox": MagicMock()},
        )

    def test_round_trip(self, tmp_path):
        settings = tmp_path / "mopac.json"
        save = self._save_fn(settings)
        load = self._load_fn(settings)

        writer = MagicMock()
        writer.presets_data = {"My Preset": {"keywords": "PM6 PRECISE", "nosym": True}}
        save(writer)
        assert settings.exists()

        reader = MagicMock()
        load(reader)
        assert reader.presets_data == writer.presets_data
        reader.update_preset_combo.assert_called_once()

    def test_load_missing_file_gives_empty(self, tmp_path):
        load = self._load_fn(tmp_path / "absent.json")
        reader = MagicMock()
        load(reader)
        assert reader.presets_data == {}

    def test_load_corrupt_json_gives_empty(self, tmp_path):
        settings = tmp_path / "corrupt.json"
        settings.write_text("{not json!", encoding="utf-8")
        load = self._load_fn(settings)
        reader = MagicMock()
        load(reader)
        assert reader.presets_data == {}


# ---------------------------------------------------------------------------
# GAMESS Input Generator — generate_content
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# PySCF Input Generator — generate_content
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Psi4 Input Generator — generate_content + update_auto_reference
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Legacy run(mw) guard — all six generators warn when no molecule is loaded
# ---------------------------------------------------------------------------


class TestRunNoMoleculeGuard:
    @pytest.mark.parametrize(
        "path",
        [MOPAC_PATH, GAMESS_PATH, PYSCF_PATH, PSI4_PATH, NWCHEM_PATH],
        ids=lambda p: p.stem,
    )
    def test_run_warns_and_returns_without_molecule(self, path):
        with mock_optional_imports():
            mod = load_plugin(path)
            mw = MagicMock()
            mw.current_mol = None
            mod.QMessageBox.warning.reset_mock()
            mod.run(mw)
            mod.QMessageBox.warning.assert_called_once()
            # No dialog should have been constructed (exec never reached)
