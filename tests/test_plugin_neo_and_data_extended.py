"""
Second-wave extended tests:
  - Gaussian Input Generator Neo (Link0/route/coords assembly, Route Builder
    preview, route parsing, charge/mult validation, save_file %chk rewrite)
  - Paste from ChemDraw          (reconstruct_from_flat_text MolBlock output)
  - XYZ Editor                   (XYZ serialization, mol signature, handlers)
  - Animated XYZ Giffer          (frame stepping, play/pause, update scheduling)
  - Structural Updater           (trigger dispatch, apply-mode state, finalize)
  - Encrypted Project            (derive_key, export/import flow, salt header)

All dialog/widget classes inherit from mocked Qt bases, so methods are
extracted standalone via AST and invoked with fake `self` objects.
"""

from __future__ import annotations

import ast
import base64
import contextlib
import json
import logging
import os
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as real_numpy
import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

GAUSSIAN_PATH = (
    PLUGINS_DIR / "Gaussian_Input_Generator_Neo" / "gaussian_input_generator_neo.py"
)
CHEMDRAW_PATH = PLUGINS_DIR / "Paste_from_ChemDraw" / "paste_chemdraw.py"
XYZ_EDITOR_PATH = PLUGINS_DIR / "XYZ_Editor" / "xyz_editor.py"
GIFFER_PATH = PLUGINS_DIR / "Animated_XYZ_Giffer" / "animated_xyz_giffer.py"
STRUCTURAL_PATH = PLUGINS_DIR / "Structural_Updater" / "structural_updater.py"
ENCRYPTED_PATH = PLUGINS_DIR / "Encrypted_Project" / "encrypted_project.py"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def mocks_with_real_numpy():
    """
    Like mock_optional_imports(), but numpy resolves to the real package.

    ALL already-imported ``numpy.*`` submodules must be restored, not just the
    top-level package — see test_plugin_io_export_extended.py for the rationale
    (a MagicMock answer for numpy._core._methods corrupts numpy process-wide).
    """
    real_mods = {
        k: v
        for k, v in sys.modules.items()
        if k == "numpy" or k.startswith("numpy.")
    }
    with mock_optional_imports():
        sys.modules.update(real_mods)
        try:
            yield
        finally:
            for k in real_mods:
                sys.modules.pop(k, None)


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


class FakeSpin:
    def __init__(self, val):
        self._val = val
        self.set_calls = []

    def value(self):
        return self._val

    def setValue(self, v):
        self._val = v
        self.set_calls.append(v)


class FakeCombo:
    def __init__(self, text="", index=0):
        self._text = text
        self._index = index
        self.set_text_calls = []
        self.set_index_calls = []

    def currentText(self):
        return self._text

    def currentIndex(self):
        return self._index

    def setCurrentText(self, t):
        self._text = t
        self.set_text_calls.append(t)

    def setCurrentIndex(self, i):
        self._index = i
        self.set_index_calls.append(i)


class FakeCheck:
    def __init__(self, checked=False):
        self._checked = checked

    def isChecked(self):
        return self._checked

    def setChecked(self, v):
        self._checked = v


class FakeLine:
    def __init__(self, text=""):
        self._text = text

    def text(self):
        return self._text

    def setText(self, t):
        self._text = t


class FakeText:
    def __init__(self, text=""):
        self._text = text

    def toPlainText(self):
        return self._text

    def setText(self, t):
        self._text = t


class FakeLabel:
    def __init__(self):
        self._text = ""

    def setText(self, t):
        self._text = t

    def text(self):
        return self._text


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


def water_mol():
    return FakeMol(
        [FakeAtom("O", 8), FakeAtom("H", 1), FakeAtom("H", 1)],
        [(0.0, 0.0, 0.117), (0.0, 0.757, -0.469), (0.0, -0.757, -0.469)],
    )


# ===========================================================================
# Gaussian Input Generator Neo — generate_input_content / get_coords_string
# ===========================================================================

_gauss_generate = _extract_method_as_fn(
    GAUSSIAN_PATH, "GaussianSetupDialog", "generate_input_content"
)
_gauss_coords = _extract_method_as_fn(
    GAUSSIAN_PATH, "GaussianSetupDialog", "get_coords_string"
)


def _gauss_self(
    nproc=4,
    mem_val=4,
    mem_unit="GB",
    chk="",
    route="#P B3LYP/6-31G(d) Opt Freq",
    title="Generated by MoleditPy Plugin",
    charge=0,
    mult=1,
    tail="",
    mol=None,
):
    fake = SimpleNamespace()
    fake.nproc_spin = FakeSpin(nproc)
    fake.mem_spin = FakeSpin(mem_val)
    fake.mem_unit = FakeCombo(mem_unit)
    fake.chk_edit = FakeLine(chk)
    fake.keywords_edit = FakeLine(route)
    fake.title_edit = FakeLine(title)
    fake.charge_spin = FakeSpin(charge)
    fake.mult_spin = FakeSpin(mult)
    fake.tail_edit = FakeText(tail)
    fake.mol = mol
    fake.get_coords_string = lambda: _gauss_coords(fake)
    return fake


class TestGaussianGenerateContent:
    def test_link0_section(self):
        out = _gauss_generate(_gauss_self(nproc=8, mem_val=16, mem_unit="GB"))
        lines = out.splitlines()
        assert lines[0] == "%nprocshared=8"
        assert lines[1] == "%mem=16GB"

    def test_chk_auto_generated_from_hint(self):
        out = _gauss_generate(_gauss_self(), filename_hint="myjob")
        assert "%chk=myjob.chk" in out.splitlines()

    def test_custom_chk_gets_suffix(self):
        out = _gauss_generate(_gauss_self(chk="custom"))
        assert "%chk=custom.chk" in out.splitlines()

    def test_custom_chk_suffix_not_doubled(self):
        out = _gauss_generate(_gauss_self(chk="custom.chk"))
        assert "%chk=custom.chk" in out.splitlines()
        assert "custom.chk.chk" not in out

    def test_route_hash_prefixed_when_missing(self):
        out = _gauss_generate(_gauss_self(route="B3LYP/6-31G(d) Opt"))
        assert "# B3LYP/6-31G(d) Opt" in out.splitlines()

    def test_route_kept_when_already_prefixed(self):
        out = _gauss_generate(_gauss_self(route="#P HF/STO-3G SP"))
        assert "#P HF/STO-3G SP" in out.splitlines()

    def test_empty_title_defaults(self):
        out = _gauss_generate(_gauss_self(title="   "))
        assert "Gaussian Job" in out.splitlines()

    def test_charge_mult_line(self):
        out = _gauss_generate(_gauss_self(charge=-2, mult=3))
        assert "-2 3" in out.splitlines()

    def test_coords_included(self):
        out = _gauss_generate(_gauss_self(mol=water_mol()))
        assert "O" in out
        assert "0.757000" in out

    def test_tail_appended(self):
        out = _gauss_generate(_gauss_self(tail="B 1 2 F"))
        assert "B 1 2 F" in out

    def test_blank_tail_omitted(self):
        out = _gauss_generate(_gauss_self(tail="   \n  "))
        assert "   \n" not in out

    def test_wfn_line_appended_for_output_wfn(self):
        out = _gauss_generate(
            _gauss_self(route="# B3LYP/6-31G(d) SP Output=WFN"),
            filename_hint="job1",
        )
        assert "job1.wfn" in out.splitlines()

    def test_no_wfn_line_without_keyword(self):
        out = _gauss_generate(_gauss_self(), filename_hint="job1")
        assert "job1.wfn" not in out

    def test_ends_with_blank_line(self):
        out = _gauss_generate(_gauss_self())
        assert out.endswith("\n")

    def test_blank_line_after_route_and_title(self):
        out = _gauss_generate(_gauss_self(route="# HF/STO-3G SP", title="T"))
        lines = out.splitlines()
        r = lines.index("# HF/STO-3G SP")
        assert lines[r + 1] == ""
        assert lines[r + 2] == "T"
        assert lines[r + 3] == ""


class TestGaussianCoordsString:
    def test_no_mol_returns_empty(self):
        assert _gauss_coords(SimpleNamespace(mol=None)) == ""

    def test_standard_symbols(self):
        s = _gauss_coords(SimpleNamespace(mol=water_mol()))
        lines = s.splitlines()
        assert len(lines) == 3
        assert lines[0].startswith("O ")
        assert lines[1].startswith("H ")

    def test_custom_symbol_property_wins(self):
        mol = FakeMol(
            [FakeAtom("C", props={"custom_symbol": "C13"})], [(1.0, 2.0, 3.0)]
        )
        s = _gauss_coords(SimpleNamespace(mol=mol))
        assert s.startswith("C13")

    def test_coordinate_formatting(self):
        mol = FakeMol([FakeAtom("N", 7)], [(1.5, -2.25, 0.125)])
        s = _gauss_coords(SimpleNamespace(mol=mol))
        assert "1.500000" in s
        assert "-2.250000" in s

    def test_error_reported_as_string(self):
        broken = SimpleNamespace(mol=MagicMock())
        broken.mol.GetConformer.side_effect = RuntimeError("no conformer")
        s = _gauss_coords(broken)
        assert s.startswith("Error extracting coordinates")


# ===========================================================================
# Gaussian Neo — calc_initial_charge_mult / validate_charge_mult
# ===========================================================================


class TestGaussianChargeMult:
    def _calc(self, mol, formal_charge=0, chem_raises=False):
        chem = MagicMock()
        if chem_raises:
            chem.GetFormalCharge.side_effect = RuntimeError("boom")
        else:
            chem.GetFormalCharge.return_value = formal_charge
        fn = _extract_method_as_fn(
            GAUSSIAN_PATH,
            "GaussianSetupDialog",
            "calc_initial_charge_mult",
            extra_globals={"Chem": chem},
        )
        fake = SimpleNamespace(
            mol=mol,
            charge_spin=FakeSpin(0),
            mult_spin=FakeSpin(1),
            validate_charge_mult=MagicMock(),
        )
        fn(fake)
        return fake

    def test_no_mol_is_noop(self):
        fake = self._calc(None)
        assert fake.charge_spin.set_calls == []

    def test_charge_from_formal_charge(self):
        fake = self._calc(water_mol(), formal_charge=-1)
        assert fake.charge_spin.set_calls == [-1]

    def test_mult_radicals_plus_one(self):
        mol = FakeMol(
            [FakeAtom("C", radicals=2), FakeAtom("O", 8)], [(0, 0, 0), (1, 0, 0)]
        )
        fake = self._calc(mol)
        assert fake.mult_spin.set_calls == [3]

    def test_closed_shell_singlet(self):
        fake = self._calc(water_mol())
        assert fake.mult_spin.set_calls == [1]
        fake.validate_charge_mult.assert_called_once()

    def _validate(self, charge, mult, protons):
        qpalette = MagicMock()
        qt = MagicMock()
        qcolor = MagicMock()
        fn = _extract_method_as_fn(
            GAUSSIAN_PATH,
            "GaussianSetupDialog",
            "validate_charge_mult",
            extra_globals={"QPalette": qpalette, "QColor": qcolor, "Qt": qt},
        )
        atoms = [FakeAtom("X", num=protons)]
        mol = FakeMol(atoms, [(0, 0, 0)])
        fake = SimpleNamespace(
            mol=mol,
            charge_spin=MagicMock(),
            mult_spin=MagicMock(),
            default_palette="DEFAULT",
            update_preview=MagicMock(),
        )
        fake.charge_spin.value.return_value = charge
        fake.mult_spin.value.return_value = mult
        fn(fake)
        return fake

    def test_even_electrons_odd_mult_valid(self):
        fake = self._validate(charge=0, mult=1, protons=8)
        fake.charge_spin.setPalette.assert_called_once_with("DEFAULT")

    def test_even_electrons_even_mult_invalid(self):
        fake = self._validate(charge=0, mult=2, protons=8)
        args = fake.charge_spin.setPalette.call_args[0]
        assert args[0] != "DEFAULT"

    def test_odd_electrons_even_mult_valid(self):
        fake = self._validate(charge=0, mult=2, protons=7)
        fake.charge_spin.setPalette.assert_called_once_with("DEFAULT")

    def test_always_refreshes_preview(self):
        fake = self._validate(charge=0, mult=1, protons=8)
        fake.update_preview.assert_called_once()


# ===========================================================================
# Gaussian Neo — save_file %chk smart rewrite
# ===========================================================================


class TestGaussianSaveFile:
    def _save(self, tmp_path, content, fname="run1.gjf"):
        target = tmp_path / fname
        qfd = MagicMock()
        qfd.getSaveFileName.return_value = (str(target), "")
        qmb = MagicMock()
        fn = _extract_method_as_fn(
            GAUSSIAN_PATH,
            "GaussianSetupDialog",
            "save_file",
            extra_globals={"QFileDialog": qfd, "QMessageBox": qmb},
        )
        fake = SimpleNamespace(preview_text=FakeText(content))
        fn(fake)
        return target, qmb

    def test_existing_chk_rewritten_to_filename(self, tmp_path):
        content = "%nprocshared=4\n%chk=old.chk\n# HF SP\n\nT\n\n0 1\n"
        target, _ = self._save(tmp_path, content)
        text = target.read_text(encoding="utf-8")
        assert "%chk=run1.chk" in text
        assert "%chk=old.chk" not in text

    def test_missing_chk_inserted_at_top(self, tmp_path):
        content = "%nprocshared=4\n# HF SP\n\nT\n\n0 1\n"
        target, _ = self._save(tmp_path, content)
        text = target.read_text(encoding="utf-8")
        assert text.splitlines()[0] == "%chk=run1.chk"

    def test_wfn_filename_rewritten(self, tmp_path):
        content = "%chk=job.chk\n# SP Output=WFN\n\nT\n\n0 1\nH 0 0 0\n\njob.wfn\n"
        target, _ = self._save(tmp_path, content)
        text = target.read_text(encoding="utf-8")
        assert "run1.wfn" in text
        assert "job.wfn" not in text

    def test_success_dialog_shown(self, tmp_path):
        _, qmb = self._save(tmp_path, "%chk=a.chk\n# SP\n")
        qmb.information.assert_called_once()

    def test_cancel_writes_nothing(self, tmp_path):
        qfd = MagicMock()
        qfd.getSaveFileName.return_value = ("", "")
        fn = _extract_method_as_fn(
            GAUSSIAN_PATH,
            "GaussianSetupDialog",
            "save_file",
            extra_globals={"QFileDialog": qfd, "QMessageBox": MagicMock()},
        )
        fn(SimpleNamespace(preview_text=FakeText("x")))
        assert list(tmp_path.iterdir()) == []


# ===========================================================================
# Gaussian Neo — RouteBuilderDialog update_preview / parse_route
# ===========================================================================

_route_preview = _extract_method_as_fn(
    GAUSSIAN_PATH, "RouteBuilderDialog", "update_preview"
)


def _route_self(
    output_idx=1,
    method="B3LYP",
    basis="6-31G(d)",
    job_idx=0,
    tight=False,
    verytight=False,
    calcfc=False,
    maxcycles=False,
    raman=False,
    anharm=False,
    solv="None",
    solvent="Water",
    dispersion=False,
    pop="None",
    density=False,
    symmetry="Default",
    grid="Default",
    td=False,
    td_nstates=6,
    wfn=False,
):
    fake = SimpleNamespace(ui_ready=True)
    fake.output_level = FakeCombo(index=output_idx)
    fake.method_name = FakeCombo(method)
    fake.basis_set = FakeCombo(basis)
    fake.job_type = FakeCombo(index=job_idx)
    fake.opt_tight = FakeCheck(tight)
    fake.opt_verytight = FakeCheck(verytight)
    fake.opt_calcfc = FakeCheck(calcfc)
    fake.opt_maxcycles = FakeCheck(maxcycles)
    fake.freq_raman = FakeCheck(raman)
    fake.freq_anharm = FakeCheck(anharm)
    fake.solv_model = FakeCombo(solv)
    fake.solvent = FakeCombo(solvent)
    fake.dispersion = FakeCheck(dispersion)
    fake.pop_analysis = FakeCombo(pop)
    fake.density_chk = FakeCheck(density)
    fake.symmetry_combo = FakeCombo(symmetry)
    fake.grid_combo = FakeCombo(grid)
    fake.td_chk = FakeCheck(td)
    fake.td_nstates = FakeSpin(td_nstates)
    fake.wfn_chk = FakeCheck(wfn)
    fake.preview_label = FakeLabel()
    return fake


class TestRouteBuilderPreview:
    def _route(self, **kw):
        fake = _route_self(**kw)
        _route_preview(fake)
        return fake.preview_str

    def test_not_ready_skips(self):
        fake = _route_self()
        fake.ui_ready = False
        _route_preview(fake)
        assert not hasattr(fake, "preview_str")

    def test_default_opt_freq(self):
        assert self._route() == "# B3LYP/6-31G(d) Opt Freq"

    def test_output_level_prefixes(self):
        assert self._route(output_idx=0).startswith("#P ")
        assert self._route(output_idx=2).startswith("#T ")

    def test_opt_options_combined(self):
        r = self._route(job_idx=1, tight=True, calcfc=True)
        assert "Opt=(Tight, CalcFC)" in r

    def test_freq_options_combined(self):
        r = self._route(job_idx=2, raman=True, anharm=True)
        assert "Freq=(Raman, Anharmonic)" in r

    def test_sp_job(self):
        assert "SP" in self._route(job_idx=3)

    def test_scrf_with_solvent(self):
        r = self._route(solv="SMD", solvent="Acetonitrile")
        assert "SCRF=(SMD, Solvent=Acetonitrile)" in r

    def test_dispersion_appended(self):
        assert "EmpiricalDispersion=GD3BJ" in self._route(dispersion=True)

    def test_pop_nbo(self):
        assert "Pop=NBO" in self._route(pop="NBO (Pop=NBO)")

    def test_nosymm(self):
        assert "NoSymm" in self._route(symmetry="None (NoSymm)")

    def test_grid(self):
        assert "Integral(UltraFine)" in self._route(grid="UltraFine")

    def test_td_nstates(self):
        assert "TD=(NStates=10)" in self._route(td=True, td_nstates=10)

    def test_wfn_output(self):
        assert "Output=WFN" in self._route(wfn=True)


class TestRouteBuilderParse:
    def _parse(self, route):
        fn = _extract_method_as_fn(GAUSSIAN_PATH, "RouteBuilderDialog", "parse_route")
        fake = SimpleNamespace(
            job_type=FakeCombo(),
            method_type=FakeCombo(),
            method_name=FakeCombo(),
            basis_set=FakeCombo(),
            solv_model=FakeCombo(),
            wfn_chk=FakeCheck(False),
        )
        fn(fake, route)
        return fake

    def test_opt_freq_selects_first_job(self):
        fake = self._parse("# B3LYP/6-31G(d) Opt Freq")
        assert fake.job_type.set_index_calls == [0]

    def test_freq_only(self):
        fake = self._parse("# B3LYP/6-31G(d) Freq")
        assert fake.job_type.set_index_calls == [2]

    def test_method_and_basis_detected(self):
        fake = self._parse("# WB97XD/def2TZVP Opt")
        assert "WB97XD" in fake.method_name.set_text_calls
        assert "def2TZVP" in fake.basis_set.set_text_calls

    def test_smd_solvation(self):
        fake = self._parse("# B3LYP/6-31G(d) SP SCRF=(SMD, Solvent=Water)")
        assert fake.solv_model.set_text_calls == ["SMD"]

    def test_wfn_checkbox(self):
        fake = self._parse("# HF/STO-3G SP Output=WFN")
        assert fake.wfn_chk.isChecked()


class TestGaussianEntryPoints:
    def test_run_plugin_warns_without_molecule(self):
        with mock_optional_imports():
            mod = load_plugin(GAUSSIAN_PATH)
            ctx = make_context()
            ctx.current_mol = None
            mod.QMessageBox.warning.reset_mock()
            mod.run_plugin(ctx)
            mod.QMessageBox.warning.assert_called_once()

    def test_run_plugin_warns_for_empty_molecule(self):
        with mock_optional_imports():
            mod = load_plugin(GAUSSIAN_PATH)
            ctx = make_context()
            mol = MagicMock()
            mol.GetNumAtoms.return_value = 0
            ctx.current_mol = mol
            mod.QMessageBox.warning.reset_mock()
            mod.run_plugin(ctx)
            mod.QMessageBox.warning.assert_called_once()

    def test_initialize_registers_export_action(self):
        with mock_optional_imports():
            mod = load_plugin(GAUSSIAN_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_export_action.assert_called_once()
            assert "Gaussian" in ctx.add_export_action.call_args[0][0]


# ===========================================================================
# Paste from ChemDraw — reconstruct_from_flat_text MolBlock output
# ===========================================================================


def _reconstruct(flat_text):
    """Run reconstruct_from_flat_text and capture the MolBlock it builds."""
    with mock_optional_imports():
        mod = load_plugin(CHEMDRAW_PATH)
        captured = {}

        def capture(block):
            captured["block"] = block
            return "MOL_SENTINEL"

        mod.Chem.MolFromMolBlock = capture
        result = mod.reconstruct_from_flat_text(flat_text)
        return result, captured.get("block")


def _flat(n_atoms, n_bonds, body_tokens):
    header = f"molecule ChemDraw 2D  {n_atoms} {n_bonds} 0 0 0 0 0 0 0 0999 V2000"
    return header + " " + " ".join(body_tokens)


class TestChemDrawReconstruct:
    def test_no_v2000_returns_none(self):
        result, block = _reconstruct("just some random text")
        assert result is None
        assert block is None

    def test_short_header_returns_none(self):
        result, block = _reconstruct("1 1 V2000 rest")
        assert result is None
        assert block is None

    def test_non_numeric_counts_return_none(self):
        text = "a b c d e f g h i j k l V2000 stuff"
        result, block = _reconstruct(text)
        assert result is None

    def test_counts_line_format(self):
        tokens = ["0.0", "0.0", "0.0", "C"] + ["0"] * 12
        _, block = _reconstruct(_flat(1, 0, tokens))
        assert "  1  0  0  0  0  0  0  0  0  0999 V2000" in block

    def test_atom_line_formatting(self):
        tokens = ["1.5", "-2.25", "0.0", "N"] + ["0"] * 12
        _, block = _reconstruct(_flat(1, 0, tokens))
        atom_line = block.splitlines()[4]
        assert atom_line.startswith("    1.5000   -2.2500    0.0000 N")

    def test_phantom_f_token_skipped(self):
        tokens = ["F", "1.0", "2.0", "3.0", "O"] + ["0"] * 12
        _, block = _reconstruct(_flat(1, 0, tokens))
        atom_line = block.splitlines()[4]
        assert " O" in atom_line
        assert "1.0000" in atom_line

    def test_bond_line_built(self):
        tokens = (
            ["0.0", "0.0", "0.0", "C"]
            + ["0"] * 12
            + ["1.0", "0.0", "0.0", "C"]
            + ["0"] * 12
            + ["1", "2", "2", "0"]
        )
        _, block = _reconstruct(_flat(2, 1, tokens))
        assert "  1  2  2  0  0  0  0" in block

    def test_m_chg_preserved_before_end(self):
        tokens = (
            ["0.0", "0.0", "0.0", "N"]
            + ["0"] * 12
            + ["M", "CHG", "1", "1", "1", "M", "END"]
        )
        _, block = _reconstruct(_flat(1, 0, tokens))
        lines = block.splitlines()
        chg_idx = next(i for i, l in enumerate(lines) if l.startswith("M  CHG"))
        assert lines[chg_idx + 1] == "M  END"

    def test_block_ends_with_m_end(self):
        tokens = ["0.0", "0.0", "0.0", "C"] + ["0"] * 12
        _, block = _reconstruct(_flat(1, 0, tokens))
        assert block.splitlines()[-1] == "M  END"

    def test_control_characters_sanitized(self):
        tokens = ["0.0", "0.0", "0.0", "C"] + ["0"] * 12
        dirty = "\x01\x16" + _flat(1, 0, tokens) + "\x02"
        result, block = _reconstruct(dirty)
        assert result == "MOL_SENTINEL"

    def test_returns_parser_result(self):
        tokens = ["0.0", "0.0", "0.0", "C"] + ["0"] * 12
        result, _ = _reconstruct(_flat(1, 0, tokens))
        assert result == "MOL_SENTINEL"

    def test_two_atoms_two_lines(self):
        tokens = (
            ["0.0", "0.0", "0.0", "C"]
            + ["0"] * 12
            + ["1.4", "0.0", "0.0", "O"]
            + ["0"] * 12
        )
        _, block = _reconstruct(_flat(2, 0, tokens))
        lines = block.splitlines()
        assert " C" in lines[4]
        assert " O" in lines[5]


# ===========================================================================
# XYZ Editor — serialization, mol signature, persistence handlers
# ===========================================================================


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


# ===========================================================================
# Animated XYZ Giffer — frame stepping, play/pause, update scheduling
# ===========================================================================

_g_next = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "next_frame")
_g_prev = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "prev_frame")
_g_toggle = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "toggle_play")
_g_fps = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "set_fps")
_g_sched = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "schedule_update")
_g_status = _extract_method_as_fn(
    GIFFER_PATH, "AnimatedXYZPlayer", "update_status_silent"
)


def _giffer_self(n_frames=5, current=0, loop=True, playing=False, fps=10):
    fake = SimpleNamespace()
    fake.frames = [{"coords": []} for _ in range(n_frames)]
    fake.current_frame_idx = current
    fake.target_frame_idx = current
    fake.chk_loop = FakeCheck(loop)
    fake.is_playing = playing
    fake.fps = fps
    fake.schedule_update = MagicMock()
    fake.toggle_play = MagicMock()
    fake.timer = MagicMock()
    fake.btn_play = MagicMock()
    fake.context = MagicMock()
    fake.base_mol = "BASE"
    fake.last_display_mol = None
    fake.is_updating_view = False
    fake.pending_update = False
    fake.do_effective_update = MagicMock()
    fake.lbl_status = FakeLabel()
    fake.slider = MagicMock()
    return fake


class TestGifferFrameStepping:
    def test_next_advances(self):
        fake = _giffer_self(current=1)
        _g_next(fake)
        assert fake.target_frame_idx == 2
        fake.schedule_update.assert_called_once()

    def test_next_wraps_with_loop(self):
        fake = _giffer_self(current=4, loop=True)
        _g_next(fake)
        assert fake.target_frame_idx == 0

    def test_next_at_end_no_loop_stops_playback(self):
        fake = _giffer_self(current=4, loop=False, playing=True)
        _g_next(fake)
        fake.toggle_play.assert_called_once()
        fake.schedule_update.assert_not_called()

    def test_next_at_end_no_loop_not_playing_noop(self):
        fake = _giffer_self(current=4, loop=False, playing=False)
        _g_next(fake)
        fake.toggle_play.assert_not_called()
        fake.schedule_update.assert_not_called()

    def test_prev_steps_back(self):
        fake = _giffer_self(current=3)
        _g_prev(fake)
        assert fake.target_frame_idx == 2

    def test_prev_wraps_with_loop(self):
        fake = _giffer_self(current=0, loop=True)
        _g_prev(fake)
        assert fake.target_frame_idx == 4

    def test_prev_at_start_no_loop_noop(self):
        fake = _giffer_self(current=0, loop=False)
        _g_prev(fake)
        fake.schedule_update.assert_not_called()


class TestGifferPlayback:
    def test_play_starts_timer_with_fps_interval(self):
        fake = _giffer_self(playing=False, fps=20)
        _g_toggle(fake)
        assert fake.is_playing is True
        fake.timer.start.assert_called_once_with(50)
        fake.btn_play.setText.assert_called_with("Pause")

    def test_play_from_last_frame_restarts(self):
        fake = _giffer_self(current=4, playing=False)
        _g_toggle(fake)
        assert fake.target_frame_idx == 0
        fake.schedule_update.assert_called_once()

    def test_pause_stops_timer_and_pushes_molecule(self):
        fake = _giffer_self(playing=True)
        fake.last_display_mol = "FRAME_MOL"
        _g_toggle(fake)
        assert fake.is_playing is False
        fake.timer.stop.assert_called_once()
        assert fake.context.current_molecule == "FRAME_MOL"

    def test_pause_falls_back_to_base_mol(self):
        fake = _giffer_self(playing=True)
        fake.last_display_mol = None
        _g_toggle(fake)
        assert fake.context.current_molecule == "BASE"

    def test_set_fps_restarts_timer_only_when_playing(self):
        fake = _giffer_self(playing=True)
        _g_fps(fake, 25)
        assert fake.fps == 25
        fake.timer.start.assert_called_once_with(40)

        fake2 = _giffer_self(playing=False)
        _g_fps(fake2, 25)
        fake2.timer.start.assert_not_called()


class TestGifferScheduling:
    def test_schedule_runs_when_idle(self):
        fake = _giffer_self()
        _g_sched(fake)
        fake.do_effective_update.assert_called_once()
        assert fake.is_updating_view is True  # cleared by do_effective_update IRL

    def test_schedule_defers_when_busy(self):
        fake = _giffer_self()
        fake.is_updating_view = True
        _g_sched(fake)
        fake.do_effective_update.assert_not_called()
        assert fake.pending_update is True

    def test_status_label_and_slider_sync(self):
        fake = _giffer_self(n_frames=8, current=2)
        _g_status(fake)
        assert fake.lbl_status.text() == "Frame: 3 / 8"
        fake.slider.blockSignals.assert_any_call(True)
        fake.slider.setValue.assert_called_once_with(2)
        fake.slider.blockSignals.assert_any_call(False)
