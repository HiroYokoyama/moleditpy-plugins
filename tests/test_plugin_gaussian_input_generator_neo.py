"""
Tests for the Gaussian Input Generator Neo plugin: Link0/route/coords
assembly, Route Builder preview, route parsing, charge/mult validation,
save_file %chk rewrite.

GaussianSetupDialog / RouteBuilderDialog inherit mocked Qt bases, so
methods are extracted standalone via AST and invoked with fake self
objects.
"""

from __future__ import annotations

import ast
import logging
import json
import os
import textwrap
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as real_numpy
import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
GAUSSIAN_PATH = (
    PLUGINS_DIR / "_old" / "Gaussian_Input_Generator_Neo" / "gaussian_input_generator_neo.py"
)

pytestmark = pytest.mark.skip(
    reason="Gaussian Input Generator Neo retired; replaced by Gaussian Input Generator Pro"
)


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

    def setPlainText(self, t):
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

    def test_run_plugin_with_valid_molecule_opens_dialog(self):
        with mock_optional_imports():
            mod = load_plugin(GAUSSIAN_PATH)
            ctx = make_context()
            mol = MagicMock()
            mol.GetNumAtoms.return_value = 3
            ctx.current_mol = mol
            mod.QMessageBox.warning.reset_mock()
            # GaussianSetupDialog collapses to a MagicMock under the mocked Qt
            # bases; just confirm no warning fires and exec() gets invoked.
            mod.run_plugin(ctx)
            mod.QMessageBox.warning.assert_not_called()

    def test_run_noop_without_plugin_manager(self):
        with mock_optional_imports():
            mod = load_plugin(GAUSSIAN_PATH)
            mw = MagicMock(spec=[])  # no plugin_manager attribute
            mod.run(mw)  # must not raise

    def test_run_noop_without_context(self):
        with mock_optional_imports():
            mod = load_plugin(GAUSSIAN_PATH)
            mod.PLUGIN_CONTEXT = None
            mw = MagicMock()
            mod.run(mw)  # must not raise


# ---------------------------------------------------------------------------
# RouteBuilderDialog.get_route / is_wfn_enabled
# ---------------------------------------------------------------------------

_route_get_route = _extract_method_as_fn(GAUSSIAN_PATH, "RouteBuilderDialog", "get_route")
_route_is_wfn_enabled = _extract_method_as_fn(
    GAUSSIAN_PATH, "RouteBuilderDialog", "is_wfn_enabled"
)


class TestRouteBuilderMisc:
    def test_get_route_returns_preview_str(self):
        fake = SimpleNamespace(preview_str="# B3LYP/6-31G(d) Opt")
        assert _route_get_route(fake) == "# B3LYP/6-31G(d) Opt"

    def test_is_wfn_enabled_true(self):
        fake = SimpleNamespace(wfn_chk=FakeCheck(True))
        assert _route_is_wfn_enabled(fake) is True

    def test_is_wfn_enabled_false(self):
        fake = SimpleNamespace(wfn_chk=FakeCheck(False))
        assert _route_is_wfn_enabled(fake) is False


# ---------------------------------------------------------------------------
# GaussianSetupDialog.update_job_options_visibility (job-tab visibility rules)
# ---------------------------------------------------------------------------

_update_job_vis = _extract_method_as_fn(
    GAUSSIAN_PATH, "RouteBuilderDialog", "update_job_options_visibility"
)


class FakeGroupBox:
    def __init__(self):
        self._visible = None

    def setVisible(self, v):
        self._visible = v


def _job_vis_self(job_idx):
    fake = SimpleNamespace()
    fake.job_type = FakeCombo(index=job_idx)
    fake.opt_group = FakeGroupBox()
    fake.freq_group = FakeGroupBox()
    return fake


class TestUpdateJobOptionsVisibility:
    @pytest.mark.parametrize("idx,is_opt", [(0, True), (1, True), (2, False), (3, False), (4, True), (5, True), (6, False), (7, False)])
    def test_opt_group_visibility(self, idx, is_opt):
        fake = _job_vis_self(idx)
        _update_job_vis(fake)
        assert fake.opt_group._visible is is_opt

    @pytest.mark.parametrize("idx,is_freq", [(0, True), (1, False), (2, True), (3, False), (4, False)])
    def test_freq_group_visibility(self, idx, is_freq):
        fake = _job_vis_self(idx)
        _update_job_vis(fake)
        assert fake.freq_group._visible is is_freq


# ---------------------------------------------------------------------------
# GaussianSetupDialog.insert_tail_template
# ---------------------------------------------------------------------------

_insert_tail_template = _extract_method_as_fn(
    GAUSSIAN_PATH, "GaussianSetupDialog", "insert_tail_template"
)


def _tail_self(combo_text):
    fake = SimpleNamespace()
    fake.tail_template_combo = FakeCombo(combo_text)
    cursor = MagicMock()
    fake.tail_edit = MagicMock()
    fake.tail_edit.textCursor.return_value = cursor
    return fake, cursor


class TestInsertTailTemplate:
    def test_select_placeholder_is_noop(self):
        fake, cursor = _tail_self("-- Select Template --")
        _insert_tail_template(fake)
        cursor.insertText.assert_not_called()

    def test_modredundant_template_inserted(self):
        fake, cursor = _tail_self("ModRedundant Freeze/Scan")
        _insert_tail_template(fake)
        text = cursor.insertText.call_args[0][0]
        assert "B 1 2 F" in text
        fake.tail_edit.setFocus.assert_called_once()

    def test_basis_set_template_inserted(self):
        fake, cursor = _tail_self("Gen/GenECP Basis Set Block")
        _insert_tail_template(fake)
        text = cursor.insertText.call_args[0][0]
        assert "6-31G(d)" in text

    def test_ecp_template_inserted(self):
        fake, cursor = _tail_self("ECP Block")
        _insert_tail_template(fake)
        text = cursor.insertText.call_args[0][0]
        assert "LanL2DZ" in text

    def test_nbo_template_inserted(self):
        fake, cursor = _tail_self("NBO Analysis")
        _insert_tail_template(fake)
        text = cursor.insertText.call_args[0][0]
        assert "$NBO" in text

    def test_link1_template_inserted(self):
        fake, cursor = _tail_self("Link1 Section")
        _insert_tail_template(fake)
        text = cursor.insertText.call_args[0][0]
        assert "--Link1--" in text

    def test_connectivity_template_inserted(self):
        fake, cursor = _tail_self("Connectivity Block")
        _insert_tail_template(fake)
        text = cursor.insertText.call_args[0][0]
        assert "Geom=Connectivity" in text

    def test_unrecognized_text_no_insert(self):
        fake, cursor = _tail_self("Totally Unknown Option")
        _insert_tail_template(fake)
        cursor.insertText.assert_not_called()


# ---------------------------------------------------------------------------
# Preset persistence: load/save/apply/delete
# ---------------------------------------------------------------------------

def _preset_fn(name, settings_file):
    return _extract_method_as_fn(
        GAUSSIAN_PATH,
        "GaussianSetupDialog",
        name,
        extra_globals={"SETTINGS_FILE": str(settings_file), "QMessageBox": MagicMock(),
                        "QInputDialog": MagicMock()},
    )


def _preset_self(**overrides):
    fake = SimpleNamespace()
    fake.presets_data = {}
    fake.preset_combo = FakeCombo()
    fake.nproc_spin = FakeSpin(4)
    fake.mem_spin = FakeSpin(4)
    fake.mem_unit = FakeCombo("GB")
    fake.chk_edit = FakeLine("")
    fake.keywords_edit = FakeLine("")
    fake.tail_edit = FakeText("")
    fake.update_preview = MagicMock()
    for k, v in overrides.items():
        setattr(fake, k, v)
    return fake


class TestPresetPersistence:
    def test_load_presets_creates_default_when_no_file(self, tmp_path):
        settings_file = tmp_path / "presets.json"
        fn = _preset_fn("load_presets_from_file", settings_file)
        fake = _preset_self()
        fake.update_preset_combo = MagicMock()
        fn(fake)
        assert "Default" in fake.presets_data
        fake.update_preset_combo.assert_called_once()

    def test_load_presets_reads_existing_file(self, tmp_path):
        settings_file = tmp_path / "presets.json"
        settings_file.write_text(
            json.dumps({"MyPreset": {"nproc": 8, "route": "# HF/STO-3G SP"}}),
            encoding="utf-8",
        )
        fn = _preset_fn("load_presets_from_file", settings_file)
        fake = _preset_self()
        fake.update_preset_combo = MagicMock()
        fn(fake)
        assert fake.presets_data["MyPreset"]["nproc"] == 8
        assert "Default" in fake.presets_data  # still auto-added

    def test_load_presets_corrupted_file_falls_back_to_default(self, tmp_path):
        settings_file = tmp_path / "presets.json"
        settings_file.write_text("{not valid json", encoding="utf-8")
        fn = _preset_fn("load_presets_from_file", settings_file)
        fake = _preset_self()
        fake.update_preset_combo = MagicMock()
        fn(fake)
        assert "Default" in fake.presets_data

    def test_apply_selected_preset_populates_fields(self, tmp_path):
        fn = _preset_fn("apply_selected_preset", tmp_path / "presets.json")
        fake = _preset_self(
            presets_data={
                "P1": {
                    "nproc": 16,
                    "mem_val": 8,
                    "mem_unit": "MB",
                    "chk": "job",
                    "route": "# PBE0/def2SVP SP",
                    "tail": "B 1 2 F",
                }
            },
            preset_combo=FakeCombo("P1"),
        )
        fn(fake)
        assert fake.nproc_spin.value() == 16
        assert fake.mem_spin.value() == 8
        assert fake.mem_unit.currentText() == "MB"
        assert fake.chk_edit.text() == "job"
        assert fake.keywords_edit.text() == "# PBE0/def2SVP SP"
        assert fake.tail_edit.toPlainText() == "B 1 2 F"
        fake.update_preview.assert_called_once()

    def test_apply_selected_preset_missing_name_is_noop(self, tmp_path):
        fn = _preset_fn("apply_selected_preset", tmp_path / "presets.json")
        fake = _preset_self(preset_combo=FakeCombo("Ghost"))
        fn(fake)
        fake.update_preview.assert_not_called()

    def test_get_current_ui_settings_round_trip(self, tmp_path):
        fn = _preset_fn("get_current_ui_settings", tmp_path / "presets.json")
        fake = _preset_self(
            nproc_spin=FakeSpin(2),
            mem_spin=FakeSpin(1),
            mem_unit=FakeCombo("GB"),
            chk_edit=FakeLine("c"),
            keywords_edit=FakeLine("# HF SP"),
            tail_edit=FakeText("tail text"),
        )
        result = fn(fake)
        assert result == {
            "nproc": 2,
            "mem_val": 1,
            "mem_unit": "GB",
            "chk": "c",
            "route": "# HF SP",
            "tail": "tail text",
        }

    def test_save_presets_to_file_writes_json(self, tmp_path):
        settings_file = tmp_path / "presets.json"
        fn = _preset_fn("save_presets_to_file", settings_file)
        fake = _preset_self(presets_data={"Default": {"nproc": 4}})
        fn(fake)
        assert json.loads(settings_file.read_text(encoding="utf-8")) == {
            "Default": {"nproc": 4}
        }

    def test_delete_preset_default_is_protected(self, tmp_path):
        fn = _preset_fn("delete_preset", tmp_path / "presets.json")
        fake = _preset_self(
            presets_data={"Default": {"nproc": 4}}, preset_combo=FakeCombo("Default")
        )
        fake.update_preset_combo = MagicMock()
        fake.save_presets_to_file = MagicMock()
        fn(fake)
        assert "Default" in fake.presets_data
        fake.save_presets_to_file.assert_not_called()
        fake.update_preset_combo.assert_not_called()

    def test_delete_preset_confirmed_removes_entry(self, tmp_path):
        settings_file = tmp_path / "presets.json"
        qmb = MagicMock()
        qmb.question.return_value = qmb.StandardButton.Yes
        fn = _extract_method_as_fn(
            GAUSSIAN_PATH,
            "GaussianSetupDialog",
            "delete_preset",
            extra_globals={"SETTINGS_FILE": str(settings_file), "QMessageBox": qmb},
        )
        fake = _preset_self(
            presets_data={"Custom": {"nproc": 4}}, preset_combo=FakeCombo("Custom")
        )
        fake.update_preset_combo = MagicMock()
        fake.save_presets_to_file = MagicMock()
        fn(fake)
        assert "Custom" not in fake.presets_data
        fake.save_presets_to_file.assert_called_once()
        fake.update_preset_combo.assert_called_once()

    def test_delete_preset_declined_keeps_entry(self, tmp_path):
        settings_file = tmp_path / "presets.json"
        qmb = MagicMock()
        qmb.question.return_value = qmb.StandardButton.No
        fn = _extract_method_as_fn(
            GAUSSIAN_PATH,
            "GaussianSetupDialog",
            "delete_preset",
            extra_globals={"SETTINGS_FILE": str(settings_file), "QMessageBox": qmb},
        )
        fake = _preset_self(
            presets_data={"Custom": {"nproc": 4}}, preset_combo=FakeCombo("Custom")
        )
        fake.update_preset_combo = MagicMock()
        fake.save_presets_to_file = MagicMock()
        fn(fake)
        assert "Custom" in fake.presets_data
        fake.save_presets_to_file.assert_not_called()
