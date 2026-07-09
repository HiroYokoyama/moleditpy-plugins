"""
Tests for the NWChem Input Generator plugin: scf-block spin-state lines,
generate_content() input-file text generation, calc_initial_charge_mult(),
check_heavy_atoms(), and save_file() start-line rewriting.
"""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import FakeMol, extract_function, load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
NWCHEM_PATH = PLUGINS_DIR / "NWChem_Input_Generator" / "nwchem_input_generator.py"


class _Stub:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class _Text:
    def __init__(self, s):
        self._s = s

    def text(self):
        return self._s

    def setText(self, s):
        self._s = s


class _Combo:
    def __init__(self, text):
        self._t = text

    def currentText(self):
        return self._t


class _Spin:
    def __init__(self, v):
        self._v = v

    def value(self):
        return self._v

    def setValue(self, v):
        self._v = v


class TestNwchemScfSpinLines:
    def _fn(self):
        with mock_optional_imports():
            return load_plugin(NWCHEM_PATH)._scf_spin_lines

    def test_singlet_uses_rhf(self):
        assert self._fn()(1) == ["  rhf", "  singlet"]

    @pytest.mark.parametrize(
        "mult,keyword",
        [
            (2, "doublet"),
            (3, "triplet"),
            (4, "quartet"),
            (5, "quintet"),
            (6, "sextet"),
            (7, "septet"),
            (8, "octet"),
        ],
    )
    def test_open_shell_uses_uhf_keyword(self, mult, keyword):
        assert self._fn()(mult) == ["  uhf", f"  {keyword}"]

    def test_unnamed_high_multiplicity_falls_back_to_comment(self):
        lines = self._fn()(9)
        assert lines[0] == "  uhf"
        assert "9" in lines[1] and lines[1].lstrip().startswith("#")


# ---------------------------------------------------------------------------
# generate_content()
# ---------------------------------------------------------------------------


def _make_self(mol=None, module="dft", functional="b3lyp", task="optimize",
               basis="6-31G*", charge=0, mult=1, title="NWChem Job", start="[filename]"):
    return _Stub(
        mol=mol,
        title_edit=_Text(title),
        start_edit=_Text(start),
        module_combo=_Combo(module),
        functional_combo=_Combo(functional),
        task_combo=_Combo(task),
        basis_combo=_Combo(basis),
        charge_spin=_Spin(charge),
        mult_spin=_Spin(mult),
    )


class TestGenerateContent:
    def _fn(self):
        with mock_optional_imports():
            mod = load_plugin(NWCHEM_PATH)
        return extract_function(NWCHEM_PATH, "NwchemSetupDialog", "generate_content",
                                 extra_globals={"_scf_spin_lines": mod._scf_spin_lines})

    def test_dft_module_has_xc_and_mult_lines(self):
        fn = self._fn()
        content = fn(_make_self(module="dft", functional="pbe0", mult=1))
        assert "  xc pbe0" in content.splitlines()
        assert "  mult 1" in content.splitlines()

    def test_scf_module_uses_spin_lines(self):
        fn = self._fn()
        content = fn(_make_self(module="scf", mult=2))
        lines = content.splitlines()
        assert "  uhf" in lines
        assert "  doublet" in lines

    def test_nonzero_charge_emits_charge_line(self):
        fn = self._fn()
        content = fn(_make_self(charge=-1))
        assert "charge -1" in content.splitlines()

    def test_zero_charge_omits_charge_line(self):
        fn = self._fn()
        content = fn(_make_self(charge=0))
        assert not any(l.startswith("charge ") for l in content.splitlines())

    def test_none_mol_produces_empty_geometry_block(self):
        fn = self._fn()
        content = fn(_make_self(mol=None))
        lines = content.splitlines()
        i = lines.index("geometry units angstroms noautosym")
        assert lines[i + 1] == "end"

    def test_mol_atoms_written_with_coords(self):
        fn = self._fn()
        mol = FakeMol(["O", "H"], [(0.0, 0.0, 0.0), (0.96, 0.0, 0.0)])
        content = fn(_make_self(mol=mol))
        assert any(l.strip().startswith("O") and "0.000000" in l for l in content.splitlines())
        assert any(l.strip().startswith("H") and "0.960000" in l for l in content.splitlines())

    def test_filename_hint_overrides_start_edit(self):
        fn = self._fn()
        content = fn(_make_self(start="ignored"), "hinted_name")
        assert "start hinted_name" in content.splitlines()

    def test_task_line_uses_module_and_task_op(self):
        fn = self._fn()
        content = fn(_make_self(module="mp2", task="freq"))
        assert content.splitlines()[-1] == "task mp2 freq"

    def test_basis_library_line_present(self):
        fn = self._fn()
        content = fn(_make_self(basis="cc-pvtz"))
        assert '  * library "cc-pvtz"' in content.splitlines()


# ---------------------------------------------------------------------------
# calc_initial_charge_mult()
# ---------------------------------------------------------------------------


class TestCalcInitialChargeMult:
    def _fn(self):
        with mock_optional_imports():
            mod = load_plugin(NWCHEM_PATH)
        return extract_function(NWCHEM_PATH, "NwchemSetupDialog", "calc_initial_charge_mult",
                                 extra_globals={"Chem": mod.Chem, "logging": mod.logging})

    def test_none_mol_returns_without_raising(self):
        fn = self._fn()
        self_stub = _Stub(mol=None, charge_spin=_Spin(0), mult_spin=_Spin(1))
        fn(self_stub)  # must not raise
        assert self_stub.charge_spin.value() == 0

    def test_neutral_even_electron_molecule_is_singlet(self, monkeypatch):
        fn = self._fn()
        mol = FakeMol(["O", "H", "H"], [(0, 0, 0), (1, 0, 0), (0, 1, 0)])
        for atom, z in zip(mol.atoms, [8, 1, 1]):
            atom.GetAtomicNum = lambda z=z: z
        with mock_optional_imports():
            reloaded = load_plugin(NWCHEM_PATH)
        monkeypatch.setattr(reloaded.Chem, "GetFormalCharge", lambda m: 0)
        fn = extract_function(NWCHEM_PATH, "NwchemSetupDialog", "calc_initial_charge_mult",
                               extra_globals={"Chem": reloaded.Chem, "logging": reloaded.logging})
        self_stub = _Stub(mol=mol, charge_spin=_Spin(-99), mult_spin=_Spin(-99))
        fn(self_stub)
        assert self_stub.charge_spin.value() == 0
        assert self_stub.mult_spin.value() == 1

    def test_exception_is_swallowed(self):
        fn = self._fn()
        bad_mol = _Stub()  # no GetAtoms/Chem support -> triggers AttributeError internally
        self_stub = _Stub(mol=bad_mol, charge_spin=_Spin(0), mult_spin=_Spin(1))
        fn(self_stub)  # must not raise


# ---------------------------------------------------------------------------
# check_heavy_atoms()
# ---------------------------------------------------------------------------


class TestCheckHeavyAtoms:
    def _fn(self):
        return extract_function(NWCHEM_PATH, "NwchemSetupDialog", "check_heavy_atoms")

    def test_no_mol_no_warning_set(self):
        fn = self._fn()
        label = MagicMock()
        self_stub = _Stub(mol=None, ecp_warning=label)
        fn(self_stub)
        label.setText.assert_not_called()

    def test_heavy_atom_triggers_warning(self):
        fn = self._fn()
        mol = FakeMol(["I", "H"], [(0, 0, 0), (1, 0, 0)])
        mol.atoms[0].GetAtomicNum = lambda: 53  # Iodine, > 36
        mol.atoms[1].GetAtomicNum = lambda: 1
        label = MagicMock()
        self_stub = _Stub(mol=mol, ecp_warning=label)
        fn(self_stub)
        label.setText.assert_called_once()
        assert "ECP" in label.setText.call_args[0][0]

    def test_light_atoms_clear_warning(self):
        fn = self._fn()
        mol = FakeMol(["C", "H"], [(0, 0, 0), (1, 0, 0)])
        mol.atoms[0].GetAtomicNum = lambda: 6
        mol.atoms[1].GetAtomicNum = lambda: 1
        label = MagicMock()
        self_stub = _Stub(mol=mol, ecp_warning=label)
        fn(self_stub)
        label.setText.assert_called_once_with("")


# ---------------------------------------------------------------------------
# save_file()
# ---------------------------------------------------------------------------


class TestSaveFile:
    def _fn(self, path):
        qfd = MagicMock()
        qfd.getSaveFileName.return_value = (path, "NWChem Input (*.nw)")
        qmb = MagicMock()
        fn = extract_function(
            NWCHEM_PATH,
            "NwchemSetupDialog",
            "save_file",
            extra_globals={"QFileDialog": qfd, "QMessageBox": qmb, "os": os},
        )
        return fn, qfd, qmb

    def test_replaces_existing_start_line(self, tmp_path):
        path = str(tmp_path / "myjob.nw")
        fn, qfd, qmb = self._fn(path)
        preview = _Stub(toPlainText=lambda: 'title "x"\nstart old_name\ngeometry\nend')
        self_stub = _Stub(filename=None, preview_text=preview)
        fn(self_stub)
        content = Path(path).read_text(encoding="utf-8")
        assert "start myjob" in content
        assert "start old_name" not in content
        qmb.information.assert_called_once()

    def test_prepends_start_line_when_missing(self, tmp_path):
        path = str(tmp_path / "job2.nw")
        fn, qfd, qmb = self._fn(path)
        preview = _Stub(toPlainText=lambda: 'title "x"\ngeometry\nend')
        self_stub = _Stub(filename=None, preview_text=preview)
        fn(self_stub)
        content = Path(path).read_text(encoding="utf-8")
        assert content.startswith("start job2\n")

    def test_no_path_selected_does_nothing(self, tmp_path):
        fn, qfd, qmb = self._fn("")
        preview = _Stub(toPlainText=lambda: "content")
        self_stub = _Stub(filename=None, preview_text=preview)
        fn(self_stub)  # must not raise
        qmb.information.assert_not_called()

    def test_write_failure_shows_critical_not_raise(self):
        fn, qfd, qmb = self._fn("/nonexistent_dir_xyz/job.nw")
        preview = _Stub(toPlainText=lambda: "start x\ncontent")
        self_stub = _Stub(filename=None, preview_text=preview)
        fn(self_stub)  # must not raise
        qmb.critical.assert_called_once()
