"""
Tests for the Charge Editor plugin: initialize registration, formal-charge and
radical commits, total-charge / multiplicity summary, nearest-atom picking, and
molecule signature.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from conftest import (
    extract_function,
    load_plugin,
    make_context,
    mock_optional_imports,
)

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CHARGE_EDITOR_PATH = PLUGINS_DIR / "Charge_Editor" / "charge_editor.py"


class TestChargeEditorInitialize:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(CHARGE_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()
            assert "Charge Editor" in ctx.add_menu_action.call_args[0][0]

    def test_initialize_registers_document_reset_handler(self):
        with mock_optional_imports():
            mod = load_plugin(CHARGE_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_document_reset_handler.assert_called_once()

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(CHARGE_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_reset_reloads_open_window(self):
        with mock_optional_imports():
            mod = load_plugin(CHARGE_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            reset = ctx.register_document_reset_handler.call_args[0][0]
            win = MagicMock()
            ctx.get_window.return_value = win
            reset()
            win.load_molecule.assert_called_once()


# ---------------------------------------------------------------------------
# fake atoms / molecule
# ---------------------------------------------------------------------------


class _Atom:
    def __init__(self, idx, symbol="C", charge=0, radical=0):
        self._idx = idx
        self._symbol = symbol
        self._charge = charge
        self._radical = radical
        self.no_implicit = False

    def GetIdx(self):
        return self._idx

    def GetSymbol(self):
        return self._symbol

    def GetFormalCharge(self):
        return self._charge

    def SetFormalCharge(self, v):
        self._charge = v

    def GetNumRadicalElectrons(self):
        return self._radical

    def SetNumRadicalElectrons(self, v):
        self._radical = v

    def SetNoImplicit(self, v):
        self.no_implicit = v


class _Mol:
    def __init__(self, atoms):
        self._atoms = atoms

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetAtoms(self):
        return self._atoms

    def GetAtomWithIdx(self, i):
        return self._atoms[i]


# ---------------------------------------------------------------------------
# update_summary
# ---------------------------------------------------------------------------


class TestUpdateSummary:
    def _fn(self):
        return extract_function(CHARGE_EDITOR_PATH, "ChargeEditorWindow", "update_summary")

    def _self(self, mol):
        captured = {}
        label = SimpleNamespace(setText=lambda s: captured.__setitem__("text", s))
        return SimpleNamespace(
            context=SimpleNamespace(current_molecule=mol),
            summary_label=label,
        ), captured

    def test_neutral_singlet(self):
        fn = self._fn()
        self_, cap = self._self(_Mol([_Atom(0), _Atom(1)]))
        fn(self_)
        assert "Total charge: +0" in cap["text"]
        assert "Multiplicity: 1" in cap["text"]

    def test_anion_with_radical(self):
        fn = self._fn()
        mol = _Mol([_Atom(0, charge=-1), _Atom(1, radical=1)])
        self_, cap = self._self(mol)
        fn(self_)
        assert "Total charge: -1" in cap["text"]
        assert "Multiplicity: 2" in cap["text"]
        assert "1 unpaired" in cap["text"]

    def test_no_molecule(self):
        fn = self._fn()
        self_, cap = self._self(None)
        fn(self_)
        assert "Total charge: 0" in cap["text"]


# ---------------------------------------------------------------------------
# get_mol_signature
# ---------------------------------------------------------------------------


class TestSignature:
    def _fn(self):
        return extract_function(
            CHARGE_EDITOR_PATH, "ChargeEditorWindow", "get_mol_signature"
        )

    def test_none(self):
        assert self._fn()(SimpleNamespace(), None) is None

    def test_charge_change_changes_signature(self):
        fn = self._fn()
        m1 = _Mol([_Atom(0, charge=0)])
        m2 = _Mol([_Atom(0, charge=1)])
        assert fn(SimpleNamespace(), m1)[1:] != fn(SimpleNamespace(), m2)[1:]

    def test_radical_change_changes_signature(self):
        fn = self._fn()
        m1 = _Mol([_Atom(0, radical=0)])
        m2 = _Mol([_Atom(0, radical=2)])
        assert fn(SimpleNamespace(), m1)[1:] != fn(SimpleNamespace(), m2)[1:]


# ---------------------------------------------------------------------------
# on_charge_changed / on_radical_changed / clear_all_*
# ---------------------------------------------------------------------------


class _FakeRW:
    def __init__(self, base):
        self._atoms = base._atoms

    def GetAtomWithIdx(self, i):
        return self._atoms[i]

    def GetAtoms(self):
        return self._atoms


def _edit_fn(name):
    return extract_function(
        CHARGE_EDITOR_PATH,
        "ChargeEditorWindow",
        name,
        extra_globals={"Chem": SimpleNamespace(RWMol=_FakeRW), "QMessageBox": MagicMock()},
    )


def _edit_self(mol):
    committed = []
    return (
        SimpleNamespace(
            context=SimpleNamespace(
                current_molecule=mol, show_status_message=MagicMock()
            ),
            _commit=lambda rw, msg: committed.append((rw, msg)),
        ),
        committed,
    )


class TestEdits:
    def test_charge_change_sets_and_commits(self):
        fn = _edit_fn("on_charge_changed")
        mol = _Mol([_Atom(0), _Atom(1)])
        self_, committed = _edit_self(mol)
        fn(self_, 1, -1)
        assert mol._atoms[1].GetFormalCharge() == -1
        assert len(committed) == 1

    def test_radical_change_sets_and_commits(self):
        fn = _edit_fn("on_radical_changed")
        mol = _Mol([_Atom(0)])
        self_, committed = _edit_self(mol)
        fn(self_, 0, 2)
        assert mol._atoms[0].GetNumRadicalElectrons() == 2
        assert len(committed) == 1

    def test_no_molecule_no_commit(self):
        fn = _edit_fn("on_charge_changed")
        self_, committed = _edit_self(None)
        fn(self_, 0, 1)
        assert committed == []

    def test_clear_all_charges(self):
        fn = _edit_fn("clear_all_charges")
        mol = _Mol([_Atom(0, charge=1), _Atom(1, charge=-1)])
        self_, committed = _edit_self(mol)
        fn(self_)
        assert all(a.GetFormalCharge() == 0 for a in mol._atoms)
        assert len(committed) == 1

    def test_clear_all_radicals(self):
        fn = _edit_fn("clear_all_radicals")
        mol = _Mol([_Atom(0, radical=1), _Atom(1, radical=2)])
        self_, committed = _edit_self(mol)
        fn(self_)
        assert all(a.GetNumRadicalElectrons() == 0 for a in mol._atoms)
        assert len(committed) == 1


# ---------------------------------------------------------------------------
# _nearest_atom_to_point
# ---------------------------------------------------------------------------


class _Pos:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class _ConfMol:
    def __init__(self, coords):
        self._coords = coords
        self._atoms = [_Atom(i) for i in range(len(coords))]

    def GetConformer(self):
        return SimpleNamespace(GetAtomPosition=lambda i: _Pos(*self._coords[i]))

    def GetAtoms(self):
        return self._atoms


class TestNearestAtom:
    def _fn(self):
        return extract_function(
            CHARGE_EDITOR_PATH, "ChargeEditorWindow", "_nearest_atom_to_point"
        )

    def test_picks_closest(self):
        fn = self._fn()
        mol = _ConfMol([(0.0, 0.0, 0.0), (5.0, 0.0, 0.0), (10.0, 0.0, 0.0)])
        assert fn(SimpleNamespace(), mol, (4.6, 0.0, 0.0)) == 1
        assert fn(SimpleNamespace(), mol, (0.1, 0.1, 0.0)) == 0


# ---------------------------------------------------------------------------
# _highlight_radius (VDW-scaled selection halo)
# ---------------------------------------------------------------------------


class _NumAtom:
    def __init__(self, atomic_num):
        self._num = atomic_num

    def GetAtomicNum(self):
        return self._num


class _NumMol:
    def __init__(self, atomic_nums):
        self._atoms = [_NumAtom(n) for n in atomic_nums]

    def GetAtomWithIdx(self, i):
        return self._atoms[i]


_RVDW = {1: 1.2, 6: 1.7, 16: 1.8}


def _highlight_fn(rvdw=None):
    pt = SimpleNamespace(GetRvdw=rvdw or _RVDW.__getitem__)
    return extract_function(
        CHARGE_EDITOR_PATH,
        "ChargeEditorWindow",
        "_highlight_radius",
        extra_globals={
            "Chem": SimpleNamespace(GetPeriodicTable=lambda: pt),
            "logging": __import__("logging"),
        },
    )


class TestHighlightRadius:
    def test_scales_with_vdw_radius(self):
        fn = _highlight_fn()
        mol = _NumMol([6, 1, 16])
        assert abs(fn(mol, 0) - 1.7 * 1.2 * 0.3) < 1e-9  # C
        assert abs(fn(mol, 1) - 1.2 * 1.2 * 0.3) < 1e-9  # H
        assert fn(mol, 1) < fn(mol, 0) < fn(mol, 2)  # H < C < S

    def test_ghost_atom_uses_fixed_fallback(self):
        fn = _highlight_fn()
        mol = _NumMol([0])  # dummy atom, atomic number 0
        assert abs(fn(mol, 0) - 0.45) < 1e-9

    def test_rvdw_failure_uses_fixed_fallback(self):
        def boom(_n):
            raise RuntimeError("no such element")

        fn = _highlight_fn(rvdw=boom)
        mol = _NumMol([6])
        assert abs(fn(mol, 0) - 0.45) < 1e-9
