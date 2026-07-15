"""
Tests for the Bond Editor plugin: initialize registration, bond-type label
mapping, ring detection (_moving_side), bond-length math, add/delete guards,
and molecule signature.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as real_numpy

from conftest import (
    extract_function,
    load_plugin,
    make_context,
    mock_optional_imports,
)

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
BOND_EDITOR_PATH = PLUGINS_DIR / "Bond_Editor" / "bond_editor.py"


class TestBondEditorInitialize:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(BOND_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()
            assert "Bond Editor" in ctx.add_menu_action.call_args[0][0]

    def test_initialize_registers_document_reset_handler(self):
        with mock_optional_imports():
            mod = load_plugin(BOND_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_document_reset_handler.assert_called_once()

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(BOND_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_reset_reloads_open_window(self):
        with mock_optional_imports():
            mod = load_plugin(BOND_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            reset = ctx.register_document_reset_handler.call_args[0][0]
            win = MagicMock()
            ctx.get_window.return_value = win
            reset()
            win.load_molecule.assert_called_once()

    def test_bond_type_label_round_trip(self):
        with mock_optional_imports():
            mod = load_plugin(BOND_EDITOR_PATH)
            # with mocked rdkit, BondType attributes are memoized MagicMocks
            assert mod.bond_type_from_label("Double") is mod.Chem.BondType.DOUBLE
            assert mod.bond_type_from_label("unknown") is mod.Chem.BondType.SINGLE


class TestLabelFromBondType:
    def _fn(self):
        return extract_function(BOND_EDITOR_PATH, None, "label_from_bond_type")

    def test_known_types(self):
        fn = self._fn()
        assert fn("BondType.SINGLE") == "Single"
        assert fn("BondType.DOUBLE") == "Double"
        assert fn("BondType.TRIPLE") == "Triple"
        assert fn("BondType.AROMATIC") == "Aromatic"

    def test_unknown_type_defaults_to_single(self):
        fn = self._fn()
        assert fn("BondType.DATIVE") == "Single"

    def test_bare_name_without_prefix(self):
        fn = self._fn()
        assert fn("DOUBLE") == "Double"


# ---------------------------------------------------------------------------
# fake bond graph
# ---------------------------------------------------------------------------


class _GraphMol:
    """Adjacency-list molecule exposing GetAtomWithIdx().GetNeighbors()."""

    def __init__(self, n, edges):
        self._adj = {i: set() for i in range(n)}
        for a, b in edges:
            self._adj[a].add(b)
            self._adj[b].add(a)

    def GetAtomWithIdx(self, i):
        return SimpleNamespace(
            GetNeighbors=lambda: [
                SimpleNamespace(GetIdx=lambda j=j: j) for j in sorted(self._adj[i])
            ]
        )


class TestMovingSide:
    def _fn(self):
        return extract_function(BOND_EDITOR_PATH, "BondEditorWindow", "_moving_side")

    def test_chain_returns_end_side(self):
        fn = self._fn()
        mol = _GraphMol(4, [(0, 1), (1, 2), (2, 3)])
        assert fn(SimpleNamespace(), mol, 1, 2) == {2, 3}

    def test_terminal_bond_returns_single_atom(self):
        fn = self._fn()
        mol = _GraphMol(3, [(0, 1), (1, 2)])
        assert fn(SimpleNamespace(), mol, 1, 2) == {2}

    def test_ring_bond_returns_none(self):
        fn = self._fn()
        mol = _GraphMol(3, [(0, 1), (1, 2), (2, 0)])
        assert fn(SimpleNamespace(), mol, 0, 1) is None

    def test_branch_on_end_side_included(self):
        fn = self._fn()
        # 0-1, 1-2, 2-3, 2-4 : moving side of bond 1-2 is {2,3,4}
        mol = _GraphMol(5, [(0, 1), (1, 2), (2, 3), (2, 4)])
        assert fn(SimpleNamespace(), mol, 1, 2) == {2, 3, 4}

    def test_ring_further_away_still_none(self):
        fn = self._fn()
        # bond 0-1 where 1 is in a ring back to 0 via 2,3
        mol = _GraphMol(4, [(0, 1), (1, 2), (2, 3), (3, 0)])
        assert fn(SimpleNamespace(), mol, 0, 1) is None

    def test_disconnected_fragment_not_included(self):
        fn = self._fn()
        mol = _GraphMol(4, [(0, 1), (2, 3)])
        assert fn(SimpleNamespace(), mol, 0, 1) == {1}


# ---------------------------------------------------------------------------
# set_bond_length
# ---------------------------------------------------------------------------


class _FakePos:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class _FakeLenConf:
    def __init__(self, coords):
        self.coords = {i: list(c) for i, c in enumerate(coords)}

    def GetAtomPosition(self, i):
        return _FakePos(*self.coords[i])

    def SetAtomPosition(self, i, p):
        self.coords[i] = [p.x, p.y, p.z]


class _FakeLenMol:
    def __init__(self, coords):
        self.conf = _FakeLenConf(coords)

    def GetNumConformers(self):
        return 1

    def GetConformer(self):
        return self.conf


class _FakeLenRW:
    def __init__(self, base):
        self.conf = base.conf

    def GetConformer(self):
        return self.conf


def _set_bond_length_fn():
    chem_ns = SimpleNamespace(RWMol=_FakeLenRW)
    return extract_function(
        BOND_EDITOR_PATH,
        "BondEditorWindow",
        "set_bond_length",
        extra_globals={
            "Chem": chem_ns,
            "np": real_numpy,
            "Point3D": _FakePos,
            "QMessageBox": MagicMock(),
        },
    )


def _len_self(mol, pair, side):
    committed = []
    return (
        SimpleNamespace(
            context=SimpleNamespace(
                current_molecule=mol, show_status_message=MagicMock()
            ),
            _row_bond_atoms=lambda row: pair,
            _moving_side=lambda m, b, e: side,
            _commit=lambda rw, msg: committed.append((rw, msg)),
            load_molecule=MagicMock(),
        ),
        committed,
    )


class TestSetBondLength:
    def test_stretches_end_atom_along_axis(self):
        fn = _set_bond_length_fn()
        mol = _FakeLenMol([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
        self_, committed = _len_self(mol, (0, 1), {1})
        fn(self_, 0, 1.5)
        assert len(committed) == 1
        assert committed[0][0].conf.coords[1] == [1.5, 0.0, 0.0]

    def test_moves_whole_side_fragment(self):
        fn = _set_bond_length_fn()
        mol = _FakeLenMol([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)])
        self_, committed = _len_self(mol, (0, 1), {1, 2})
        fn(self_, 0, 1.0)
        conf = committed[0][0].conf
        assert conf.coords[1] == [1.0, 0.0, 0.0]
        assert conf.coords[2] == [2.0, 0.0, 0.0]  # rigidly translated
        assert conf.coords[0] == [0.0, 0.0, 0.0]  # anchor side untouched

    def test_diagonal_direction_preserved(self):
        fn = _set_bond_length_fn()
        mol = _FakeLenMol([(0.0, 0.0, 0.0), (1.0, 1.0, 0.0)])
        self_, committed = _len_self(mol, (0, 1), {1})
        fn(self_, 0, 2.0 * 2.0**0.5)
        x, y, z = committed[0][0].conf.coords[1]
        assert abs(x - 2.0) < 1e-9 and abs(y - 2.0) < 1e-9 and abs(z) < 1e-9

    def test_ring_bond_refused_with_message(self):
        fn = _set_bond_length_fn()
        mol = _FakeLenMol([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
        self_, committed = _len_self(mol, (0, 1), None)
        fn(self_, 0, 1.5)
        assert committed == []
        self_.context.show_status_message.assert_called_once()
        self_.load_molecule.assert_called_once()

    def test_coincident_atoms_refused(self):
        fn = _set_bond_length_fn()
        mol = _FakeLenMol([(1.0, 1.0, 1.0), (1.0, 1.0, 1.0)])
        self_, committed = _len_self(mol, (0, 1), {1})
        fn(self_, 0, 1.5)
        assert committed == []
        self_.context.show_status_message.assert_called_once()


# ---------------------------------------------------------------------------
# add_bond / delete_selected_bonds guards
# ---------------------------------------------------------------------------


class _Idx:
    def __init__(self, row):
        self._row = row

    def row(self):
        return self._row


class _FakeAddRW:
    def __init__(self, base):
        self.added = []

    def AddBond(self, a1, a2, t):
        self.added.append((a1, a2, t))


def _add_bond_fn():
    chem_ns = SimpleNamespace(RWMol=_FakeAddRW)
    return extract_function(
        BOND_EDITOR_PATH,
        "BondEditorWindow",
        "add_bond",
        extra_globals={
            "Chem": chem_ns,
            "QMessageBox": MagicMock(),
            "bond_type_from_label": lambda label: label.upper(),
        },
    )


def _add_bond_self(mol, a1, a2, type_label="Single"):
    committed = []
    return (
        SimpleNamespace(
            context=SimpleNamespace(
                current_molecule=mol, show_status_message=MagicMock()
            ),
            spin_a1=SimpleNamespace(value=lambda: a1),
            spin_a2=SimpleNamespace(value=lambda: a2),
            add_type_combo=SimpleNamespace(currentText=lambda: type_label),
            _commit=lambda rw, msg: committed.append((rw, msg)),
        ),
        committed,
    )


def _bond_mol(n, existing_pairs=()):
    pairs = {frozenset(p) for p in existing_pairs}
    return SimpleNamespace(
        GetNumAtoms=lambda: n,
        GetBondBetweenAtoms=lambda a, b: (
            "bond" if frozenset((a, b)) in pairs else None
        ),
    )


class TestAddBond:
    def test_adds_bond_and_commits(self):
        fn = _add_bond_fn()
        self_, committed = _add_bond_self(_bond_mol(3), 0, 2, "Double")
        fn(self_)
        assert committed[0][0].added == [(0, 2, "DOUBLE")]

    def test_self_bond_refused(self):
        fn = _add_bond_fn()
        self_, committed = _add_bond_self(_bond_mol(3), 1, 1)
        fn(self_)
        assert committed == []
        self_.context.show_status_message.assert_called_once()

    def test_existing_bond_refused(self):
        fn = _add_bond_fn()
        self_, committed = _add_bond_self(_bond_mol(3, [(0, 1)]), 1, 0)
        fn(self_)
        assert committed == []

    def test_no_molecule_refused(self):
        fn = _add_bond_fn()
        self_, committed = _add_bond_self(None, 0, 1)
        fn(self_)
        assert committed == []


class _FakeDelRW:
    def __init__(self, base):
        self.removed = []

    def RemoveBond(self, a1, a2):
        self.removed.append((a1, a2))


def _delete_bonds_fn():
    chem_ns = SimpleNamespace(RWMol=_FakeDelRW)
    return extract_function(
        BOND_EDITOR_PATH,
        "BondEditorWindow",
        "delete_selected_bonds",
        extra_globals={"Chem": chem_ns, "QMessageBox": MagicMock()},
    )


def _delete_self(mol, selected_rows, row_pairs):
    committed = []
    return (
        SimpleNamespace(
            context=SimpleNamespace(
                current_molecule=mol, show_status_message=MagicMock()
            ),
            table=SimpleNamespace(
                selectedIndexes=lambda: [_Idx(r) for r in selected_rows]
            ),
            _row_bond_atoms=lambda row: row_pairs.get(row),
            _commit=lambda rw, msg: committed.append((rw, msg)),
        ),
        committed,
    )


class TestDeleteSelectedBonds:
    def test_removes_selected_pairs(self):
        fn = _delete_bonds_fn()
        self_, committed = _delete_self(
            "mol", [0, 2], {0: (0, 1), 1: (1, 2), 2: (2, 3)}
        )
        fn(self_)
        assert committed[0][0].removed == [(0, 1), (2, 3)]

    def test_no_selection_shows_message(self):
        fn = _delete_bonds_fn()
        self_, committed = _delete_self("mol", [], {})
        fn(self_)
        assert committed == []
        self_.context.show_status_message.assert_called_once()

    def test_duplicate_rows_deduplicated(self):
        fn = _delete_bonds_fn()
        self_, committed = _delete_self("mol", [1, 1], {1: (1, 2)})
        fn(self_)
        assert committed[0][0].removed == [(1, 2)]


# ---------------------------------------------------------------------------
# _row_bond_atoms / get_mol_signature
# ---------------------------------------------------------------------------


class _CellTable:
    def __init__(self, cells):
        self._cells = cells

    def item(self, row, col):
        text = self._cells.get((row, col))
        if text is None:
            return None
        return SimpleNamespace(text=lambda: text)


class TestRowBondAtoms:
    def _fn(self):
        return extract_function(BOND_EDITOR_PATH, "BondEditorWindow", "_row_bond_atoms")

    def _self(self, cells):
        s = SimpleNamespace(table=_CellTable(cells))
        s.COL_A1, s.COL_A2 = 1, 2
        return s

    def test_parses_labels(self):
        fn = self._fn()
        self_ = self._self({(0, 1): "3 (C)", (0, 2): "12 (Ag*)"})
        assert fn(self_, 0) == (3, 12)

    def test_missing_item_returns_none(self):
        fn = self._fn()
        assert fn(self._self({}), 0) is None

    def test_garbage_text_returns_none(self):
        fn = self._fn()
        self_ = self._self({(0, 1): "abc", (0, 2): "1 (C)"})
        assert fn(self_, 0) is None


class _SigBond:
    def __init__(self, b, e, t):
        self._b, self._e, self._t = b, e, t

    def GetBeginAtomIdx(self):
        return self._b

    def GetEndAtomIdx(self):
        return self._e

    def GetBondType(self):
        return self._t


class _SigMol:
    def __init__(self, n_atoms, bonds, coords):
        self._n = n_atoms
        self._bonds = bonds
        self._coords = coords

    def GetNumAtoms(self):
        return self._n

    def GetNumBonds(self):
        return len(self._bonds)

    def GetBonds(self):
        return self._bonds

    def GetNumConformers(self):
        return 1

    def GetConformer(self):
        return SimpleNamespace(
            GetPositions=lambda: real_numpy.array(self._coords, dtype=float)
        )


class TestBondSignature:
    def _fn(self):
        return extract_function(
            BOND_EDITOR_PATH,
            "BondEditorWindow",
            "get_mol_signature",
            extra_globals={"np": real_numpy},
        )

    def test_none_mol(self):
        assert self._fn()(SimpleNamespace(), None) is None

    def test_bond_type_change_changes_signature(self):
        fn = self._fn()
        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
        m1 = _SigMol(2, [_SigBond(0, 1, "SINGLE")], coords)
        m2 = _SigMol(2, [_SigBond(0, 1, "DOUBLE")], coords)
        s1, s2 = fn(SimpleNamespace(), m1), fn(SimpleNamespace(), m2)
        # ignore the id(mol) element (position 0)
        assert s1[1:] != s2[1:]

    def test_same_content_same_tail(self):
        fn = self._fn()
        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
        m1 = _SigMol(2, [_SigBond(0, 1, "SINGLE")], coords)
        m2 = _SigMol(2, [_SigBond(0, 1, "SINGLE")], coords)
        assert fn(SimpleNamespace(), m1)[1:] == fn(SimpleNamespace(), m2)[1:]
