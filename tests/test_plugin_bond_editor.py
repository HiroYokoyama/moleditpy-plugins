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


def _add_bond_self(mol, type_label="Single"):
    committed = []
    return (
        SimpleNamespace(
            context=SimpleNamespace(
                current_molecule=mol, show_status_message=MagicMock()
            ),
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
        self_, committed = _add_bond_self(_bond_mol(3), "Double")
        fn(self_, 0, 2)
        assert committed[0][0].added == [(0, 2, "DOUBLE")]

    def test_self_bond_refused(self):
        fn = _add_bond_fn()
        self_, committed = _add_bond_self(_bond_mol(3))
        fn(self_, 1, 1)
        assert committed == []
        self_.context.show_status_message.assert_called_once()

    def test_existing_bond_refused(self):
        fn = _add_bond_fn()
        self_, committed = _add_bond_self(_bond_mol(3, [(0, 1)]))
        fn(self_, 1, 0)
        assert committed == []

    def test_no_molecule_refused(self):
        fn = _add_bond_fn()
        self_, committed = _add_bond_self(None)
        fn(self_, 0, 1)
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


class TestModeOverlay:
    def _fn(self):
        import logging as real_logging

        return extract_function(
            BOND_EDITOR_PATH,
            "BondEditorWindow",
            "_update_mode_overlay",
            extra_globals={"logging": real_logging},
        )

    @staticmethod
    def _self(mode, first_idx=None, plotter=None, type_label="Single"):
        return SimpleNamespace(
            context=SimpleNamespace(plotter=plotter),
            click_mode_combo=SimpleNamespace(currentText=lambda: mode),
            add_type_combo=SimpleNamespace(currentText=lambda: type_label),
            _first_pick_idx=first_idx,
        )

    def test_select_mode_removes_label(self):
        plotter = MagicMock()
        self._fn()(self._self("Select bond", plotter=plotter))
        plotter.remove_actor.assert_called_once_with("bond_editor_mode_label")
        plotter.add_text.assert_not_called()
        plotter.render.assert_called_once()

    def test_create_bond_first_click_prompt_includes_type(self):
        plotter = MagicMock()
        self._fn()(self._self("Create bond", plotter=plotter, type_label="Double"))
        text = plotter.add_text.call_args[0][0]
        assert "first atom" in text
        assert "double" in text
        assert plotter.add_text.call_args[1]["name"] == "bond_editor_mode_label"

    def test_create_bond_second_click_shows_picked_atom(self):
        plotter = MagicMock()
        self._fn()(self._self("Create bond", first_idx=7, plotter=plotter))
        text = plotter.add_text.call_args[0][0]
        assert "atom 7" in text
        assert "second atom" in text

    def test_no_plotter_is_noop(self):
        self._fn()(self._self("Create bond", plotter=None))


class TestPickedAtomLabels:
    def _fn(self):
        import logging as real_logging

        return extract_function(
            BOND_EDITOR_PATH,
            "BondEditorWindow",
            "_update_picked_atom_labels",
            extra_globals={"logging": real_logging},
        )

    @staticmethod
    def _mol(n=2):
        return SimpleNamespace(
            GetNumAtoms=lambda: n,
            GetNumConformers=lambda: 1,
            GetConformer=lambda: SimpleNamespace(
                GetAtomPosition=lambda i: SimpleNamespace(x=float(i), y=0.0, z=0.0)
            ),
        )

    def _self(self, mol, picked, plotter):
        return SimpleNamespace(
            context=SimpleNamespace(plotter=plotter, current_mol=mol),
            _picked_atoms=picked,
            _atom_label=lambda m, i: f"{i} (C)",
        )

    def test_picked_atom_gets_label(self):
        plotter = MagicMock()
        self._fn()(self._self(self._mol(3), {"Atom 1": 2}, plotter))
        args, kwargs = plotter.add_point_labels.call_args
        assert args[0] == [[2.0, 0.0, 0.0]]
        assert args[1] == ["Atom 1: 2 (C)"]
        assert kwargs["name"] == "bond_editor_atom_labels"
        assert kwargs["pickable"] is False
        assert kwargs["reset_camera"] is False

    def test_no_pick_removes_labels(self):
        plotter = MagicMock()
        self._fn()(self._self(self._mol(2), {}, plotter))
        plotter.remove_actor.assert_called_once_with("bond_editor_atom_labels")
        plotter.add_point_labels.assert_not_called()

    def test_stale_index_ignored(self):
        plotter = MagicMock()
        self._fn()(self._self(self._mol(2), {"Atom 1": 9}, plotter))
        plotter.remove_actor.assert_called_once_with("bond_editor_atom_labels")
        plotter.add_point_labels.assert_not_called()

    def test_no_molecule_removes_labels(self):
        plotter = MagicMock()
        self._fn()(self._self(None, {"Atom 1": 0}, plotter))
        plotter.remove_actor.assert_called_once_with("bond_editor_atom_labels")
        plotter.add_point_labels.assert_not_called()

    def test_no_plotter_is_noop(self):
        self._fn()(self._self(self._mol(1), {"Atom 1": 0}, None))


class TestCreateBondPick:
    def _fn(self):
        return extract_function(
            BOND_EDITOR_PATH, "BondEditorWindow", "_create_bond_pick"
        )

    @staticmethod
    def _self(first_idx=None):
        return SimpleNamespace(
            _first_pick_idx=first_idx,
            _picked_atoms={},
            context=SimpleNamespace(show_status_message=MagicMock()),
            add_bond=MagicMock(),
            _update_mode_overlay=MagicMock(),
            _update_picked_atom_labels=MagicMock(),
        )

    def test_first_click_stores_pick(self):
        self_ = self._self()
        self._fn()(self_, 3)
        assert self_._first_pick_idx == 3
        assert self_._picked_atoms == {"Atom 1": 3}
        self_.add_bond.assert_not_called()
        self_._update_picked_atom_labels.assert_called_once()

    def test_second_click_creates_bond(self):
        self_ = self._self(first_idx=3)
        self_._picked_atoms = {"Atom 1": 3}
        self._fn()(self_, 5)
        self_.add_bond.assert_called_once_with(3, 5)
        assert self_._first_pick_idx is None
        assert self_._picked_atoms == {}

    def test_same_atom_click_cancels(self):
        self_ = self._self(first_idx=3)
        self_._picked_atoms = {"Atom 1": 3}
        self._fn()(self_, 3)
        self_.add_bond.assert_not_called()
        assert self_._first_pick_idx is None
        assert self_._picked_atoms == {}
