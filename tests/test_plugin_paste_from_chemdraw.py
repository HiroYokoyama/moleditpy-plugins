"""
Tests for the Paste from ChemDraw plugin: reconstruct_from_flat_text
(module-level function).
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CHEMDRAW_PATH = PLUGINS_DIR / "Paste_from_ChemDraw" / "paste_chemdraw.py"

with mock_optional_imports():
    _chemdraw = load_plugin(CHEMDRAW_PATH)


class TestReconstructFromFlatText:
    def test_empty_string_returns_none(self):
        assert _chemdraw.reconstruct_from_flat_text("") is None

    def test_no_v2000_marker_returns_none(self):
        assert _chemdraw.reconstruct_from_flat_text("some ordinary text") is None

    def test_too_few_header_tokens_returns_none(self):
        # V2000 present but fewer than 10 header tokens before it
        result = _chemdraw.reconstruct_from_flat_text("only three tokens V2000 rest")
        assert result is None

    def test_non_printable_chars_dont_raise(self):
        # Control characters should be sanitised; no V2000 → None
        result = _chemdraw.reconstruct_from_flat_text("\x01\x02\x16garbage")
        assert result is None




# ---------------------------------------------------------------------------
# reconstruct_from_flat_text MolBlock output
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# run() drawing loop: atom items must come from scene.atom_items, not
# scene.data.atoms[id]["item"] (regression test for KeyError: 'item')
# ---------------------------------------------------------------------------


class _FakeQPointF:
    def __init__(self, x, y):
        self._x, self._y = x, y

    def x(self):
        return self._x

    def y(self):
        return self._y


class _FakeAtom:
    def __init__(self, idx, symbol, charge=0):
        self.idx, self.symbol, self.charge = idx, symbol, charge

    def GetSymbol(self):
        return self.symbol

    def GetFormalCharge(self):
        return self.charge


class _FakeBond:
    def __init__(self, begin_idx, end_idx, order=1.0):
        self.begin_idx, self.end_idx, self.order = begin_idx, end_idx, order

    def GetBeginAtomIdx(self):
        return self.begin_idx

    def GetEndAtomIdx(self):
        return self.end_idx

    def GetBondTypeAsDouble(self):
        return self.order

    def GetBondDir(self):
        return None


class _FakeConf:
    def __init__(self, coords):
        self._coords = coords

    def GetAtomPosition(self, idx):
        return self._coords[idx]


class _FakeMol:
    """Two carbons, single bond -- just enough to drive run()'s draw loop."""

    def __init__(self):
        self._atoms = [_FakeAtom(0, "C"), _FakeAtom(1, "C")]
        self._bonds = [_FakeBond(0, 1, 1.0)]
        self._conf = _FakeConf({0: type("P", (), {"x": 0.0, "y": 0.0})(),
                                 1: type("P", (), {"x": 1.4, "y": 0.0})()})

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetAtomWithIdx(self, i):
        return self._atoms[i]

    def GetBonds(self):
        return list(self._bonds)

    def GetNumConformers(self):
        return 1

    def GetConformer(self):
        return self._conf


class _FakeAtomItem:
    """Mirrors the real AtomItem: has .atom_id, nothing else needed here."""

    def __init__(self, atom_id):
        self.atom_id = atom_id


class _FakeSceneData:
    """Mirrors MolecularData: atom dicts have no "item" key."""

    def __init__(self):
        self.atoms = {}


class _FakeScene:
    """Mirrors the real molecule scene's create_atom/create_bond contract."""

    def __init__(self):
        self.data = _FakeSceneData()
        self.atom_items = {}
        self.bonds_created = []
        self._next_id = 0

    def create_atom(self, symbol, pos, charge=0):
        atom_id = self._next_id
        self._next_id += 1
        self.data.atoms[atom_id] = {
            "symbol": symbol,
            "pos": (pos.x(), pos.y()),
            "charge": charge,
            "radical": 0,
        }
        self.atom_items[atom_id] = _FakeAtomItem(atom_id)
        return atom_id

    def create_bond(self, atom1, atom2, bond_order=None, bond_stereo=None):
        self.bonds_created.append((atom1.atom_id, atom2.atom_id, bond_order, bond_stereo))


class TestChemDrawDrawingLoop:
    def test_atom_items_come_from_scene_atom_items_not_data_dict(self):
        """
        Regression test: run() must resolve pasted-atom Qt items via
        scene.atom_items[id], never scene.data.atoms[id]["item"] (that key
        does not exist on the real MolecularData atom dict and raised
        KeyError: 'item' in production).
        """
        scene = _FakeScene()
        mol = _FakeMol()

        # Exercise exactly the code path run() uses after a mol is parsed:
        # build rdkit_idx_to_item via scene.atom_items, then create bonds.
        rdkit_idx_to_item = {}
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            a_id = scene.create_atom(
                atom.GetSymbol(), _FakeQPointF(0, 0), charge=atom.GetFormalCharge()
            )
            # This is the line that was buggy: scene.data.atoms[a_id]["item"]
            # would raise KeyError('item') against a real-shaped data dict.
            rdkit_idx_to_item[i] = scene.atom_items[a_id]

        for bond in mol.GetBonds():
            idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            atom1, atom2 = rdkit_idx_to_item[idx1], rdkit_idx_to_item[idx2]
            scene.create_bond(atom1, atom2, bond_order=int(bond.GetBondTypeAsDouble()))

        assert scene.bonds_created == [(0, 1, 1, None)]
        assert "item" not in scene.data.atoms[0]

    def test_source_uses_atom_items_not_data_atoms_item_key(self):
        """
        Static guard: the specific KeyError-causing expression must not
        reappear in the source even if the surrounding code is refactored.
        """
        source = CHEMDRAW_PATH.read_text(encoding="utf-8")
        assert 'data.atoms[a_id]["item"]' not in source
        assert "context.scene.atom_items[a_id]" in source
