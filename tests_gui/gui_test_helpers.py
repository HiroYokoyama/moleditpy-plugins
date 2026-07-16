"""
Shared fakes/helpers used by more than one plugin's GUI test file.

Currently used by:
  - test_gui_plugin_bond_editor.py
  - test_gui_plugin_charge_editor.py
"""

from __future__ import annotations

from unittest.mock import MagicMock


def ctx_no_mol() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.plotter = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


class FakeAtom:
    def __init__(self, idx, symbol="C", charge=0, radicals=0, hs=0):
        self._idx = idx
        self._symbol = symbol
        self._charge = charge
        self._radicals = radicals
        self._hs = hs

    def GetIdx(self):
        return self._idx

    def GetSymbol(self):
        return self._symbol

    def HasProp(self, name):
        return False

    def GetFormalCharge(self):
        return self._charge

    def GetNumRadicalElectrons(self):
        return self._radicals

    def GetTotalNumHs(self):
        return self._hs


class FakeBond:
    def __init__(self, idx, begin, end, bond_type="SINGLE"):
        self._idx = idx
        self._begin = begin
        self._end = end
        self._type = bond_type

    def GetIdx(self):
        return self._idx

    def GetBeginAtomIdx(self):
        return self._begin

    def GetEndAtomIdx(self):
        return self._end

    def GetBondType(self):
        return self._type


class FakeMol:
    """No conformer — keeps the plugins away from mocked numpy."""

    def __init__(self, atoms, bonds=()):
        self._atoms = list(atoms)
        self._bonds = list(bonds)

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetNumBonds(self):
        return len(self._bonds)

    def GetNumConformers(self):
        return 0

    def GetAtoms(self):
        return list(self._atoms)

    def GetBonds(self):
        return list(self._bonds)

    def GetAtomWithIdx(self, idx):
        return self._atoms[idx]
