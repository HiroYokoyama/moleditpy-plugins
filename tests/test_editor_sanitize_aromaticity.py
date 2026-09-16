"""Regression tests for the three table editors' sanitize fallback.

An edit that breaks a ring makes ``Chem.SanitizeMol`` raise. The editors used
to swallow that and fall back to ``UpdatePropertyCache(strict=False)``, which
does not clear ``IsAromatic`` — so atoms stayed flagged aromatic while no
longer being in a ring. The molecule committed looking valid (the 3D view and
the editor tables read fine) and only failed later, when Optimize kekulizes it
in ``Chem.MolToMolBlock`` or MMFF atom typing.

These tests exercise the real ``sanitize_or_clear_aromaticity`` from each
plugin against real RDKit, then push the result through the same calls the host
makes on Optimize.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from conftest import extract_function

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")
from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
PLUGINS = {
    "xyz": ROOT / "plugins" / "XYZ_Editor" / "xyz_editor.py",
    "bond": ROOT / "plugins" / "Bond_Editor" / "bond_editor.py",
    "charge": ROOT / "plugins" / "Charge_Editor" / "charge_editor.py",
}
FUNC = "sanitize_or_clear_aromaticity"

HELPERS = {
    name: extract_function(
        path,
        None,
        FUNC,
        extra_globals={"Chem": Chem, "logging": logging, "PLUGIN_NAME": path.stem},
    )
    for name, path in PLUGINS.items()
}
ALL = pytest.mark.parametrize("helper_name", sorted(HELPERS))


def _naphthalene():
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccc2ccccc2c1"))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    return mol


def _stale_aromatics(mol):
    """Atoms still flagged aromatic that are no longer in any ring."""
    return [
        a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic() and not a.IsInRing()
    ]


def _optimize_as_host_does(mol):
    """The host's Optimize path: MolToMolBlock -> MolFromMolBlock -> MMFF."""
    mol_block = Chem.MolToMolBlock(mol, includeStereo=True)
    reparsed = Chem.MolFromMolBlock(mol_block, removeHs=False)
    assert reparsed is not None, "Failed to parse MOL block."
    return reparsed


# Edits that break a ring and so make the first SanitizeMol raise.
RING_BREAKING = {
    "delete_ring_bond": lambda rw: rw.RemoveBond(0, 1),
    "delete_ring_atom": lambda rw: rw.RemoveAtom(0),
}


@ALL
@pytest.mark.parametrize("edit_name", sorted(RING_BREAKING))
def test_ring_breaking_edit_leaves_no_stale_aromatic_flag(helper_name, edit_name):
    rw = Chem.RWMol(_naphthalene())
    RING_BREAKING[edit_name](rw)
    HELPERS[helper_name](rw)
    assert _stale_aromatics(rw.GetMol()) == []


@ALL
@pytest.mark.parametrize("edit_name", sorted(RING_BREAKING))
def test_ring_breaking_edit_still_optimizes(helper_name, edit_name):
    rw = Chem.RWMol(_naphthalene())
    RING_BREAKING[edit_name](rw)
    HELPERS[helper_name](rw)
    # Before the fix this raised AtomKekulizeException("non-ring atom 0
    # marked aromatic") inside MolToMolBlock, on the GUI thread.
    _optimize_as_host_does(rw.GetMol())


@ALL
def test_unkekulizable_edit_still_optimizes(helper_name):
    """Retyping a ring carbon leaves a ring RDKit cannot kekulize."""
    rw = Chem.RWMol(_naphthalene())
    rw.GetAtomWithIdx(0).SetAtomicNum(7)
    HELPERS[helper_name](rw)
    mol = _optimize_as_host_does(rw.GetMol())
    assert AllChem.MMFFGetMoleculeProperties(mol) is not None


@ALL
def test_intact_aromatic_ring_keeps_its_flags(helper_name):
    """A clean molecule must sanitize on the first pass, untouched."""
    mol = _naphthalene()
    aromatic_before = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    rw = Chem.RWMol(mol)
    HELPERS[helper_name](rw)
    out = rw.GetMol()
    assert sum(1 for a in out.GetAtoms() if a.GetIsAromatic()) == aromatic_before
    _optimize_as_host_does(out)


@ALL
def test_surviving_ring_stays_aromatic_after_the_other_is_broken(helper_name):
    """Clearing flags is a retry step, not a permanent de-aromatization."""
    rw = Chem.RWMol(_naphthalene())
    rw.RemoveBond(0, 1)
    HELPERS[helper_name](rw)
    out = rw.GetMol()
    assert sum(1 for a in out.GetAtoms() if a.GetIsAromatic()) == 6


@ALL
def test_hopeless_molecule_does_not_raise(helper_name):
    """A valence RDKit cannot fix at all must still fall through quietly."""
    rw = Chem.RWMol(Chem.MolFromSmiles("C", sanitize=False))
    rw.GetAtomWithIdx(0).SetNumExplicitHs(9)
    HELPERS[helper_name](rw)  # must not raise
