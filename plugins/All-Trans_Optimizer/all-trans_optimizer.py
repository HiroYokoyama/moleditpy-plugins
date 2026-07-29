from PyQt6.QtWidgets import QMessageBox
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

PLUGIN_VERSION = "2026.07.30"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Convert non-cyclic chain torsions (including heteroatoms such as O, N, S) "
    "to an all-trans conformation."
)
PLUGIN_NAME = "All-Trans Optimizer"

# Backbone atoms allowed in the chain being straightened.
# C, N, O, F(as terminal only), P, S, Si, Cl/Br/I are common in organic backbones.
# The central (rotated) atoms are restricted to elements that actually form chains;
# terminal reference atoms may be any heavy atom.
_CENTRAL_ELEMENTS = "#6,#7,#8,#15,#16,#14"

# Central atom: a chain-forming element that is NOT part of a double/triple bond
# (excludes carbonyls, imines, nitro, etc. so we only rotate genuine single bonds).
_CENTRAL_ATOM = f"[{_CENTRAL_ELEMENTS};!$([*]=*);!$([*]#*)]"

# terminal - central -!@ central - terminal
# The central bond must be a single, acyclic bond (a rotatable backbone bond).
_ALL_TRANS_SMARTS = f"[!#1]-{_CENTRAL_ATOM}-;!@{_CENTRAL_ATOM}-[!#1]"


def _main_chain_atoms(mol):
    """Atom indices on the longest heavy-atom path (the backbone).

    Uses the topological distance matrix rather than all-pairs shortest
    paths, so this stays O(n^2) on large molecules.
    """
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if len(heavy) < 2:
        return set(heavy)

    dmat = Chem.GetDistanceMatrix(mol)
    best_pair, best_dist = (heavy[0], heavy[1]), -1.0
    for pos, i in enumerate(heavy):
        for j in heavy[pos + 1 :]:
            d = dmat[i][j]
            # Disconnected atoms get a huge sentinel distance, not a real path.
            if d > best_dist and d < mol.GetNumAtoms():
                best_dist, best_pair = d, (i, j)

    return set(Chem.rdmolops.GetShortestPath(mol, best_pair[0], best_pair[1]))


def _select_torsions(matches, backbone=None):
    """Keep exactly one torsion quartet per rotatable central bond.

    A branched backbone atom produces several matches that share the same
    central bond (idx2-idx3). Applying all of them would make each
    SetDihedralDeg undo the previous one for that bond, so we keep only one
    quartet per central bond.

    When ``backbone`` is given, the quartet whose reference atoms both lie on
    it wins. Taking the first match instead straightened a substituent and
    left the main chain bent -- 3-ethylnonane kept a 62 degree backbone
    dihedral. Ties keep the earlier match, so the result stays deterministic.
    """
    best = {}
    order = []
    for match in matches:
        idx1, idx2, idx3, idx4 = match
        bond_key = (idx2, idx3) if idx2 <= idx3 else (idx3, idx2)
        if backbone is None:
            score = 0
        else:
            score = (idx1 in backbone) + (idx4 in backbone)
        if bond_key not in best:
            best[bond_key] = (score, match)
            order.append(bond_key)
        elif score > best[bond_key][0]:
            best[bond_key] = (score, match)
    return [best[key][1] for key in order]


def run_plugin(context):
    """Straighten the acyclic chains of the current molecule to all-trans.

    Supports carbon chains as well as chains containing heteroatoms
    (ethers, alcohols, amines, thioethers, phosphates, ...).
    """
    mw = context.get_main_window()
    mol = context.current_mol

    if not mol:
        QMessageBox.warning(mw, PLUGIN_NAME, "No molecule loaded.")
        return

    try:
        # Require existing 3D coordinates; do not fabricate them.
        if mol.GetNumConformers() == 0:
            QMessageBox.warning(mw, PLUGIN_NAME, "Molecule has no 3D coordinates.")
            return

        conf = mol.GetConformer()

        patt = Chem.MolFromSmarts(_ALL_TRANS_SMARTS)
        matches = mol.GetSubstructMatches(patt)
        torsions = _select_torsions(matches, _main_chain_atoms(mol))

        if not torsions:
            QMessageBox.information(
                mw, PLUGIN_NAME, "No rotatable chain torsions found."
            )
            return

        count = 0
        for idx1, idx2, idx3, idx4 in torsions:
            # Set the backbone dihedral to 180 degrees (anti / trans).
            rdMolTransforms.SetDihedralDeg(conf, idx1, idx2, idx3, idx4, 180.0)
            count += 1

        # Push updated molecule back so 3D view redraws with new coordinates
        context.current_mol = mol
        context.refresh_3d_view()

        # Push undo state via V3 API
        context.push_undo_checkpoint()

        context.show_status_message(f"Applied All-Trans to {count} torsions.")
        QMessageBox.information(
            mw, PLUGIN_NAME, f"Applied All-Trans to {count} torsions."
        )

    except Exception as e:
        QMessageBox.critical(mw, PLUGIN_NAME, f"Error: {str(e)}")


def initialize(context):
    context.add_menu_action("3D Edit/All-Trans Optimizer", lambda: run_plugin(context))
