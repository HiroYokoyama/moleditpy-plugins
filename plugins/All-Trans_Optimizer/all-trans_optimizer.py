from PyQt6.QtWidgets import QMessageBox
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

PLUGIN_VERSION = "2026.07.11"
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


def _select_torsions(matches):
    """Keep exactly one torsion quartet per rotatable central bond.

    A branched backbone atom produces several matches that share the same
    central bond (idx2-idx3). Applying all of them would make each
    SetDihedralDeg undo the previous one for that bond, so we keep only the
    first quartet seen for each central bond. Order is preserved for
    determinism.
    """
    seen_bonds = set()
    selected = []
    for match in matches:
        idx1, idx2, idx3, idx4 = match
        bond_key = (idx2, idx3) if idx2 <= idx3 else (idx3, idx2)
        if bond_key in seen_bonds:
            continue
        seen_bonds.add(bond_key)
        selected.append(match)
    return selected


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
        torsions = _select_torsions(matches)

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
    context.add_menu_action(
        "3D Edit/All-Trans Optimizer", lambda: run_plugin(context)
    )
