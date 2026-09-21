"""Interactive one-attachment functional-group replacement tool."""
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QComboBox, QPushButton, QLabel, QMessageBox
from rdkit import Chem

PLUGIN_NAME = "Functional Group Toolbox"
PLUGIN_VERSION = "2026.09.21"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_SUPPORTED_PYTHON_VERSION = ">=3.9, <3.15"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Replace a selected atom with a common functional group while retaining its neighbours."
PLUGIN_DEPENDENCIES = ["rdkit", "PyQt6"]
GROUPS = {
    "Methyl": "[*:1]C", "Ethyl": "[*:1]CC", "Hydroxyl": "[*:1]O",
    "Amino": "[*:1]N", "Methoxy": "[*:1]OC", "Acetyl": "[*:1]C(=O)C",
    "Carboxyl": "[*:1]C(=O)O", "Cyano": "[*:1]C#N", "Nitro": "[*:1][N+](=O)[O-]",
    "Phenyl": "[*:1]c1ccccc1", "tert-Butyl": "[*:1]C(C)(C)C",
    "Fluoro": "[*:1]F", "Chloro": "[*:1]Cl", "Bromo": "[*:1]Br", "Iodo": "[*:1]I",
}

class FunctionalGroupToolbox(QWidget):
    def __init__(self, context):
        super().__init__(context.get_main_window())
        self.context = context
        self.setWindowTitle(PLUGIN_NAME)
        self.resize(300, 150)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Atom index to replace:"))
        self.atom_combo = QComboBox()
        layout.addWidget(self.atom_combo)
        layout.addWidget(QLabel("Functional group:"))
        self.group_combo = QComboBox()
        self.group_combo.addItems(GROUPS)
        layout.addWidget(self.group_combo)
        button = QPushButton("Replace selected atom")
        button.clicked.connect(self.replace_atom)
        layout.addWidget(button)
        self.refresh_atoms()
        context.register_window("functional_group_toolbox", self)

    def refresh_atoms(self):
        self.atom_combo.clear()
        mol = self.context.current_mol
        if mol:
            self.atom_combo.addItems([f"{a.GetIdx()}: {a.GetSymbol()}" for a in mol.GetAtoms()])

    def replace_atom(self):
        mol = self.context.current_mol
        target = self.atom_combo.currentIndex()
        if mol is None or target < 0:
            return
        fragment = Chem.MolFromSmiles(GROUPS[self.group_combo.currentText()])
        if fragment is None:
            return
        dummy = next((a.GetIdx() for a in fragment.GetAtoms() if a.GetAtomicNum() == 0), None)
        if dummy is None:
            return
        attachment_atoms = list(fragment.GetAtomWithIdx(dummy).GetNeighbors())
        if len(attachment_atoms) != 1:
            QMessageBox.warning(self, "Replacement failed", "The selected group must have exactly one attachment point.")
            return
        attach = attachment_atoms[0].GetIdx()
        try:
            rw = Chem.RWMol(mol)
            target_atom = rw.GetAtomWithIdx(target)
            source_atom = fragment.GetAtomWithIdx(attach)
            target_atom.SetAtomicNum(source_atom.GetAtomicNum())
            target_atom.SetFormalCharge(source_atom.GetFormalCharge())
            target_atom.SetIsAromatic(source_atom.GetIsAromatic())
            mapping = {attach: target}
            for atom in fragment.GetAtoms():
                if atom.GetIdx() != dummy and atom.GetIdx() != attach:
                    mapping[atom.GetIdx()] = rw.AddAtom(Chem.Atom(atom))
            for bond in fragment.GetBonds():
                a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if dummy in (a, b):
                    continue
                if rw.GetBondBetweenAtoms(mapping[a], mapping[b]) is None:
                    rw.AddBond(mapping[a], mapping[b], bond.GetBondType())
            Chem.SanitizeMol(rw)
            new_mol = rw.GetMol()
            self.context.current_molecule = new_mol
            self.context.push_undo_checkpoint()
            refresh = getattr(self.context, "refresh_3d_view", None)
            if callable(refresh):
                refresh()
            self.refresh_atoms()
            self.context.show_status_message("Functional group replacement applied.")
        except (RuntimeError, ValueError, AttributeError) as exc:
            QMessageBox.critical(self, "Replacement failed", str(exc))


def initialize(context):
    def show():
        window = context.get_window("functional_group_toolbox")
        if window is None:
            window = FunctionalGroupToolbox(context)
        window.show(); window.raise_(); window.activateWindow()
    context.add_menu_action("Edit/Functional Group Toolbox...", show)
