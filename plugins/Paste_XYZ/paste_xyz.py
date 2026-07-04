import logging
import traceback

from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QTextEdit,
    QPushButton,
    QLabel,
    QLineEdit,
    QDialogButtonBox,
    QMessageBox,
)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    Chem = None

PLUGIN_NAME = "Paste XYZ"
PLUGIN_VERSION = "2026.07.04"
PLUGIN_CATEGORY = "File"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Paste XYZ coordinates from the clipboard to create a new molecule."


def initialize(context):
    context.add_menu_action("File/Paste XYZ...", lambda: run_plugin(context))


# ---------------------------------------------------------------------------
# Pure helpers (module-level for testability)
# ---------------------------------------------------------------------------

def parse_xyz_lines(text: str) -> list:
    """Return [(symbol, x, y, z), ...] from XYZ text, ignoring headers/invalid lines."""
    atoms = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        symbol = parts[0]
        if not symbol[0].isalpha():
            continue
        atoms.append((symbol, x, y, z))
    return atoms


def _build_rwmol(atoms_data: list) -> "Chem.RWMol":
    mol = Chem.RWMol()
    for i, (symbol, _, __, ___) in enumerate(atoms_data):
        atom = Chem.Atom(symbol)
        atom.SetIntProp("xyz_unique_id", i)
        mol.AddAtom(atom)
    conf = Chem.Conformer(len(atoms_data))
    for i, (_, x, y, z) in enumerate(atoms_data):
        conf.SetAtomPosition(i, AllChem.rdGeometry.Point3D(x, y, z))
    mol.AddConformer(conf)
    return mol


def _apply_bonds(rwmol: "Chem.RWMol", charge: int, mw) -> "Chem.Mol":
    """Try rdDetermineBonds; fall back to distance-based estimation."""
    try:
        from rdkit.Chem import rdDetermineBonds
        candidate = Chem.RWMol(rwmol)
        rdDetermineBonds.DetermineBonds(candidate, charge=charge)
        result = candidate.GetMol()
        result.SetIntProp("_xyz_charge", int(charge))
        return result
    except Exception:
        pass

    if hasattr(mw, "io_manager") and hasattr(mw.io_manager, "estimate_bonds_from_distances"):
        try:
            mw.io_manager.estimate_bonds_from_distances(rwmol)
        except Exception as exc:
            logging.warning("[paste_xyz] estimate_bonds_from_distances: %s", exc)
    result = rwmol.GetMol()
    result.SetIntProp("_xyz_charge", int(charge))
    return result


# ---------------------------------------------------------------------------
# Charge prompt dialog
# ---------------------------------------------------------------------------

def _prompt_charge(mw) -> tuple:
    """Return (charge: int, ok: bool, skip: bool)."""
    dialog = QDialog(mw)
    dialog.setWindowTitle("Import XYZ — Charge")
    layout = QVBoxLayout(dialog)
    layout.addWidget(QLabel("Enter total molecular charge:"))
    line_edit = QLineEdit(dialog)
    line_edit.setText("0")
    layout.addWidget(line_edit)

    btn_box = QDialogButtonBox(
        QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
        parent=dialog,
    )
    skip_btn = QPushButton("Skip chemistry", dialog)
    hl = QHBoxLayout()
    hl.addWidget(btn_box)
    hl.addWidget(skip_btn)
    layout.addLayout(hl)

    result = {"accepted": False, "skip": False}

    def _on_ok():
        result["accepted"] = True
        dialog.accept()

    def _on_skip():
        result["skip"] = True
        dialog.accept()

    btn_box.accepted.connect(_on_ok)
    btn_box.rejected.connect(dialog.reject)
    skip_btn.clicked.connect(_on_skip)

    if dialog.exec() != QDialog.DialogCode.Accepted:
        return 0, False, False
    if result["skip"]:
        return 0, True, True
    try:
        return int(line_edit.text().strip()), True, False
    except ValueError:
        return 0, True, False


# ---------------------------------------------------------------------------
# Main dialog
# ---------------------------------------------------------------------------

class PasteXYZDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Paste XYZ")
        self.resize(600, 400)
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("Paste XYZ coordinates below (atom count / comment lines are ignored)."))

        self.text_edit = QTextEdit()
        self.text_edit.setPlaceholderText(
            "Paste XYZ data here...\nExample:\nC  0.000  0.000  0.000\nH  1.089  0.000  0.000\n..."
        )
        layout.addWidget(self.text_edit)

        btn_layout = QHBoxLayout()
        self.load_btn = QPushButton("Load")
        self.cancel_btn = QPushButton("Cancel")
        btn_layout.addStretch()
        btn_layout.addWidget(self.load_btn)
        btn_layout.addWidget(self.cancel_btn)
        layout.addLayout(btn_layout)

        self.load_btn.clicked.connect(self.accept)
        self.cancel_btn.clicked.connect(self.reject)

    def get_data(self) -> str:
        return self.text_edit.toPlainText()


# ---------------------------------------------------------------------------
# Plugin action
# ---------------------------------------------------------------------------

def _resolve_mol_with_charge(rwmol, context, mw) -> "Chem.Mol | None":
    """Run bond determination, prompting for charge when needed."""
    always_ask = context.get_setting("always_ask_charge", False)

    if not always_ask:
        try:
            return _apply_bonds(rwmol, 0, mw)
        except Exception:
            pass

    while True:
        charge, ok, skip = _prompt_charge(mw)
        if not ok:
            return None
        if skip:
            if hasattr(mw, "io_manager") and hasattr(mw.io_manager, "estimate_bonds_from_distances"):
                try:
                    mw.io_manager.estimate_bonds_from_distances(rwmol)
                except Exception as exc:
                    logging.warning("[paste_xyz] estimate_bonds_from_distances: %s", exc)
            result = rwmol.GetMol()
            result.SetIntProp("_xyz_skip_checks", 1)
            return result
        try:
            return _apply_bonds(rwmol, charge, mw)
        except Exception:
            context.show_status_message("Bond determination failed for that charge — try again.", 3000)


def run_plugin(context) -> None:
    mw = context.get_main_window()
    if Chem is None:
        QMessageBox.critical(mw, "Error", "RDKit is not available.")
        return

    dialog = PasteXYZDialog(mw)
    if dialog.exec() != QDialog.DialogCode.Accepted:
        return

    xyz_text = dialog.get_data()
    if not xyz_text.strip():
        return

    try:
        atoms_data = parse_xyz_lines(xyz_text)
        if not atoms_data:
            QMessageBox.warning(
                mw,
                "Paste XYZ",
                "No valid coordinate lines found.\nExpected format: Symbol X Y Z",
            )
            return

        context.clear_canvas(push_to_undo=False)
        rwmol = _build_rwmol(atoms_data)

        skip_checks = context.get_setting("skip_chemistry_checks", False)
        if skip_checks:
            if hasattr(mw, "io_manager") and hasattr(mw.io_manager, "estimate_bonds_from_distances"):
                try:
                    mw.io_manager.estimate_bonds_from_distances(rwmol)
                except Exception as exc:
                    logging.warning("[paste_xyz] estimate_bonds_from_distances: %s", exc)
            final_mol = rwmol.GetMol()
            final_mol.SetIntProp("_xyz_skip_checks", 1)
        else:
            final_mol = _resolve_mol_with_charge(rwmol, context, mw)

        if final_mol is None:
            return

        context.current_molecule = final_mol
        context.push_undo_checkpoint()
        context.check_chemistry_problems()
        context.refresh_ui()
        context.enter_3d_mode()
        context.fit_3d_view()
        context.show_status_message(f"Pasted {len(atoms_data)} atoms from clipboard.")

    except Exception as exc:
        logging.exception("Paste XYZ: failed to parse or load: %s", exc)
        QMessageBox.critical(mw, "Error", f"Failed to parse or load data:\n{exc}")
