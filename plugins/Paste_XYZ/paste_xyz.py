
import traceback
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QTextEdit, QPushButton, QHBoxLayout, 
    QMessageBox, QLabel, QLineEdit, QDialogButtonBox
)
import logging
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    Chem = None

PLUGIN_NAME = "Paste XYZ"
PLUGIN_VERSION = "2026.04.11"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Allows pasting XYZ coordinates directly from the clipboard to create a new molecule."

class PasteXYZDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Paste XYZ")
        self.resize(600, 400)
        self.layout = QVBoxLayout(self)

        info_label = QLabel("Paste XYZ coordinates below (Header/Atom count is IGNORED).")
        self.layout.addWidget(info_label)

        self.text_edit = QTextEdit()
        self.text_edit.setPlaceholderText("Paste XYZ data here...\nExample:\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\n...")
        self.layout.addWidget(self.text_edit)

        btn_layout = QHBoxLayout()
        self.load_btn = QPushButton("Load")
        self.cancel_btn = QPushButton("Cancel")
        btn_layout.addStretch()
        btn_layout.addWidget(self.load_btn)
        btn_layout.addWidget(self.cancel_btn)
        self.layout.addLayout(btn_layout)

        self.load_btn.clicked.connect(self.accept)
        self.cancel_btn.clicked.connect(self.reject)

    def get_data(self):
        return self.text_edit.toPlainText()

def run_plugin(context):
    mw = context.get_main_window()
    if Chem is None:
        QMessageBox.critical(mw, "Error", "RDKit is not available.")
        return

    dialog = PasteXYZDialog(mw)
    if dialog.exec() == QDialog.DialogCode.Accepted:
        xyz_text = dialog.get_data()
        if not xyz_text.strip():
            return

        try:
            # Robust Parsing Logic: Scan for 'Symbol X Y Z' lines
            atoms_data = []
            lines = xyz_text.splitlines()
            
            for line in lines:
                parts = line.split()
                if len(parts) >= 4:
                    # Check if last 3 are standard floats (coordinates)
                    try:
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        symbol = parts[0]
                        # Verify symbol is likely an element (alpha)
                        if not symbol[0].isalpha():
                                continue
                        
                        atoms_data.append((symbol, x, y, z))
                    except ValueError:
                        continue

            if not atoms_data:
                QMessageBox.warning(mw, "Paste XYZ", "No valid coordinate lines found.\nExpected format: Symbol X Y Z")
                return

            # Clean workspace
            # [DIRECT ACCESS] to main window for state reset
            if hasattr(mw, 'clear_all'):
                mw.edit_actions_manager.clear_all()
            
            # Create RWMol and add atoms/conformers
            mol = Chem.RWMol()
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                atom = Chem.Atom(symbol)
                atom.SetIntProp("xyz_unique_id", i)
                mol.AddAtom(atom)
            
            conf = Chem.Conformer(len(atoms_data))
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                conf.SetAtomPosition(i, AllChem.rdGeometry.Point3D(x, y, z))
            mol.AddConformer(conf)

            # Helper: prompt for charge
            def prompt_for_charge():
                dialog = QDialog(mw)
                dialog.setWindowTitle("Import XYZ Charge")
                layout = QVBoxLayout(dialog)
                label = QLabel("Enter total molecular charge:")
                line_edit = QLineEdit(dialog)
                line_edit.setText("0")
                btn_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, parent=dialog)
                skip_btn = QPushButton("Skip chemistry", dialog)
                hl = QHBoxLayout()
                hl.addWidget(btn_box)
                hl.addWidget(skip_btn)
                layout.addWidget(label)
                layout.addWidget(line_edit)
                layout.addLayout(hl)
                
                result = {"accepted": False, "skip": False}
                def on_ok():
                    result["accepted"] = True
                    dialog.accept()
                def on_skip():
                    result["skip"] = True
                    dialog.accept()

                btn_box.accepted.connect(on_ok)
                btn_box.rejected.connect(dialog.reject)
                skip_btn.clicked.connect(on_skip)

                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return None, False, False
                if result["skip"]:
                    return 0, True, True
                
                try:
                    return int(line_edit.text().strip()), True, False
                except:
                    return 0, True, False

            # Inner helper: process with charge
            def _process_with_charge(charge_val):
                used_rd_determine = False
                mol_to_finalize = None
                try:
                    from rdkit.Chem import rdDetermineBonds
                    mol_candidate = Chem.RWMol(mol)
                    rdDetermineBonds.DetermineBonds(mol_candidate, charge=charge_val)
                    mol_to_finalize = mol_candidate
                    used_rd_determine = True
                except:
                    used_rd_determine = False
                    mol_to_finalize = mol
                
                if not used_rd_determine:
                    if hasattr(mw, 'io_manager') and hasattr(mw.io_manager, 'estimate_bonds_from_distances'):
                        mw.io_manager.estimate_bonds_from_distances(mol_to_finalize)

                candidate_mol = mol_to_finalize.GetMol()
                candidate_mol.SetIntProp("_xyz_charge", int(charge_val))
                
                if hasattr(mw, 'edit_actions_manager') and hasattr(mw.edit_actions_manager, '_apply_chem_check_and_set_flags'):
                    mw.edit_actions_manager._apply_chem_check_and_set_flags(candidate_mol, source_desc='PasteXYZ')

                return candidate_mol

            # Main Logic Loop
            final_mol = None
            settings = getattr(mw.init_manager, 'settings', {}) if hasattr(mw, 'init_manager') else {}
            always_ask = bool(settings.get('always_ask_charge', False))
            skip_checks_global = bool(settings.get('skip_chemistry_checks', False))

            if skip_checks_global:
                if hasattr(mw, 'io_manager') and hasattr(mw.io_manager, 'estimate_bonds_from_distances'):
                    try: mw.io_manager.estimate_bonds_from_distances(mol)
                    except Exception as _e:
                        logging.warning("[paste_xyz.py:178] silenced: %s", _e)
                final_mol = mol.GetMol()
                final_mol.SetIntProp("_xyz_skip_checks", 1)
            else:
                try:
                    if not always_ask:
                        try:
                            final_mol = _process_with_charge(0)
                        except:
                            while True:
                                charge_val, ok, skip_flag = prompt_for_charge()
                                if not ok: return
                                if skip_flag:
                                    if hasattr(mw, 'io_manager') and hasattr(mw.io_manager, 'estimate_bonds_from_distances'):
                                        try: mw.io_manager.estimate_bonds_from_distances(mol)
                                        except Exception as _e:
                                            logging.warning("[paste_xyz.py:193] silenced: %s", _e)
                                    final_mol = mol.GetMol()
                                    final_mol.SetIntProp("_xyz_skip_checks", 1)
                                    break
                                try:
                                    final_mol = _process_with_charge(charge_val)
                                    break
                                except:
                                    context.show_status_message("Bond determination failed for that charge.")
                    else:
                        while True:
                            charge_val, ok, skip_flag = prompt_for_charge()
                            if not ok: return
                            if skip_flag:
                                if hasattr(mw, 'io_manager') and hasattr(mw.io_manager, 'estimate_bonds_from_distances'):
                                    try: mw.io_manager.estimate_bonds_from_distances(mol)
                                    except Exception as _e:
                                        logging.warning("[paste_xyz.py:209] silenced: %s", _e)
                                final_mol = mol.GetMol()
                                final_mol.SetIntProp("_xyz_skip_checks", 1)
                                break
                            try:
                                final_mol = _process_with_charge(charge_val)
                                break
                            except:
                                context.show_status_message("Bond determination failed.")
                except Exception as _e:
                    logging.warning("[paste_xyz.py:218] silenced: %s", _e)

            if final_mol is None:
                if hasattr(mw, 'io_manager') and hasattr(mw.io_manager, 'estimate_bonds_from_distances'):
                    try: mw.io_manager.estimate_bonds_from_distances(mol)
                    except Exception as _e:
                        logging.warning("[paste_xyz.py:224] silenced: %s", _e)
                final_mol = mol.GetMol()

            if final_mol:
                rw_mol = Chem.RWMol(final_mol)
                context.current_mol = rw_mol
                context.push_undo_checkpoint()
                context.show_status_message(f"Pasted {len(atoms_data)} atoms from clipboard.")
                
                if hasattr(mw, 'ui_manager') and hasattr(mw.ui_manager, '_enter_3d_viewer_ui_mode'):
                    try: mw.ui_manager._enter_3d_viewer_ui_mode()
                    except Exception as _e:
                        logging.warning("[paste_xyz.py:235] silenced: %s", _e)
                # Minimize the 2D editor after loading XYZ into 3D-first workflow.
                try:
                    if hasattr(mw, "ui_manager") and hasattr(mw.ui_manager, "minimize_2d_panel"):
                        mw.ui_manager.minimize_2d_panel()
                    elif hasattr(mw, "init_manager") and hasattr(mw.init_manager, "splitter"):
                        splitter = mw.init_manager.splitter
                        if splitter and splitter.count() > 1:
                            splitter.setSizes([0, 1000])
                except Exception as _e:
                    logging.warning("[paste_xyz.py:244] silenced: %s", _e)

        except Exception as e:
            traceback.print_exc()
            QMessageBox.critical(mw, "Error", f"Failed to parse or load data:\n{e}")

def run(mw):
    if not hasattr(mw, 'plugin_manager'):
        return

    from moleditpy.plugins.plugin_interface import PluginContext
    context = PluginContext(mw.plugin_manager, PLUGIN_NAME)
    run_plugin(context)
