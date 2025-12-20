
import traceback
import io
import contextlib
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QTextEdit, QPushButton, QHBoxLayout, 
    QMessageBox, QLabel, QLineEdit, QDialogButtonBox, QInputDialog
)
from PyQt6.QtCore import Qt
try:
    from rdkit import Chem
    from rdkit.Chem import rdGeometry, AllChem
except ImportError:
    Chem = None

PLUGIN_NAME = "Paste XYZ"
__version__="2025.12.20"
__author__="HiroYokoyama"

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

def run(main_window):
    if Chem is None:
        QMessageBox.critical(main_window, "Error", "RDKit is not available.")
        return

    dialog = PasteXYZDialog(main_window)
    if dialog.exec() == QDialog.Accepted:
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
                             # Maybe it's 'Index Symbol X Y Z' or something else?
                             # Strict requirement: "Header/Atom count ignored", assume Paste is pure or standard XYZ
                             # If line starts with a number, it's likely count or index, SKIP.
                             # Standard XYZ atom lines start with Symbol.
                             continue
                        
                        atoms_data.append((symbol, x, y, z))
                    except ValueError:
                        # Not coordinates, skip (likely header or title)
                        continue

            if not atoms_data:
                QMessageBox.warning(main_window, "Paste XYZ", "No valid coordinate lines found.\nExpected format: Symbol X Y Z")
                return

            # Clean workspace
            main_window.clear_all()
            if hasattr(main_window, 'plotter'):
                main_window.plotter.clear()
            
            # Create RWMol and add atoms/conformers
            mol = Chem.RWMol()
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                atom = Chem.Atom(symbol)
                atom.SetIntProp("xyz_unique_id", i)
                mol.AddAtom(atom)
            
            conf = Chem.Conformer(len(atoms_data))
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                conf.SetAtomPosition(i, rdGeometry.Point3D(x, y, z))
            mol.AddConformer(conf)

            # --- Chemistry Check Logic (Adapted from main_window_molecular_parsers.py) ---
            
            # Helper: prompt for charge
            def prompt_for_charge():
                try:
                    dialog = QDialog(main_window)
                    dialog.setWindowTitle("Import XYZ Charge")
                    layout = QVBoxLayout(dialog)
                    label = QLabel("Enter total molecular charge:")
                    line_edit = QLineEdit(dialog)
                    line_edit.setText("")
                    btn_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dialog)
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
                    def on_cancel():
                        dialog.reject()
                    def on_skip():
                        result["skip"] = True
                        dialog.accept()

                    try:
                        btn_box.button(QDialogButtonBox.Ok).clicked.connect(on_ok)
                        btn_box.button(QDialogButtonBox.Cancel).clicked.connect(on_cancel)
                    except Exception:
                        btn_box.accepted.connect(on_ok)
                        btn_box.rejected.connect(on_cancel)
                    skip_btn.clicked.connect(on_skip)

                    if dialog.exec() != QDialog.Accepted:
                        return None, False, False
                    if result["skip"]:
                        return 0, True, True
                    if not result["accepted"]:
                        return None, False, False
                    
                    charge_text = line_edit.text()
                except Exception:
                    # Fallback to simple input
                    try:
                        charge_text, ok = QInputDialog.getText(main_window, "Import XYZ Charge", "Enter total molecular charge:", text="0")
                        if not ok: return None, False, False
                    except Exception:
                        return 0, True, False

                try:
                    return int(str(charge_text).strip()), True, False
                except Exception:
                    try:
                        return int(float(str(charge_text).strip())), True, False
                    except Exception:
                        return 0, True, False

            # Inner helper: process with charge
            def _process_with_charge(charge_val):
                buf = io.StringIO()
                used_rd_determine = False
                mol_to_finalize = None
                with contextlib.redirect_stderr(buf):
                    try:
                        from rdkit.Chem import rdDetermineBonds
                        try:
                            # Try to copy mol
                            try:
                                mol_candidate = Chem.RWMol(Chem.Mol(mol))
                            except Exception:
                                mol_candidate = Chem.RWMol(mol)
                            
                            rdDetermineBonds.DetermineBonds(mol_candidate, charge=charge_val)
                            mol_to_finalize = mol_candidate
                            used_rd_determine = True
                        except Exception:
                            # DetermineBonds failed
                            raise RuntimeError("DetermineBondsFailed")
                    except RuntimeError:
                        raise
                    except Exception:
                        used_rd_determine = False
                        mol_to_finalize = mol
                    
                    if not used_rd_determine:
                        # Fallback to distance based
                        if hasattr(main_window, 'estimate_bonds_from_distances'):
                             main_window.estimate_bonds_from_distances(mol_to_finalize)

                    try:
                         candidate_mol = mol_to_finalize.GetMol()
                    except Exception:
                         candidate_mol = None

                    if candidate_mol is None:
                         # Salvage
                         try:
                             candidate_mol = mol.GetMol()
                         except Exception:
                             candidate_mol = None
                    
                    if candidate_mol is None:
                        raise ValueError("Failed to create valid molecule object")

                    # Attach charge prop
                    try:
                        candidate_mol.SetIntProp("_xyz_charge", int(charge_val))
                    except Exception:
                         pass
                    
                    # Apply chem check flags (if available in main_window)
                    if hasattr(main_window, '_apply_chem_check_and_set_flags'):
                        main_window._apply_chem_check_and_set_flags(candidate_mol, source_desc='PasteXYZ')

                    return candidate_mol

            # Main Logic Loop
            final_mol = None
            
            # Check settings if available (defaulting to safe behavior)
            settings = getattr(main_window, 'settings', {})
            always_ask = bool(settings.get('always_ask_charge', False))
            skip_checks_global = bool(settings.get('skip_chemistry_checks', False))

            if skip_checks_global:
                 # Skip path
                 if hasattr(main_window, 'estimate_bonds_from_distances'):
                      try: main_window.estimate_bonds_from_distances(mol)
                      except: pass
                 try:
                     final_mol = mol.GetMol()
                     final_mol.SetIntProp("_xyz_skip_checks", 1)
                     # Disable optimization for this mol
                     main_window.current_mol = final_mol
                     main_window.is_xyz_derived = True
                 except: pass
            else:
                # Normal path
                try:
                    if not always_ask:
                        try:
                            final_mol = _process_with_charge(0)
                        except RuntimeError:
                            # DetermineBonds failed for 0, loop prompt
                            while True:
                                charge_val, ok, skip_flag = prompt_for_charge()
                                if not ok: return # User cancel
                                if skip_flag:
                                    # User skipped
                                    if hasattr(main_window, 'estimate_bonds_from_distances'):
                                        try: main_window.estimate_bonds_from_distances(mol)
                                        except: pass
                                    try:
                                        final_mol = mol.GetMol()
                                        final_mol.SetIntProp("_xyz_skip_checks", 1)
                                    except: pass
                                    break
                                
                                try:
                                    final_mol = _process_with_charge(charge_val)
                                    break
                                except RuntimeError:
                                    main_window.statusBar().showMessage("DetermineBonds failed for that charge...")
                                    continue
                                except Exception:
                                    # Other error, try salvage if skip checks allowed? No, here we just continue or break
                                    continue
                    else:
                        # Always ask
                        while True:
                            charge_val, ok, skip_flag = prompt_for_charge()
                            if not ok: return
                            if skip_flag:
                                if hasattr(main_window, 'estimate_bonds_from_distances'):
                                    try: main_window.estimate_bonds_from_distances(mol)
                                    except: pass
                                try:
                                    final_mol = mol.GetMol()
                                    final_mol.SetIntProp("_xyz_skip_checks", 1)
                                except: pass
                                break
                            try:
                                final_mol = _process_with_charge(charge_val)
                                break
                            except RuntimeError:
                                main_window.statusBar().showMessage("DetermineBonds failed...")
                                continue
                            except Exception:
                                continue
                except Exception:
                     # Any other unhandled fallback
                     pass

            # If failed to get final_mol, fallback to raw
            if final_mol is None:
                if hasattr(main_window, 'estimate_bonds_from_distances'):
                    try: main_window.estimate_bonds_from_distances(mol)
                    except: pass
                try: final_mol = mol.GetMol()
                except: pass

            if final_mol:
                # We need RWMol for editing
                rw_mol = Chem.RWMol(final_mol)
                main_window.current_mol = rw_mol
                main_window.current_file_path = None
                main_window.has_unsaved_changes = True

                main_window.update_window_title()
                main_window.reset_undo_stack()
                main_window.draw_molecule_3d(rw_mol)
                main_window.fit_to_view()
                main_window.statusBar().showMessage(f"Loaded {len(atoms_data)} atoms from clipboard data.")
                
                # Enter 3D only mode as requested
                if hasattr(main_window, '_enter_3d_viewer_ui_mode'):
                    try:
                        main_window._enter_3d_viewer_ui_mode()
                    except Exception as e:
                        print(f"Could not switch to 3D mode: {e}")

        except Exception as e:
            traceback.print_exc()
            QMessageBox.critical(main_window, "Error", f"Failed to parse or load data:\n{e}")
