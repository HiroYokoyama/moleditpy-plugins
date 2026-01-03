# -*- coding: utf-8 -*-
import os
import json
from PyQt6.QtWidgets import (QMessageBox, QDialog, QVBoxLayout, QLabel, 
                             QLineEdit, QSpinBox, QPushButton, QGroupBox, 
                             QHBoxLayout, QComboBox, QTextEdit, QFileDialog, QFormLayout, QCheckBox, QInputDialog, QSizePolicy)
from PyQt6.QtCore import Qt
from rdkit import Chem

PLUGIN_NAME = "PySCF Input Generator"
PLUGIN_VERSION = "2026.01.03"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Generate Python scripts for PySCF quantum chemistry calculations."
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "pyscf_input_generator.json")

class PyscfSetupDialog(QDialog):
    def __init__(self, parent=None, mol=None, filename=None):
        super().__init__(parent)
        self.setWindowTitle(PLUGIN_NAME)
        self.resize(600, 700)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.WindowMinMaxButtonsHint)
        self.setSizeGripEnabled(True)
        self.mol = mol
        self.filename = filename
        self.setup_ui()
        self.load_presets_from_file()
        self.calc_initial_charge_mult()

    def setup_ui(self):
        layout = QVBoxLayout()

        # --- Preset Management ---
        preset_group = QGroupBox("User Presets")
        preset_layout = QHBoxLayout()
        
        self.preset_combo = QComboBox()
        self.preset_combo.currentIndexChanged.connect(self.apply_selected_preset)
        preset_layout.addWidget(QLabel("Preset:"))
        preset_layout.addWidget(self.preset_combo, 1)

        self.btn_save_preset = QPushButton("Save New...")
        self.btn_save_preset.clicked.connect(self.save_preset_dialog)
        preset_layout.addWidget(self.btn_save_preset)

        self.btn_del_preset = QPushButton("Delete")
        self.btn_del_preset.clicked.connect(self.delete_preset)
        preset_layout.addWidget(self.btn_del_preset)

        preset_group.setLayout(preset_layout)
        layout.addWidget(preset_group)

        # --- Method ---
        method_group = QGroupBox("Calculation Method")
        method_layout = QFormLayout()
        
        self.category_combo = QComboBox()
        self.category_combo.addItems(["Hartree-Fock", "DFT", "Post-HF (MP2, CCSD)"])
        self.category_combo.currentIndexChanged.connect(self.update_method_options)
        method_layout.addRow("Category:", self.category_combo)
        
        self.functional_combo = QComboBox() # For DFT
        self.functional_combo.addItems(["b3lyp", "pbe", "blyp", "svwn", "m06", "wb97x-d"])
        self.functional_combo.setEnabled(False)
        self.functional_combo.currentIndexChanged.connect(self.update_preview)
        method_layout.addRow("Functional (DFT):", self.functional_combo)
        
        self.post_hf_combo = QComboBox() # For Post-HF
        self.post_hf_combo.addItems(["MP2", "CCSD", "CCSD(T)"])
        self.post_hf_combo.setEnabled(False)
        self.post_hf_combo.currentIndexChanged.connect(self.update_preview)
        method_layout.addRow("Post-HF Method:", self.post_hf_combo)
        
        self.post_hf_ref_combo = QComboBox() # Reference for Post-HF
        self.post_hf_ref_combo.addItems(["RHF", "UHF", "ROHF"])
        self.post_hf_ref_combo.setEnabled(False)
        self.post_hf_ref_combo.currentIndexChanged.connect(self.update_preview)
        method_layout.addRow("Post-HF Reference:", self.post_hf_ref_combo)
        
        method_group.setLayout(method_layout)
        layout.addWidget(method_group)

        # --- Basis & Molecule ---
        mol_group = QGroupBox("Molecule & Basis")
        mol_layout = QFormLayout()
        
        self.basis_combo = QComboBox()
        self.basis_combo.addItems(["sto-3g", "6-31g", "6-31g*", "cc-pvdz", "cc-pvtz", "def2-svp", "def2-tzvp"])
        self.basis_combo.setCurrentText("def2-svp")
        self.basis_combo.currentIndexChanged.connect(self.update_preview)
        mol_layout.addRow("Basis Set:", self.basis_combo)
        
        # Charge/Spin
        cs_layout = QHBoxLayout()
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.valueChanged.connect(self.update_preview)
        
        self.mult_spin = QSpinBox()
        self.mult_spin.setRange(1, 10)
        self.mult_spin.valueChanged.connect(self.update_preview)
        
        cs_layout.addWidget(QLabel("Charge:"))
        cs_layout.addWidget(self.charge_spin)
        cs_layout.addWidget(QLabel("Multiplicity (2S+1):"))
        cs_layout.addWidget(self.mult_spin)
        
        mol_layout.addRow(cs_layout)
        
        self.symmetry_check = QCheckBox("Use Symmetry")
        self.symmetry_check.setChecked(True)
        self.symmetry_check.toggled.connect(self.update_preview)
        mol_layout.addRow(self.symmetry_check)
        
        mol_group.setLayout(mol_layout)
        layout.addWidget(mol_group)

        # --- Preview ---
        preview_group = QGroupBox("Preview (Python Code)")
        preview_layout = QVBoxLayout()
        self.preview_text = QTextEdit()
        self.preview_text.setReadOnly(False) # Editable
        self.preview_text.setMinimumHeight(200) # Adjusted height
        self.preview_text.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        preview_layout.addWidget(self.preview_text, 1)
        
        btn_refresh = QPushButton("Reset/Refresh Preview")
        btn_refresh.clicked.connect(self.update_preview)
        preview_layout.addWidget(btn_refresh)
        
        preview_group.setLayout(preview_layout)
        layout.addWidget(preview_group, 1)

        # --- Buttons ---
        btn_layout = QHBoxLayout()
        self.btn_save = QPushButton("Save Python Script...")
        self.btn_save.clicked.connect(self.save_file)
        self.btn_save.setStyleSheet("font-weight: bold; padding: 5px;")
        btn_layout.addWidget(self.btn_save)
        layout.addLayout(btn_layout)
        
        self.setLayout(layout)
        self.update_preview()

    def calc_initial_charge_mult(self):
        if not self.mol: return
        try:
            charge = Chem.GetFormalCharge(self.mol)
            self.charge_spin.setValue(charge)
            
            total_electrons = 0
            for atom in self.mol.GetAtoms():
                total_electrons += atom.GetAtomicNum()
            total_electrons -= charge
            
            mult = 1 if total_electrons % 2 == 0 else 2
            self.mult_spin.setValue(mult)
        except:
            pass

    def update_method_options(self):
        cat = self.category_combo.currentText()
        self.functional_combo.setEnabled(cat == "DFT")
        self.post_hf_combo.setEnabled("Post-HF" in cat)
        self.post_hf_ref_combo.setEnabled("Post-HF" in cat)
        self.update_preview()

    def generate_content(self, filename_hint="job"):
        lines = []
        lines.append("from pyscf import gto, scf, dft, mp, cc")
        lines.append("")
        
        # Molecule
        lines.append("mol = gto.M(")
        lines.append("    atom='''")
        if self.mol:
            conf = self.mol.GetConformer()
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                atom = self.mol.GetAtomWithIdx(i)
                lines.append(f"    {atom.GetSymbol()} {pos.x:.8f} {pos.y:.8f} {pos.z:.8f}")
        lines.append("    ''',")
        
        lines.append(f"    basis='{self.basis_combo.currentText()}',")
        lines.append(f"    charge={self.charge_spin.value()},")
        
        mult = self.mult_spin.value()
        spin = mult - 1
        lines.append(f"    spin={spin},  # 2S = Multiplicity - 1")
        
        sym = "True" if self.symmetry_check.isChecked() else "False"
        lines.append(f"    symmetry={sym}")
        lines.append(")")
        lines.append("mol.build()")
        lines.append("")
        
        # Method
        cat = self.category_combo.currentText()
        
        if cat == "Hartree-Fock":
            if mult == 1:
                lines.append("mf = scf.RHF(mol)")
            else:
                lines.append("mf = scf.ROHF(mol)") 
        elif cat == "DFT":
            func = self.functional_combo.currentText()
            if mult == 1:
                lines.append("mf = dft.RKS(mol)")
            else:
                lines.append("mf = dft.ROKS(mol)")
            lines.append(f"mf.xc = '{func}'")
        elif "Post-HF" in cat:
            ref_type = self.post_hf_ref_combo.currentText()
            if ref_type == "RHF": lines.append("mf = scf.RHF(mol)")
            elif ref_type == "UHF": lines.append("mf = scf.UHF(mol)")
            elif ref_type == "ROHF": lines.append("mf = scf.ROHF(mol)")
        
        # Checkpoint file logic
        # Default hint is 'job', but will be updated on save
        chk_name = f"{filename_hint}.chk"
        lines.append(f"mf.chkfile = '{chk_name}'")
        
        lines.append("mf.kernel()")
        lines.append("")

        if "Post-HF" in cat:
            phf = self.post_hf_combo.currentText()
            if phf == "MP2":
                lines.append("pt = mp.MP2(mf)")
                lines.append("pt.kernel()")
            elif phf == "CCSD":
                lines.append("mycc = cc.CCSD(mf)")
                lines.append("mycc.kernel()")
            elif phf == "CCSD(T)":
                lines.append("from pyscf.cc import ccsd_t")
                lines.append("mycc = cc.CCSD(mf)")
                lines.append("mycc.kernel()")
                lines.append("e_t = ccsd_t.kernel(mycc, mycc.ao2mo())")
                lines.append("print('CCSD(T) correction =', e_t)")
                lines.append("print('Total CCSD(T) energy =', mycc.e_tot + e_t)")

        return "\n".join(lines)

    def update_preview(self):
        # Use a placeholder when just previewing
        self.preview_text.setText(self.generate_content(filename_hint="[filename]"))

    def save_file(self):
        default_name = "pyscf_job.py"
        if self.filename:
            base = os.path.splitext(os.path.basename(self.filename))[0]
            default_name = f"{base}_pyscf.py"
            
        path, _ = QFileDialog.getSaveFileName(self, "Save PySCF Input", default_name, "Python Script (*.py)")
        if path:
            try:
                # Re-generate content with correct chkfile name if standard pattern is mostly kept
                # Or try to replace it in the current text if possible. 
                # Since we allow editing, a full regen might overwrite user edits.
                # Heuristic: If line "mf.chkfile =" exists, replace it. If not, append it?
                # User request: "chk file name space will be inputted in saving"
                
                base_name = os.path.splitext(os.path.basename(path))[0]
                chk_line = f"mf.chkfile = '{base_name}.chk'"
                
                content = self.preview_text.toPlainText()
                
                # Simple replacement regex-like
                import re
                if "mf.chkfile =" in content:
                    content = re.sub(r"mf\.chkfile\s*=\s*['\"].*['\"]", chk_line, content)
                else:
                    # Insert before kernel
                    if "mf.kernel()" in content:
                        content = content.replace("mf.kernel()", f"{chk_line}\nmf.kernel()")
                    else:
                        content += f"\n{chk_line}\n"
                
                with open(path, "w", encoding="utf-8") as f:
                    f.write(content)
                    
                QMessageBox.information(self, "Success", f"File saved to {path}")
                # Do not close automatically
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    # --- Preset Management ---
    def load_presets_from_file(self):
        self.presets_data = {}
        if os.path.exists(SETTINGS_FILE):
            try:
                with open(SETTINGS_FILE, 'r', encoding='utf-8') as f:
                    self.presets_data = json.load(f)
            except Exception as e:
                print(f"Error loading presets: {e}")
        self.update_preset_combo()

    def update_preset_combo(self):
        current = self.preset_combo.currentText()
        self.preset_combo.blockSignals(True)
        self.preset_combo.clear()
        self.preset_combo.addItems(sorted(self.presets_data.keys()))
        if current in self.presets_data:
            self.preset_combo.setCurrentText(current)
        self.preset_combo.blockSignals(False)
        self.apply_selected_preset()

    def apply_selected_preset(self):
        name = self.preset_combo.currentText()
        if name not in self.presets_data: return
        data = self.presets_data[name]
        
        # Block signals to prevent redundant updates
        self.category_combo.blockSignals(True)
        self.functional_combo.blockSignals(True)
        self.post_hf_combo.blockSignals(True)
        self.post_hf_ref_combo.blockSignals(True)
        self.basis_combo.blockSignals(True)
        self.symmetry_check.blockSignals(True)
        
        self.category_combo.setCurrentText(data.get("category", "Hartree-Fock"))
        self.update_method_options() # Force update enabling
        self.functional_combo.setCurrentText(data.get("functional", "b3lyp"))
        self.post_hf_combo.setCurrentText(data.get("post_hf", "MP2"))
        self.post_hf_ref_combo.setCurrentText(data.get("post_hf_ref", "RHF"))
        self.basis_combo.setCurrentText(data.get("basis", "def2-svp"))
        self.symmetry_check.setChecked(data.get("symmetry", True))
        
        self.category_combo.blockSignals(False)
        self.functional_combo.blockSignals(False)
        self.post_hf_combo.blockSignals(False)
        self.post_hf_ref_combo.blockSignals(False)
        self.basis_combo.blockSignals(False)
        self.symmetry_check.blockSignals(False)
        
        self.update_preview()

    def save_preset_dialog(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:")
        if ok and name:
            self.presets_data[name] = {
                "category": self.category_combo.currentText(),
                "functional": self.functional_combo.currentText(),
                "post_hf": self.post_hf_combo.currentText(),
                "post_hf_ref": self.post_hf_ref_combo.currentText(),
                "basis": self.basis_combo.currentText(),
                "symmetry": self.symmetry_check.isChecked()
            }
            self.save_presets_to_file()
            self.update_preset_combo()
            self.preset_combo.setCurrentText(name)

    def delete_preset(self):
        name = self.preset_combo.currentText()
        if not name: return
        confirm = QMessageBox.question(self, "Confirm", f"Delete preset '{name}'?")
        if confirm == QMessageBox.StandardButton.Yes:
            del self.presets_data[name]
            self.save_presets_to_file()
            self.update_preset_combo()

    def save_presets_to_file(self):
        try:
            with open(SETTINGS_FILE, 'w', encoding='utf-8') as f:
                json.dump(self.presets_data, f, indent=4)
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to save presets: {e}")

def run(mw):
    mol = getattr(mw, 'current_mol', None)
    if not mol:
        QMessageBox.warning(mw, PLUGIN_NAME, "No molecule loaded.")
        return
    filename = getattr(mw, 'current_file_path', None)
    dialog = PyscfSetupDialog(parent=mw, mol=mol, filename=filename)
    dialog.exec()
