# -*- coding: utf-8 -*-
import json
import os
from PyQt6.QtWidgets import (QMessageBox, QDialog, QVBoxLayout, QLabel, 
                             QLineEdit, QSpinBox, QPushButton, QGroupBox, 
                             QHBoxLayout, QComboBox, QTextEdit, QFileDialog, QFormLayout, QInputDialog, QSizePolicy, QCheckBox)
from PyQt6.QtCore import Qt
from rdkit import Chem

PLUGIN_NAME = "Psi4 Input Generator"
PLUGIN_VERSION = "2026.01.03"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Generate Psi4 input files for quantum chemistry calculations."
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "psi4_input_generator.json")

class Psi4SetupDialog(QDialog):
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

        # --- Resources ---
        res_group = QGroupBox("Resources")
        res_layout = QHBoxLayout()
        self.mem_spin = QSpinBox()
        self.mem_spin.setRange(1, 1000)
        self.mem_spin.setValue(2)
        self.mem_spin.setSuffix(" GB")
        self.mem_spin.valueChanged.connect(self.update_preview)
        
        self.thread_spin = QSpinBox()
        self.thread_spin.setRange(1, 64)
        self.thread_spin.setValue(4)
        self.thread_spin.valueChanged.connect(self.update_preview)
        
        res_layout.addWidget(QLabel("Memory:"))
        res_layout.addWidget(self.mem_spin)
        res_layout.addWidget(QLabel("Threads:"))
        res_layout.addWidget(self.thread_spin)
        res_group.setLayout(res_layout)
        layout.addWidget(res_group)

        # --- Method ---
        kw_group = QGroupBox("Calculation Settings")
        kw_layout = QFormLayout()
        
        self.task_combo = QComboBox()
        self.task_combo.addItems(["energy", "optimize", "frequency", "gradient"])
        self.task_combo.currentIndexChanged.connect(self.update_preview)
        kw_layout.addRow("Task:", self.task_combo)
        
        self.method_combo = QComboBox()
        self.method_combo.addItems(["scf", "b3lyp", "pbe", "mp2", "ccsd(t)", "wB97X-D", "m06-2x"])
        self.method_combo.setCurrentText("b3lyp")
        self.method_combo.currentIndexChanged.connect(self.update_preview)
        kw_layout.addRow("Method:", self.method_combo)
        
        self.basis_combo = QComboBox()
        self.basis_combo.addItems(["cc-pvdz", "cc-pvtz", "aug-cc-pvdz", "def2-svp", "def2-tzvp", "6-31G*", "6-311G**"])
        self.basis_combo.setCurrentText("def2-svp")
        self.basis_combo.currentIndexChanged.connect(self.update_preview)
        kw_layout.addRow("Basis Set:", self.basis_combo)
        
        self.ref_combo = QComboBox()
        self.ref_combo.addItems(["rhf", "uhf", "rohf", "rks", "uks"])
        self.ref_combo.setCurrentText("rks")
        self.ref_combo.currentIndexChanged.connect(self.update_preview)
        kw_layout.addRow("Reference:", self.ref_combo)

        kw_group.setLayout(kw_layout)
        layout.addWidget(kw_group)

        # --- Molecular State ---
        mol_group = QGroupBox("Molecular State")
        mol_layout = QHBoxLayout()
        
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.valueChanged.connect(self.update_auto_reference)
        
        self.mult_spin = QSpinBox()
        self.mult_spin.setRange(1, 10)
        self.mult_spin.valueChanged.connect(self.update_auto_reference)

        mol_layout.addWidget(QLabel("Charge:"))
        mol_layout.addWidget(self.charge_spin)
        mol_layout.addWidget(QLabel("Multiplicity:"))
        mol_layout.addWidget(self.mult_spin)
        
        mol_group.setLayout(mol_layout)
        layout.addWidget(mol_group)

        # --- Preview ---
        preview_group = QGroupBox("Preview")
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
        self.btn_save = QPushButton("Save Input File...")
        self.btn_save.clicked.connect(self.save_file)
        self.btn_save.setStyleSheet("font-weight: bold; padding: 5px;")
        btn_layout.addWidget(self.btn_save)
        layout.addLayout(btn_layout)
        
        # Connect signals
        self.method_combo.currentIndexChanged.connect(self.update_auto_reference)

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
            
            # Auto update reference triggered by signals or explicit call
            self.update_auto_reference()
        except:
            pass

    def update_auto_reference(self):
        method = self.method_combo.currentText().lower()
        mult = self.mult_spin.value()
        
        # Determine if DFT or WF
        is_dft = method not in ["scf", "mp2", "ccsd", "ccsd(t)"]
        
        # Block signals to prevent update_preview from being called multiple times
        self.ref_combo.blockSignals(True)
        if is_dft:
            # RKS or UKS
            if mult == 1:
                self.ref_combo.setCurrentText("rks")
            else:
                self.ref_combo.setCurrentText("uks")
        else:
            # RHF, UHF, ROHF
            if mult == 1:
                self.ref_combo.setCurrentText("rhf")
            else:
                self.ref_combo.setCurrentText("uhf") # Default to UHF for open shell WF
        self.ref_combo.blockSignals(False)
                
        self.update_preview()

    def generate_content(self, filename_hint="job"):
        lines = []
        
        # Memory & Threads
        lines.append(f"memory {self.mem_spin.value()} GB")
        lines.append(f"set_num_threads({self.thread_spin.value()})")
        
        # Output file logic (useful if running python style or to set scratch)
        # Psi4 input style (psi4 -i input.dat) usually sets output name from shell. 
        # But setting output_file inside allows explicit control.
        lines.append(f"psi4.core.set_output_file('{filename_hint}.out', False)")
        lines.append("")
        
        # Molecule
        lines.append("molecule {")
        lines.append(f"{self.charge_spin.value()} {self.mult_spin.value()}") # Charge Mult
        
        if self.mol:
            conf = self.mol.GetConformer()
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                atom = self.mol.GetAtomWithIdx(i)
                lines.append(f"  {atom.GetSymbol(): <3} {pos.x: >10.6f} {pos.y: >10.6f} {pos.z: >10.6f}")
        
        lines.append("}")
        lines.append("")
        
        # Set
        lines.append("set {")
        lines.append(f"  basis {self.basis_combo.currentText()}")
        lines.append(f"  reference {self.ref_combo.currentText()}")
        lines.append("}")
        lines.append("")
        
        # Task
        task = self.task_combo.currentText()
        method = self.method_combo.currentText()
        lines.append(f"{task}('{method}')")
        
        return "\n".join(lines)

    def update_preview(self):
        self.preview_text.setText(self.generate_content(filename_hint="[filename]"))

    def save_file(self):
        default_name = "psi4_job.dat"
        if self.filename:
            base = os.path.splitext(os.path.basename(self.filename))[0]
            default_name = f"{base}.dat"
            
        path, _ = QFileDialog.getSaveFileName(self, "Save PSI4 Input", default_name, "PSI4 Input (*.dat *.in)")
        if path:
            try:
                # Update output filename in content if present
                base_name = os.path.splitext(os.path.basename(path))[0]
                out_line = f"psi4.core.set_output_file('{base_name}.out', False)"
                
                content = self.preview_text.toPlainText()
                
                import re
                if "psi4.core.set_output_file" in content:
                    content = re.sub(r"psi4\.core\.set_output_file\s*\(.*?\)", out_line, content)
                else:
                    # Prepend if missing? Or just insert early
                    content = f"{out_line}\n\n{content}"
                
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
        if not name or name not in self.presets_data: return
        data = self.presets_data[name]
        
        # Block signals for dependent combos
        self.method_combo.blockSignals(True)
        self.ref_combo.blockSignals(True)
        self.basis_combo.blockSignals(True)
        self.task_combo.blockSignals(True)
        self.mem_spin.blockSignals(True)
        self.thread_spin.blockSignals(True)
        
        self.method_combo.setCurrentText(data.get("method", "b3lyp"))
        self.ref_combo.setCurrentText(data.get("reference", "rks"))
        self.basis_combo.setCurrentText(data.get("basis", "def2-svp"))
        self.task_combo.setCurrentText(data.get("task", "energy"))
        self.mem_spin.setValue(data.get("memory", 2))
        self.thread_spin.setValue(data.get("threads", 4))
        
        self.method_combo.blockSignals(False)
        self.ref_combo.blockSignals(False)
        self.basis_combo.blockSignals(False)
        self.task_combo.blockSignals(False)
        self.mem_spin.blockSignals(False)
        self.thread_spin.blockSignals(False)
        
        self.update_preview()

    def save_preset_dialog(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:")
        if ok and name:
            self.presets_data[name] = {
                "method": self.method_combo.currentText(),
                "reference": self.ref_combo.currentText(),
                "basis": self.basis_combo.currentText(),
                "task": self.task_combo.currentText(),
                "memory": self.mem_spin.value(),
                "threads": self.thread_spin.value()
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
    dialog = Psi4SetupDialog(parent=mw, mol=mol, filename=filename)
    dialog.exec()
