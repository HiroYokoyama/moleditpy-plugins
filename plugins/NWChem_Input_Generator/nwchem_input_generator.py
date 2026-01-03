# -*- coding: utf-8 -*-
import os
import json
from PyQt6.QtWidgets import (QMessageBox, QDialog, QVBoxLayout, QLabel, 
                             QLineEdit, QSpinBox, QPushButton, QGroupBox, 
                             QHBoxLayout, QComboBox, QTextEdit, QFileDialog, QFormLayout, QInputDialog, QSizePolicy, QCheckBox)
from PyQt6.QtCore import Qt
from rdkit import Chem

PLUGIN_NAME = "NWChem Input Generator"
PLUGIN_VERSION = "2026.01.03"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Generate NWChem input files for quantum chemistry calculations."
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "nwchem_input_generator.json")

class NwchemSetupDialog(QDialog):
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

        # --- General ---
        gen_group = QGroupBox("General")
        gen_layout = QFormLayout()
        
        self.title_edit = QLineEdit("NWChem Job")
        self.title_edit.textChanged.connect(self.update_preview)
        gen_layout.addRow("Title:", self.title_edit)
        
        self.start_edit = QLineEdit("[filename]")
        self.start_edit.textChanged.connect(self.update_preview)
        gen_layout.addRow("Start/Scratch Name:", self.start_edit)
        
        gen_group.setLayout(gen_layout)
        layout.addWidget(gen_group)

        # --- Method ---
        kw_group = QGroupBox("Calculation Settings")
        kw_layout = QFormLayout()
        
        self.module_combo = QComboBox()
        self.module_combo.addItems(["dft", "scf", "mp2", "ccsd", "tce"])
        self.module_combo.setCurrentText("dft")
        self.module_combo.currentIndexChanged.connect(self.update_ui_state)
        kw_layout.addRow("Theory Module:", self.module_combo)
        
        self.functional_combo = QComboBox()
        self.functional_combo.addItems(["b3lyp", "pbe0", "pbe96", "blyp", "lda", "tpss", "m06", "m06-2x"])
        self.functional_combo.setCurrentText("b3lyp")
        self.functional_combo.currentIndexChanged.connect(self.update_preview)
        kw_layout.addRow("DFT Functional:", self.functional_combo)
        
        self.task_combo = QComboBox()
        self.task_combo.addItems(["energy", "optimize", "freq", "gradient", "saddle"])
        self.task_combo.setCurrentText("optimize")
        self.task_combo.currentIndexChanged.connect(self.update_preview)
        kw_layout.addRow("Task Operation:", self.task_combo)
        
        self.basis_combo = QComboBox()
        self.basis_combo.addItems(["6-31G*", "6-311G**", "cc-pvdz", "cc-pvtz", "def2-svp", "def2-tzvp"])
        self.basis_combo.setCurrentText("6-31G*")
        self.basis_combo.currentIndexChanged.connect(self.update_preview)
        kw_layout.addRow("Basis Set:", self.basis_combo)
        
        self.ecp_warning = QLabel("")
        self.ecp_warning.setStyleSheet("color: red; font-weight: bold;")
        kw_layout.addRow("", self.ecp_warning)

        kw_group.setLayout(kw_layout)
        layout.addWidget(kw_group)

        # --- Molecular State ---
        mol_group = QGroupBox("Molecular State")
        mol_layout = QHBoxLayout()
        
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.valueChanged.connect(self.update_preview)
        
        self.mult_spin = QSpinBox()
        self.mult_spin.setRange(1, 10)
        self.mult_spin.valueChanged.connect(self.update_preview)

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

        # Signals
        self.title_edit.textChanged.connect(self.update_preview)
        self.start_edit.textChanged.connect(self.update_preview)
        self.functional_combo.currentIndexChanged.connect(self.update_preview)

        self.setLayout(layout)
        self.update_ui_state()
        self.check_heavy_atoms()
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

    def check_heavy_atoms(self):
        # ECP warning if Z > 36 (Kr) or typical transition metals if basis suggests problem
        # Simple check: Any atom > Ar (18) might need consideration, but definitly > Zn(30) or Kr(36)
        if not self.mol: return
        has_heavy = False
        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() > 36: # heavier than Kr
                has_heavy = True
                break
        
        if has_heavy:
            self.ecp_warning.setText("Warning: Heavy atoms detected. Consider using ECPs (e.g. LANL2DZ) manually.")
        else:
            self.ecp_warning.setText("")

    def update_ui_state(self):
        mod = self.module_combo.currentText()
        self.functional_combo.setEnabled(mod == "dft")
        self.update_preview()

    def generate_content(self, filename_hint=None):
        lines = []
        
        lines.append(f'title "{self.title_edit.text()}"')
        
        start_name = filename_hint if filename_hint else self.start_edit.text()
        lines.append(f'start {start_name}')
        lines.append("")
        
        # Geometry
        lines.append("geometry units angstroms noautosym")
        if self.mol:
            conf = self.mol.GetConformer()
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                atom = self.mol.GetAtomWithIdx(i)
                lines.append(f"  {atom.GetSymbol(): <3} {pos.x: >10.6f} {pos.y: >10.6f} {pos.z: >10.6f}")
        lines.append("end")
        lines.append("")
        
        # Basis
        lines.append("basis")
        lines.append(f"  * library \"{self.basis_combo.currentText()}\"")
        lines.append("end")
        lines.append("")
        
        # Charge
        charge = self.charge_spin.value()
        if charge != 0:
            lines.append(f"charge {charge}")
            lines.append("")
            
        # Module block
        module = self.module_combo.currentText()
        mult = self.mult_spin.value()
        
        lines.append(module)
        if module == "dft":
            lines.append(f"  xc {self.functional_combo.currentText()}")
            lines.append(f"  mult {mult}")
        elif module == "scf":
            if mult == 1:
                lines.append("  rhf")
                lines.append("  singlet")
            elif mult == 2:
                lines.append("  uhf")
                lines.append("  doublet")
            elif mult == 3:
                lines.append("  uhf")
                lines.append("  triplet")
            else:
                lines.append("  uhf")
                # High multiplicity, see previous turn comments
                kw_map = {4: "quartet", 5: "quintet", 6: "sextet", 7: "septet", 8: "octet"}
                if mult in kw_map:
                    lines.append(f"  {kw_map[mult]}")
                else:
                    lines.append(f"  # Generic multiplicity {mult}")
        
        lines.append("end")
        lines.append("")
        
        # Task
        task_op = self.task_combo.currentText()
        lines.append(f"task {module} {task_op}")
        
        return "\n".join(lines)

    def update_preview(self):
        # Use current start_edit value for preview
        self.preview_text.setText(self.generate_content())

    def save_file(self):
        default_name = "nwchem_job.nw"
        if self.filename:
            base = os.path.splitext(os.path.basename(self.filename))[0]
            default_name = f"{base}.nw"
            
        path, _ = QFileDialog.getSaveFileName(self, "Save NWChem Input", default_name, "NWChem Input (*.nw)")
        if path:
            try:
                # Update 'start' directive
                base_name = os.path.splitext(os.path.basename(path))[0]
                start_line = f"start {base_name}"
                
                content = self.preview_text.toPlainText()
                
                import re
                if re.search(r"^start\s+.*$", content, re.MULTILINE):
                    content = re.sub(r"^start\s+.*$", start_line, content, flags=re.MULTILINE)
                else:
                    content = f"{start_line}\n{content}"
                
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
        
        # Block signals
        self.module_combo.blockSignals(True)
        self.functional_combo.blockSignals(True)
        self.task_combo.blockSignals(True)
        self.basis_combo.blockSignals(True)
        self.title_edit.blockSignals(True)
        self.start_edit.blockSignals(True)
        
        self.module_combo.setCurrentText(data.get("module", "dft"))
        self.update_ui_state() # Enable/disable
        self.functional_combo.setCurrentText(data.get("functional", "b3lyp"))
        self.task_combo.setCurrentText(data.get("task", "optimize"))
        self.basis_combo.setCurrentText(data.get("basis", "6-31G*"))
        self.title_edit.setText(data.get("title", ""))
        self.start_edit.setText(data.get("start", ""))
        
        self.module_combo.blockSignals(False)
        self.functional_combo.blockSignals(False)
        self.task_combo.blockSignals(False)
        self.basis_combo.blockSignals(False)
        
        self.update_preview()
        self.title_edit.blockSignals(False)
        self.start_edit.blockSignals(False)
        
        self.update_preview()

    def save_preset_dialog(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:")
        if ok and name:
            self.presets_data[name] = {
                "module": self.module_combo.currentText(),
                "functional": self.functional_combo.currentText(),
                "task": self.task_combo.currentText(),
                "basis": self.basis_combo.currentText(),
                "title": self.title_edit.text(),
                "start": self.start_edit.text()
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
    dialog = NwchemSetupDialog(parent=mw, mol=mol, filename=filename)
    dialog.exec()
