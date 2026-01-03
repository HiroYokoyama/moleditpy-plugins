# -*- coding: utf-8 -*-
import os
import json
from PyQt6.QtWidgets import (QMessageBox, QDialog, QVBoxLayout, QLabel, 
                             QLineEdit, QSpinBox, QPushButton, QGroupBox, 
                             QHBoxLayout, QComboBox, QTextEdit, QFileDialog, QFormLayout, QInputDialog, QSizePolicy, QCheckBox)
from PyQt6.QtCore import Qt
from rdkit import Chem

PLUGIN_NAME = "GAMESS Input Generator"
PLUGIN_VERSION = "2026.01.03"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Generate GAMESS input files for ab initio quantum chemistry."
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "gamess_input_generator.json")

class GamessSetupDialog(QDialog):
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

        # --- Control Options ---
        ctrl_group = QGroupBox("Control Options ($CONTRL)")
        ctrl_layout = QFormLayout()
        
        self.run_type = QComboBox()
        self.run_type.addItems(["ENERGY", "OPTIMIZE", "HESSIAN", "SADPOINT", "GRADIENT"])
        self.run_type.setCurrentText("OPTIMIZE")
        self.run_type.currentIndexChanged.connect(self.update_preview)
        ctrl_layout.addRow("Run Type (RUNTYP):", self.run_type)
        
        self.scf_type = QComboBox()
        self.scf_type.addItems(["RHF", "UHF", "ROHF", "GVB", "MCSCF"])
        self.scf_type.currentIndexChanged.connect(self.update_preview)
        ctrl_layout.addRow("SCF Type (SCFTYP):", self.scf_type)
        
        self.dft_type = QComboBox()
        self.dft_type.addItems(["None (Hartree-Fock)", "B3LYP", "PBE0", "PBE", "BLYP", "SVWN", "M06", "M06-2X"])
        self.dft_type.currentIndexChanged.connect(self.update_preview)
        ctrl_layout.addRow("DFT Functional (DFTTYP):", self.dft_type)
        
        self.chk_nosym = QCheckBox("Disable Symmetry (NOSYM=1)")
        self.chk_nosym.setToolTip("Important to prevent 'symmetry unique' errors for arbitrary coordinates.")
        self.chk_nosym.setChecked(True) # Strongly recommended
        self.chk_nosym.stateChanged.connect(self.update_preview)
        ctrl_layout.addRow(self.chk_nosym)

        ctrl_group.setLayout(ctrl_layout)
        layout.addWidget(ctrl_group)

        # --- Basis Set ---
        basis_group = QGroupBox("Basis Set ($BASIS)")
        basis_layout = QFormLayout()
        
        self.basis_gbasis = QComboBox()
        self.basis_gbasis.addItems(["N31", "N21", "STO", "TZV", "MINI", "MIDI"])
        self.basis_gbasis.setCurrentText("N31")
        self.basis_gbasis.currentIndexChanged.connect(self.update_basis_options)
        basis_layout.addRow("Group Basis (GBASIS):", self.basis_gbasis)
        
        self.basis_ngauss = QSpinBox()
        self.basis_ngauss.setValue(6)
        self.basis_ngauss.valueChanged.connect(self.update_preview)
        basis_layout.addRow("Num Gaussians (NGAUSS):", self.basis_ngauss)
        
        self.basis_ndfunc = QSpinBox()
        self.basis_ndfunc.setValue(1)
        self.basis_ndfunc.valueChanged.connect(self.update_preview)
        basis_layout.addRow("Polarization d-func (NDFUNC):", self.basis_ndfunc)
        
        self.basis_npfunc = QSpinBox()
        self.basis_npfunc.setValue(0)
        self.basis_npfunc.valueChanged.connect(self.update_preview)
        basis_layout.addRow("Polarization p-func (NPFUNC):", self.basis_npfunc)
        
        basis_group.setLayout(basis_layout)
        layout.addWidget(basis_group)

        # --- System ---
        sys_group = QGroupBox("System ($SYSTEM)")
        sys_layout = QHBoxLayout()
        sys_layout.addWidget(QLabel("Memory (MWORDS):"))
        self.mem_spin = QSpinBox()
        self.mem_spin.setRange(1, 9999)
        self.mem_spin.setValue(100)
        self.mem_spin.valueChanged.connect(self.update_preview)
        sys_layout.addWidget(self.mem_spin)
        sys_group.setLayout(sys_layout)
        layout.addWidget(sys_group)

        # --- Molecular State ---
        mol_group = QGroupBox("Molecular State")
        mol_layout = QHBoxLayout()
        
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.valueChanged.connect(self.update_preview)
        
        self.mult_spin = QSpinBox()
        self.mult_spin.setRange(1, 10)
        self.mult_spin.valueChanged.connect(self.update_preview)

        mol_layout.addWidget(QLabel("Charge (ICHARG):"))
        mol_layout.addWidget(self.charge_spin)
        mol_layout.addWidget(QLabel("Multiplicity (MULT):"))
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
        
        self.setLayout(layout)
        self.update_basis_options() # Initial Set
        self.update_preview()

    def update_basis_options(self):
        gbasis = self.basis_gbasis.currentText()
        if gbasis == "STO":
            self.basis_ngauss.setValue(3)
            self.basis_ngauss.setEnabled(True)
        elif gbasis in ["N21", "N31"]:
            self.basis_ngauss.setValue(6)
            self.basis_ngauss.setEnabled(True)
        elif gbasis == "TZV":
            self.basis_ngauss.setValue(1) 
            pass
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
            
            # Suggest SCFTYP based on mult
            if mult == 1:
                self.scf_type.setCurrentText("RHF")
            else:
                self.scf_type.setCurrentText("ROHF") # Safer default for open shell
        except:
            pass

    def generate_content(self):
        lines = []
        
        # $CONTRL
        runtyp = self.run_type.currentText()
        scftyp = self.scf_type.currentText()
        dfttyp = self.dft_type.currentText()
        icharg = self.charge_spin.value()
        mult = self.mult_spin.value()
        
        contrl = f" $CONTRL SCFTYP={scftyp} RUNTYP={runtyp} ICHARG={icharg} MULT={mult} COORD=CART"
        if "None" not in dfttyp:
            contrl += f" DFTTYP={dfttyp}"
        
        if self.chk_nosym.isChecked():
            contrl += " NOSYM=1"
        
        contrl += " $END"
        lines.append(contrl)
        
        # $SYSTEM
        mwords = self.mem_spin.value()
        lines.append(f" $SYSTEM MWORDS={mwords} $END")
        
        # $BASIS
        gbasis = self.basis_gbasis.currentText()
        ngauss = self.basis_ngauss.value()
        ndfunc = self.basis_ndfunc.value()
        npfunc = self.basis_npfunc.value()
        
        basis = f" $BASIS GBASIS={gbasis} NGAUSS={ngauss}"
        if ndfunc > 0: basis += f" NDFUNC={ndfunc}"
        if npfunc > 0: basis += f" NPFUNC={npfunc}"
        basis += " $END"
        lines.append(basis)
        
        # $DATA
        lines.append(" $DATA")
        lines.append("GAMESS Job generated by MoleditPy")
        lines.append("C1") # Point Group
        
        if self.mol:
            conf = self.mol.GetConformer()
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                atom = self.mol.GetAtomWithIdx(i)
                an = atom.GetAtomicNum()
                symbol = atom.GetSymbol()
                # GAMESS Data format: Name Z X Y Z
                # Name can be Symbol+Number or just Symbol
                lines.append(f"{symbol: <4} {an: >3} {pos.x: >12.8f} {pos.y: >12.8f} {pos.z: >12.8f}")
                
        lines.append(" $END")
        
        return "\n".join(lines)

    def update_preview(self):
        self.preview_text.setText(self.generate_content())

    def save_file(self):
        # Use content from preview text to allow manual edits
        content = self.preview_text.toPlainText()
        default_name = "gamess_job.inp"
        if self.filename:
            base = os.path.splitext(os.path.basename(self.filename))[0]
            default_name = f"{base}.inp"
            
        path, _ = QFileDialog.getSaveFileName(self, "Save GAMESS Input", default_name, "GAMESS Input (*.inp *.src)")
        if path:
            try:
                with open(path, "w", encoding="utf-8") as f:
                    f.write(content)
                QMessageBox.information(self, "Success", f"File saved to {path}")
                # Do not close dialog automatically
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
        self.run_type.blockSignals(True)
        self.scf_type.blockSignals(True)
        self.dft_type.blockSignals(True)
        self.chk_nosym.blockSignals(True)
        self.basis_gbasis.blockSignals(True)
        self.basis_ngauss.blockSignals(True)
        self.basis_ndfunc.blockSignals(True)
        self.basis_npfunc.blockSignals(True)
        self.mem_spin.blockSignals(True)
        
        self.run_type.setCurrentText(data.get("run_type", "OPTIMIZE"))
        self.scf_type.setCurrentText(data.get("scf_type", "RHF"))
        self.dft_type.setCurrentText(data.get("dft_type", "None (Hartree-Fock)"))
        self.chk_nosym.setChecked(data.get("nosym", True))
        self.basis_gbasis.setCurrentText(data.get("gbasis", "N31"))
        self.basis_ngauss.setValue(data.get("ngauss", 6))
        self.basis_ndfunc.setValue(data.get("ndfunc", 1))
        self.basis_npfunc.setValue(data.get("npfunc", 0))
        
        self.mem_spin.setValue(data.get("memory", 100))
        
        self.run_type.blockSignals(False)
        self.scf_type.blockSignals(False)
        self.dft_type.blockSignals(False)
        self.chk_nosym.blockSignals(False)
        self.basis_gbasis.blockSignals(False)
        self.basis_ngauss.blockSignals(False)
        self.basis_ndfunc.blockSignals(False)
        self.basis_npfunc.blockSignals(False)
        self.mem_spin.blockSignals(False)
        
        self.update_preview()

    def save_preset_dialog(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:")
        if ok and name:
            self.presets_data[name] = {
                "run_type": self.run_type.currentText(),
                "scf_type": self.scf_type.currentText(),
                "dft_type": self.dft_type.currentText(),
                "nosym": self.chk_nosym.isChecked(),
                "gbasis": self.basis_gbasis.currentText(),
                "ngauss": self.basis_ngauss.value(),
                "ndfunc": self.basis_ndfunc.value(),
                "npfunc": self.basis_npfunc.value(),
                "memory": self.mem_spin.value()
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
    dialog = GamessSetupDialog(parent=mw, mol=mol, filename=filename)
    dialog.exec()
