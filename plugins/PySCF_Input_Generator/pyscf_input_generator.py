# -*- coding: utf-8 -*-
import os
import json
import logging
import re
from PySide6.QtWidgets import (
    QMessageBox,
    QDialog,
    QVBoxLayout,
    QLabel,
    QSpinBox,
    QDoubleSpinBox,
    QPushButton,
    QGroupBox,
    QHBoxLayout,
    QComboBox,
    QTextEdit,
    QFileDialog,
    QFormLayout,
    QCheckBox,
    QInputDialog,
    QSizePolicy,
)
from PySide6.QtCore import Qt
from rdkit import Chem

PLUGIN_NAME = "PySCF Input Generator"
PLUGIN_VERSION = "2026.10.19"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=5.0.0, <6.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Generates Python scripts for PySCF quantum chemistry calculations."
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "pyscf_input_generator.json")

FUNCTIONALS = [
    "b3lyp", "pbe", "bp86", "blyp", "pw91", "pbe0", "b3pw91", "hse06",
    "tpss", "scan", "r2scan", "m06-l", "m06", "m06-2x", "tpssh", "m11",
    "wb97x-v", "wb97m-v", "cam-b3lyp", "lrc-wpbe", "lrc-wpbeh", "lda", "svwn",
]
BASES = [
    "sto-3g", "3-21g", "6-31g", "6-31+g", "6-31g*", "6-31+g*", "6-31++g**",
    "6-311g", "6-311g*", "6-311+g*", "6-311++g**", "6-311++g(2d,2p)",
    "6-311++g(3df,3pd)", "cc-pvdz", "cc-pvtz", "cc-pvqz", "cc-pv5z",
    "aug-cc-pvdz", "aug-cc-pvtz", "aug-cc-pvqz", "def2-svp", "def2-tzvp",
    "def2-qzvp", "def2-tzvpp", "def2-qzvpp", "ma-def2-svp", "ma-def2-tzvp",
    "ma-def2-qzvp", "lanl2dz",
]
XC_ALIASES = {"m11": "hyb_mgga_x_m11,mgga_c_m11"}
DISPERSION = {"None": None, "D3(BJ)": "d3bj", "D3(zero)": "d3zero", "D4": "d4"}
SOLVENTS = {
    "None (Vacuum)": None, "Water": 78.2, "Ethanol": 24.5, "Methanol": 32.7,
    "Acetone": 20.7, "THF": 7.58, "Chloroform": 4.81,
    "Dichloromethane": 8.93, "Toluene": 2.38, "Benzene": 2.27,
}


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

        self.functional_combo = QComboBox()  # For DFT
        self.functional_combo.addItems(FUNCTIONALS)
        self.functional_combo.setEnabled(False)
        self.functional_combo.currentIndexChanged.connect(self.update_preview)
        method_layout.addRow("Functional (DFT):", self.functional_combo)

        self.post_hf_combo = QComboBox()  # For Post-HF
        self.post_hf_combo.addItems(["MP2", "CCSD", "CCSD(T)"])
        self.post_hf_combo.setEnabled(False)
        self.post_hf_combo.currentIndexChanged.connect(self.update_preview)
        method_layout.addRow("Post-HF Method:", self.post_hf_combo)

        self.post_hf_ref_combo = QComboBox()  # Reference for Post-HF
        self.post_hf_ref_combo.addItems(["RHF", "UHF", "ROHF"])
        self.post_hf_ref_combo.setEnabled(False)
        self.post_hf_ref_combo.currentIndexChanged.connect(self.update_preview)
        method_layout.addRow("Post-HF Reference:", self.post_hf_ref_combo)

        self.reference_combo = QComboBox()
        self.reference_combo.addItems(["RHF", "UHF", "ROHF", "RKS", "UKS", "ROKS"])
        self.reference_combo.currentIndexChanged.connect(self.update_preview)
        method_layout.addRow("HF/DFT Reference:", self.reference_combo)

        method_group.setLayout(method_layout)
        layout.addWidget(method_group)

        # --- Basis & Molecule ---
        mol_group = QGroupBox("Molecule & Basis")
        mol_layout = QFormLayout()

        self.basis_combo = QComboBox()
        self.basis_combo.addItems(BASES)
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
        self.symmetry_check.setChecked(False)
        self.symmetry_check.toggled.connect(self.update_preview)
        mol_layout.addRow(self.symmetry_check)

        mol_group.setLayout(mol_layout)
        layout.addWidget(mol_group)

        options_group = QGroupBox("Calculation Options")
        options_layout = QFormLayout()
        self.memory_spin = QSpinBox()
        self.memory_spin.setRange(100, 1000000)
        self.memory_spin.setValue(4000)
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(0, 1024)
        self.cycles_spin = QSpinBox()
        self.cycles_spin.setRange(1, 2000)
        self.cycles_spin.setValue(100)
        self.grid_spin = QSpinBox()
        self.grid_spin.setRange(0, 9)
        self.grid_spin.setValue(3)
        self.tolerance_spin = QDoubleSpinBox()
        self.tolerance_spin.setDecimals(12)
        self.tolerance_spin.setRange(1e-12, 1)
        self.tolerance_spin.setValue(1e-9)
        self.numeric_options = {
            "memory": self.memory_spin, "threads": self.threads_spin,
            "cycles": self.cycles_spin, "grid": self.grid_spin,
            "tolerance": self.tolerance_spin,
        }
        for label, widget in (
            ("Memory (MB):", self.memory_spin), ("Threads (0 = auto):", self.threads_spin),
            ("SCF cycles:", self.cycles_spin), ("DFT grid:", self.grid_spin),
            ("SCF tolerance:", self.tolerance_spin),
        ):
            widget.valueChanged.connect(self.update_preview)
            options_layout.addRow(label, widget)
        self.solvent_combo = QComboBox()
        self.solvent_combo.addItems(list(SOLVENTS))
        self.dispersion_combo = QComboBox()
        self.dispersion_combo.addItems(list(DISPERSION))
        for label, widget in (("Solvent (ddCOSMO):", self.solvent_combo),
                              ("Dispersion:", self.dispersion_combo)):
            widget.currentIndexChanged.connect(self.update_preview)
            options_layout.addRow(label, widget)
        options_group.setLayout(options_layout)
        layout.addWidget(options_group)

        # --- Preview ---
        preview_group = QGroupBox("Preview (Python Code)")
        preview_layout = QVBoxLayout()
        self.preview_text = QTextEdit()
        self.preview_text.setReadOnly(False)  # Editable
        self.preview_text.setMinimumHeight(200)  # Adjusted height
        self.preview_text.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
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
        if not self.mol:
            return
        try:
            charge = Chem.GetFormalCharge(self.mol)
            self.charge_spin.setValue(charge)

            total_electrons = 0
            unpaired = 0
            for atom in self.mol.GetAtoms():
                total_electrons += atom.GetAtomicNum()
                unpaired += atom.GetNumRadicalElectrons()
            total_electrons -= charge

            # Explicit radicals set the multiplicity; parity only settles what they
            # leave open. Deriving from parity alone wrote triplet O2 as a singlet.
            if (total_electrons - unpaired) % 2 != 0:
                unpaired = unpaired + 1 if unpaired == 0 else unpaired - 1
            mult = unpaired + 1
            self.mult_spin.setValue(mult)
        except Exception as _e:
            logging.warning("[pyscf_input_generator.py:157] silenced: %s", _e)

    def update_method_options(self):
        cat = self.category_combo.currentText()
        self.functional_combo.setEnabled(cat == "DFT")
        self.post_hf_combo.setEnabled("Post-HF" in cat)
        self.post_hf_ref_combo.setEnabled("Post-HF" in cat)
        self.reference_combo.setEnabled("Post-HF" not in cat)
        reference = self.reference_combo.currentText()
        if cat == "DFT" and reference.endswith("HF"):
            self.reference_combo.setCurrentText(reference.replace("HF", "KS"))
        elif cat == "Hartree-Fock" and reference.endswith("KS"):
            self.reference_combo.setCurrentText(reference.replace("KS", "HF"))
        self.update_preview()

    def generate_content(self, filename_hint="job"):
        cat = self.category_combo.currentText()
        functional = self.functional_combo.currentText()
        dispersion = DISPERSION[self.dispersion_combo.currentText()]
        solvent_eps = SOLVENTS[self.solvent_combo.currentText()]
        if cat == "DFT" and functional not in FUNCTIONALS:
            raise ValueError(f"Unsupported PySCF functional: {functional}")
        if cat == "DFT" and functional.endswith("-v") and dispersion:
            raise ValueError("VV10 functionals already include dispersion; select None.")
        if "Post-HF" in cat and dispersion:
            raise ValueError("Dispersion corrections are available for HF and DFT jobs.")
        lines = []
        lines.append("from pyscf import gto, scf, dft, mp, cc")
        lines.append("from pyscf import lib")
        if solvent_eps is not None:
            lines.append("from pyscf import solvent")
        if dispersion:
            lines.append("from pyscf import dispersion  # requires pyscf-dispersion")
        lines.append("")
        if self.threads_spin.value():
            lines.append(f"lib.num_threads({self.threads_spin.value()})")

        # Molecule
        lines.append("mol = gto.M(")
        lines.append("    atom='''")
        if self.mol:
            conf = self.mol.GetConformer()
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                atom = self.mol.GetAtomWithIdx(i)
                lines.append(
                    f"    {atom.GetSymbol()} {pos.x:.8f} {pos.y:.8f} {pos.z:.8f}"
                )
        lines.append("    ''',")

        lines.append(f"    basis='{self.basis_combo.currentText()}',")
        basis = self.basis_combo.currentText()
        ecp = {}
        if self.mol:
            for atom in self.mol.GetAtoms():
                if basis == "lanl2dz" and atom.GetAtomicNum() > 10:
                    ecp[atom.GetSymbol()] = "lanl2dz"
                elif (basis.startswith("def2-") or basis.startswith("ma-def2-")) and atom.GetAtomicNum() > 36:
                    ecp[atom.GetSymbol()] = "def2-svp"
        if ecp:
            lines.append(f"    ecp={ecp!r},")
        lines.append(f"    max_memory={self.memory_spin.value()},")
        lines.append(f"    charge={self.charge_spin.value()},")

        mult = self.mult_spin.value()
        spin = mult - 1
        lines.append(f"    spin={spin},  # 2S = Multiplicity - 1")

        sym = "True" if self.symmetry_check.isChecked() else "False"
        lines.append(f"    symmetry={sym}")
        lines.append(")")
        lines.append("")

        # Method
        cat = self.category_combo.currentText()

        reference = (self.post_hf_ref_combo.currentText() if "Post-HF" in cat
                     else self.reference_combo.currentText())
        allowed = ("RKS", "UKS", "ROKS") if cat == "DFT" else ("RHF", "UHF", "ROHF")
        if reference not in allowed:
            raise ValueError(f"{cat} requires one of {', '.join(allowed)}.")
        if mult > 1 and reference in ("RHF", "RKS"):
            reference = "UHF" if reference == "RHF" else "UKS"
        namespace = "dft" if cat == "DFT" else "scf"
        lines.append(f"mf = {namespace}.{reference}(mol)")
        if cat == "DFT":
            lines.append(f"mf.xc = {XC_ALIASES.get(functional, functional)!r}")
            lines.append(f"mf.grids.level = {self.grid_spin.value()}")
            if self.grid_spin.value() >= 4:
                lines.append("mf.grids.prune = None")
        if solvent_eps is not None:
            lines.append("mf = solvent.ddCOSMO(mf)")
            lines.append(f"mf.with_solvent.eps = {solvent_eps}")
        if dispersion:
            lines.append(f"mf.disp = {dispersion!r}")
        lines.append(f"mf.max_cycle = {self.cycles_spin.value()}")
        lines.append(f"mf.conv_tol = {self.tolerance_spin.value()!r}")

        # Checkpoint file logic
        # Default hint is 'job', but will be updated on save
        chk_name = f"{filename_hint}.chk"
        lines.append(f"mf.chkfile = {chk_name!r}")

        lines.append("mf.kernel()")
        lines.extend(["if not mf.converged:", "    raise RuntimeError('SCF did not converge')"])
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
                lines.append("mycc = cc.CCSD(mf)")
                lines.append("mycc.kernel()")
                # mycc.ccsd_t() dispatches to the restricted or unrestricted
                # (T) code by itself; the manual pyscf.cc.ccsd_t.kernel(...)
                # call is RHF-only and crashes for UHF/ROHF references.
                lines.append("e_t = mycc.ccsd_t()")
                lines.append("print('CCSD(T) correction =', e_t)")
                lines.append("print('Total CCSD(T) energy =', mycc.e_tot + e_t)")

        return "\n".join(lines)

    def update_preview(self):
        self._generation_error = None
        try:
            self.preview_text.setText(self.generate_content(filename_hint="[filename]"))
        except ValueError as exc:
            self._generation_error = str(exc)
            self.preview_text.setText(f"# Invalid input: {exc}")

    def save_file(self):
        if self._generation_error:
            QMessageBox.warning(self, "Invalid PySCF input", self._generation_error)
            return
        default_name = "pyscf_job.py"
        if self.filename:
            base = os.path.splitext(os.path.basename(self.filename))[0]
            default_name = f"{base}_pyscf.py"

        path, _ = QFileDialog.getSaveFileName(
            self, "Save PySCF Input", default_name, "Python Script (*.py)"
        )
        if path:
            try:
                # Re-generate content with correct chkfile name if standard pattern is mostly kept
                # Or try to replace it in the current text if possible.
                # Since we allow editing, a full regen might overwrite user edits.
                # Heuristic: If line "mf.chkfile =" exists, replace it. If not, append it?
                # User request: "chk file name space will be inputted in saving"

                base_name = os.path.splitext(os.path.basename(path))[0]
                chk_line = f"mf.chkfile = {(base_name + '.chk')!r}"

                content = self.preview_text.toPlainText()

                chk_pattern = r"(?m)^([ \t]*)mf\.chkfile\s*=[^\n]*$"
                if re.search(chk_pattern, content):
                    content = re.sub(
                        chk_pattern, lambda match: match[1] + chk_line, content
                    )
                else:
                    # Insert before kernel
                    if "mf.kernel()" in content:
                        content = content.replace(
                            "mf.kernel()", f"{chk_line}\nmf.kernel()"
                        )
                    else:
                        content += f"\n{chk_line}\n"

                with open(path, "w", encoding="utf-8", newline="\n") as f:
                    f.write(content)

                QMessageBox.information(self, "Success", f"File saved to {path}")
                # Do not close automatically
            except OSError as e:
                QMessageBox.critical(self, "Error", str(e))

    # --- Preset Management ---
    def load_presets_from_file(self):
        self.presets_data = {}
        if os.path.exists(SETTINGS_FILE):
            try:
                with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
                    self.presets_data = json.load(f)
            except (OSError, TypeError, ValueError) as e:
                logging.warning("Error loading presets: %s", e)
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
        if name not in self.presets_data:
            return
        data = self.presets_data[name]
        if data.get("functional", "b3lyp") not in FUNCTIONALS:
            QMessageBox.warning(self, "Unsupported preset", "Select a supported PySCF functional before using this preset.")
            return
        option_widgets = list(self.numeric_options.values()) + [
            self.reference_combo, self.solvent_combo, self.dispersion_combo,
        ]
        for widget in option_widgets:
            widget.blockSignals(True)

        # Block signals to prevent redundant updates
        self.category_combo.blockSignals(True)
        self.functional_combo.blockSignals(True)
        self.post_hf_combo.blockSignals(True)
        self.post_hf_ref_combo.blockSignals(True)
        self.basis_combo.blockSignals(True)
        self.symmetry_check.blockSignals(True)

        self.category_combo.setCurrentText(data.get("category", "Hartree-Fock"))
        self.update_method_options()  # Force update enabling
        self.functional_combo.setCurrentText(data.get("functional", "b3lyp"))
        self.post_hf_combo.setCurrentText(data.get("post_hf", "MP2"))
        self.post_hf_ref_combo.setCurrentText(data.get("post_hf_ref", "RHF"))
        self.basis_combo.setCurrentText(data.get("basis", "def2-svp"))
        self.symmetry_check.setChecked(data.get("symmetry", False))
        self.reference_combo.setCurrentText(data.get("reference", "RKS" if data.get("category") == "DFT" else "RHF"))
        for key, widget in self.numeric_options.items():
            if key in data:
                widget.setValue(data[key])
        self.solvent_combo.setCurrentText(data.get("solvent", "None (Vacuum)"))
        self.dispersion_combo.setCurrentText(data.get("dispersion", "None"))
        for widget in option_widgets:
            widget.blockSignals(False)

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
                "symmetry": self.symmetry_check.isChecked(),
                "reference": self.reference_combo.currentText(),
                "solvent": self.solvent_combo.currentText(),
                "dispersion": self.dispersion_combo.currentText(),
                **{key: widget.value() for key, widget in self.numeric_options.items()},
            }
            self.save_presets_to_file()
            self.update_preset_combo()
            self.preset_combo.setCurrentText(name)

    def delete_preset(self):
        name = self.preset_combo.currentText()
        if not name:
            return
        confirm = QMessageBox.question(self, "Confirm", f"Delete preset '{name}'?")
        if confirm == QMessageBox.StandardButton.Yes:
            del self.presets_data[name]
            self.save_presets_to_file()
            self.update_preset_combo()

    def save_presets_to_file(self):
        try:
            with open(SETTINGS_FILE, "w", encoding="utf-8") as f:
                json.dump(self.presets_data, f, indent=4)
        except (OSError, TypeError, ValueError) as e:
            QMessageBox.warning(self, "Error", f"Failed to save presets: {e}")


def run(mw):
    mol = getattr(mw, "current_mol", None)
    if not mol:
        QMessageBox.warning(mw, PLUGIN_NAME, "No molecule loaded.")
        return
    filename = getattr(mw, "current_file_path", None)
    dialog = PyscfSetupDialog(parent=mw, mol=mol, filename=filename)
    dialog.exec()
