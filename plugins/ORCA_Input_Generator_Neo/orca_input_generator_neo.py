# -*- coding: utf-8 -*-
import os
from PyQt6.QtWidgets import (QMessageBox, QDialog, QVBoxLayout, QLabel, 
                             QLineEdit, QSpinBox, QPushButton, QFileDialog, 
                             QFormLayout, QGroupBox, QHBoxLayout, QComboBox, QTextEdit, 
                             QTabWidget, QCheckBox, QWidget, QScrollArea, QMenu)
from PyQt6.QtGui import QPalette, QColor, QAction
from PyQt6.QtCore import Qt
from rdkit import Chem
import json

__version__="2025.12.18"
__author__="HiroYokoyama"
PLUGIN_NAME = "ORCA Input Generator Neo"
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "orca_input_generator_neo.json")

class OrcaKeywordBuilderDialog(QDialog):
    """
    Dialog to construct the ORCA Job Route line.
    """
    def __init__(self, parent=None, current_route=""):
        super().__init__(parent)
        self.setWindowTitle("ORCA Keyword Builder")
        self.resize(700, 600)
        self.ui_ready = False
        self.current_route = current_route
        self.setup_ui()
        self.parse_route(current_route)

    def setup_ui(self):
        layout = QVBoxLayout()
        
        self.tabs = QTabWidget()
        
        # --- Tab 1: Method & Basis ---
        self.tab_method = QWidget()
        self.setup_method_tab()
        self.tabs.addTab(self.tab_method, "Method/Basis")
        
        # --- Tab 2: Job Type ---
        self.tab_job = QWidget()
        self.setup_job_tab()
        self.tabs.addTab(self.tab_job, "Job Type")
        
        # --- Tab 3: Solvation & Disp ---
        self.tab_solvation = QWidget()
        self.setup_solvation_tab()
        self.tabs.addTab(self.tab_solvation, "Solvation/Dispersion")
        
        # --- Tab 4: Properties ---
        self.tab_props = QWidget()
        self.setup_props_tab()
        self.tabs.addTab(self.tab_props, "Properties")

        layout.addWidget(self.tabs)

        # --- Preview ---
        preview_group = QGroupBox("Keyword Preview")
        preview_layout = QVBoxLayout()
        self.preview_label = QLabel()
        self.preview_label.setWordWrap(True)
        self.preview_label.setStyleSheet("font-weight: bold; color: blue; font-size: 14px;")
        preview_layout.addWidget(self.preview_label)
        preview_group.setLayout(preview_layout)
        layout.addWidget(preview_group)

        # --- Buttons ---
        btn_layout = QHBoxLayout()
        self.btn_ok = QPushButton("Apply to Job")
        self.btn_ok.clicked.connect(self.accept)
        self.btn_cancel = QPushButton("Cancel")
        self.btn_cancel.clicked.connect(self.reject)
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_ok)
        btn_layout.addWidget(self.btn_cancel)
        layout.addLayout(btn_layout)

        self.setLayout(layout)
        
        self.connect_signals()
        self.ui_ready = True
        self.update_ui_state() # Initial UI state update
        self.update_preview()

    def setup_method_tab(self):
        layout = QFormLayout()

        # Method Type
        self.method_type = QComboBox()
        self.method_type.addItems([
            "DFT (GGA/Hybrid/Meta)", 
            "DFT (Range-Separated)", 
            "DFT (Double Hybrid)",
            "Wavefunction (HF/MP2)",
            "Wavefunction (Coupled Cluster)",
            "Wavefunction (Multireference)",
            "Semi-Empirical"
        ])
        self.method_type.currentIndexChanged.connect(self.update_method_list)
        layout.addRow("Method Type:", self.method_type)
        
        # Method Name
        self.method_name = QComboBox()
        self.update_method_list()
        layout.addRow("Method:", self.method_name)
        
        # Basis Set
        self.basis_set = QComboBox()
        basis_groups = [
            "--- Karlsruhe (Def2) ---",
            "def2-SV(P)", "def2-SVP", "def2-TZVP", "def2-QZVP",
            "ma-def2-SVP", "ma-def2-TZVP", "ma-def2-QZVP",
            "def2-TZVPP", "def2-QZVPP",
            "--- Dunning (cc-pV) ---",
            "cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z",
            "aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ", "aug-cc-pV5Z",
            "--- Pople ---",
            "6-31G", "6-31G*", "6-311G", "6-311G*", "6-311G**", 
            "6-31+G*", "6-311+G**", "6-31++G**",
            "--- Jensen (pc) ---",
            "pc-0", "pc-1", "pc-2", "pc-3", "aug-pc-1", "aug-pc-2",
            "--- Other ---",
            "EPR-II", "EPR-III", "IGLO-II", "IGLO-III"
        ]
        self.basis_set.addItems(basis_groups)
        self.basis_set.setCurrentText("def2-SVP")
        layout.addRow("Basis Set:", self.basis_set)
        
        # Auxiliary Basis (RI/RIJCOSX)
        self.aux_basis = QComboBox()
        self.aux_basis.addItems([
            "Auto (Def2/J, etc)", "None", "Def2/J", "Def2/JK", 
            "AutoAux", "NoAux"
        ])
        layout.addRow("Aux Basis (RI):", self.aux_basis)

        self.tab_method.setLayout(layout)

    def update_method_list(self):
        mtype = self.method_type.currentText()
        self.method_name.blockSignals(True)
        self.method_name.clear()
        
        if "GGA/Hybrid" in mtype:
             self.method_name.addItems([
                 "B3LYP", "PBE0", "PBE", "BP86", "BLYP", "PW91", 
                 "TPSSh", "TPSS", "SCAN", "M06", "M06-2X", "M06-HF", "M06-L",
                 "X3LYP", "O3LYP", "B3PW91", "BH&HLYP"
             ])
        elif "Range-Separated" in mtype:
             self.method_name.addItems([
                 "wB97X-D3", "wB97X-V", "wB97M-V", "wB97X", "wB97",
                 "CAM-B3LYP", "LC-wPBE", "wPBE", "wPBEh"
             ])
        elif "Double Hybrid" in mtype:
             self.method_name.addItems([
                 "B2PLYP", "mPW2PLYP", "B2GP-PLYP", "DSD-PBEP86", "DSD-BLYP",
                 "PTPSS-D3", "PWPB95"
             ])
        elif "Wavefunction (HF/MP2)" in mtype:
             self.method_name.addItems([
                 "HF", "UHF", "ROHF", 
                 "MP2", "RI-MP2", "SCS-MP2", "OO-RI-MP2"
             ])
        elif "Coupled Cluster" in mtype:
             self.method_name.addItems([
                 "DLPNO-CCSD(T)", "DLPNO-CCSD(T1)", "DLPNO-CCSD", 
                 "CCSD(T)", "CCSD", "QCISD(T)", "CEPA/1"
             ])
        elif "Multireference" in mtype:
             self.method_name.addItems([
                 "CASSCF", "NEVPT2", "DLPNO-NEVPT2", "MRCI"
             ])
        elif "Semi-Empirical" in mtype:
            self.method_name.addItems(["XT", "GFN1-xTB", "GFN2-xTB", "PM3", "AM1", "ZINDO/S", "PM6", "MNDO"])
        else:
            self.method_name.addItem("Custom")
            
        self.method_name.blockSignals(False)
        self.update_ui_state() # Update visibility when method type changes
        self.update_preview()

    def setup_job_tab(self):
        layout = QVBoxLayout()
        
        self.job_type = QComboBox()
        self.job_type.addItems([
            "Optimization + Freq (Opt Freq)", 
            "Optimization Only (Opt)", 
            "Frequency Only (Freq)", 
            "Single Point Energy (SP)",
            "Scan (Relaxed Surface)",
            "Transition State Opt (OptTS)",
            "Gradient",
            "Hessian"
        ])
        layout.addWidget(QLabel("Job Task:"))
        layout.addWidget(self.job_type)
        self.job_type.currentIndexChanged.connect(self.update_ui_state)
        
        # Opt Options
        self.opt_group = QGroupBox("Optimization Options")
        opt_layout = QHBoxLayout()
        self.opt_tight = QCheckBox("TightOpt")
        self.opt_loose = QCheckBox("LooseOpt")
        self.opt_calcfc = QCheckBox("CalcFC")
        self.opt_ts_mode = QCheckBox("CalcHess (for TS)")
        opt_layout.addWidget(self.opt_tight)
        opt_layout.addWidget(self.opt_loose)
        opt_layout.addWidget(self.opt_calcfc)
        opt_layout.addWidget(self.opt_ts_mode)
        self.opt_group.setLayout(opt_layout)
        layout.addWidget(self.opt_group)
        
        # Freq Options
        self.freq_group = QGroupBox("Freq Options")
        freq_layout = QHBoxLayout()
        self.freq_num = QCheckBox("NumFreq")
        self.freq_raman = QCheckBox("Raman")
        freq_layout.addWidget(self.freq_num)
        freq_layout.addWidget(self.freq_raman)
        self.freq_group.setLayout(freq_layout)
        layout.addWidget(self.freq_group)
        
        layout.addStretch()
        self.tab_job.setLayout(layout)

    def setup_solvation_tab(self):
        layout = QFormLayout()
        
        self.solv_model = QComboBox()
        self.solv_model.addItems(["None", "CPCM", "SMD", "IEFPCM", "CPC(Water) (Short)"])
        self.solv_model.currentIndexChanged.connect(self.update_ui_state)
        layout.addRow("Solvation Model:", self.solv_model)
        
        self.solvent = QComboBox()
        solvents = [
            "Water", "Acetonitrile", "Methanol", "Ethanol", 
            "Chloroform", "Dichloromethane", "Toluene", 
            "THF", "DMSO", "Cyclohexane", "Benzene", "Acetone",
            "CCl4", "DMF", "HMPA", "Pyridine"
        ]
        self.solvent.addItems(solvents)
        layout.addRow("Solvent:", self.solvent)
        
        layout.addRow(QLabel(" "))
        
        self.dispersion = QComboBox()
        self.dispersion.addItems(["None", "D3BJ", "D3Zero", "D4", "D2", "NL"])
        layout.addRow("Dispersion Correction:", self.dispersion)
        
        self.tab_solvation.setLayout(layout)

    def setup_props_tab(self):
        layout = QFormLayout()
        
        self.rijcosx = QCheckBox("RIJCOSX / RI approximation")
        self.rijcosx.setChecked(True)
        layout.addRow(self.rijcosx)
        
        self.grid_combo = QComboBox()
        self.grid_combo.addItems(["Default", "DefGrid2", "DefGrid3", "Grid4", "Grid5", "Grid6", "NoGrid"])
        layout.addRow("Grid:", self.grid_combo)

        self.scf_conv = QComboBox()
        self.scf_conv.addItems(["Default", "LooseSCF", "TightSCF", "VeryTightSCF"])
        layout.addRow("SCF Convergence:", self.scf_conv)

        # NBO
        self.pop_nbo = QCheckBox("NBO Analysis (! NBO)")
        layout.addRow(self.pop_nbo)

        self.tab_props.setLayout(layout)

    def connect_signals(self):
        widgets = [
            self.method_type, self.method_name, self.basis_set, self.aux_basis,
            self.job_type, self.opt_tight, self.opt_loose, self.opt_calcfc, self.opt_ts_mode,
            self.freq_num, self.freq_raman,
            self.solv_model, self.solvent, self.dispersion,
            self.rijcosx, self.grid_combo, self.scf_conv, self.pop_nbo
        ]
        for w in widgets:
            if isinstance(w, QComboBox):
                w.currentIndexChanged.connect(self.update_preview)
            elif isinstance(w, QCheckBox):
                w.toggled.connect(self.update_preview)
            elif isinstance(w, QSpinBox):
                w.valueChanged.connect(self.update_preview)

    def update_ui_state(self):
        """Update usability of widgets based on current selection."""
        if not getattr(self, 'ui_ready', False): return

        # 1. Method Dependent
        mtype = self.method_type.currentText()
        is_semi = "Semi-Empirical" in mtype
        
        # Disable Basis Set & Aux Basis for Semi-Empirical
        self.basis_set.setEnabled(not is_semi)
        self.aux_basis.setEnabled(not is_semi)
        
        # Handling RI / RIJCOSX
        if is_semi:
            self.rijcosx.setEnabled(False)
            self.rijcosx.setChecked(False)
            self.rijcosx.setText("RIJCOSX (N/A)")
        else:
            self.rijcosx.setEnabled(True)
            if "Wavefunction" in mtype:
                self.rijcosx.setText("RI Approximation (! RI ...)")
            else:
                self.rijcosx.setText("RIJCOSX (Speed up Hybrid DFT)")

        # 2. Solvation
        solv = self.solv_model.currentText()
        is_solvated = solv != "None"
        self.solvent.setEnabled(is_solvated)
        if "CPC(Water)" in solv:
             self.solvent.setEnabled(False) # Water is implied

        # 3. Job Type
        job_txt = self.job_type.currentText()
        is_opt = "Opt" in job_txt or "Scan" in job_txt
        is_freq = "Freq" in job_txt
        
        self.opt_group.setVisible(is_opt)
        self.freq_group.setVisible(is_freq)

        # 4. TD-DFT (Removed from Route Builder, handled via blocks)

    def update_preview(self):
        if not getattr(self, 'ui_ready', False):
            return

        route_parts = ["!"]
        
        # Method / Basis
        method = self.method_name.currentText()
        basis = self.basis_set.currentText()
        
        # 3c methods usually don't need basis set
        mtype = self.method_type.currentText()
        if "Semi-Empirical" in mtype:
            route_parts.append(method)
        elif "3c" in method:
            route_parts.append(method)
        else:
            route_parts.append(method)
            route_parts.append(basis)
            
            # RIJCOSX / RI
            if self.rijcosx.isEnabled() and self.rijcosx.isChecked():
                if "Wavefunction" in mtype:
                     route_parts.append("RI")
                else:
                     route_parts.append("RIJCOSX")
                     
                aux = self.aux_basis.currentText()
                if "Def2/J" in aux: route_parts.append("Def2/J")
                elif "Def2/JK" in aux: route_parts.append("Def2/JK")

        # Job Type
        job_txt = self.job_type.currentText()
        if "Opt Freq" in job_txt: route_parts.extend(["Opt", "Freq"])
        elif "Opt Only" in job_txt: route_parts.append("Opt")
        elif "OptTS" in job_txt: route_parts.append("OptTS")
        elif "Freq Only" in job_txt: route_parts.append("Freq")
        elif "Scan" in job_txt: route_parts.append("Scan")
        elif "Gradient" in job_txt: route_parts.append("Gradient")
        elif "Hessian" in job_txt: route_parts.append("Hessian")
        elif "SP" in job_txt: pass # No keyword
        
        # Opt Options
        if self.opt_group.isVisible():
            if self.opt_tight.isChecked(): route_parts.append("TightOpt")
            if self.opt_loose.isChecked(): route_parts.append("LooseOpt")
            if self.opt_calcfc.isChecked(): route_parts.append("CalcFC")
            if self.opt_ts_mode.isChecked(): route_parts.append("CalcHess")
        
        # Freq Options
        if self.freq_group.isVisible():
            if self.freq_num.isChecked(): route_parts.append("NumFreq")
        
        # Solvation
        solv = self.solv_model.currentText()
        if solv != "None":
            if "CPC(Water)" in solv:
                 route_parts.append("CPC(Water)")
            else:
                solvent = self.solvent.currentText()
                if "CPCM" == solv:
                    route_parts.append(f"CPCM({solvent})")
                elif "SMD" == solv:
                    route_parts.append(f"CPCM({solvent})")
                    route_parts.append("SMD")
                elif "IEFPCM" == solv:
                    route_parts.append(f"CPCM({solvent})") 

        # Dispersion
        disp = self.dispersion.currentText()
        if disp != "None":
            route_parts.append(disp)

        # SCF / Grid
        scf = self.scf_conv.currentText()
        if scf != "Default": route_parts.append(scf)
        
        grid = self.grid_combo.currentText()
        if grid != "Default": route_parts.append(grid)
        
        # NBO
        if self.pop_nbo.isChecked(): route_parts.append("NBO")

        self.preview_str = " ".join(route_parts)
        self.preview_label.setText(self.preview_str)

    def get_route(self):
        return self.preview_str
        
    def parse_route(self, route):
        pass

class OrcaSetupDialogNeo(QDialog):
    """
    ORCA Input Generator Neo
    """
    def __init__(self, parent=None, mol=None):
        super().__init__(parent)
        self.setWindowTitle(PLUGIN_NAME)
        self.resize(600, 700)
        self.mol = mol
        self.setup_ui()
        self.load_presets_from_file()
        self.calc_initial_charge_mult()

    def setup_ui(self):
        main_layout = QVBoxLayout()

        # --- 0. Preset Management ---
        preset_group = QGroupBox("Preset Management")
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
        main_layout.addWidget(preset_group)

        # --- 1. Resource Configuration ---
        res_group = QGroupBox("Resources (%pal, %maxcore)")
        res_layout = QFormLayout()

        # NProcs
        self.nproc_spin = QSpinBox()
        self.nproc_spin.setRange(1, 128)
        self.nproc_spin.setValue(4)
        res_layout.addRow("Number of Processors (nprocs):", self.nproc_spin)

        # Memory
        self.mem_spin = QSpinBox()
        self.mem_spin.setRange(100, 999999)
        self.mem_spin.setValue(2000) 
        self.mem_spin.setSuffix(" MB")
        res_layout.addRow("Memory per Core (MaxCore):", self.mem_spin)

        res_group.setLayout(res_layout)
        main_layout.addWidget(res_group)

        # --- 2. Simple Input Line ---
        kw_group = QGroupBox("Simple Input Line (!)")
        kw_layout = QVBoxLayout()
        
        kw_h_layout = QHBoxLayout()
        self.keywords_edit = QLineEdit("! B3LYP def2-SVP Opt Freq")
        self.btn_route = QPushButton("Keyword Builder...")
        self.btn_route.clicked.connect(self.open_keyword_builder)
        kw_h_layout.addWidget(self.keywords_edit)
        kw_h_layout.addWidget(self.btn_route)
        
        kw_layout.addWidget(QLabel("Keywords (starts with !):"))
        kw_layout.addLayout(kw_h_layout)
        
        # Comment
        self.comment_edit = QLineEdit("Generated by MoleditPy ORCA Neo")
        kw_layout.addWidget(QLabel("Comment (# ...):"))
        kw_layout.addWidget(self.comment_edit)

        kw_group.setLayout(kw_layout)
        main_layout.addWidget(kw_group)

        # --- 3. Molecular State ---
        mol_group = QGroupBox("Molecular Specification")
        mol_layout = QHBoxLayout()
        
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.valueChanged.connect(self.validate_charge_mult)
        
        self.mult_spin = QSpinBox()
        self.mult_spin.setRange(1, 10)
        self.mult_spin.valueChanged.connect(self.validate_charge_mult)

        mol_layout.addWidget(QLabel("Charge:"))
        mol_layout.addWidget(self.charge_spin)
        mol_layout.addWidget(QLabel("Multiplicity:"))
        mol_layout.addWidget(self.mult_spin)

        self.default_palette = self.charge_spin.palette()
        
        mol_group.setLayout(mol_layout)
        main_layout.addWidget(mol_group)

        # --- 4. Advanced/Blocks ---
        adv_group = QGroupBox("Advanced Blocks (Pre/Post Coordinates)")
        adv_layout = QVBoxLayout()
        
        # Block Helper
        blk_h_layout = QHBoxLayout()
        blk_h_layout.addWidget(QLabel("Block Helper:"))
        
        self.block_combo = QComboBox()
        self.block_combo.addItems([
             "Select Block to Insert...",
             "%scf ... end",
             "%output ... end",
             "%geom ... end",
             "%elprop ... end",
             "%plots ... end",
             "%tddft ... end",
             "%cis ... end", 
             "%mrci ... end",
             "%casscf ... end",
             "%rr ... end"
        ])
        blk_h_layout.addWidget(self.block_combo, 1)
        
        self.btn_insert_block = QPushButton("Insert")
        self.btn_insert_block.clicked.connect(self.insert_block_template)
        blk_h_layout.addWidget(self.btn_insert_block)
        
        adv_layout.addLayout(blk_h_layout)
        
        self.adv_edit = QTextEdit()
        self.adv_edit.setPlaceholderText("Example:\n%scf\n MaxIter 100\nend\n\n%output\n Print[P_Hirshfeld] 1\nend")
        adv_layout.addWidget(self.adv_edit)

        adv_group.setLayout(adv_layout)
        main_layout.addWidget(adv_group)

        # --- Save Button ---
        self.save_btn = QPushButton("Save Input File...")
        self.save_btn.clicked.connect(self.save_file)
        self.save_btn.setStyleSheet("font-weight: bold; padding: 5px;")
        main_layout.addWidget(self.save_btn)
        
        self.setLayout(main_layout)

    def insert_block_template(self):
        txt = self.block_combo.currentText()
        if "Select" in txt: return
        
        template = ""
        if "%scf" in txt:
             template = "%scf\n MaxIter 125\nend\n"
        elif "%output" in txt:
             template = "%output\n Print[P_Hirshfeld] 1\nend\n"
        elif "%geom" in txt:
             template = "%geom\n MaxIter 100\nend\n"
        elif "%elprop" in txt:
             template = "%elprop\n Dipole True\n Quadrupole True\nend\n"
        elif "%plots" in txt:
             template = "%plots\n Format Gaussian_Cube\nend\n"
        elif "%tddft" in txt:
             template = (
                 "%tddft\n"
                 "  NRoots 10       # Number of excited states\n"
                 "  MaxDim 10       # Max dimension of expansion space\n"
                 "  TDA true        # Tamm-Dancoff Approximation (true/false)\n"
                 "  IRoot 1         # State of interest for gradient properties\n"
                 "  Triplets true   # Calculate triplet states\n"
                 "  Singlets true   # Calculate singlet states\n"
                 "  DoQuad true     # Compute quadrupole intensities\n"
                 "end\n"
             )
        elif "%cis" in txt:
             template = "%cis\n NRoots 10\nend\n"
        elif "%mrci" in txt:
             template = "%mrci\n NewBlocks 1 1\nend\n"
        elif "%casscf" in txt:
             template = "%casscf\n Nel 2\n Norb 2\n Mult 1\nend\n"
        
        cursor = self.adv_edit.textCursor()
        cursor.insertText(template)

    def open_keyword_builder(self):
        dialog = OrcaKeywordBuilderDialog(self, self.keywords_edit.text())
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.keywords_edit.setText(dialog.get_route())

    def get_coords_lines(self):
        if not self.mol: return []
        lines = []
        try:
            conf = self.mol.GetConformer()
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                atom = self.mol.GetAtomWithIdx(i)
                lines.append(f"  {atom.GetSymbol(): <4} {pos.x: >12.6f} {pos.y: >12.6f} {pos.z: >12.6f}")
        except Exception as e:
            return [f"# Error: {e}"]
        return lines

    def save_file(self):
        coord_lines = self.get_coords_lines()
        if any("Error" in line for line in coord_lines):
            QMessageBox.critical(self, "Error", "\n".join(coord_lines))
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save ORCA Input", "", "ORCA Input (*.inp);;All Files (*)"
        )

        if file_path:
            try:
                with open(file_path, 'w', encoding='utf-8') as f:
                    comment = self.comment_edit.text().strip()
                    if comment:
                         f.write(f"# {comment}\n")
                    
                    # Resources
                    nprocs = self.nproc_spin.value()
                    if nprocs > 1:
                        f.write(f"%pal nprocs {nprocs} end\n")
                    f.write(f"%maxcore {self.mem_spin.value()}\n\n")
                    
                    # Keywords
                    kw = self.keywords_edit.text().strip()
                    if not kw.startswith("!"): kw = "! " + kw
                    f.write(f"{kw}\n\n")
                    
                    # Advanced Blocks
                    adv = self.adv_edit.toPlainText().strip()
                    if adv:
                        f.write(f"{adv}\n\n")
                    
                    # Coordinates
                    f.write(f"* xyz {self.charge_spin.value()} {self.mult_spin.value()}\n")
                    for l in coord_lines:
                        f.write(l + "\n")
                    f.write("*\n")

                QMessageBox.information(self, "Success", f"File saved:\n{file_path}")
                self.accept()
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save file:\n{e}")

    # --- Preset Management (Similar to Gaussian Neo) ---
    def load_presets_from_file(self):
        self.presets_data = {}
        if os.path.exists(SETTINGS_FILE):
            try:
                with open(SETTINGS_FILE, 'r', encoding='utf-8') as f:
                    self.presets_data = json.load(f)
            except Exception as e:
                print(f"Error loading presets: {e}")
        
        if "Default" not in self.presets_data:
            self.presets_data["Default"] = {
                "nproc": 4, "maxcore": 2000, 
                "route": "! B3LYP def2-SVP Opt Freq", "adv": ""
            }
        
        self.update_preset_combo()

    def update_preset_combo(self):
        current = self.preset_combo.currentText()
        self.preset_combo.blockSignals(True)
        self.preset_combo.clear()
        self.preset_combo.addItems(sorted(self.presets_data.keys()))
        
        if current in self.presets_data:
            self.preset_combo.setCurrentText(current)
        elif "Default" in self.presets_data:
            self.preset_combo.setCurrentText("Default")
        
        self.preset_combo.blockSignals(False)
        self.apply_selected_preset()

    def apply_selected_preset(self):
        name = self.preset_combo.currentText()
        if name not in self.presets_data: return
        data = self.presets_data[name]
        
        self.nproc_spin.setValue(data.get("nproc", 4))
        self.mem_spin.setValue(data.get("maxcore", 2000))
        self.keywords_edit.setText(data.get("route", "! B3LYP def2-SVP Opt Freq"))
        self.adv_edit.setPlainText(data.get("adv", ""))

    def save_preset_dialog(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:")
        if ok and name:
            self.presets_data[name] = {
                "nproc": self.nproc_spin.value(),
                "maxcore": self.mem_spin.value(),
                "route": self.keywords_edit.text(),
                "adv": self.adv_edit.toPlainText()
            }
            self.save_presets_to_file()
            self.update_preset_combo()
            self.preset_combo.setCurrentText(name)

    def delete_preset(self):
        name = self.preset_combo.currentText()
        if name == "Default":
             QMessageBox.warning(self, "Warning", "Cannot delete Default preset.")
             return
        
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

    # --- Charge/Mult Logic ---
    def calc_initial_charge_mult(self):
        if not self.mol: return
        try:
            try: charge = Chem.GetFormalCharge(self.mol)
            except: charge = 0
            
            num_radical = sum(atom.GetNumRadicalElectrons() for atom in self.mol.GetAtoms())
            mult = num_radical + 1
            
            self.charge_spin.setValue(int(charge))
            self.mult_spin.setValue(int(mult))
            self.validate_charge_mult()
        except Exception: pass

    def validate_charge_mult(self):
        if not self.mol: return
        try:
            charge = self.charge_spin.value()
            mult = self.mult_spin.value()
            
            total_protons = sum(atom.GetAtomicNum() for atom in self.mol.GetAtoms())
            total_electrons = total_protons - charge
            
            is_valid = (total_electrons % 2 == 0 and mult % 2 != 0) or \
                       (total_electrons % 2 != 0 and mult % 2 == 0)
                       
            if is_valid:
                self.charge_spin.setPalette(self.default_palette)
                self.mult_spin.setPalette(self.default_palette)
            else:
                p = self.charge_spin.palette()
                p.setColor(QPalette.ColorRole.Base, QColor("#FFDDDD"))
                p.setColor(QPalette.ColorRole.Text, Qt.GlobalColor.red)
                self.charge_spin.setPalette(p)
                self.mult_spin.setPalette(p)
        except: pass

from PyQt6.QtWidgets import QInputDialog 

def run(main_window):
    mol = getattr(main_window, 'current_mol', None)
    if not mol:
        QMessageBox.warning(main_window, PLUGIN_NAME, "No molecule loaded.")
        return
    dialog = OrcaSetupDialogNeo(parent=main_window, mol=mol)
    dialog.exec()
