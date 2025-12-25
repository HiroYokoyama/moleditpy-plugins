# -*- coding: utf-8 -*-
import os
from PyQt6.QtWidgets import (QMessageBox, QDialog, QVBoxLayout, QLabel, 
                             QLineEdit, QSpinBox, QPushButton, QFileDialog, 
                             QFormLayout, QGroupBox, QHBoxLayout, QComboBox, QTextEdit, 
                             QTabWidget, QCheckBox, QWidget, QScrollArea, QMenu)
from PyQt6.QtGui import QPalette, QColor, QAction, QFont
from PyQt6.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import json

__version__="2025.12.25"
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
                 "TPSSh", "TPSS", "SCAN", "r2SCAN-3c", "B97-3c", "PBEh-3c", # 3c系を追加
                 "M06", "M06-2X", "M06-HF", "M06-L",
                 "X3LYP", "O3LYP", "B3PW91", "BHandHLYP" # BH&HLYPをBHandHLYPに修正
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
                 "HF", "HF-3c", "UHF", "ROHF", # HF-3cを追加
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
    def __init__(self, parent=None, mol=None, filename=None):
        super().__init__(parent)
        self.setWindowTitle(PLUGIN_NAME)
        self.resize(600, 700)
        self.mol = mol
        self.filename = filename
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

        # --- 3b. Coordinate Format ---
        coord_group = QGroupBox("Coordinate Format")
        coord_layout = QHBoxLayout()
        self.coord_format_combo = QComboBox()
        self.coord_format_combo.addItems(["Cartesian (XYZ)", "Internal (Z-Matrix / * int)", "Internal (Z-Matrix / * gzmt)"])
        self.coord_format_combo.setCurrentIndex(0)
        self.coord_format_combo.setEnabled(True)
        coord_layout.addWidget(QLabel("Format:"))
        coord_layout.addWidget(self.coord_format_combo)
        coord_group.setLayout(coord_layout)
        main_layout.addWidget(coord_group)

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
             "%eprnmr ... end (Post)",
        ])
        blk_h_layout.addWidget(self.block_combo, 1)
        
        self.btn_insert_block = QPushButton("Insert")
        self.btn_insert_block.clicked.connect(self.insert_block_template)
        blk_h_layout.addWidget(self.btn_insert_block)
        
        adv_layout.addLayout(blk_h_layout)
        
        # Tabs for Pre/Post
        self.adv_tabs = QTabWidget()
        
        self.adv_edit = QTextEdit()
        self.adv_edit.setPlaceholderText("Pre-Coordinate Blocks\nExample:\n%scf\n MaxIter 100\nend")
        self.adv_tabs.addTab(self.adv_edit, "Pre-Coordinate")
        
        self.post_adv_edit = QTextEdit()
        self.post_adv_edit.setPlaceholderText("Post-Coordinate Blocks")
        self.adv_tabs.addTab(self.post_adv_edit, "Post-Coordinate")
        
        adv_layout.addWidget(self.adv_tabs)

        adv_group.setLayout(adv_layout)
        main_layout.addWidget(adv_group)

        # --- Save/Preview Buttons ---
        btn_box = QHBoxLayout()
        
        self.btn_preview = QPushButton("Preview Input")
        self.btn_preview.clicked.connect(self.preview_file)
        btn_box.addWidget(self.btn_preview)
        
        self.save_btn = QPushButton("Save Input File...")
        self.save_btn.clicked.connect(self.save_file)
        self.save_btn.setStyleSheet("font-weight: bold; padding: 5px;")
        btn_box.addWidget(self.save_btn)
        
        main_layout.addLayout(btn_box)
        
        self.setLayout(main_layout)

    def insert_block_template(self):
        txt = self.block_combo.currentText()
        if "Select" in txt: return
        
        template = ""
        if "%scf" in txt:
             template = "%scf\n MaxIter 125\nend\n"
        elif "%output" in txt:
             template = "%output\n  PrintLevel     Normal\nend\n"
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
        elif "%eprnmr" in txt:
             template = "%eprnmr\n     nuclei = all h {shift, ssall}\nend\n"
             # Switch to Post-Coordinate tab automatically
             self.adv_tabs.setCurrentWidget(self.post_adv_edit)
        
        current_widget = self.adv_tabs.currentWidget()
        if isinstance(current_widget, QTextEdit):
            cursor = current_widget.textCursor()
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

    def _build_zmatrix_data(self):
        """Helper to build Z-Matrix connectivity and values."""
        if not self.mol: return None
        try:
            atoms = list(self.mol.GetAtoms())
            conf = self.mol.GetConformer()
            
            def get_dist(i, j): return rdMolTransforms.GetBondLength(conf, i, j)
            def get_angle(i, j, k): return rdMolTransforms.GetAngleDeg(conf, i, j, k)
            def get_dihedral(i, j, k, l): return rdMolTransforms.GetDihedralDeg(conf, i, j, k, l)
            
            defined = []
            z_data = [] # List of dicts for each atom
            
            for i, atom in enumerate(atoms):
                symbol = atom.GetSymbol()
                
                # Atom 0
                if i == 0:
                    z_data.append({"symbol": symbol, "refs": []})
                    defined.append(i)
                    continue
                
                # Find neighbors in defined set
                current_idx = atom.GetIdx()
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                candidates = [n for n in neighbors if n in defined]
                if not candidates: candidates = defined[:] # Fallback
                
                refs = []
                # Ref 1 (Distance)
                if candidates: refs.append(candidates[-1])
                else: refs.append(0)
                
                # Ref 2 (Angle)
                candidates_2 = [x for x in defined if x != refs[0]]
                if candidates_2: refs.append(candidates_2[-1])
                
                # Ref 3 (Dihedral)
                candidates_3 = [x for x in defined if x not in refs]
                if candidates_3: refs.append(candidates_3[-1])
                
                # Calculate Values
                row = {"symbol": symbol, "refs": [], "values": []}
                
                if len(refs) >= 1:
                    row["refs"].append(refs[0]) # 0-based index for calculation
                    row["values"].append(get_dist(i, refs[0]))
                
                if len(refs) >= 2:
                    row["refs"].append(refs[1])
                    row["values"].append(get_angle(i, refs[0], refs[1]))
                
                if len(refs) >= 3:
                    row["refs"].append(refs[2])
                    row["values"].append(get_dihedral(i, refs[0], refs[1], refs[2]))
                
                z_data.append(row)
                defined.append(i)
                
            return z_data
        except Exception as e:
            raise e

    def get_zmatrix_standard_lines(self):
        """
        Generates * int style lines.
        Format: Symbol Ref1 Ref2 Ref3 R Angle Dihed
        Refs are 1-based integers. Missing refs/values are 0/0.0.
        """
        try:
            data = self._build_zmatrix_data()
            if not data: return []
            
            lines = []
            for i, row in enumerate(data):
                symbol = row["symbol"]
                
                # Prepare 3 refs (1-based) and 3 values
                refs_out = [0, 0, 0]
                vals_out = [0.0, 0.0, 0.0]
                
                current_refs = row["refs"]
                current_vals = row.get("values", [])
                
                for k in range(min(3, len(current_refs))):
                    refs_out[k] = current_refs[k] + 1 # Convert to 1-based
                    vals_out[k] = current_vals[k]
                    
                line = f"  {symbol: <3} {refs_out[0]: >3} {refs_out[1]: >3} {refs_out[2]: >3} " \
                       f"{vals_out[0]: >10.6f} {vals_out[1]: >10.6f} {vals_out[2]: >10.6f}"
                lines.append(line)
            return lines
        except Exception as e:
            return [f"# Error generating Standard Z-Matrix: {e}"]

    def get_zmatrix_gzmt_lines(self):
        """
        Generates * gzmt style lines (Compact).
        Format: Symbol Ref1 R Ref2 Angle Ref3 Dihed
        Refs are 1-based.
        """
        try:
            data = self._build_zmatrix_data()
            if not data: return []
            
            lines = []
            for i, row in enumerate(data):
                symbol = row["symbol"]
                line = f"  {symbol: <3}"
                
                current_refs = row["refs"]
                current_vals = row.get("values", [])
                
                # Z-Matrix logic:
                # Atom 1: Symbol
                # Atom 2: Symbol Ref1 R
                # Atom 3: Symbol Ref1 R Ref2 A
                # Atom 4: Symbol Ref1 R Ref2 A Ref3 D
                
                if i == 0:
                    pass
                else:
                    count = len(current_refs)
                    if count >= 1:
                        line += f"  {current_refs[0] + 1: >3} {current_vals[0]: .6f}"
                    if count >= 2:
                        line += f"  {current_refs[1] + 1: >3} {current_vals[1]: .6f}"
                    if count >= 3:
                        line += f"  {current_refs[2] + 1: >3} {current_vals[2]: .6f}"
                
                lines.append(line)
            return lines
        except Exception as e:
            return [f"# Error generating GZMT Z-Matrix: {e}"]


    def generate_input_content(self):
        """Generates the full content of the input file as a string."""
        content = []
        
        comment = self.comment_edit.text().strip()
        if comment:
                content.append(f"# {comment}")
        
        # Resources
        nprocs = self.nproc_spin.value()
        if nprocs > 1:
            content.append(f"%pal nprocs {nprocs} end")
        content.append(f"%maxcore {self.mem_spin.value()}\n")
        
        # Keywords
        kw = self.keywords_edit.text().strip()
        if not kw.startswith("!"): kw = "! " + kw
        content.append(f"{kw}\n")
        
        # Advanced Blocks
        adv = self.adv_edit.toPlainText().strip()
        if adv:
            content.append(f"{adv}\n")
        
        # Coordinates
        is_cartesian = "Cartesian" in self.coord_format_combo.currentText()
        coord_lines = self.get_coords_lines()
        
        if is_cartesian:
            content.append(f"* xyz {self.charge_spin.value()} {self.mult_spin.value()}")
            content.extend(coord_lines)
            content.append("*")
        else:
            # Z-Matrix
            is_gzmt = "gzmt" in self.coord_format_combo.currentText()
            
            if is_gzmt:
                zmat_lines = self.get_zmatrix_gzmt_lines()
                header = "* gzmt"
            else:
                zmat_lines = self.get_zmatrix_standard_lines()
                header = "* int"
                
            if any("Error" in line for line in zmat_lines):
                    content.append(f"# ERROR: Z-Matrix generation failed.")
                    content.append(f"* xyz {self.charge_spin.value()} {self.mult_spin.value()}")
                    content.extend(self.get_coords_lines())
                    content.append("*")
            else:
                    content.append(f"{header} {self.charge_spin.value()} {self.mult_spin.value()}")
                    content.extend(zmat_lines)
                    content.append("*")
                    
        # Post-Coordinate Blocks
        adv_post = self.post_adv_edit.toPlainText().strip()
        if adv_post:
            content.append(f"\n{adv_post}")
            
        return "\n".join(content)

    def preview_file(self):
        content = self.generate_input_content()
        
        d = QDialog(self)
        d.setWindowTitle("Preview Input - ORCA Neo")
        d.resize(600, 500)
        l = QVBoxLayout()
        t = QTextEdit()
        t.setPlainText(content)
        t.setReadOnly(True)
        t.setFont(QFont("Courier New", 10))
        l.addWidget(t)
        
        btn = QPushButton("Close")
        btn.clicked.connect(d.accept)
        l.addWidget(btn)
        d.setLayout(l)
        d.exec()

    def save_file(self):
        # Validation Check (e.g. Z-Matrix Error)
        # We can re-check here or just trust generate_input_content
        
        default_name = ""
        if self.filename:
            base, _ = os.path.splitext(self.filename)
            default_name = base + ".inp"
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save ORCA Input", default_name, "ORCA Input (*.inp);;All Files (*)"
        )

        if file_path:
            try:
                content = self.generate_input_content()
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(content)

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
                "route": "! B3LYP def2-SVP Opt Freq", "adv": "", "adv_post": ""
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
        self.post_adv_edit.setPlainText(data.get("adv_post", ""))

    def save_preset_dialog(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:")
        if ok and name:
            self.presets_data[name] = {
                "nproc": self.nproc_spin.value(),
                "maxcore": self.mem_spin.value(),
                "route": self.keywords_edit.text(),
                "adv": self.adv_edit.toPlainText(),
                "adv_post": self.post_adv_edit.toPlainText()
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

def run(mw):
    mol = getattr(mw, 'current_mol', None)
    if not mol:
        QMessageBox.warning(mw, PLUGIN_NAME, "No molecule loaded.")
        return
        
    filename = getattr(mw, 'current_file_path', None)
    dialog = OrcaSetupDialogNeo(parent=mw, mol=mol, filename=filename)
    dialog.exec()

# initialize removed as it only registered the menu action

