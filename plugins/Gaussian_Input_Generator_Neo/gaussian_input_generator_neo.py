# -*- coding: utf-8 -*-
import os
from PyQt6.QtWidgets import (QMessageBox, QDialog, QVBoxLayout, QLabel, 
                             QLineEdit, QSpinBox, QPushButton, QFileDialog, 
                             QFormLayout, QGroupBox, QHBoxLayout, QComboBox, QTextEdit, 
                             QInputDialog, QTabWidget, QCheckBox, QRadioButton, QButtonGroup,
                             QWidget, QScrollArea, QSizePolicy)
from PyQt6.QtGui import QPalette, QColor, QFont
from PyQt6.QtCore import Qt
from rdkit import Chem
import json

PLUGIN_NAME = "Gaussian Input Generator Neo"
PLUGIN_VERSION = "2026.01.03"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Advanced Gaussian Input Generator with Preview and Presets"
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "gaussian_input_generator_neo.json")

class RouteBuilderDialog(QDialog):
    """
    Dialog to construct the Gaussian Job Route line.
    """
    def __init__(self, parent=None, current_route=""):
        super().__init__(parent)
        self.setWindowTitle("Route Builder")
        self.resize(600, 500)
        self.ui_ready = False  # Flag to prevent premature updates
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
        preview_group = QGroupBox("Route Preview")
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
        
        # Connect signals to update preview
        self.connect_signals()
        
        self.ui_ready = True  # Enable updates
        self.update_preview()

    def setup_method_tab(self):
        layout = QFormLayout()

        # Output Level
        self.output_level = QComboBox()
        self.output_level.addItems(["Additional Output (#P)", "Standard Output (#)", "Terse Output (#T)"])
        # Default to # (Standard)
        self.output_level.setCurrentIndex(1) 
        self.output_level.currentIndexChanged.connect(self.update_preview)
        layout.addRow("Print Level:", self.output_level)
        
        # Method Type
        self.method_type = QComboBox()
        self.method_type.addItems(["DFT", "Double Hybrid", "MP2", "Hartree-Fock", "Semi-Empirical", "Other"])
        self.method_type.currentIndexChanged.connect(self.update_method_list)
        layout.addRow("Method Type:", self.method_type)
        
        # Method Name
        self.method_name = QComboBox()
        # Initial population
        # Note: update_method_list calls update_preview, which calls currentText()
        # We need to make sure update_preview guards against uninitialized widgets
        self.update_method_list() 
        layout.addRow("Method:", self.method_name)
        
        # Basis Set
        self.basis_set = QComboBox()
        basis_groups = [
            "6-31G(d)", "6-31+G(d,p)", "6-311G(d,p)", "6-311+G(d,p)", # Pople
            "cc-pVDZ", "cc-pVTZ", "aug-cc-pVDZ", "aug-cc-pVTZ", # Dunning
            "def2SVP", "def2TZVP", "def2QZVP", # Karlsruhe
            "LanL2DZ", "SDD", "Gen", "GenECP" # ECP/Gen
        ]
        self.basis_set.addItems(basis_groups)
        layout.addRow("Basis Set:", self.basis_set)
        
        self.tab_method.setLayout(layout)

    def update_method_list(self):
        mtype = self.method_type.currentText()
        self.method_name.blockSignals(True)
        self.method_name.clear()
        
        if mtype == "DFT":
            self.method_name.addItems([
                "B3LYP", "WB97XD", "M062X", "PBE0", "CAM-B3LYP", "APFD", "B97D3", "TPSSTPSS",
                "MN15", "MN15L", "BHandHLYP"
            ])
        elif mtype == "Double Hybrid":
            self.method_name.addItems(["B2PLYP"])
        elif mtype == "MP2":
            self.method_name.addItems(["MP2", "MP3", "MP4", "CCSD", "CCSD(T)"])
        elif mtype == "Hartree-Fock":
            self.method_name.addItems(["HF", "ROHF", "UHF"])
        elif mtype == "Semi-Empirical":
            self.method_name.addItems(["AM1", "PM6", "PM7", "ZINDO"])
        else:
            self.method_name.addItem("Custom")
            
        self.method_name.blockSignals(False)
        self.update_preview()

    def setup_job_tab(self):
        layout = QVBoxLayout()
        
        self.job_type = QComboBox()
        self.job_type.addItems([
            "Optimization + Freq (Opt Freq)", 
            "Optimization Only (Opt)", 
            "Frequency Only (Freq)", 
            "Single Point Energy (SP)",
            "Scan (ModRedundant)",
            "IRC",
            "Stability Analysis (Stable)",
            "Volume"
        ])
        layout.addWidget(QLabel("Job Task:"))
        layout.addWidget(self.job_type)
        self.job_type.currentIndexChanged.connect(self.update_job_options_visibility)
        
        # Opt Options
        self.opt_group = QGroupBox("Optimization Options")
        opt_layout = QHBoxLayout()
        self.opt_tight = QCheckBox("Tight")
        self.opt_verytight = QCheckBox("VeryTight")
        self.opt_calcfc = QCheckBox("CalcFC")
        self.opt_maxcycles = QCheckBox("MaxCycles=100")
        opt_layout.addWidget(self.opt_tight)
        opt_layout.addWidget(self.opt_verytight)
        opt_layout.addWidget(self.opt_calcfc)
        opt_layout.addWidget(self.opt_maxcycles)
        self.opt_group.setLayout(opt_layout)
        layout.addWidget(self.opt_group)
        
        # Freq Options
        self.freq_group = QGroupBox("Freq Options")
        freq_layout = QHBoxLayout()
        self.freq_raman = QCheckBox("Raman")
        self.freq_anharm = QCheckBox("Anharmonic")
        freq_layout.addWidget(self.freq_raman)
        freq_layout.addWidget(self.freq_anharm)
        self.freq_group.setLayout(freq_layout)
        layout.addWidget(self.freq_group)
        
        layout.addStretch()
        self.tab_job.setLayout(layout)
        
        # Initial visibility update
        self.update_job_options_visibility()

    def update_job_options_visibility(self):
        job_idx = self.job_type.currentIndex()
        txt = self.job_type.currentText()
        
        # Opt options: Show for Opt, Opt+Freq, Scan, IRC, Stable, Volume (some imply optimization)
        # Strictly speaking: Opt tasks. 
        # 0: Opt+Freq, 1: Opt, 4: Scan (ModRedundant), 5: IRC
        is_opt = job_idx in [0, 1, 4, 5] 
        self.opt_group.setVisible(is_opt)
        
        # Freq options: Show for Opt+Freq, Freq
        # 0: Opt+Freq, 2: Freq
        is_freq = job_idx in [0, 2]
        self.freq_group.setVisible(is_freq)

    def setup_solvation_tab(self):
        layout = QFormLayout()
        
        self.solv_model = QComboBox()
        self.solv_model.addItems(["None", "PCM", "CPCM", "SMD", "IEFPCM"])
        layout.addRow("Solvation Model:", self.solv_model)
        
        self.solvent = QComboBox()
        solvents = [
            "Water", "Acetonitrile", "Methanol", "Ethanol", 
            "Chloroform", "Dichloromethane", "Toluene", 
            "Tetrahydrofuran", "DimethylSulfoxide", "Cyclohexane", 
            "Benzene", "Acetone"
        ]
        self.solvent.addItems(solvents)
        layout.addRow("Solvent:", self.solvent)
        
        layout.addRow(QLabel(" ")) # Spacer
        
        self.dispersion = QCheckBox("EmpiricalDispersion=GD3BJ")
        layout.addRow("Dispersion:", self.dispersion)
        
        self.tab_solvation.setLayout(layout)

    def setup_props_tab(self):
        layout = QFormLayout()
        
        self.pop_analysis = QComboBox()
        self.pop_analysis.addItems(["None", "NBO (Pop=NBO)", "Hirshfeld", "MK (Merz-Kollman)", "Regular (Pop=Reg)"])
        layout.addRow("Population Analysis:", self.pop_analysis)
        
        self.density_chk = QCheckBox("Density=Current")
        layout.addRow(QLabel("Density:"), self.density_chk)
        
        self.symmetry_combo = QComboBox()
        self.symmetry_combo.addItems(["Default", "Loose", "None (NoSymm)"])
        layout.addRow("Symmetry:", self.symmetry_combo)
        
        self.grid_combo = QComboBox()
        self.grid_combo.addItems(["Default", "FineGrid", "UltraFine", "SuperFine"])
        layout.addRow("Integration Grid:", self.grid_combo)

        # TD-DFT Section
        td_group = QGroupBox("Excited States (TD-DFT)")
        td_layout = QHBoxLayout()
        self.td_chk = QCheckBox("Enable TD")
        self.td_nstates = QSpinBox()
        self.td_nstates.setValue(6)
        self.td_nstates.setPrefix("NStates=")
        td_layout.addWidget(self.td_chk)
        td_layout.addWidget(self.td_nstates)
        td_group.setLayout(td_layout)
        layout.addRow(td_group)
        
        # WFN Output
        self.wfn_chk = QCheckBox("Output Wavefunction File (.wfn)")
        layout.addRow("WFN Output:", self.wfn_chk)
        self.wfn_chk.toggled.connect(self.update_preview)

        self.tab_props.setLayout(layout)

    def connect_signals(self):
        # Connect all widgets to update_preview
        widgets = [
            self.method_type, self.method_name, self.basis_set,
            self.job_type, self.opt_tight, self.opt_verytight, self.opt_calcfc, self.opt_maxcycles,
            self.freq_raman, self.freq_anharm,
            self.solv_model, self.solvent, self.dispersion,
            self.pop_analysis, self.density_chk, self.symmetry_combo, self.grid_combo,
            self.td_chk, self.td_nstates, self.wfn_chk
        ]
        for w in widgets:
            if isinstance(w, QComboBox):
                w.currentIndexChanged.connect(self.update_preview)
            elif isinstance(w, QCheckBox):
                w.toggled.connect(self.update_preview)
            elif isinstance(w, QSpinBox):
                w.valueChanged.connect(self.update_preview)

    def update_preview(self):
        if not getattr(self, 'ui_ready', False):
            return

        # Output Level
        lvl_map = {0: "#P", 1: "#", 2: "#T"}
        prefix = lvl_map.get(self.output_level.currentIndex(), "#P")
        route_parts = [prefix]
        
        # Method / Basis
        method = self.method_name.currentText()
        basis = self.basis_set.currentText()
        route_parts.append(f"{method}/{basis}")
        
        # Job Type
        job_idx = self.job_type.currentIndex()
        if job_idx == 0: route_parts.extend(["Opt", "Freq"])
        elif job_idx == 1: route_parts.append("Opt")
        elif job_idx == 2: route_parts.append("Freq")
        elif job_idx == 3: route_parts.append("SP")
        elif job_idx == 4: route_parts.append("Scan") # ModRedundant usually implies Opt=ModRedundant but simple Scan is obsolete
        elif job_idx == 5: route_parts.append("IRC")
        elif job_idx == 6: route_parts.append("Stable")
        elif job_idx == 7: route_parts.append("Volume")
        
        # Opt Options
        opt_opts = []
        if self.opt_tight.isChecked(): opt_opts.append("Tight")
        if self.opt_verytight.isChecked(): opt_opts.append("VeryTight")
        if self.opt_calcfc.isChecked(): opt_opts.append("CalcFC")
        if self.opt_maxcycles.isChecked(): opt_opts.append("MaxCycles=100")
        
        # If Opt is in route and we have options, combine them
        if "Opt" in route_parts and opt_opts:
             idx = route_parts.index("Opt")
             route_parts[idx] = f"Opt=({', '.join(opt_opts)})"
        
        # Freq Options
        freq_opts = []
        if self.freq_raman.isChecked(): freq_opts.append("Raman")
        if self.freq_anharm.isChecked(): freq_opts.append("Anharmonic")
        
        if "Freq" in route_parts and freq_opts:
             idx = route_parts.index("Freq")
             route_parts[idx] = f"Freq=({', '.join(freq_opts)})"

        # Solvation
        solv = self.solv_model.currentText()
        if solv != "None":
            solvent = self.solvent.currentText()
            scrf_opts = [f"{solv}", f"Solvent={solvent}"]
            route_parts.append(f"SCRF=({', '.join(scrf_opts)})")
            
        # Dispersion
        if self.dispersion.isChecked():
            route_parts.append("EmpiricalDispersion=GD3BJ")
            
        # Properties
        pop = self.pop_analysis.currentText()
        if "NBO" in pop: route_parts.append("Pop=NBO")
        elif "Hirshfeld" in pop: route_parts.append("Pop=Hirshfeld")
        elif "MK" in pop: route_parts.append("Pop=MK")
        elif "Reg" in pop: route_parts.append("Pop=Reg")
        
        if self.density_chk.isChecked():
            route_parts.append("Density=Current")
            
        sym = self.symmetry_combo.currentText()
        if "Loose" in sym: route_parts.append("Symmetry=Loose")
        elif "None" in sym: route_parts.append("NoSymm")
        
        grid = self.grid_combo.currentText()
        if grid != "Default":
             route_parts.append(f"Integral({grid})")
             
        # TD-DFT
        if self.td_chk.isChecked():
            route_parts.append(f"TD=(NStates={self.td_nstates.value()})")
            
        # WFN
        if self.wfn_chk.isChecked():
            route_parts.append("Output=WFN")
            
        self.preview_str = " ".join(route_parts)
        self.preview_label.setText(self.preview_str)

    def get_route(self):
        return self.preview_str
    
    def is_wfn_enabled(self):
        # Helper for main dialog to know if we need to append filename
        return self.wfn_chk.isChecked()

    def parse_route(self, route):
        # Very basic parsing to try and populate fields from an existing string
        # This is hard to do perfectly, so we might just best-effort match.
        route = route.upper()
        
        # Check Job Type
        if "OPT" in route and "FREQ" in route: self.job_type.setCurrentIndex(0)
        elif "OPT" in route: self.job_type.setCurrentIndex(1)
        elif "FREQ" in route: self.job_type.setCurrentIndex(2)
        elif "SP" in route: self.job_type.setCurrentIndex(3)
        elif "SCAN" in route: self.job_type.setCurrentIndex(4)
        elif "IRC" in route: self.job_type.setCurrentIndex(5)
        elif "STABLE" in route: self.job_type.setCurrentIndex(6)
        elif "VOLUME" in route: self.job_type.setCurrentIndex(7)
        
        # Check Method (Just a few examples)
        if "WB97XD" in route: 
            self.method_type.setCurrentText("DFT")
            self.method_name.setCurrentText("WB97XD")
        elif "B3LYP" in route:
            self.method_type.setCurrentText("DFT")
            self.method_name.setCurrentText("B3LYP")
        elif "M062X" in route:
            self.method_type.setCurrentText("DFT")
            self.method_name.setCurrentText("M062X")
        elif "HF" in route:
            self.method_type.setCurrentText("Hartree-Fock")
            self.method_name.setCurrentText("HF")
            
        # Check Basis
        if "6-31G(D)" in route: self.basis_set.setCurrentText("6-31G(d)")
        elif "CC-PVTZ" in route: self.basis_set.setCurrentText("cc-pVTZ")
        elif "DEF2TZVP" in route: self.basis_set.setCurrentText("def2TZVP")
        elif "GEN" in route: self.basis_set.setCurrentText("Gen")
        
        # Solvation
        if "SCRF" in route:
            if "SMD" in route: self.solv_model.setCurrentText("SMD")
            elif "CPCM" in route: self.solv_model.setCurrentText("CPCM")
            else: self.solv_model.setCurrentText("PCM")
            
        # WFN
        if "OUTPUT=WFN" in route or "OUTPUT=WFX" in route:
            self.wfn_chk.setChecked(True)

class GaussianSetupDialog(QDialog):
    """
    Gaussianインプット作成のための高機能設定ダイアログ
    Link 0, Route, Title, Charge/Mult, および追記データに対応
    """
    def __init__(self, parent=None, mol=None):
        super().__init__(parent)
        self.setWindowTitle("Gaussian Input Setup")
        self.resize(600, 700)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.WindowMinMaxButtonsHint)
        self.setSizeGripEnabled(True)
        self.mol = mol
        self.force_wfn_line = False # Track if WFN line is needed
        self.ui_ready = False # [FIXED]
        self.setup_ui()
        self.load_presets_from_file()

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

        # --- 1. Link 0 Section (System Resources) ---
        link0_group = QGroupBox("Link 0 Commands (System Resources)")
        link0_layout = QFormLayout()

        # Memory
        mem_layout = QHBoxLayout()
        self.mem_spin = QSpinBox()
        self.mem_spin.setRange(1, 999999)
        self.mem_spin.setValue(4) # Default
        self.mem_spin.valueChanged.connect(self.update_preview)
        self.mem_unit = QComboBox()
        self.mem_unit.addItems(["GB", "MB", "MW"])
        self.mem_unit.currentIndexChanged.connect(self.update_preview)
        mem_layout.addWidget(self.mem_spin)
        mem_layout.addWidget(self.mem_unit)
        link0_layout.addRow("Memory (%mem):", mem_layout)

        # Processors
        self.nproc_spin = QSpinBox()
        self.nproc_spin.setRange(1, 128)
        self.nproc_spin.setValue(4) # Default
        self.nproc_spin.valueChanged.connect(self.update_preview)
        link0_layout.addRow("Processors (%nprocshared):", self.nproc_spin)
        
        # Checkpoint (Optional override)
        self.chk_edit = QLineEdit()
        self.chk_edit.setPlaceholderText("Auto-generated from filename if empty")
        self.chk_edit.textChanged.connect(self.update_preview)
        link0_layout.addRow("Checkpoint (%chk):", self.chk_edit)

        link0_group.setLayout(link0_layout)
        main_layout.addWidget(link0_group)

        # --- 2. Job Specification (Route & Title) ---
        job_group = QGroupBox("Job Specification")
        job_layout = QFormLayout()

        # Route section
        route_layout = QHBoxLayout()
        self.keywords_edit = QLineEdit("#P B3LYP/6-31G(d) Opt Freq")
        self.keywords_edit.textChanged.connect(self.update_preview)
        self.btn_route_builder = QPushButton("Route Builder...")
        self.btn_route_builder.clicked.connect(self.open_route_builder)
        route_layout.addWidget(self.keywords_edit)
        route_layout.addWidget(self.btn_route_builder)
        job_layout.addRow("Route Section (#):", route_layout)

        # Title
        self.title_edit = QLineEdit("Generated by MoleditPy Plugin")
        self.title_edit.textChanged.connect(self.update_preview)
        job_layout.addRow("Title:", self.title_edit)

        job_group.setLayout(job_layout)
        main_layout.addWidget(job_group)

        # --- 3. Molecular State ---
        mol_group = QGroupBox("Molecular Specification")
        mol_layout = QHBoxLayout()
        
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.valueChanged.connect(self.validate_charge_mult)
        
        # Save default palette for restoration
        self.default_palette = self.charge_spin.palette()

        # Multiplicity
        self.mult_spin = QSpinBox()
        self.mult_spin.setRange(1, 10)
        self.mult_spin.valueChanged.connect(self.validate_charge_mult)

        mol_layout.addWidget(QLabel("Charge:"))
        mol_layout.addWidget(self.charge_spin)
        mol_layout.addWidget(QLabel("Multiplicity:"))
        mol_layout.addWidget(self.mult_spin)

        mol_group.setLayout(mol_layout)
        main_layout.addWidget(mol_group)
        
        # Calculate initial values
        self.calc_initial_charge_mult()

        # --- 4. Additional Input (The "Everything" section) ---
        tail_group = QGroupBox("Additional Input (Post-Coordinates)")
        tail_layout = QVBoxLayout()
        
        help_label = QLabel("For ModRedundant, Basis Sets (Gen), or Link1.\nAppended after the molecule specification.")
        help_label.setStyleSheet("color: gray; font-size: 10px;")
        tail_layout.addWidget(help_label)

        self.tail_edit = QTextEdit()
        self.tail_edit.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.tail_edit.setFixedHeight(100)
        self.tail_edit.setPlaceholderText("Example:\nD 1 2 3 4 S 10 10.0\n\nOr basis set data for Gen/GenECP...")
        self.tail_edit.textChanged.connect(self.update_preview)
        tail_layout.addWidget(self.tail_edit)

        tail_group.setLayout(tail_layout)
        main_layout.addWidget(tail_group)

        # --- Preview ---
        preview_group = QGroupBox("Preview")
        preview_layout = QVBoxLayout()
        self.preview_text = QTextEdit()
        self.preview_text.setReadOnly(False) # Editable
        self.preview_text.setFont(QFont("Courier New", 10))
        self.preview_text.setMinimumHeight(200) # Adjusted height
        self.preview_text.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        preview_layout.addWidget(self.preview_text, 1)
        
        btn_refresh = QPushButton("Reset/Refresh Preview")
        btn_refresh.clicked.connect(self.update_preview)
        preview_layout.addWidget(btn_refresh)
        
        preview_group.setLayout(preview_layout)
        main_layout.addWidget(preview_group, 1)

        # --- Save/Preview Buttons ---
        btn_box = QHBoxLayout()
        
        self.save_btn = QPushButton("Save Input File...")
        self.save_btn.clicked.connect(self.save_file)
        self.save_btn.setStyleSheet("font-weight: bold; padding: 5px;")
        btn_box.addWidget(self.save_btn)
        
        main_layout.addLayout(btn_box)

        self.setLayout(main_layout)
        
        # Initial preview
        self.ui_ready = True # [FIXED]
        self.update_preview()

    def get_coords_string(self):
        """
        分子オブジェクトから原子座標文字列を生成する
        """
        if not self.mol:
            return ""
        
        lines = []
        try:
            # RDKitの場合の座標取得
            conf = self.mol.GetConformer()
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                atom = self.mol.GetAtomWithIdx(i)
                symbol = atom.GetSymbol()
                # Symbol X Y Z formatted
                lines.append(f"{symbol: <4} {pos.x: >12.6f} {pos.y: >12.6f} {pos.z: >12.6f}")
        except Exception as e:
            return f"Error extracting coordinates: {str(e)}"
            
        return "\n".join(lines)

    def generate_input_content(self, filename_hint="gaussian_job"):
        """Generates the full content of the input file as a string."""
        lines = []
        
        # --- Link 0 Section ---
        lines.append(f"%nprocshared={self.nproc_spin.value()}")
        lines.append(f"%mem={self.mem_spin.value()}{self.mem_unit.currentText()}")
        
        chk_name = self.chk_edit.text().strip()
        if not chk_name:
             # Use hint if empty
             chk_name = f"{filename_hint}.chk"
             
        if chk_name and not chk_name.endswith(".chk"): chk_name += ".chk"
        if chk_name: lines.append(f"%chk={chk_name}")

        # --- Route Section ---
        route_line = self.keywords_edit.text().strip()
        if not route_line.startswith("#"):
            route_line = "# " + route_line
        lines.append(f"{route_line}")
        lines.append("") # Blank line required
        
        # --- Title Section ---
        title_line = self.title_edit.text().strip()
        if not title_line:
            title_line = "Gaussian Job"
        lines.append(f"{title_line}")
        lines.append("") # Blank line required
        
        # --- Charge & Multiplicity ---
        lines.append(f"{self.charge_spin.value()} {self.mult_spin.value()}")
        
        # --- Molecule Specification ---
        coords_block = self.get_coords_string()
        if "Error" in coords_block:
            return f"# Error generating coordinates: {coords_block}"
        lines.append(coords_block)
        lines.append("") # Must end coord block with newline
        
        # --- Additional Input (Tail) ---
        tail_content = self.tail_edit.toPlainText()
        if tail_content.strip():
            # Ensure blank line before tail if needed (already added above)
            lines.append(tail_content)
            # Ensure tail ends with newline (managed by join)
        
        # WFN Output Line ?
        # If output=wfn is in route, we need a filename line at the end
        if "OUTPUT=WFN" in route_line.upper() or "OUTPUT=WFX" in route_line.upper():
            wfn_file = f"{filename_hint}.wfn"
            lines.append(wfn_file)
            lines.append("")
        
        # Final blank line if not already
        if lines[-1] != "": lines.append("")
        
        return "\n".join(lines)

    def update_preview(self):
        if not getattr(self, 'ui_ready', False): # [FIXED]
            return 
        # We don't know filename yet, use placeholder
        self.preview_text.setText(self.generate_input_content(filename_hint="[filename]"))

    def preview_file(self):
        # Legacy stub or redirect
        self.update_preview()

    def save_file(self):
        # ファイル保存ダイアログ
        default_name = "gaussian_job.gjf"
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Gaussian Input", default_name, "Gaussian Input (*.gjf *.com );;All Files (*)"
        )

        if file_path:
            try:
                # ファイル名に基づいてchk名を決定
                filename_base = os.path.splitext(os.path.basename(file_path))[0]
                
                # Check for %chk in content
                content = self.preview_text.toPlainText()
                
                # Smart update of %chk
                import re
                chk_line = f"%chk={filename_base}.chk"
                
                # Update CHK
                if re.search(r"^%chk=.*$", content, re.MULTILINE):
                     content = re.sub(r"^%chk=.*$", chk_line, content, flags=re.MULTILINE)
                else:
                     # Insert at top (after other link0 if possible, or just top)
                     content = f"{chk_line}\n{content}"
                     
                # Smart update of WFN filename?
                # If content ends with "job.wfn", replace it with "filename.wfn"?
                if "job.wfn" in content:
                    content = content.replace("job.wfn", f"{filename_base}.wfn")

                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(content)

                QMessageBox.information(self, "Success", f"File saved:\n{file_path}")
                # Do not close dialog automatically
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save file:\n{str(e)}")

    def open_route_builder(self):
        """Open the Route Builder dialog."""
        dialog = RouteBuilderDialog(self, self.keywords_edit.text())
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.keywords_edit.setText(dialog.get_route())
            # Check if WFN was enabled in builder
            if dialog.is_wfn_enabled():
                self.force_wfn_line = True
            else:
                self.force_wfn_line = False
            self.update_preview()

    def calc_initial_charge_mult(self):
        """
        Calculate default charge and multiplicity from the molecule.
        """
        if not self.mol:
            return

        try:
            # Logic taken from ORCA xyz2inp GUI
            
            # Charge
            # Use Chem.GetFormalCharge for consistency
            try:
                charge = Chem.GetFormalCharge(self.mol)
            except:
                # Fallback if Chem.GetFormalCharge fails
                charge = 0
                if hasattr(self.mol, "GetFormalCharge"):
                   charge = self.mol.GetFormalCharge()

            # Multiplicity
            # RDKit keeps track of radical electrons on atoms
            num_radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in self.mol.GetAtoms())
            mult = num_radical_electrons + 1
            
            # Note: The ORCA plugin logic is simple (radicals + 1).
            # It does not force even-electron systems to be doublets explicitly if no radicals are set.
            # But the user said "orca... works well", so we trust this primary logic.
            
            # Set values
            self.charge_spin.setValue(int(charge))
            self.mult_spin.setValue(int(mult))
            self.validate_charge_mult()
            
        except Exception as e:
            print(f"Error calculating charge/mult: {e}")

    def validate_charge_mult(self):
        """
        Validate if the charge/multiplicity combination is chemically consistent.
        Turns input red if invalid.
        """
        if not self.mol:
            return

        try:
            charge = self.charge_spin.value()
            mult = self.mult_spin.value()
            
            total_protons = 0
            for i in range(self.mol.GetNumAtoms()):
                atom = self.mol.GetAtomWithIdx(i)
                total_protons += atom.GetAtomicNum()
                
            total_electrons = total_protons - charge
            
            # Check consistency
            # If electrons are even, multiplicity must be odd (1, 3, 5)
            # If electrons are odd, multiplicity must be even (2, 4, 6)
            is_valid = False
            if total_electrons % 2 == 0:
                if mult % 2 != 0: is_valid = True
            else:
                if mult % 2 == 0: is_valid = True
                
            if is_valid:
                self.charge_spin.setPalette(self.default_palette)
                self.mult_spin.setPalette(self.default_palette)
            else:
                p = self.charge_spin.palette()
                p.setColor(QPalette.ColorRole.Base, QColor("#FFDDDD"))
                p.setColor(QPalette.ColorRole.Text, Qt.GlobalColor.red)
                self.charge_spin.setPalette(p)
                self.mult_spin.setPalette(p)
            
            self.update_preview()
            
        except Exception as e:
            print(f"Validation error: {e}")

    # --- Preset Management Methods ---

    def load_presets_from_file(self):
        """Load presets from JSON file and populate combo box."""
        self.presets_data = {}
        if os.path.exists(SETTINGS_FILE):
            try:
                with open(SETTINGS_FILE, 'r', encoding='utf-8') as f:
                    self.presets_data = json.load(f)
            except Exception as e:
                print(f"Error loading presets: {e}")
        
        # Ensure Default exists if empty
        if "Default" not in self.presets_data:
            self.presets_data["Default"] = {
                "nproc": 4, "mem_val": 4, "mem_unit": "GB", "chk": "",
                "route": "# B3LYP/6-31G(d) Opt Freq", 
                "tail": ""
            }
        
        self.update_preset_combo()

    def update_preset_combo(self):
        """Update the combobox items from self.presets_data."""
        current = self.preset_combo.currentText()
        self.preset_combo.blockSignals(True)
        self.preset_combo.clear()
        self.preset_combo.addItems(sorted(self.presets_data.keys()))
        
        # Restore selection or select Default
        if current in self.presets_data:
            self.preset_combo.setCurrentText(current)
        elif "Default" in self.presets_data:
            self.preset_combo.setCurrentText("Default")
        self.preset_combo.blockSignals(False)
        
        # Apply the current selection
        self.apply_selected_preset()

    def apply_selected_preset(self):
        """Apply the settings from the selected preset to the UI."""
        name = self.preset_combo.currentText()
        if name not in self.presets_data:
            return
            
        data = self.presets_data[name]
        
        # Apply values
        self.nproc_spin.setValue(data.get("nproc", 4))
        self.mem_spin.setValue(data.get("mem_val", 4))
        self.mem_unit.setCurrentText(data.get("mem_unit", "GB"))
        self.chk_edit.setText(data.get("chk", ""))
        self.keywords_edit.setText(data.get("route", ""))
        # DO NOT set Title, Charge, Mult (Molecule specific)
        self.tail_edit.setPlainText(data.get("tail", ""))
        self.update_preview()

    def get_current_ui_settings(self):
        """Return a dict of current UI settings."""
        return {
            "nproc": self.nproc_spin.value(),
            "mem_val": self.mem_spin.value(),
            "mem_unit": self.mem_unit.currentText(),
            "chk": self.chk_edit.text(),
            "route": self.keywords_edit.text(),
            # DO NOT save Title, Charge, Mult
            "tail": self.tail_edit.toPlainText()
        }

    def save_presets_to_file(self):
        """Save self.presets_data to JSON file."""
        try:
            with open(SETTINGS_FILE, 'w', encoding='utf-8') as f:
                json.dump(self.presets_data, f, indent=4)
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to save presets:\n{e}")

    def save_preset_dialog(self):
        """Show dialog to save current settings as a new preset."""
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:", text=self.preset_combo.currentText())
        if ok and name:
            name = name.strip()
            if not name:
                return
            
            # Confirm overwrite
            if name in self.presets_data:
                ret = QMessageBox.question(self, "Overwrite", f"Preset '{name}' already exists. Overwrite?",
                                           QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
                if ret != QMessageBox.StandardButton.Yes:
                    return

            # Save
            self.presets_data[name] = self.get_current_ui_settings()
            self.save_presets_to_file()
            self.update_preset_combo()
            self.preset_combo.setCurrentText(name) # Select the new one

    def delete_preset(self):
        """Delete the currently selected preset."""
        name = self.preset_combo.currentText()
        if name == "Default":
            QMessageBox.warning(self, "Warning", "Cannot delete 'Default' preset.")
            return

        ret = QMessageBox.question(self, "Delete Preset", f"Are you sure you want to delete '{name}'?",
                                   QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if ret == QMessageBox.StandardButton.Yes:
            if name in self.presets_data:
                del self.presets_data[name]
                self.save_presets_to_file()
                self.update_preset_combo()

def run(mw):
    mol = getattr(mw, 'current_mol', None)
    
    if not mol or mol.GetNumAtoms() == 0:
        QMessageBox.warning(mw, PLUGIN_NAME, "No molecule loaded or molecule is empty.")
        return

    dialog = GaussianSetupDialog(parent=mw, mol=mol)
    dialog.exec()

def initialize(context):
    def show_dialog():
        mw = context.get_main_window()
        run(mw)
    context.add_export_action("Gaussian Input...", show_dialog)