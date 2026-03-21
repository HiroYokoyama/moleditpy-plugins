from PyQt6 import QtCore
# -*- coding: utf-8 -*-
import os
from PyQt6.QtWidgets import (
    QCompleter, QMessageBox, QDialog, QVBoxLayout, QLabel, 
    QLineEdit, QSpinBox, QPushButton, QFileDialog, 
    QFormLayout, QGroupBox, QHBoxLayout, QComboBox, QTextEdit, 
    QTabWidget, QCheckBox, QWidget, QScrollArea, QMenu, QSizePolicy,
    QInputDialog)
from PyQt6.QtGui import QPalette, QColor, QAction, QFont, QSyntaxHighlighter, QTextCharFormat
from PyQt6.QtCore import Qt, QRegularExpression
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import json

PLUGIN_NAME = "ORCA Input Generator Neo"
PLUGIN_VERSION = "2026.03.21"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Advanced ORCA Input Generator with Preview and Presets"
SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "orca_input_generator_neo.json")

# ORCA Methods (Cleaned)
ALL_ORCA_METHODS = [
    # 密度汎関数理論 (DFT)
    "B2GP-PLYP", "B2PLYP", "B3LYP", "B3PW91", "B97-3c", "B97M-rV", "B97M-V", "BHandHLYP", "BLYP", "BP", 
    "BP86", "CAM-B3LYP", "DSD-BLYP", "DSD-PBEP86", "LC-wPBE", "M06", "M06-2X", "M06-HF", "M06-L", 
    "mPW2PLYP", "O3LYP", "PBE", "PBE0", "PBEh-3c", "PTPSS-D3", "PW91", "PWPB95", "r2SCAN-3c", "SCAN", 
    "TPSS", "TPSSh", "wB97", "wB97M-V", "wB97X", "wB97X-D3", "wB97X-V", "wPBE", "wPBEh", "X3LYP",
    # 波動関数理論 (HF / Post-HF)
    "ACPF", "AQCC", "BD", "CASSCF", "CCSD", "CCSD(T)", "CCSD(T)-F12", "CCSD(T)-F12/RI", "CCSD(T)-F12D/RI", 
    "CCSD-F12", "CCSD-F12/RI", "CCSD-F12D/RI", "CEPA/1", "CEPA/2", "CEPA/3", "DLPNO-CCSD", "DLPNO-CCSD(T)", 
    "DLPNO-CCSD(T)-F12", "DLPNO-CCSD(T)-F12/D", "DLPNO-CCSD(T1)", "DLPNO-CCSD(T1)-F12", "DLPNO-CCSD(T1)-F12/D", 
    "DLPNO-CCSD-F12", "DLPNO-CCSD-F12/D", "DLPNO-MP2", "DLPNO-NEVPT2", "DLPNO-SCS-MP2", "DLPNO-SOS-MP2", 
    "F12-DLPNO-MP2", "F12-MP2", "F12-RI-MP2", "F12/D-DLPNO-MP2", "F12/D-RI-MP2", "HF", "MP2", "MP3", "MRCI", 
    "NCEPA/1", "NCPF/1", "NEVPT2", "OO-RI-MP2", "OO-RI-SCS-MP2", "OO-RI-SOS-MP2", "QCISD", "QCISD(T)", 
    "QCISD(T)-F12", "QCISD(T)-F12/RI", "QCISD-F12", "QCISD-F12/RI", "RI-CEPA/1-F12", "RI-MP2", "RI-SCS-MP2", 
    "RI-SOS-MP2", "ROHF", "SCS-MP2", "SCS-MP3", "UHF",
    # 半経験的手法・xTB・力場
    "AM1", "GFN-FF", "GFN0-XTB", "GFN1-XTB", "GFN1-xTB", "GFN2-XTB", "GFN2-xTB", "MNDO", "NATIVE-GFN-XTB", 
    "NATIVE-GFN1-XTB", "NATIVE-GFN2-XTB", "NATIVE-spGFN-XTB", "NATIVE-spGFN1-XTB", "NATIVE-spGFN2-XTB", 
    "NATIVE-spXTB", "NATIVE-spXTB1", "NATIVE-spXTB2", "NATIVE-XTB", "NATIVE-XTB1", "NATIVE-XTB2", "PM3", 
    "PM6", "XTB", "XTB0", "XTB1", "XTB2", "XTBFF", "ZINDO/1", "ZINDO/2", "ZINDO/S", "ZNDDO/1", "ZNDDO/2"
]

# ORCA Basis Sets (Cleaned)
ALL_ORCA_BASIS_SETS = [
    "3-21G", "3-21GSP", "4-22GSP", "6-31++G(2d2p)", "6-31++G(2df2p)", "6-31++G(2df2pd)", "6-31++G(2dp)", 
    "6-31++G(dp)", "6-31++G**", "6-31+G(2d)", "6-31+G(2d2p)", "6-31+G(2df)", "6-31+G(2df2p)", "6-31+G(2df2pd)", 
    "6-31+G(2dp)", "6-31+G(d)", "6-31+G(dp)", "6-31+G*", "6-31+G**", "6-311++G(2d2p)", "6-311++G(2df2p)", 
    "6-311++G(2df2pd)", "6-311++G(2dp)", "6-311++G(3df3pd)", "6-311++G(dp)", "6-311++G**", "6-311+G(2d)", 
    "6-311+G(2d2p)", "6-311+G(2df)", "6-311+G(2df2p)", "6-311+G(2df2pd)", "6-311+G(2dp)", "6-311+G(3df)", 
    "6-311+G(3df2p)", "6-311+G(3df3pd)", "6-311+G(d)", "6-311+G(dp)", "6-311+G*", "6-311+G**", "6-311G", 
    "6-311G(2d)", "6-311G(2d2p)", "6-311G(2df)", "6-311G(2df2p)", "6-311G(2df2pd)", "6-311G(2dp)", "6-311G(3df)", 
    "6-311G(3df3pd)", "6-311G(d)", "6-311G(dp)", "6-311G*", "6-311G**", "6-31G", "6-31G(2d)", "6-31G(2d2p)", 
    "6-31G(2df)", "6-31G(2df2p)", "6-31G(2df2pd)", "6-31G(2dp)", "6-31G(d)", "6-31G(dp)", "6-31G*", "6-31G**", 
    "AHGBS-5", "AHGBS-7", "AHGBS-9", "AHGBSP1-5", "AHGBSP1-7", "AHGBSP1-9", "AHGBSP2-5", "AHGBSP2-7", 
    "AHGBSP2-9", "AHGBSP3-5", "AHGBSP3-7", "AHGBSP3-9", "ANO-pV5Z", "ANO-pV6Z", "ANO-pVDZ", "ANO-pVQZ", 
    "ANO-pVTZ", "ANO-RCC-DZP", "ANO-RCC-Full", "ANO-RCC-QZP", "ANO-RCC-TZP", "ANO-SZ", "apr-cc-pV(Q+d)Z", 
    "aug-ANO-pV5Z", "aug-ANO-pVDZ", "aug-ANO-pVQZ", "aug-ANO-pVTZ", "aug-cc-pCV5Z", "aug-cc-pCV5Z-PP", 
    "aug-cc-pCV6Z", "aug-cc-pCVDZ", "aug-cc-pCVDZ-PP", "aug-cc-pCVQZ", "aug-cc-pCVQZ-PP", "aug-cc-pCVTZ", 
    "aug-cc-pCVTZ-PP", "aug-cc-pV5(+d)Z", "aug-cc-pV5Z", "aug-cc-pV5Z-DK", "aug-cc-pV5Z-PP", "aug-cc-pV6(+d)Z", 
    "aug-cc-pV6Z", "aug-cc-pVD(+d)Z", "aug-cc-pVDZ", "aug-cc-pVDZ-DK", "aug-cc-pVDZ-PP", "aug-cc-pVQ(+d)Z", 
    "aug-cc-pVQZ", "aug-cc-pVQZ-DK", "aug-cc-pVQZ-PP", "aug-cc-pVT(+d)Z", "aug-cc-pVTZ", "aug-cc-pVTZ-DK", 
    "aug-cc-pVTZ-J", "aug-cc-pVTZ-PP", "aug-cc-pwCV5Z", "aug-cc-pwCV5Z-DK", "aug-cc-pwCV5Z-PP", "aug-cc-pwCVDZ", 
    "aug-cc-pwCVDZ-DK", "aug-cc-pwCVDZ-PP", "aug-cc-pwCVQZ", "aug-cc-pwCVQZ-DK", "aug-cc-pwCVQZ-PP", 
    "aug-cc-pwCVTZ", "aug-cc-pwCVTZ-DK", "aug-cc-pwCVTZ-PP", "aug-pc-0", "aug-pc-1", "aug-pc-2", "aug-pc-3", 
    "aug-pc-4", "aug-pcH-1", "aug-pcH-2", "aug-pcH-3", "aug-pcH-4", "aug-pcJ-0", "aug-pcJ-1", "aug-pcJ-2", 
    "aug-pcJ-3", "aug-pcJ-4", "aug-pcseg-0", "aug-pcseg-1", "aug-pcseg-2", "aug-pcseg-3", "aug-pcseg-4", 
    "aug-pcSseg-0", "aug-pcSseg-1", "aug-pcSseg-2", "aug-pcSseg-3", "aug-pcSseg-4", "aug-pcX-1", "aug-pcX-2", 
    "aug-pcX-3", "aug-pcX-4", "cc-pCV5Z", "cc-pCV5Z-PP", "cc-pCV6Z", "cc-pCVDZ", "cc-pCVDZ-F12", "cc-pCVDZ-PP", 
    "cc-pCVQZ", "cc-pCVQZ-F12", "cc-pCVQZ-PP", "cc-pCVTZ", "cc-pCVTZ-F12", "cc-pCVTZ-PP", "cc-pV5(+d)Z", 
    "cc-pV5Z", "cc-pV5Z-DK", "cc-pV5Z-PP", "cc-pV6Z", "cc-pVD(+d)Z", "cc-pVDZ", "cc-pVDZ-DK", "cc-pVDZ-DK3", 
    "cc-pVDZ-F12", "cc-pVDZ-PP", "cc-pVDZ-PP-F12", "cc-pVQ(+d)Z", "cc-pVQZ", "cc-pVQZ-DK", "cc-pVQZ-DK3", 
    "cc-pVQZ-F12", "cc-pVQZ-PP", "cc-pVQZ-PP-F12", "cc-pVT(+d)Z", "cc-pVTZ", "cc-pVTZ-DK", "cc-pVTZ-DK3", 
    "cc-pVTZ-F12", "cc-pVTZ-PP", "cc-pVTZ-PP-F12", "cc-pwCV5Z", "cc-pwCV5Z-DK", "cc-pwCV5Z-PP", "cc-pwCVDZ", 
    "cc-pwCVDZ-DK", "cc-pwCVDZ-DK3", "cc-pwCVDZ-PP", "cc-pwCVQZ", "cc-pwCVQZ-DK", "cc-pwCVQZ-DK3", 
    "cc-pwCVQZ-PP", "cc-pwCVTZ", "cc-pwCVTZ-DK", "cc-pwCVTZ-DK3", "cc-pwCVTZ-PP", "CRENBL", "D95", "D95p", 
    "def-SV(P)", "def-SVP", "def-TZVP", "def-TZVPP", "def2-mSVP", "def2-mTZVP", "def2-mTZVPP", "def2-QZVP", 
    "def2-QZVPD", "def2-QZVPP", "def2-QZVPPD", "def2-SV(P)", "def2-SVP", "def2-SVPD", "def2-TZVP", 
    "def2-TZVP(-f)", "def2-TZVPD", "def2-TZVPP", "def2-TZVPPD", "dhf-QZVP", "dhf-QZVP-2c", "dhf-QZVPP", 
    "dhf-QZVPP-2c", "dhf-SV(P)", "dhf-SVP", "dhf-SVP-2c", "dhf-TZVP", "dhf-TZVP-2c", "dhf-TZVPP", 
    "dhf-TZVPP-2c", "DKH-def2-QZVPP", "DKH-def2-SV(P)", "DKH-def2-SVP", "DKH-def2-TZVP", "DKH-def2-TZVP(-f)", 
    "DKH-def2-TZVPP", "DKH-QZVP", "DKH-QZVPP", "DKH-SV(P)", "DKH-SVP", "DKH-TZV(P)", "DKH-TZVP", "DKH-TZVPP", 
    "EPR-II", "EPR-III", "haV(5+d)Z", "haV(Q+d)Z", "haV(T+d)Z", "HGBS-5", "HGBS-7", "HGBS-9", "HGBSP1-5", 
    "HGBSP1-7", "HGBSP1-9", "HGBSP2-5", "HGBSP2-7", "HGBSP2-9", "HGBSP3-5", "HGBSP3-7", "HGBSP3-9", "IGLO-II", 
    "IGLO-III", "jul-cc-pV(D+d)Z", "jul-cc-pV(Q+d)Z", "jul-cc-pV(T+d)Z", "jun-cc-pV(D+d)Z", "jun-cc-pV(Q+d)Z", 
    "jun-cc-pV(T+d)Z", "LANL08", "LANL08(f)", "LANL2DZ", "LANL2TZ", "LANL2TZ(f)", "m6-31G", "m6-31G*", 
    "ma-def-TZVP", "ma-def2-mSVP", "ma-def2-QZVP", "ma-def2-QZVPP", "ma-def2-SV(P)", "ma-def2-SVP", "ma-def2-TZVP", 
    "ma-def2-TZVP(-f)", "ma-def2-TZVPP", "ma-DKH-def2-QZVPP", "ma-DKH-def2-SV(P)", "ma-DKH-def2-SVP", 
    "ma-DKH-def2-TZVP", "ma-DKH-def2-TZVP(-f)", "ma-DKH-def2-TZVPP", "ma-ZORA-def2-QZVPP", "ma-ZORA-def2-SV(P)", 
    "ma-ZORA-def2-SVP", "ma-ZORA-def2-TZVP", "ma-ZORA-def2-TZVP(-f)", "ma-ZORA-def2-TZVPP", "maug-cc-pV(D+d)Z", 
    "maug-cc-pV(Q+d)Z", "maug-cc-pV(T+d)Z", "may-cc-pV(Q+d)Z", "may-cc-pV(T+d)Z", "MIDI", "MINI", "MINIS", 
    "MINIX", "old-DKH-SV(P)", "old-DKH-SVP", "old-DKH-TZV(P)", "old-DKH-TZVP", "old-DKH-TZVPP", "old-SV", 
    "old-SV(P)", "old-SVP", "old-TZV", "old-TZV(P)", "old-TZVP", "old-TZVPP", "old-ZORA-SV(P)", "old-ZORA-SVP", 
    "old-ZORA-TZV(P)", "old-ZORA-TZVP", "old-ZORA-TZVPP", "Partridge-1", "Partridge-2", "Partridge-3", 
    "Partridge-4", "pc-0", "pc-1", "pc-2", "pc-3", "pc-4", "pcH-1", "pcH-2", "pcH-3", "pcH-4", "pcJ-0", 
    "pcJ-1", "pcJ-2", "pcJ-3", "pcJ-4", "pcseg-0", "pcseg-1", "pcseg-2", "pcseg-3", "pcseg-4", "pcSseg-0", 
    "pcSseg-1", "pcSseg-2", "pcSseg-3", "pcSseg-4", "pcX-1", "pcX-2", "pcX-3", "pcX-4", "QZVP", "QZVPP", 
    "Sapporo-DKH3-DZP-2012", "Sapporo-DKH3-QZP-2012", "Sapporo-DKH3-TZP-2012", "Sapporo-DZP-2012", 
    "Sapporo-QZP-2012", "Sapporo-TZP-2012", "SARC-DKH-SVP", "SARC-DKH-TZVP", "SARC-DKH-TZVPP", "SARC-ZORA-SVP", 
    "SARC-ZORA-TZVP", "SARC-ZORA-TZVPP", "SARC2-DKH-QZV", "SARC2-DKH-QZVP", "SARC2-ZORA-QZV", "SARC2-ZORA-QZVP", 
    "saug-ANO-pV5Z", "saug-ANO-pVDZ", "saug-ANO-pVQZ", "saug-ANO-pVTZ", "STO-3G", "SV", "SV(P)", "SVP", "TZV", 
    "TZV(P)", "TZVP", "TZVPP", "UGBS", "vDZP", "W1-DZ", "W1-mtsmall", "W1-Opt", "W1-QZ", "W1-TZ", "Wachters+f", 
    "x2c-QZVPall", "x2c-QZVPall-2c", "x2c-QZVPall-2c-s", "x2c-QZVPall-s", "x2c-QZVPPall", "x2c-QZVPPall-2c", 
    "x2c-QZVPPall-2c-s", "x2c-QZVPPall-s", "x2c-SV(P)all", "x2c-SV(P)all-2c", "x2c-SV(P)all-s", "x2c-SVPall", 
    "x2c-SVPall-2c", "x2c-SVPall-s", "x2c-TZVPall", "x2c-TZVPall-2c", "x2c-TZVPall-s", "x2c-TZVPPall", 
    "x2c-TZVPPall-2c", "x2c-TZVPPall-s", "ZORA-def2-QZVPP", "ZORA-def2-SV(P)", "ZORA-def2-SVP", "ZORA-def2-TZVP", 
    "ZORA-def2-TZVP(-f)", "ZORA-def2-TZVPP", "ZORA-QZVP", "ZORA-QZVPP", "ZORA-SV(P)", "ZORA-SVP", "ZORA-TZV(P)", 
    "ZORA-TZVP", "ZORA-TZVPP"
]

class OrcaSyntaxHighlighter(QSyntaxHighlighter):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.rules = []

        # Keywords (!)
        keyword_format = QTextCharFormat()
        keyword_format.setForeground(QColor("#D32F2F")) # Red
        keyword_format.setFontWeight(QFont.Weight.Bold)
        self.rules.append((QRegularExpression("^!.*"), keyword_format))

        # Blocks (%)
        block_format = QTextCharFormat()
        block_format.setForeground(QColor("#1976D2")) # Blue
        block_format.setFontWeight(QFont.Weight.Bold)
        self.rules.append((QRegularExpression("^%.*"), block_format))
        self.rules.append((QRegularExpression("^end.*"), block_format))

        # Coordinates/Title (*)
        coord_format = QTextCharFormat()
        coord_format.setForeground(QColor("#388E3C")) # Green
        coord_format.setFontWeight(QFont.Weight.Bold)
        self.rules.append((QRegularExpression("^\\*.*"), coord_format))

        # Comments (#)
        comment_format = QTextCharFormat()
        comment_format.setForeground(QColor("#757575")) # Grey
        self.rules.append((QRegularExpression("#.*"), comment_format))

    def highlightBlock(self, text):
        for pattern, format in self.rules:
            match_iterator = pattern.globalMatch(text)
            while match_iterator.hasNext():
                match = match_iterator.next()
                self.setFormat(match.capturedStart(), match.capturedLength(), format)

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
        
        # --- Tab 5: TD-DFT ---
        self.tab_tddft = QWidget()
        self.setup_tddft_tab()
        self.tabs.addTab(self.tab_tddft, "TD-DFT")

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

    def get_inferred_category(self, text):
        if not text: return "All Methods"
        text_upper = text.upper()
        
        if text_upper in [m.upper() for m in ['CAM-B3LYP', 'LC-wPBE', 'wB97', 'wB97M-V', 'wB97X', 'wB97X-D3', 'wB97X-V', 'wPBE', 'wPBEh']]:
            return "DFT (Range-Separated)"
        elif text_upper in [m.upper() for m in ['B2GP-PLYP', 'B2PLYP', 'DSD-BLYP', 'DSD-PBEP86', 'mPW2PLYP', 'PTPSS-D3', 'PWPB95']]:
            return "DFT (Double Hybrid)"
        elif text_upper in [m.upper() for m in ['B3LYP', 'B3PW91', 'B97-3c', 'B97M-rV', 'B97M-V', 'BHandHLYP', 'BLYP', 'BP', 'BP86', 'M06', 'M06-2X', 'M06-HF', 'M06-L', 'O3LYP', 'PBE', 'PBE0', 'PBEh-3c', 'PW91', 'r2SCAN-3c', 'SCAN', 'TPSS', 'TPSSh', 'X3LYP']]:
            return "DFT (GGA/Hybrid/Meta)"
        elif text_upper in [m.upper() for m in ['HF', 'MP2', 'MP3', 'DLPNO-MP2', 'DLPNO-SCS-MP2', 'DLPNO-SOS-MP2', 'F12-DLPNO-MP2', 'F12-MP2', 'F12-RI-MP2', 'F12/D-DLPNO-MP2', 'F12/D-RI-MP2', 'OO-RI-MP2', 'OO-RI-SCS-MP2', 'OO-RI-SOS-MP2', 'RI-MP2', 'RI-SCS-MP2', 'RI-SOS-MP2', 'ROHF', 'SCS-MP2', 'SCS-MP3', 'UHF']]:
            return "Wavefunction (HF/MP2)"
        elif text_upper in [m.upper() for m in ['ACPF', 'AQCC', 'BD', 'CCSD', 'CCSD(T)', 'CCSD(T)-F12', 'CCSD(T)-F12/RI', 'CCSD(T)-F12D/RI', 'CCSD-F12', 'CCSD-F12/RI', 'CCSD-F12D/RI', 'CEPA/1', 'CEPA/2', 'CEPA/3', 'DLPNO-CCSD', 'DLPNO-CCSD(T)', 'DLPNO-CCSD(T)-F12', 'DLPNO-CCSD(T)-F12/D', 'DLPNO-CCSD(T1)', 'DLPNO-CCSD(T1)-F12', 'DLPNO-CCSD(T1)-F12/D', 'DLPNO-CCSD-F12', 'DLPNO-CCSD-F12/D', 'NCEPA/1', 'NCPF/1', 'QCISD', 'QCISD(T)', 'QCISD(T)-F12', 'QCISD(T)-F12/RI', 'QCISD-F12', 'QCISD-F12/RI', 'RI-CEPA/1-F12']]:
            return "Wavefunction (Coupled Cluster)"
        elif text_upper in [m.upper() for m in ['CASSCF', 'NEVPT2', 'DLPNO-NEVPT2', 'MRCI']]:
            return "Wavefunction (Multireference)"
        elif text_upper in [m.upper() for m in ['AM1', 'GFN-FF', 'GFN0-XTB', 'GFN1-XTB', 'GFN1-xTB', 'GFN2-XTB', 'GFN2-xTB', 'MNDO', 'NATIVE-GFN-XTB', 'NATIVE-GFN1-XTB', 'NATIVE-GFN2-XTB', 'NATIVE-spGFN-XTB', 'NATIVE-spGFN1-XTB', 'NATIVE-spGFN2-XTB', 'NATIVE-spXTB', 'NATIVE-spXTB1', 'NATIVE-spXTB2', 'NATIVE-XTB', 'NATIVE-XTB1', 'NATIVE-XTB2', 'PM3', 'PM6', 'XTB', 'XTB0', 'XTB1', 'XTB2', 'XTBFF', 'ZINDO/1', 'ZINDO/2', 'ZINDO/S', 'ZNDDO/1', 'ZNDDO/2']]:
            return "Semi-Empirical"
            
        return "All Methods"

    def setup_method_tab(self):
        layout = QFormLayout()

        self.method_type = QComboBox()
        self.method_type.addItems([
            "DFT (GGA/Hybrid/Meta)", 
            "DFT (Range-Separated)", 
            "DFT (Double Hybrid)",
            "Wavefunction (HF/MP2)",
            "Wavefunction (Coupled Cluster)",
            "Wavefunction (Multireference)",
            "Semi-Empirical",
            "All Methods"
        ])
        self.method_type.currentIndexChanged.connect(self.update_method_list)
        layout.addRow("Method Type:", self.method_type)
        
        self.method_name = QComboBox()
        self.method_name.setEditable(True)
        m_completer = QCompleter(ALL_ORCA_METHODS, self)
        m_completer.setCaseSensitivity(QtCore.Qt.CaseSensitivity.CaseInsensitive)
        m_completer.setFilterMode(QtCore.Qt.MatchFlag.MatchContains)
        self.method_name.setCompleter(m_completer)
        self.update_method_list()
        layout.addRow("Method:", self.method_name)
        
        self.basis_set = QComboBox()
        self.basis_set.setEditable(True)
        self.basis_set.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        basis_groups = [
            "--- Karlsruhe (Def2) ---",
            "def2-SV(P)", "def2-SVP", "def2-TZVP", "def2-QZVP",
            "def2-SVPD", "def2-TZVPD", "def2-QZVPD",
            "ma-def2-SVP", "ma-def2-TZVP", "ma-def2-QZVP",
            "def2-TZVPP", "def2-QZVPP", "def2-TZVPPD", "def2-QVPPD",
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
        b_completer = QCompleter(ALL_ORCA_BASIS_SETS, self)
        b_completer.setCaseSensitivity(QtCore.Qt.CaseSensitivity.CaseInsensitive)
        b_completer.setFilterMode(QtCore.Qt.MatchFlag.MatchContains)
        self.basis_set.setCompleter(b_completer)
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
        current_text = self.method_name.currentText()
        
        self.method_name.blockSignals(True)
        self.method_name.clear()
        
        if "GGA/Hybrid" in mtype:
             self.method_name.addItems(['B3LYP', 'B3PW91', 'B97-3c', 'B97M-rV', 'B97M-V', 'BHandHLYP', 'BLYP', 'BP', 'BP86', 'M06', 'M06-2X', 'M06-HF', 'M06-L', 'O3LYP', 'PBE', 'PBE0', 'PBEh-3c', 'PW91', 'r2SCAN-3c', 'SCAN', 'TPSS', 'TPSSh', 'X3LYP'])
        elif "Range-Separated" in mtype:
             self.method_name.addItems(['CAM-B3LYP', 'LC-wPBE', 'wB97', 'wB97M-V', 'wB97X', 'wB97X-D3', 'wB97X-V', 'wPBE', 'wPBEh'])
        elif "Double Hybrid" in mtype:
             self.method_name.addItems(['B2GP-PLYP', 'B2PLYP', 'DSD-BLYP', 'DSD-PBEP86', 'mPW2PLYP', 'PTPSS-D3', 'PWPB95'])
        elif "HF/MP2" in mtype:
             self.method_name.addItems(['HF', 'MP2', 'MP3', 'DLPNO-MP2', 'DLPNO-SCS-MP2', 'DLPNO-SOS-MP2', 'F12-DLPNO-MP2', 'F12-MP2', 'F12-RI-MP2', 'F12/D-DLPNO-MP2', 'F12/D-RI-MP2', 'OO-RI-MP2', 'OO-RI-SCS-MP2', 'OO-RI-SOS-MP2', 'RI-MP2', 'RI-SCS-MP2', 'RI-SOS-MP2', 'ROHF', 'SCS-MP2', 'SCS-MP3', 'UHF'])
        elif "Coupled Cluster" in mtype:
             self.method_name.addItems(['ACPF', 'AQCC', 'BD', 'CCSD', 'CCSD(T)', 'CCSD(T)-F12', 'CCSD(T)-F12/RI', 'CCSD(T)-F12D/RI', 'CCSD-F12', 'CCSD-F12/RI', 'CCSD-F12D/RI', 'CEPA/1', 'CEPA/2', 'CEPA/3', 'DLPNO-CCSD', 'DLPNO-CCSD(T)', 'DLPNO-CCSD(T)-F12', 'DLPNO-CCSD(T)-F12/D', 'DLPNO-CCSD(T1)', 'DLPNO-CCSD(T1)-F12', 'DLPNO-CCSD(T1)-F12/D', 'DLPNO-CCSD-F12', 'DLPNO-CCSD-F12/D', 'NCEPA/1', 'NCPF/1', 'QCISD', 'QCISD(T)', 'QCISD(T)-F12', 'QCISD(T)-F12/RI', 'QCISD-F12', 'QCISD-F12/RI', 'RI-CEPA/1-F12'])
        elif "Multireference" in mtype:
             self.method_name.addItems(['CASSCF', 'NEVPT2', 'DLPNO-NEVPT2', 'MRCI'])
        elif "Semi-Empirical" in mtype:
             self.method_name.addItems(['AM1', 'GFN-FF', 'GFN0-XTB', 'GFN1-XTB', 'GFN1-xTB', 'GFN2-XTB', 'GFN2-xTB', 'MNDO', 'NATIVE-GFN-XTB', 'NATIVE-GFN1-XTB', 'NATIVE-GFN2-XTB', 'NATIVE-spGFN-XTB', 'NATIVE-spGFN1-XTB', 'NATIVE-spGFN2-XTB', 'NATIVE-spXTB', 'NATIVE-spXTB1', 'NATIVE-spXTB2', 'NATIVE-XTB', 'NATIVE-XTB1', 'NATIVE-XTB2', 'PM3', 'PM6', 'XTB', 'XTB0', 'XTB1', 'XTB2', 'XTBFF', 'ZINDO/1', 'ZINDO/2', 'ZINDO/S', 'ZNDDO/1', 'ZNDDO/2'])
        elif mtype == "All Methods":
             self.method_name.addItems(ALL_ORCA_METHODS)
             
        if mtype == "All Methods":
            if current_text:
                self.method_name.setCurrentText(current_text)
        else:
            if self.method_name.count() > 0:
                 self.method_name.setCurrentIndex(0)
             
        self.method_name.blockSignals(False)
        self.update_ui_state()
        self.update_preview()

    def setup_job_tab(self):
        layout = QVBoxLayout()
        
        self.job_type = QComboBox()
        self.job_type.addItems([
            "Optimization + Freq (Opt Freq)", 
            "Optimization Only (Opt)", 
            "Frequency Only (Freq)", 
            "Single Point Energy (SP)",
            "NMR",
            "Scan (Relaxed Surface)",
            "Transition State Opt (OptTS)",
            "Gradient",
            "Hessian"
        ])
        layout.addWidget(QLabel("Job Task:"))
        layout.addWidget(self.job_type)
        self.job_type.currentIndexChanged.connect(self.update_ui_state)
        
        # SCF Options (Moved to below Task)
        self.scf_group = QGroupBox("SCF Convergence")
        scf_layout = QGridLayout()
        
        self.scf_sloppy = QCheckBox("Sloppy")
        self.scf_loose = QCheckBox("Loose")
        self.scf_normal = QCheckBox("Normal")
        self.scf_strong = QCheckBox("Strong")
        self.scf_tight = QCheckBox("Tight")
        self.scf_verytight = QCheckBox("VeryTight")
        self.scf_extreme = QCheckBox("Extreme")
        
        # Row 1
        scf_layout.addWidget(self.scf_sloppy, 0, 0)
        scf_layout.addWidget(self.scf_loose, 0, 1)
        scf_layout.addWidget(self.scf_normal, 0, 2)
        scf_layout.addWidget(self.scf_strong, 0, 3)
        # Row 2
        scf_layout.addWidget(self.scf_tight, 1, 0)
        scf_layout.addWidget(self.scf_verytight, 1, 1)
        scf_layout.addWidget(self.scf_extreme, 1, 2)
        
        self.scf_group.setLayout(scf_layout)
        layout.addWidget(self.scf_group)
        
        # Opt Options
        self.opt_group = QGroupBox("Optimization Options")
        opt_layout = QHBoxLayout()
        self.opt_tight = QCheckBox("TightOpt")
        self.opt_verytight = QCheckBox("VeryTightOpt")
        self.opt_loose = QCheckBox("LooseOpt")
        self.opt_cart = QCheckBox("COpt (Cartesian)")
        self.opt_calcfc = QCheckBox("CalcFC")
        self.opt_ts_mode = QCheckBox("CalcHess (for TS)")
        
        opt_layout.addWidget(self.opt_tight)
        opt_layout.addWidget(self.opt_verytight)
        opt_layout.addWidget(self.opt_loose)
        opt_layout.addWidget(self.opt_cart)
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
        self.grid_combo.addItems(["Default", "defgrid1", "defgrid2", "defgrid3", "Grid4", "Grid5", "Grid6", "NoGrid"])
        self.grid_combo.setCurrentText("Default")
        layout.addRow("Grid:", self.grid_combo)


        # NBO
        self.pop_nbo = QCheckBox("NBO Analysis (! NBO)")
        layout.addRow(self.pop_nbo)

        self.tab_props.setLayout(layout)

    def setup_tddft_tab(self):
        layout = QFormLayout()
        
        self.tddft_enable = QCheckBox("Enable TD-DFT (%)")
        self.tddft_enable.toggled.connect(self.update_preview)
        layout.addRow(self.tddft_enable)
        
        self.tddft_nroots = QSpinBox()
        self.tddft_nroots.setRange(1, 500)
        self.tddft_nroots.setValue(10)
        layout.addRow("Number of Roots (NRoots):", self.tddft_nroots)
        
        self.tddft_singlets = QCheckBox("Singlets")
        self.tddft_singlets.setChecked(True)
        self.tddft_triplets = QCheckBox("Triplets")
        self.tddft_triplets.setChecked(False)
        
        s_layout = QHBoxLayout()
        s_layout.addWidget(self.tddft_singlets)
        s_layout.addWidget(self.tddft_triplets)
        layout.addRow("States:", s_layout)
        
        self.tddft_tda = QCheckBox("Use TDA (Tamm-Dancoff Approximation)")
        self.tddft_tda.setChecked(True)
        layout.addRow(self.tddft_tda)
        
        self.tddft_iroot = QSpinBox()
        self.tddft_iroot.setRange(1, 500)
        self.tddft_iroot.setValue(1)
        layout.addRow("Root for Polar./Grad. (IRoot):", self.tddft_iroot)

        self.tab_tddft.setLayout(layout)

    def connect_signals(self):
        widgets = [
            self.method_type, self.method_name, self.basis_set, self.aux_basis,
            self.job_type, self.opt_tight, self.opt_verytight, self.opt_loose, 
            self.opt_cart, self.opt_calcfc, self.opt_ts_mode,
            self.freq_num, self.freq_raman,
            self.solv_model, self.solvent, self.dispersion,
            self.solv_model, self.solvent, self.dispersion,
            self.rijcosx, self.grid_combo, 
            self.scf_sloppy, self.scf_loose, self.scf_normal, self.scf_strong, 
            self.scf_tight, self.scf_verytight, self.scf_extreme,
            self.pop_nbo, self.tddft_enable, self.tddft_nroots, self.tddft_singlets, 
            self.tddft_triplets, self.tddft_tda, self.tddft_iroot
        ]
        for w in widgets:
            if isinstance(w, QComboBox):
                w.currentIndexChanged.connect(self.update_preview)
                if w.isEditable():
                    w.currentTextChanged.connect(self.update_preview)
            elif isinstance(w, QCheckBox):
                w.toggled.connect(self.update_preview)
                # Mutual exclusivity for SCF
                if w in [self.scf_sloppy, self.scf_loose, self.scf_normal, self.scf_strong, 
                           self.scf_tight, self.scf_verytight, self.scf_extreme]:
                    w.clicked.connect(self.enforce_scf_mutual_exclusion)
                # Mutual exclusivity for Opt
                if w in [self.opt_tight, self.opt_verytight, self.opt_loose]:
                    w.clicked.connect(self.enforce_opt_mutual_exclusion)
            elif isinstance(w, QSpinBox):
                w.valueChanged.connect(self.update_preview)

    def update_ui_state(self):
        """Update usability of widgets based on current selection."""
        if not getattr(self, 'ui_ready', False): return

        # 1. Method Dependent
        method_text = self.method_name.currentText()
        mtype = self.get_inferred_category(method_text)
        is_semi = "Semi-Empirical" in mtype
        is_3c = "3C" in method_text.upper()
        no_basis = is_semi or is_3c
        
        # Disable Basis Set & Aux Basis for Semi-Empirical and 3c
        self.basis_set.setEnabled(not no_basis)
        self.aux_basis.setEnabled(not no_basis)
        
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

    def enforce_scf_mutual_exclusion(self):
        ctx = self.sender()
        if not ctx.isChecked(): return
        scf_boxes = [
            self.scf_sloppy, self.scf_loose, self.scf_normal, self.scf_strong, 
            self.scf_tight, self.scf_verytight, self.scf_extreme
        ]
        for cb in scf_boxes:
            if cb != ctx:
                cb.blockSignals(True)
                cb.setChecked(False)
                cb.blockSignals(False)
        self.update_preview()

    def enforce_opt_mutual_exclusion(self):
        ctx = self.sender()
        if not ctx.isChecked(): return
        for cb in [self.opt_tight, self.opt_verytight, self.opt_loose]:
            if cb != ctx:
                cb.blockSignals(True)
                cb.setChecked(False)
                cb.blockSignals(False)
        self.update_preview()

    def update_preview(self):
        if not getattr(self, 'ui_ready', False):
            return
        self.update_ui_state()

        route_parts = ["!"]
        
        # Method / Basis
        method = self.method_name.currentText()
        basis = self.basis_set.currentText()
        
        # 3c methods usually don't need basis set
        mtype = self.get_inferred_category(self.method_name.currentText())
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
        elif "NMR" in job_txt: route_parts.append("NMR")
        elif "SP" in job_txt: pass # No keyword
        
        # Opt Options
        if self.opt_group.isVisible():
            if self.opt_tight.isChecked(): route_parts.append("TightOpt")
            if self.opt_verytight.isChecked(): route_parts.append("VeryTightOpt")
            if self.opt_loose.isChecked(): route_parts.append("LooseOpt")
            if self.opt_cart.isChecked(): route_parts.append("COpt")
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
        if self.scf_sloppy.isChecked(): route_parts.append("SloppySCF")
        elif self.scf_loose.isChecked(): route_parts.append("LooseSCF")
        elif self.scf_normal.isChecked(): route_parts.append("NormalSCF")
        elif self.scf_strong.isChecked(): route_parts.append("StrongSCF")
        elif self.scf_tight.isChecked(): route_parts.append("TightSCF")
        elif self.scf_verytight.isChecked(): route_parts.append("VeryTightSCF")
        elif self.scf_extreme.isChecked(): route_parts.append("ExtremeSCF")
        
        grid = self.grid_combo.currentText()
        if grid != "Default": route_parts.append(grid)
        
        # NBO
        if self.pop_nbo.isChecked(): route_parts.append("NBO")

        self.preview_str = " ".join(route_parts)
        
        # Add %tddft block if enabled
        if self.tddft_enable.isChecked():
            block = (
                f"\n\n%tddft\n"
                f"  NRoots {self.tddft_nroots.value()}\n"
                f"  Singlets {'true' if self.tddft_singlets.isChecked() else 'false'}\n"
                f"  Triplets {'true' if self.tddft_triplets.isChecked() else 'false'}\n"
                f"  TDA {'true' if self.tddft_tda.isChecked() else 'false'}\n"
                f"  IRoot {self.tddft_iroot.value()}\n"
                f"end"
            )
            self.preview_str += block

        self.preview_label.setText(self.preview_str)

    def get_route(self):
        return self.preview_str
        
    def parse_route(self, route):
        if not route: return
        self.ui_ready = False 
        
        # Normalize route
        cleaned_route = route.strip()
        if cleaned_route.startswith("!"):
            cleaned_route = cleaned_route[1:].strip()
            
        tokens = cleaned_route.split()
        if not tokens: 
            self.ui_ready = True
            return

        method_list_upper = [m.upper() for m in ALL_ORCA_METHODS]
        basis_list_upper = [b.upper() for b in ALL_ORCA_BASIS_SETS]
        
        found_method = False
        found_basis = False
        
        for t in tokens:
            tu = t.upper()
            
            # 1. Method & Basis (Priority)
            if not found_method and tu in method_list_upper:
                idx = method_list_upper.index(tu)
                self.method_name.setCurrentText(ALL_ORCA_METHODS[idx])
                found_method = True
                continue
            elif not found_basis and tu in basis_list_upper:
                idx = basis_list_upper.index(tu)
                self.basis_set.setCurrentText(ALL_ORCA_BASIS_SETS[idx])
                found_basis = True
                continue
            
            # 2. Job Types
            if tu == "OPT":
                 if self.job_type.currentText() == "Frequency Only (Freq)":
                      self.job_type.setCurrentText("Optimization + Freq (Opt Freq)")
                 else:
                      self.job_type.setCurrentText("Optimization Only (Opt)")
            elif tu == "FREQ": 
                if self.job_type.currentText() == "Optimization Only (Opt)":
                    self.job_type.setCurrentText("Optimization + Freq (Opt Freq)")
                else:
                    self.job_type.setCurrentText("Frequency Only (Freq)")
            elif tu == "OPTTS": self.job_type.setCurrentText("Transition State Opt (OptTS)")
            elif tu == "SCAN": self.job_type.setCurrentText("Scan (Relaxed Surface)")
            elif tu == "NMR": self.job_type.setCurrentText("NMR")
            elif tu in ["GRADIENT", "HESSIAN"]:
                 # Direct match for Gradient/Hessian
                 for i in range(self.job_type.count()):
                     if tu in self.job_type.itemText(i).upper():
                         self.job_type.setCurrentIndex(i)
                         break
            
            # 3. Opt Options
            if tu == "TIGHTOPT": self.opt_tight.setChecked(True)
            elif tu == "VERYTIGHTOPT": self.opt_verytight.setChecked(True)
            elif tu == "LOOSEOPT": self.opt_loose.setChecked(True)
            elif tu == "COPT": self.opt_cart.setChecked(True)
            elif tu == "CALCFC": self.opt_calcfc.setChecked(True)
            elif tu == "CALCHESS": self.opt_ts_mode.setChecked(True)
            
            # 4. Freq Options
            if tu == "NUMFREQ": self.freq_num.setChecked(True)
            if tu == "RAMAN": self.freq_raman.setChecked(True)
            
            # 5. Solvation
            if "(" in tu and any(x in tu for x in ["CPCM", "SMD", "IEFPCM"]):
                s_model = tu.split("(")[0]
                s_name = t.split("(")[1].split(")")[0]
                
                if s_model == "CPCM": self.solv_model.setCurrentText("CPCM")
                elif s_model == "SMD": self.solv_model.setCurrentText("SMD")
                elif s_model == "IEFPCM": self.solv_model.setCurrentText("IEFPCM")
                
                for i in range(self.solvent.count()):
                    if self.solvent.itemText(i).upper() == s_name.upper():
                        self.solvent.setCurrentIndex(i)
                        break
            elif tu == "SMD": self.solv_model.setCurrentText("SMD")
            elif tu == "CPC(WATER)": self.solv_model.setCurrentText("CPC(Water) (Short)")
            
            # 6. Dispersion
            if tu in ["D3BJ", "D3ZERO", "D4", "D2", "NL"]:
                self.dispersion.setCurrentText(tu)
            
            # 7. RI / RIJCOSX
            if tu in ["RIJCOSX", "RI"]:
                self.rijcosx.setChecked(True)
            
            # 8. Aux Basis
            if tu == "DEF2/J": self.aux_basis.setCurrentText("Def2/J")
            elif tu == "DEF2/JK": self.aux_basis.setCurrentText("Def2/JK")
            elif tu == "AUTOMX": self.aux_basis.setCurrentText("AutoAux")
            
            # 9. SCF / Grid
            if tu == "SLOPPYSCF": self.scf_sloppy.setChecked(True)
            elif tu == "LOOSESCF": self.scf_loose.setChecked(True)
            elif tu == "NORMALSCF": self.scf_normal.setChecked(True)
            elif tu == "STRONGSCF": self.scf_strong.setChecked(True)
            elif tu == "TIGHTSCF": self.scf_tight.setChecked(True)
            elif tu == "VERYTIGHTSCF": self.scf_verytight.setChecked(True)
            elif tu == "EXTREMESCF": self.scf_extreme.setChecked(True)
            
            for combo in [self.grid_combo]:
                for i in range(combo.count()):
                    if combo.itemText(i).lower() == tu.lower():
                        combo.setCurrentIndex(i)
                        break
            
            # 10. NBO
            if tu == "NBO": self.pop_nbo.setChecked(True)

        self.ui_ready = True
        self.update_ui_state()
        self.update_preview()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_Escape:
            focused = self.focusWidget()
            if isinstance(focused, (QLineEdit, QComboBox, QSpinBox, QTextEdit)):
                focused.clearFocus()
                return
        super().keyPressEvent(event)

class OrcaSetupDialogNeo(QDialog):
    """
    ORCA Input Generator Neo
    """
    def __init__(self, parent=None, mol=None, filename=None):
        super().__init__(parent)
        self.setWindowTitle(PLUGIN_NAME)
        self.resize(1100, 800)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.WindowMinMaxButtonsHint)
        self.setSizeGripEnabled(True)
        self.mol = mol
        self.filename = filename
        self.ui_ready = False
        self.setup_ui()
        self.ui_ready = True
        self.load_presets_from_file()
        self.calc_initial_charge_mult()

    def setup_ui(self):
        main_layout = QVBoxLayout()
        
        # --- Horizontal Split Layer ---
        content_layout = QHBoxLayout()
        
        # --- Right Side: Settings (Scrollable) ---
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QScrollArea.Shape.NoFrame)
        
        settings_container = QWidget()
        settings_layout = QVBoxLayout(settings_container)

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
        settings_layout.addWidget(preset_group)

        # --- 1. Resource Configuration ---
        res_group = QGroupBox("Resources (%pal, %maxcore)")
        res_layout = QFormLayout()
        res_h1 = QHBoxLayout()
        self.nproc_spin = QSpinBox()
        self.nproc_spin.setRange(1, 128)
        self.nproc_spin.setValue(4)
        self.nproc_spin.valueChanged.connect(self.update_preview)
        btn_auto_proc = QPushButton("Auto")
        btn_auto_proc.setFixedWidth(50)
        btn_auto_proc.clicked.connect(self.auto_detect_nproc)
        res_h1.addWidget(self.nproc_spin)
        res_h1.addWidget(btn_auto_proc)
        res_layout.addRow("Processors:", res_h1)

        res_h2 = QHBoxLayout()
        self.mem_spin = QSpinBox()
        self.mem_spin.setRange(100, 999999)
        self.mem_spin.setValue(2000) 
        self.mem_spin.setSuffix(" MB")
        self.mem_spin.valueChanged.connect(self.update_preview)
        btn_auto_mem = QPushButton("Auto")
        btn_auto_mem.setFixedWidth(50)
        btn_auto_mem.clicked.connect(self.auto_detect_mem)
        res_h2.addWidget(self.mem_spin)
        res_h2.addWidget(btn_auto_mem)
        res_layout.addRow("MaxCore:", res_h2)
        res_group.setLayout(res_layout)
        settings_layout.addWidget(res_group)

        # --- 2. Simple Input Line ---
        kw_group = QGroupBox("Simple Input Line (!)")
        kw_layout = QVBoxLayout()
        kw_h_layout = QHBoxLayout()
        self.keywords_edit = QTextEdit("! B3LYP def2-SVP RIJCOSX Def2/J Opt Freq")
        self.keywords_edit.setFixedHeight(70)
        self.keywords_edit.textChanged.connect(self.update_preview)
        self.kw_highlighter = OrcaSyntaxHighlighter(self.keywords_edit.document())
        self.btn_route = QPushButton("Builder...")
        self.btn_route.clicked.connect(self.open_keyword_builder)
        kw_h_layout.addWidget(self.keywords_edit)
        kw_h_layout.addWidget(self.btn_route)
        kw_layout.addWidget(QLabel("Keywords:"))
        kw_layout.addLayout(kw_h_layout)
        self.comment_edit = QLineEdit("Generated by MoleditPy")
        self.comment_edit.textChanged.connect(self.update_preview)
        kw_layout.addWidget(QLabel("Comment:"))
        kw_layout.addWidget(self.comment_edit)
        kw_group.setLayout(kw_layout)
        settings_layout.addWidget(kw_group)

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
        mol_layout.addWidget(QLabel("Mult:"))
        mol_layout.addWidget(self.mult_spin)
        self.default_palette = self.charge_spin.palette()
        mol_group.setLayout(mol_layout)
        settings_layout.addWidget(mol_group)

        # --- 3b. Coordinate Format ---
        coord_group = QGroupBox("Coordinate Format")
        coord_layout = QHBoxLayout()
        self.coord_format_combo = QComboBox()
        self.coord_format_combo.addItems(["Cartesian (XYZ)", "Internal (* int)", "Internal (* gzmt)"])
        self.coord_format_combo.currentIndexChanged.connect(self.update_preview)
        coord_layout.addWidget(self.coord_format_combo)
        coord_group.setLayout(coord_layout)
        settings_layout.addWidget(coord_group)

        # --- 4. Advanced/Blocks ---
        adv_group = QGroupBox("Advanced Blocks")
        adv_layout = QVBoxLayout()
        blk_h_layout = QHBoxLayout()
        self.block_combo = QComboBox()
        self.block_combo.addItems([
            "Select Block to Insert...",
            "%output (Basis/MOs)",
            "%eprnmr (NMR/J-coupling)",
            "%scf ... end",
            "%geom ... end",
            "%elprop ... end",
            "%plots ... end",
            "%tddft ... end",
            "%cis ... end", 
            "%mrci ... end",
            "%casscf ... end"
        ])
        blk_h_layout.addWidget(self.block_combo, 1)
        self.btn_insert_block = QPushButton("Insert")
        self.btn_insert_block.clicked.connect(self.insert_block_template)
        blk_h_layout.addWidget(self.btn_insert_block)
        adv_layout.addLayout(blk_h_layout)
        self.adv_tabs = QTabWidget()
        self.adv_tabs.setFixedHeight(120)
        self.adv_edit = QTextEdit()
        self.adv_edit.textChanged.connect(self.update_preview)
        self.adv_tabs.addTab(self.adv_edit, "Pre-Coord")
        self.post_adv_edit = QTextEdit()
        self.post_adv_edit.textChanged.connect(self.update_preview)
        self.adv_tabs.addTab(self.post_adv_edit, "Post-Coord")
        adv_layout.addWidget(self.adv_tabs)
        adv_group.setLayout(adv_layout)
        settings_layout.addWidget(adv_group)

        settings_layout.addStretch()
        scroll_area.setWidget(settings_container)
        content_layout.addWidget(scroll_area, 3) # Left side settings

        # --- Right Side: Preview ---
        preview_group = QGroupBox("Input Preview")
        preview_layout = QVBoxLayout()
        self.preview_text = QTextEdit()
        self.preview_text.setReadOnly(False)
        self.preview_text.setFont(QFont("Courier New", 10))
        self.preview_text.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.highlighter = OrcaSyntaxHighlighter(self.preview_text.document())
        preview_layout.addWidget(self.preview_text)
        btn_refresh = QPushButton("Reset/Refresh Preview")
        btn_refresh.clicked.connect(self.update_preview)
        preview_layout.addWidget(btn_refresh)
        preview_group.setLayout(preview_layout)
        content_layout.addWidget(preview_group, 4) # Right side preview
        
        main_layout.addLayout(content_layout)

        # --- Save Button ---
        self.save_btn = QPushButton("Save ORCA Input File...")
        self.save_btn.clicked.connect(self.save_file)
        self.save_btn.setStyleSheet("font-weight: bold; padding: 8px; font-size: 14px;")
        main_layout.addWidget(self.save_btn)
        
        self.setLayout(main_layout)
        self.update_preview()
        
    def update_preview(self):
        if not getattr(self, 'ui_ready', False):
            return
        self.preview_text.setText(self.generate_input_content())

    def auto_detect_nproc(self):
        try:
            import psutil
            cpus = psutil.cpu_count(logical=False)
        except:
            cpus = os.cpu_count()
            
        if cpus:
            self.nproc_spin.setValue(cpus)

    def auto_detect_mem(self):
        # Default fallback
        total_mb = 8000
        try:
            import psutil
            # Total system memory
            total_mb = psutil.virtual_memory().total // (1024 * 1024)
        except:
             pass
        
        # Use roughly 75-80% of total memory to be safe, divided by nprocs
        nprocs = self.nproc_spin.value()
        rec_core = int((total_mb * 0.80) / nprocs)
        self.mem_spin.setValue(max(500, rec_core))

    def preview_file(self):
        # Legacy stub
        self.update_preview()

    def save_file(self):
        # 1. 座標データの取得 (Check for error first)
        coord_style = self.coord_format_combo.currentText()
        if "XYZ" in coord_style:
            lines = self.get_coords_lines()
        elif "gzmt" in coord_style:
             lines = self.get_zmatrix_gzmt_lines()
        else:
             lines = self.get_zmatrix_standard_lines()

        if any("Error" in l for l in lines):
            err = "\\n".join(l for l in lines if "Error" in l)
            QMessageBox.critical(self, "Error", f"Coordinate Generation Failed:\\n{err}")
            return

        # 2. ファイル保存ダイアログ
        default_name = "orca_job.inp"
        if self.filename:
            base = os.path.splitext(os.path.basename(self.filename))[0]
            default_name = f"{base}.inp"
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save ORCA Input", default_name, "ORCA Input (*.inp);;All Files (*)"
        )

        if file_path:
            try:
                # Use content from preview (editable)
                content = self.preview_text.toPlainText()
                
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(content)

                QMessageBox.information(self, "Success", f"File saved:\\n{file_path}")
                #QMessageBox.information(self, "Success", f"File saved:\n{file_path}")
                # Do not close automatically
                # self.accept()
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save file:\n{str(e)}")

    def insert_block_template(self):
        txt = self.block_combo.currentText()
        if "Select" in txt: return
        
        template = ""
        if "%scf" in txt:
             template = "%scf\n MaxIter 125\nend\n"
        elif "%output" in txt:
             template = "%output\n  Print[P_Basis] 2  # Required for Basis Set parsing\n  Print[P_Mos] 1    # Ensure MO coefficients are printed\nend\n"
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
                 "  Triplets false  # Calculate triplet states\n"
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
             template = "%eprnmr\n  NUCLEI = ALL H {SHIFT, SSALL} # Required for J-coupling (nmrsim)\nend\n"
             # Switch to Post-Coordinate tab automatically
             self.adv_tabs.setCurrentWidget(self.post_adv_edit)
        
        current_widget = self.adv_tabs.currentWidget()
        if isinstance(current_widget, QTextEdit):
            cursor = current_widget.textCursor()
            cursor.insertText(template)

    def open_keyword_builder(self):
        dialog = OrcaKeywordBuilderDialog(self, self.keywords_edit.toPlainText())
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.keywords_edit.setPlainText(dialog.get_route())

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
        res_part = []
        nprocs = self.nproc_spin.value()
        if nprocs > 1:
            res_part.append(f"%pal nprocs {nprocs} end")
        res_part.append(f"%maxcore {self.mem_spin.value()}")
        content.append("\n".join(res_part))
        
        # Keywords
        kw = self.keywords_edit.toPlainText().strip()
        if not kw.startswith("!"): kw = "! " + kw
        content.append(f"\n{kw}")
        
        # Advanced Blocks
        adv = self.adv_edit.toPlainText().strip()
        if adv:
            content.append(f"\n{adv}")
        
        # Coordinates
        is_cartesian = "Cartesian" in self.coord_format_combo.currentText()
        coord_lines = self.get_coords_lines()
        
        content.append("") # Blank line before coordinates
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
                "route": "! B3LYP def2-SVP RIJCOSX Def2/J Opt Freq", "adv": "", "adv_post": ""
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
        self.keywords_edit.setText(data.get("route", "! B3LYP def2-SVP RIJCOSX Def2/J Opt Freq"))
        self.adv_edit.setPlainText(data.get("adv", ""))
        self.post_adv_edit.setPlainText(data.get("adv_post", ""))
            
        self.update_preview()

    def save_preset_dialog(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:")
        if ok and name:
            self.presets_data[name] = {
                "nproc": self.nproc_spin.value(),
                "maxcore": self.mem_spin.value(),
                "route": self.keywords_edit.toPlainText(),
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
            
            self.update_preview()
            
        except: pass

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_Escape:
            focused = self.focusWidget()
            if isinstance(focused, (QLineEdit, QComboBox, QSpinBox, QTextEdit)):
                focused.clearFocus()
                return
        super().keyPressEvent(event)

def run(mw):
    mol = getattr(mw, 'current_mol', None)
    if not mol:
        QMessageBox.warning(mw, PLUGIN_NAME, "No molecule loaded.")
        return
        
    filename = getattr(mw, 'current_file_path', None)
    dialog = OrcaSetupDialogNeo(parent=mw, mol=mol, filename=filename)
    dialog.exec()

def initialize(context):
    def show_dialog():
        mw = context.get_main_window()
        run(mw)
    context.add_export_action("ORCA Input...", show_dialog)

