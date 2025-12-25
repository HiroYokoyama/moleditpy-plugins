import random
import math
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, 
    QMessageBox, QLabel, QProgressBar, QSpinBox,
    QGroupBox, QFormLayout, QComboBox
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

PLUGIN_NAME = "Complex Molecule Untangler"
__version__="2025.12.25"
__author__="HiroYokoyama"


class UntangleWorker(QThread):
    """
    回転可能な結合をランダムに回して、衝突（エネルギー）が減る配置を探すスレッド
    """
    progress = pyqtSignal(int)
    finished = pyqtSignal(object, str)
    
    def __init__(self, mol, max_iter=500, force_field="MMFF94"):
        super().__init__()
        self.mol = mol
        self.max_iter = max_iter
        self.force_field = force_field

    def run(self):
        try:
            # 元の分子をコピーして操作
            work_mol = Chem.Mol(self.mol)
            
            # フォースフィールドのセットアップ
            ff = None
            if self.force_field == "MMFF94":
                try:
                    props = AllChem.MMFFGetMoleculeProperties(work_mol)
                    if props:
                        ff = AllChem.MMFFGetMoleculeForceField(work_mol, props)
                except:
                    pass
            elif self.force_field == "UFF":
                try:
                    ff = AllChem.UFFGetMoleculeForceField(work_mol)
                except:
                    pass

            # フォースフィールド構築失敗時のフォールバックなどは今回は厳密にしない（エラー通知）
            if not ff:
                msg = f"Could not setup Force Field ({self.force_field})."
                if self.force_field == "MMFF94":
                    msg += "\nTry using UFF if MMFF94 parameters are missing."
                self.finished.emit(None, msg)
                return

            # 初期エネルギー（衝突具合）
            current_energy = ff.CalcEnergy()
            
            # 回転可能な結合（二面角）を探索
            # SMARTSパターン: アミド結合などを除外した、厳密な回転可能結合
            rotatable_smarts = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
            matches = work_mol.GetSubstructMatches(rotatable_smarts)
            
            if not matches:
                self.finished.emit(None, "No rotatable bonds found.")
                return

            # 結合ごとの4原子インデックス(i, j, k, l)のリストを作成
            dihedrals = []
            for (j, k) in matches:
                # j-k が回転軸。それぞれの隣接原子 i, l を探す
                atom_j = work_mol.GetAtomWithIdx(j)
                atom_k = work_mol.GetAtomWithIdx(k)
                
                # jの隣接原子でk以外
                nbrs_j = [n.GetIdx() for n in atom_j.GetNeighbors() if n.GetIdx() != k]
                # kの隣接原子でj以外
                nbrs_k = [n.GetIdx() for n in atom_k.GetNeighbors() if n.GetIdx() != j]
                
                if nbrs_j and nbrs_k:
                    dihedrals.append((nbrs_j[0], j, k, nbrs_k[0]))

            if not dihedrals:
                self.finished.emit(None, "Could not define dihedrals.")
                return

            # --- メインループ: ランダム回転による衝突回避 ---
            conf = work_mol.GetConformer()
            
            for i in range(self.max_iter):
                # ランダムに1つ選択
                i_idx, j_idx, k_idx, l_idx = random.choice(dihedrals)
                
                # 現在の角度を保存
                old_angle = rdMolTransforms.GetDihedralDeg(conf, i_idx, j_idx, k_idx, l_idx)
                
                # ランダムに回転 (-180度 〜 +180度 の範囲で新しい角度を決定)
                new_angle = random.uniform(-180, 180)
                
                # 回転適用
                rdMolTransforms.SetDihedralDeg(conf, i_idx, j_idx, k_idx, l_idx, new_angle)
                
                # 判定
                # 座標が変わったのでFFを更新する必要があるか？ -> RDKitのFFは座標更新を追跡しない場合があるが、
                # CalcEnergyは現在のCoordsを使うはず。ただしInitializeが必要な場合も。
                # RDKit通常の使用法では座標を変えたらそのままCalcEnergyで反映される。
                new_energy = ff.CalcEnergy()
                
                if new_energy < current_energy:
                    # 改善した（衝突が減った） -> 採用
                    current_energy = new_energy
                else:
                    # 悪化した（ぶつかった） -> 元に戻す
                    rdMolTransforms.SetDihedralDeg(conf, i_idx, j_idx, k_idx, l_idx, old_angle)

                # 進捗通知
                self.progress.emit(i + 1)

            # 最後に軽く整列（微調整）して仕上げ
            # 選択されたFFで最適化
            if self.force_field == "MMFF94":
                try:
                    AllChem.MMFFOptimizeMolecule(work_mol, maxIters=50)
                except: 
                    pass
            elif self.force_field == "UFF":
                try:
                    AllChem.UFFOptimizeMolecule(work_mol, maxIters=50)
                except:
                    pass

            self.finished.emit(work_mol, f"Processed {len(matches)} bonds.\nFinal Score ({self.force_field}): {current_energy:.2f}")

        except Exception as e:
            self.finished.emit(None, str(e))

class UntanglerDialog(QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setWindowTitle("Complex Molecule Untangler")
        self.resize(320, 350)
        
        self.worker = None
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # --- Header Section ---
        title_label = QLabel("Complex Molecule Untangler")
        title_font = QFont()
        title_font.setBold(True)
        title_font.setPointSize(11)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title_label)

        sub_label = QLabel("Resolves steric clashes by randomly rotating single bonds (Monte Carlo).")
        sub_label.setWordWrap(True)
        sub_label.setStyleSheet("color: #555; margin-bottom: 10px;")
        layout.addWidget(sub_label)

        # --- Configuration Group ---
        config_group = QGroupBox("Configuration")
        form_layout = QFormLayout()

        # Force Field Selection
        self.combo_ff = QComboBox()
        self.combo_ff.addItems(["MMFF94", "UFF"])
        self.combo_ff.setToolTip("Select the Force Field for energy calculation.")
        form_layout.addRow("Force Field:", self.combo_ff)
        
        # Set default based on main window setting
        default_method = getattr(self.main_window, "optimization_method", "MMFF_RDKIT")
        if default_method:
            default_method = default_method.upper()
            if "UFF" in default_method:
                self.combo_ff.setCurrentText("UFF")
            else:
                self.combo_ff.setCurrentText("MMFF94")

        # Max Iterations
        self.spin_iter = QSpinBox()
        self.spin_iter.setRange(100, 10000)
        self.spin_iter.setValue(500)
        self.spin_iter.setSingleStep(100)
        self.spin_iter.setToolTip("Number of random rotation attempts.")
        form_layout.addRow("Max Iterations:", self.spin_iter)

        config_group.setLayout(form_layout)
        layout.addWidget(config_group)

        # --- Progress Area ---
        self.pbar = QProgressBar()
        self.pbar.setValue(0)
        self.pbar.setTextVisible(False)
        layout.addWidget(self.pbar)

        # --- Action Area ---
        self.btn_run = QPushButton("Untangle Molecule")
        self.btn_run.setMinimumHeight(40)
        font_btn = QFont()
        font_btn.setBold(True)
        self.btn_run.setFont(font_btn)
        self.btn_run.clicked.connect(self.run_untangle)
        layout.addWidget(self.btn_run)

    def run_untangle(self):
        mol = getattr(self.main_window, "current_mol", None)
        if not mol:
            QMessageBox.warning(self, PLUGIN_NAME, "No molecule loaded.\nPlease load a molecule first.")
            return

        self.btn_run.setEnabled(False)
        self.btn_run.setText("Processing...")
        
        max_iter = self.spin_iter.value()
        ff_choice = self.combo_ff.currentText()
        
        self.pbar.setRange(0, max_iter)
        self.pbar.setValue(0)
        self.pbar.setTextVisible(True)

        self.worker = UntangleWorker(mol, max_iter=max_iter, force_field=ff_choice)
        self.worker.progress.connect(self.pbar.setValue)
        self.worker.finished.connect(self.on_finished)
        self.worker.start()

    def on_finished(self, new_mol, msg):
        self.btn_run.setEnabled(True)
        self.btn_run.setText("Untangle Molecule")
        self.pbar.setTextVisible(False)
        self.pbar.setValue(0)

        if new_mol:
            self.main_window.current_mol = new_mol
            
            # ビュー更新
            if hasattr(self.main_window, "draw_molecule_3d"):
                self.main_window.draw_molecule_3d(new_mol)
            elif hasattr(self.main_window, "update_view"):
                self.main_window.update_view()
            elif hasattr(self.main_window, "gl_widget"):
                getattr(self.main_window.gl_widget, "update", lambda: None)()
            
            # Push undo state using the newly applied molecule
            if hasattr(self.main_window, "push_undo_state"):
                self.main_window.push_undo_state()
                
            QMessageBox.information(self, PLUGIN_NAME, f"Untangling Complete!\n{msg}")
        else:
            QMessageBox.warning(self, PLUGIN_NAME, f"Error: {msg}")

def run(mw):
    if hasattr(mw, "_untangler_dialog") and mw._untangler_dialog.isVisible():
        mw._untangler_dialog.raise_()
        mw._untangler_dialog.activateWindow()
        return

    dialog = UntanglerDialog(mw, parent=mw)
    mw._untangler_dialog = dialog
    dialog.show()

# initialize removed as it only registered the menu action