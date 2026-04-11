from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTableWidget, 
    QTableWidgetItem, QPushButton, QMessageBox, QLabel, QHeaderView, QAbstractItemView,
    QApplication, QComboBox, QCheckBox
)
from PyQt6.QtCore import Qt
from rdkit.Chem import AllChem
import copy

PLUGIN_NAME = "Conformational Search"
__version__="2026.04.04"
__author__="HiroYokoyama"

class ConformerSearchDialog(QDialog):
    def __init__(self, context, parent=None):
        super().__init__(parent)
        self.context = context
        self.main_window = context.get_main_window()
        self.setWindowTitle("Conformational Search & Preview")
        self.resize(400, 500)
        
        # メインウィンドウの分子への参照
        self.target_mol = context.current_mol
        
        # 計算用の一時的な分子（オリジナルを汚染しないため）
        self.temp_mol = None
        # 生成された配座データのリスト [(Energy, ConformerID), ...]
        self.conformer_data = []
        # 全ての計算結果（未フィルタ）
        self.results_raw = []
        
        # Original coordinates for restoration on cancel
        self.original_coords = []
        if self.target_mol:
            conf = self.target_mol.GetConformer()
            self.original_coords = [conf.GetAtomPosition(i) for i in range(self.target_mol.GetNumAtoms())]
        
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # 説明ラベル
        self.lbl_info = QLabel("Click 'Run Search' to generate conformers.\nSelect a row to preview.")
        layout.addWidget(self.lbl_info)

        # Force Field Selection
        hbox_ff = QHBoxLayout()
        hbox_ff.addWidget(QLabel("Force Field:"))
        self.combo_ff = QComboBox()
        self.combo_ff.addItems(["MMFF94", "UFF"])
        hbox_ff.addWidget(self.combo_ff)
        hbox_ff.addStretch()
        layout.addLayout(hbox_ff)

        # Show All Checkbox (Moving under FF box)
        self.cb_show_all = QCheckBox("Show All Conformers (ignore energy redundancy)")
        self.cb_show_all.setChecked(False)
        self.cb_show_all.toggled.connect(self.apply_filter_and_update)
        layout.addWidget(self.cb_show_all)
        
        # Set default based on main window setting
        default_method = self.context.get_setting("optimization_method", "MMFF_RDKIT")
        if default_method:
            default_method = default_method.upper()
            if "UFF" in default_method:
                self.combo_ff.setCurrentText("UFF")
            else:
                self.combo_ff.setCurrentText("MMFF94")

        # 結果表示用テーブル
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Rank", "Energy (kcal/mol)"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.currentItemChanged.connect(self.preview_conformer)
        layout.addWidget(self.table)

        # ボタンエリア
        btn_layout = QHBoxLayout()
        self.btn_run = QPushButton("Run Search")
        self.btn_run.clicked.connect(self.run_search)
        
        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.accept) # 閉じる（現在のプレビュー状態で確定）

        btn_layout.addWidget(self.btn_run)
        btn_layout.addWidget(self.btn_close)
        layout.addLayout(btn_layout)

    def accept(self):
        # Push undo state when closing the dialog (confirming the selection)
        if self.target_mol:
            self.context.push_undo_checkpoint()
        super().accept()

    def reject(self):
        # Restore original coordinates if user cancels/closes without 'Accept'
        if self.target_mol and self.original_coords:
            conf = self.target_mol.GetConformer()
            for i, pos in enumerate(self.original_coords):
                conf.SetAtomPosition(i, pos)
            self.context.refresh_3d_view()
        super().reject()

    def run_search(self):
        if not self.target_mol:
            return

        self.btn_run.setEnabled(False)
        self.lbl_info.setText("Running conformational search... please wait.")
        QApplication.processEvents()

        try:
            # 計算用に分子を複製（水素が付加されていることを推奨）
            mol_calc = copy.deepcopy(self.target_mol)
            
            # 1. 配座生成 (ETKDGv3)
            params = AllChem.ETKDGv3()
            params.useSmallRingTorsions = True
            cids = AllChem.EmbedMultipleConfs(mol_calc, numConfs=30, params=params)

            if not cids:
                QMessageBox.warning(self, PLUGIN_NAME, "Failed to generate conformers.")
                self.lbl_info.setText("Failed.")
                self.btn_run.setEnabled(True)
                return

            # 2. 構造最適化とエネルギー計算
            results = []
            selected_ff = self.combo_ff.currentText()
            
            for i, cid in enumerate(cids):
                energy = None
                
                if selected_ff == "MMFF94":
                     # MMFF94 Optimize
                    if AllChem.MMFFOptimizeMolecule(mol_calc, confId=cid) != -1:
                         # Calculate Energy
                         prop = AllChem.MMFFGetMoleculeProperties(mol_calc)
                         if prop:
                             ff = AllChem.MMFFGetMoleculeForceField(mol_calc, prop, confId=cid)
                             if ff:
                                 energy = ff.CalcEnergy()
                
                elif selected_ff == "UFF":
                     # UFF Optimize
                     if AllChem.UFFOptimizeMolecule(mol_calc, confId=cid) != -1:
                         # Calculate Energy
                         ff = AllChem.UFFGetMoleculeForceField(mol_calc, confId=cid)
                         if ff:
                             energy = ff.CalcEnergy()

                if energy is not None:
                    results.append((energy, cid))
                
                # UIの応答性を維持
                if i % 5 == 0:
                    QApplication.processEvents()

            if not results:
                 QMessageBox.warning(self, PLUGIN_NAME, f"Optimization failed with {selected_ff}.")
                 self.btn_run.setEnabled(True)
                 return

            # エネルギーが低い順にソート
            results.sort(key=lambda x: x[0])
            self.results_raw = results

            # データを保持 & 更新
            self.temp_mol = mol_calc
            self.apply_filter_and_update()

        except Exception as e:
            QMessageBox.critical(self, PLUGIN_NAME, f"Error during search: {str(e)}")
            self.lbl_info.setText("Error occurred.")
        finally:
            self.btn_run.setEnabled(True)

    def apply_filter_and_update(self):
        """現在のフィルタ設定に基づいてデータを抽出し、テーブルを更新する"""
        if not self.results_raw:
            return

        if self.cb_show_all.isChecked():
            self.conformer_data = self.results_raw
        else:
            # エネルギーの重複を排除
            filtered = []
            ENERGY_THRESHOLD = 0.0001
            for energy, cid in self.results_raw:
                is_redundant = False
                for existing_energy, _ in filtered:
                    if abs(energy - existing_energy) < ENERGY_THRESHOLD:
                        is_redundant = True
                        break
                if not is_redundant:
                    filtered.append((energy, cid))
            self.conformer_data = filtered

        self.update_table()
        self.lbl_info.setText(f"Showing {len(self.conformer_data)} conformers (Total found: {len(self.results_raw)}).")

    def update_table(self):
        self.table.setRowCount(0)
        # base_energy = self.conformer_data[0][0] if self.conformer_data else 0

        for rank, (energy, cid) in enumerate(self.conformer_data):
            row_idx = self.table.rowCount()
            self.table.insertRow(row_idx)
            
            # Rank
            self.table.setItem(row_idx, 0, QTableWidgetItem(str(rank + 1)))
            
            # Energy
            energy_str = f"{energy:.4f}"
            self.table.setItem(row_idx, 1, QTableWidgetItem(energy_str))
            
            # 隠しデータとしてConformer IDを持たせる
            self.table.item(row_idx, 0).setData(Qt.ItemDataRole.UserRole, cid)

    def preview_conformer(self, current, previous):
        """リスト選択時にメインウィンドウの表示を更新"""
        if not self.temp_mol or not self.target_mol or not current:
            return
        
        row = current.row()
        # Rankカラム(0)にCIDを埋め込んでいるので取得
        cid = self.table.item(row, 0).data(Qt.ItemDataRole.UserRole)
        
        # 選択された配座の座標を取得
        source_conf = self.temp_mol.GetConformer(cid)
        target_conf = self.target_mol.GetConformer() # 現在の表示用Conformer
        
        # Safety check: Atom count must match
        if self.temp_mol.GetNumAtoms() != self.target_mol.GetNumAtoms():
             self.lbl_info.setText("<font color='red'>Error: Molecule changed in main window. Restart search.</font>")
             return

        # 座標のコピー
        for i in range(self.target_mol.GetNumAtoms()):
            pos = source_conf.GetAtomPosition(i)
            target_conf.SetAtomPosition(i, pos)
            
        # ビューの更新 (重要: V3ではdraw_molecule_3dを呼ぶことで座標更新を反映させる)
        self.context.current_mol = self.target_mol
        self.context.refresh_3d_view()

def run_plugin(context):
    mw = context.get_main_window()
    if not context.current_mol:
        QMessageBox.warning(mw, PLUGIN_NAME, "No molecule loaded.")
        return
        
    win = context.get_window("main_panel")
    if win:
        win.show()
        win.raise_()
        win.activateWindow()
        return

    dialog = ConformerSearchDialog(context, parent=mw)
    context.register_window("main_panel", dialog)
    dialog.show()

_launch_fn = None

def initialize(context):
    """Register the plugin in the 3D Edit menu."""
    global _launch_fn
    _launch_fn = lambda: run_plugin(context)
    context.add_menu_action("3D Edit/Conformational Search...", _launch_fn)

def run(mw):
    if hasattr(mw, 'host'):
        mw = mw.host
    if _launch_fn:
        _launch_fn()
