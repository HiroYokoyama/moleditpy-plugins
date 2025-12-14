from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTableWidget, 
    QTableWidgetItem, QPushButton, QMessageBox, QLabel, QHeaderView, QAbstractItemView,
    QApplication, QComboBox
)
from PyQt6.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import AllChem
import copy

PLUGIN_NAME = "Conformational Search"

class ConformerSearchDialog(QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setWindowTitle("Conformational Search & Preview")
        self.resize(400, 500)
        
        # メインウィンドウの分子への参照
        self.target_mol = getattr(main_window, "current_mol", None)
        
        # 計算用の一時的な分子（オリジナルを汚染しないため）
        self.temp_mol = None
        # 生成された配座データのリスト [(Energy, ConformerID), ...]
        self.conformer_data = []

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
        
        # Set default based on main window setting
        default_method = getattr(self.main_window, "optimization_method", "MMFF_RDKIT")
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
        self.table.itemClicked.connect(self.preview_conformer)
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
            
            # データを保持
            self.temp_mol = mol_calc
            self.conformer_data = results

            # テーブル更新
            self.update_table()
            self.lbl_info.setText(f"Found {len(results)} conformers ({selected_ff}).")

        except Exception as e:
            QMessageBox.critical(self, PLUGIN_NAME, f"Error during search: {str(e)}")
            self.lbl_info.setText("Error occurred.")
        finally:
            self.btn_run.setEnabled(True)

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

    def preview_conformer(self, item):
        """リスト選択時にメインウィンドウの表示を更新"""
        if not self.temp_mol or not self.target_mol:
            return
        
        row = item.row()
        # Rankカラム(0)にCIDを埋め込んでいるので取得
        cid = self.table.item(row, 0).data(Qt.ItemDataRole.UserRole)
        
        # 選択された配座の座標を取得
        source_conf = self.temp_mol.GetConformer(cid)
        target_conf = self.target_mol.GetConformer() # 現在の表示用Conformer
        
        # 座標のコピー
        for i in range(self.target_mol.GetNumAtoms()):
            pos = source_conf.GetAtomPosition(i)
            target_conf.SetAtomPosition(i, pos)
            
        # ビューの更新（ユーザー提供コードのロジックに従う）
        if hasattr(self.main_window, "draw_molecule_3d"):
            self.main_window.draw_molecule_3d(self.target_mol)
        elif hasattr(self.main_window, "update_view"):
            self.main_window.update_view()
        elif hasattr(self.main_window, "gl_widget"):
            # GLWidgetのリフレッシュ
            getattr(self.main_window.gl_widget, "update", lambda: None)()

def run(main_window):
    """
    プラグインのエントリーポイント
    """
    mol = getattr(main_window, "current_mol", None)
    if not mol:
        QMessageBox.warning(main_window, PLUGIN_NAME, "No molecule loaded.")
        return
        
    # 既存のダイアログがあればアクティブにする
    if hasattr(main_window, "_conformer_search_dialog") and main_window._conformer_search_dialog.isVisible():
        main_window._conformer_search_dialog.raise_()
        main_window._conformer_search_dialog.activateWindow()
        return

    dialog = ConformerSearchDialog(main_window, parent=main_window)
    # 参照を保持してGCを防ぐ
    main_window._conformer_search_dialog = dialog
    dialog.show() # モーダルではなくModeless（非ブロック）で表示