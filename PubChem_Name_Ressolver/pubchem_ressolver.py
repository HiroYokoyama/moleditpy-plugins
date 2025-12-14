import sys
import requests # API通信に必要 (pip install requests)
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTableWidget, 
    QTableWidgetItem, QPushButton, QMessageBox, QLabel, QHeaderView, 
    QAbstractItemView, QApplication, QLineEdit, QComboBox
)
from PyQt6.QtCore import Qt, QPointF
from rdkit import Chem
from rdkit.Chem import AllChem

PLUGIN_NAME = "PubChem Name Resolver"

class MoleculeResolverDialog(QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setWindowTitle("PubChem Name Resolver")
        self.resize(500, 600)
        
        # 取得した候補データのリスト
        # 各要素は辞書: {'name': str, 'smiles': str, 'formula': str}
        self.candidates_data = []
        
        # 生成されたRDKit分子オブジェクト（一時保存）
        self.generated_mol = None

        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # --- 入力エリア ---
        input_layout = QHBoxLayout()
        
        self.combo_type = QComboBox()
        self.combo_type.addItems(["Auto (Name/CAS)", "SMILES"])
        input_layout.addWidget(self.combo_type)

        self.line_input = QLineEdit()
        self.line_input.setPlaceholderText("Enter Name or SMILES...")
        self.line_input.returnPressed.connect(self.run_search) # Enterキーで検索
        input_layout.addWidget(self.line_input)
        
        self.btn_search = QPushButton("Search Online")
        self.btn_search.clicked.connect(self.run_search)
        input_layout.addWidget(self.btn_search)
        
        layout.addLayout(input_layout)

        # --- 説明ラベル ---
        self.lbl_info = QLabel("Enter a chemical identifier and click Search.")
        layout.addWidget(self.lbl_info)

        # --- 結果表示用テーブル ---
        self.table = QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["Name/Synonym", "Formula", "SMILES"])
        # ヘッダー調整
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        # 選択変更時にロードボタンを有効化するなどの処理を入れるならここ
        layout.addWidget(self.table)

        # --- ボタンエリア ---
        btn_layout = QHBoxLayout()
        
        self.btn_load = QPushButton("Load to 2D Editor")
        self.btn_load.clicked.connect(self.load_molecule)
        
        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.close)

        btn_layout.addWidget(self.btn_load)
        btn_layout.addWidget(self.btn_close)
        layout.addLayout(btn_layout)

    def run_search(self):
        query = self.line_input.text().strip()
        if not query:
            return

        self.lbl_info.setText("Searching PubChem... please wait.")
        self.btn_search.setEnabled(False)
        self.table.setRowCount(0)
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        QApplication.processEvents()

        results = []
        error_msg = None
        network_error = False

        try:
            search_type = self.combo_type.currentText()

            if search_type == "SMILES":
                # SMILESの場合は直接リストに追加（検証含む）
                mol = Chem.MolFromSmiles(query)
                if mol:
                    results.append({
                        'name': 'User Input SMILES',
                        'smiles': query,
                        'formula': Chem.rdMolDescriptors.CalcMolFormula(mol)
                    })
                else:
                    raise ValueError("Invalid SMILES string.")
            
            else: # Auto (PubChem API)
                # PUG REST APIを使用して検索
                # プロパティとしてSMILESと分子式を取得
                # CanonicalSMILESも取得してフォールバックに使用
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/property/IsomericSMILES,CanonicalSMILES,MolecularFormula,Title/JSON"
                
                response = requests.get(url, timeout=10)
                
                if response.status_code == 200:
                    data = response.json()
                    props = data.get('PropertyTable', {}).get('Properties', [])
                    for p in props:
                        # SMILESの取得（Isomericを優先、なければCanonical、それでもなければキー検索）
                        smiles = p.get('IsomericSMILES')
                        if not smiles:
                            smiles = p.get('CanonicalSMILES')
                        
                        # まだ取得できていない場合、キー名に'SMILES'が含まれるものを探す
                        if not smiles:
                            for k in p.keys():
                                if 'SMILES' in k:
                                    smiles = p[k]
                                    if smiles:
                                        break
                        
                        if not smiles:
                            smiles = "" # 見つからない場合

                        results.append({
                            'name': p.get('Title', query),
                            'smiles': smiles,
                            'formula': p.get('MolecularFormula', '')
                        })
                else:
                    # found nothing or error status
                    pass

        except requests.exceptions.RequestException:
            network_error = True
        except Exception as e:
            error_msg = str(e)
        finally:
            QApplication.restoreOverrideCursor()
            self.btn_search.setEnabled(True)

        # UI updates after cursor restore
        if network_error:
            QMessageBox.critical(self, PLUGIN_NAME, "Network error. Please check your internet connection.")
            self.lbl_info.setText("Network error.")
        elif error_msg:
            QMessageBox.critical(self, PLUGIN_NAME, f"Error: {error_msg}")
            self.lbl_info.setText("Error occurred.")
        elif not results and 'response' in locals() and response.status_code != 200:
            self.lbl_info.setText("Not found in PubChem search.")
        else:
            # データを保持
            self.candidates_data = results
            self.update_table()
            
            if results:
                self.lbl_info.setText(f"Found {len(results)} candidates. Select one and click Load to 2D Editor.")
            else:
                self.lbl_info.setText("No results found.")

    def update_table(self):
        self.table.setRowCount(0)
        
        for i, data in enumerate(self.candidates_data):
            row_idx = self.table.rowCount()
            self.table.insertRow(row_idx)
            
            self.table.setItem(row_idx, 0, QTableWidgetItem(str(data['name'])))
            self.table.setItem(row_idx, 1, QTableWidgetItem(str(data['formula'])))
            self.table.setItem(row_idx, 2, QTableWidgetItem(str(data['smiles'])))

    def load_molecule(self):
        """選択された行のSMILESから2D構造を生成し、メインウィンドウのエディタに入れる"""
        selected_items = self.table.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, PLUGIN_NAME, "Please select a molecule from the list.")
            return
            
        row = selected_items[0].row()
        smiles = self.candidates_data[row]['smiles']
        name = self.candidates_data[row]['name']
        
        if not smiles:
             QMessageBox.warning(self, PLUGIN_NAME, "No SMILES data available for this entry.")
             return

        self.lbl_info.setText("Loading into 2D Editor...")
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        QApplication.processEvents()

        success = False
        error_msg = None

        try:
            # メインウィンドウのSMILES読み込み機能を使用
            if hasattr(self.main_window, "load_from_smiles"):
                self.main_window.load_from_smiles(smiles)
                success = True
            else:
                error_msg = "Main window does not support 'load_from_smiles'."

        except Exception as e:
            error_msg = str(e)
            
        finally:
            QApplication.restoreOverrideCursor()

        if success:
             self.lbl_info.setText(f"Loaded: {name}")
             QMessageBox.information(self, PLUGIN_NAME, f"Successfully loaded: {name}")
             self.accept() # ダイアログを閉じる
        elif error_msg:
             QMessageBox.critical(self, PLUGIN_NAME, f"Error: {error_msg}")
             self.lbl_info.setText("Load failed.")

def run(main_window):
    """
    プラグインのエントリーポイント
    """
    # 既存のダイアログがあればアクティブにする
    if hasattr(main_window, "_molecule_resolver_dialog") and main_window._molecule_resolver_dialog.isVisible():
        main_window._molecule_resolver_dialog.raise_()
        main_window._molecule_resolver_dialog.activateWindow()
        return

    dialog = MoleculeResolverDialog(main_window, parent=main_window)
    main_window._molecule_resolver_dialog = dialog
    dialog.show()