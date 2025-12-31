# --- Plugin Metadata ---
PLUGIN_NAME = "Symmetry Analyzer"
PLUGIN_VERSION = "2025.12.31"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Analyzes molecular symmetry (point group) and symmetrizes structures."

import sys
import numpy as np
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, 
    QDoubleSpinBox, QPushButton, QListWidget, 
    QTextEdit, QGroupBox, QMessageBox, QSplitter
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal

# --- RDKit Imports ---
try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
except ImportError:
    pass # MoleditPy環境なら入っているはず

# --- Pymatgen Imports ---
try:
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False

class SymmetryAnalysisWorker(QThread):
    """Background worker for symmetry analysis to prevent UI freezing"""
    finished = pyqtSignal(dict, bool) # processed_data, found_any
    
    def __init__(self, mol_pmg, max_tol):
        super().__init__()
        self.mol_pmg = mol_pmg
        self.max_tol = max_tol
        
    def run(self):
        import numpy as np
        tolerances = np.arange(0.1, self.max_tol + 0.05, 0.05)
        
        group_data = {}
        found_any = False
        
        for tol in tolerances:
            try:
                tol_val = float(tol)
                # Heavy calculation here
                analyzer = PointGroupAnalyzer(self.mol_pmg, tolerance=tol_val)
                sym = analyzer.sch_symbol
                
                if sym not in group_data:
                    group_data[sym] = {
                        'analyzer': analyzer,
                        'tols': [tol_val]
                    }
                else:
                    group_data[sym]['tols'].append(tol_val)
                found_any = True
                
            except Exception:
                pass 
        
        self.finished.emit(group_data, found_any)

class SymmetryAnalysisPlugin(QWidget):
    """
    MoleditPy Plugin: Molecular Symmetry Analyzer & Symmetrizer
    
    [機能]
    1. 分子の点群 (Point Group) の判定
    2. 対称操作 (回転・鏡映など) のリスト表示
    3. 許容誤差 (Tolerance) の調整
    4. 射影演算子法による構造の対照化 (Symmetrization)
    """
    
    def __init__(self, interface):
        """
        :param interface: MoleditPyのメインインターフェースオブジェクト
                          (get_molecule(), update_view() 等を持つと想定)
        """
        # 親ウィンドウを設定することで、メインウィンドウの前面に表示され、かつ他のアプリの後ろに行くようになる
        # interfaceがQWidget継承なら親にする
        parent = interface if hasattr(interface, 'windowTitle') else None 
        super().__init__(parent)
        self.setWindowFlags(Qt.WindowType.Window) # 独立したウィンドウとして振る舞う
        
        self.interface = interface
        self.analyzer = None     # pymatgen PointGroupAnalyzer instance
        self.symmetry_ops = []   # Detected symmetry operations
        self.worker = None       # QThread instance
        
        # 選択色を薄くするスタイルシート (Light Sky Blue)
        self.setStyleSheet("""
            QListWidget::item:selected {
                background-color: #87CEFA;
                color: black;
            }
        """)
        
        self.init_ui()

        if hasattr(self.interface, 'update_view'):
            self.interface.update_view()
        elif hasattr(self.interface, 'canvas'):
            self.interface.canvas.update() # 例

    def init_ui(self):
        main_layout = QVBoxLayout(self)

        # --- 1. Top Settings (Analyze & Symmetrize only) ---
        settings_group = QGroupBox("Actions")
        settings_layout = QGridLayout()
        
        # Max Tolerance Input
        tol_label = QLabel("Max Tol (Å):")
        self.max_tol_spin = QDoubleSpinBox()
        self.max_tol_spin.setRange(0.1, 10.0)
        self.max_tol_spin.setSingleStep(0.1)
        self.max_tol_spin.setValue(1.5) # Default
        
        settings_layout.addWidget(tol_label, 0, 0)
        settings_layout.addWidget(self.max_tol_spin, 0, 1)
        
        # Analyze Button
        self.calc_btn = QPushButton("Analyze (Scan)")
        self.calc_btn.setToolTip("Scan tolerances to find likely point groups.")
        self.calc_btn.clicked.connect(self.analyze_symmetry)

        # Symmetrize Button
        self.sym_btn = QPushButton("Symmetrize Detected")
        self.sym_btn.setToolTip("Symmetrize structure to match the selected point group.")
        self.sym_btn.clicked.connect(self.symmetrize_structure)
        self.sym_btn.setEnabled(False) 
        
        settings_layout.addWidget(self.calc_btn, 1, 0)
        settings_layout.addWidget(self.sym_btn, 1, 1)
        
        settings_group.setLayout(settings_layout)
        main_layout.addWidget(settings_group)

        # --- 2. Results Area (Splitter) ---
        splitter = QSplitter(Qt.Orientation.Vertical)
        
        # A. Likely Groups List
        groups_container = QWidget()
        groups_layout = QVBoxLayout(groups_container)
        groups_layout.setContentsMargins(0, 0, 0, 0)
        groups_layout.addWidget(QLabel("1. Likely Point Groups (Select one):"))
        
        self.groups_list = QListWidget()
        self.groups_list.setSelectionMode(QListWidget.SelectionMode.SingleSelection)
        self.groups_list.itemClicked.connect(self.on_group_selected)
        groups_layout.addWidget(self.groups_list)
        splitter.addWidget(groups_container)

        # B. Operations List
        ops_container = QWidget()
        ops_layout = QVBoxLayout(ops_container)
        ops_layout.setContentsMargins(0, 0, 0, 0)
        ops_layout.addWidget(QLabel("2. Symmetry Operations:"))
        
        self.ops_list = QListWidget()
        self.ops_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection) # 複数選択可(Ctrl/Shift)
        self.ops_list.setAlternatingRowColors(True)
        self.ops_list.itemSelectionChanged.connect(self.on_op_selection_changed)
        ops_layout.addWidget(self.ops_list)
        splitter.addWidget(ops_container)
        
        # C. Matrix Details
        details_container = QWidget()
        details_layout = QVBoxLayout(details_container)
        details_layout.setContentsMargins(0, 0, 0, 0)
        details_layout.addWidget(QLabel("3. Operation Details:"))
        
        self.op_details = QTextEdit()
        self.op_details.setReadOnly(True)
        self.op_details.setPlaceholderText("Select an operation above to view matrix details.")
        details_layout.addWidget(self.op_details)
        splitter.addWidget(details_container)

        # Set initial sizes for splitter (optional)
        splitter.setSizes([150, 200, 150])

        main_layout.addWidget(splitter)
        
        # Check Dependency
        if not HAS_PYMATGEN:
            self.calc_btn.setEnabled(False)
            QMessageBox.critical(self, "Dependency Error", 
                "This plugin requires 'pymatgen'.\nPlease install it via: pip install pymatgen")

        # Data storage
        self.group_data = {} # Symbol -> {'analyzer': obj, 'tols': [float]}

    def get_pymatgen_molecule(self):
        """MoleditPy(RDKit)の分子をpymatgen形式に変換"""
        # Manual says: mw.current_mol holds the RDKit object
        rd_mol = getattr(self.interface, 'current_mol', None)
        # Fallback for Mock or if attribute is missing
        if rd_mol is None and hasattr(self.interface, 'get_molecule'):
             rd_mol = self.interface.get_molecule()
             
        if rd_mol is None:
            return None

        try:
            conf = rd_mol.GetConformer()
        except ValueError:
            return None # Conformerがない場合

        species = [atom.GetAtomicNum() for atom in rd_mol.GetAtoms()]
        coords = [list(conf.GetAtomPosition(i)) for i in range(rd_mol.GetNumAtoms())]
        
        return Molecule(species, coords)

    def analyze_symmetry(self):
        """許容誤差を変えながら点群をスキャンし、UIにリストアップする"""
        mol_pmg = self.get_pymatgen_molecule()
        if mol_pmg is None:
            QMessageBox.warning(self, "Error", "No molecule to analyze.")
            return

        self.groups_list.clear()
        self.ops_list.clear()
        self.op_details.clear()
        self.sym_btn.setEnabled(False)
        self.group_data = {}
        
        # UI controls update
        self.calc_btn.setEnabled(False)
        self.calc_btn.setText("Scanning...")
        
        max_tol = self.max_tol_spin.value()
        
        # Start Worker Thread
        self.worker = SymmetryAnalysisWorker(mol_pmg, max_tol)
        self.worker.finished.connect(self.on_analysis_finished)
        self.worker.start()

    def on_analysis_finished(self, group_data, found_any):
        """スレッド完了後の処理"""
        self.calc_btn.setEnabled(True)
        self.calc_btn.setText("Analyze (Scan)")
        self.group_data = group_data
        
        if not found_any:
            self.groups_list.addItem("No point groups found.")
            return
            
        # ユーザー要望: Range(Tolerance)が小さい(=厳密に合致している)順に表示
        # ソートキー: (最小許容誤差, -最大許容誤差) -> 小さい誤差で見つかり、かつ範囲が広いものを優先
        sorted_keys = sorted(self.group_data.keys(), 
                             key=lambda k: (min(self.group_data[k]['tols']), -max(self.group_data[k]['tols'])))
        
        for sym in sorted_keys:
            data = self.group_data[sym]
            min_t = min(data['tols'])
            max_t = max(data['tols'])
            
            # リスト表示： "Td (Tol: 0.10 - 2.00)" のようにシンプルに
            item_text = f"{sym}  (Tol: {min_t:.2f} - {max_t:.2f} Å)"
            self.groups_list.addItem(item_text)
            
        QMessageBox.information(self, "Done", 
            f"Found {len(self.group_data)} potential point groups.\n"
            "Sorted by strictness (smaller tolerance first).")

    def _get_op_sort_key(self, op):
        """対称操作のソートキーを生成 (クラスごとにグループ化)"""
        m = op.rotation_matrix
        det = np.linalg.det(m)
        trace = np.trace(m)
        
        # Priority:
        # 1. Identity (0)
        # 2. Proper Rotation (1, order)
        # 3. Reflection (2)
        # 4. Inversion (3)
        # 5. Improper (4)

        if np.allclose(m, np.eye(3)):
            return (0, 0)
            
        if np.allclose(m, -np.eye(3)):
            return (3, 0)

        is_proper = np.isclose(det, 1.0)
        if is_proper:
             val = (trace - 1) / 2.0
             val = np.clip(val, -1.0, 1.0)
             angle = np.degrees(np.arccos(val))
             if angle < 1.0: order = 1
             else: order = int(round(360.0 / angle))
             return (1, order)
        else:
            if np.isclose(trace, 1.0):
                return (2, 0) # Reflection
            return (4, 0) # Improper

    def on_group_selected(self, item):
        """グループが選択されたら、そのオペレーションを表示"""
        text = item.text()
        # "Td  (Range...)" から "Td" を取り出す
        sym = text.split()[0]
        
        if sym in self.group_data:
            self.analyzer = self.group_data[sym]['analyzer']
            ops = self.analyzer.get_symmetry_operations()
            # クラス順にソート (Identity -> Rotation -> Reflection -> Inversion -> Improper)
            ops.sort(key=self._get_op_sort_key)
            self.symmetry_ops = ops
            
            self.update_ops_list()
            self.sym_btn.setEnabled(True)
            self.op_details.clear()
            
    def update_ops_list(self):
        """リストウィジェットの更新"""
        self.ops_list.clear()
        self.op_details.clear()
        
        for i, op in enumerate(self.symmetry_ops):
            label = self._get_op_label(op, i)
            self.ops_list.addItem(label)

    def _get_op_label(self, op, i):
        """回転行列から操作へのラベルを生成"""
        m = op.rotation_matrix
        det = np.linalg.det(m)
        trace = np.trace(m)
        
        if np.allclose(m, np.eye(3)):
            return f"#{i+1}: Identity (E)"
            
        if np.allclose(m, -np.eye(3)):
            return f"#{i+1}: Inversion (i)"
            
        if np.isclose(det, 1.0):
            # Proper rotation
            val = (trace - 1) / 2.0
            val = np.clip(val, -1.0, 1.0)
            angle = np.degrees(np.arccos(val))
            if angle > 1.0:
                order = round(360.0 / angle)
                return f"#{i+1}: Rotation (C{order})"
            else:
                 return f"#{i+1}: Rotation (Unknown)"
        else:
            # Improper rotation / Reflection
            # Reflection (sigma) usually has trace 1
            if np.isclose(trace, 1.0):
                return f"#{i+1}: Reflection (sigma)"
            
            # Improper axis Sn
            return f"#{i+1}: Improper Rotation (Sn)"

    def on_op_selection_changed(self):
        """操作の選択状態が変わったときの処理 (複数選択対応)"""
        items = self.ops_list.selectedItems()
        
        ops_to_show = []
        for item in items:
            row = self.ops_list.row(item)
            if 0 <= row < len(self.symmetry_ops):
                ops_to_show.append(self.symmetry_ops[row])
        
        # 3D可視化の更新
        self.visualize_ops(ops_to_show)
        
        # テキスト詳細の更新
        if len(ops_to_show) == 0:
            self.op_details.clear()
        elif len(ops_to_show) == 1:
            # 1つだけ選ばれているなら詳細を表示
            self._display_single_op_details(ops_to_show[0])
        else:
            self.op_details.setText(f"{len(ops_to_show)} operations selected.\n"
                                    "See 3D view for visualization.")

    def _display_single_op_details(self, op):
        """単一操作の詳細テキスト表示"""
        mat_str = np.array2string(op.rotation_matrix, precision=3, suppress_small=True)
        trans_str = np.array2string(op.translation_vector, precision=3, suppress_small=True)
        
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            xyz_str = op.as_xyz_str()
        
        text = (f"--- Rotation Matrix ---\n{mat_str}\n\n"
                f"--- Translation Vector ---\n{trans_str}\n\n"
                f"Type: {xyz_str}")
        self.op_details.setText(text)

    def visualize_ops(self, ops_list):
        """指定された複数の対称操作を3Dビューに描画する"""
        if not hasattr(self.interface, 'plotter'):
            return
            
        import numpy as np
        try:
            import pyvista as pv
        except ImportError:
            return

        plotter = self.interface.plotter
        
        # 以前の表示をクリア
        if not hasattr(self, 'vis_actors'):
            self.vis_actors = []
            
        for actor in self.vis_actors:
            plotter.remove_actor(actor)
        self.vis_actors = []

        if not ops_list:
            plotter.render()
            return

        # 分子の重心(COM)を計算 (共通)
        rd_mol = getattr(self.interface, 'current_mol', None)
        if rd_mol is None and hasattr(self.interface, 'get_molecule'):
            rd_mol = self.interface.get_molecule()
            
        if rd_mol:
            conf = rd_mol.GetConformer()
            positions = [conf.GetAtomPosition(i) for i in range(rd_mol.GetNumAtoms())]
            coords = np.array([[p.x, p.y, p.z] for p in positions])
            com = np.mean(coords, axis=0) # 重心
            # 分子の大きさ (最大半径) を計算
            dists = np.linalg.norm(coords - com, axis=1)
            mol_radius = np.max(dists) if len(dists) > 0 else 2.0
        else:
            com = np.array([0., 0., 0.])
            mol_radius = 2.0
            
        # 少し余裕を持たせる
        scale = mol_radius

        # 各操作を描画
        for op in ops_list:
            self._add_op_visualization(plotter, op, com, pv, scale)
            
        plotter.render()

    def _add_op_visualization(self, plotter, op, com, pv, scale=2.0):
        """単一の操作をシーンに追加 (内部ヘルパー)"""
        m = op.rotation_matrix
        det = np.linalg.det(m)
        trace = np.trace(m)
        
        # 1. Identity
        if np.allclose(m, np.eye(3)):
            return

        # 2. Inversion
        if np.allclose(m, -np.eye(3)):
            r = 0.25 # Fixed size (User request: half of 0.5)
            sphere = pv.Sphere(radius=r, center=com)
            actor = plotter.add_mesh(sphere, color="orange", opacity=0.8, name="sym_inversion")
            self.vis_actors.append(actor)
            return

        # 3. Rotation / Reflection
        eigvals, eigvecs = np.linalg.eig(m)
        axis_idx = np.where(np.isclose(eigvals, 1.0))[0]
        normal_idx = np.where(np.isclose(eigvals, -1.0))[0]
        
        # A. Proper Rotation
        if np.isclose(det, 1.0) and len(axis_idx) > 0:
            axis = np.real(eigvecs[:, axis_idx[0]])
            if np.linalg.norm(axis) < 1e-3: return 
            axis = axis / np.linalg.norm(axis)
            
            length = scale * 1.5
            start = com - axis * length
            end = com + axis * length
            line = pv.Line(start, end)
            actor_line = plotter.add_mesh(line, color="cyan", line_width=4)
            self.vis_actors.append(actor_line)

        # B. Mirror Reflection
        elif len(normal_idx) > 0 and np.isclose(trace, 1.0):
            normal = np.real(eigvecs[:, normal_idx[0]])
            if np.linalg.norm(normal) < 1e-3: return
            normal = normal / np.linalg.norm(normal)
            
            size = scale * 1.5
            disk = pv.Disc(center=com, inner=0, outer=size, normal=normal, c_res=30)
            actor_plane = plotter.add_mesh(disk, color="magenta", opacity=0.3)
            self.vis_actors.append(actor_plane)



    def symmetrize_structure(self):
        """構造の対照化 (Symmetrization)"""
        if self.analyzer is None:
            return

        mol_pmg = self.get_pymatgen_molecule()
        if mol_pmg is None:
            return
            
        coords = mol_pmg.cart_coords
        ops = self.analyzer.get_symmetry_operations()
        
        if not ops:
            return

        # 射影演算子法 (Permutation対応版)
        # 単純に op.operate(coords) を平均すると、原子が移動している場合(置換)
        # 全ての原子が重心に集まってしまう(Collapse)ため、
        # 「どの原子がどこに移ったか」を追跡して平均をとる必要がある。
        
        new_coords = np.zeros_like(coords)
        n_ops = len(ops)
        n_atoms = len(coords)
        
        for i in range(n_atoms):
            target_sum = np.zeros(3)
            current_pos = coords[i]
            
            for op in ops:
                # 1. 原子iを操作opで移動させた位置を計算
                rotated_pos = op.operate(current_pos)
                
                # 2. その位置に最も近い原子jを探す (Mapping/Permutation)
                #    歪があるため完全には一致しないが、最も近いものが対応する原子
                dists = np.linalg.norm(coords - rotated_pos, axis=1)
                j = np.argmin(dists)
                
                # 3. 対応する原子jの位置を、操作opの逆で引き戻した位置を加算
                #    「もし原子jが原子iの対称移動先なら、逆操作で戻せばiの理想位置になるはず」
                target_sum += op.inverse.operate(coords[j])
                
            new_coords[i] = target_sum / n_ops
        
        # RDKit側に反映
        self.update_rdkit_coords(new_coords)
        
        QMessageBox.information(self, "Symmetrized", 
            f"Structure symmetrized to {self.analyzer.sch_symbol}.\n"
            f"(Averaged over {len(ops)} operations)")

    def update_rdkit_coords(self, new_coords):
        """計算された座標をRDKitオブジェクトに戻し、ビューを更新"""
        rd_mol = getattr(self.interface, 'current_mol', None)
        if rd_mol is None and hasattr(self.interface, 'get_molecule'):
             rd_mol = self.interface.get_molecule()

        if rd_mol is None:
            return
            
        # Undo stateを保存 (Manual Section 4)
        if hasattr(self.interface, 'push_undo_state'):
            self.interface.push_undo_state()

        conf = rd_mol.GetConformer()
        for i in range(rd_mol.GetNumAtoms()):
            # floatキャスト (numpy.float64 は RDKit C++ API で弾かれることがあるため)
            x, y, z = map(float, new_coords[i])
            conf.SetAtomPosition(i, Point3D(x, y, z))
            
        # 3Dビューの更新 (Manual Section 4: mw.draw_molecule_3d(mol))
        if hasattr(self.interface, 'draw_molecule_3d'):
             self.interface.draw_molecule_3d(rd_mol)
             
        # その他の更新シグナル (Legacy support or generic)
        if hasattr(self.interface, 'update_view'):
            self.interface.update_view()
        elif hasattr(self.interface, 'canvas'):
            self.interface.canvas.update()

    def closeEvent(self, event):
        """ウィンドウが閉じられるときの処理 (クリーンアップ)"""
        # 1. 3D可視化の消去
        if hasattr(self.interface, 'plotter') and hasattr(self, 'vis_actors'):
            plotter = self.interface.plotter
            for actor in self.vis_actors:
                plotter.remove_actor(actor)
            self.vis_actors = []
            plotter.render()
            
        # 2. UIと内部データのリセット
        self.groups_list.clear()
        self.ops_list.clear()
        self.op_details.clear()
        self.sym_btn.setEnabled(False)
        self.group_data = {}
        self.symmetry_ops = []
        self.analyzer = None
        
        super().closeEvent(event)

def run(interface):
    """
    MoleditPy Plugin Entry Point
    :param interface: The main application window or interface object
    """
    # 既存のウィンドウがあれば閉じるなどの処理が必要かもしれませんが、
    # ここではシンプルに新しいウィンドウを作成して表示します。
    
    # interface (MainWindow) を親として渡すと、メインウィンドウと一緒に最小化/終了されます
    # ガベージコレクションされないように参照を保持します
    if not hasattr(interface, 'symmetry_plugin_window'):
        interface.symmetry_plugin_window = SymmetryAnalysisPlugin(interface)
        interface.symmetry_plugin_window.setWindowTitle("Symmetry Analyzer")
        interface.symmetry_plugin_window.resize(400, 600)
    
    interface.symmetry_plugin_window.show()
    interface.symmetry_plugin_window.raise_()

if __name__ == "__main__":
    from PyQt6.QtWidgets import QApplication
    
    # Mock Interface for testing
    class MockInterface:
        def __init__(self):
            self.mol = None
            self._create_sample_molecule()
            
        def _create_sample_molecule(self):
            try:
                # Create Methane (CH4) with 3D coordinates
                from rdkit import Chem
                from rdkit.Chem import AllChem
                
                m = Chem.MolFromSmiles('C')
                m = Chem.AddHs(m)
                AllChem.EmbedMolecule(m, randomSeed=42) # Generate 3D coords
                AllChem.MMFFOptimizeMolecule(m) # Optimize
                self.mol = m
                print("Mock: Created sample molecule (CH4) - Optimized.")
            except ImportError:
                print("Mock: RDKit not found or error creating molecule.")
                
        def get_molecule(self):
            return self.mol
            
        def update_view(self):
            print("Mock: View updated.")

    app = QApplication(sys.argv)
    
    # Check if dependencies are met for the mock
    try:
        import rdkit
        import pymatgen
        print("Dependencies found.")
    except ImportError as e:
        print(f"Warning: Missing dependencies for full functionality: {e}")

    interface = MockInterface()
    window = SymmetryAnalysisPlugin(interface)
    window.setWindowTitle("Symmetry Analyzer (Standalone Test)")
    window.resize(400, 600)
    window.show()
    
    sys.exit(app.exec())

