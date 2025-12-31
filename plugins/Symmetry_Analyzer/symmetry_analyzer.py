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
            QListWidget::item {
                padding: 1px;
                margin: 0px;
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
        groups_box = QGroupBox("1. Likely Point Groups (Select one)")
        groups_layout = QVBoxLayout()
        
        self.groups_list = QListWidget()
        self.groups_list.setSelectionMode(QListWidget.SelectionMode.SingleSelection)
        self.groups_list.itemSelectionChanged.connect(self.on_group_selected)
        groups_layout.addWidget(self.groups_list)
        groups_box.setLayout(groups_layout)
        splitter.addWidget(groups_box)

        # B. Operations List
        ops_box = QGroupBox("2. Symmetry Operations")
        ops_layout = QVBoxLayout()
        
        # Selected Group Display
        self.selected_group_label = QLabel("Symmetry Group: -")
        self.selected_group_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.selected_group_label.setStyleSheet("QLabel { font-size: 12pt; color: #2c3e50; margin: 2px; }")
        ops_layout.addWidget(self.selected_group_label)
        
        self.ops_list = QListWidget()
        self.ops_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection) # 複数選択可(Ctrl/Shift)
        self.ops_list.setAlternatingRowColors(True)
        self.ops_list.itemSelectionChanged.connect(self.on_op_selection_changed)
        ops_layout.addWidget(self.ops_list)
        ops_box.setLayout(ops_layout)
        splitter.addWidget(ops_box)
        
        # C. Matrix Details
        details_box = QGroupBox("3. Operation Details")
        details_layout = QVBoxLayout()
        
        self.op_details = QTextEdit()
        self.op_details.setReadOnly(True)
        self.op_details.setPlaceholderText("Select an operation above to view matrix details.")
        details_layout.addWidget(self.op_details)
        details_box.setLayout(details_layout)
        splitter.addWidget(details_box)

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
        self.selected_group_label.setText("Symmetry Group: -") # Reset to placeholder
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

        # 1つ目をデフォルトで選択 (Auto-select first item)
        if self.groups_list.count() > 0:
            self.groups_list.setCurrentRow(0)
            
        QMessageBox.information(self, "Done", 
            f"Found {len(self.group_data)} potential point groups.\n"
            "Sorted by strictness (smaller tolerance first).")

    def _get_op_sort_key(self, op):
        """
        対称操作のソート順を決定するキーを生成
        順序: Identity -> Rotation(Cn) -> Reflection(sigma) -> Inversion(i) -> Improper(Sn)
        """
        m = op.rotation_matrix
        det = np.linalg.det(m)
        trace = np.trace(m)
        tol = 1e-2

        # 1. Identity (Priority: 0)
        if np.allclose(m, np.eye(3), atol=tol):
            return (0, 0)

        # 2. Proper Rotation (Priority: 1)
        if np.isclose(det, 1.0, atol=tol):
            val = np.clip((trace - 1) / 2.0, -1.0, 1.0)
            angle = np.degrees(np.arccos(val))
            order = round(360.0 / angle) if angle > 1.0 else 1
            return (1, -order)

        # 3. Reflection (Priority: 2)
        if np.isclose(trace, 1.0, atol=tol):
            return (2, 0)

        # 4. Inversion (Priority: 3)
        if np.isclose(trace, -3.0, atol=tol):
            return (3, 0)
            
        # 5. Improper Rotation (Priority: 4)
        val = np.clip((trace + 1) / 2.0, -1.0, 1.0)
        angle = np.degrees(np.arccos(val))
        order = round(360.0 / angle) if angle > 1.0 else 1
        
        # S2 (Inversion) check for safety -> Priority 3
        if order == 2:
            return (3, 0)

        return (4, -order)

    def on_group_selected(self):
        """グループが選択されたら、そのオペレーションを表示"""
        item = self.groups_list.currentItem()
        if not item:
            return
            
        text = item.text()
        # "Td  (Range...)" から "Td" を取り出す
        sym = text.split()[0]
        
        # Display nicely formatted symbol
        s = self._format_symmetry_symbol(sym)
        if s.startswith("<html>"):
             # Keep HTML structure valid
             inner = s.replace("<html>", "").replace("</html>", "")
             self.selected_group_label.setText(f"<html>Symmetry Group: {inner}</html>")
        else:
             self.selected_group_label.setText(f"Symmetry Group: {s}")
        
        if sym in self.group_data:
            self.analyzer = self.group_data[sym]['analyzer']
            ops = self.analyzer.get_symmetry_operations()
            # クラス順にソート (Identity -> Rotation -> Reflection -> Inversion -> Improper)
            ops.sort(key=self._get_op_sort_key)
            self.symmetry_ops = ops
            
            self.update_ops_list()
            self.sym_btn.setEnabled(True)
            self.op_details.clear()
            
    def _format_symmetry_symbol(self, sym):
        """SchoenfliesシンボルをHTML形式に整形 (イタリック体 + 下付き文字)"""
        import re
        # C2v -> C, 2v
        # D3h -> D, 3h
        # Td  -> T, d
        # C*v -> C, *v -> C, ∞v
        
        match = re.match(r"^([A-Z]+)(.*)$", sym)
        if match:
            main = match.group(1)
            sub = match.group(2)
            
            # 無限の処理
            sub = sub.replace('*', '∞') 
            
            # 視認性を高めるためのフォント調整等はスタイルシートで行っていますが、
            # ここでは構造的なマークアップを提供します。
            return f"<html><i>{main}</i><sub>{sub}</sub></html>"
        return sym
            
    def update_ops_list(self):
        """リストウィジェットの更新"""
        self.ops_list.clear()
        self.op_details.clear()
        
        for i, op in enumerate(self.symmetry_ops):
            label = self._get_op_label(op, i)
            self.ops_list.addItem(label)

    def _get_op_label(self, op, i):
        """
        回転行列から操作へのラベルを生成 (Unicode表記)
        順序やロジックはソートキーと整合させています。
        Order: Identity -> Rotation -> Reflection -> Inversion -> Improper
        """
        m = op.rotation_matrix
        det = np.linalg.det(m)
        trace = np.trace(m)
        tol = 1e-2
        
        # Helper for subscripts (下付き文字生成)
        def to_sub(n):
            return str(n).translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))

        # 1. Identity (恒等操作)
        if np.allclose(m, np.eye(3), atol=tol):
            return f"#{i+1}: E (Identity)"
            
        # 2. Proper Rotation (回転操作)
        # Det = +1
        if np.isclose(det, 1.0, atol=tol):
            val = (trace - 1) / 2.0
            val = np.clip(val, -1.0, 1.0)
            angle = np.degrees(np.arccos(val))
            
            # 角度がほとんど0ならIdentity
            if angle < 1.0: 
                return f"#{i+1}: E (Identity)"
            
            order = round(360.0 / angle)
            return f"#{i+1}: C{to_sub(order)} (Rotation)"
            
        else:
            # Improper Operations (回映・鏡映・反転) 
            # Det = -1
            
            # 3. Reflection (鏡映: Trace = 1)
            if np.isclose(trace, 1.0, atol=tol):
                return f"#{i+1}: σ (Reflection)"

            # 4. Inversion (反転: Trace = -3)
            if np.isclose(trace, -3.0, atol=tol):
                return f"#{i+1}: i (Inversion)"
            
            # 5. Improper Rotation (回映)
            val = (trace + 1) / 2.0
            val = np.clip(val, -1.0, 1.0)
            angle = np.degrees(np.arccos(val))
            
            if angle < 1.0:
                return f"#{i+1}: σ (Reflection)"

            order = round(360.0 / angle)
            
            # S2 は Inversion と等価
            if order == 2:
                 return f"#{i+1}: i (Inversion)"
                 
            return f"#{i+1}: S{to_sub(order)} (Improper Rotation)"

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
        if np.isclose(trace, -3.0, atol=1e-2): 
            r = 0.25 
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
        """構造の対照化 (Symmetrization) - 一対一対応保証版"""
        if self.analyzer is None:
            return

        mol_pmg = self.get_pymatgen_molecule()
        if mol_pmg is None:
            return
            
        # 1. 重心補正
        original_coords = mol_pmg.cart_coords
        center_of_mass = np.mean(original_coords, axis=0)
        centered_coords = original_coords - center_of_mass
        
        species = mol_pmg.species
        ops = self.analyzer.get_symmetry_operations()
        
        if not ops:
            return

        # 割り当て問題（Kuhn-Munkres法）を解くためのライブラリを試行
        try:
            from scipy.optimize import linear_sum_assignment
            HAS_SCIPY = True
        except ImportError:
            HAS_SCIPY = False

        new_coords = np.zeros_like(centered_coords)
        n_ops = len(ops)
        n_atoms = len(centered_coords)
        mapping_error = False

        for op in ops:
            # 対称操作を適用した全座標を一度に計算
            rotated_coords = np.array([op.operate(p) for p in centered_coords])
            
            # 各原子 i (操作後) と 各原子 j (元) の距離行列を作成
            # ただし、元素種が異なる組み合わせには大きなペナルティ（無限大）を与える
            cost_matrix = np.zeros((n_atoms, n_atoms))
            for i in range(n_atoms):
                diff = centered_coords - rotated_coords[i]
                dists = np.linalg.norm(diff, axis=1)
                
                # 同一元素種でない場合は距離を無限大にする
                mask = [s != species[i] for s in species]
                dists[mask] = 1e9 
                cost_matrix[i] = dists

            # 一対一の最適マッピングを計算
            if HAS_SCIPY:
                # scipyがある場合はハンガリアン法で最適化
                row_ind, col_ind = linear_sum_assignment(cost_matrix)
            else:
                # scipyがない場合のフォールバック: 簡易的な最近接（重複排除付き）
                row_ind = np.arange(n_atoms)
                col_ind = np.full(n_atoms, -1, dtype=int)
                used_j = set()
                for i in range(n_atoms):
                    sorted_j = np.argsort(cost_matrix[i])
                    for j in sorted_j:
                        if j not in used_j and cost_matrix[i, j] < 1.5:
                            col_ind[i] = j
                            used_j.add(j)
                            break
                    if col_ind[i] == -1: # マッピング失敗
                        col_ind[i] = i # フォールバック
                        mapping_error = True

            # マッピング結果に基づき、逆操作で戻して加算
            for i, j in zip(row_ind, col_ind):
                if cost_matrix[i, j] > 1.5: mapping_error = True
                new_coords[i] += op.inverse.operate(centered_coords[j])

        # 全操作で平均化
        new_coords /= n_ops
        
        if mapping_error:
            QMessageBox.warning(self, "Warning", 
                "Some atoms could not be mapped cleanly.\n"
                "The structure might be too distorted for this point group.")

        # 2. 座標を元の重心位置に戻して反映
        final_coords = new_coords + center_of_mass
        self.update_rdkit_coords(final_coords)
        
        QMessageBox.information(self, "Symmetrized", 
            f"Structure symmetrized to {self.analyzer.sch_symbol}.\n"
            f"(Averaged over {len(ops)} operations with 1-to-1 mapping)")

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
        self.selected_group_label.clear()
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

