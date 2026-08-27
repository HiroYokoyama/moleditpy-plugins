from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QMessageBox,
    QLabel,
    QHeaderView,
    QAbstractItemView,
    QComboBox,
    QCheckBox,
    QSpinBox,
    QProgressBar,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from rdkit import Chem
from rdkit.Chem import AllChem
import copy

PLUGIN_NAME = "Conformational Search"
PLUGIN_VERSION = "2026.08.27"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Perform conformational search using RDKit ETKDG."


class ConformerSearchWorker(QThread):
    """Embeds and optimizes conformers off the GUI thread, stoppable at any point."""

    progress = pyqtSignal(int, int, str)
    failed = pyqtSignal(str)
    completed = pyqtSignal(object, list, bool)

    # A single EmbedMultipleConfs() call cannot be interrupted, so embedding is
    # done in small batches with a stop check between them.
    EMBED_BATCH = 5

    def __init__(self, mol, num_confs, force_field, parent=None):
        super().__init__(parent)
        self.mol = mol
        self.num_confs = num_confs
        self.force_field = force_field
        self._stop = False

    def stop(self):
        self._stop = True

    def run(self):
        try:
            mol_calc = self.mol

            # ETKDG embeds from the graph, so it follows the chiral tags rather
            # than the coordinates on screen.  Re-perceive stereochemistry from
            # the displayed 3D geometry first: without this, a molecule whose
            # tags are unset yields a mix of both enantiomers, and one whose
            # tags are stale after a 3D edit yields conformers that are all the
            # mirror image of what the user is looking at.
            try:
                Chem.AssignStereochemistryFrom3D(mol_calc)
            except (ValueError, RuntimeError):
                pass

            params = AllChem.ETKDGv3()
            params.useSmallRingTorsions = True
            # Batches accumulate. clearConfs is not a keyword of the params
            # overload of EmbedMultipleConfs; it has to be set on params.
            params.clearConfs = False

            cids = []
            while len(cids) < self.num_confs and not self._stop:
                batch = min(self.EMBED_BATCH, self.num_confs - len(cids))
                new_cids = AllChem.EmbedMultipleConfs(
                    mol_calc, numConfs=batch, params=params
                )
                if not new_cids:
                    break
                cids.extend(new_cids)
                self.progress.emit(
                    0, self.num_confs, f"Embedding {len(cids)}/{self.num_confs}..."
                )

            if not cids:
                if self._stop:
                    self.completed.emit(mol_calc, [], True)
                else:
                    self.failed.emit("Failed to generate conformers.")
                return

            results = []
            total = len(cids)
            for i, cid in enumerate(cids):
                if self._stop:
                    break
                energy = self._optimize(mol_calc, cid)
                if energy is not None:
                    results.append((energy, cid))
                self.progress.emit(
                    i + 1, total, f"Optimizing {i + 1}/{total} ({self.force_field})..."
                )

            if not results and not self._stop:
                self.failed.emit(f"Optimization failed with {self.force_field}.")
                return

            results.sort(key=lambda x: x[0])
            self.completed.emit(mol_calc, results, self._stop)

        except Exception as e:
            self.failed.emit(f"Error during search: {str(e)}")

    def _optimize(self, mol_calc, cid):
        """Optimize one conformer in place; return its energy, or None on failure."""
        if self.force_field == "MMFF94":
            if AllChem.MMFFOptimizeMolecule(mol_calc, confId=cid) == -1:
                return None
            prop = AllChem.MMFFGetMoleculeProperties(mol_calc)
            if not prop:
                return None
            ff = AllChem.MMFFGetMoleculeForceField(mol_calc, prop, confId=cid)
        else:
            if AllChem.UFFOptimizeMolecule(mol_calc, confId=cid) == -1:
                return None
            ff = AllChem.UFFGetMoleculeForceField(mol_calc, confId=cid)
        if not ff:
            return None
        # CalcEnergy() reads the snapshot taken at the first call; rebind first.
        ff.Initialize()
        return ff.CalcEnergy()


class ConformerSearchDialog(QDialog):
    def __init__(self, context, parent=None):
        super().__init__(parent)
        self.context = context
        self.main_window = context.get_main_window()
        self.setWindowTitle("Conformational Search & Preview")
        self.resize(400, 560)

        # メインウィンドウの分子への参照
        self.target_mol = context.current_mol

        # 計算用の一時的な分子（オリジナルを汚染しないため）
        self.temp_mol = None
        # 生成された配座データのリスト [(Energy, ConformerID), ...]
        self.conformer_data = []
        # 全ての計算結果（未フィルタ）
        self.results_raw = []
        self.worker = None

        # Original coordinates for restoration on cancel
        self.original_coords = []
        if self.target_mol:
            conf = self.target_mol.GetConformer()
            self.original_coords = [
                conf.GetAtomPosition(i) for i in range(self.target_mol.GetNumAtoms())
            ]

        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # 説明ラベル
        self.lbl_info = QLabel(
            "Click 'Run Search' to generate conformers.\nSelect a row to preview."
        )
        layout.addWidget(self.lbl_info)

        # Force Field Selection
        hbox_ff = QHBoxLayout()
        hbox_ff.addWidget(QLabel("Force Field:"))
        self.combo_ff = QComboBox()
        self.combo_ff.addItems(["MMFF94", "UFF"])
        hbox_ff.addWidget(self.combo_ff)
        hbox_ff.addSpacing(12)
        hbox_ff.addWidget(QLabel("Conformers:"))
        self.spin_confs = QSpinBox()
        self.spin_confs.setRange(1, 1000)
        self.spin_confs.setValue(30)
        hbox_ff.addWidget(self.spin_confs)
        hbox_ff.addStretch()
        layout.addLayout(hbox_ff)

        # Show All Checkbox (Moving under FF box)
        self.cb_show_all = QCheckBox("Show All Conformers (ignore energy redundancy)")
        self.cb_show_all.setChecked(False)
        self.cb_show_all.toggled.connect(self.apply_filter_and_update)
        layout.addWidget(self.cb_show_all)

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

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
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.currentItemChanged.connect(self.preview_conformer)
        layout.addWidget(self.table)

        # ボタンエリア
        btn_layout = QHBoxLayout()
        self.btn_run = QPushButton("Run Search")
        self.btn_run.clicked.connect(self.run_search)

        self.btn_stop = QPushButton("Stop")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.stop_search)

        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(
            self.accept
        )  # 閉じる（現在のプレビュー状態で確定）

        btn_layout.addWidget(self.btn_run)
        btn_layout.addWidget(self.btn_stop)
        btn_layout.addWidget(self.btn_close)
        layout.addLayout(btn_layout)

    def accept(self):
        self._shutdown_worker()
        # Commit current conformer state to the main window
        if self.target_mol:
            self.context.push_undo_checkpoint()
        super().accept()
        # Unregister so next open creates a fresh dialog with clean state
        self.context.register_window("main_panel", None)

    def reject(self):
        self._shutdown_worker()
        # Restore original coordinates if user cancels/closes without 'Accept'
        if self.target_mol and self.original_coords:
            conf = self.target_mol.GetConformer()
            for i, pos in enumerate(self.original_coords):
                conf.SetAtomPosition(i, pos)
            self.context.refresh_3d_view()
        super().reject()
        # Unregister so next open creates a fresh dialog with clean state
        self.context.register_window("main_panel", None)

    def closeEvent(self, event):
        # X button: behave like the Close button (accept), not cancel
        self.accept()
        event.ignore()

    def _shutdown_worker(self):
        """Stop a running search and wait for the thread before tearing down."""
        worker = self.worker
        if worker is None:
            return
        worker.stop()
        worker.wait()
        self.worker = None

    def _set_running(self, running):
        self.btn_run.setEnabled(not running)
        self.btn_stop.setEnabled(running)
        self.combo_ff.setEnabled(not running)
        self.spin_confs.setEnabled(not running)
        self.progress.setVisible(running)
        if running:
            self.progress.setValue(0)

    def run_search(self):
        if self.worker is not None:
            return
        # Always pick up the molecule currently loaded in the main window
        current_mol = self.context.current_mol
        if not current_mol:
            QMessageBox.warning(self, PLUGIN_NAME, "No molecule loaded.")
            return
        # If the molecule changed since last run, refresh target and original_coords
        if current_mol is not self.target_mol:
            self.target_mol = current_mol
            conf = self.target_mol.GetConformer()
            self.original_coords = [
                conf.GetAtomPosition(i) for i in range(self.target_mol.GetNumAtoms())
            ]

        self._set_running(True)
        self.lbl_info.setText("Running conformational search...")

        # 計算用に分子を複製（水素が付加されていることを推奨）
        # Copied on the GUI thread so the worker never touches the molecule the
        # main window is drawing.
        mol_calc = copy.deepcopy(self.target_mol)

        self.worker = ConformerSearchWorker(
            mol_calc, self.spin_confs.value(), self.combo_ff.currentText(), parent=self
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.failed.connect(self.on_failed)
        self.worker.completed.connect(self.on_completed)
        self.worker.start()

    def stop_search(self):
        if self.worker is None:
            return
        self.worker.stop()
        self.btn_stop.setEnabled(False)
        self.lbl_info.setText("Stopping... (keeping conformers found so far)")

    def on_progress(self, done, total, message):
        if total > 0:
            self.progress.setValue(int(done * 100 / total))
        self.lbl_info.setText(message)

    def on_failed(self, message):
        self._set_running(False)
        self.worker = None
        QMessageBox.warning(self, PLUGIN_NAME, message)
        self.lbl_info.setText("Failed.")

    def on_completed(self, mol_calc, results, was_stopped):
        self._set_running(False)
        self.worker = None
        if not results:
            self.lbl_info.setText("Stopped before any conformer was optimized.")
            return
        self.temp_mol = mol_calc
        self.results_raw = results
        self.apply_filter_and_update()
        if was_stopped:
            self.lbl_info.setText(f"Stopped. {self.lbl_info.text()}")

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
        self.lbl_info.setText(
            f"Showing {len(self.conformer_data)} conformers (Total found: {len(self.results_raw)})."
        )

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
        target_conf = self.target_mol.GetConformer()  # 現在の表示用Conformer

        # Safety check: Atom count must match
        if self.temp_mol.GetNumAtoms() != self.target_mol.GetNumAtoms():
            self.lbl_info.setText(
                "<font color='red'>Error: Molecule changed in main window. Restart search.</font>"
            )
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
    if hasattr(mw, "host"):
        mw = mw.host
    if _launch_fn:
        _launch_fn()
