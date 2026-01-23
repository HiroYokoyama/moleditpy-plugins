
import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QListWidget, QListWidgetItem,
    QLabel, QComboBox, QGroupBox, QRadioButton, QDockWidget, QMessageBox,
    QColorDialog, QFormLayout, QTableWidget, QTableWidgetItem, QHeaderView,
    QFileDialog, QMenuBar, QCheckBox, QToolButton, QMenu, QProgressDialog
)
from PyQt6.QtCore import Qt, QTimer, QThread, pyqtSignal
from PyQt6.QtGui import QColor, QIcon, QAction
from rdkit import Chem
from rdkit.Chem import AllChem
import copy
import sys

class AlignmentWorker(QThread):
    progress = pyqtSignal(int, int, str) # current, total, status_message
    finished_signal = pyqtSignal(list)   # list of result dicts
    error_signal = pyqtSignal(str)

    def __init__(self, ref_mol, target_entries, method, ignore_hs):
        super().__init__()
        self.ref_mol_raw = ref_mol
        self.ref_work = copy.deepcopy(ref_mol)
        
        self.targets = []
        for idx, mol in target_entries:
             self.targets.append((idx, copy.deepcopy(mol)))
             
        self.method = method
        self.ignore_hs = ignore_hs
        
        # ★追加: 組み合わせ計算の上限回数（これを超えたら打ち切る）
        self.MAX_COMBINATIONS = 100000

    def run(self):
        results = []
        total_steps = len(self.targets)
        
        try:
            # Pre-calculation for Reference
            if self.ignore_hs and self.method == "Substructure (MCS)":
                ref_calc = Chem.RemoveHs(self.ref_work)
            else:
                ref_calc = self.ref_work

            for i, (orig_idx, probe_mol) in enumerate(self.targets):
                if self.isInterruptionRequested():
                    break
                    
                self.progress.emit(i, total_steps, f"Aligning Mol {orig_idx + 1}...")
                
                probe_work = copy.deepcopy(probe_mol)
                
                result_entry = {
                    'index': orig_idx,
                    'mol': probe_work,
                    'rms': -1.0
                }

                # Pre-calc for Probe
                if self.ignore_hs and self.method == "Substructure (MCS)":
                    probe_calc = Chem.RemoveHs(probe_work)
                else:
                    probe_calc = probe_work

                best_rms = float('inf')
                best_transform = None

                # --- Method A: Atom IDs ---
                if self.method == "Atom IDs":
                    if self.ignore_hs:
                        probe_heavy = [a.GetIdx() for a in probe_work.GetAtoms() if a.GetAtomicNum() != 1]
                        ref_heavy = [a.GetIdx() for a in self.ref_work.GetAtoms() if a.GetAtomicNum() != 1]
                        
                        if len(probe_heavy) == len(ref_heavy):
                            atom_map = list(zip(probe_heavy, ref_heavy))
                            try:
                                rms = AllChem.AlignMol(probe_work, self.ref_work, atomMap=atom_map, reflect=False)
                                best_rms = rms
                                result_entry['mol'] = probe_work
                            except RuntimeError:
                                pass
                    else:
                        if probe_work.GetNumAtoms() == self.ref_work.GetNumAtoms():
                            atom_map = [(k, k) for k in range(probe_work.GetNumAtoms())]
                            try:
                                rms = AllChem.AlignMol(probe_work, self.ref_work, atomMap=atom_map, reflect=False)
                                best_rms = rms
                                result_entry['mol'] = probe_work
                            except RuntimeError:
                                pass

                # --- Method B: MCS (修正箇所) ---
                elif self.method == "Substructure (MCS)":
                    from rdkit.Chem import rdFMCS
                    # タイムアウトを少し短く設定 (5秒は長い場合があるため適宜調整)
                    res = rdFMCS.FindMCS(
                        [ref_calc, probe_calc],
                        matchValences=True, ringMatchesRingOnly=True, completeRingsOnly=True, timeout=5
                    )
                    
                    if res.numAtoms > 0:
                        patt = Chem.MolFromSmarts(res.smartsString)
                        # uniquify=False は対称性を考慮するが、数が爆発する原因
                        ref_matches = ref_calc.GetSubstructMatches(patt, uniquify=False)
                        probe_matches = probe_calc.GetSubstructMatches(patt, uniquify=False)
                        
                        if ref_matches and probe_matches:
                            total_combos = len(ref_matches) * len(probe_matches)
                            
                            # ★ユーザーへの警告: 組み合わせが多すぎる場合
                            if total_combos > self.MAX_COMBINATIONS:
                                print(f"[WARNING] Mol {orig_idx}: Too many MCS matches ({total_combos} combinations).")
                                print(f"          Search will be limited to {self.MAX_COMBINATIONS} iterations.")
                                print(f"          Consider using 'Atom IDs' fitting method for better performance.")
                            
                            combo_count = 0
                            abort_search = False # 打ち切りフラグ
                            
                            for ref_match in ref_matches:
                                if abort_search: break 
                                
                                for probe_match in probe_matches:
                                    if self.isInterruptionRequested():
                                        abort_search = True
                                        break
                                    
                                    combo_count += 1
                                    
                                    # ★追加: 計算回数上限チェック
                                    if combo_count > self.MAX_COMBINATIONS:
                                        #print(f"Align stop: Exceeded max combinations ({self.MAX_COMBINATIONS}) for Mol {orig_idx}")
                                        abort_search = True
                                        break

                                    # 進捗更新（頻度を落として負荷軽減）
                                    if combo_count % 50 == 0:
                                         # Calculate percentage: if total < 2000, use total; otherwise cap at 2000
                                         max_value = min(total_combos, self.MAX_COMBINATIONS)
                                         progress_pct = int(combo_count / max_value * 100)
                                         self.progress.emit(i, total_steps, f"Aligning Mol {orig_idx+1}: {progress_pct}% ({combo_count}/{max_value})")

                                    atom_map = list(zip(probe_match, ref_match))
                                    probe_temp = copy.deepcopy(probe_calc)
                                    try:
                                        rms = AllChem.AlignMol(probe_temp, ref_calc, atomMap=atom_map, reflect=False)
                                    except RuntimeError:
                                        continue
                                        
                                    if rms < best_rms:
                                        best_rms = rms
                                        
                                        # ベストアラインメントが見つかった場合の変換適用ロジック
                                        if self.ignore_hs:
                                            p_heavy = [a.GetIdx() for a in probe_work.GetAtoms() if a.GetAtomicNum()!=1]
                                            r_heavy = [a.GetIdx() for a in self.ref_work.GetAtoms() if a.GetAtomicNum()!=1]
                                            full_map = [(p_heavy[p], r_heavy[r]) for p, r in atom_map]
                                            
                                            probe_final = copy.deepcopy(probe_work)
                                            try:
                                                AllChem.AlignMol(probe_final, self.ref_work, atomMap=full_map, reflect=False)
                                                best_transform = probe_final
                                            except:
                                                pass
                                        else:
                                            best_transform = probe_temp

                            if best_transform is not None:
                                result_entry['mol'] = best_transform
                                result_entry['rms'] = best_rms
                
                # Store valid result
                if best_rms != float('inf'):
                    result_entry['rms'] = best_rms
                
                results.append(result_entry)
            
            self.finished_signal.emit(results)
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.error_signal.emit(str(e))

PLUGIN_NAME = "Molecule Comparator"
PLUGIN_VERSION = "2026.01.23"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Compare multiple molecules in 3D, calculate RMSD, and align them."

# Default Palette
DEFAULT_COLORS = [
    "#FF0000", # Red
    "#00FF00", # Green
    "#0000FF", # Blue
    "#FFFF00", # Yellow
    "#00FFFF", # Cyan
    "#FF00FF", # Magenta
    "#FFA500", # Orange
    "#800080", # Purple
]

class MoleculeComparator(QWidget):
    def __init__(self, mw, ctrl=None):
        super().__init__(mw) # Parent to mw to ensure close on exit
        self.mw = mw
        self.ctrl = ctrl  # Plugin3DController for API access
        self.setWindowTitle("Molecule Comparator")
        self.molecules = [] # List of dicts: {'name': str, 'mol': Mol, 'color': str, 'scope': str, 'rms': float}
        
        self.setWindowFlags(Qt.WindowType.Window)
        self.setAcceptDrops(True)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # Molecule List
        header_layout = QHBoxLayout()
        header_layout.addWidget(QLabel("Loaded Molecules:"))
        self.btn_set_ref = QPushButton("Set as Ref")
        self.btn_set_ref.clicked.connect(self.set_as_reference)
        # self.btn_set_ref.setEnabled(False) # Enable/Disable based on selection?
        header_layout.addWidget(self.btn_set_ref)
        layout.addLayout(header_layout)
        
        self.mol_list = QListWidget()
        layout.addWidget(self.mol_list)

        # Controls
        btn_layout = QHBoxLayout()
        self.btn_add_current = QPushButton("Add Current")
        self.btn_add_current.clicked.connect(self.add_current_molecule)
        self.btn_load_file = QPushButton("Load File...")
        self.btn_load_file.clicked.connect(self.load_from_file)
        self.btn_remove = QPushButton("Remove")
        self.btn_remove.clicked.connect(self.remove_molecule)
        
        self.btn_redraw = QPushButton("Redraw")
        self.btn_redraw.clicked.connect(self.redraw_visualization)

        btn_layout.addWidget(self.btn_add_current)
        btn_layout.addWidget(self.btn_load_file)
        btn_layout.addWidget(self.btn_remove)
        btn_layout.addWidget(self.btn_redraw)
        layout.addLayout(btn_layout)

        # Global Visualization Settings
        style_layout = QHBoxLayout()
        style_layout.addWidget(QLabel("Global Style:"))
        self.combo_style = QComboBox()
        self.combo_style.addItems(["CPK", "Ball and Stick", "Sticks", "Wireframe"])
        self.combo_style.currentTextChanged.connect(self.change_style)
        self.combo_style.setCurrentIndex(2) # Default to Sticks
        style_layout.addWidget(self.combo_style)
        layout.addLayout(style_layout)
        
        # Wireframe Lighting Option
        self.check_wireframe_lighting = QCheckBox("Enable Lighting in Wireframe")
        self.check_wireframe_lighting.setChecked(False)
        self.check_wireframe_lighting.setToolTip("Enable or disable lighting when Wireframe style is selected")
        self.check_wireframe_lighting.stateChanged.connect(self.update_wireframe_lighting)
        layout.addWidget(self.check_wireframe_lighting)

        # Color Settings (Per Molecule)
        color_group = QGroupBox("Color Settings")
        color_layout = QFormLayout(color_group)
        
        self.btn_color = QPushButton("Pick Color...")
        self.btn_color.clicked.connect(self.change_color)
        self.combo_scope = QComboBox()
        self.combo_scope.addItems(["Carbon Only", "All Atoms"])
        self.combo_scope.currentIndexChanged.connect(self.change_scope)
        
        color_layout.addRow("Color:", self.btn_color)
        color_layout.addRow("Scope:", self.combo_scope)
        layout.addWidget(color_group)

        # Alignment Settings
        align_group = QGroupBox("Alignment & RMSD")
        align_layout = QFormLayout(align_group)
        
        self.combo_align_method = QComboBox()
        self.combo_align_method.addItems(["Substructure (MCS)", "Atom IDs"])
        self.combo_align_method.setCurrentIndex(0)  # Default to MCS
        
        self.check_ignore_hs = QCheckBox("Ignore Hydrogens")
        self.check_ignore_hs.setChecked(False) 
        self.check_ignore_hs.setToolTip("Exclude hydrogen atoms from RMSD calculation")
        
        self.btn_align = QPushButton("Align & Calculate RMSD")
        self.btn_align.clicked.connect(self.run_alignment)
        
        align_layout.addRow("Method:", self.combo_align_method)
        align_layout.addRow("", self.check_ignore_hs)
        align_layout.addRow(self.btn_align)
        layout.addWidget(align_group)
        
        # Results
        self.table_results = QTableWidget()
        self.table_results.setColumnCount(2)
        self.table_results.setHorizontalHeaderLabels(["Molecule", "RMSD (Å)"])
        self.table_results.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table_results.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)  # Make non-editable
        self.table_results.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows) # Select full row
        self.table_results.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.table_results.customContextMenuRequested.connect(self.show_results_context_menu)
        layout.addWidget(QLabel("Results:"))
        layout.addWidget(self.table_results)
        
        # Result action buttons
        result_btn_layout = QHBoxLayout()
        self.btn_copy_results = QPushButton("Copy Results")
        self.btn_copy_results.clicked.connect(self.copy_results_to_clipboard)
        self.btn_save_results = QPushButton("Save Results...")
        self.btn_save_results.clicked.connect(self.save_results_to_file)
        self.btn_export_png = QPushButton("Export as PNG...")
        self.btn_export_png.clicked.connect(self.export_as_png)
        result_btn_layout.addWidget(self.btn_copy_results)
        result_btn_layout.addWidget(self.btn_save_results)
        result_btn_layout.addWidget(self.btn_export_png)
        result_btn_layout.addStretch()
        layout.addLayout(result_btn_layout)

        # Close Button
        self.btn_close = QPushButton("Close Plugin")
        self.btn_close.clicked.connect(self.close_plugin)
        layout.addWidget(self.btn_close)

        # Events
        self.mol_list.currentRowChanged.connect(self.on_selection_changed)

    def show_results_context_menu(self, pos):
        item = self.table_results.itemAt(pos)
        if not item:
            return
            
        row = item.row()
        if row < 0 or row >= len(self.molecules):
            return

        from PyQt6.QtWidgets import QMenu
        menu = QMenu(self)
        
        action_mol = QAction("Export Coordinates as MOL...", self)
        action_mol.triggered.connect(lambda: self.export_molecule_coords(row, "mol"))
        menu.addAction(action_mol)
        
        action_xyz = QAction("Export Coordinates as XYZ...", self)
        action_xyz.triggered.connect(lambda: self.export_molecule_coords(row, "xyz"))
        menu.addAction(action_xyz)
        
        menu.exec(self.table_results.mapToGlobal(pos))

    def export_molecule_coords(self, row, format_type):
        if row < 0 or row >= len(self.molecules):
            return
            
        entry = self.molecules[row]
        mol = entry['mol']
        name = entry['name']
        
        default_ext = f".{format_type}"
        file_filter = f"{format_type.upper()} File (*{default_ext});;All Files (*)"
        
        default_name = f"{name}_aligned{default_ext}"
        
        file_path, _ = QFileDialog.getSaveFileName(
            self.mw, f"Export {name} as {format_type.upper()}", default_name, file_filter
        )
        
        if not file_path:
            return
            
        try:
            if format_type == "mol":
                print(Chem.MolToMolBlock(mol), file=open(file_path, 'w'))
            elif format_type == "xyz":
                print(Chem.MolToXYZBlock(mol), file=open(file_path, 'w'))
                
            QMessageBox.information(self.mw, "Exported", f"Saved to:\n{file_path}")
        except Exception as e:
            QMessageBox.critical(self.mw, "Error", f"Failed to export:\n{str(e)}")

    def closeEvent(self, event):
        self.cleanup_and_close()
        super().closeEvent(event)

    def close_plugin(self):
        self.cleanup_and_close()
        self.hide()

    def cleanup_and_close(self):
        # Restore original state
        self.mw._plugin_color_overrides = {}
        if self.mw.current_mol:
            self.mw.main_window_view_3d.draw_molecule_3d(self.mw.current_mol)
        else:
            self.mw.plotter.clear()
            
        self.exit_3d_only_mode()

    def exit_3d_only_mode(self):
         # Placeholder if not defined in snippet, but cleanup calls it
         pass

    def add_current_molecule(self):
        mol = self.mw.current_mol
        if not mol:
            QMessageBox.warning(self.mw, "Error", "No molecule loaded in Main Window.")
            return
            
        # Create a copy to preserve state
        mol_copy = Chem.Mol(mol)
        
        # Determine name
        name = ""
        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        
        if not name:
            name = f"Mol {len(self.molecules) + 1}"
            
        # Determine Color (cycle)
        color = DEFAULT_COLORS[len(self.molecules) % len(DEFAULT_COLORS)]
        
        entry = {
            'name': name,
            'mol': mol_copy,
            'color': color,
            'scope': "Carbon Only", # Default per user pref seems to be customization, start with Carbon
            'rms': None
        }
        
        self.molecules.append(entry)
        self.update_list()
        self.update_visualization()
        self.update_wireframe_lighting()
        self.reset_view()

    def reset_view(self):
        # Placeholder for view reset if needed
        pass

    def load_from_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self.mw, "Open Molecule File", "", 
            "Molecule Files (*.mol *.sdf *.pdb *.xyz);;All Files (*)"
        )
        
        if not file_path:
            return
            
        self.add_molecule_from_path(file_path)

    def add_molecule_from_path(self, file_path):
        try:
            ext = os.path.splitext(file_path)[1].lower()
            mol = None
            if ext in ['.mol', '.sdf']:
                mol = Chem.MolFromMolFile(file_path, removeHs=False)
            elif ext == '.pdb':
                mol = Chem.MolFromPDBFile(file_path, removeHs=False)
            
            if not mol:
                QMessageBox.warning(self.mw, "Error", f"Failed to load molecule from {file_path}")
                return

            # Sanitize to ensure proper 3D rendering properties
            try:
                Chem.SanitizeMol(mol)
            except:
                pass

            # Determine name
            name = os.path.basename(file_path)
            
            # Determine Color (cycle)
            color = DEFAULT_COLORS[len(self.molecules) % len(DEFAULT_COLORS)]
            
            entry = {
                'name': name,
                'mol': mol,
                'color': color,
                'scope': "Carbon Only", 
                'rms': None
            }
            
            self.molecules.append(entry)
            self.update_list()
            self.update_visualization()
            self.update_wireframe_lighting()
            self.reset_view()
            
        except Exception as e:
            QMessageBox.critical(self.mw, "Error", f"An error occurred loading the file:\n{str(e)}")

    def set_as_reference(self):
        row = self.mol_list.currentRow()
        if row <= 0:
            return
            
        # Move selected to top
        item = self.molecules.pop(row)
        self.molecules.insert(0, item)
        self.update_list()
        self.update_visualization()
        self.update_wireframe_lighting()
        self.reset_view()

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            if os.path.isfile(file_path):
                self.add_molecule_from_path(file_path)

    def remove_molecule(self):
        row = self.mol_list.currentRow()
        if row >= 0:
            self.molecules.pop(row)
            self.update_list()
            self.update_visualization()
            self.update_wireframe_lighting()
            self.reset_view()

    def update_list(self):
        self.mol_list.clear()
        for i, entry in enumerate(self.molecules):
            ref_mark = " (Ref)" if i == 0 else ""
            self.mol_list.addItem(f"{entry['name']}{ref_mark}")
        
        # Select last if added
        if self.molecules:
            self.mol_list.setCurrentRow(len(self.molecules) - 1)

    def on_selection_changed(self, row):
        if row < 0 or row >= len(self.molecules):
            self.btn_color.setEnabled(False)
            self.combo_scope.setEnabled(False)
            return
            
        self.btn_color.setEnabled(True)
        self.combo_scope.setEnabled(True)
        
        entry = self.molecules[row]
        # Update UI to reflect entry state
        # Set button style with automatic text color based on background brightness
        text_color = self._get_text_color_for_background(entry['color'])
        self.btn_color.setStyleSheet(f"background-color: {entry['color']}; color: {text_color};")
        
        scope_idx = self.combo_scope.findText(entry['scope'])
        if scope_idx >= 0:
            self.combo_scope.blockSignals(True)
            self.combo_scope.setCurrentIndex(scope_idx)
            self.combo_scope.blockSignals(False)

    def change_color(self):
        row = self.mol_list.currentRow()
        if row < 0: return
        
        entry = self.molecules[row]
        curr_color = QColor(entry['color'])
        new_color = QColorDialog.getColor(curr_color, self.mw, "Select Molecule Color")
        
        if new_color.isValid():
            entry['color'] = new_color.name()
            text_color = self._get_text_color_for_background(entry['color'])
            self.btn_color.setStyleSheet(f"background-color: {entry['color']}; color: {text_color};")
            self.update_visualization()
            self.update_wireframe_lighting()

    def _get_text_color_for_background(self, color_hex):
        """Calculate the best text color (black or white) for a given background color.
        
        Uses relative luminance formula to determine readability.
        """
        try:
            # Remove '#' if present
            hex_color = color_hex.lstrip('#')
            
            # Convert to RGB
            r = int(hex_color[0:2], 16) / 255.0
            g = int(hex_color[2:4], 16) / 255.0
            b = int(hex_color[4:6], 16) / 255.0
            
            # Calculate relative luminance (ITU-R BT.709)
            luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
            
            # Return white text for dark backgrounds, black for light backgrounds
            return 'white' if luminance < 0.5 else 'black'
        except:
            # Fallback to black if parsing fails
            return 'black'

    def redraw_visualization(self):
        """Redraw visualization and apply lighting settings."""
        self.update_visualization()
        self.update_wireframe_lighting()

    def change_scope(self):
        row = self.mol_list.currentRow()
        if row < 0: return
        
        entry = self.molecules[row]
        entry['scope'] = self.combo_scope.currentText()
        self.update_visualization()
        self.update_wireframe_lighting()

    def run_alignment(self):
        if len(self.molecules) < 2:
            return

        ref_entry = self.molecules[0]
        # Reset reference RMSD
        ref_entry['rms'] = 0.0

        # Prepare Data for Worker
        target_entries = [] # List of (index, mol)
        for i in range(1, len(self.molecules)):
            target_entries.append((i, self.molecules[i]['mol']))
        
        if not target_entries:
            return

        method = self.combo_align_method.currentText()
        ignore_hs = self.check_ignore_hs.isChecked()

        # Setup Progress Dialog
        self.progress_dialog = QProgressDialog("Initializing alignment...", "Stop", 0, len(target_entries), self.mw)
        self.progress_dialog.setWindowTitle("Alignment Progress")
        self.progress_dialog.setWindowModality(Qt.WindowModality.WindowModal)
        self.progress_dialog.setMinimumDuration(0) # Show immediately
        self.progress_dialog.setValue(0)

        # Setup Worker
        self.worker = AlignmentWorker(ref_entry['mol'], target_entries, method, ignore_hs)
        
        # Connect Signals
        self.worker.progress.connect(self.update_progress)
        self.worker.finished_signal.connect(self.on_alignment_finished)
        self.worker.error_signal.connect(self.on_alignment_error)
        self.progress_dialog.canceled.connect(self.worker.requestInterruption)
        
        # Start
        self.worker.start()

    def update_progress(self, current, total, message):
        if self.progress_dialog.wasCanceled():
            return
        self.progress_dialog.setMaximum(total)
        self.progress_dialog.setValue(current)
        self.progress_dialog.setLabelText(message)

    def on_alignment_finished(self, results):
        # Update Data
        for res in results:
            idx = res['index']
            if idx < len(self.molecules):
                self.molecules[idx]['mol'] = res['mol']
                self.molecules[idx]['rms'] = res['rms']
        
        # self.progress_dialog.display() # Ensure events processed
        self.progress_dialog.close()
        
        self.update_results_table()
        self.update_visualization()
        self.update_wireframe_lighting()
        
        # Determine success count
        success_count = sum(1 for r in results if r['rms'] != -1.0)
        if success_count < len(results):
             # Maybe show a small distinct message if some failed/stopped?
             pass

    def on_alignment_error(self, message):
        self.progress_dialog.close()
        QMessageBox.critical(self.mw, "Alignment Error", f"An error occurred during alignment:\n{message}")

    def update_results_table(self):
        self.table_results.setRowCount(len(self.molecules))
        for i, entry in enumerate(self.molecules):
            self.table_results.setItem(i, 0, QTableWidgetItem(entry['name']))
            rms_val = entry['rms']
            rms_str = f"{rms_val:.4f}" if rms_val is not None and rms_val >= 0 else ("N/A" if rms_val == -1.0 else "-")
            self.table_results.setItem(i, 1, QTableWidgetItem(rms_str))

    def copy_results_to_clipboard(self):
        """Copy results table to clipboard as tab-separated text."""
        if not self.molecules:
            return
        
        # Create tab-separated text
        lines = ["Molecule\tRMSD (Å)"]
        for entry in self.molecules:
            name = entry['name']
            rms_val = entry['rms']
            rms_str = f"{rms_val:.4f}" if rms_val is not None and rms_val >= 0 else ("N/A" if rms_val == -1.0 else "-")
            lines.append(f"{name}\t{rms_str}")
        
        text = "\n".join(lines)
        
        # Copy to clipboard
        from PyQt6.QtWidgets import QApplication
        clipboard = QApplication.clipboard()
        clipboard.setText(text)
        
        # Show confirmation
        self.mw.statusBar().showMessage("Results copied to clipboard", 3000)

    def save_results_to_file(self):
        """Save results table to a CSV file."""
        if not self.molecules:
            QMessageBox.warning(self.mw, "No Data", "No results to save.")
            return
        
        # Open file dialog
        file_path, _ = QFileDialog.getSaveFileName(
            self.mw, "Save Results", "", 
            "CSV Files (*.csv);;Text Files (*.txt);;All Files (*)"
        )
        
        if not file_path:
            return
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                # Write header
                f.write("Molecule,RMSD (Å)\n")
                
                # Write data
                for entry in self.molecules:
                    name = entry['name']
                    rms_val = entry['rms']
                    rms_str = f"{rms_val:.4f}" if rms_val is not None and rms_val >= 0 else ("N/A" if rms_val == -1.0 else "-")
                    # Escape commas in name if present
                    if ',' in name:
                        name = f'"{name}"'
                    f.write(f"{name},{rms_str}\n")
            
            # Show confirmation
            QMessageBox.information(self.mw, "Saved", f"Results saved to:\n{file_path}")
        except Exception as e:
            QMessageBox.critical(self.mw, "Error", f"Failed to save file:\n{str(e)}")

    def export_as_png(self):
        """Export the current 3D visualization as a PNG image."""
        if not self.molecules:
            QMessageBox.warning(self.mw, "No Data", "No molecules to export.")
            return
        
        # Create custom dialog for export options
        from PyQt6.QtWidgets import QDialog, QVBoxLayout, QCheckBox, QDialogButtonBox
        
        dialog = QDialog(self.mw)
        dialog.setWindowTitle("Export Settings")
        layout = QVBoxLayout(dialog)
        
        # Transparent background option
        check_transparent = QCheckBox("Transparent Background")
        check_transparent.setChecked(True)  # Default to transparent
        layout.addWidget(check_transparent)
        
        # Dialog buttons
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        
        # Show dialog
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        
        transparent = check_transparent.isChecked()
        
        # Open file dialog with default filename
        file_path, _ = QFileDialog.getSaveFileName(
            self.mw, "Export as PNG", "molecule_comparison.png", 
            "PNG Files (*.png);;All Files (*)"
        )
        
        if not file_path:
            return
        
        # Ensure .png extension
        if not file_path.lower().endswith('.png'):
            file_path += '.png'
        
        try:
            # Save screenshot using PyVista
            self.mw.plotter.screenshot(file_path, transparent_background=transparent)
            
            # Show confirmation
            QMessageBox.information(self.mw, "Exported", f"Image saved to:\n{file_path}")
        except Exception as e:
            QMessageBox.critical(self.mw, "Error", f"Failed to export image:\n{str(e)}")

    def update_visualization(self):
        # Combine all molecules
        if not self.molecules:
            # Clear colors using API if available
            if self.ctrl:
                # No direct "clear all" in API, so just draw None
                pass
            self.mw.main_window_view_3d.draw_molecule_3d(None)
            return

        combined_mol = Chem.Mol()
        
        current_atom_offset = 0
        
        for entry in self.molecules:
            mol = entry['mol']
            
            # Combine
            if combined_mol.GetNumAtoms() == 0:
                combined_mol = Chem.Mol(mol)
            else:
                combined_mol = Chem.CombineMols(combined_mol, mol)
            
            current_atom_offset += mol.GetNumAtoms()

        # Draw the combined molecule first
        self.mw.main_window_view_3d.draw_molecule_3d(combined_mol)
        
        # Then apply colors using the API
        if self.ctrl:
            current_atom_offset = 0
            for entry in self.molecules:
                mol = entry['mol']
                color_hex = entry['color']
                scope = entry['scope']
                
                for atom in mol.GetAtoms():
                    # Global index in the combined molecule
                    global_idx = current_atom_offset + atom.GetIdx()
                    
                    apply_color = False
                    if scope == "All Atoms":
                        apply_color = True
                    elif scope == "Carbon Only":
                        if atom.GetAtomicNum() == 6:
                            apply_color = True
                    
                    if apply_color:
                        try:
                            self.ctrl.set_atom_color(global_idx, color_hex)
                        except Exception as e:
                            print(f"Failed to set atom color: {e}")
                
                current_atom_offset += mol.GetNumAtoms()
            
            # Trigger redraw after colors are set
            self.mw.main_window_view_3d.draw_molecule_3d(combined_mol)
        else:
            # Fallback to direct access if controller not available (legacy)
            color_overrides = {}
            current_atom_offset = 0
            for entry in self.molecules:
                mol = entry['mol']
                color_hex = entry['color']
                scope = entry['scope']
                
                for atom in mol.GetAtoms():
                    global_idx = current_atom_offset + atom.GetIdx()
                    
                    apply_color = False
                    if scope == "All Atoms":
                        apply_color = True
                    elif scope == "Carbon Only":
                        if atom.GetAtomicNum() == 6:
                            apply_color = True
                    
                    if apply_color:
                        color_overrides[global_idx] = color_hex
                
                current_atom_offset += mol.GetNumAtoms()
            
            self.mw._plugin_color_overrides = color_overrides
            self.mw.main_window_view_3d.draw_molecule_3d(combined_mol)
        
    def _find_style_tool_button(self):
        """Helper to find the Main Window's 3D Style QToolButton."""
        try:
            # It's likely in a toolbar
            for btn in self.mw.findChildren(QToolButton):
                if "3D Style" in btn.text():
                    return btn
        except Exception:
            pass
        return None

    def _find_style_actions(self):
        """Helper to find style-related QActions in the Main Window."""
        actions = []
        try:
            # Search all actions in main window
            for action in self.mw.findChildren(QAction):
                # We check the text of the action
                txt = action.text().replace('&', '') 
                
                # Exclude specific settings dialog opener if we want to keep it enabled
                # But here we want to FIND it for potential disabling or checking?
                # User asked to keep "3D View Settings" enabled (toolbar only setting is wrong).
                # But "Style" actions inside the button should be identified.
                
                if ("Ball" in txt and "Stick" in txt) or ("CPK" in txt) or ("Wireframe" in txt) or ("Stick" in txt):
                    # Filter out unrelated actions if strictly necessary, but names are usually specific
                    actions.append(action)
        except Exception:
            pass
        return actions

    def change_style(self, style_name):
        # Update Main Window Style using explicit mapping
        mapping = {
            "CPK": "cpk",
            "Ball and Stick": "ball_and_stick",
            "Sticks": "stick",
            "Wireframe": "wireframe"
        }
        style_id = mapping.get(style_name, "ball_and_stick")
            
        # Update View 3D component
        if hasattr(self.mw, 'main_window_view_3d') and hasattr(self.mw.main_window_view_3d, 'current_3d_style'):
             self.mw.main_window_view_3d.current_3d_style = style_id
        
        # Also try direct mw attribute if exists (some versions might use mixin)
        if hasattr(self.mw, 'current_3d_style'):
            self.mw.current_3d_style = style_id

        # Push to settings
        if hasattr(self.mw, 'settings'):
             self.mw.settings['default_3d_style'] = style_id
        
        # Sync with Main Window Actions (Visual Checkmark)
        try:
            for action in self._find_style_actions():
                txt = action.text().replace('&', '')
                match = False
                if style_name == "CPK" and "CPK" in txt: match = True
                elif style_name == "Wireframe" and "Wireframe" in txt: match = True
                elif (style_name == "Sticks" or style_name == "Stick") and "Stick" in txt and "Ball" not in txt: match = True
                elif (style_name == "Ball and Stick" or style_name == "Ball & Stick") and "Ball" in txt and "Stick" in txt: match = True
                
                if match:
                    # Block signals to prevent recursion if the action triggers style change again
                    was_blocked = action.blockSignals(True)
                    action.setChecked(True)
                    action.blockSignals(was_blocked)
                else:
                    # Uncheck others if they are in an exclusive group (actions usually manage this, but being explicit helps)
                    # But be careful not to uncheck unrelated things.
                    # QActionGroup handles exclusivity automatically for these.
                    pass
        except Exception:
            pass

        # Ensure Main Window applies any side-effects of style change
        if hasattr(self.mw, 'apply_3d_settings'):
            try:
                self.mw.apply_3d_settings(redraw=False)
            except Exception:
                pass

        self.update_visualization()
        
        # Apply wireframe lighting setting AFTER visualization update
        # so it's not reset by redraw
        self.update_wireframe_lighting()

    def update_wireframe_lighting(self):
        """Control lighting based on current style and wireframe lighting checkbox."""
        try:
            current_style = self.combo_style.currentText()
            
            if hasattr(self.mw, 'plotter'):
                # If Wireframe mode, apply checkbox setting
                if current_style == "Wireframe":
                    enable_lighting = self.check_wireframe_lighting.isChecked()
                else:
                    # For other styles, always enable lighting
                    enable_lighting = True
                
                # Control lighting by setting actor properties
                # This is the correct way to enable/disable lighting in VTK/PyVista
                actors = self.mw.plotter.renderer.GetActors()
                if actors:
                    actors.InitTraversal()
                    for i in range(actors.GetNumberOfItems()):
                        actor = actors.GetNextItem()
                        if actor and hasattr(actor, 'GetProperty'):
                            prop = actor.GetProperty()
                            if prop:
                                # Enable/disable lighting for this actor
                                prop.SetLighting(enable_lighting)
                
                # Render to apply changes
                self.mw.plotter.render()
        except Exception as e:
            # Silently handle if plotter doesn't support lighting control
            pass

    def reset_view(self):
        # Use a timer to ensure the view is reset AFTER the visualization update is fully rendered/processed.
        # NOTE: mw.fit_to_view() is for 2D. We must use plotter methods for 3D.
        def _do_reset():
            if hasattr(self.mw, 'plotter'):
                try:
                    self.mw.plotter.reset_camera()
                    self.mw.plotter.render()
                except Exception:
                    pass
        
        QTimer.singleShot(100, _do_reset)

    def enter_3d_only_mode(self):
        # Save current splitter state if not already collapsed (approximately)
        if hasattr(self.mw, 'splitter'):
            current_sizes = self.mw.splitter.sizes()
            # Assuming [left, right] or similar. If > 0, left pane is visible.
            if len(current_sizes) >= 2 and current_sizes[0] > 0:
                self.saved_splitter_sizes = current_sizes
                # Collapse left (index 0)
                # We can try setting strict 0.
                self.mw.splitter.setSizes([0, 10000])

        # Disable 3D Drag (Edit Mode)
        if hasattr(self.mw, 'toggle_3d_edit_mode'):
            self.mw.toggle_3d_edit_mode(False)
        if hasattr(self.mw, 'edit_3d_action'):
            self.mw.edit_3d_action.setEnabled(False)

        # Disable Main Window Style Button (The Toolbar Button)
        btn = self._find_style_tool_button()
        if btn:
            btn.setEnabled(False)
            
        # Disable specific Style Actions (Just in case they are accessible elsewhere)
        for action in self._find_style_actions():
            action.setEnabled(False)

    def exit_3d_only_mode(self):
        # Restore splitter state if saved
        if hasattr(self, 'saved_splitter_sizes') and self.saved_splitter_sizes:
            if hasattr(self.mw, 'splitter'):
                self.mw.splitter.setSizes(self.saved_splitter_sizes)
            self.saved_splitter_sizes = None
            
        # Re-enable 3D Drag Action
        if hasattr(self.mw, 'edit_3d_action'):
            self.mw.edit_3d_action.setEnabled(True)

        # Re-enable Main Window Style Button
        btn = self._find_style_tool_button()
        if btn:
            btn.setEnabled(True)
            
        # Re-enable Style Actions
        for action in self._find_style_actions():
            action.setEnabled(True)

def run(mw, ctrl=None):
    # Check if already exists to prevent duplicates
    if not hasattr(mw, 'molecule_comparator_window'):
        # Pass mw and ctrl to ensure API access
        # Task: Fix RMSD Calculation Consistency
        win = MoleculeComparator(mw, ctrl)
        mw.molecule_comparator_window = win
    
    win = mw.molecule_comparator_window
    if win.isVisible():
        # Toggle behavior: if open, close it (which triggers clean up)
        win.close()
    else:
        win.show()
        win.raise_()
        win.enter_3d_only_mode()

def initialize(context):
    def run_plugin():
        mw = context.get_main_window()
        ctrl = context.get_3d_controller()
        run(mw, ctrl)
