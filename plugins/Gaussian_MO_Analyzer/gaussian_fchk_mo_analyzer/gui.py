import os
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton,
                             QListWidget, QLabel, QDoubleSpinBox, QProgressBar, QSpinBox, 
                             QGroupBox, QMessageBox, QWidget, QAbstractItemView, QListWidgetItem, QFormLayout, QCheckBox)
from PyQt6.QtGui import QColor, QBrush
from PyQt6.QtCore import Qt, QThread, pyqtSignal

# Relative imports
from .analyzer import FCHKReader, BasisSetEngine
from .vis import CubeVisualizer
from .worker import CalcWorker

class OrbitalWidget(QWidget):
    def __init__(self, parent, context, fchk_path):
        super().__init__(parent)
        # Make this a Tool window so it stays on top of parent but doesn't block (Modeless)
        # and doesn't sit on top of other apps (unlike WindowStaysOnTopHint)
        # self.setWindowFlags(Qt.WindowType.Tool)
        
        self.context = context
        self.fchk_path = fchk_path
        self.setWindowTitle(f"Gaussian MO Analyzer - {os.path.basename(fchk_path)}")
        self.resize(500, 600)
        
        # Init Backend
        self.reader = FCHKReader(fchk_path)
        self.engine = BasisSetEngine(self.reader)
        
        # Load Molecule to Main Window
        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds
            
            xyz = self.reader.get_xyz_block()
            if xyz:
                raw_mol = Chem.MolFromXYZBlock(xyz)
                if raw_mol:
                    mol = Chem.RWMol(raw_mol)
                    
                    # Get charge from FCHK for better bond order guess
                    charge = 0
                    try:
                        c_arr = self.reader.get("Charge")
                        if c_arr is not None:
                            charge = int(c_arr[0])
                    except:
                        pass
                        
                    # Connectivity and Bond Order
                    try:
                        rdDetermineBonds.DetermineConnectivity(mol)
                        rdDetermineBonds.DetermineBondOrders(mol, charge=charge)
                    except Exception as e:
                        print(f"rdDetermineBonds failed: {e}")
                    
                    context.current_molecule = mol
                    
                    mw = context.get_main_window()
                    
                    if hasattr(mw, '_enter_3d_viewer_ui_mode'):
                        mw._enter_3d_viewer_ui_mode()
                        
                    mw.plotter.reset_camera()
        except Exception as e:
            print(f"Failed to load molecule structure: {e}")
        
        # Layout
        layout = QVBoxLayout(self)
        
        # Settings
        grp = QGroupBox("Generation Settings")
        form = QFormLayout(grp)
        
        self.spin_points = QSpinBox()
        self.spin_points.setRange(10, 500)
        self.spin_points.setValue(40)
        self.spin_points.setSingleStep(5)
        form.addRow("Grid Points (x,y,z):", self.spin_points)
        
        self.spin_gap = QDoubleSpinBox()
        self.spin_gap.setRange(1.0, 20.0) 
        self.spin_gap.setValue(4.0) 
        form.addRow("Margin (Bohr):", self.spin_gap)
        
        layout.addWidget(grp)
        
        # Isovalue Control
        iso_layout = QHBoxLayout()
        iso_layout.addWidget(QLabel("Isovalue:"))
        self.spin_iso = QDoubleSpinBox()
        self.spin_iso.setRange(0.001, 0.5)
        self.spin_iso.setSingleStep(0.005)
        self.spin_iso.setDecimals(3)
        self.spin_iso.setValue(0.02)
        self.spin_iso.valueChanged.connect(self.on_iso_changed)
        iso_layout.addWidget(self.spin_iso)
        layout.addLayout(iso_layout)

        # Opacity Control
        op_layout = QHBoxLayout()
        op_layout.addWidget(QLabel("Opacity:"))
        self.spin_opacity = QDoubleSpinBox()
        self.spin_opacity.setRange(0.0, 1.0)
        self.spin_opacity.setSingleStep(0.1)
        self.spin_opacity.setDecimals(1)
        self.spin_opacity.setValue(0.4)
        self.spin_opacity.valueChanged.connect(self.on_iso_changed)
        op_layout.addWidget(self.spin_opacity)
        layout.addLayout(op_layout)
        
        # Shading
        self.check_smooth = QCheckBox("Smooth Shading")
        self.check_smooth.setChecked(True)
        self.check_smooth.toggled.connect(self.on_iso_changed)
        layout.addWidget(self.check_smooth)
        
        # Mode Selection (Disabled/Commented out per user request)
        
    # List
        self.list_widget = QListWidget()
        self.list_widget.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.list_widget.itemDoubleClicked.connect(self.on_double_click)
        self.list_widget.itemClicked.connect(self.on_single_click)
        layout.addWidget(self.list_widget)
        
        # Populate List
        self.populate_list()
        
        # Progress
        self.pbar = QProgressBar()
        self.pbar.setVisible(False)
        layout.addWidget(self.pbar)
        
        # Buttons
        btn_layout = QHBoxLayout()
        self.btn_gen = QPushButton("Generate & Visualize")
        self.btn_gen.clicked.connect(self.generate_selected)
        btn_layout.addWidget(self.btn_gen)
        layout.addLayout(btn_layout)
        
        self.worker = None
        self.last_cube_path = None
        self.generation_queue = []

    def on_iso_changed(self, val):
        if self.last_cube_path and os.path.exists(self.last_cube_path):
            self.visualize(self.last_cube_path)

    def get_cube_path(self, mo_idx):
        # mo_idx is 1-based index
        base_dir = os.path.dirname(self.fchk_path)
        base_name = os.path.splitext(os.path.basename(self.fchk_path))[0]
        out_dir = os.path.join(base_dir, f"{base_name}_cubes")
        # Ensure dir exists? Maybe not here, only when writing.
        # But for checking existence, we need the path.
        out_name = f"{base_name}_MO{mo_idx}.cube"
        return os.path.join(out_dir, out_name)

    def populate_list(self):
        self.list_widget.clear()
        # mode = self.combo_mode.currentIndex() # 0=MO, 1=AO
        mode = 0 # Force MO
        
        if mode == 0:
            # MO Mode
            n_basis = self.engine.n_basis
            coeffs = self.reader.get("Alpha MO coefficients")
            if not n_basis or not len(coeffs):
                self.list_widget.addItem("Error reading MOs")
                return
                
            n_mos = len(coeffs) // n_basis
            energies = self.reader.get("Alpha Orbital Energies")
            n_alpha = self.reader.get("Number of alpha electrons")
            if n_alpha is not None and len(n_alpha) > 0:
                n_occ = n_alpha[0]
            else:
                n_el_list = self.reader.get("Number of electrons")
                n_el = n_el_list[0] if n_el_list is not None else 0
                n_occ = int(np.ceil(n_el / 2))
            
            # Populate in Reverse Order (High Energy Top)
            for i in range(n_mos - 1, -1, -1):
                mo_idx = i + 1
                label = f"MO {mo_idx:3d}"
                if energies is not None and i < len(energies):
                    label += f"   E={energies[i]:.5f} Ha"
                
                # HOMO/LUMO Labels
                homo_idx = n_occ - 1
                diff = i - homo_idx
                
                start_marker = "  <-- "
                
                if diff == 0:
                    label += f"{start_marker}HOMO"
                elif diff == 1:
                    label += f"{start_marker}LUMO"
                elif diff < 0 and diff >= -99:
                    label += f"{start_marker}HOMO{diff}"
                elif diff > 1 and diff <= 99:
                    label += f"{start_marker}LUMO+{diff-1}" # LUMO+1 is diff=2
                
                item = QListWidgetItem(label)
                # Store actual MO Index (1-based)
                item.setData(Qt.ItemDataRole.UserRole, mo_idx)
                
                # Check existence
                path = self.get_cube_path(mo_idx)
                if os.path.exists(path):
                    item.setForeground(QBrush(QColor("darkgreen")))
                
                self.list_widget.addItem(item)
                
            # Scroll to HOMO (or slightly above for context)
            # Find the item with HOMO
            for r in range(self.list_widget.count()):
                item = self.list_widget.item(r)
                if "HOMO" in item.text() and "-" not in item.text().split("HOMO")[-1]: # Exact HOMO match or check UserRole
                     idx = item.data(Qt.ItemDataRole.UserRole)
                     if idx == n_occ: # Verify
                         self.list_widget.setCurrentRow(r)
                         self.list_widget.scrollToItem(item, QAbstractItemView.ScrollHint.PositionAtCenter)
                         break

    def on_single_click(self, item):
        mo_idx = item.data(Qt.ItemDataRole.UserRole)
        path = self.get_cube_path(mo_idx)
        if os.path.exists(path):
            self.last_cube_path = path # Update current reference
            self.visualize(path)

    def on_double_click(self, item):
        self.generate_selected()

    def generate_selected(self):
        items = self.list_widget.selectedItems()
        if not items: return
        
        self.generation_queue = []
        for item in items:
            mo_idx = item.data(Qt.ItemDataRole.UserRole)
            self.generation_queue.append(mo_idx)
            
        if not self.generation_queue: return
        
        # Lock UI
        self.btn_gen.setEnabled(False)
        self.list_widget.setEnabled(False)
        self.pbar.setValue(0)
        self.pbar.setVisible(True)
        
        self.process_generation_queue()

    def process_generation_queue(self):
        if not self.generation_queue:
            # All done
            self.pbar.setVisible(False)
            self.btn_gen.setEnabled(True)
            self.list_widget.setEnabled(True)
            return

        mo_idx = self.generation_queue[0] # Peek/Process first. Don't pop yet if we want to use worker's mo_idx
        # Actually it's cleaner to pop now.
        self.current_gen_mo_idx = self.generation_queue.pop(0)
        
        # Setup Path
        path = self.get_cube_path(self.current_gen_mo_idx)
        self.last_cube_path = path

        # Check Existence
        if os.path.exists(path):
            # Already exists, just skip calculation but trigger visualization of last one?
            # Or visualize all? Standard behavior: Visualize last one or all?
            # Let's visualize and move to next
            self.update_list_item(self.current_gen_mo_idx)
            # If it's the last one (or single selection), visualize it
            # But continuous visualization might be heavy.
            # Let's visualize only if queue is empty or one by one?
            # User might want to see them pop up.
            self.visualize(path)
            self.process_generation_queue()
            return
            
        # Calc
        mode_str = "MO"
        row_idx = self.current_gen_mo_idx - 1 # 0-based for engine
        
        atoms = self.reader.get("Current cartesian coordinates").reshape(-1, 3)
        coeffs = self.reader.get("Alpha MO coefficients")
        
        self.pbar.setValue(0)
        
        self.worker = CalcWorker(
            self.engine, row_idx, 
            self.spin_points.value(), 
            self.spin_gap.value(),
            atoms, coeffs, path,
            mode=mode_str
        )
        self.worker.progress_sig.connect(self.pbar.setValue)
        self.worker.finished_sig.connect(self.on_finished)
        self.worker.start()

    def update_list_item(self, mo_idx):
        # Update color
        for r in range(self.list_widget.count()):
            item = self.list_widget.item(r)
            if item.data(Qt.ItemDataRole.UserRole) == mo_idx:
                item.setForeground(QBrush(QColor("darkgreen")))
                break

    def on_finished(self, success, result):
        if success:
            # Mark as generated
            self.update_list_item(self.current_gen_mo_idx)
            
            # Visualize
            self.visualize(result)
        else:
            QMessageBox.warning(self, "Error", f"Calculation failed for MO {self.current_gen_mo_idx}: {result}")
            # Continue queue even if fail? Yes.
        
        # Next
        self.process_generation_queue()

    def visualize(self, cube_path):
        mw = self.context.get_main_window()
        
        # Get Isovalue
        val = self.spin_iso.value()
        op_val = self.spin_opacity.value()
        is_smooth = self.check_smooth.isChecked()

        # Use our visualizer
        vis = CubeVisualizer(mw)
        if vis.load_file(cube_path):
            vis.show_iso(isovalue=val, opacity=op_val, smooth_shading=is_smooth)
            # Only reset camera if it's the first time or user requests it?
            # Creating a new actor usually keeps camera if not reset.
            # But let's reset to ensure visibility.
            # mw.plotter.reset_camera() 
            mw.plotter.render()
            mw.statusBar().showMessage(f"Visualizing {os.path.basename(cube_path)} (Iso={val}, Opacity={op_val}, Smooth={is_smooth})")
        else:
            QMessageBox.warning(self, "Error", "Failed to load generated Cube file.")
