import os
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton,
                             QListWidget, QLabel, QDoubleSpinBox, QProgressBar, QSpinBox, 
                             QGroupBox, QMessageBox, QWidget, QAbstractItemView, QListWidgetItem, QFormLayout)
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
        
        # Mode Selection (Disabled/Commented out per user request)
        # ...
        
        # List
        self.list_widget = QListWidget()
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
        self.btn_gen.clicked.connect(self.generate_current)
        btn_layout.addWidget(self.btn_gen)
        layout.addLayout(btn_layout)
        
        self.worker = None
        self.last_cube_path = None

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
            
            for i in range(n_mos):
                mo_idx = i + 1
                label = f"MO {mo_idx}"
                if energies is not None and i < len(energies):
                    label += f"   E={energies[i]:.5f} Ha"
                if i < n_occ:
                    if i == n_occ - 1:
                        label += "  <-- HOMO"
                elif i == n_occ:
                    label += "  <-- LUMO"
                
                item = QListWidgetItem(label)
                
                # Check existence
                path = self.get_cube_path(mo_idx)
                if os.path.exists(path):
                    item.setForeground(QBrush(QColor("darkgreen")))
                
                self.list_widget.addItem(item)
                
            # Scroll to HOMO
            # n_occ is number of occupied orbitals.
            # Indices are 0-based. So HOMO index is n_occ - 1.
            # Example: 21 occupied. HOMO is index 20 (MO 21).
            homo_idx = n_occ - 1
            if homo_idx >= 0 and homo_idx < self.list_widget.count():
                 self.list_widget.setCurrentRow(homo_idx)
                 self.list_widget.scrollToItem(self.list_widget.item(homo_idx), 
                                               QAbstractItemView.ScrollHint.PositionAtCenter)

    def on_single_click(self, item):
        row = self.list_widget.row(item)
        mo_idx = row + 1
        path = self.get_cube_path(mo_idx)
        if os.path.exists(path):
            self.visualize(path)

    def on_double_click(self, item):
        self.generate_current()

    def generate_current(self):
        row = self.list_widget.currentRow()
        if row < 0: return
        
        # Setup File Path
        base_dir = os.path.dirname(self.fchk_path)
        base_name = os.path.splitext(os.path.basename(self.fchk_path))[0]
        
        # Create output directory
        out_dir = os.path.join(base_dir, f"{base_name}_cubes")
        os.makedirs(out_dir, exist_ok=True)
        
        out_name = f"{base_name}_MO{row+1}.cube"
        out_path = os.path.join(out_dir, out_name)
        
        self.last_cube_path = out_path
        
        # Check Existence
        if os.path.exists(out_path):
            print(f"Using existing file: {out_path}")
            self.visualize(out_path)
            return

        # Lock UI
        self.btn_gen.setEnabled(False)
        self.list_widget.setEnabled(False)
        self.pbar.setValue(0)
        self.pbar.setVisible(True)
        
        mode_str = "MO"
        
        # Start Worker
        atoms = self.reader.get("Current cartesian coordinates").reshape(-1, 3)
        coeffs = self.reader.get("Alpha MO coefficients")
        
        self.worker = CalcWorker(
            self.engine, row, 
            self.spin_points.value(), 
            self.spin_gap.value(),
            atoms, coeffs, out_path,
            mode=mode_str
        )
        self.worker.progress_sig.connect(self.pbar.setValue)
        self.worker.finished_sig.connect(self.on_finished)
        self.worker.start()

    def on_finished(self, success, result):
        self.pbar.setVisible(False)
        self.btn_gen.setEnabled(True)
        self.list_widget.setEnabled(True)
        
        if success:
            # Mark as generated
            if self.worker:
                row = self.worker.mo_idx # This is 0-based index from list logic: row passed to worker
                if row >= 0 and row < self.list_widget.count():
                    item = self.list_widget.item(row)
                    item.setForeground(QBrush(QColor("darkgreen")))

            # Visualize
            self.visualize(result)
        else:
            QMessageBox.critical(self, "Error", f"Calculation failed: {result}")

    def visualize(self, cube_path):
        mw = self.context.get_main_window()
        
        # Get Isovalue
        val = self.spin_iso.value()
        op_val = self.spin_opacity.value()

        # Use our visualizer
        vis = CubeVisualizer(mw)
        if vis.load_file(cube_path):
            vis.show_iso(isovalue=val, opacity=op_val)
            # Only reset camera if it's the first time or user requests it?
            # Creating a new actor usually keeps camera if not reset.
            # But let's reset to ensure visibility.
            # mw.plotter.reset_camera() 
            mw.plotter.render()
            mw.statusBar().showMessage(f"Visualizing {os.path.basename(cube_path)} (Iso={val}, Opacity={op_val})")
        else:
            QMessageBox.warning(self, "Error", "Failed to load generated Cube file.")
