import os
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, 
                             QListWidget, QLabel, QDoubleSpinBox, QProgressBar, 
                             QGroupBox, QMessageBox, QWidget, QAbstractItemView)
from PyQt6.QtCore import Qt, QThread, pyqtSignal

# Relative imports
from .analyzer import FCHKReader, BasisSetEngine, CubeWriter
from .vis import CubeVisualizer

class CalcWorker(QThread):
    progress_sig = pyqtSignal(int)
    finished_sig = pyqtSignal(bool, str) # success, message/path

    def __init__(self, engine, mo_idx, resolution, margin, atoms, mo_coeffs, output_path, mode="MO"):
        super().__init__()
        self.engine = engine
        self.mo_idx = mo_idx
        self.resolution = resolution
        self.margin = margin
        self.atoms = atoms
        self.mo_coeffs = mo_coeffs
        self.output_path = output_path
        self.mode = mode
        self._is_cancelled = False

    def run(self):
        try:
            # Grid Definition
            min_c = np.min(self.atoms, axis=0) - self.margin
            max_c = np.max(self.atoms, axis=0) + self.margin
            span = max_c - min_c
            
            nx = int(np.ceil(span[0] / self.resolution))
            ny = int(np.ceil(span[1] / self.resolution))
            nz = int(np.ceil(span[2] / self.resolution))
            
            # Generate points
            x = np.linspace(min_c[0], min_c[0] + (nx-1)*self.resolution, nx)
            y = np.linspace(min_c[1], min_c[1] + (ny-1)*self.resolution, ny)
            z = np.linspace(min_c[2], min_c[2] + (nz-1)*self.resolution, nz)
            
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
            grid_points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
            
            n_total = grid_points.shape[0]
            chunk_size = 50000
            result_flat = np.zeros(n_total)
            
            for i, start in enumerate(range(0, n_total, chunk_size)):
                if self._is_cancelled: return
                
                end = min(start + chunk_size, n_total)
                chunk_pts = grid_points[start:end]
                
                # Evaluate
                if self.mode == "MO":
                    val = self.engine.evaluate_mo_on_grid(self.mo_idx, chunk_pts, self.mo_coeffs)
                else:
                    # self.mo_idx is basis_idx in this context
                    val = self.engine.evaluate_basis_on_grid(self.mo_idx, chunk_pts)
                    
                result_flat[start:end] = val
                
                # Progress
                pct = int((end / n_total) * 100)
                self.progress_sig.emit(pct)
            
            # Write
            vol_data = result_flat.reshape(nx, ny, nz)
            vectors = np.identity(3) * self.resolution
            
            atom_nos = self.engine.fchk.get("Atomic numbers")
            
            label = f"MO {self.mo_idx+1}" if self.mode == "MO" else f"Basis {self.mo_idx+1}"
            CubeWriter.write(self.output_path, self.atoms, atom_nos, min_c, vectors, vol_data, comment=label)
            
            self.finished_sig.emit(True, self.output_path)
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.finished_sig.emit(False, str(e))

class OrbitalDialog(QDialog):
    def __init__(self, parent, context, fchk_path):
        super().__init__(parent)
        # Make this a Tool window so it stays on top of parent but doesn't block (Modeless)
        # and doesn't sit on top of other apps (unlike WindowStaysOnTopHint)
        self.setWindowFlags(Qt.WindowType.Tool)
        
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
                    
                    # Push undo state so user can revert
                    mw = context.get_main_window()
                    mw.push_undo_state()
                    mw.plotter.reset_camera()
        except Exception as e:
            print(f"Failed to load molecule structure: {e}")
        
        # Layout
        layout = QVBoxLayout(self)
        
        # Settings
        grp = QGroupBox("Generation Settings")
        form = QHBoxLayout(grp)
        
        form.addWidget(QLabel("Resolution (Bohr):"))
        self.spin_res = QDoubleSpinBox()
        self.spin_res.setRange(0.05, 1.0)
        self.spin_res.setValue(0.2)
        self.spin_res.setSingleStep(0.05)
        form.addWidget(self.spin_res)
        
        form.addWidget(QLabel("Margin (Bohr):"))
        self.spin_gap = QDoubleSpinBox()
        self.spin_gap.setRange(1.0, 20.0) 
        self.spin_gap.setValue(4.0) 
        form.addWidget(self.spin_gap)
        
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
        
        # Mode Selection (Disabled/Commented out per user request)
        # ...
        
        # List
        self.list_widget = QListWidget()
        self.list_widget.itemDoubleClicked.connect(self.on_double_click)
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
                self.list_widget.addItem(label)
                
            # Scroll to HOMO
            # n_occ is number of occupied orbitals.
            # Indices are 0-based. So HOMO index is n_occ - 1.
            # Example: 21 occupied. HOMO is index 20 (MO 21).
            homo_idx = n_occ - 1
            if homo_idx >= 0 and homo_idx < self.list_widget.count():
                 self.list_widget.setCurrentRow(homo_idx)
                 self.list_widget.scrollToItem(self.list_widget.item(homo_idx), 
                                               QAbstractItemView.ScrollHint.PositionAtCenter)

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
            self.spin_res.value(), 
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
            # Visualize
            self.visualize(result)
        else:
            QMessageBox.critical(self, "Error", f"Calculation failed: {result}")

    def visualize(self, cube_path):
        mw = self.context.get_main_window()
        
        # Get Isovalue
        val = self.spin_iso.value()

        # Use our visualizer
        vis = CubeVisualizer(mw)
        if vis.load_file(cube_path):
            vis.show_iso(isovalue=val)
            # Only reset camera if it's the first time or user requests it?
            # Creating a new actor usually keeps camera if not reset.
            # But let's reset to ensure visibility.
            # mw.plotter.reset_camera() 
            mw.plotter.render()
            mw.statusBar().showMessage(f"Visualizing {os.path.basename(cube_path)} (Iso={val})")
        else:
            QMessageBox.warning(self, "Error", "Failed to load generated Cube file.")
