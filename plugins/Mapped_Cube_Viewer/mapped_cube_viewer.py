
import os
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QFileDialog, QDockWidget, QWidget, QVBoxLayout, 
                             QSlider, QLabel, QHBoxLayout, QPushButton, QMessageBox, QDoubleSpinBox, QProgressBar, QComboBox, QCheckBox)
from PyQt6.QtCore import Qt, QThread, pyqtSignal

# --- Dependency Management ---
try:
    from rdkit import Chem
    from rdkit import Geometry
except ImportError:
    Chem = None
    Geometry = None

# Periodic Table for Radii
try:
    from moleditpy.modules.constants import pt
except ImportError:
    try:
        from modules.constants import pt
    except ImportError:
        try:
            from rdkit import Chem
            pt = Chem.GetPeriodicTable()
        except:
            pt = None
            
__version__="2025.12.17"
__author__="HiroYokoyama"
PLUGIN_NAME = "Mapped Cube Viewer"

# --- Core Logic: Read Cube ---
def read_cube(filename):
    """
    Robust Cube file parser.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    if len(lines) < 6:
        raise ValueError("File content too short.")

    # Parse Header
    # L3: N_atoms, Origin
    try:
        tokens = lines[2].split()
        n_atoms_raw = int(tokens[0])
        n_atoms = abs(n_atoms_raw)
        origin = np.array([float(x) for x in tokens[1:4]])
        
        # L4: NX, Xvec
        tokens = lines[3].split()
        nx = int(tokens[0])
        x_vec = np.array([float(x) for x in tokens[1:4]])
        
        # L5: NY, Yvec
        tokens = lines[4].split()
        ny = int(tokens[0])
        y_vec = np.array([float(x) for x in tokens[1:4]])
        
        # L6: NZ, Zvec
        tokens = lines[5].split()
        nz = int(tokens[0])
        z_vec = np.array([float(x) for x in tokens[1:4]])
    except Exception as e:
        raise ValueError(f"Header parsing failed: {e}")

    # Unit Conversion (Bohr to Angstrom)
    # If N < 0, units are Bohr? Actually checking the sign of NX/NY/NZ is standard.
    # Gaussian Cube spec: If N < 0, the input file is Angstrom? No.
    # Standard: Distance units are Bohr.
    # However, some sources say: If NX > 0, units are Bohr. If NX < 0, units are Angstrom.
    # Let's check typical spec. 
    # Usually Cube files are Bohr.
    
    BOHR_TO_ANG = 0.529177
    
    # Check simple heuristic: if step size is very small (< 0.1), it might be Angstrom? 
    # Or if Origin is large?
    # Standard implementation assumes Bohr unless flagged.
    # Let's assume Bohr for now, consistent with previous code.
    
    # Fix negative counts logic
    if nx < 0: nx = abs(nx)
    if ny < 0: ny = abs(ny)
    if nz < 0: nz = abs(nz)
    
    # Scale spatial vectors
    origin *= BOHR_TO_ANG
    x_vec *= BOHR_TO_ANG
    y_vec *= BOHR_TO_ANG
    z_vec *= BOHR_TO_ANG
    
    # Parse Atoms
    atoms = []
    line_idx = 6
    for _ in range(n_atoms):
        tokens = lines[line_idx].split()
        line_idx += 1
        atomic_num = int(tokens[0])
        # charge = float(tokens[1]) # Skip
        pos = np.array([float(tokens[2]), float(tokens[3]), float(tokens[4])])
        pos *= BOHR_TO_ANG
        atoms.append((atomic_num, pos))
        
    # Skip extra orbital info if present (n_atoms_raw < 0)
    if n_atoms_raw < 0:
        tokens = lines[line_idx].split()
        # Usually checking number of orbitals
        try:
             n_mo = int(tokens[0])
             # Skip lines = n_mo
             line_idx += 1 # The line with n_mo
             line_idx += n_mo # Skip header lines for MOs?
             # Actually often it's just one line of indices?
             # Let's just consume until we see data pattern?
             # Safe fallback: look for line with many floats
        except: 
            pass
        # Logic for "skip header" is tricky without robust parser.
        # But commonly just +1 line.
        # Let's look for the start of volumetric data.
        
    # Volumetric Data
    # Flatten remaining lines
    data = []
    for l in lines[line_idx:]:
        parts = l.split()
        for p in parts:
            try:
                data.append(float(p))
            except: pass
            
    data_np = np.array(data)
    expected = nx * ny * nz
    
    # If we parsed too few, that's bad. Too many? Truncate.
    if len(data_np) < expected:
        # Maybe we missed the header skipping or it was really short
        # Try to recover: take last N values
        pass 
    if len(data_np) > expected:
        data_np = data_np[:expected] # Truncate header noise if any

    vol_data = data_np.reshape((nx, ny, nz)) # Cube is usually X outer, then Y, then Z inner?
    # Actually Cube is: loops X, then Y, then Z.
    # So flatten order is C-like?
    # Z indices change fastest.
    # vol_data[ix, iy, iz]
    
    # Create PyVista Grid
    # Meshgrid
    # Note: numpy meshgrid with indexing='ij' gives [nx, ny, nz]
    
    grid = pv.ImageData() 
    # StructuredGrid is better for non-orthogonal, but ImageData is faster if orthogonal.
    # Are Cube axes always orthogonal? usually YES. 
    # But vectors are provided. If they have off-diagonal, we need StructuredGrid.
    
    is_orthogonal = True
    # check orthogonality
    if (x_vec[1]!=0 or x_vec[2]!=0 or 
        y_vec[0]!=0 or y_vec[2]!=0 or 
        z_vec[0]!=0 or z_vec[1]!=0):
        is_orthogonal = False

    if is_orthogonal:
        grid = pv.ImageData()
        grid.dimensions = (nx, ny, nz)
        grid.origin = origin
        # spacing requires scalar if orthogonal
        grid.spacing = (x_vec[0], y_vec[1], z_vec[2])
        grid.point_data["property_values"] = vol_data.flatten(order='F') # Z fastest
    else:
        # Structured Grid for rotated cubes
        x_range = np.arange(nx)
        y_range = np.arange(ny)
        z_range = np.arange(nz)
        gx, gy, gz = np.meshgrid(x_range, y_range, z_range, indexing='ij')
        
        # Flatten for calculation
        gx, gy, gz = gx.flatten(), gy.flatten(), gz.flatten()
        points = (origin + 
                  np.outer(gx, x_vec) + 
                  np.outer(gy, y_vec) + 
                  np.outer(gz, z_vec))
        
        grid = pv.StructuredGrid()
        grid.dimensions = [nx, ny, nz]
        grid.points = points
        grid.point_data["property_values"] = vol_data.flatten(order='F')

    return {"atoms": atoms}, grid

from PyQt6.QtWidgets import (QDialog, QLineEdit, QFormLayout, QDialogButtonBox)

class MappedCubeSetupDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Mapped Cube Setup")
        self.surface_file = None
        self.property_file = None
        
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        form = QFormLayout()
        
        # File 1: Surface
        self.le_surf = QLineEdit()
        btn_surf = QPushButton("Browse")
        btn_surf.clicked.connect(self.browse_surf)
        h1 = QHBoxLayout()
        h1.addWidget(self.le_surf)
        h1.addWidget(btn_surf)
        form.addRow("Surface (Density):", h1)
        
        # File 2: Property
        self.le_prop = QLineEdit()
        btn_prop = QPushButton("Browse")
        btn_prop.clicked.connect(self.browse_prop)
        h2 = QHBoxLayout()
        h2.addWidget(self.le_prop)
        h2.addWidget(btn_prop)
        form.addRow("Property (Color):", h2)
        
        layout.addLayout(form)
        
        # Buttons
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        
    def browse_surf(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Surface Cube", "", "Cube (*.cube *.cub)")
        if f: self.le_surf.setText(f)
            
    def browse_prop(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Property Cube", "", "Cube (*.cube *.cub)")
        if f: self.le_prop.setText(f)
        
    def accept(self):
        self.surface_file = self.le_surf.text()
        self.property_file = self.le_prop.text()
        if not os.path.exists(self.surface_file):
            QMessageBox.warning(self, "Error", "Surface file not found.")
            return
        if not os.path.exists(self.property_file):
            QMessageBox.warning(self, "Error", "Property file not found.")
            return
        super().accept()

# --- Dock Widget ---
class MappedWidget(QWidget):
    def __init__(self, mw, dock, grid_surf, grid_prop):
        super().__init__()
        self.mw = mw
        self.dock = dock
        self.grid_surf = grid_surf
        self.grid_prop = grid_prop
        self.actor = None
        
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        
        # Controls
        # Isovalue
        self.layout.addWidget(QLabel("<b>Surface Isovalue:</b>"))
        self.iso_spin = QDoubleSpinBox()
        self.iso_spin.setRange(-1000, 1000)
        self.iso_spin.setDecimals(5)
        self.iso_spin.setSingleStep(0.001)
        
        # Default ISO
        vals = self.grid_surf.point_data.get("property_values", np.array([0]))
        default_iso = 0.00200 if vals.max() > 0.1 else vals.mean()
        self.iso_spin.setValue(default_iso)
        self.iso_spin.setKeyboardTracking(False)
        self.iso_spin.valueChanged.connect(lambda: self.update_mesh(auto_fit=False))
        self.layout.addWidget(self.iso_spin)
        
        # Color Range
        self.layout.addWidget(QLabel("<b>Property Color Range:</b>"))
        
        # Colormap Selector
        cbox = QHBoxLayout()
        cbox.addWidget(QLabel("Style:"))
        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems([
            'jet_r', 'jet', 
            'rainbow_r', 'rainbow', 
            'bwr_r', 'bwr', 
            'seismic_r', 'seismic', 
            'coolwarm_r', 'coolwarm', 
            'viridis_r', 'viridis', 
            'plasma_r', 'plasma', 
            'magma_r', 'magma', 
            'inferno_r', 'inferno'
        ])
        self.cmap_combo.setCurrentText('jet_r')
        self.cmap_combo.currentTextChanged.connect(lambda: self.update_mesh(auto_fit=False))
        cbox.addWidget(self.cmap_combo)
        self.layout.addLayout(cbox)
        
        rbox = QHBoxLayout()
        
        pvals = self.grid_prop.point_data.get("property_values", np.array([-0.1, 0.1]))
        vmin, vmax = pvals.min(), pvals.max()
        
        self.min_spin = QDoubleSpinBox()
        self.min_spin.setRange(-1e20, 1e20)
        self.min_spin.setDecimals(6)
        self.min_spin.setSingleStep(0.01)
        self.min_spin.setKeyboardTracking(False)
        self.min_spin.setValue(vmin)
        self.min_spin.valueChanged.connect(lambda: self.update_mesh(auto_fit=False))
        
        self.max_spin = QDoubleSpinBox()
        self.max_spin.setRange(-1e20, 1e20)
        self.max_spin.setDecimals(6)
        self.max_spin.setSingleStep(0.01)
        self.max_spin.setKeyboardTracking(False)
        self.max_spin.setValue(vmax)
        self.max_spin.valueChanged.connect(lambda: self.update_mesh(auto_fit=False))
        
        rbox.addWidget(QLabel("Min:"))
        rbox.addWidget(self.min_spin)
        rbox.addWidget(QLabel("Max:"))
        rbox.addWidget(self.max_spin)
        self.layout.addLayout(rbox)
        
        btn_fit = QPushButton("Fit Range to Surface")
        btn_fit.clicked.connect(lambda: self.update_mesh(auto_fit=True))
        self.layout.addWidget(btn_fit)
        
        # Opacity
        self.layout.addWidget(QLabel("<b>Opacity:</b>"))
        self.opacity_spin = QDoubleSpinBox()
        self.opacity_spin.setRange(0, 1)
        self.opacity_spin.setSingleStep(0.1)
        self.opacity_spin.setValue(0.4)
        self.opacity_spin.valueChanged.connect(lambda: self.update_mesh(auto_fit=False))
        self.layout.addWidget(self.opacity_spin)
        
        # Exports
        ebox = QHBoxLayout()
        btn_view = QPushButton("Export View")
        btn_view.clicked.connect(self.export_view)
        ebox.addWidget(btn_view)
        
        btn_cbar = QPushButton("Export Color Bar")
        btn_cbar.clicked.connect(self.export_colorbar)
        ebox.addWidget(btn_cbar)
        self.layout.addLayout(ebox)
        
        # Options
        self.check_transparent = QCheckBox("Transparent Background")
        self.check_transparent.setChecked(True)
        self.layout.addWidget(self.check_transparent)
        
        # Close
        self.layout.addStretch()
        btn = QPushButton("Close")
        btn.clicked.connect(self.close_plugin)
        self.layout.addWidget(btn)
        
        # Initial Draw
        self.update_mesh(auto_fit=True)

    def export_view(self):
        f, _ = QFileDialog.getSaveFileName(self, "Save View", "mapped_view.png", "PNG Image (*.png)")
        if f:
            try:
                self.mw.plotter.screenshot(f, transparent_background=self.check_transparent.isChecked())
                QMessageBox.information(self, "Success", f"Saved View to {f}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Save failed: {e}")

    def export_colorbar(self):
        f, _ = QFileDialog.getSaveFileName(self, "Save Color Bar", "mapped_colorbar.png", "PNG Image (*.png)")
        if f:
            try:
                # Create off-screen plotter for just the bar
                # Use a wide layout (Horizontal) - Even larger and wider
                pl = pv.Plotter(off_screen=True, window_size=(2000, 400))
                
                # Checkbox sensitive background
                if self.check_transparent.isChecked():
                    pl.set_background(None) # Transparent
                else:
                    pl.set_background('white')
                
                # Dummy mesh to attach mapper
                vmin = self.min_spin.value()
                vmax = self.max_spin.value()
                cmap = self.cmap_combo.currentText()
                
                # Box with scalars covering range
                mesh = pv.Box()
                mesh.point_data['scalars'] = np.linspace(vmin, vmax, mesh.n_points)
                
                # Add mesh with scalar bar args directly
                pl.add_mesh(mesh, 
                            scalars='scalars', 
                            cmap=cmap, 
                            clim=[vmin, vmax],
                            show_scalar_bar=True, # Show it via add_mesh
                            opacity=0.0, # Invisible mesh
                            scalar_bar_args={
                                "title": "",
                                "vertical": False,
                                "n_labels": 5,
                                "italic": False,
                                "bold": False,
                                "title_font_size": 1,
                                "label_font_size": 50, # Slightly reduced but still large
                                "color": "black",
                                "height": 0.4,      # THICK bar (Restored)
                                "width": 0.8,       # 80% width
                                "position_x": 0.1,  # Centered (10% margin)
                                "position_y": 0.3
                            }
                )
                
                # FIX ALIGNMENT & TICKS
                if hasattr(pl, 'scalar_bar'):
                    try:
                        sb = pl.scalar_bar
                        sb.GetLabelTextProperty().SetJustificationToCentered()
                    except: pass
                
                pl.screenshot(f, transparent_background=self.check_transparent.isChecked())
                pl.close()
                QMessageBox.information(self, "Success", f"Saved Color Bar to {f}")
            except Exception as e:
                import traceback
                traceback.print_exc()
                QMessageBox.critical(self, "Error", f"Save failed: {e}")

    def update_mesh(self, auto_fit=False):
        try:
            iso_val = self.iso_spin.value()
            opacity = self.opacity_spin.value()
            
            # 1. Generate Surface
            iso = self.grid_surf.contour([iso_val], scalars="property_values")
            if iso.n_points == 0:
                pass
                return
            
            # 2. Sample Property
            mapped = iso.sample(self.grid_prop)
            
            # 3. Auto Fit
            clim = [self.min_spin.value(), self.max_spin.value()]
            if auto_fit:
                mvals = mapped.point_data.get("property_values")
                if mvals is not None and len(mvals) > 0:
                    vmin, vmax = mvals.min(), mvals.max()
                    if vmax - vmin < 1e-9: vmax += 0.001
                    clim = [vmin, vmax]
                    
                    self.min_spin.blockSignals(True); self.min_spin.setValue(vmin); self.min_spin.blockSignals(False)
                    self.max_spin.blockSignals(True); self.max_spin.setValue(vmax); self.max_spin.blockSignals(False)
            
            # 4. Render
            if self.actor:
                self.mw.plotter.remove_actor(self.actor)
                
            self.actor = self.mw.plotter.add_mesh(
                mapped,
                scalars="property_values",
                cmap=self.cmap_combo.currentText(),
                clim=clim,
                smooth_shading=True,
                opacity=opacity,
                name="mapped_mesh_dock",
                scalar_bar_args={'title': ''}
            )
            

            
            self.mw.plotter.render()
            
        except Exception as e:
            print(f"Update Error: {e}")

    def close_plugin(self):
        try:
            self.mw.plotter.remove_actor(self.actor)
            self.mw.plotter.render()
        except: pass
        if self.dock:
            self.dock.close()
        try:
            self.mw.clear_all()
        except: pass
        self.close()

# --- Entry Point ---
def run(mw):
    # 1. Setup Dialog
    dlg = MappedCubeSetupDialog(mw)
    if dlg.exec() != QDialog.DialogCode.Accepted:
        return
        
    s_file = dlg.surface_file
    p_file = dlg.property_file
    
    if not s_file or not p_file: return
    
    try:
        from PyQt6.QtWidgets import QApplication, QProgressDialog
        
        # 2. Read Files
        progress = QProgressDialog("Reading Cube Files...", "Cancel", 0, 2, mw)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.show()
        

        meta1, grid_surf = read_cube(s_file)
        progress.setValue(1)
        QApplication.processEvents()
        

        meta2, grid_prop = read_cube(p_file)
        progress.setValue(2)
        QApplication.processEvents()
        progress.close()
        
        # 3. Visualize Atoms
        atoms = meta1['atoms']
        if Chem:
            s = f"{len(atoms)}\n\n"
            for n, p in atoms:
                sym = pt.GetElementSymbol(n) if pt else "C"
                s += f"{sym} {p[0]} {p[1]} {p[2]}\n"
            mol = Chem.MolFromXYZBlock(s)
            try: 
                 from rdkit.Chem import rdDetermineBonds
                 rdDetermineBonds.DetermineConnectivity(mol)
            except: pass
            
            mw.current_mol = mol
            if hasattr(mw, 'draw_molecule_3d'):
                mw.draw_molecule_3d(mol)

        # 4. Enter 3D Mode (Hides other panels)
        if hasattr(mw, 'main_window_ui_manager'):
            try: mw.main_window_ui_manager._enter_3d_viewer_ui_mode()
            except: pass

        # 5. Launch Dock
        # Close old docks
        for d in mw.findChildren(QDockWidget):
            if d.windowTitle() == "Mapped Viewer":
                d.close()

        dock = QDockWidget("Mapped Viewer", mw)
        dock.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea)
        widget = MappedWidget(mw, dock, grid_surf, grid_prop)
        dock.setWidget(widget)
        mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        QMessageBox.critical(mw, "Error", f"Failed:\n{e}")

