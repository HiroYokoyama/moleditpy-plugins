import os
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QFileDialog, QDockWidget, QWidget, QVBoxLayout, 
                             QSlider, QLabel, QHBoxLayout, QPushButton, QMessageBox, QDoubleSpinBox, QProgressBar, QComboBox, QCheckBox, QDialog, QLineEdit, QFormLayout, QDialogButtonBox)
from PyQt6.QtCore import Qt

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
            
__version__="2025.12.21" # Fixed version
__author__="HiroYokoyama"
PLUGIN_NAME = "Mapped Cube Viewer"

# --- Core Logic: Robust Parser from cube_viewer.py ---

def parse_cube_data(filename):
    """
    Robust Cube file parser copied from cube_viewer.py
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    if len(lines) < 6:
        raise ValueError("File too short to be a Cube file.")

    # --- Header Parsing ---
    tokens = lines[2].split()
    n_atoms_raw = int(tokens[0])
    n_atoms = abs(n_atoms_raw)
    origin_raw = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])])

    def parse_vec(line):
        t = line.split()
        return int(t[0]), np.array([float(t[1]), float(t[2]), float(t[3])])

    nx, x_vec_raw = parse_vec(lines[3])
    ny, y_vec_raw = parse_vec(lines[4])
    nz, z_vec_raw = parse_vec(lines[5])
    
    is_angstrom_header = (nx < 0 or ny < 0 or nz < 0)
    nx, ny, nz = abs(nx), abs(ny), abs(nz)

    # --- Atoms Parsing ---
    atoms = []
    current_line = 6
    if n_atoms_raw < 0:
        try:
            parts = lines[current_line].split()
            if len(parts) != 5: 
                 current_line += 1
        except:
             current_line += 1

    for _ in range(n_atoms):
        line = lines[current_line].split()
        current_line += 1
        atomic_num = int(line[0])
        try:
            x, y, z = float(line[2]), float(line[3]), float(line[4])
        except:
            x, y, z = 0.0, 0.0, 0.0
        atoms.append((atomic_num, np.array([x, y, z])))

    # --- Volumetric Data Parsing (Skip Metadata) ---
    while current_line < len(lines):
        line_content = lines[current_line].strip()
        parts = line_content.split()
        if not parts:
            current_line += 1
            continue
        if len(parts) < 6: # Heuristic: Data lines usually have 6 cols
            current_line += 1
            continue
        try:
            float(parts[0])
        except ValueError:
            current_line += 1
            continue
        break

    full_str = " ".join(lines[current_line:])
    try:
        data_values = np.fromstring(full_str, sep=' ')
    except:
        data_values = np.array([])
    
    expected_size = nx * ny * nz
    actual_size = len(data_values)
    
    if actual_size > expected_size:
        excess = actual_size - expected_size
        data_values = data_values[excess:] # Trim from start if header noise included
    elif actual_size < expected_size:
        pad = np.zeros(expected_size - actual_size)
        data_values = np.concatenate((data_values, pad))
    
    return {
        "atoms": atoms,
        "origin": origin_raw,
        "x_vec": x_vec_raw,
        "y_vec": y_vec_raw,
        "z_vec": z_vec_raw,
        "dims": (nx, ny, nz),
        "data_flat": data_values,
        "is_angstrom_header": is_angstrom_header
    }

def build_grid_from_meta(meta):
    """
    Reconstructs the PyVista grid.
    Correctly handles standard Cube (Z-fast) mapping.
    """
    nx, ny, nz = meta['dims']
    origin = meta['origin'].copy()
    x_vec = meta['x_vec'].copy()
    y_vec = meta['y_vec'].copy()
    z_vec = meta['z_vec'].copy()
    atoms = []
    
    BOHR_TO_ANGSTROM = 0.529177210903
    convert_to_ang = True
    if meta['is_angstrom_header']:
        convert_to_ang = False
            
    if convert_to_ang:
        origin *= BOHR_TO_ANGSTROM
        x_vec *= BOHR_TO_ANGSTROM
        y_vec *= BOHR_TO_ANGSTROM
        z_vec *= BOHR_TO_ANGSTROM
        
    for anum, pos in meta['atoms']:
        p = pos.copy()
        if convert_to_ang:
            p *= BOHR_TO_ANGSTROM
        atoms.append((anum, p))

    # Grid Points Generation
    x_range = np.arange(nx)
    y_range = np.arange(ny)
    z_range = np.arange(nz)
    
    gx, gy, gz = np.meshgrid(x_range, y_range, z_range, indexing='ij')
    
    # Flatten using Fortran order (X-fast) for VTK Structure
    # !!! IMPORTANT FIX: Must match Flatten order of data !!!
    gx_f = gx.flatten(order='F')
    gy_f = gy.flatten(order='F')
    gz_f = gz.flatten(order='F')
    
    points = (origin + 
              np.outer(gx_f, x_vec) + 
              np.outer(gy_f, y_vec) + 
              np.outer(gz_f, z_vec))
    
    grid = pv.StructuredGrid()
    grid.points = points
    grid.dimensions = [nx, ny, nz]
    
    # Data Mapping
    raw_data = meta['data_flat']
    
    # Standard Cube: Z-Fast (C-order reshape)
    # Then flatten F-order to match the X-fast points we just generated
    vol_3d = raw_data.reshape((nx, ny, nz), order='C')
    
    # Key name "property_values" for Mapped Viewer
    grid.point_data["property_values"] = vol_3d.flatten(order='F')
        
    return {"atoms": atoms}, grid

def read_cube(filename):
    meta = parse_cube_data(filename)
    return build_grid_from_meta(meta)

# --- Dialog & Widget (Identical Logic, just ensured imports) ---

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
        self.layout.addWidget(QLabel("<b>Surface Isovalue:</b>"))
        self.iso_spin = QDoubleSpinBox()
        self.iso_spin.setRange(-1000, 1000)
        self.iso_spin.setDecimals(5)
        self.iso_spin.setSingleStep(0.001)
        
        # Default ISO
        vals = self.grid_surf.point_data.get("property_values", np.array([0]))
        # Robust default calculation
        if len(vals) > 0:
            if vals.max() > 0.1:
                default_iso = 0.00200
            else:
                default_iso = float(np.mean(vals))
        else:
            default_iso = 0.002
            
        self.iso_spin.setValue(default_iso)
        self.iso_spin.setKeyboardTracking(False)
        self.iso_spin.valueChanged.connect(lambda: self.update_mesh(auto_fit=False))
        self.layout.addWidget(self.iso_spin)
        
        # Color Range
        self.layout.addWidget(QLabel("<b>Property Color Range:</b>"))
        
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
        if len(pvals) > 0:
            vmin, vmax = pvals.min(), pvals.max()
        else:
            vmin, vmax = -0.1, 0.1
            
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
        
        self.layout.addStretch()
        btn = QPushButton("Close")
        btn.clicked.connect(self.close_plugin)
        self.layout.addWidget(btn)
        
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
                pl = pv.Plotter(off_screen=True, window_size=(2000, 400))
                if self.check_transparent.isChecked():
                    pl.set_background(None)
                else:
                    pl.set_background('white')
                
                vmin = self.min_spin.value()
                vmax = self.max_spin.value()
                cmap = self.cmap_combo.currentText()
                
                mesh = pv.Box()
                mesh.point_data['scalars'] = np.linspace(vmin, vmax, mesh.n_points)
                
                pl.add_mesh(mesh, 
                            scalars='scalars', 
                            cmap=cmap, 
                            clim=[vmin, vmax],
                            show_scalar_bar=True,
                            opacity=0.0,
                            scalar_bar_args={
                                "title": "",
                                "vertical": False,
                                "n_labels": 5,
                                "italic": False,
                                "bold": False,
                                "title_font_size": 1,
                                "label_font_size": 50,
                                "color": "black",
                                "height": 0.4,
                                "width": 0.8,
                                "position_x": 0.1,
                                "position_y": 0.3
                            }
                )
                
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
            
            iso = self.grid_surf.contour([iso_val], scalars="property_values")
            if iso.n_points == 0:
                return
            
            mapped = iso.sample(self.grid_prop)
            
            clim = [self.min_spin.value(), self.max_spin.value()]
            if auto_fit:
                mvals = mapped.point_data.get("property_values")
                if mvals is not None and len(mvals) > 0:
                    vmin, vmax = mvals.min(), mvals.max()
                    if vmax - vmin < 1e-9: vmax += 0.001
                    clim = [vmin, vmax]
                    
                    self.min_spin.blockSignals(True); self.min_spin.setValue(vmin); self.min_spin.blockSignals(False)
                    self.max_spin.blockSignals(True); self.max_spin.setValue(vmax); self.max_spin.blockSignals(False)
            
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
    dlg = MappedCubeSetupDialog(mw)
    if dlg.exec() != QDialog.DialogCode.Accepted:
        return
        
    s_file = dlg.surface_file
    p_file = dlg.property_file
    
    if not s_file or not p_file: return
    
    try:
        from PyQt6.QtWidgets import QApplication, QProgressDialog
        
        progress = QProgressDialog("Reading Cube Files...", "Cancel", 0, 2, mw)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.show()
        
        # Using the robust read_cube logic
        meta1, grid_surf = read_cube(s_file)
        progress.setValue(1)
        QApplication.processEvents()
        
        meta2, grid_prop = read_cube(p_file)
        progress.setValue(2)
        QApplication.processEvents()
        progress.close()
        
        # Visualize Atoms
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

        if hasattr(mw, 'main_window_ui_manager'):
            try: mw.main_window_ui_manager._enter_3d_viewer_ui_mode()
            except: pass

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
