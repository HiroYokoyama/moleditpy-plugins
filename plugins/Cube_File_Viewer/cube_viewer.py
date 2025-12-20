
import os
import tempfile
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QFileDialog, QDockWidget, QWidget, QVBoxLayout, 
                             QSlider, QLabel, QHBoxLayout, QPushButton, QMessageBox, 
                             QDoubleSpinBox, QColorDialog, QInputDialog, QDialog, 
                             QFormLayout, QDialogButtonBox, QSpinBox, QCheckBox)
from PyQt6.QtGui import QColor
from PyQt6.QtCore import Qt

# RDKit imports for molecule construction
try:
    from rdkit import Chem
    from rdkit import Geometry
    try:
        from rdkit.Chem import rdDetermineBonds
    except ImportError:
        rdDetermineBonds = None
except ImportError:
    Chem = None
    Geometry = None
    rdDetermineBonds = None
    
__version__="2025.12.20"
__author__="HiroYokoyama"
PLUGIN_NAME = "Cube File Viewer"

def parse_cube_data(filename):
    """
    Parses a Gaussian Cube file and returns raw data structures.
    Adapted from test.py with robust header handling.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    if len(lines) < 6:
        raise ValueError("File too short to be a Cube file.")

    # --- Header Parsing ---
    # Line 3: Natoms, Origin
    tokens = lines[2].split()
    n_atoms_raw = int(tokens[0])
    n_atoms = abs(n_atoms_raw)
    origin_raw = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])])

    # Lines 4-6: NX, NY, NZ and vectors
    def parse_vec(line):
        t = line.split()
        return int(t[0]), np.array([float(t[1]), float(t[2]), float(t[3])])

    nx, x_vec_raw = parse_vec(lines[3])
    ny, y_vec_raw = parse_vec(lines[4])
    nz, z_vec_raw = parse_vec(lines[5])
    
    # Auto-detect units based on sign of NX/NY/NZ (Gaussian standard)
    is_angstrom_header = (nx < 0 or ny < 0 or nz < 0)
    
    nx, ny, nz = abs(nx), abs(ny), abs(nz)

    # --- Atoms Parsing ---
    atoms = []
    current_line = 6
    if n_atoms_raw < 0:
        # Smart check: if next line does not look like an atom line, skip it.
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

    # --- Volumetric Data Parsing ---
    # Read rest of file
    full_str = " ".join(lines[current_line:])
    try:
        data_values = np.fromstring(full_str, sep=' ')
    except:
        data_values = np.array([])
    
    expected_size = nx * ny * nz
    actual_size = len(data_values)
    
    # FIX: Trim from START if excess > 0 (The header values are at the start)
    # This logic fixes the shift issue test.py had.
    if actual_size > expected_size:
        excess = actual_size - expected_size
        data_values = data_values[excess:]
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
    
    # Units Handling (replicated logic)
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

    # Grid Points Generation (Matches test.py logic for F-order consistency)
    x_range = np.arange(nx)
    y_range = np.arange(ny)
    z_range = np.arange(nz)
    
    gx, gy, gz = np.meshgrid(x_range, y_range, z_range, indexing='ij')
    
    # Flatten using Fortran order (X-fast) for VTK Structure
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
    # This logic matches test.py's implementation which user preferred
    vol_3d = raw_data.reshape((nx, ny, nz), order='C')
    grid.point_data["values"] = vol_3d.flatten(order='F')
        
    return {"atoms": atoms}, grid

def read_cube(filename):
    meta = parse_cube_data(filename)
    return build_grid_from_meta(meta)

class CubeViewerWidget(QWidget):
    def __init__(self, parent_window, dock_widget, grid):
        super().__init__(parent_window)
        self.mw = parent_window
        self.dock = dock_widget 
        self.grid = grid
        self.iso_actor_p = None
        self.iso_actor_n = None
        
        self.plotter = self.mw.plotter
        
        # Initial Colors
        self.color_p = "#0000FF" # Blue
        self.color_n = "#FF0000" # Red

        self.init_ui()
        self.update_iso() 

    def init_ui(self):
        layout = QVBoxLayout()
        
        # --- Isovalue Controls ---
        ctrl_layout = QHBoxLayout()
        ctrl_layout.addWidget(QLabel("Isovalue:"))
        
        # Spinbox: Fixed max 1.0 as requested
        self.max_val = 1.0
        
        self.spin = QDoubleSpinBox()
        self.spin.setRange(0.0001, self.max_val) 
        self.spin.setSingleStep(0.002)
        self.spin.setDecimals(4)
        self.spin.setValue(0.05)
        
        self.spin.valueChanged.connect(self.on_spin_changed)
        ctrl_layout.addWidget(self.spin)
        
        # Slider
        self.slider_max_int = 1000
        
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(0, self.slider_max_int)
        
        # Default 0.05
        default_val = 0.05
        self.slider.setValue(int((default_val / self.max_val) * self.slider_max_int))
        
        self.slider.valueChanged.connect(self.on_slider_changed)
        ctrl_layout.addWidget(self.slider)

        layout.addLayout(ctrl_layout)
        
        # --- Color Controls ---
        # Fixed size for color buttons to align them nicely
        btn_size = (50, 24)
        
        # Positive
        pos_color_layout = QHBoxLayout()
        pos_color_layout.addWidget(QLabel("Pos Color:"))
        self.btn_color_p = QPushButton()
        self.btn_color_p.setFixedSize(*btn_size)
        self.btn_color_p.setStyleSheet(f"background-color: {self.color_p}; border: 1px solid gray;")
        self.btn_color_p.clicked.connect(self.choose_color_p)
        pos_color_layout.addWidget(self.btn_color_p)
        pos_color_layout.addStretch() # Align left
        layout.addLayout(pos_color_layout)

        # Negative
        neg_color_layout = QHBoxLayout()
        neg_color_layout.addWidget(QLabel("Neg Color:"))
        self.btn_color_n = QPushButton()
        self.btn_color_n.setFixedSize(*btn_size)
        self.btn_color_n.setStyleSheet(f"background-color: {self.color_n}; border: 1px solid gray;")
        self.btn_color_n.clicked.connect(self.choose_color_n)
        neg_color_layout.addWidget(self.btn_color_n)
        neg_color_layout.addStretch() # Align left
        layout.addLayout(neg_color_layout)

        # Complementary Checkbox (Next line)
        comp_layout = QHBoxLayout()
        self.check_comp_color = QCheckBox("Use complementary color for neg")
        self.check_comp_color.toggled.connect(self.on_comp_color_toggled)
        comp_layout.addWidget(self.check_comp_color)
        comp_layout.addStretch()
        layout.addLayout(comp_layout)

        # --- Opacity Controls ---
        opacity_layout = QHBoxLayout()
        # Static label (no numeric value displayed)
        self.opacity_label = QLabel("Opacity:")
        opacity_layout.addWidget(self.opacity_label)
        
        # Numeric opacity input: place before the slider to match isovalue controls
        self.opacity_spin = QDoubleSpinBox()
        self.opacity_spin.setRange(0.0, 1.0)
        self.opacity_spin.setSingleStep(0.01)
        self.opacity_spin.setDecimals(2)
        self.opacity_spin.setValue(0.4)
        self.opacity_spin.valueChanged.connect(self.on_opacity_spin_changed)
        opacity_layout.addWidget(self.opacity_spin)

        self.opacity_slider = QSlider(Qt.Orientation.Horizontal)
        self.opacity_slider.setRange(0, 100)
        self.opacity_slider.setValue(40) # Match spinbox
        self.opacity_slider.valueChanged.connect(self.on_opacity_changed)
        opacity_layout.addWidget(self.opacity_slider)
        
        layout.addLayout(opacity_layout)

        close_btn = QPushButton("Close Plugin")
        close_btn.clicked.connect(self.close_plugin)
        layout.addWidget(close_btn)
        
        layout.addStretch()
        self.setLayout(layout)

    def on_comp_color_toggled(self, checked):
        self.btn_color_n.setEnabled(not checked)
        if checked:
            self.update_complementary_color()

    def update_complementary_color(self):
        # Convert hex/name to QColor
        col_p = QColor(self.color_p)
        if not col_p.isValid():
            return
            
        h = col_p.hue()
        s = col_p.saturation()
        v = col_p.value()
        
        # Complementary = hue + 180 degrees
        if h != -1: # -1 means grayscale/achromatic
             new_h = (h + 180) % 360
        else:
             new_h = h # Keep achromatic
        
        col_n = QColor.fromHsv(new_h, s, v)
        self.color_n = col_n.name()
        self.btn_color_n.setStyleSheet(f"background-color: {self.color_n}; border: 1px solid gray;")
        self.update_iso()

    def choose_color_p(self):
        c = QColorDialog.getColor(initial=QColor(self.color_p), title="Select Positive Lobe Color")
        if c.isValid():
            self.color_p = c.name()
            self.btn_color_p.setStyleSheet(f"background-color: {self.color_p}; border: 1px solid gray;")
            
            if self.check_comp_color.isChecked():
                self.update_complementary_color()
            else:
                self.update_iso()

    def choose_color_n(self):
        c = QColorDialog.getColor(initial=QColor(self.color_n), title="Select Negative Lobe Color")
        if c.isValid():
            self.color_n = c.name()
            self.btn_color_n.setStyleSheet(f"background-color: {self.color_n}; border: 1px solid gray;")
            self.update_iso()

    def update_iso(self):
        val = self.spin.value()
        # Prefer numeric spinbox value (kept in sync with slider)
        opacity = self.opacity_spin.value() if hasattr(self, 'opacity_spin') else self.opacity_slider.value() / 100.0
        
        try:
            # Cleanup previous
            if self.iso_actor_p: self.plotter.remove_actor(self.iso_actor_p)
            if self.iso_actor_n: self.plotter.remove_actor(self.iso_actor_n)
            
            # Additional safety cleanup by name
            self.plotter.remove_actor("cube_iso_p")
            self.plotter.remove_actor("cube_iso_n")
            
            # Using full grid
            using_grid = self.grid

            # Positive lobe
            iso_p = using_grid.contour(isosurfaces=[val])
            if iso_p.n_points > 0:
                self.iso_actor_p = self.plotter.add_mesh(iso_p, color=self.color_p, opacity=opacity, name="cube_iso_p", reset_camera=False)
                
            # Negative lobe
            iso_n = using_grid.contour(isosurfaces=[-val])
            if iso_n.n_points > 0:
                self.iso_actor_n = self.plotter.add_mesh(iso_n, color=self.color_n, opacity=opacity, name="cube_iso_n", reset_camera=False)
                
            self.plotter.render()
            
        except Exception as e:
            print(f"Iso update error: {e}")
            import traceback
            traceback.print_exc()

    def on_opacity_changed(self, val):
        opacity = val / 100.0
        # Sync numeric spinbox without re-triggering signals
        if hasattr(self, 'opacity_spin'):
            self.opacity_spin.blockSignals(True)
            self.opacity_spin.setValue(opacity)
            self.opacity_spin.blockSignals(False)
        self.update_iso()

    def on_opacity_spin_changed(self, val):
        # val is between 0.0 and 1.0
        # Sync slider (0-100)
        if hasattr(self, 'opacity_slider'):
            int_val = int(round(val * 100))
            self.opacity_slider.blockSignals(True)
            self.opacity_slider.setValue(int_val)
            self.opacity_slider.blockSignals(False)
        self.update_iso()

    def on_slider_changed(self, val):
        float_val = (val / self.slider_max_int) * self.max_val
        self.spin.blockSignals(True)
        self.spin.setValue(float_val)
        self.spin.blockSignals(False)
        self.update_iso()

    def on_spin_changed(self, val):
        if self.max_val > 0:
            int_val = int((val / self.max_val) * self.slider_max_int)
            self.slider.blockSignals(True)
            self.slider.setValue(int_val)
            self.slider.blockSignals(False)
        self.update_iso()

    def close_plugin(self):
        try:
             # Full cleanup
             self.mw.plotter.clear()
             self.mw.current_mol = None
             self.mw.current_file_path = None
             self.mw.plotter.render()
             
             # Restore UI state
             if hasattr(self.mw, 'restore_ui_for_editing'):
                 self.mw.restore_ui_for_editing()
        except: pass
        
        if self.dock:
            self.mw.removeDockWidget(self.dock)
            self.dock.deleteLater()
            self.dock = None
        
        self.deleteLater()

class ChargeDialog(QDialog):
    def __init__(self, parent=None, current_charge=0):
        super().__init__(parent)
        self.setWindowTitle("Bond Connectivity Error")
        self.result_action = "cancel" # retry, skip, cancel
        self.charge = current_charge
        
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Could not determine connectivity with charge: " + str(self.charge)))
        layout.addWidget(QLabel("Please specificy correct charge or skip chemistry check."))
        
        form = QFormLayout()
        self.spin = QSpinBox()
        self.spin.setRange(-20, 20)
        self.spin.setValue(self.charge)
        form.addRow("Charge:", self.spin)
        layout.addLayout(form)
        
        btns = QHBoxLayout()
        retry_btn = QPushButton("Retry")
        retry_btn.clicked.connect(self.on_retry)
        
        skip_btn = QPushButton("Skip Chemistry")
        skip_btn.clicked.connect(self.on_skip)
        
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        
        btns.addWidget(retry_btn)
        btns.addWidget(skip_btn)
        btns.addWidget(cancel_btn)
        layout.addLayout(btns)
        
        self.setLayout(layout)
        
    def on_retry(self):
        self.charge = self.spin.value()
        self.result_action = "retry"
        self.accept()
        
    def on_skip(self):
        self.result_action = "skip"
        self.accept()

def open_cube_viewer(main_window, fname):
    """Core logic to open cube viewer with a specific file."""
    if Chem is None:
        QMessageBox.critical(main_window, "Error", "RDKit is required for this plugin.")
        return

    # Close existing docks
    docks_to_close = []
    for dock in main_window.findChildren(QDockWidget):
        if dock.windowTitle() == "Cube Viewer":
            docks_to_close.append(dock)
    
    for dock in docks_to_close:
        try:
            widget = dock.widget()
            if hasattr(widget, 'close_plugin'):
                widget.close_plugin()
            else:
                main_window.removeDockWidget(dock)
                dock.deleteLater()
        except:
             pass

    try:
        if not fname: # Should not happen if called correctly
            return

        try:
            meta, grid = read_cube(fname)
        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.critical(main_window, "Error", f"Failed to parse Cube file:\n{e}")
            return
        
        if hasattr(main_window, 'plotter'):
            main_window.plotter.clear()
        
        if hasattr(main_window, 'main_window_ui_manager'):
            main_window.main_window_ui_manager._enter_3d_viewer_ui_mode()
        
        # Create Molecule (XYZ)
        atoms = meta['atoms']
        xyz_lines = [f"{len(atoms)}", "Generated by Cube Plugin"]
        pt = Chem.GetPeriodicTable()
        for atomic_num, pos in atoms:
            symbol = pt.GetElementSymbol(atomic_num)
            xyz_lines.append(f"{symbol} {pos[0]:.4f} {pos[1]:.4f} {pos[2]:.4f}")
        
        xyz_content = "\n".join(xyz_lines)
        
        mol = Chem.MolFromXYZBlock(xyz_content)
        
        # Use rdDetermineBonds if available and no bonds found
        current_charge = 0
        if mol is not None and mol.GetNumBonds() == 0 and rdDetermineBonds is not None:
            while True:
                # Re-create fresh molecule for this attempt to avoid dirty state on retry
                mol = Chem.MolFromXYZBlock(xyz_content)
                
                try:
                    rdDetermineBonds.DetermineConnectivity(mol, charge=current_charge)
                    rdDetermineBonds.DetermineBondOrders(mol, charge=current_charge)
                    break # Success
                except Exception as e:
                    # Show Dialog
                    dlg = ChargeDialog(main_window, current_charge)
                    if dlg.exec() == QDialog.DialogCode.Accepted:
                        if dlg.result_action == "retry":
                            current_charge = dlg.charge
                            continue # Loop again with new charge
                        elif dlg.result_action == "skip":
                            # Ensure we have a clean unbonded mol
                            mol = Chem.MolFromXYZBlock(xyz_content)
                            rdDetermineBonds.DetermineConnectivity(mol, charge=0)
                            break # Break loop, accept no bonds/bad connectivity
                    else:
                        return # User Cancelled plugin load
        
        if mol is None:
            # Fallback
            mol = Chem.RWMol()
            conf = Chem.Conformer()
            for i, (atomic_num, pos) in enumerate(atoms):
                idx = mol.AddAtom(Chem.Atom(atomic_num))
                conf.SetAtomPosition(idx, Geometry.Point3D(pos[0], pos[1], pos[2]))
            mol.AddConformer(conf)

        # Set current molecular data in main window for consistency
        main_window.current_mol = mol
        main_window.current_file_path = fname

        # Draw
        if hasattr(main_window, 'draw_molecule_3d'):
             main_window.draw_molecule_3d(mol)
        elif hasattr(main_window, 'main_window_view_3d'):
             main_window.main_window_view_3d.draw_molecule_3d(mol)
             
        # Report bonds
        nb = mol.GetNumBonds() if mol else 0
        if hasattr(main_window, 'statusBar'):
            main_window.statusBar().showMessage(f"Loaded Cube. Atoms: {len(atoms)}, Bonds: {nb}")

        # Setup Dock
        dock = QDockWidget("Cube Viewer", main_window)
        dock.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        
        viewer = CubeViewerWidget(main_window, dock, grid)
        dock.setWidget(viewer)
        
        main_window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        
        main_window.plotter.reset_camera()

    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Plugin Error: {e}")

def run(main_window):
    if Chem is None:
        QMessageBox.critical(main_window, "Error", "RDKit is required for this plugin.")
        return

    fname, _ = QFileDialog.getOpenFileName(main_window, "Open Gaussian Cube File", "", "Cube Files (*.cube *.cub);;All Files (*)")
    if fname:
        open_cube_viewer(main_window, fname)

def autorun(main_window):
    # --- 修正: 現在のメソッドをローカル変数としてキャプチャする ---
    # こうすることで、他のプラグインが設定したラッパーも含めてチェーンできる
    original_dragEnterEvent = main_window.dragEnterEvent
    original_dropEvent = main_window.dropEvent

    def custom_dragEnterEvent(event):
        # Cubeファイルが含まれているかチェック
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                if url.isLocalFile():
                    file_path = url.toLocalFile()
                    if file_path.lower().endswith(('.cube', '.cub')):
                        event.acceptProposedAction()
                        return
        
        # フォールバック: キャプチャしておいたメソッドを呼ぶ
        original_dragEnterEvent(event)

    def custom_dropEvent(event):
        urls = event.mimeData().urls()
        cube_file = None
        for url in urls:
            if url.isLocalFile():
                fpath = url.toLocalFile()
                if fpath.lower().endswith(('.cube', '.cub')):
                    cube_file = fpath
                    break
        
        if cube_file:
            open_cube_viewer(main_window, cube_file)
            event.acceptProposedAction()
        else:
            # フォールバック: キャプチャしておいたメソッドを呼ぶ
            original_dropEvent(event)

    # 上書き適用
    # インスタンスへの関数代入なので、self引数は自動付与されないため
    # custom_dragEnterEvent(event) という定義で正しい
    main_window.dragEnterEvent = custom_dragEnterEvent
    main_window.dropEvent = custom_dropEvent
    
    # 2. Check Command Line Arguments
    import sys
    # Simple check for CLI args
    for arg in sys.argv[1:]:
        if arg.lower().endswith(('.cube', '.cub')) and os.path.exists(arg):
            # Use QTimer to delay slightly to ensure UI is ready
            from PyQt6.QtCore import QTimer
            QTimer.singleShot(100, lambda f=arg: open_cube_viewer(main_window, f))
            break

