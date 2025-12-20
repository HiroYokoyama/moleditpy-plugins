import os
import re
import numpy as np
import traceback
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QListWidget, QSlider, QCheckBox, QFileDialog, QMessageBox,
                             QDockWidget, QWidget, QSplitter, QApplication, QTreeWidget, QTreeWidgetItem, QHeaderView)
from PyQt6.QtCore import Qt, QTimer, pyqtSignal

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
except ImportError:
    Chem = None

PLUGIN_NAME = "Gaussian Freq Analyzer"
__version__="2025.12.20"
__author__="HiroYokoyama"

class FCHKParser:
    def __init__(self):
        self.filename = ""
        self.atoms = [] # List of atomic numbers
        self.coords = [] # List of (x, y, z) in Angstrom
        self.frequencies = [] # List of frequencies in cm^-1
        self.vib_modes = [] # List of displacement vectors (natoms, 3)
        self.charge = 0
        self.multiplicity = 1

    def parse(self, filename):
        self.filename = filename
        self.atoms = []
        self.coords = []
        self.frequencies = []
        self.intensities = [] # IR Intensities
        self.vib_modes = [] # list of list of (dx, dy, dz) 
        
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        data = {}
        current_section = None
        current_values = []
        
        # Helper to process accumulated values
        def process_section(section, values):
            if section == "Atomic numbers":
                # Integers
                return [int(v) for v in values]
            elif section == "Current cartesian coordinates":
                # Floats, Bohr -> Angstrom? No, FCHK usually stores in Bohr.
                # Use standard conversion 0.529177
                return [float(v) for v in values]
            elif section == "Vib-E2":
                return [float(v) for v in values]
            elif section == "Vib-Modes":
                return [float(v) for v in values]
            elif section == "Dipole Derivatives":
                return [float(v) for v in values]
            elif section == "IR Inten":
                return [float(v) for v in values]
            elif section == "Charge":
                return [int(v) for v in values]
            elif section == "Multiplicity":
                return [int(v) for v in values]
            return values

        # Basic FCHK parsing
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line:
                i += 1
                continue
            
            # Section header? 
            # Pattern: Name (starts upper) ... Type (I/R) ... N= (optional) ... Value/Count
            # Example: "Vib-Modes              R   N=          27"
            # Or: "Charge                 I                0"
            # We look for Capital Start, then some space, then I or R, then optional N=
            
            # Simple check: line starts with upper char
            if line[0].isupper():
                # Check for Type and N= signature more loosely but reliably
                # If "   I   " or "   R   " is present (at least 3 spaces before/after or N=)
                # Let's try Regex for robustness
                # Matches: Start with Word, spaces, I or R, spaces, (N= xxxxx)?
                match = re.search(r'^([A-Za-z0-9\-\s]+?)\s+([IR])\s+(?:N=\s+(\d+)|([0-9\-]+))', line)
                
                # Actually, simpler: check for "   I   " or "   R   " or "   I " at specific columns?
                # FCHK is fixed format mostly. 
                # Name: 0-40. Type: 43. 
                # But let's trust "   I" or "   R" presence for now if regex is too complex to get right blindly.
                # However, the previous "R   N=" failed. Maybe it was "R    N=".
                # Let's use a regex that handles variable whitespace.
                if re.search(r'\s+[IR]\s+(N=)?\s+', line):
                     # Store previous
                    if current_section:
                        data[current_section] = process_section(current_section, current_values)
                        
                    # Extract section name
                    parts = line.split()
                    # Name is usually first 1-N tokens until I/R
                    # But Name can have spaces "Atomic numbers".
                    # split by "   I" or "   R" is safest if exists.
                    
                    if "   I" in line:
                         current_section = line.split("   I")[0].strip()
                    elif "   R" in line:
                         current_section = line.split("   R")[0].strip()
                    else:
                         # Fallback: scan tokens for I or R
                         label_parts = []
                         for p in parts:
                             if p == 'I' or p == 'R':
                                 break
                             label_parts.append(p)
                         current_section = " ".join(label_parts)
                    
                    current_values = []
                    i += 1
                    continue
            
            # Data line (if not header)
                
            # Data line
            # Accumulated values
            # FCHK data lines are space separated
            tokens = line.split()
            current_values.extend(tokens)
            i += 1
            
        # Last section
        if current_section:
            data[current_section] = process_section(current_section, current_values)
            
        # Extract specific data
        if "Atomic numbers" in data:
            self.atoms = data["Atomic numbers"]
        
        BOHR_TO_ANG = 0.529177210903
        
        if "Current cartesian coordinates" in data:
            raw_coords = data["Current cartesian coordinates"]
            # Convert to list of tuples (x,y,z)
            # FCHK coords are X1, Y1, Z1, X2... in Bohr
            coords_ang = []
            for j in range(0, len(raw_coords), 3):
                if j+2 < len(raw_coords):
                    x = raw_coords[j] * BOHR_TO_ANG
                    y = raw_coords[j+1] * BOHR_TO_ANG
                    z = raw_coords[j+2] * BOHR_TO_ANG
                    coords_ang.append((x, y, z))
            self.coords = coords_ang
            
        # Parse Vib-E2 section properly
        # Vib-E2 contains blocks: [Frequencies, Reduced Masses, Force Constants, IR Intensities, ...]
        # Each block has n_modes values
        if "Vib-E2" in data:
            raw_e2 = data["Vib-E2"]
            
            # Determine n_modes from Vib-Modes section (safest approach)
            if "Vib-Modes" in data:
                n_atoms = len(self.atoms) if self.atoms else 0
                if n_atoms > 0:
                    # Vib-Modes size = 3 * n_atoms * n_modes
                    n_modes = len(data["Vib-Modes"]) // (3 * n_atoms)
                else:
                    n_modes = 0
            else:
                # Fallback: estimate from 3N-6 (or 3N-5 for linear)
                n_atoms = len(self.atoms)
                n_modes = max(1, 3 * n_atoms - 6) if n_atoms > 2 else 1
            
            if n_modes > 0 and len(raw_e2) >= n_modes:
                # Block 1: Frequencies (0 to n_modes)
                self.frequencies = raw_e2[0:n_modes]
                
                # Block 4: IR Intensities (3*n_modes to 4*n_modes)
                # Already in km/mol units - NO conversion needed!
                if len(raw_e2) >= 4 * n_modes:
                    ir_start = 3 * n_modes
                    ir_end = 4 * n_modes
                    self.intensities = raw_e2[ir_start:ir_end]
        
        # Get actual masses used by Gaussian
        self.masses = []
        if "Real atomic weights" in data:
            self.masses = [float(m) for m in data["Real atomic weights"]]
        elif "Atomic numbers" in data:
            # Fallback: RDKit average atomic weight (slight difference)
            if Chem:
                from rdkit.Chem import GetPeriodicTable
                pt = GetPeriodicTable()
                self.masses = [pt.GetAtomicWeight(int(z)) for z in data["Atomic numbers"]]

        # Fallback: Check for separate IR Inten section (uncommon but possible)
        # Gaussian's conversion factor: Atomic Units -> km/mol
        CONV_FACTOR = 974.868

        # Only look for separate IR section if we didn't get it from Vib-E2
        if not self.intensities or len(self.intensities) == 0:
            ir_key = None
            for key in data.keys():
                if key.strip().lower() in ["ir inten", "ir intensities"]:
                    ir_key = key
                    break
            
            if ir_key:
                raw_vals = [float(val) for val in data[ir_key]]
                # Always convert from a.u. to km/mol using Gaussian's factor
                self.intensities = [v * CONV_FACTOR for v in raw_vals]

        if "Dipole Derivatives" in data:
            self.dipole_derivs = data["Dipole Derivatives"]
            
        if "Vib-Modes" in data:
            # Modes are stored as X1, Y1, Z1... for mode 1, then mode 2...
            # Size should be 3*N_atoms * N_freqs
            raw_modes = data["Vib-Modes"]
            n_atoms = len(self.atoms)
            n_modes = len(self.frequencies)
            mode_len = n_atoms * 3
            
            parsed_modes = []
            for m in range(n_modes):
                start = m * mode_len
                end = start + mode_len
                if end <= len(raw_modes):
                    mode_vec = raw_modes[start:end]
                    # Structure as list of (dx, dy, dz)
                    vecs = []
                    for k in range(0, len(mode_vec), 3):
                        dx = mode_vec[k]
                        dy = mode_vec[k+1]
                        dz = mode_vec[k+2]
                        vecs.append((dx, dy, dz))
                    parsed_modes.append(vecs)
            self.vib_modes = parsed_modes
            
            # Consistency Fix:
            # Sync frequencies, modes, and intensities to avoid mismatches
            if len(self.frequencies) != len(self.vib_modes):
                n_valid = min(len(self.frequencies), len(self.vib_modes))
                self.frequencies = self.frequencies[:n_valid]
                self.vib_modes = self.vib_modes[:n_valid]
                if self.intensities and len(self.intensities) > n_valid:
                    self.intensities = self.intensities[:n_valid]



        if "Charge" in data and len(data["Charge"]) > 0:
            self.charge = data["Charge"][0]
            
        if "Multiplicity" in data and len(data["Multiplicity"]) > 0:
            self.multiplicity = data["Multiplicity"][0]


from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QListWidget, QSlider, QCheckBox, QFileDialog, QMessageBox,
                             QDockWidget, QWidget, QSplitter, QFormLayout, QDialogButtonBox, QSpinBox, QDoubleSpinBox)
from PyQt6.QtGui import QImage, QPainter, QPen, QColor, QFont, QPaintEvent
try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

# ... (Previous imports safely handled by original file or above)


class GaussianFCHKFreqAnalyzer(QWidget):
    def __init__(self, main_window, dock_widget=None):
        super().__init__(main_window)
        self.mw = main_window
        self.dock = dock_widget
        self.setAcceptDrops(True)
        
        self.parser = None
        self.base_mol = None 
        self.timer = QTimer()
        self.timer.timeout.connect(self.animate_frame)
        self.animation_step = 0
        self.is_playing = False
        self.vector_actor = None
        
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Info Label
        self.lbl_info = QLabel("Drop .fchk file here or click Open")
        self.lbl_info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbl_info.setStyleSheet("border: 2px dashed #AAA; padding: 20px; color: #555;")
        layout.addWidget(self.lbl_info)
        
        # Open Button
        btn_open = QPushButton("Open FCHK File")
        btn_open.clicked.connect(self.open_file_dialog)
        layout.addWidget(btn_open)
        
        # Frequency List
        # Frequency List
        layout.addWidget(QLabel("Vibrational Frequencies:"))
        self.list_freq = QTreeWidget()
        self.list_freq.setColumnCount(3)
        self.list_freq.setHeaderLabels(["No.", "Frequency (cm⁻¹)", "Intensity (km/mol)"])
        self.list_freq.header().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.list_freq.header().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.list_freq.headerItem().setTextAlignment(0, Qt.AlignmentFlag.AlignCenter)  # Center No. header
        self.list_freq.headerItem().setTextAlignment(1, Qt.AlignmentFlag.AlignCenter)  # Center Frequency header
        self.list_freq.headerItem().setTextAlignment(2, Qt.AlignmentFlag.AlignCenter)  # Center Intensity header
        self.list_freq.currentItemChanged.connect(self.on_freq_selected)
        layout.addWidget(self.list_freq)
        
        # Spectrum Button
        btn_spectrum = QPushButton("Show Spectrum")
        btn_spectrum.clicked.connect(self.show_spectrum)
        layout.addWidget(btn_spectrum)

        # Animation Controls
        anim_layout = QVBoxLayout()
        
        # Vector Controls
        vec_layout = QHBoxLayout()
        self.chk_vectors = QCheckBox("Show Vectors")
        self.chk_vectors.setChecked(False)
        self.chk_vectors.stateChanged.connect(lambda state: self.update_vectors())
        vec_layout.addWidget(self.chk_vectors)
        
        vec_layout.addWidget(QLabel("Scale:"))
        self.spin_vec_scale = QDoubleSpinBox()
        self.spin_vec_scale.setRange(0.1, 200.0)
        self.spin_vec_scale.setSingleStep(1.0)
        self.spin_vec_scale.setValue(2.0)
        self.spin_vec_scale.valueChanged.connect(lambda val: self.update_vectors())
        vec_layout.addWidget(self.spin_vec_scale)
        vec_layout.addStretch()
        
        anim_layout.addLayout(vec_layout)
        
        # Amplitude
        amp_layout = QHBoxLayout()
        amp_layout.addWidget(QLabel("Amplitude:"))
        self.slider_amp = QSlider(Qt.Orientation.Horizontal)
        self.slider_amp.setRange(1, 20)
        self.slider_amp.setValue(5)
        
        self.lbl_amp_val = QLabel("5")
        self.slider_amp.valueChanged.connect(lambda v: self.lbl_amp_val.setText(str(v)))
        
        amp_layout.addWidget(self.slider_amp)
        amp_layout.addWidget(self.lbl_amp_val)
        anim_layout.addLayout(amp_layout)
        
        # FPS (Speed)
        speed_layout = QHBoxLayout()
        speed_layout.addWidget(QLabel("FPS:"))
        self.slider_speed = QSlider(Qt.Orientation.Horizontal)
        self.slider_speed.setRange(1, 60)
        self.slider_speed.setValue(20) 
        
        self.lbl_speed_val = QLabel("20")
        self.slider_speed.valueChanged.connect(lambda v: self.lbl_speed_val.setText(str(v)))
        
        self.slider_speed.valueChanged.connect(self.update_timer_interval)
        speed_layout.addWidget(self.slider_speed)
        speed_layout.addWidget(self.lbl_speed_val)
        anim_layout.addLayout(speed_layout)
        
        # Buttons
        btn_layout = QHBoxLayout()
        self.btn_play = QPushButton("Play")
        self.btn_play.clicked.connect(self.toggle_play)
        self.btn_stop = QPushButton("Stop")
        self.btn_stop.clicked.connect(self.stop_play)
        
        btn_layout.addWidget(self.btn_play)
        btn_layout.addWidget(self.btn_stop)
        anim_layout.addLayout(btn_layout)
        
        layout.addLayout(anim_layout)
        
        # Export GIF Button
        self.btn_gif = QPushButton("Export GIF")
        self.btn_gif.clicked.connect(self.save_as_gif)
        self.btn_gif.setEnabled(HAS_PIL)
        layout.addWidget(self.btn_gif)
        
        # Metadata Info
        self.lbl_meta = QLabel("")
        layout.addWidget(self.lbl_meta)
        
        layout.addStretch()
        

        # Close Button
        btn_close = QPushButton("Close")
        def close_action():
            self.stop_play()
            self.reset_geometry()
            self.remove_vectors()
            if self.dock:
                self.dock.close()
            else:
                self.close()
        btn_close.clicked.connect(close_action)
        layout.addWidget(btn_close)
        
        self.setLayout(layout)
        
    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                fname = url.toLocalFile().lower()
                if fname.endswith(".fchk") or fname.endswith(".fck"):
                    event.acceptProposedAction()
                    return
        event.ignore()
        
    def dropEvent(self, event):
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            if file_path.lower().endswith((".fchk", ".fck")):
                self.load_file(file_path)
                event.acceptProposedAction()
                break
                
    def open_file_dialog(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Open Gaussian FCHK", "", "FCHK Files (*.fchk *.fck)")
        if fname:
            self.load_file(fname)
            
    def load_file(self, filename):
        self.parser = FCHKParser()
        try:
            self.parser.parse(filename)
            self.update_ui_after_load()
            self.lbl_info.setText(os.path.basename(filename))
            self.lbl_info.setStyleSheet("border: 2px solid #4CAF50; padding: 10px; color: #4CAF50;")
            
            # Update Main Window Context
            if hasattr(self.mw, 'current_file_path'):
                self.mw.current_file_path = filename
                if hasattr(self.mw, 'update_window_title'):
                    self.mw.update_window_title()
                else:
                    self.mw.setWindowTitle(f"{os.path.basename(filename)} - MoleditPy")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to parse FCHK:\n{e}")
            traceback.print_exc()

    def update_ui_after_load(self):
        self.list_freq.clear()
        # Ensure freqs match modes logic handled in parser consistency check
        if hasattr(self.parser, 'frequencies'):
            for i, freq in enumerate(self.parser.frequencies):
                item = QTreeWidgetItem()
                item.setText(0, str(i + 1))  # Mode number
                item.setText(1, f"{freq:.2f}")
                
                if hasattr(self.parser, 'intensities') and i < len(self.parser.intensities):
                    inten = self.parser.intensities[i]
                    # Higher precision to see small values and compare with LOG
                    item.setText(2, f"{inten:.4f}")
                else:
                    item.setText(2, "-")
                item.setTextAlignment(0, Qt.AlignmentFlag.AlignCenter)  # Center mode number
                item.setTextAlignment(1, Qt.AlignmentFlag.AlignCenter)  # Center frequency
                item.setTextAlignment(2, Qt.AlignmentFlag.AlignCenter)  # Center intensity
                self.list_freq.addTopLevelItem(item)
            
        self.lbl_meta.setText(f"Charge: {self.parser.charge}, Multiplicity: {self.parser.multiplicity}\nAtoms: {len(self.parser.atoms)}")
        
        # Load molecule into main window
        if len(self.parser.atoms) > 0 and Chem:
            self.create_base_molecule()

    def create_base_molecule(self):
        if not self.parser: return
        
        mol = Chem.RWMol()
        
        # Mapping atomic number to symbol?
        # Needed: FCHK gives Atomic Numbers
        pt = Chem.GetPeriodicTable()
        
        for ans in self.parser.atoms:
            sym = pt.GetElementSymbol(int(ans))
            atom = Chem.Atom(sym)
            mol.AddAtom(atom)
            
        conf = Chem.Conformer(len(self.parser.atoms))
        for idx, (x, y, z) in enumerate(self.parser.coords):
            conf.SetAtomPosition(idx, Point3D(x, y, z))
        mol.AddConformer(conf)
        
        if hasattr(self.mw, 'estimate_bonds_from_distances'):
            self.mw.estimate_bonds_from_distances(mol)
            
        self.base_mol = mol.GetMol()
        self.mw.current_mol = self.base_mol
        
        if hasattr(self.mw, '_enter_3d_viewer_ui_mode'):
            self.mw._enter_3d_viewer_ui_mode()
            
        self.mw.draw_molecule_3d(self.base_mol)
        if hasattr(self.mw, 'plotter'):
            self.mw.plotter.reset_camera()

    def on_freq_selected(self, current, previous):
        if self.is_playing:
            # Transition smoothly to new selected mode
            pass
        else:
            self.update_vectors()

    def toggle_play(self):
        curr = self.list_freq.currentItem()
        if not curr:
            return
        # Row index logic
        row = self.list_freq.indexOfTopLevelItem(curr)
        if row < 0: return
            
        if self.is_playing:
            # Pause logic
            self.is_playing = False
            self.timer.stop()
            self.btn_play.setText("Play")
            # Do NOT reset geometry
            return
            
        self.is_playing = True
        self.btn_play.setText("Pause")
        self.timer.start(50) 
        self.update_timer_interval()
        
    def stop_play(self):
        self.is_playing = False
        self.timer.stop()
        self.btn_play.setText("Play")

        self.reset_geometry()
        QApplication.processEvents()

    def update_timer_interval(self):
        fps = self.slider_speed.value()
        if fps <= 0: fps = 1
        interval = 1000 / fps
        self.timer.setInterval(int(interval))

    def reset_geometry(self):
        if not self.base_mol or not self.parser: return
        conf = self.base_mol.GetConformer()
        for idx, (x, y, z) in enumerate(self.parser.coords):
            conf.SetAtomPosition(idx, Point3D(x, y, z))
        self.mw.draw_molecule_3d(self.base_mol)
        
        self.update_vectors()
        
        if hasattr(self.mw, 'plotter'):
            self.mw.plotter.render()

    def animate_frame(self):
        if not self.parser or not self.base_mol:
            self.stop_play()
            return
            
        curr = self.list_freq.currentItem()
        if not curr: return
        row = self.list_freq.indexOfTopLevelItem(curr)
        if row < 0 or row >= len(self.parser.vib_modes):
            return
            
        mode_vecs = self.parser.vib_modes[row]
        
        self.animation_step += 1
        # Use 20 steps per cycle
        cycle_pos = (self.animation_step % 20) / 20.0
        phase = cycle_pos * 2 * np.pi
        
        scale = self.slider_amp.value() / 20.0 
        factor = np.sin(phase) * scale
        
        self.apply_displacement(mode_vecs, factor)
        self.mw.draw_molecule_3d(self.base_mol)
        self.update_vectors(mode_vecs=mode_vecs, scale_factor=factor)
        
    def apply_displacement(self, mode_vecs, factor):
        conf = self.base_mol.GetConformer()
        base_coords = self.parser.coords
        for idx, (bx, by, bz) in enumerate(base_coords):
            if idx < len(mode_vecs):
                dx, dy, dz = mode_vecs[idx]
                nx = bx + dx * factor
                ny = by + dy * factor
                nz = bz + dz * factor
                conf.SetAtomPosition(idx, Point3D(nx, ny, nz))

    def remove_vectors(self):
        if self.vector_actor and hasattr(self.mw, 'plotter'):
            try:
                self.mw.plotter.remove_actor(self.vector_actor)
            except: pass
        self.vector_actor = None

    def update_vectors(self, mode_vecs=None, scale_factor=0.0):
        # Clean up existing vectors
        self.remove_vectors()
        
        if not self.chk_vectors.isChecked():
            return
            
        if not self.parser or not self.base_mol or not hasattr(self.mw, 'plotter'):
            return

        # Get current frequency
        curr = self.list_freq.currentItem()
        if not curr: return
        row = self.list_freq.indexOfTopLevelItem(curr)
        if row < 0 or row >= len(self.parser.vib_modes): return
        
        # Get vectors if not provided
        if mode_vecs is None:
            mode_vecs = self.parser.vib_modes[row]
            
        # Current coords from molecule conformer
        conf = self.base_mol.GetConformer()
        coords = []
        vectors = []
        
        # Amplitude for vector length scaling
        # Now decoupled from animation amplitude
        vis_scale = self.spin_vec_scale.value()
        
        for idx in range(len(mode_vecs)):
             pos = conf.GetAtomPosition(idx)
             coords.append([pos.x, pos.y, pos.z])
             
             dx, dy, dz = mode_vecs[idx]
             vectors.append([dx, dy, dz])
             
        if not coords: return
        
        coords = np.array(coords)
        vectors = np.array(vectors)
        
        try:
            self.vector_actor = self.mw.plotter.add_arrows(coords, vectors, mag=vis_scale, color='lightgreen', show_scalar_bar=False)
        except Exception as e:
            print(f"Error adding arrows: {e}")

    def save_as_gif(self):
        if not self.parser or not self.base_mol: return
        
        # Pause to configure
        was_playing = self.is_playing
        if self.is_playing:
            self.toggle_play() # Pause

        curr = self.list_freq.currentItem()
        if not curr:
            QMessageBox.warning(self, "Select Frequency", "Please select a frequency to export.")
            return
        row = self.list_freq.indexOfTopLevelItem(curr)

        dialog = QDialog(self)
        dialog.setWindowTitle("Export GIF Settings")
        form = QFormLayout(dialog)
        
        # Calculate current FPS
        # Slider value is now FPS directly
        current_fps = self.slider_speed.value()
        
        spin_fps = QSpinBox()
        spin_fps.setRange(1, 60)
        spin_fps.setValue(current_fps)
        
        chk_transparent = QCheckBox()
        chk_transparent.setChecked(True)
        
        chk_hq = QCheckBox()
        chk_hq.setChecked(True)
        
        form.addRow("FPS:", spin_fps)
        form.addRow("Transparent Background:", chk_transparent)
        form.addRow("High Quality (Adaptive):", chk_hq)
        
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        form.addRow(btns)
        
        if dialog.exec() != QDialog.DialogCode.Accepted:
            if was_playing: self.toggle_play()
            return # Cancel

        target_fps = spin_fps.value()
        use_transparent = chk_transparent.isChecked()
        use_hq = chk_hq.isChecked()
        
        file_path, _ = QFileDialog.getSaveFileName(self, "Save GIF", "", "GIF Files (*.gif)")
        if not file_path:
            if was_playing: self.toggle_play()
            return
            
        if not file_path.lower().endswith('.gif'):
             file_path += '.gif'
             
        # Generate Frames
        # 1 Cycle = 20 steps
        images = []
        mode_vecs = self.parser.vib_modes[row]
        
        import copy
        # Store current geometry to restore later
        self.reset_geometry() # align to base
        
        try:
            for i in range(20):
                cycle_pos = i / 20.0
                phase = cycle_pos * 2 * np.pi
                scale = self.slider_amp.value() / 20.0
                factor = np.sin(phase) * scale # Calculate factor here
                self.apply_displacement(mode_vecs, factor)
                self.mw.draw_molecule_3d(self.base_mol)
                self.update_vectors(mode_vecs, factor)
                self.mw.plotter.render()
                
                img_array = self.mw.plotter.screenshot(transparent_background=use_transparent, return_img=True)
                if img_array is not None:
                     img = Image.fromarray(img_array)
                     images.append(img)
            
            if images:
                duration_ms = int(1000 / target_fps)
                
                processed_images = []
                for img in images:
                    if use_hq:
                         if use_transparent:
                            # Alpha preservation with adaptive palette wrapper
                            alpha = img.split()[3]
                            img_rgb = img.convert("RGB")
                            # Quantize to 255 colors to leave room for transparency
                            img_p = img_rgb.convert('P', palette=Image.Palette.ADAPTIVE, colors=255)
                            # Set simple transparency
                            mask = Image.eval(alpha, lambda a: 255 if a <= 128 else 0)
                            img_p.paste(255, mask)
                            img_p.info['transparency'] = 255
                            processed_images.append(img_p)
                         else:
                            processed_images.append(img.convert("P", palette=Image.Palette.ADAPTIVE, colors=256))
                    else:
                        if use_transparent:
                            img = img.convert("RGBA")
                            processed_images.append(img)
                        else:
                            processed_images.append(img.convert("RGB"))

                processed_images[0].save(file_path, save_all=True, append_images=processed_images[1:], duration=duration_ms, loop=0, disposal=2)
                QMessageBox.information(self, "Success", f"Saved GIF to:\n{file_path}")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save GIF: {e}")
            traceback.print_exc()
        finally:
            self.reset_geometry()
            # Restore play state if needed, or leave paused
            # User might want to inspect
            if was_playing:
                 self.toggle_play()
    
    def on_dock_visibility_changed(self, visible):
        if not visible and self.is_playing:
            self.stop_play()

    def close_plugin(self):
        self.stop_play()
        self.remove_vectors()
        self.animation_step = 0
    def show_spectrum(self):
        if not self.parser or not self.parser.frequencies:
            return
        
        # Filter out low frequencies (translations/rotations)
        freqs = []
        intensities = []
        
        parser_intensities = self.parser.intensities if hasattr(self.parser, 'intensities') and self.parser.intensities else [1.0]*len(self.parser.frequencies)
        
        for i, freq in enumerate(self.parser.frequencies):
            # Use abs() to preserve imaginary frequencies (negative values)
            # Only exclude low-frequency modes (translations/rotations)
            if abs(freq) > 10.0:
                freqs.append(freq)
                if i < len(parser_intensities):
                    intensities.append(parser_intensities[i])
                else:
                    intensities.append(1.0)
        
        dlg = SpectrumDialog(freqs, intensities, self)
        dlg.exec()


class SpectrumDialog(QDialog):
    def __init__(self, freqs, intensities, parent=None):
        super().__init__(parent)
        self.setWindowTitle("IR Spectrum")
        self.resize(800, 600)
        
        self.freqs = np.array(freqs)
        self.intensities = np.array(intensities)
        self.scaling_factor = 1.0
        
        # Layout
        layout = QVBoxLayout(self)
        
        # Plot Area
        self.plot_widget = SpectrumPlotWidget(self.freqs, self.intensities)
        layout.addWidget(self.plot_widget)
        
        # Controls
        controls = QHBoxLayout()
        
        # Scaling Factor
        controls.addWidget(QLabel("Scaling Factor:"))
        from PyQt6.QtWidgets import QDoubleSpinBox
        self.spin_scale = QDoubleSpinBox()
        self.spin_scale.setRange(0.5, 1.5)
        self.spin_scale.setValue(1.0)
        self.spin_scale.setSingleStep(0.01)
        self.spin_scale.setDecimals(3)
        self.spin_scale.valueChanged.connect(self.on_scaling_changed)
        controls.addWidget(self.spin_scale)
        
        controls.addWidget(QLabel("FWHM (cm⁻¹):"))
        self.spin_fwhm = QSpinBox()
        self.spin_fwhm.setRange(1, 500)
        self.spin_fwhm.setValue(50)
        self.spin_fwhm.valueChanged.connect(self.on_fwhm_changed)
        controls.addWidget(self.spin_fwhm)
        
        # Axis Range
        controls.addWidget(QLabel("Min WN:"))
        self.spin_min = QSpinBox()
        self.spin_min.setRange(0, 5000)
        self.spin_min.setValue(0)
        self.spin_min.setSingleStep(100)
        self.spin_min.valueChanged.connect(self.on_range_changed)
        controls.addWidget(self.spin_min)

        controls.addWidget(QLabel("Max WN:"))
        self.spin_max = QSpinBox()
        self.spin_max.setRange(0, 5000)
        self.spin_max.setValue(4000)
        self.spin_max.setSingleStep(100)
        self.spin_max.valueChanged.connect(self.on_range_changed)
        controls.addWidget(self.spin_max)
        
        btn_csv = QPushButton("Export CSV")
        btn_csv.clicked.connect(self.export_csv)
        controls.addWidget(btn_csv)
        
        btn_png = QPushButton("Export Image")
        btn_png.clicked.connect(self.export_png)
        controls.addWidget(btn_png)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        controls.addWidget(btn_close)
        
        layout.addLayout(controls)
        self.setLayout(layout)
        
        # Initial Plot
        self.on_range_changed()

    def on_scaling_changed(self, val):
        self.scaling_factor = val
        scaled_freqs = self.freqs * self.scaling_factor
        self.plot_widget.set_frequencies(scaled_freqs)
    
    def on_fwhm_changed(self, val):
        self.plot_widget.set_fwhm(val)

    def on_range_changed(self):
        mn = self.spin_min.value()
        mx = self.spin_max.value()
        if mx > mn:
             self.plot_widget.set_range(mn, mx)

    def export_csv(self):
        fname, _ = QFileDialog.getSaveFileName(self, "Save Spectrum Data", "", "CSV Files (*.csv)")
        if fname:
            if not fname.lower().endswith('.csv'): fname += '.csv'
            try:
                x, y = self.plot_widget.get_curve_data()
                with open(fname, 'w') as f:
                    f.write("Frequency,Intensity\n")
                    for xi, yi in zip(x, y):
                        f.write(f"{xi:.2f},{yi:.4f}\n")
                QMessageBox.information(self, "Success", "Saved CSV.")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def export_png(self):
        fname, _ = QFileDialog.getSaveFileName(self, "Save Spectrum Image", "", "PNG Files (*.png)")
        if fname:
            if not fname.lower().endswith('.png'): fname += '.png'
            try:
                # Capture the widget
                pixmap = self.plot_widget.grab()
                pixmap.save(fname)
                QMessageBox.information(self, "Success", "Saved Image.")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

class SpectrumPlotWidget(QWidget):
    def __init__(self, freqs, intensities, parent=None):
        super().__init__(parent)
        self.freqs = freqs
        self.intensities = intensities
        self.fwhm = 80.0
        self.curve_x = []
        self.curve_y = []

        self.setAutoFillBackground(True)
        self.setStyleSheet("background-color: white;")
        
        self.min_x = 0.0
        self.max_x = 4000.0

    def set_fwhm(self, val):
        self.fwhm = val
        self.recalc_curve()
        self.update()

    def set_frequencies(self, freqs):
        self.freqs = freqs
        self.recalc_curve()
        self.update()

    def set_range(self, mn, mx):
        self.min_x = float(mn)
        self.max_x = float(mx)
        self.recalc_curve()
        self.update()

    def get_curve_data(self):
        return self.curve_x, self.curve_y

    def recalc_curve(self):
        # Determine range
        if len(self.freqs) == 0: return
        
        # X resolution based on custom range
        self.curve_x = np.linspace(self.min_x, self.max_x, 1000)
        self.curve_y = np.zeros_like(self.curve_x)
        
        # Sum Gaussians weighted by intensity
        sigma = self.fwhm / 2.35482 
        
        for f, i in zip(self.freqs, self.intensities):
            self.curve_y += i * np.exp(-(self.curve_x - f)**2 / (2 * sigma**2))
            
        # Do NOT normalize - preserve actual intensity values


    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        w = self.width()
        h = self.height()
        
        # Margins
        margin_l = 50
        margin_r = 20
        margin_t = 20
        margin_b = 60  # Increased from 40 for better spacing
        
        plot_w = w - margin_l - margin_r
        plot_h = h - margin_t - margin_b
        
        if len(self.curve_x) == 0:
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No Data")
            return
            
        # Calculate max_y from both curve AND raw intensities to prevent stick normalization
        # Add 10% padding for better visibility
        max_curve = np.max(self.curve_y) if len(self.curve_y) > 0 else 1.0
        max_intensity = np.max(self.intensities) if len(self.intensities) > 0 else 1.0
        max_y = max(max_curve, max_intensity) * 1.1  # 1.1x for padding
        if max_y == 0: max_y = 1.0
        min_x = np.min(self.curve_x)
        max_x = np.max(self.curve_x)
        range_x = max_x - min_x
        if range_x == 0: range_x = 100
        
        # Helper to transform
        # Inverted X: Max (left) -> Min (right)
        # Inverted Y: 0 (top) -> Max (bottom) - peaks point down
        
        def to_screen(x, y):
            # X: margin_l corresponds to max_x, w - margin_r corresponds to min_x
            sx = margin_l + (max_x - x) / range_x * plot_w
            
            # Y: margin_t corresponds to 0, h - margin_b corresponds to max_y
            sy = margin_t + (y / max_y) * plot_h
            return sx, sy

        # Draw Axes
        painter.setPen(QPen(Qt.GlobalColor.black, 2))
        painter.drawLine(margin_l, margin_t, w - margin_r, margin_t) # X-axis at top (Baseline)
        painter.drawLine(margin_l, h - margin_b, w - margin_r, h - margin_b) # X-axis at bottom
        
        painter.drawLine(margin_l, h - margin_b, margin_l, margin_t) # Y (Left)
        painter.drawLine(w - margin_r, h - margin_b, w - margin_r, margin_t) # Y (Right border)

        # Draw Ticks / Labels (Simplified)
        font = painter.font()
        font.setPointSize(12)  # Increased from 8
        painter.setFont(font)
        
        # X Ticks (approx 5)
        # Inverted: Left is Max, Right is Min
        n_ticks = 5
        for i in range(n_ticks + 1):
            # val goes from max_x to min_x
            val = max_x - (range_x * i / n_ticks)
            px, py = to_screen(val, 0)
            # Label at bottom
            painter.drawText(int(px)-20, h - margin_b + 5, 40, 20, Qt.AlignmentFlag.AlignCenter, f"{int(val)}")
            painter.drawLine(int(px), h - margin_b, int(px), h - margin_b + 5)

        # Labels
        font.setPointSize(14)  # Increased from 10
        font.setBold(True)
        painter.setFont(font)
        painter.drawText(0, h-25, w, 20, Qt.AlignmentFlag.AlignCenter, "Wavenumber (cm⁻¹)")
        
        # Draw baseline at y=0
        baseline_x_start, baseline_y = to_screen(max_x, 0)
        baseline_x_end, _ = to_screen(min_x, 0)
        painter.setPen(QPen(QColor(150, 150, 150), 1, Qt.PenStyle.DashLine))
        painter.drawLine(int(baseline_x_start), int(baseline_y), int(baseline_x_end), int(baseline_y))
        
        # Draw Curve
        painter.setPen(QPen(Qt.GlobalColor.blue, 2))
        path_points = []
        for x, y in zip(self.curve_x, self.curve_y):
            sx, sy = to_screen(x, y)
            path_points.append( (sx, sy) )
            
        if len(path_points) > 1:
            from PyQt6.QtCore import QPointF
            qpoints = [QPointF(x, y) for x, y in path_points]
            painter.drawPolyline(qpoints)
            
        # Draw Sticks (Bars) for original frequencies
        painter.setPen(QPen(QColor(255, 0, 0, 100), 1))
        for f, i in zip(self.freqs, self.intensities):
             sx, sy = to_screen(f, i)
             px_base, py_base = to_screen(f, 0)
             painter.drawLine(int(sx), int(py_base), int(sx), int(sy))
        # If we need to remove self from layout or close dock?
        # Done by caller usually.

    # def closeEvent(self, event):
    #     self.close_plugin()
    #     super().closeEvent(event)


    
    def on_dock_visibility_changed(self, visible):
        if not visible:
            self.stop_play()
            self.animation_step = 0
            if self.base_mol:
                self.reset_geometry()
                # Force redraw
                if hasattr(self.mw, 'plotter'):
                    self.mw.plotter.render()

def load_from_file(main_window, fname):
    # Check for existing dock
    dock = None
    analyzer = None
    
    # Check existing docks
    for d in main_window.findChildren(QDockWidget):
        if d.windowTitle() == "Gaussian Freq Analyzer":
            dock = d
            analyzer = d.widget()
            break
            
    if not dock:
        dock = QDockWidget("Gaussian Freq Analyzer", main_window)
        dock.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        analyzer = GaussianFCHKFreqAnalyzer(main_window, dock)
        dock.setWidget(analyzer)
        main_window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        
        # Connect visibility change
        dock.visibilityChanged.connect(analyzer.on_dock_visibility_changed)
    
    dock.show()
    dock.raise_()
    
    if analyzer:
        analyzer.load_file(fname)

def run(main_window):
    # Smart Open Logic
    if hasattr(main_window, 'current_file_path') and main_window.current_file_path:
        fpath = main_window.current_file_path.lower()
        if fpath.endswith((".fchk", ".fck")):
             load_from_file(main_window, main_window.current_file_path)
             return

    # Just show empty or open dialog?
    # Logic: Open File Dialog first
    fname, _ = QFileDialog.getOpenFileName(main_window, "Open Gaussian FCHK", "", "Gaussian FCHK (*.fchk *.fck);;All Files (*)")
    if not fname:
        return
    load_from_file(main_window, fname)

def autorun(main_window):
    try:
        # Check if already patched
        if not hasattr(main_window, '_custom_drop_handlers'):
            main_window._custom_drop_handlers = {}
            
            main_window._original_dragEnterEvent = main_window.dragEnterEvent
            main_window._original_dropEvent = main_window.dropEvent
            
            # --- Monkey patch for Command Line Argument Opening ---
            if not hasattr(main_window, '_custom_loader_handlers'):
                main_window._custom_loader_handlers = {}
                main_window._original_load_command_line_file = main_window.load_command_line_file

                def dynamic_load_command_line_file(self, file_path):
                    if not file_path: return
                    
                    lower_path = file_path.lower()
                    handled = False
                    
                    # Check custom handlers (cli handlers)
                    # Actually we can reuse drop handlers map if we want, but keeping separate is cleaner?
                    # The plan said "register handler". 
                    # Let's check _custom_loader_handlers. 
                    # Wait, below I see usage of _custom_loader_handlers.
                    if hasattr(self, '_custom_loader_handlers'):
                        for ext, handler in self._custom_loader_handlers.items():
                            if lower_path.endswith(ext):
                                try:
                                    handler(self, file_path)
                                    handled = True
                                    break
                                except Exception as e:
                                    print(f"Error loading {file_path} via plugin: {e}")
                    
                    if not handled:
                        # Fallback to original
                        if hasattr(self, '_original_load_command_line_file'):
                            self._original_load_command_line_file(file_path)

                import types
                main_window.load_command_line_file = types.MethodType(dynamic_load_command_line_file, main_window)

            
            
            def dynamic_dragEnterEvent(self, event):
                has_handler = False
                if event.mimeData().hasUrls():
                    for url in event.mimeData().urls():
                        if url.isLocalFile():
                            fpath = url.toLocalFile().lower()
                            for ext in self._custom_drop_handlers:
                                if fpath.endswith(ext):
                                    has_handler = True
                                    break
                        if has_handler: break
                
                if has_handler:
                    event.acceptProposedAction()
                else:
                    if hasattr(self, '_original_dragEnterEvent'):
                        self._original_dragEnterEvent(event)
            
            def dynamic_dropEvent(self, event):
                handled = False
                if event.mimeData().hasUrls():
                    for url in event.mimeData().urls():
                        if url.isLocalFile():
                            fpath = url.toLocalFile()
                            lower_path = fpath.lower()
                            for ext, handler in self._custom_drop_handlers.items():
                                if lower_path.endswith(ext):
                                    try:
                                        handler(self, fpath)
                                        handled = True
                                        event.acceptProposedAction()
                                    except Exception as e:
                                        print(f"Error in drop handler for {ext}: {e}")
                                    break
                        if handled: break
                
                if not handled:
                    if hasattr(self, '_original_dropEvent'):
                        self._original_dropEvent(event)

            import types
            main_window.dragEnterEvent = types.MethodType(dynamic_dragEnterEvent, main_window)
            main_window.dropEvent = types.MethodType(dynamic_dropEvent, main_window)
            
        main_window._custom_drop_handlers['.fchk'] = load_from_file
        main_window._custom_drop_handlers['.fck'] = load_from_file
        # Register CLI handler too
        if hasattr(main_window, '_custom_loader_handlers'):
             main_window._custom_loader_handlers['.fchk'] = load_from_file
             main_window._custom_loader_handlers['.fck'] = load_from_file
        
    except Exception as e:
        print(f"Failed to register FCHK DnD: {e}")
