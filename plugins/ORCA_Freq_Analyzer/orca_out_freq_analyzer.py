import os
import numpy as np
import traceback
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QListWidget, QSlider, QCheckBox, QFileDialog, QMessageBox,
                             QDockWidget, QWidget, QFormLayout, QDialogButtonBox, QSpinBox, QApplication, QTreeWidget, QTreeWidgetItem, QHeaderView, QDoubleSpinBox)
from PyQt6.QtGui import QImage, QPainter, QPen, QColor, QFont, QPaintEvent
try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
from PyQt6.QtCore import Qt, QTimer, pyqtSignal

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
except ImportError:
    Chem = None

PLUGIN_NAME = "ORCA Freq Analyzer"
__version__="2025.12.20"
__author__="HiroYokoyama"

class OrcaParser:
    def __init__(self):
        self.filename = ""
        self.atoms = []
        self.coords = []
        self.frequencies = []
        self.intensities = [] # IR Intensities
        self.vib_modes = [] # list of list of (dx, dy, dz)
        self.charge = 0
        self.multiplicity = 1

    def parse(self, filename):
        self.filename = filename
        self.atoms = []
        self.coords = []
        self.frequencies = []
        self.vib_modes = []
        self.charge = 0
        self.multiplicity = 1
        
        with open(filename, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
            
        # 1. basic properties: charge, mult
        # Look for "Total Charge           Charge          ....   0"
        # OR In input block: "* xyz   0   1"
        # scan lines
        
        # 2. Geometry: Look for "CARTESIAN COORDINATES (ANGSTROEM)" or similar.
        # Use the LAST occurrence to get optimized geometry.
        
        coord_start_lines = []
        freq_start_line = -1
        modes_start_line = -1
        
        for i, line in enumerate(lines):
            line_s = line.strip()
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                coord_start_lines.append(i)
            elif "CARTESIAN COORDINATES (A.U.)" in line:
                pass # ignore AU for now unless needed
            elif "VIBRATIONAL FREQUENCIES" in line:
                freq_start_line = i
            elif "NORMAL MODES" in line:
                modes_start_line = i
            elif "Total Charge" in line and "Charge" in line:
                 # e.g. "Total Charge           Charge          ....   0"
                 parts = line.split()
                 try:
                     self.charge = int(parts[-1])
                 except: pass
            elif "Multiplicity" in line and "Mult" in line:
                 parts = line.split()
                 try:
                     self.multiplicity = int(parts[-1])
                 except: pass
        
        # Parse Geometry (Last one)
        if coord_start_lines:
            start = coord_start_lines[-1]
            try:
                # format:
                # CARTESIAN COORDINATES (ANGSTROEM)
                # ---------------------------------
                #   O      0.000000    0.000000    0.000000
                #   H      0.000000    0.759337    0.596043
                
                # Check where data starts. Usually start+2
                curr = start + 2
                while curr < len(lines):
                    l = lines[curr].strip()
                    if not l:
                        curr += 1
                        # If multiple empty lines, maybe end
                        if curr < len(lines) and not lines[curr].strip():
                            break
                        continue
                        
                    parts = l.split()
                    if len(parts) >= 4:
                        sym = parts[0]
                        # Check if sym is element
                        if not sym[0].isalpha(): 
                             break # End of block
                             
                        try:
                            x = float(parts[1])
                            y = float(parts[2])
                            z = float(parts[3])
                            
                            # Valid atom
                            # Convert Sym to atomic num? RDKit needs num or Valid symbol
                            # We can keep symbol or convert
                            # Let's keep symbol for now, convert to num for uniformity if needed
                            # RDKit Atom(symbol) works
                            self.atoms.append(sym)
                            self.coords.append((x, y, z))
                        except ValueError:
                            pass
                    else:
                        break
                    curr += 1
            except Exception as e:
                print(f"Error parsing coords: {e}")

        # Parse Frequencies
        # -----------------------
        # VIBRATIONAL FREQUENCIES
        # -----------------------
        # ...
        #    0:       0.00 cm**-1
        #    6:    1709.03 cm**-1
        
        if freq_start_line > 0:
            curr = freq_start_line + 4 # skip header
            while curr < len(lines):
                l = lines[curr].strip()
                if not l: 
                    curr += 1
                    continue
                if "NORMAL MODES" in l: # safe guard
                    break
                # Format: "   6:    1709.03 cm**-1"
                if ":" in l and "cm**-1" in l:
                    parts = l.split()
                    # parts example: ['6:', '1709.03', 'cm**-1']
                    # sometimes: '0:', '0.00', 'cm**-1'
                    try:
                        val_str = parts[1]
                        val = float(val_str)
                        if val != 0.0: # Only non-zero? FCHK has all. User wants vibrations.
                             # But let's verify if user wants 0 modes (translation/rotation). Usually not.
                             # But keeping index alignment with Normal Modes is crucial.
                             # ORCA Normal Modes output block usually includes all 3*N modes.
                             pass
                        self.frequencies.append(val)
                    except: pass
                elif "NORMAL MODES" in l or "-----" in l:
                    if len(self.frequencies) > 0: # If we parsed some, stop
                         break
                curr += 1

        # Parse IR Spectrum for Intensities
        # Store as a dictionary: mode_id -> intensity
        intensity_map = {}  # {mode_idx: intensity}
        
        ir_start = -1
        for i, line in enumerate(lines):
             if "IR SPECTRUM" in line:
                 ir_start = i  # Keep updating to get the LAST occurrence
                 # Don't break - continue to find last one
                 
        if ir_start > 0:
            curr = ir_start + 6 # Skip headers and dashed line
            # Expected format:
            #  Mode   freq       eps      Int      T**2         TX        TY        TZ
            #        cm**-1   L/(mol*cm) km/mol    a.u.
            # ----------------------------------------------------------------------------
            #   6:   1709.03   0.015725   79.47  0.002871  (-0.001018 -0.053574 -0.000000)
            
            while curr < len(lines):
                l = lines[curr].strip()
                if not l:
                    curr += 1
                    continue
                if "-----" in l or "The first frequency" in l or "*" in l:
                    break
                    
                # line: "   6:   1709.03   0.015725   79.47  0.002871 ..."
                if ":" in l:
                    parts = l.split()
                    # parts[0] -> "6:"
                    # parts[1] -> Freq
                    # parts[2] -> eps
                    # parts[3] -> Int (km/mol)
                    # parts[4] -> T**2
                    try:
                         mode_id_str = parts[0].rstrip(':')
                         mode_id = int(mode_id_str)
                         if len(parts) >= 4:
                             inten = float(parts[3])
                             intensity_map[mode_id] = inten
                    except: 
                        pass
                curr += 1

        # Parse Normal Modes
        # ------------
        # NORMAL MODES
        # ------------
        # ...
        #                   0          1          2          3          4          5
        #       0       0.000000   0.000000   ...
        # ...
        #                   6          7          8
        #       0       0.001345  -0.000914   0.069964
        
        if modes_start_line > 0 and len(self.atoms) > 0:
            # We need to reconstruct full vectors for each mode.
            # ORCA output is blocked by columns (modes).
            # Rows are atoms * 3 (X, Y, Z coordinates).
            
            n_atoms = len(self.atoms)
            n_coords = n_atoms * 3
            # Initialize empty modes
            # We don't know exactly how many modes total yet, but freq list gives a hint.
            # Let's verify number of frequencies to initialize
            
            # Since parsing columnar data is tricky line-by-line without knowing total columns,
            # we will create a dictionary mode_index -> [vector]
            
            mode_data = {} # {mode_idx: [val0, val1, ...]}
            
            curr = modes_start_line + 7 # skip headers (approx)
            # Find first line of numbers
            while curr < len(lines) and not lines[curr].strip():
                curr += 1
                
            # Now we are at column headers: "                  0          1          2  ..."
            while curr < len(lines):
                header = lines[curr].strip()
                if not header:
                    curr += 1
                    continue
                if "IR SPECTRUM" in header or "--------" in header:
                    break
                    
                # Parse column indices
                try:
                    cols = [int(c) for c in header.split()]
                except ValueError:
                    # Maybe not a header line, skip
                    curr += 1
                    continue
                
                # Check next line: "      0       0.000000   0.000000 ..."
                # This corresponds to coordinate index 0 (Atom 0, X)
                
                start_data = curr + 1
                # Read 3*N lines for coordinates
                for row_idx in range(n_coords):
                    if start_data + row_idx >= len(lines): break
                    row_line = lines[start_data + row_idx].strip()
                    row_parts = row_line.split()
                    
                    # First part is coordinate index (0, 1, 2...)
                    # Remaining parts are values for the columns
                    values = row_parts[1:]
                    
                    if len(values) != len(cols):
                        # Mismatched data?
                        continue
                        
                    for c_idx, val_str in enumerate(values):
                        mode_idx = cols[c_idx]
                        val = float(val_str)
                        if mode_idx not in mode_data:
                            mode_data[mode_idx] = []
                        mode_data[mode_idx].append(val)
                
                curr = start_data + n_coords
                
            # Convert mode_data to self.vib_modes aligned with self.frequencies?
            # Frequencies list has indices 0...N.
            # We only care about modes that match frequencies we found.
            # Note ORCA freq list usually filters out 0.00 translations if "VIBRATIONAL FREQUENCIES" is used, 
            # BUT the list index in freq block ("   6: ...") matches the mode index.
            
            # Let's align them.
            # `self.frequencies` is just a flat list of values found. We need pairs (index, freq).
            # The freq parsing above was simplistic. Let's re-parse freq to get IDs.
            
            # Re-scanning frequencies for ID mapping:
            freq_map = {} # id -> freq
            if freq_start_line > 0:
                curr = freq_start_line + 4
                while curr < len(lines):
                    l = lines[curr].strip()
                    if ":" in l and "cm**-1" in l:
                        parts = l.split(':')
                        try:
                            mid = int(parts[0].strip())
                            freq_val = float(parts[1].split()[0])
                            freq_map[mid] = freq_val
                        except: pass
                    if "NORMAL MODES" in l: break
                    curr += 1
            
            # Now build final list
            sorted_mids = sorted(mode_data.keys())
            
            # Only keep modes that have frequencies (or all?)
            # Usually we only want non-zero modes (vibrations).
            # ORCA often prints 0-5 as translations/rotations.
            
            # Let's store pairs: (freq, vector)
            # Filter: only if freq > 10.0 cm-1 (avoid translations)
            
            self.final_modes = [] # item: {'freq': f, 'intensity': I, 'vector': [(x,y,z),...]}
            
            for mid in sorted_mids:
                freq = freq_map.get(mid, 0.0)
                # Use abs() to preserve imaginary frequencies (negative values)
                # Only exclude low-frequency modes (translations/rotations)
                if abs(freq) < 10.0: continue
                
                raw_vec = mode_data[mid]
                if len(raw_vec) != n_coords: continue
                
                # format to list of (dx,dy,dz)
                vec_formatted = []
                for k in range(0, len(raw_vec), 3):
                    vec_formatted.append((raw_vec[k], raw_vec[k+1], raw_vec[k+2]))
                
                # Get intensity for this mode from intensity_map
                intensity = intensity_map.get(mid, 0.0)
                    
                self.final_modes.append({'freq': freq, 'intensity': intensity, 'vector': vec_formatted})


class OrcaOutFreqAnalyzer(QWidget):
    def __init__(self, main_window, dock_widget=None):
        super().__init__(main_window)
        self.mw = main_window
        self.dock = dock_widget # Store reference
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
        self.lbl_info = QLabel("Drop .out file here or click Open")
        self.lbl_info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbl_info.setStyleSheet("border: 2px dashed #AAA; padding: 20px; color: #555;")
        layout.addWidget(self.lbl_info)
        
        # Open Button
        btn_open = QPushButton("Open ORCA Output")
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
                if fname.endswith(".out") or fname.endswith(".log"):
                    event.acceptProposedAction()
                    return
        event.ignore()
        
    def dropEvent(self, event):
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            if file_path.lower().endswith((".out", ".log")):
                self.load_file(file_path)
                event.acceptProposedAction()
                break
                
    def open_file_dialog(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Open ORCA Out", "", "Output Files (*.out *.log)")
        if fname:
            self.load_file(fname)
            
    def load_file(self, filename):
        # Validation before parsing
        if not is_valid_orca_file(filename):
            QMessageBox.critical(self, "Invalid File", "The selected file does not appear to be a valid ORCA output file.\n(Missing 'ORCA' header)")
            return

        self.parser = OrcaParser()
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
            QMessageBox.critical(self, "Error", f"Failed to parse Output:\n{e}")
            traceback.print_exc()

    def update_ui_after_load(self):
        self.list_freq.clear()
        if hasattr(self.parser, 'final_modes'):
            for i, mode in enumerate(self.parser.final_modes):
                freq = mode['freq']
                item = QTreeWidgetItem()
                item.setText(0, str(i + 1))  # Mode number
                item.setText(1, f"{freq:.2f}")
                
                # Get intensity from the mode dictionary
                inten = mode.get('intensity', 0.0)
                if inten > 0:
                     item.setText(2, f"{inten:.2f}")
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
        
        for sym in self.parser.atoms:
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
            self.stop_play()
        else:
            self.update_vectors()

    def toggle_play(self):
        curr = self.list_freq.currentItem()
        if not curr:
            return
            
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
        self.mw.draw_molecule_3d(self.base_mol)
        
        self.update_vectors()
        
        if hasattr(self.mw, 'plotter'):
            self.mw.plotter.render()

    def animate_frame(self):
        if not self.parser or not self.base_mol:
            self.stop_play()
            return
            
        curr = self.list_freq.currentItem()
        if not curr:
             self.stop_play()
             return
        row = self.list_freq.indexOfTopLevelItem(curr)
        if row < 0 or row >= len(self.parser.final_modes):
            return
            
        mode_data = self.parser.final_modes[row]
        mode_vecs = mode_data['vector']
        
        self.animation_step += 1
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
        if row < 0 or row >= len(self.parser.final_modes): return
        
        # Get vectors if not provided
        if mode_vecs is None:
            mode_data = self.parser.final_modes[row]
            mode_vecs = mode_data['vector']
            
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
        
        was_playing = self.is_playing
        if self.is_playing:
            self.toggle_play()

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
            return

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
        images = []
        mode_data = self.parser.final_modes[row]
        mode_vecs = mode_data['vector']
        
        self.reset_geometry()
        
        try:
            for i in range(20):
                cycle_pos = i / 20.0
                phase = cycle_pos * 2 * np.pi
                scale = self.slider_amp.value() / 20.0
                factor = np.sin(phase) * scale
                
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
                            processed_images.append(img.convert("RGBA"))
                        else:
                            processed_images.append(img.convert("RGB"))

                processed_images[0].save(file_path, save_all=True, append_images=processed_images[1:], duration=duration_ms, loop=0, disposal=2)
                QMessageBox.information(self, "Success", f"Saved GIF to:\n{file_path}")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save GIF: {e}")
            traceback.print_exc()
        finally:
            self.reset_geometry()
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
        # Use final_modes which has proper freq-intensity alignment
        # and already filters out low frequencies (translations/rotations)
        
        if not hasattr(self.parser, 'final_modes') or not self.parser.final_modes:
             QMessageBox.warning(self, "No Data", "No vibrational frequencies available.")
             return
        
        freqs = []
        intensities = []
        
        for mode in self.parser.final_modes:
            freqs.append(mode['freq'])
            intensities.append(mode.get('intensity', 0.0))
        
        dlg = SpectrumDialog(freqs, intensities, self)
        dlg.exec()


class SpectrumDialog(QDialog):
    def __init__(self, freqs, intensities, parent=None):
        super().__init__(parent)
        self.setWindowTitle("IR Spectrum")
        self.resize(800, 600)
        
        self.freqs = np.array(freqs)
        self.intensities = np.array(intensities)
        
        # Layout
        layout = QVBoxLayout(self)
        
        # Plot Area
        self.plot_widget = SpectrumPlotWidget(self.freqs, self.intensities)
        layout.addWidget(self.plot_widget)
        
        # Controls
        controls = QHBoxLayout()
        
        controls.addWidget(QLabel("Gaussian Broadening (FWHM, cm⁻¹):"))
        self.spin_fwhm = QSpinBox()
        self.spin_fwhm.setRange(1, 500)
        self.spin_fwhm.setValue(50)
        self.spin_fwhm.valueChanged.connect(self.on_fwhm_changed)
        controls.addWidget(self.spin_fwhm)
        
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

    def set_range(self, mn, mx):
        self.min_x = float(mn)
        self.max_x = float(mx)
        self.recalc_curve()
        self.update()

    def get_curve_data(self):
        return self.curve_x, self.curve_y

    def recalc_curve(self):
        if len(self.freqs) == 0: return
        
        # X resolution based on custom range
        self.curve_x = np.linspace(self.min_x, self.max_x, 1000)
        self.curve_y = np.zeros_like(self.curve_x)
        
        # Sum Gaussians
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
        # Inverted Y: 0 (top) -> Max (bottom) 
        
        def to_screen(x, y):
            # X: margin_l corresponds to max_x, w - margin_r corresponds to min_x
            # formula: x_ratio = (max_x - x) / range_x
            # sx = margin_l + x_ratio * plot_w
            sx = margin_l + (max_x - x) / range_x * plot_w
            
            # Y: margin_t corresponds to 0, h - margin_b corresponds to max_y
            # formula: y_ratio = y / max_y
            # sy = margin_t + y_ratio * plot_h
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
    dock = None
    analyzer = None
    
    # Check existing docks
    for d in main_window.findChildren(QDockWidget):
        if d.windowTitle() == "ORCA Output Freq Analyzer":
            dock = d
            analyzer = d.widget()
            break
            
    if not dock:
        dock = QDockWidget("ORCA Output Freq Analyzer", main_window)
        dock.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        analyzer = OrcaOutFreqAnalyzer(main_window, dock)
        dock.setWidget(analyzer)
        main_window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        
        # Connect visibility change
        dock.visibilityChanged.connect(analyzer.on_dock_visibility_changed)
    
    dock.show()
    dock.raise_()
    
    if analyzer:
        analyzer.load_file(fname)

def is_valid_orca_file(filepath):
    try:
        with open(filepath, 'r', errors='ignore') as f:
            # Check first 500 lines for "ORCA" keyword to be safe
            for _ in range(500):
                line = f.readline()
                if not line: break
                if "ORCA" in line or "O   R   C   A" in line:
                    return True
        return False
    except:
        return False

def run(main_window):
    # Smart Open Logic
    if hasattr(main_window, 'current_file_path') and main_window.current_file_path:
        fpath = main_window.current_file_path.lower()
        if fpath.endswith((".out", ".log")):
             # Validate content before auto-loading
             if is_valid_orca_file(main_window.current_file_path):
                 load_from_file(main_window, main_window.current_file_path)
                 return
                 
    fname, _ = QFileDialog.getOpenFileName(main_window, "Open ORCA Output", "", "ORCA Output (*.out *.log);;All Files (*)")
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
            
        main_window._custom_drop_handlers['.out'] = load_from_file
        main_window._custom_drop_handlers['.log'] = load_from_file
        # Register CLI handler too
        if hasattr(main_window, '_custom_loader_handlers'):
             main_window._custom_loader_handlers['.out'] = load_from_file
             main_window._custom_loader_handlers['.log'] = load_from_file
        
    except Exception as e:
        print(f"Failed to register ORCA DnD: {e}")

