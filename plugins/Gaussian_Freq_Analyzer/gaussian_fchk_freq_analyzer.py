import os
import re
import math
import numpy as np

from PySide6.QtGui import QPainter, QPen, QColor
import logging

try:
    from PIL import Image

    HAS_PIL = True
except ImportError:
    HAS_PIL = False


from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QSlider,
    QCheckBox,
    QFileDialog,
    QMessageBox,
    QDockWidget,
    QWidget,
    QApplication,
    QTreeWidget,
    QTreeWidgetItem,
    QHeaderView,
    QFormLayout,
    QDialogButtonBox,
    QSpinBox,
    QDoubleSpinBox,
)
from PySide6.QtCore import Qt, QTimer, QPointF

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from rdkit.Chem import rdDetermineBonds
except ImportError:
    Chem = None
    rdDetermineBonds = None

PLUGIN_NAME = "Gaussian Freq Analyzer"
PLUGIN_VERSION = "2026.10.19"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=5.0.0, <6.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Analyzes vibrational frequencies from Gaussian FCHK files. View IR spectrum, animate normal modes, and export GIF animations."
PLUGIN_CONTEXT = None


class FCHKParser:
    def __init__(self):
        self.filename = ""
        self.atoms = []  # List of atomic numbers
        self.coords = []  # List of (x, y, z) in Angstrom
        self.frequencies = []  # List of frequencies in cm^-1
        self.vib_modes = []  # List of displacement vectors (natoms, 3)
        self.charge = 0
        self.multiplicity = 1

    def parse(self, filename):
        self.charge = 0
        self.multiplicity = 1
        self.masses = []
        self.dipole_derivs = []
        self.filename = filename
        self.atoms = []
        self.coords = []
        self.frequencies = []
        self.intensities = []  # IR Intensities
        self.vib_modes = []  # list of list of (dx, dy, dz)

        with open(filename, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        data = {}
        wanted = {"Atomic numbers", "Current cartesian coordinates", "Vib-E2", "Vib-Modes", "Dipole Derivatives", "IR Inten", "IR intensities", "Real atomic weights", "Charge", "Multiplicity", "Number of Normal Modes"}
        wanted_labels = {name.casefold() for name in wanted}
        header_pattern = re.compile(r"^(.+?)\s+([IRCLH])\s+(?:N=\s*(\d+)|(.+))$")
        numeric_pattern = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?")
        i = 0
        while i < len(lines):
            match = header_pattern.match(lines[i].strip())
            i += 1
            if not match:
                continue
            label, dtype, count_text, scalar = match.groups()
            if label.casefold() not in wanted_labels:
                continue
            expected_type = "I" if label in {"Atomic numbers", "Charge", "Multiplicity", "Number of Normal Modes"} else "R"
            if dtype != expected_type:
                raise ValueError(f"Invalid FCHK record type: {label}")
            count = int(count_text) if count_text is not None else 1
            if count < 0:
                raise ValueError(f"Invalid FCHK count for {label}")
            values = []
            text = scalar if scalar is not None else ""
            while True:
                if text.strip():
                    # Like the MO Analyzer, accept packed Fortran numbers and D exponents.
                    tokens = numeric_pattern.findall(text)
                    if numeric_pattern.sub("", text).strip():
                        raise ValueError(f"Malformed FCHK numeric record: {label}")
                    for token in tokens:
                        value = int(token) if dtype == "I" else float(token.replace("D", "E").replace("d", "e"))
                        if not math.isfinite(value):
                            raise ValueError(f"Nonfinite FCHK value: {label}")
                        values.append(value)
                if len(values) >= count:
                    break
                if i >= len(lines) or header_pattern.match(lines[i].strip()):
                    raise ValueError(f"Incomplete FCHK record: {label}")
                text = lines[i].strip()
                i += 1
            if len(values) != count:
                raise ValueError(f"FCHK count mismatch: {label}")
            next_line = i
            while next_line < len(lines) and not lines[next_line].strip():
                next_line += 1
            if next_line < len(lines):
                extra = lines[next_line].strip()
                if numeric_pattern.findall(extra) and not numeric_pattern.sub("", extra).strip():
                    raise ValueError(f"Excess FCHK values after declared count: {label}")
            data[label] = values

        # Extract specific data
        if "Atomic numbers" in data:
            self.atoms = data["Atomic numbers"]

        BOHR_TO_ANG = 0.529177210903

        if "Current cartesian coordinates" in data:
            raw_coords = data["Current cartesian coordinates"]
            if len(raw_coords) != 3 * len(self.atoms):
                raise ValueError("FCHK coordinates do not match the atom count")
            # Convert to list of tuples (x,y,z)
            # FCHK coords are X1, Y1, Z1, X2... in Bohr
            coords_ang = []
            for j in range(0, len(raw_coords), 3):
                if j + 2 < len(raw_coords):
                    x = raw_coords[j] * BOHR_TO_ANG
                    y = raw_coords[j + 1] * BOHR_TO_ANG
                    z = raw_coords[j + 2] * BOHR_TO_ANG
                    coords_ang.append((x, y, z))
            self.coords = coords_ang

        if any(z < 1 or z > 118 for z in self.atoms):
            raise ValueError("Invalid FCHK atomic number")
        if self.atoms and len(self.coords) != len(self.atoms):
            raise ValueError("FCHK is missing atom coordinates")
        if "Vib-Modes" in data and (not self.atoms or len(data["Vib-Modes"]) % (3 * len(self.atoms))):
            raise ValueError("FCHK normal modes do not match the atom count")

        # Parse Vib-E2 section properly
        # Vib-E2 contains blocks: [Frequencies, Reduced Masses, Force Constants, IR Intensities, ...]
        # Each block has n_modes values
        if "Vib-E2" in data:
            raw_e2 = data["Vib-E2"]

            # Determine n_modes from Vib-Modes section (safest approach)
            if "Number of Normal Modes" in data:
                n_modes = data["Number of Normal Modes"][0]
                if n_modes < 0:
                    raise ValueError("Invalid number of FCHK normal modes")
                if "Vib-Modes" in data and len(data["Vib-Modes"]) != 3 * len(self.atoms) * n_modes:
                    raise ValueError("FCHK normal-mode count mismatch")
            elif "Vib-Modes" in data:
                n_atoms = len(self.atoms) if self.atoms else 0
                if n_atoms > 0:
                    # Vib-Modes size = 3 * n_atoms * n_modes
                    n_modes = len(data["Vib-Modes"]) // (3 * n_atoms)
                else:
                    n_modes = 0
            else:
                # Fallback: estimate from 3N-6 (or 3N-5 for linear)
                n_atoms = len(self.atoms)
                if n_atoms < 2:
                    n_modes = 0
                else:
                    origin = self.coords[0]
                    vectors = [tuple(point[j] - origin[j] for j in range(3)) for point in self.coords[1:]]
                    axis = max(vectors, key=lambda v: sum(x * x for x in v))
                    axis_norm = sum(x * x for x in axis)
                    linear = True
                    for vector in vectors:
                        cross = (axis[1] * vector[2] - axis[2] * vector[1],
                                 axis[2] * vector[0] - axis[0] * vector[2],
                                 axis[0] * vector[1] - axis[1] * vector[0])
                        cross_norm = sum(value * value for value in cross)
                        vector_norm = sum(value * value for value in vector)
                        if cross_norm > 1e-12 * axis_norm * vector_norm:
                            linear = False
                            break
                    n_modes = 3 * n_atoms - (5 if linear else 6)

            if n_modes > 0 and len(raw_e2) % n_modes:
                raise ValueError("FCHK frequency property blocks do not match the mode count")
            if n_modes > 0 and len(raw_e2) < n_modes:
                raise ValueError("Incomplete FCHK frequency block")
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
                pt = Chem.GetPeriodicTable()
                self.masses = [
                    pt.GetAtomicWeight(int(z)) for z in data["Atomic numbers"]
                ]

        if self.masses and (len(self.masses) != len(self.atoms) or any(mass <= 0 for mass in self.masses)):
            raise ValueError("FCHK atomic masses do not match the atoms")

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
                        dy = mode_vec[k + 1]
                        dz = mode_vec[k + 2]
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


class GaussianFCHKFreqAnalyzer(QWidget):
    def __init__(self, context, dock_widget=None):
        super().__init__(context.get_main_window())
        self.context = context
        # [DIRECT ACCESS] to main window for legacy plotter/view access
        self.mw = context.get_main_window()
        self.dock = dock_widget
        self.setAcceptDrops(True)

        self.parser = None
        self.base_mol = None
        self._bound_mol = None
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.animate_frame)
        self._structure_watch_timer = QTimer(self)
        self._structure_watch_timer.setInterval(250)
        self._structure_watch_timer.timeout.connect(self._check_structure_still_bound)
        self.animation_step = 0
        self.is_playing = False
        self.vector_actor = None

        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # Info Label
        self.lbl_info = QLabel("Drop .fchk file here or click Open")
        self.lbl_info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbl_info.setStyleSheet(
            "border: 2px dashed #AAA; padding: 20px; color: #555;"
        )
        layout.addWidget(self.lbl_info)

        # Open Button
        btn_open = QPushButton("Open FCHK File")
        btn_open.clicked.connect(self.open_file_dialog)
        layout.addWidget(btn_open)

        # Scaling Factor
        scale_layout = QHBoxLayout()
        scale_layout.addWidget(QLabel("Scaling Factor:"))
        self.spin_sf = QDoubleSpinBox()
        self.spin_sf.setRange(0.0, 2.0)
        self.spin_sf.setSingleStep(0.001)
        self.spin_sf.setDecimals(4)
        self.spin_sf.setValue(1.0)  # Default
        self.spin_sf.valueChanged.connect(self.update_list_and_spectrum_values)
        scale_layout.addWidget(self.spin_sf)
        scale_layout.addStretch()
        layout.addLayout(scale_layout)

        # Frequency List

        # Frequency List
        # Frequency List
        layout.addWidget(QLabel("Vibrational Frequencies:"))
        self.list_freq = QTreeWidget()
        self.list_freq.setColumnCount(3)
        self.list_freq.setHeaderLabels(
            ["No.", "Frequency (cm⁻¹)", "Intensity (km/mol)"]
        )
        self.list_freq.header().setSectionResizeMode(
            0, QHeaderView.ResizeMode.ResizeToContents
        )
        self.list_freq.header().setSectionResizeMode(
            1, QHeaderView.ResizeMode.ResizeToContents
        )
        self.list_freq.headerItem().setTextAlignment(
            0, Qt.AlignmentFlag.AlignCenter
        )  # Center No. header
        self.list_freq.headerItem().setTextAlignment(
            1, Qt.AlignmentFlag.AlignCenter
        )  # Center Frequency header
        self.list_freq.headerItem().setTextAlignment(
            2, Qt.AlignmentFlag.AlignCenter
        )  # Center Intensity header
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
        self.slider_speed.valueChanged.connect(
            lambda v: self.lbl_speed_val.setText(str(v))
        )

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
            # Restore 2D editing UI when leaving the 3D viewer mode
            try:
                mw = (
                    self.mw
                    if getattr(self, "mw", None) is not None
                    else self.context.get_main_window()
                )
                if hasattr(mw, "ui_manager") and hasattr(
                    mw.ui_manager, "restore_ui_for_editing"
                ):
                    mw.ui_manager.restore_ui_for_editing()
            except Exception as _e:
                logging.warning("[gaussian_fchk_freq_analyzer.py:428] silenced: %s", _e)
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
                if (
                    fname.endswith(".fchk")
                    or fname.endswith(".fck")
                    or fname.endswith(".fch")
                ):
                    event.acceptProposedAction()
                    return
        event.ignore()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            if file_path.lower().endswith((".fchk", ".fck", ".fch")):
                self.load_file(file_path)
                event.acceptProposedAction()
                break

    def open_file_dialog(self):
        fname, _ = QFileDialog.getOpenFileName(
            self, "Open Gaussian FCHK", "", "FCHK Files (*.fchk *.fck *.fch)"
        )
        if fname:
            self.load_file(fname)

    def load_file(self, filename):
        candidate = FCHKParser()
        try:
            candidate.parse(filename)
            self.stop_play()
            self.remove_vectors()
            self.parser = candidate
            self.animation_step = 0
            self.update_ui_after_load()
            self.lbl_info.setText(os.path.basename(filename))
            self.lbl_info.setStyleSheet(
                "border: 2px solid #4CAF50; padding: 10px; color: #4CAF50;"
            )

            # Reset SF if needed? Or keep user setting? User usually keeps preferred SF.
            # self.spin_sf.setValue(1.0)

            # Update Main Window Context
            if hasattr(self.mw, "init_manager") and hasattr(
                self.mw.init_manager, "current_file_path"
            ):
                self.mw.init_manager.current_file_path = filename
                self.context.refresh_ui()

        except Exception as e:
            logging.exception("Failed to parse FCHK: %s", e)
            QMessageBox.critical(self, "Error", f"Failed to parse FCHK:\n{e}")

    def update_ui_after_load(self):
        self.list_freq.clear()
        # Ensure freqs match modes logic handled in parser consistency check
        if hasattr(self.parser, "frequencies"):
            for i, freq in enumerate(self.parser.frequencies):
                item = QTreeWidgetItem()
                item.setText(0, str(i + 1))  # Mode number
                item.setText(1, f"{freq:.2f}")

                if hasattr(self.parser, "intensities") and i < len(
                    self.parser.intensities
                ):
                    inten = self.parser.intensities[i]
                    # Higher precision to see small values and compare with LOG
                    item.setText(2, f"{inten:.4f}")
                else:
                    item.setText(2, "-")
                item.setTextAlignment(
                    0, Qt.AlignmentFlag.AlignCenter
                )  # Center mode number
                item.setTextAlignment(
                    1, Qt.AlignmentFlag.AlignCenter
                )  # Center frequency
                item.setTextAlignment(
                    2, Qt.AlignmentFlag.AlignCenter
                )  # Center intensity
                self.list_freq.addTopLevelItem(item)

        self.lbl_meta.setText(
            f"Charge: {self.parser.charge}, Multiplicity: {self.parser.multiplicity}\nAtoms: {len(self.parser.atoms)}"
        )

        # Load molecule into main window
        if len(self.parser.atoms) > 0 and Chem:
            self.create_base_molecule()

    def update_list_and_spectrum_values(self):
        if not self.parser or not self.parser.frequencies:
            return

        sf = self.spin_sf.value()

        # Update List items (assume order matches parser.frequencies)
        # We can iterate top level items
        root = self.list_freq.invisibleRootItem()
        for i in range(root.childCount()):
            item = root.child(i)
            # Index is in col 0 text "1", "2"... but assume linear
            idx = int(item.text(0)) - 1
            if 0 <= idx < len(self.parser.frequencies):
                raw_freq = self.parser.frequencies[idx]
                scaled_freq = raw_freq * sf
                item.setText(1, f"{scaled_freq:.2f}")

    def create_base_molecule(self):
        if not self.parser:
            return

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

        try:
            rdDetermineBonds.DetermineConnectivity(mol)
            rdDetermineBonds.DetermineBondOrders(mol, charge=int(self.parser.charge))
        except (ImportError, RuntimeError, AttributeError, TypeError, ValueError):
            if hasattr(self.mw, "io_manager") and hasattr(
                self.mw.io_manager, "estimate_bonds_from_distances"
            ):
                self.mw.io_manager.estimate_bonds_from_distances(mol)

        self.base_mol = mol.GetMol()
        self.context.enter_3d_mode()
        self.context.current_molecule = self.base_mol
        self._bound_mol = self.base_mol
        self._structure_watch_timer.start()
        self.context.reset_3d_camera()

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
        if row < 0:
            return

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
        if fps <= 0:
            fps = 1
        interval = 1000 / fps
        self.timer.setInterval(int(interval))

    def _check_structure_still_bound(self):
        if self._bound_mol is None:
            return True
        if self.context.current_molecule is self._bound_mol:
            return True
        self.timer.stop()
        self._structure_watch_timer.stop()
        self.is_playing = False
        self.btn_play.setText("Play")
        self.remove_vectors()
        self.base_mol = None
        self._bound_mol = None
        return False

    def closeEvent(self, event):
        self.close_plugin()
        self._structure_watch_timer.stop()
        super().closeEvent(event)

    def reset_geometry(self):
        if not self._check_structure_still_bound():
            return
        if not self.base_mol or not self.parser:
            return
        conf = self.base_mol.GetConformer()
        for idx, (x, y, z) in enumerate(self.parser.coords):
            conf.SetAtomPosition(idx, Point3D(x, y, z))
        self.context.draw_molecule_3d(self.base_mol)

        self.update_vectors()

        self.context.plotter.render()

    def animate_frame(self):
        if not self._check_structure_still_bound():
            return
        if not self.parser or not self.base_mol:
            self.stop_play()
            return

        curr = self.list_freq.currentItem()
        if not curr:
            return
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
        self.context.draw_molecule_3d(self.base_mol)
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
        if self.vector_actor:
            try:
                self.context.plotter.remove_actor(self.vector_actor)
            except (RuntimeError, AttributeError, KeyError, TypeError, ValueError) as _e:
                logging.warning("[gaussian_fchk_freq_analyzer.py:659] silenced: %s", _e)
        self.vector_actor = None

    def update_vectors(self, mode_vecs=None, scale_factor=0.0):
        if not self._check_structure_still_bound():
            return
        # Clean up existing vectors
        self.remove_vectors()

        if not self.chk_vectors.isChecked():
            return

        if not self.parser or not self.base_mol:
            return

        # Get current frequency
        curr = self.list_freq.currentItem()
        if not curr:
            return
        row = self.list_freq.indexOfTopLevelItem(curr)
        if row < 0 or row >= len(self.parser.vib_modes):
            return

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

        if not coords:
            return

        coords = np.array(coords)
        vectors = np.array(vectors)

        try:
            self.vector_actor = self.context.plotter.add_arrows(
                coords,
                vectors,
                mag=vis_scale,
                color="lightgreen",
                show_scalar_bar=False,
            )
        except (RuntimeError, AttributeError, KeyError, TypeError, ValueError) as e:
            logging.warning("Error adding arrows: %s", e)

    def save_as_gif(self):
        if not self.parser or not self.base_mol:
            return

        # Pause to configure
        was_playing = self.is_playing
        if self.is_playing:
            self.toggle_play()  # Pause

        curr = self.list_freq.currentItem()
        if not curr:
            QMessageBox.warning(
                self, "Select Frequency", "Please select a frequency to export."
            )
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

        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        form.addRow(btns)

        if dialog.exec() != QDialog.DialogCode.Accepted:
            if was_playing:
                self.toggle_play()
            return  # Cancel

        target_fps = spin_fps.value()
        use_transparent = chk_transparent.isChecked()
        use_hq = chk_hq.isChecked()

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save GIF", "", "GIF Files (*.gif)"
        )
        if not file_path:
            if was_playing:
                self.toggle_play()
            return

        if not file_path.lower().endswith(".gif"):
            file_path += ".gif"

        # Generate Frames
        # 1 Cycle = 20 steps
        images = []
        mode_vecs = self.parser.vib_modes[row]

        # Store current geometry to restore later
        self.reset_geometry()  # align to base

        try:
            for i in range(20):
                cycle_pos = i / 20.0
                phase = cycle_pos * 2 * np.pi
                scale = self.slider_amp.value() / 20.0
                factor = np.sin(phase) * scale  # Calculate factor here
                self.apply_displacement(mode_vecs, factor)
                self.context.draw_molecule_3d(self.base_mol)
                self.update_vectors(mode_vecs, factor)
                self.context.plotter.render()

                img_array = self.context.plotter.screenshot(
                    transparent_background=use_transparent, return_img=True
                )
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
                            img_p = img_rgb.convert(
                                "P", palette=Image.Palette.ADAPTIVE, colors=255
                            )
                            # Set simple transparency
                            mask = Image.eval(alpha, lambda a: 255 if a <= 128 else 0)
                            img_p.paste(255, mask)
                            img_p.info["transparency"] = 255
                            processed_images.append(img_p)
                        else:
                            processed_images.append(
                                img.convert(
                                    "P", palette=Image.Palette.ADAPTIVE, colors=256
                                )
                            )
                    else:
                        if use_transparent:
                            img = img.convert("RGBA")
                            processed_images.append(img)
                        else:
                            processed_images.append(img.convert("RGB"))

                processed_images[0].save(
                    file_path,
                    save_all=True,
                    append_images=processed_images[1:],
                    duration=duration_ms,
                    loop=0,
                    disposal=2,
                )
                QMessageBox.information(self, "Success", f"Saved GIF to:\n{file_path}")

        except Exception as e:
            logging.exception("Failed to save GIF: %s", e)
            QMessageBox.critical(self, "Error", f"Failed to save GIF: {e}")
        finally:
            self.reset_geometry()
            # Restore play state if needed, or leave paused
            # User might want to inspect
            if was_playing:
                self.toggle_play()

    def on_dock_visibility_changed(self, visible):
        if not visible:
            self.close_plugin()
        elif self._bound_mol is not None:
            self._structure_watch_timer.start()
            self.update_vectors()

    def close_plugin(self):
        self._structure_watch_timer.stop()
        self.stop_play()
        self.remove_vectors()
        self.animation_step = 0

    def show_spectrum(self):
        if not self.parser or not self.parser.frequencies:
            return

        # Filter out low frequencies (translations/rotations)
        freqs = []
        intensities = []

        parser_intensities = self.parser.intensities

        sf = self.spin_sf.value()

        for i, freq in enumerate(self.parser.frequencies):
            # Use abs() to preserve imaginary frequencies (negative values)
            # Only exclude low-frequency modes (translations/rotations)
            if freq != 0 and i < len(parser_intensities):
                freqs.append(freq * sf)  # Apply Scale

                # Check Intensity
                val = 0.0
                if i < len(parser_intensities):
                    val = parser_intensities[i]
                intensities.append(val)

        if not freqs:
            QMessageBox.information(self, "Info", "No frequencies with IR intensities are available.")
            return

        # Show Plot
        dlg = SpectrumDialog(
            freqs,
            intensities,
            title=f"IR Spectrum (SF: {sf:.4f}) - {os.path.basename(self.parser.filename)}",
            parent=self,
        )
        dlg.exec()


class SpectrumDialog(QDialog):
    def __init__(self, freqs, intensities, title="IR Spectrum", parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
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

        controls.addWidget(QLabel("FWHM (cm⁻¹):"))
        self.spin_fwhm = QSpinBox()
        self.spin_fwhm.setRange(1, 500)
        self.spin_fwhm.setValue(50)
        self.spin_fwhm.valueChanged.connect(self.on_fwhm_changed)
        controls.addWidget(self.spin_fwhm)

        # Axis Range
        controls.addWidget(QLabel("Min WN:"))
        self.spin_min = QSpinBox()
        self.spin_min.setRange(0, 10000)
        self.spin_min.setValue(0)
        self.spin_min.setSingleStep(100)
        self.spin_min.valueChanged.connect(self.on_range_changed)
        controls.addWidget(self.spin_min)

        controls.addWidget(QLabel("Max WN:"))
        self.spin_max = QSpinBox()
        self.spin_max.setRange(0, 10000)
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

        # Initial Plot.  The spin boxes were populated before their signals
        # were connected, so nothing has pushed their values into the plot
        # widget yet — without this the first curve (and any CSV exported from
        # it) is broadened with the widget's own default rather than the FWHM
        # the dialog is displaying.
        self.plot_widget.set_fwhm(float(self.spin_fwhm.value()))
        self.on_range_changed()

    def on_fwhm_changed(self, val):
        self.plot_widget.set_fwhm(val)

    def on_range_changed(self):
        mn = self.spin_min.value()
        mx = self.spin_max.value()
        if mx > mn:
            self.plot_widget.set_range(mn, mx)

    def export_csv(self):
        fname, _ = QFileDialog.getSaveFileName(
            self, "Save Spectrum Data", "", "CSV Files (*.csv)"
        )
        if fname:
            if not fname.lower().endswith(".csv"):
                fname += ".csv"
            try:
                x, y = self.plot_widget.get_curve_data()
                with open(fname, "w", encoding="utf-8") as f:
                    f.write("Frequency,Intensity\n")
                    for xi, yi in zip(x, y):
                        f.write(f"{xi:.2f},{yi:.4f}\n")
                QMessageBox.information(self, "Success", "Saved CSV.")
            except (OSError, RuntimeError, AttributeError, TypeError, ValueError) as e:
                QMessageBox.critical(self, "Error", str(e))

    def export_png(self):
        fname, _ = QFileDialog.getSaveFileName(
            self, "Save Spectrum Image", "", "PNG Files (*.png)"
        )
        if fname:
            if not fname.lower().endswith(".png"):
                fname += ".png"
            try:
                # Capture the widget
                pixmap = self.plot_widget.grab()
                pixmap.save(fname)
                QMessageBox.information(self, "Success", "Saved Image.")
            except (RuntimeError, AttributeError) as e:
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
        if len(self.freqs) == 0:
            return

        # X resolution based on custom range
        self.curve_x = np.linspace(self.min_x, self.max_x, 1000)
        self.curve_y = np.zeros_like(self.curve_x)

        # Sum Gaussians weighted by intensity
        sigma = self.fwhm / 2.35482

        for f, i in zip(self.freqs, self.intensities):
            self.curve_y += i * np.exp(-((self.curve_x - f) ** 2) / (2 * sigma**2))

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
        if max_y == 0:
            max_y = 1.0
        min_x = np.min(self.curve_x)
        max_x = np.max(self.curve_x)
        range_x = max_x - min_x
        if range_x == 0:
            range_x = 100

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
        painter.drawLine(
            margin_l, margin_t, w - margin_r, margin_t
        )  # X-axis at top (Baseline)
        painter.drawLine(
            margin_l, h - margin_b, w - margin_r, h - margin_b
        )  # X-axis at bottom

        painter.drawLine(margin_l, h - margin_b, margin_l, margin_t)  # Y (Left)
        painter.drawLine(
            w - margin_r, h - margin_b, w - margin_r, margin_t
        )  # Y (Right border)

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
            painter.drawText(
                int(px) - 20,
                h - margin_b + 5,
                40,
                20,
                Qt.AlignmentFlag.AlignCenter,
                f"{int(val)}",
            )
            painter.drawLine(int(px), h - margin_b, int(px), h - margin_b + 5)

        # Labels
        font.setPointSize(14)  # Increased from 10
        font.setBold(True)
        painter.setFont(font)
        painter.drawText(
            0, h - 25, w, 20, Qt.AlignmentFlag.AlignCenter, "Wavenumber (cm⁻¹)"
        )

        # Draw baseline at y=0
        baseline_x_start, baseline_y = to_screen(max_x, 0)
        baseline_x_end, _ = to_screen(min_x, 0)
        painter.setPen(QPen(QColor(150, 150, 150), 1, Qt.PenStyle.DashLine))
        painter.drawLine(
            int(baseline_x_start), int(baseline_y), int(baseline_x_end), int(baseline_y)
        )

        # Draw Curve
        painter.setPen(QPen(Qt.GlobalColor.blue, 2))
        path_points = []
        for x, y in zip(self.curve_x, self.curve_y):
            sx, sy = to_screen(x, y)
            path_points.append((sx, sy))

        if len(path_points) > 1:
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
                self.context.plotter.render()


def load_from_file(context, fname):
    main_window = context.get_main_window()
    # FIRST: Close conflicting plugin (ORCA)
    for d in main_window.findChildren(QDockWidget):
        if d.windowTitle() == "ORCA Output Freq Analyzer":
            d.close()

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
        dock.setAllowedAreas(
            Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea
        )
        analyzer = GaussianFCHKFreqAnalyzer(context, dock)
        dock.setWidget(analyzer)

        # [DIRECT ACCESS] for dock injection
        main_window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)

        # Register for management
        context.register_window("main_panel", dock)

        # Connect visibility change
        if hasattr(analyzer, "on_dock_visibility_changed"):
            dock.visibilityChanged.connect(analyzer.on_dock_visibility_changed)

    dock.show()
    dock.raise_()

    if analyzer:
        analyzer.load_file(fname)


def run_plugin(context):
    """Entry point for manual launch via menu."""
    mw = context.get_main_window()
    # Smart Open Logic
    if (
        hasattr(mw, "init_manager")
        and hasattr(mw.init_manager, "current_file_path")
        and mw.init_manager.current_file_path
    ):
        fpath = mw.init_manager.current_file_path.lower()
        if fpath.endswith((".fchk", ".fck")):
            load_from_file(context, mw.init_manager.current_file_path)
            return

    fname, _ = QFileDialog.getOpenFileName(
        mw, "Open Gaussian FCHK", "", "Gaussian FCHK (*.fchk *.fck);;All Files (*)"
    )
    if fname:
        load_from_file(context, fname)


def initialize(context):
    """Initialize the Gaussian Freq Analyzer plugin."""
    global PLUGIN_CONTEXT
    PLUGIN_CONTEXT = context
    # 1. Register Menu Action
    # context.add_menu_action("File/Gaussian Freq Analyzer...", lambda: run_plugin(context))

    # 2. Register File Openers
    def load_wrapper(fname):
        load_from_file(context, fname)

    context.register_file_opener(".fchk", load_wrapper)
    context.register_file_opener(".fck", load_wrapper)
    context.register_file_opener(".fch", load_wrapper)

    # 3. Register Drop Handler
    def drop_handler(file_path):
        if file_path.lower().endswith((".fchk", ".fck", ".fch")):
            load_from_file(context, file_path)
            return True
        return False

    context.register_drop_handler(drop_handler, priority=10)


def run(mw):
    if not hasattr(mw, "plugin_manager"):
        return

    context = PLUGIN_CONTEXT
    if not context:
        return

    run_plugin(context)
