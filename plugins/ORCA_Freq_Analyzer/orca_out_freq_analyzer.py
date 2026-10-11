import os
import re
import math
import numpy as np
from PyQt6.QtWidgets import (
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
    QFormLayout,
    QDialogButtonBox,
    QSpinBox,
    QApplication,
    QTreeWidget,
    QTreeWidgetItem,
    QHeaderView,
    QDoubleSpinBox,
)
from PyQt6.QtGui import QPainter, QPen, QColor

try:
    from PIL import Image

    HAS_PIL = True
except ImportError:
    HAS_PIL = False
from PyQt6.QtCore import Qt, QTimer, QPointF
import logging

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
except ImportError:
    Chem = None

PLUGIN_NAME = "ORCA Freq Analyzer"
PLUGIN_DESCRIPTION = "Analyzes vibrational frequencies from ORCA output files. View IR spectrum, animate normal modes, and export GIF animations."
PLUGIN_VERSION = "2026.10.13"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_CONTEXT = None


def _find_intensity_column(lines, section_start, wanted, default):
    """Index of the first column named in *wanted* within a spectrum header.

    ORCA's header token positions line up with the data rows only once inline
    unit tokens are dropped: the ORCA 4 IR header reads
    "Mode freq (cm**-1) T**2 ..." while its rows carry no unit column, so
    counting "(cm**-1)" shifts every lookup one place right.

    Falls back to *default* when no header is recognized.
    """
    end = min(section_start + 6, len(lines))
    for i in range(section_start + 1, end):
        tokens = lines[i].split()
        if not tokens or tokens[0].lower() != "mode":
            continue
        lowered = [t.lower().strip(":") for t in tokens if not t.startswith("(")]
        for name in wanted:
            if name in lowered:
                return lowered.index(name)
    return default


class OrcaParser:
    def __init__(self):
        self.filename = ""
        self.atoms = []
        self.coords = []
        self.frequencies = []
        self.intensities = []  # IR Intensities
        self.vib_modes = []  # list of list of (dx, dy, dz)
        self.charge = 0
        self.multiplicity = 1
        self.final_modes = []

    def parse(self, filename):
        self.__init__()
        self.filename = filename
        with open(filename, "r", encoding="utf-8", errors="strict") as source:
            lines = source.readlines()

        def number(token):
            value = float(token.replace("D", "E").replace("d", "e"))
            if not math.isfinite(value):
                raise ValueError("ORCA frequency data contains nonfinite values")
            return value

        freq_headers = [i for i, line in enumerate(lines) if "VIBRATIONAL FREQUENCIES" in line.upper()]
        freq_start = freq_headers[-1] if freq_headers else -1
        coord_headers = [i for i, line in enumerate(lines) if "CARTESIAN COORDINATES (ANGSTROEM)" in line.upper() and (freq_start < 0 or i < freq_start)]
        for line in lines[:freq_start if freq_start >= 0 else len(lines)]:
            if "Total Charge" in line and "Charge" in line:
                self.charge = int(line.split()[-1])
            elif "Multiplicity" in line and "Mult" in line:
                self.multiplicity = int(line.split()[-1])
        if coord_headers:
            started = False
            for line in lines[coord_headers[-1] + 1:]:
                parts = line.split()
                if not parts or set(line.strip()) == {"-"}:
                    if started:
                        break
                    continue
                if len(parts) != 4:
                    break
                self.atoms.append(parts[0])
                self.coords.append(tuple(number(value) for value in parts[1:]))
                started = True
        if freq_start < 0:
            return

        # Keep tables attached to the last frequency calculation, as in the
        # ORCA Result Analyzer; earlier jobs cannot supply missing intensities.
        after = lines[freq_start + 1:]
        freq_map = {}
        for line in after:
            if "NORMAL MODES" in line.upper() or "IR SPECTRUM" in line.upper():
                break
            match = re.match(r"\s*(\d+)\s*:\s*(\S+)\s+cm(?:\*\*)?-1", line)
            if match:
                mode_id = int(match.group(1))
                if mode_id in freq_map:
                    raise ValueError("Duplicate ORCA frequency mode ID")
                freq_map[mode_id] = number(match.group(2))
        self.frequencies = list(freq_map.values())
        modes_headers = [i for i in range(freq_start + 1, len(lines)) if "NORMAL MODES" in lines[i].upper()]
        ir_headers = [i for i in range(freq_start + 1, len(lines)) if "IR SPECTRUM" in lines[i].upper()]
        intensity_map = {}
        if ir_headers:
            start = ir_headers[-1]
            column = _find_intensity_column(lines, start, ("int", "t**2"), 3)
            started = False
            for line in lines[start + 1:]:
                parts = line.split()
                if parts and re.fullmatch(r"\d+:", parts[0]):
                    if len(parts) <= column:
                        raise ValueError("Incomplete ORCA IR row")
                    intensity_map[int(parts[0][:-1])] = number(parts[column])
                    started = True
                elif started and (not parts or "-----" in line or "The first frequency" in line):
                    break
        if not modes_headers or not self.atoms:
            return

        n_coords = 3 * len(self.atoms)
        mode_data = {}
        curr = modes_headers[-1] + 1
        while curr < len(lines):
            text = lines[curr].strip()
            if "IR SPECTRUM" in text or "RAMAN SPECTRUM" in text:
                break
            parts = text.split()
            if not parts or not all(part.isdigit() for part in parts):
                curr += 1
                continue
            cols = [int(part) for part in parts]
            if len(set(cols)) != len(cols) or any(mode in mode_data for mode in cols):
                raise ValueError("Duplicate ORCA normal-mode column")
            vectors = {mode: [] for mode in cols}
            for row in range(n_coords):
                curr += 1
                if curr >= len(lines):
                    raise ValueError("Incomplete ORCA normal-mode block")
                fields = lines[curr].split()
                if len(fields) != len(cols) + 1 or fields[0] != str(row):
                    raise ValueError("ORCA normal-mode coordinate rows are incomplete or out of order")
                for mode, value in zip(cols, fields[1:]):
                    vectors[mode].append(number(value))
            mode_data.update(vectors)
            curr += 1

        first_vibration = None
        for line in after:
            match = re.search(r"The first frequency considered to be a vibration is\s+(\d+)", line)
            if match:
                first_vibration = int(match.group(1))
        if mode_data and any(mode not in mode_data and frequency != 0 for mode, frequency in freq_map.items()):
            raise ValueError("ORCA frequency is missing its normal-mode vector")
        for mode, vector in sorted(mode_data.items()):
            if mode not in freq_map:
                raise ValueError("ORCA normal mode has no matching frequency")
            frequency = freq_map[mode]
            if (first_vibration is not None and mode < first_vibration) or frequency == 0:
                continue
            self.final_modes.append(dict(freq=frequency, intensity=intensity_map.get(mode), vector=[tuple(vector[i:i + 3]) for i in range(0, n_coords, 3)]))


class OrcaOutFreqAnalyzer(QWidget):
    def __init__(self, main_window, dock_widget=None):
        super().__init__(main_window)
        self.mw = main_window
        self.dock = dock_widget  # Store reference
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
        self.lbl_info = QLabel("Drop .out file here or click Open")
        self.lbl_info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbl_info.setStyleSheet(
            "border: 2px dashed #AAA; padding: 20px; color: #555;"
        )
        layout.addWidget(self.lbl_info)

        # Open Button
        btn_open = QPushButton("Open ORCA Output")
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
        fname, _ = QFileDialog.getOpenFileName(
            self, "Open ORCA Out", "", "Output Files (*.out *.log)"
        )
        if fname:
            self.load_file(fname)

    def load_file(self, filename):
        # Validation before parsing
        if not is_valid_orca_file(filename):
            QMessageBox.critical(
                self,
                "Invalid File",
                "The selected file does not appear to be a valid ORCA output file.\n(Missing 'ORCA' header)",
            )
            return

        candidate = OrcaParser()
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

            # Reset SF if needed? Or keep user setting?
            # self.spin_sf.setValue(1.0)

        except (OSError, ValueError, RuntimeError, AttributeError) as e:
            logging.exception("Failed to parse Output: %s", e)
            QMessageBox.critical(self, "Error", f"Failed to parse Output:\n{e}")

    def update_ui_after_load(self):
        self.list_freq.clear()
        if hasattr(self.parser, "final_modes"):
            for i, mode in enumerate(self.parser.final_modes):
                freq = mode["freq"]
                item = QTreeWidgetItem()
                item.setText(0, str(i + 1))  # Mode number
                item.setText(1, f"{freq:.2f}")

                # Get intensity from the mode dictionary
                inten = mode.get("intensity", None)
                if inten is not None:
                    item.setText(2, f"{inten:.2f}")
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
        if not self.parser or not self.parser.final_modes:
            return

        sf = self.spin_sf.value()

        # Update List values
        # Top level iteration
        root = self.list_freq.invisibleRootItem()
        for i in range(root.childCount()):
            item = root.child(i)
            # Find corresponding mode freq
            idx = int(item.text(0)) - 1
            if 0 <= idx < len(self.parser.final_modes):
                mode = self.parser.final_modes[idx]
                scaled = mode["freq"] * sf
                item.setText(1, f"{scaled:.2f}")

    def create_base_molecule(self):
        if not self.parser:
            return

        mol = Chem.RWMol()

        for sym in self.parser.atoms:
            atom = Chem.Atom(sym)
            mol.AddAtom(atom)

        conf = Chem.Conformer(len(self.parser.atoms))
        for idx, (x, y, z) in enumerate(self.parser.coords):
            conf.SetAtomPosition(idx, Point3D(x, y, z))
        mol.AddConformer(conf)

        if hasattr(self.mw, "io_manager") and hasattr(
            self.mw.io_manager, "estimate_bonds_from_distances"
        ):
            self.mw.io_manager.estimate_bonds_from_distances(mol)

        self.base_mol = mol.GetMol()
        PLUGIN_CONTEXT.enter_3d_mode()
        PLUGIN_CONTEXT.current_mol = self.base_mol
        PLUGIN_CONTEXT.reset_3d_camera()

    def on_freq_selected(self, current, previous):
        if self.is_playing:
            # Transition smoothly: Reset geometry base for next frame calculation?
            # animate_frame uses base coords + displacement, so it handles switch naturally.
            # Just ensure vectors update if needed (animate_frame updates them too)
            pass
        else:
            self.reset_geometry()
            self.update_vectors()

    def toggle_play(self):
        curr = self.list_freq.currentItem()
        if not curr:
            return

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
        if PLUGIN_CONTEXT.current_molecule is self._bound_mol:
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
        PLUGIN_CONTEXT.draw_molecule_3d(self.base_mol)

        self.update_vectors()

        PLUGIN_CONTEXT.plotter.render()

    def animate_frame(self):
        if not self._check_structure_still_bound():
            return
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
        mode_vecs = mode_data["vector"]

        self.animation_step += 1
        cycle_pos = (self.animation_step % 20) / 20.0
        phase = cycle_pos * 2 * np.pi

        scale = self.slider_amp.value() / 20.0
        factor = np.sin(phase) * scale

        self.apply_displacement(mode_vecs, factor)
        PLUGIN_CONTEXT.draw_molecule_3d(self.base_mol)
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
                PLUGIN_CONTEXT.plotter.remove_actor(self.vector_actor)
            except (RuntimeError, AttributeError, KeyError, ValueError) as _e:
                logging.warning("[orca_out_freq_analyzer.py:847] silenced: %s", _e)
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
        if row < 0 or row >= len(self.parser.final_modes):
            return

        # Get vectors if not provided
        if mode_vecs is None:
            mode_data = self.parser.final_modes[row]
            mode_vecs = mode_data["vector"]

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
            self.vector_actor = PLUGIN_CONTEXT.plotter.add_arrows(
                coords,
                vectors,
                mag=vis_scale,
                color="lightgreen",
                show_scalar_bar=False,
            )
        except (RuntimeError, AttributeError, KeyError, ValueError) as e:
            logging.warning("Error adding arrows: %s", e)

    def save_as_gif(self):
        if not self.parser or not self.base_mol:
            return

        was_playing = self.is_playing
        if self.is_playing:
            self.toggle_play()

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
            return

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
        images = []
        mode_data = self.parser.final_modes[row]
        mode_vecs = mode_data["vector"]

        self.reset_geometry()

        try:
            for i in range(20):
                cycle_pos = i / 20.0
                phase = cycle_pos * 2 * np.pi
                scale = self.slider_amp.value() / 20.0
                factor = np.sin(phase) * scale

                self.apply_displacement(mode_vecs, factor)
                PLUGIN_CONTEXT.draw_molecule_3d(self.base_mol)
                self.update_vectors(mode_vecs, factor)
                PLUGIN_CONTEXT.plotter.render()

                img_array = PLUGIN_CONTEXT.plotter.screenshot(
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
                            processed_images.append(img.convert("RGBA"))
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
        # Use final_modes which has proper freq-intensity alignment
        # and already filters out low frequencies (translations/rotations)

        if not hasattr(self.parser, "final_modes") or not self.parser.final_modes:
            QMessageBox.warning(
                self, "No Data", "No vibrational frequencies available."
            )
            return

        freqs = []
        intensities = []

        sf = self.spin_sf.value()

        for mode in self.parser.final_modes:
            if mode.get("intensity") is None:
                continue
            freqs.append(mode["freq"] * sf)
            intensities.append(mode.get("intensity", 0.0))

        if not freqs:
            QMessageBox.warning(self, "No IR Data", "No IR intensities are available for these modes.")
            return

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
            except (OSError, RuntimeError, AttributeError, ValueError) as e:
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

    def set_range(self, mn, mx):
        self.min_x = float(mn)
        self.max_x = float(mx)
        self.recalc_curve()
        self.update()

    def get_curve_data(self):
        return self.curve_x, self.curve_y

    def recalc_curve(self):
        if len(self.freqs) == 0:
            return

        # X resolution based on custom range
        self.curve_x = np.linspace(self.min_x, self.max_x, 1000)
        self.curve_y = np.zeros_like(self.curve_x)

        # Sum Gaussians
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


def load_from_file(main_window, fname):
    pass

    # FIRST: Close conflicting plugin (Gaussian)
    for d in main_window.findChildren(QDockWidget):
        if d.windowTitle() == "Gaussian Freq Analyzer":
            d.close()

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
        dock.setAllowedAreas(
            Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea
        )
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
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            # Check first 500 lines for "ORCA" keyword to be safe
            for _ in range(500):
                line = f.readline()
                if not line:
                    break
                if "ORCA" in line or "O   R   C   A" in line:
                    return True
        return False
    except (OSError, ValueError):
        return False


def get_host_file_path(mw):
    """The file the host currently has open.

    It lives on ``init_manager``; ``mw.current_file_path`` is not a proxy, so
    reading it always came back empty.
    """
    getter = getattr(mw, "get_current_file_path", None)
    if callable(getter):
        return getter()
    im = getattr(mw, "init_manager", None)
    return getattr(im, "current_file_path", None) if im is not None else None


def run(mw):
    # Smart Open Logic
    current = get_host_file_path(mw)
    if current:
        fpath = current.lower()
        if fpath.endswith((".out", ".log")):
            if is_valid_orca_file(current):
                load_from_file(mw, current)
                return

    fname, _ = QFileDialog.getOpenFileName(
        mw, "Open ORCA Output", "", "ORCA Output (*.out *.log);;All Files (*)"
    )
    if not fname:
        return

    if is_valid_orca_file(fname):
        load_from_file(mw, fname)
    else:
        QMessageBox.warning(
            mw, "Invalid File", "This does not appear to be a valid ORCA output file."
        )


def initialize(context):
    global PLUGIN_CONTEXT
    PLUGIN_CONTEXT = context
    mw = context.get_main_window()

    def load_wrapper(fname):
        # Validate ORCA file before loading
        if is_valid_orca_file(fname):
            load_from_file(mw, fname)
        else:
            print(f"Skipping invalid ORCA file: {fname}")

    # 1. Register File Openers
    # Note: .log can be generic, but our valid check helps
    context.register_file_opener(".out", load_wrapper)
    context.register_file_opener(".log", load_wrapper)

    # 2. Register Drop Handler
    def drop_handler(file_path):
        fpath = file_path.lower()
        if fpath.endswith((".out", ".log")):
            if is_valid_orca_file(file_path):
                load_from_file(mw, file_path)
                return True
        return False

    if hasattr(context, "register_drop_handler"):
        context.register_drop_handler(drop_handler, priority=10)
