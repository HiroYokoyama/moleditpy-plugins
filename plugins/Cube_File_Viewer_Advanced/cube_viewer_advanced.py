
import os
import json
import tempfile
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QFileDialog, QDockWidget, QWidget, QVBoxLayout, 
                             QSlider, QLabel, QHBoxLayout, QPushButton, QMessageBox, 
                             QDoubleSpinBox, QColorDialog, QInputDialog, QDialog, 
                             QFormLayout, QDialogButtonBox, QSpinBox, QCheckBox, QComboBox, QLineEdit)
from PyQt6.QtGui import QColor
from PyQt6.QtCore import Qt, QTimer, QCoreApplication

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
    
__version__="2026.02.04"
__author__="HiroYokoyama"
PLUGIN_NAME = "Cube File Viewer Advanced"
PLUGIN_DESCRIPTION = "Advanced 3D visualization for Gaussian Cube files with PBR, SSAO, and other effects."

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
    
    # Skip metadata lines (e.g. "1 150") before data starts
    while current_line < len(lines):
        line_content = lines[current_line].strip()
        parts = line_content.split()
        
        # Skip empty lines
        if not parts:
            current_line += 1
            continue
            
        # Skip short lines (metadata often has few columns, data usually has 6)
        if len(parts) < 6:
            current_line += 1
            continue
            
        # Check if start is numeric
        try:
            float(parts[0])
        except ValueError:
            current_line += 1
            continue
            
        # If we get here, it's likely data
        break

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

class FlexibleDoubleSpinBox(QDoubleSpinBox):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDecimals(10) # Set extremely high to allow typing, but we format display manually

    def textFromValue(self, value):
        # Return simplest string representation (stripping 0s)
        # Using string formatting to avoid scientific notation if possible
        text = f"{value:.10f}".rstrip('0').rstrip('.')
        if text == "": text = "0"
        return text

class CubeViewerWidget(QWidget):
    def __init__(self, parent_window, dock_widget, grid, data_max=1.0):
        super().__init__(parent_window)
        self.mw = parent_window
        self.dock = dock_widget 
        self.grid = grid
        # Ensure we have a reasonable positive max value
        self.data_max = max(abs(float(data_max)), 1e-6)
        
        self.iso_actor_p = None
        self.iso_actor_n = None
        self.iso_actor_p_sil = None
        self.iso_actor_n_sil = None
        
        self.plotter = self.mw.plotter
        
        # Initial Colors
        self.color_p = (0, 0, 255)  # Default blue
        self.color_n = (255, 0, 0)  # Default orange
        # Alternative: Blue/Red scheme: (0, 100, 255) / (255, 100, 0)

        # Advanced Rendering State
        self.use_pbr = False
        self.metallic = 0.5
        self.roughness = 0.5
        self.use_ssao = False
        self.use_depth_peeling = False
        self.use_silhouette = False
        self.env_texture_path = ""
        
        # New Graphics Options
        self.use_aa = False
        self.use_edl = False
        self.edl_strength = 0.2
        self.use_shadows = False
        self.light_intensity = 1.0
        

        # Atom PBR specific - REMOVED (Handled by different plugin)

        
        # Presets {name: settings_dict}
        self.presets = {}
        self.default_preset_names = set()  # Track which presets are defaults (read-only)
        self._init_default_presets()

        self.init_ui()
        
        # Populate preset dropdown with defaults
        self.update_preset_combo()
        
        # DELAY INITIAL UPDATE: 
        # When opening from command line, the window might not be fully ready.
        # A small delay ensures the plotter is ready for actors.
        QTimer.singleShot(100, self.initial_update)
        
        # Ensure settings are saved when application quits (Main Window X button)
        QCoreApplication.instance().aboutToQuit.connect(self.save_settings)

    def initial_update(self):
        self.update_iso()
        self.plotter.reset_camera()
        self.plotter.render()

    def init_ui(self):
        layout = QVBoxLayout()
        
        # --- Isovalue Controls ---
        ctrl_layout = QHBoxLayout()

        preset_layout = QHBoxLayout()
        preset_layout.addWidget(QLabel("Preset:"))
        
        self.combo_presets = QComboBox()
        self.combo_presets.setEditable(False)  # Not editable - use Save dialog to name presets
        self.combo_presets.setPlaceholderText("Select preset...")
        self.combo_presets.setMinimumWidth(150)
        self.combo_presets.currentTextChanged.connect(self.on_preset_text_changed)
        self.combo_presets.activated.connect(self.on_preset_activated)
        preset_layout.addWidget(self.combo_presets)
        
        btn_save_preset = QPushButton("Save")
        btn_save_preset.setToolTip("Save current settings")
        btn_save_preset.clicked.connect(self.save_preset)
        preset_layout.addWidget(btn_save_preset)
        
        btn_del_preset = QPushButton("Del")
        btn_del_preset.setToolTip("Delete preset")
        btn_del_preset.clicked.connect(self.delete_preset)
        preset_layout.addWidget(btn_del_preset)
        
        layout.addLayout(preset_layout)
        
        ctrl_layout.addWidget(QLabel("Isovalue:"))
        
        # Spinbox: Fixed max 1.0 as requested
        # Spinbox: Dynamic max based on data
        # We allow going a bit higher than the max found in data, e.g. 1.2x, just in case
        self.max_val = self.data_max * 1.2
        
        self.spin = FlexibleDoubleSpinBox()
        self.spin.setRange(0.00001, self.max_val) 
        self.spin.setSingleStep(0.01) # User requested 0.01 step
        # self.spin.setDecimals(10) # already set in class __init__
        
        # Default value strategy: 
        # User requested hardcoded 0.05
        default_val = 0.05
        
        # If default is outside the max range, adjust it
        if default_val > self.max_val:
            default_val = self.max_val * 0.1
            
        self.spin.setValue(default_val)
        
        self.spin.valueChanged.connect(self.on_spin_changed)
        ctrl_layout.addWidget(self.spin)
        
        # Slider
        self.slider_max_int = 1000
        
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(0, self.slider_max_int)
        
        # Default 0.05
        # Default set from spinbox value
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
        
        # --- Style Controls ---
        style_layout = QHBoxLayout()
        style_layout.addWidget(QLabel("Style:"))
        self.combo_style = QComboBox()
        self.combo_style.addItems(["Surface", "Smoothed Surface", "Wireframe", "Points"])
        self.combo_style.currentIndexChanged.connect(self.update_iso)
        style_layout.addWidget(self.combo_style)
        
        self.check_smooth = QCheckBox("Smooth Shading")
        self.check_smooth.setChecked(True)
        self.check_smooth.toggled.connect(self.update_iso)
        style_layout.addWidget(self.check_smooth)
        
        layout.addLayout(style_layout)

        # --- Advanced Rendering ---
        self.init_advanced_ui(layout)

        # --- Reset Button (Outside Groups) ---
        btn_reset = QPushButton("Reset All Settings")
        btn_reset.setToolTip("Reset all viewer settings to defaults")
        btn_reset.clicked.connect(self.on_reset_all)
        layout.addWidget(btn_reset)

        layout.addStretch()
        
        close_btn = QPushButton("Close Plugin")
        close_btn.clicked.connect(self.close_plugin)
        layout.addWidget(close_btn)
        self.setLayout(layout)
        
        # Load settings after UI init
        self.load_settings()

    def init_advanced_ui(self, layout):
        from PyQt6.QtWidgets import QGroupBox, QGridLayout, QTabWidget, QWidget, QVBoxLayout, QComboBox, QMessageBox, QInputDialog

        group = QGroupBox("Advanced Rendering")
        # Main Layout of Group is just a VBox holding the TabWidget
        group_layout = QVBoxLayout()
        
        tabs = QTabWidget()
        
        # --- TAB 1: SCENE (Global Effects) ---
        tab_scene = QWidget()
        scene_layout = QGridLayout()
        
        # Texture
        scene_layout.addWidget(QLabel("Env Texture:"), 0, 0)
        tex_layout = QHBoxLayout()
        self.line_env_path = QLineEdit()
        self.line_env_path.setPlaceholderText("Path to .png/.jpg texture...")
        self.line_env_path.returnPressed.connect(self.on_texture_path_entered)
        tex_layout.addWidget(self.line_env_path)
        btn_browse_tex = QPushButton("...")
        btn_browse_tex.setFixedWidth(30); btn_browse_tex.clicked.connect(self.on_load_env_texture)
        tex_layout.addWidget(btn_browse_tex)
        btn_clear_tex = QPushButton("X")
        btn_clear_tex.setFixedWidth(30); btn_clear_tex.clicked.connect(self.on_remove_env_texture)
        tex_layout.addWidget(btn_clear_tex)
        scene_layout.addLayout(tex_layout, 0, 1)

        # Shadows
        self.check_shadows = QCheckBox("Shadows")
        self.check_shadows.setToolTip("Enable Shadows (requires compatible lighting)")
        self.check_shadows.toggled.connect(self.on_shadows_toggled)
        scene_layout.addWidget(self.check_shadows, 1, 0)

        # Light Intensity
        scene_layout.addWidget(QLabel("Light Int:"), 2, 0)
        self.slider_light = QSlider(Qt.Orientation.Horizontal)
        self.slider_light.setRange(0, 500)
        self.slider_light.setValue(int(self.light_intensity * 100))
        self.slider_light.valueChanged.connect(self.on_light_intensity_changed)
        scene_layout.addWidget(self.slider_light, 2, 1)

        # SSAO
        self.check_ssao = QCheckBox("SSAO")
        if not hasattr(self.plotter, 'enable_ssao'): self.check_ssao.setEnabled(False)
        self.check_ssao.toggled.connect(self.on_ssao_toggled)
        scene_layout.addWidget(self.check_ssao, 3, 0)

        # AA
        self.check_aa = QCheckBox("Anti-Aliasing")
        self.check_aa.toggled.connect(self.on_aa_toggled)
        scene_layout.addWidget(self.check_aa, 3, 1)
        
        # EDL
        self.check_edl = QCheckBox("EDL (Depth)")
        self.check_edl.toggled.connect(self.on_edl_toggled)
        scene_layout.addWidget(self.check_edl, 4, 0)
        
        edl_layout = QHBoxLayout()
        edl_layout.addWidget(QLabel("Str:"))
        self.slider_edl = QSlider(Qt.Orientation.Horizontal)
        self.slider_edl.setRange(1, 100) 
        self.slider_edl.setValue(int(self.edl_strength * 100)) 
        self.slider_edl.valueChanged.connect(self.on_edl_strength_changed)
        self.slider_edl.setEnabled(False)
        edl_layout.addWidget(self.slider_edl)
        scene_layout.addLayout(edl_layout, 4, 1)

        tab_scene.setLayout(scene_layout)
        tabs.addTab(tab_scene, "Scene")

        # --- TAB 2: ORBITAL (Isosurface) ---
        tab_orb = QWidget()
        orb_layout = QGridLayout()
        
        # PBR Toggle
        self.check_pbr = QCheckBox("Enable PBR")
        self.check_pbr.toggled.connect(self.on_pbr_toggled)
        orb_layout.addWidget(self.check_pbr, 0, 0, 1, 2)
        
        # Metallic
        orb_layout.addWidget(QLabel("Metallic:"), 1, 0)
        self.slider_metallic = QSlider(Qt.Orientation.Horizontal)
        self.slider_metallic.setRange(0, 100)
        self.slider_metallic.setValue(int(self.metallic * 100))
        self.slider_metallic.valueChanged.connect(self.on_metallic_changed)
        self.slider_metallic.setEnabled(False)
        orb_layout.addWidget(self.slider_metallic, 1, 1)

        # Roughness
        orb_layout.addWidget(QLabel("Roughness:"), 2, 0)
        self.slider_roughness = QSlider(Qt.Orientation.Horizontal)
        self.slider_roughness.setRange(0, 100)
        self.slider_roughness.setValue(int(self.roughness * 100))
        self.slider_roughness.valueChanged.connect(self.on_roughness_changed)
        self.slider_roughness.setEnabled(False)
        orb_layout.addWidget(self.slider_roughness, 2, 1)
        
        # Silhouette
        self.check_silhouette = QCheckBox("Silhouette")
        self.check_silhouette.toggled.connect(self.on_silhouette_toggled)
        orb_layout.addWidget(self.check_silhouette, 3, 0)
        
        # Depth Peeling
        self.check_depth = QCheckBox("Depth Peeling")
        if not hasattr(self.plotter, 'enable_depth_peeling'): self.check_depth.setEnabled(False)
        self.check_depth.toggled.connect(self.on_depth_peeling_toggled)
        orb_layout.addWidget(self.check_depth, 3, 1)
        
        tab_orb.setLayout(orb_layout)
        tabs.addTab(tab_orb, "Orbital")


        group_layout.addWidget(tabs)
        group.setLayout(group_layout)
        layout.addWidget(group)


    def _init_default_presets(self):
        """Initialize built-in default presets that cannot be overwritten or deleted."""
        # Default preset - standard settings
        self.presets["Default"] = {
            "isovalue": 0.02,
            "color_p": [0, 0, 255],
            "color_n": [255, 0, 0],
            "opacity": 0.4,
            "metallic": 0.5,
            "roughness": 0.5,
            "use_pbr": False,
            "env_texture_path": "",
            "use_ssao": False,
            "use_depth_peeling": False,
            "use_silhouette": False,
            "use_aa": False,
            "use_edl": False,
            "edl_strength": 0.2,
            "use_shadows": False,
            "edl_strength": 0.2,
            "use_shadows": False,
            "light_intensity": 1.0,
            "smooth_shading": True,
            "style": "Surface"
        }
        
        
        # Mark these as default (read-only)
        self.default_preset_names = {"Default"}


    def get_settings_path(self):
        """Returns the path to the JSON settings file in the same directory."""
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "cube_viewer_advanced.json")

    def load_settings(self):
        """Loads settings from JSON file."""
        try:
            settings_path = self.get_settings_path()
            if os.path.exists(settings_path):
                with open(settings_path, 'r') as f:
                    settings = json.load(f)
                
                # Apply settings
                # Merge user presets with defaults (don't overwrite defaults)
                if "presets" in settings:
                    user_presets = settings["presets"]
                    # Add user presets that aren't defaults
                    for name, preset_data in user_presets.items():
                        if name not in self.default_preset_names:
                            self.presets[name] = preset_data
                    self.update_preset_combo()

                if "isovalue" in settings:
                    self.spin.setValue(float(settings["isovalue"]))
                
                if "color_p" in settings:
                    col = settings["color_p"]
                    # JSON deserializes tuples as lists, convert back
                    if isinstance(col, list):
                        col = tuple(col)
                    self.color_p = col
                    self.btn_color_p.setStyleSheet(f"background-color: rgb{self.color_p}; border: 1px solid gray;")

                if "color_n" in settings:
                    col = settings["color_n"]
                    # JSON deserializes tuples as lists, convert back
                    if isinstance(col, list):
                        col = tuple(col)
                    self.color_n = col
                    self.btn_color_n.setStyleSheet(f"background-color: rgb{self.color_n}; border: 1px solid gray;")
                
                if "use_complementary" in settings:
                    self.check_comp_color.setChecked(bool(settings["use_complementary"]))
                    # Force update of neg color button enabled state
                    self.btn_color_n.setEnabled(not self.check_comp_color.isChecked())

                if "opacity" in settings:
                    val = float(settings["opacity"])
                    self.opacity_spin.setValue(val)
                    
                if "style" in settings:
                    # Robust find
                    text = settings["style"]
                    idx = self.combo_style.findText(text)
                    if idx < 0:
                        # try case insensitive
                        for i in range(self.combo_style.count()):
                            if self.combo_style.itemText(i).lower() == text.lower():
                                idx = i
                                break
                    if idx >= 0: 
                        self.combo_style.setCurrentIndex(idx)
                    
                if "smooth_shading" in settings:
                    self.check_smooth.setChecked(bool(settings["smooth_shading"]))

                # Advanced Settings
                if "env_texture_path" in settings:
                    path = settings["env_texture_path"]
                    # Check for spaces first (VTK can't handle them)
                    if path and ' ' in path:
                        print(f"Skipping texture with spaces in path: {path}")
                        self.env_texture_path = ""
                    elif path and os.path.exists(path):
                        try:
                            self.load_env_texture(path)
                        except Exception as e:
                            print(f"Failed to restore texture from settings: {e}")
                            # Clear the invalid path
                            self.env_texture_path = ""

                if "metallic" in settings: 
                    self.metallic = float(settings["metallic"])
                    self.slider_metallic.setValue(int(self.metallic * 100))
                if "roughness" in settings: 
                    self.roughness = float(settings["roughness"])
                    self.slider_roughness.setValue(int(self.roughness * 100))
                
                # Load PBR toggle last so it checks against loaded metallic/texture state
                if "use_pbr" in settings: self.check_pbr.setChecked(bool(settings["use_pbr"]))

                if "use_ssao" in settings: self.check_ssao.setChecked(bool(settings["use_ssao"]))
                if "use_depth_peeling" in settings: self.check_depth.setChecked(bool(settings["use_depth_peeling"]))
                if "use_silhouette" in settings: self.check_silhouette.setChecked(bool(settings["use_silhouette"]))
                
                # New Graphics Settings
                if "use_aa" in settings: self.check_aa.setChecked(bool(settings["use_aa"]))
                
                if "edl_strength" in settings:
                    self.edl_strength = float(settings["edl_strength"])
                    self.slider_edl.setValue(int(self.edl_strength * 100))
                    
                if "use_edl" in settings: self.check_edl.setChecked(bool(settings["use_edl"]))
                
                if "use_edl" in settings: self.check_edl.setChecked(bool(settings["use_edl"]))



                # Update Text Field if texture loaded
                if self.env_texture_path:
                    self.line_env_path.setText(self.env_texture_path)
                
                # Redraw orbital with loaded settings
                self.update_iso()
                self.plotter.render()
                    
        except Exception as e:
            print(f"Error loading settings: {e}")

    def save_settings(self):
        """Saves current settings to JSON file."""
        try:
            settings = {
                "isovalue": self.spin.value(),
                "color_p": self.color_p,
                "color_n": self.color_n,
                "use_complementary": self.check_comp_color.isChecked(),
                "opacity": self.opacity_spin.value(),
                "style": self.combo_style.currentText(),
                "smooth_shading": self.check_smooth.isChecked(),
                "use_pbr": self.check_pbr.isChecked(),
                "metallic": self.metallic,
                "roughness": self.roughness,
                "use_ssao": self.check_ssao.isChecked(),
                "use_depth_peeling": self.check_depth.isChecked(),
                "use_silhouette": self.check_silhouette.isChecked(),
                "env_texture_path": self.env_texture_path,
                "use_aa": self.use_aa,
                "use_edl": self.use_edl,
                "edl_strength": self.edl_strength,
                "use_shadows": self.use_shadows,
                "light_intensity": self.light_intensity,
                "use_shadows": self.use_shadows,
                "light_intensity": self.light_intensity,
                "presets": self.presets
            }
            
            # Helper to sanitize data for JSON
            def sanitize(obj):
                if isinstance(obj, np.generic):
                    return obj.item()
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                if isinstance(obj, (QColor,)):  # Add other PyQt types if needed
                    return [obj.red(), obj.green(), obj.blue()]
                if isinstance(obj, dict):
                    return {k: sanitize(v) for k, v in obj.items()}
                if isinstance(obj, list):
                    return [sanitize(v) for v in obj]
                if isinstance(obj, tuple):
                    return [sanitize(v) for v in obj] # JSON doesn't support tuples
                return obj
            
            clean_settings = sanitize(settings)
            
            with open(self.get_settings_path(), 'w') as f:
                json.dump(clean_settings, f, indent=4)
                
        except Exception as e:
            print(f"Error saving settings: {e}")

    def on_comp_color_toggled(self, checked):
        self.btn_color_n.setEnabled(not checked)
        if checked:
            self.update_complementary_color()

    def update_complementary_color(self):
        # Ensure color is a tuple (might be string from old settings)
        if isinstance(self.color_p, str):
            col_temp = QColor(self.color_p)
            self.color_p = (col_temp.red(), col_temp.green(), col_temp.blue())
        
        # Unpack RGB tuple for QColor
        col_p = QColor(*self.color_p)
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
        self.color_n = (col_n.red(), col_n.green(), col_n.blue())
        self.btn_color_n.setStyleSheet(f"background-color: rgb{self.color_n}; border: 1px solid gray;")
        self.update_iso()

    def choose_color_p(self):
        # Ensure color is a tuple (might be string from old settings)
        if isinstance(self.color_p, str):
            col_temp = QColor(self.color_p)
            self.color_p = (col_temp.red(), col_temp.green(), col_temp.blue())
        
        # Unpack RGB tuple for QColor
        c = QColorDialog.getColor(initial=QColor(*self.color_p), title="Select Positive Lobe Color")
        if c.isValid():
            self.color_p = (c.red(), c.green(), c.blue())
            self.btn_color_p.setStyleSheet(f"background-color: rgb{self.color_p}; border: 1px solid gray;")
            
            if self.check_comp_color.isChecked():
                self.update_complementary_color()
            else:
                self.update_iso()

    def choose_color_n(self):
        # Ensure color is a tuple (might be string from old settings)
        if isinstance(self.color_n, str):
            col_temp = QColor(self.color_n)
            self.color_n = (col_temp.red(), col_temp.green(), col_temp.blue())
        
        # Unpack RGB tuple for QColor
        c = QColorDialog.getColor(initial=QColor(*self.color_n), title="Select Negative Lobe Color")
        if c.isValid():
            self.color_n = (c.red(), c.green(), c.blue())
            self.btn_color_n.setStyleSheet(f"background-color: rgb{self.color_n}; border: 1px solid gray;")
            self.update_iso()

    def update_iso(self):
        val = self.spin.value()
        # Prefer numeric spinbox value (kept in sync with slider)
        opacity_val = self.opacity_spin.value() if hasattr(self, 'opacity_spin') else self.opacity_slider.value() / 100.0
        
        # Style settings
        style = self.combo_style.currentText() # Keep original case for check
        smooth = self.check_smooth.isChecked()
        
        render_style = style.lower()
        do_geometric_smooth = False
        is_density_mode = False
        
        if style == "Smoothed Surface":
            render_style = "surface"
            do_geometric_smooth = True
        elif render_style in ["wireframe", "points"]:
            is_density_mode = True
            
        # Update Label to reflect usage
        if is_density_mode:
            self.opacity_label.setText("Density:")
        else:
            self.opacity_label.setText("Opacity:")
        
        try:
            # Cleanup previous
            if self.iso_actor_p: self.plotter.remove_actor(self.iso_actor_p)
            if self.iso_actor_n: self.plotter.remove_actor(self.iso_actor_n)
            
            # Explicit cleanup of silhouette actors
            if self.iso_actor_p_sil: self.plotter.remove_actor(self.iso_actor_p_sil)
            if self.iso_actor_n_sil: self.plotter.remove_actor(self.iso_actor_n_sil)
            self.iso_actor_p_sil = None
            self.iso_actor_n_sil = None
            
            # Additional safety cleanup by name
            self.plotter.remove_actor("cube_iso_p")
            self.plotter.remove_actor("cube_iso_n")
            
            # Using full grid
            using_grid = self.grid

            # Decimation / Opacity Logic
            render_opacity = opacity_val
            target_reduction = 0.0
            
            if is_density_mode:
                # In density mode, opacity slider controls density (1.0 = full, 0.0 = empty)
                # target_reduction is fraction to REMOVE.
                # If opacity_val is 1.0 (Density), reduction is 0.0
                # If opacity_val is 0.1 (Density), reduction is 0.9
                target_reduction = 1.0 - opacity_val
                # Clamp to avoid total destruction or errors
                target_reduction = max(0.0, min(0.99, target_reduction))
                
                # Force full opacity for the lines/points themselves so they are visible
                render_opacity = 1.0

            def process_mesh(iso_mesh):
                m = iso_mesh
                
                # 1. Geometric Smoothing (Before decimation usually better for shape, 
                # but if decimating heavily, smoothing after might be cleaner? 
                # Let's stick to user request order: Smoothing is a style feature, usually implies high qual.
                # Decimation is for sparse view.)
                if do_geometric_smooth:
                     m = m.smooth(n_iter=100)
                     
                # 2. Decimation (Density)
                if is_density_mode and target_reduction > 0.0:
                    m = m.decimate(target_reduction)
                    
                return m

            # Positive lobe
            iso_p = using_grid.contour(isosurfaces=[val])
            if iso_p.n_points > 0:
                iso_p = process_mesh(iso_p)
                self.iso_actor_p = self.plotter.add_mesh(
                    iso_p, 
                    color=self.color_p, 
                    opacity=render_opacity, 
                    name="cube_iso_p", 
                    reset_camera=False, 
                    style=render_style, 
                    smooth_shading=smooth,
                    pbr=self.use_pbr,
                    metallic=self.metallic,
                    roughness=self.roughness,
                    # silhouette=self.use_silhouette # Doing manual handle
                )
                if self.use_silhouette and hasattr(self.plotter, 'add_silhouette'):
                     self.iso_actor_p_sil = self.plotter.add_silhouette(iso_p, color='black', line_width=2.0)
                
            # Negative lobe
            iso_n = using_grid.contour(isosurfaces=[-val])
            if iso_n.n_points > 0:
                iso_n = process_mesh(iso_n)
                self.iso_actor_n = self.plotter.add_mesh(
                    iso_n, 
                    color=self.color_n, 
                    opacity=render_opacity, 
                    name="cube_iso_n", 
                    reset_camera=False, 
                    style=render_style, 
                    smooth_shading=smooth,
                    pbr=self.use_pbr,
                    metallic=self.metallic,
                    roughness=self.roughness,
                    # silhouette=self.use_silhouette # Doing manual handle
                )
                if self.use_silhouette and hasattr(self.plotter, 'add_silhouette'):
                     self.iso_actor_n_sil = self.plotter.add_silhouette(iso_n, color='black', line_width=2.0)
                
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

    # --- Advanced Rendering Slots ---

    def on_pbr_toggled(self, checked):
        self.use_pbr = checked
        self.slider_metallic.setEnabled(checked)
        self.slider_roughness.setEnabled(checked)
        
        # Removed MessageBox as requested by user
        if checked and not self.env_texture_path:
            # Just print to console or set status tip if needed, but avoid popup
            print("PBR enabled without environment texture. Surfaces may appear dark.")
        
        self.update_iso()


    def on_metallic_changed(self, val):
        self.metallic = val / 100.0
        self.update_iso()


    def on_roughness_changed(self, val):
        self.roughness = val / 100.0
        self.update_iso()


    def on_ssao_toggled(self, checked):
        self.use_ssao = checked
        self._disable_conflicting_effects(exclude="ssao" if checked else "")

        try:
            self._clean_render_pipeline()
            if checked:
                if hasattr(self.plotter, 'enable_ssao'):

                    self.plotter.enable_ssao()
            else:
                if hasattr(self.plotter, 'disable_ssao'):
                    self.plotter.disable_ssao()
            self.plotter.render()
        except Exception as e:
            print(f"SSAO error: {e}")

    def _clean_render_pipeline(self):
        """
        重要: パス切り替え時に発生するVTKエラーを防ぐため、
        一度レンダリングパイプラインを完全にリセットする。
        """
        if not self.plotter or not hasattr(self.plotter, 'renderer'): return
        
        try:
            # 既存のパスを強制解除 (これがReleaseGraphicsResourcesエラーを防ぐ鍵)
            if hasattr(self.plotter.renderer, 'SetPasses'):
                self.plotter.renderer.SetPasses(None)
            
            # Note: Do NOT call render() here. Rendering with None passes can cause
            # VTK to complain about missing resources or invalid state.
            # self.plotter.render() 
        except Exception as e:
            print(f"Pipeline clean error: {e}")


    def on_depth_peeling_toggled(self, checked):
        self.use_depth_peeling = checked
        self._disable_conflicting_effects(exclude="depth" if checked else "")
            
        try:
            self._clean_render_pipeline() # Clean before switching transparency mode
            if checked:
                if hasattr(self.plotter, 'enable_depth_peeling'):
                    self.plotter.enable_depth_peeling()
            else:
                if hasattr(self.plotter, 'disable_depth_peeling'):
                    self.plotter.disable_depth_peeling()
            self.plotter.render()
        except Exception as e:
            print(f"Depth Peeling error: {e}")

    def on_silhouette_toggled(self, checked):
        self.use_silhouette = checked
        self.update_iso()
             
    def on_shadows_toggled(self, checked):
        self.use_shadows = checked
        self._disable_conflicting_effects(exclude="shadows" if checked else "")

        try:
            self._clean_render_pipeline()
            if checked:
                self.plotter.enable_shadows()
            else:
                self.plotter.disable_shadows()
        except Exception:
            pass
        self.plotter.render()
            
    def on_light_intensity_changed(self, val):
        self.light_intensity = val / 100.0
        try:
            # Update all lights
            if hasattr(self.plotter, 'renderer'):
                lights = self.plotter.renderer.GetLights()
                for light in lights:
                    light.SetIntensity(self.light_intensity)
            self.plotter.render()
        except Exception:
            pass  # Silently fail for compatibility

    def on_aa_toggled(self, checked):
        self.use_aa = checked
        try:
            if checked:
                self.plotter.enable_anti_aliasing()
            else:
                self.plotter.disable_anti_aliasing()
            self.plotter.render()
        except Exception as e:
            print(f"AA Error: {e}")

    def on_edl_toggled(self, checked):
        self.use_edl = checked
        self._disable_conflicting_effects(exclude="edl" if checked else "")
            
        self.slider_edl.setEnabled(checked)
        self._clean_render_pipeline()
        self.apply_edl()

    def on_edl_strength_changed(self, val):
        self.edl_strength = val / 100.0
        if self.use_edl:
            self.apply_edl()
            
    def apply_edl(self):
        try:
            if self.use_edl:
                self.plotter.enable_eye_dome_lighting()
                
                # Access the EDL pass to set strength
                try:
                    if hasattr(self.plotter.renderer, '_edl_pass'):
                        edl_pass = self.plotter.renderer._edl_pass
                        if edl_pass and hasattr(edl_pass, 'SetEDLStrength'):
                            edl_pass.SetEDLStrength(self.edl_strength)
                except Exception:
                    pass
                
            else:
                if hasattr(self.plotter, 'disable_eye_dome_lighting'):
                    self.plotter.disable_eye_dome_lighting()
            self.plotter.render()
        except Exception as e:
            # Plotter might be closed or invalid
            pass 

    def _disable_conflicting_effects(self, exclude=""):
        """
        排他制御ヘルパー:
        - Depth Peeling (透明モード)
        - EDL / Shadows / SSAO (不透明エフェクト)
        これらは共存できないため、一方をONにしたら他方を無効化(Disable & Uncheck)する。
        """
        # 無限ループ防止
        self.blockSignals(True) 
        
        # 1. Depth PeelingがONの場合 -> 他のエフェクトをOFFにしてDisable
        if exclude == "depth" and self.use_depth_peeling:
            # EDL check OFF
            if self.use_edl:
                self.check_edl.setChecked(False)
                self.use_edl = False
                try: 
                    if hasattr(self.plotter, 'disable_eye_dome_lighting'):
                        self.plotter.disable_eye_dome_lighting()
                except: pass
            
            # Shadows check OFF
            if self.use_shadows:
                self.check_shadows.setChecked(False)
                self.use_shadows = False
                try: self.plotter.disable_shadows()
                except: pass

            # SSAO check OFF
            if self.use_ssao:
                self.check_ssao.setChecked(False)
                self.use_ssao = False
                try: 
                    if hasattr(self.plotter, 'disable_ssao'): self.plotter.disable_ssao()
                except: pass

        # 2. エフェクト(EDL/Shadows/SSAO)がONの場合 -> Depth PeelingをOFFにしてDisable
        elif exclude in ["edl", "shadows", "ssao"] and (self.use_edl or self.use_shadows or self.use_ssao):
            if self.use_depth_peeling:
                self.check_depth.setChecked(False)
                self.use_depth_peeling = False
                try:
                    if hasattr(self.plotter, 'disable_depth_peeling'):
                        self.plotter.disable_depth_peeling()
                except: pass
        

        # 3. 最後に有効/無効状態を更新 (グレーアウト処理) -> REMOVED
        # User requested to allow direct switching (auto-uncheck instead of disable).
        # We only keep the static checks (e.g. if VTK doesn't support it) which are done at init.
        
        self.blockSignals(False)

    # --- ATOM METHODS REMOVED ---

        
        # Save state of advanced rendering effects that get lost during plotter.clear()
        edl_was_enabled = self.use_edl
        shadows_were_enabled = self.use_shadows
        
        # Temporarily disable these effects to prevent FrameBufferObject cleanup errors
        if edl_was_enabled:
            try:
                self.plotter.disable_eye_dome_lighting()
            except Exception:
                pass
        
        if shadows_were_enabled:
            try:
                self.plotter.disable_shadows()
            except Exception:
                pass
        
        # Redraw molecule 3D and orbital
        # Redraw molecule using main window's draw_molecule_3d method
        if hasattr(self.mw, 'current_mol') and self.mw.current_mol:
            try:
                self.mw.draw_molecule_3d(self.mw.current_mol)
            except Exception as e:
                print(f"Error redrawing molecule: {e}")
        
        # Redraw orbital
        try:
            self.update_iso()
        except Exception as e:
            print(f"Error updating orbital: {e}")
        
        # Restore advanced rendering effects
        if shadows_were_enabled:
            try:
                self.plotter.enable_shadows()
            except Exception as e:
                print(f"Error restoring shadows: {e}")
        
        if edl_was_enabled:
            try:
                self.plotter.enable_eye_dome_lighting()
            except Exception as e:
                print(f"Error restoring EDL: {e}")
        
        # Final render to apply all effects
        try:
            self.plotter.render()
        except Exception:
            pass

    def on_texture_path_entered(self):
        path = self.line_env_path.text().strip()
        if path:
             # Just update internal state, load is triggered by Browse or manual logic if needed? 
             # Actually existing code handled path directly. Let's ensure consistency.
             self.env_texture_path = path

    # --- PRESET HANDLERS ---
    def update_preset_combo(self):
        current = self.combo_presets.currentText()
        self.combo_presets.blockSignals(True)
        self.combo_presets.clear()
        self.combo_presets.addItems(sorted(self.presets.keys()))
        self.combo_presets.setEditText(current)
        self.combo_presets.blockSignals(False)
        
    def on_preset_text_changed(self, text):
        pass # Just typing

    def on_preset_activated(self, index):
        name = self.combo_presets.currentText()
        if name in self.presets:
            # Redraw molecule and orbital before loading preset
            if hasattr(self.mw, 'current_mol') and self.mw.current_mol:
                try:
                    self.mw.draw_molecule_3d(self.mw.current_mol)
                except Exception:
                    pass
            
            try:
                self.update_iso()
            except Exception:
                pass
            
            self.load_preset_settings(self.presets[name])
            print(f"Loaded preset: {name}")

    def save_preset(self):
        # Show dialog to ask for preset name
        from PyQt6.QtWidgets import QInputDialog
        name, ok = QInputDialog.getText(
            self, 
            "Save Preset", 
            "Enter preset name:",
            text=self.combo_presets.currentText()
        )
        
        if not ok or not name.strip():
            return
        
        name = name.strip()
        
        # Prevent overwriting default presets
        if name in self.default_preset_names:
            QMessageBox.warning(self, "Cannot Overwrite Default", 
                              f"'{name}' is a built-in preset and cannot be overwritten.\nPlease use a different name.")
            return
        
        # Gather current settings
        current_data = {
            "isovalue": self.spin.value(),
            "color_p": list(self.color_p),
            "color_n": list(self.color_n),
            "opacity": self.opacity_spin.value(),
            "metallic": self.metallic,
            "roughness": self.roughness,
            "use_pbr": self.use_pbr,
            "env_texture_path": self.env_texture_path,
            "use_ssao": self.use_ssao,
            "use_depth_peeling": self.use_depth_peeling,
            "use_silhouette": self.use_silhouette,
            "use_aa": self.use_aa,
            "use_edl": self.use_edl,
            "edl_strength": self.edl_strength,
            "use_shadows": self.use_shadows,
            "light_intensity": self.light_intensity,
            "smooth_shading": self.check_smooth.isChecked(),
            "style": self.combo_style.currentText()
        }
        
        self.presets[name] = current_data
        self.update_preset_combo()
        self.combo_presets.setCurrentText(name)
        self.save_settings()
        print(f"Saved preset: {name}")
        
    def delete_preset(self):
        name = self.combo_presets.currentText().strip()
        
        # Prevent deleting default presets
        if name in self.default_preset_names:
            QMessageBox.warning(self, "Cannot Delete Default", 
                              f"'{name}' is a built-in preset and cannot be deleted.")
            return
        
        if name in self.presets:
            del self.presets[name]
            self.update_preset_combo()
            self.save_settings() # Update disk
            print(f"Deleted preset: {name}")
            
    def load_preset_settings(self, settings):
        # Apply settings dict to UI
        if "isovalue" in settings: self.spin.setValue(float(settings["isovalue"]))
        
        if "color_p" in settings:
            col = settings["color_p"]
            if isinstance(col, list): col = tuple(col)
            self.color_p = col
            self.btn_color_p.setStyleSheet(f"background-color: rgb{self.color_p}; border: 1px solid gray;")

        if "color_n" in settings:
            col = settings["color_n"]
            if isinstance(col, list): col = tuple(col)
            self.color_n = col
            self.btn_color_n.setStyleSheet(f"background-color: rgb{self.color_n}; border: 1px solid gray;")
            
        if "opacity" in settings:
            val = float(settings["opacity"])
            self.opacity_spin.setValue(val)
            self.opacity_slider.setValue(int(val * 100))
            
        if "metallic" in settings:
            self.metallic = float(settings["metallic"])
            self.slider_metallic.setValue(int(self.metallic * 100))
            
        if "roughness" in settings:
            self.roughness = float(settings["roughness"])
            self.slider_roughness.setValue(int(self.roughness * 100))
            
        if "use_pbr" in settings:
            self.use_pbr = bool(settings["use_pbr"])
            self.check_pbr.setChecked(self.use_pbr) # This triggers toggle signal -> Enable/Disable sliders

        if "env_texture_path" in settings:
            self.env_texture_path = settings["env_texture_path"]
            if self.env_texture_path:
                self.line_env_path.setText(self.env_texture_path)
                if os.path.exists(self.env_texture_path):
                    self.load_env_texture(self.env_texture_path)
            else:
                 self.line_env_path.clear()

        if "use_ssao" in settings: self.check_ssao.setChecked(bool(settings["use_ssao"]))
        if "use_depth_peeling" in settings: self.check_depth.setChecked(bool(settings["use_depth_peeling"]))
        if "use_silhouette" in settings: self.check_silhouette.setChecked(bool(settings["use_silhouette"]))
        if "use_aa" in settings: self.check_aa.setChecked(bool(settings["use_aa"]))
        
        if "edl_strength" in settings:
            self.edl_strength = float(settings["edl_strength"])
            self.slider_edl.setValue(int(self.edl_strength * 100))
            
        if "use_edl" in settings: self.check_edl.setChecked(bool(settings["use_edl"]))
        if "use_shadows" in settings: self.check_shadows.setChecked(bool(settings["use_shadows"]))
        if "light_intensity" in settings:
            self.light_intensity = float(settings["light_intensity"])
            self.slider_light.setValue(int(self.light_intensity * 100))
            
        if "smooth_shading" in settings:
            self.check_smooth.setChecked(bool(settings["smooth_shading"]))

        if "style" in settings:
            # Robust find
            text = settings["style"]
            idx = self.combo_style.findText(text)
            if idx < 0:
                # try case insensitive
                for i in range(self.combo_style.count()):
                    if self.combo_style.itemText(i).lower() == text.lower():
                        idx = i
                        break
            if idx >= 0: 
                self.combo_style.setCurrentIndex(idx)
        
        # Redraw orbital with loaded preset settings
        self.update_iso()

    def on_load_env_texture(self):
        fname, _ = QFileDialog.getOpenFileName(
            self, 
            "Load Environment Texture", 
            "", 
            "Images (*.png *.jpg *.jpeg);;All Files (*.*)"
        )
        if fname:
            # Check for spaces BEFORE attempting to load
            if ' ' in fname:
                QMessageBox.critical(
                    self,
                    "Invalid File Path",
                    "Texture file path contains spaces.\n\n"
                    "VTK cannot load files with spaces in the path.\n"
                    "Please move the file to a path without spaces and try again."
                )
                return
            self.load_env_texture(fname)

    def load_env_texture(self, path):
        """Load environment texture - PNG/JPG only."""
        progress = None
        try:
            # Validate file exists
            if not os.path.exists(path):
                raise FileNotFoundError(f"Texture file not found: {path}")
            
            # Check for spaces in path (VTK doesn't handle them well)
            if ' ' in path:
                raise ValueError("Texture file path contains spaces.\nVTK cannot load files with spaces in the path.\nPlease move the file to a path without spaces.")
            
            # Validate file extension
            ext = os.path.splitext(path)[1].lower()
            if ext not in ['.png', '.jpg', '.jpeg']:
                raise ValueError(f"Unsupported texture format: {ext}. Only PNG/JPG supported.")
            
            # Check file size to prevent freezing (limit to 10MB)
            file_size_mb = os.path.getsize(path) / (1024 * 1024)
            if file_size_mb > 10:
                reply = QMessageBox.question(
                    self, 
                    "Large Texture File",
                    f"Texture file is {file_size_mb:.1f}MB.\n"
                    f"Large textures may cause VTK to freeze.\n\n"
                    f"Continue loading?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply == QMessageBox.StandardButton.No:
                    return
            
            # Show progress dialog with cancel button
            from PyQt6.QtWidgets import QProgressDialog
            progress = QProgressDialog("Loading texture...", "Stop", 0, 0, self)
            progress.setWindowTitle("Loading Texture")
            progress.setWindowModality(Qt.WindowModality.WindowModal)
            progress.setMinimumDuration(0)
            progress.show()
            QCoreApplication.processEvents()
            
            # Check for cancellation
            if progress.wasCanceled():
                print("Texture loading cancelled by user")
                return
            
            # Load and apply texture with improved error handling
            try:
                texture = pv.read_texture(path)
            except Exception as tex_err:
                raise ValueError(f"Failed to read texture file: {tex_err}")
            
            # Check for cancellation
            if progress.wasCanceled():
                print("Texture loading cancelled by user")
                return
            
            progress.setLabelText("Applying texture to renderer...")
            QCoreApplication.processEvents()
            
            try:
                self.plotter.set_environment_texture(texture)
            except Exception as vtk_err:
                raise ValueError(f"VTK failed to apply texture (possibly incompatible format): {vtk_err}")
            
            self.env_texture_path = path
            self.line_env_path.setText(path)
            
            # Check for cancellation
            if progress.wasCanceled():
                print("Texture loading cancelled by user")
                return
            
            progress.setLabelText("Rendering...")
            QCoreApplication.processEvents()
            
            self.plotter.render()
            print(f"Loaded texture: {os.path.basename(path)}")
            
        except Exception as e:
            print(f"Failed to load texture: {e}")
            QMessageBox.critical(self, "Texture Loading Failed", 
                               f"Could not load texture:\n{os.path.basename(path)}\n\nError: {str(e)}")
            self.env_texture_path = ""
            # Re-raise exception so caller can handle it appropriately
            raise
        finally:
            # Always close progress dialog
            if progress:
                progress.close()

    def on_remove_env_texture(self):
        try:
            if hasattr(self.plotter, 'remove_environment_texture'):
                self.plotter.remove_environment_texture()
            elif hasattr(self.plotter, 'set_environment_texture'):
                # Fallback if remove_ doesn't exist but set_ does (unlikely given dir() output)
                try:
                    self.plotter.set_environment_texture(None)
                except:
                    print("Could not clear environment texture (set_environment_texture(None) failed).")
            
            self.env_texture_path = ""
            self.line_env_path.clear()
            self.plotter.render()
        except Exception as e:
            print(f"Error removing texture: {e}")

    def on_reset_all(self):
        """Resets ALL settings to defaults (not just advanced)."""
        # Block signals to prevent redundant expensive updates
        self.blockSignals(True)
        
        try:
            # 1. Reset Isovalue
            self.spin.setValue(0.02)
            
            # 2. Reset Colors
            self.color_p = (0, 0, 255) 
            self.color_n = (255, 0, 0) 
            self.btn_color_p.setStyleSheet(f"background-color: rgb{self.color_p};")
            self.btn_color_n.setStyleSheet(f"background-color: rgb{self.color_n};")
            
            # 3. Reset Opacity
            self.opacity_spin.setValue(0.4)
            self.opacity_slider.setValue(40)
            
            # 4. Reset Style
            self.combo_style.setCurrentIndex(0)  # Surface
            self.check_smooth.setChecked(True)
            
            # 5. Reset Advanced Settings
            self.check_pbr.setChecked(False)
            self.metallic = 0.5
            self.roughness = 0.5
            self.slider_metallic.setValue(50)
            self.slider_roughness.setValue(50)
            
            # 6. Reset Effects
            # --- FIX START: シグナルブロック中は手動で無効化が必要 ---
            
            # SSAO
            self.check_ssao.setChecked(False)
            if self.use_ssao:
                try: 
                    if hasattr(self.plotter, 'disable_ssao'): self.plotter.disable_ssao()
                except: pass
                self.use_ssao = False

            # Depth Peeling
            self.check_depth.setChecked(False)
            if self.use_depth_peeling:
                try:
                    if hasattr(self.plotter, 'disable_depth_peeling'): self.plotter.disable_depth_peeling()
                except: pass
                self.use_depth_peeling = False

            # Silhouette
            self.check_silhouette.setChecked(False)
            self.use_silhouette = False
            
            # EDL (Eye Dome Lighting)
            self.check_edl.setChecked(False)
            self.slider_edl.setValue(20) # 0.2
            if self.use_edl:
                try:
                    self.plotter.disable_eye_dome_lighting()
                except: pass
                self.use_edl = False

            # Anti-Aliasing
            self.check_aa.setChecked(False)
            if self.use_aa:
                try:
                    self.plotter.disable_anti_aliasing()
                except: pass
                self.use_aa = False

            # Shadows
            self.check_shadows.setChecked(False)
            self.slider_light.setValue(100) # 1.0
            if self.use_shadows:
                try:
                    self.plotter.disable_shadows()
                except: pass
                self.use_shadows = False
            
            # --- FIX END ---


             
            # 7. Clear Environment Texture
            self.on_remove_env_texture()
            
            # 8. Trigger Update
            self.blockSignals(False)
            
            # 確実に描画更新
            self.update_iso()
            self.plotter.render()

            print("All settings have been reset to defaults.")
            
        except Exception as e:
            self.blockSignals(False)
            print(f"Error resetting settings: {e}")

    def closeEvent(self, event):
        self.save_settings()
        super().closeEvent(event)

    def close_plugin(self):
        self.save_settings()
        try:
            # --- FIX START: FBO Warning対策 ---
            # Plotterをクリアする前に、レンダリングパス(Shadows, EDL)を
            # 明示的に無効化しないと、バインドされていないFBOへのアクセス警告が出る
            if self.use_edl:
                try:
                    self.plotter.disable_eye_dome_lighting()
                except: pass
            
            if self.use_shadows:
                try:
                    self.plotter.disable_shadows()
                except: pass
            
            if self.use_ssao:
                try:
                    if hasattr(self.plotter, 'disable_ssao'):
                        self.plotter.disable_ssao()
                except: pass
            # --- FIX END ---

            # Full cleanup
            self.mw.plotter.clear()
            self.mw.current_mol = None
            self.mw.current_file_path = None
            self.mw.plotter.render()
             
            # Restore UI state
            if hasattr(self.mw, 'restore_ui_for_editing'):
                self.mw.restore_ui_for_editing()
        except Exception as e: 
            print(f"Error closing plugin: {e}")
        
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
        layout.addWidget(QLabel("Please specify correct charge or skip chemistry check."))
        
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
        
        # Calculate max absolute value in data for dynamic scaling
        try:
            # Data is stored in grid.point_data["values"]
            # "make sure that the data is only in the plot data"
            flat_data = grid.point_data["values"]
            if len(flat_data) > 0:
                data_max = float(np.max(np.abs(flat_data)))
            else:
                data_max = 1.0
        except:
             data_max = 1.0

        viewer = CubeViewerWidget(main_window, dock, grid, data_max=data_max)
        dock.setWidget(viewer)
        
        main_window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        
        
        # main_window.plotter.reset_camera() # Handled in initial_update of widget

    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Plugin Error: {e}")

def run(mw):
    if Chem is None:
        QMessageBox.critical(mw, "Error", "RDKit is required for this plugin.")
        return

    fname, _ = QFileDialog.getOpenFileName(mw, "Open Gaussian Cube File", "", "Cube Files (*.cube *.cub);;All Files (*)")
    if fname:
        open_cube_viewer(mw, fname)

def initialize(context):
    """
    New Plugin System Entry Point
    """
    mw = context.get_main_window()

    def open_cube_wrapper(fname):
        open_cube_viewer(mw, fname)

    # 1. Register File Opener (Handle File > Import)
    context.register_file_opener('.cube', open_cube_wrapper)
    context.register_file_opener('.cub', open_cube_wrapper)

    # 2. Register Drop Handler (for robustness)
    # The system iterates drop handlers. Return True if we handled it.
    def drop_handler(file_path):
        if file_path.lower().endswith(('.cube', '.cub')):
            open_cube_viewer(mw, file_path)
            return True
        return False

    if hasattr(context, 'register_drop_handler'):
        context.register_drop_handler(drop_handler, priority=10)




