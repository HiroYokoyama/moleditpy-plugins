# -*- coding: utf-8 -*-
"""
Animated XYZ Player Plugin for MoleditPy

Allows loading and playing multi-frame XYZ files (e.g., MD trajectories).
"""

import os
import time
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton,
    QSlider, QLabel, QSpinBox, QFileDialog, QWidget,
    QMessageBox, QDockWidget
)
from PyQt6.QtCore import Qt, QTimer, QSize
from rdkit import Chem
from rdkit.Chem import AllChem, rdGeometry

PLUGIN_NAME = "Animated XYZ Player"

class AnimatedXYZPlayer(QDialog):
    def __init__(self, main_window):
        super().__init__(main_window)
        self.mw = main_window
        self.setWindowTitle("Animated XYZ Player")
        self.setWindowFlags(Qt.WindowType.Window) # Make it a separate window acting like a tool
        self.resize(400, 150)

        # Data
        self.frames = [] # List of list of (symbol, x, y, z)
        self.current_frame_idx = 0
        self.base_mol = None # RDKit Mol with topology
        self.is_playing = False
        self.fps = 10

        # UI Layout
        layout = QVBoxLayout(self)

        # File controls
        file_layout = QHBoxLayout()
        self.btn_load = QPushButton("Load XYZ")
        self.btn_load.clicked.connect(self.load_file)
        self.lbl_file = QLabel("No file loaded")
        file_layout.addWidget(self.btn_load)
        file_layout.addWidget(self.lbl_file)
        layout.addLayout(file_layout)

        # Status
        self.lbl_status = QLabel("Frame: 0 / 0")
        layout.addWidget(self.lbl_status)

        # Slider
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self.on_slider_changed)
        layout.addWidget(self.slider)

        # Playback controls
        ctrl_layout = QHBoxLayout()
        
        self.btn_prev = QPushButton("<<")
        self.btn_prev.clicked.connect(self.prev_frame)
        self.btn_prev.setEnabled(False)
        
        self.btn_play = QPushButton("Play")
        self.btn_play.clicked.connect(self.toggle_play)
        self.btn_play.setEnabled(False)
        
        self.btn_next = QPushButton(">>")
        self.btn_next.clicked.connect(self.next_frame)
        self.btn_next.setEnabled(False)

        ctrl_layout.addWidget(self.btn_prev)
        ctrl_layout.addWidget(self.btn_play)
        ctrl_layout.addWidget(self.btn_next)
        layout.addLayout(ctrl_layout)

        # FPS control
        fps_layout = QHBoxLayout()
        fps_layout.addWidget(QLabel("FPS:"))
        self.spin_fps = QSpinBox()
        self.spin_fps.setRange(1, 60)
        self.spin_fps.setValue(self.fps)
        self.spin_fps.valueChanged.connect(self.set_fps)
        fps_layout.addWidget(self.spin_fps)
        layout.addLayout(fps_layout)

        # Timer
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.on_timer)

    def load_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Animated XYZ", "", "XYZ Files (*.xyz);;All Files (*)"
        )
        if not file_path:
            return

        try:
            frames = self.parse_multi_frame_xyz(file_path)
            if not frames:
                QMessageBox.warning(self, "Error", "No valid frames found in XYZ file.")
                return
            
            self.frames = frames
            self.current_frame_idx = 0
            self.lbl_file.setText(os.path.basename(file_path))
            self.slider.setRange(0, len(frames) - 1)
            self.slider.setValue(0)
            self.slider.setEnabled(True)
            self.btn_prev.setEnabled(True)
            self.btn_play.setEnabled(True)
            self.btn_next.setEnabled(True)
            
            self.create_base_molecule()
            self.update_view()
            self.update_status()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file:\n{e}")

    def parse_multi_frame_xyz(self, file_path):
        """
        Parses a concatenated XYZ file.
        Returns a list of frames, where each frame is tuple of (atoms, coords).
        Wait, we just need coordinates if topology is constant.
        Let's store: [ {'symbols': [str], 'coords': [(x,y,z)]}, ... ]
        """
        frames = []
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        i = 0
        n_lines = len(lines)
        while i < n_lines:
            line = lines[i].strip()
            if not line:
                i += 1
                continue
            
            try:
                num_atoms = int(line)
            except ValueError:
                # Might be a blank line or garbage
                i += 1
                continue
            
            # Start of a frame
            # i = atom count
            # i+1 = comment
            # i+2 ... i+2+num_atoms = atoms
            
            if i + 2 + num_atoms > n_lines:
                break # Incomplete frame
            
            comment = lines[i+1].strip()
            frame_atoms = []
            frame_coords = []
            
            start_data = i + 2
            for j in range(num_atoms):
                parts = lines[start_data + j].split()
                if len(parts) >= 4:
                    sym = parts[0]
                    try:
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        frame_atoms.append(sym)
                        frame_coords.append((x, y, z))
                    except ValueError:
                        pass
            
            frames.append({
                'symbols': frame_atoms,
                'coords': frame_coords,
                'comment': comment
            })
            
            i = start_data + num_atoms
            
        return frames

    def create_base_molecule(self):
        """
        Creates the RDKit Mol object from the first frame and establishes topology.
        """
        if not self.frames:
            return

        frame0 = self.frames[0]
        mol = Chem.RWMol()
        
        # Add atoms
        for sym in frame0['symbols']:
            # Handle unknown symbols or numbers
            try:
                atom = Chem.Atom(sym)
            except:
                atom = Chem.Atom('C') # Fallback
            mol.AddAtom(atom)
            
        # Add conformer
        conf = Chem.Conformer(mol.GetNumAtoms())
        for idx, (x, y, z) in enumerate(frame0['coords']):
            conf.SetAtomPosition(idx, rdGeometry.Point3D(x, y, z))
        mol.AddConformer(conf)

        # Estimate bonds (topology)
        # We try to use the main window's helper function if available, 
        # otherwise we manually do simple distance check or leave it unconnected
        if hasattr(self.mw, 'estimate_bonds_from_distances'):
            self.mw.estimate_bonds_from_distances(mol)
        
        self.base_mol = mol.GetMol()
        
        # Set as current mol in main window so it can be drawn
        self.mw.current_mol = self.base_mol

        # Ensure 3D capabilities are on
        if hasattr(self.mw, '_enter_3d_viewer_ui_mode'):
            self.mw._enter_3d_viewer_ui_mode()
        
        # Reset camera on first load
        if hasattr(self.mw, 'plotter'):
             self.mw.plotter.reset_camera()

    def update_view(self):
        if not self.frames or self.base_mol is None:
            return
        
        if self.current_frame_idx >= len(self.frames):
            self.current_frame_idx = 0
            
        frame = self.frames[self.current_frame_idx]
        
        # Update conformer positions
        # Assuming topology (atom count/order) hasn't changed
        conf = self.base_mol.GetConformer()
        coords = frame['coords']
        
        # Safety check for atom count mismatch
        if len(coords) != self.base_mol.GetNumAtoms():
            # If atom count changes, we might need to recreate the molecule
            # For this simple plugin, we'll just ignore or warn
            # print("Atom count mismatch in animation")
            pass
        else:
             for idx, (x, y, z) in enumerate(coords):
                conf.SetAtomPosition(idx, rdGeometry.Point3D(x, y, z))
        
        # Redraw
        # This calls main_window.draw_molecule_3d which rebuilds the scene's actors
        if hasattr(self.mw, 'draw_molecule_3d'):
            self.mw.draw_molecule_3d(self.base_mol)
            # Update frame comment/title if possible?
            self.mw.statusBar().showMessage(f"Frame {self.current_frame_idx+1}/{len(self.frames)}: {frame['comment']}")

    def update_status(self):
        self.lbl_status.setText(f"Frame: {self.current_frame_idx + 1} / {len(self.frames)}")
        self.slider.blockSignals(True)
        self.slider.setValue(self.current_frame_idx)
        self.slider.blockSignals(False)

    def on_slider_changed(self, value):
        self.current_frame_idx = value
        self.update_view()
        self.update_status()

    def toggle_play(self):
        self.is_playing = not self.is_playing
        if self.is_playing:
            self.btn_play.setText("Pause")
            self.timer.start(int(1000 / self.fps))
        else:
            self.btn_play.setText("Play")
            self.timer.stop()
            # Ensure the main window knows this is the generic current molecule 
            # so the user can use File->Save As... to export the current frame.
            self.mw.current_mol = self.base_mol

    def next_frame(self):
        self.current_frame_idx = (self.current_frame_idx + 1) % len(self.frames)
        self.update_view()
        self.update_status()

    def prev_frame(self):
        self.current_frame_idx = (self.current_frame_idx - 1) % len(self.frames)
        self.update_view()
        self.update_status()

    def on_timer(self):
        self.next_frame()

    def set_fps(self, value):
        self.fps = value
        if self.is_playing:
            self.timer.start(int(1000 / self.fps))
            
    def closeEvent(self, event):
        self.timer.stop()
        
        # Clear the main window view
        try:
            if hasattr(self.mw, 'plotter'):
                self.mw.plotter.clear()
        except:
            pass

        try:
            self.mw.current_mol = None
        except:
            pass
            
        # Exit 3D mode and restore 2D editor UI
        # We try calling it directly.
        try:
            self.mw.restore_ui_for_editing()
        except Exception as e:
            print(f"Error restoring UI: {e}")

        # Force a re-render/clear of the generic 3D draw function
        try:
             self.mw.draw_molecule_3d(None)
        except:
             pass
             
        # Remove reference from main window so next run starts fresh check
        if hasattr(self.mw, '_plugin_animated_xyz_player'):
            del self.mw._plugin_animated_xyz_player

        super().closeEvent(event)

def run(main_window):
    # Always close/destroy old instance to reset variables and state
    if hasattr(main_window, '_plugin_animated_xyz_player'):
        try:
            main_window._plugin_animated_xyz_player.close()
        except:
            pass
        # Depending on if closeEvent did its job or not, strictly remove ref
        if hasattr(main_window, '_plugin_animated_xyz_player'):
             del main_window._plugin_animated_xyz_player

    # Create fresh instance
    player = AnimatedXYZPlayer(main_window)
    main_window._plugin_animated_xyz_player = player
    player.show()
    player.raise_()
    player.activateWindow()

