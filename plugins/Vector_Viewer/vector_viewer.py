
import sys
import numpy as np
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, 
    QLineEdit, QSpinBox, QDoubleSpinBox, QCheckBox, 
    QComboBox, QPushButton, QFileDialog, QMessageBox, QGroupBox,
    QColorDialog
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor

# --- Plugin Metadata ---
PLUGIN_NAME = "Vector Viewer"
PLUGIN_VERSION = "2025.12.31"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Visualizes vectors and exports PNGs."

class VectorViewerPlugin(QWidget):
    def __init__(self, interface):
        """
        :param interface: MoleditPy Main Window Interface
        """
        parent = interface if hasattr(interface, 'windowTitle') else None 
        super().__init__(parent)
        self.setWindowFlags(Qt.WindowType.Window)
        self.setWindowTitle("Vector Viewer")
        
        self.interface = interface
        self.vis_actor = None # Store the arrow actor
        self.arrow_color = QColor("green")
        
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout(self)
        
        # --- 1. Vector Input & Update ---
        vec_group = QGroupBox("Vector")
        vec_layout = QVBoxLayout()
        
        self.vec_input = QLineEdit()
        self.vec_input.setPlaceholderText("Paste vector (e.g. '1.5 0.0 -2.3')")
        self.vec_input.returnPressed.connect(self.update_visualization)
        
        self.update_btn = QPushButton("Update View")
        self.update_btn.clicked.connect(self.update_visualization)
        
        vec_layout.addWidget(QLabel("Vector (x y z), e.g. dipole moment etc.:"))
        vec_layout.addWidget(self.vec_input)
        
        self.reverse_chk = QCheckBox("Reverse Vector")
        #self.reverse_chk.toggled.connect(self.update_visualization)
        vec_layout.addWidget(self.reverse_chk)

        vec_layout.addWidget(self.update_btn)
        
        vec_group.setLayout(vec_layout)
        main_layout.addWidget(vec_group)
        
        # --- 2. Appearance Settings ---
        app_group = QGroupBox("Appearance")
        app_layout = QGridLayout()
        
        # Scale
        app_layout.addWidget(QLabel("Scale:"), 0, 0)
        self.scale_spin = QDoubleSpinBox()
        self.scale_spin.setRange(0.1, 10.0)
        self.scale_spin.setSingleStep(0.1)
        self.scale_spin.setValue(1.0)
        self.scale_spin.valueChanged.connect(self.update_visualization)
        app_layout.addWidget(self.scale_spin, 0, 1)
        
        # Resolution
        app_layout.addWidget(QLabel("Resolution:"), 1, 0)
        self.res_spin = QSpinBox()
        self.res_spin.setRange(4, 100)
        self.res_spin.setValue(20)
        self.res_spin.valueChanged.connect(self.update_visualization)
        app_layout.addWidget(self.res_spin, 1, 1)
        
        # Opacity
        app_layout.addWidget(QLabel("Opacity:"), 2, 0)
        self.opacity_spin = QDoubleSpinBox()
        self.opacity_spin.setRange(0.0, 1.0)
        self.opacity_spin.setSingleStep(0.1)
        self.opacity_spin.setValue(0.5) # Default Opaque
        self.opacity_spin.valueChanged.connect(self.update_visualization)
        app_layout.addWidget(self.opacity_spin, 2, 1)
        
        # Color Picker
        app_layout.addWidget(QLabel("Color:"), 3, 0)
        self.color_btn = QPushButton("Select Color")
        self.update_color_btn_style()
        self.color_btn.clicked.connect(self.choose_color)
        app_layout.addWidget(self.color_btn, 3, 1)
        
        app_group.setLayout(app_layout)
        main_layout.addWidget(app_group)
        
        # --- 4. Export ---
        exp_group = QGroupBox("Export")
        exp_layout = QVBoxLayout()
        
        self.trans_chk = QCheckBox("Transparent Background")
        self.trans_chk.setChecked(True)
        exp_layout.addWidget(self.trans_chk)
        
        self.exp_btn = QPushButton("Export to PNG")
        self.exp_btn.clicked.connect(self.export_image)
        exp_layout.addWidget(self.exp_btn)
        
        exp_group.setLayout(exp_layout)
        main_layout.addWidget(exp_group)
        
        main_layout.addStretch()

    def update_color_btn_style(self):
        name = self.arrow_color.name()
        # Text color contrasting
        text_col = "black" if self.arrow_color.lightness() > 128 else "white"
        self.color_btn.setStyleSheet(f"background-color: {name}; color: {text_col}; font-weight: bold;")
        self.color_btn.setText(name.upper())

    def choose_color(self):
        col = QColorDialog.getColor(self.arrow_color, self, "Choose Arrow Color")
        if col.isValid():
            self.arrow_color = col
            self.update_color_btn_style()
            self.update_visualization()

    def get_com(self):
        """Calculate Center of Mass (COM) of the current molecule."""
        rd_mol = getattr(self.interface, 'current_mol', None)
        if rd_mol is None and hasattr(self.interface, 'get_molecule'):
            rd_mol = self.interface.get_molecule()
            
        if rd_mol:
            try:
                conf = rd_mol.GetConformer()
                positions = [conf.GetAtomPosition(i) for i in range(rd_mol.GetNumAtoms())]
                coords = np.array([[p.x, p.y, p.z] for p in positions])
                if len(coords) > 0:
                    return np.mean(coords, axis=0)
            except:
                pass
        return np.array([0., 0., 0.])

    def update_visualization(self):
        """Parse vector and draw arrow in PyVista plotter."""
        if not hasattr(self.interface, 'plotter'):
            return # No plotter available
            
        text = self.vec_input.text().strip()
        if not text:
            return # Empty
            
        # Parse Vector
        try:
            # Handle comma separated or space separated
            cleaned = text.replace(',', ' ')
            parts = cleaned.split()
            if len(parts) < 3:
                raise ValueError
            vec = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
            
            if self.reverse_chk.isChecked():
                vec = -vec
        except ValueError:
            # Don't annoy user while typing, just return or subtle hint
            return

        # Settings
        scale = self.scale_spin.value()
        res = self.res_spin.value()
        opacity = self.opacity_spin.value()
        color = self.arrow_color.name()
        
        scaled_vec = vec * scale
        mag = np.linalg.norm(scaled_vec)
        
        if mag < 1e-6:
            return 
            
        com = self.get_com()
        
        # Geometry Logic: Center of arrow at COM
        # Start = COM - (Vector / 2)
        start = com - (scaled_vec / 2.0)
        
        import pyvista as pv
        plotter = self.interface.plotter
        
        # Remove old actor
        if self.vis_actor:
            plotter.remove_actor(self.vis_actor)
            self.vis_actor = None
            
        arrow = pv.Arrow(start=start, direction=scaled_vec, scale=mag, 
                         tip_length=0.2, tip_radius=0.1, shaft_radius=0.04, 
                         shaft_resolution=res, tip_resolution=res)
        
        # opacity引数は add_mesh で指定
        self.vis_actor = plotter.add_mesh(arrow, color=color, opacity=opacity)
        plotter.render()

    def export_image(self):
        """Export current view to PNG."""
        if not hasattr(self.interface, 'plotter'):
            QMessageBox.warning(self, "Error", "No 3D plotter found.")
            return

        filename, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "PNG Files (*.png)")
        if not filename:
            return
            
        transparent = self.trans_chk.isChecked()
        
        try:
            self.interface.plotter.screenshot(filename, transparent_background=transparent)
            QMessageBox.information(self, "Success", f"Saved to {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save image:\n{e}")

def run(interface):
    """Entry Point"""
    if not hasattr(interface, 'vector_viewer_window'):
        interface.vector_viewer_window = VectorViewerPlugin(interface)
        interface.vector_viewer_window.resize(300, 450)
    
    interface.vector_viewer_window.show()
    interface.vector_viewer_window.raise_()

if __name__ == "__main__":
    # Standalone Test
    from PyQt6.QtWidgets import QApplication
    
    class MockInterface:
        def __init__(self):
            self.plotter = None
            self.init_plotter()
            self.current_mol = None
            self._create_sample_molecule()
            
        def _create_sample_molecule(self):
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem
                m = Chem.MolFromSmiles('C') # Methane
                m = Chem.AddHs(m)
                AllChem.EmbedMolecule(m, randomSeed=42) 
                AllChem.MMFFOptimizeMolecule(m)
                self.current_mol = m
            except ImportError:
                print("RDKit not found")

        def init_plotter(self):
            try:
                import pyvista as pv
                self.plotter = pv.Plotter()
                self.plotter.add_text("Mock Plotter", position='upper_left')
                # Add a dummy sphere to represent molecule
                self.plotter.add_mesh(pv.Sphere(radius=0.5), color='white')
                self.plotter.show(auto_close=False, interactive_update=True)
            except ImportError:
                print("PyVista not found")
                
        def get_molecule(self):
            return self.current_mol

    app = QApplication(sys.argv)
    
    # Check PyVista
    try:
        import pyvista
    except ImportError:
        print("This tool requires 'pyvista'.")
        sys.exit(1)

    mock = MockInterface()
    plugin = VectorViewerPlugin(mock)
    plugin.show()
    
    print("Test: Paste '1 1 1' into vector input and press Enter.")
    
    sys.exit(app.exec())
