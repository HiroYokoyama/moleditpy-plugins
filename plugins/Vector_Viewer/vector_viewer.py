import numpy as np
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QGridLayout,
    QLabel,
    QLineEdit,
    QSpinBox,
    QDoubleSpinBox,
    QCheckBox,
    QPushButton,
    QFileDialog,
    QMessageBox,
    QGroupBox,
    QColorDialog,
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor
import logging

import pyvista as pv

# --- Plugin Metadata ---
PLUGIN_NAME = "Vector Viewer"
PLUGIN_VERSION = "2026.10.18"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=5.0.0a1, <6.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Visualizes vectors (e.g. dipole moment) and exports PNGs."


class VectorViewerPlugin(QWidget):
    """
    Plugin for visualizing vectors (e.g. dipole moments) in the 3D scene.
    Refactored for MoleditPy V3.0 API.
    """

    def __init__(self, context):
        """
        :param context: MoleditPy PluginContext
        """
        super().__init__(parent=context.get_main_window())
        self.setWindowFlags(Qt.WindowType.Window)
        self.setWindowTitle("Vector Viewer")

        self.context = context
        self.vis_actor = None
        self._actor_plotter = None
        self._molecule = context.current_molecule
        self._anchor = None
        self.arrow_color = QColor("green")

        self.init_ui()
        self._molecule_timer = QTimer(self)
        self._molecule_timer.setInterval(250)
        self._molecule_timer.timeout.connect(self._check_molecule)
        self._molecule_timer.start()

        # Register window for V3 lifecycle management
        self.context.register_window("main_panel", self)

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
        self.validation_label = QLabel()
        vec_layout.addWidget(self.validation_label)

        self.reverse_chk = QCheckBox("Reverse Vector")
        self.reverse_chk.toggled.connect(self.update_visualization)
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
        self.opacity_spin.setValue(0.5)  # Default Opaque
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
        self.color_btn.setStyleSheet(
            f"background-color: {name}; color: {text_col}; font-weight: bold;"
        )
        self.color_btn.setText(name.upper())

    def choose_color(self):
        col = QColorDialog.getColor(self.arrow_color, self, "Choose Arrow Color")
        if col.isValid():
            self.arrow_color = col
            self.update_color_btn_style()
            self.update_visualization()

    def get_com(self):
        """Calculate Center of Mass (COM) of the current molecule."""
        rd_mol = self.context.current_molecule
        if rd_mol and rd_mol.GetNumAtoms() > 0:
            try:
                conf = rd_mol.GetConformer()
                positions = [
                    conf.GetAtomPosition(i) for i in range(rd_mol.GetNumAtoms())
                ]
                coords = np.array([[p.x, p.y, p.z] for p in positions])
                masses = np.array([atom.GetMass() for atom in rd_mol.GetAtoms()], dtype=float)
                if len(masses) != len(coords) or not np.all(np.isfinite(coords)) or not np.all(np.isfinite(masses)) or np.any(masses <= 0):
                    raise ValueError("Molecule coordinates or atomic masses are invalid")
                return np.average(coords, axis=0, weights=masses)
            except (RuntimeError, AttributeError, IndexError, TypeError, ValueError) as _e:
                raise ValueError("The molecule needs finite 3D coordinates and atomic masses") from _e
        return np.array([0.0, 0.0, 0.0])

    def _clear_arrow(self):
        actor, self.vis_actor = self.vis_actor, None
        plotter = self._actor_plotter or self.context.plotter
        self._actor_plotter = None
        self._anchor = None
        if actor is not None and plotter is not None:
            try:
                plotter.remove_actor(actor)
                plotter.render()
            except (RuntimeError, AttributeError, KeyError, TypeError, ValueError):
                logging.exception("Could not remove Vector Viewer arrow")

    def _check_molecule(self):
        molecule = self.context.current_molecule
        if molecule is not self._molecule:
            self._molecule = molecule
            self._clear_arrow()
            self.vec_input.clear()
            self.validation_label.setText("Molecule changed. Enter its vector to display it.")
        elif self.vis_actor is not None:
            try:
                anchor = self.get_com()
                if self._anchor is not None and not np.allclose(anchor, self._anchor):
                    self.update_visualization()
            except ValueError:
                self._clear_arrow()

    def update_visualization(self):
        """Draw exactly three finite vector components, centered at molecular COM."""
        if self.vis_actor is not None and self._molecule is not self.context.current_molecule:
            self._check_molecule()
            return
        plotter = self.context.plotter
        if plotter is None:
            self._clear_arrow()
            return
        text = self.vec_input.text().strip()
        self.validation_label.clear()
        if not text:
            self._clear_arrow()
            return
        try:
            parts = text.replace(",", " ").split()
            if len(parts) != 3:
                raise ValueError("Enter exactly three finite components.")
            vec = np.array([float(part.replace("D", "E").replace("d", "e")) for part in parts])
            if not np.all(np.isfinite(vec)):
                raise ValueError("Enter exactly three finite components.")
            if self.reverse_chk.isChecked():
                vec = -vec
            with np.errstate(over="ignore", invalid="ignore"):
                scaled_vec = vec * self.scale_spin.value()
                magnitude = np.linalg.norm(scaled_vec)
            if not np.isfinite(magnitude):
                raise ValueError("The scaled vector is too large to display.")
            if magnitude < 1e-6:
                self._clear_arrow()
                return
            anchor = self.get_com()
            arrow = pv.Arrow(start=anchor - scaled_vec / 2, direction=scaled_vec,
                             scale=magnitude, tip_length=0.2, tip_radius=0.1,
                             shaft_radius=0.04, shaft_resolution=self.res_spin.value(),
                             tip_resolution=self.res_spin.value())
            self._clear_arrow()
            self.vis_actor = plotter.add_mesh(arrow, color=self.arrow_color.name(),
                                             opacity=self.opacity_spin.value(), name="vector_viewer_arrow")
            self._actor_plotter = plotter
            self._molecule = self.context.current_molecule
            self._anchor = anchor
            plotter.render()
        except (ValueError, RuntimeError, TypeError, AttributeError) as exc:
            self._clear_arrow()
            self.validation_label.setText(str(exc))
            logging.debug("Could not display vector: %s", exc)

    def closeEvent(self, event):
        self._molecule_timer.stop()
        self._clear_arrow()
        super().closeEvent(event)

    def showEvent(self, event):
        super().showEvent(event)
        self._check_molecule()
        if self.vec_input.text().strip():
            self.update_visualization()
        self._molecule_timer.start()

    def export_image(self):
        """Export current view to PNG."""
        plotter = self.context.plotter
        if plotter is None:
            QMessageBox.warning(self, "Error", "No 3D plotter found.")
            return

        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Image", "", "PNG Files (*.png)"
        )
        if not filename:
            return

        transparent = self.trans_chk.isChecked()

        try:
            plotter.screenshot(filename, transparent_background=transparent)
            self.context.show_status_message(f"Saved to {filename}")
        except (OSError, RuntimeError, AttributeError, KeyError, TypeError, ValueError) as e:
            QMessageBox.critical(self, "Error", f"Failed to save image:\n{e}")


_launch_fn = None


def initialize(context):
    """MoleditPy Plugin Entry Point (V3.0)"""
    global _launch_fn

    def launch():
        win = context.get_window("main_panel")
        if win:
            win.show()
            win.raise_()
            win.activateWindow()
            return
        win = VectorViewerPlugin(context)
        win.show()

    _launch_fn = launch
    context.show_status_message(f"{PLUGIN_NAME} Loaded.")


def run(mw):
    if hasattr(mw, "host"):
        mw = mw.host
    if _launch_fn:
        _launch_fn()
