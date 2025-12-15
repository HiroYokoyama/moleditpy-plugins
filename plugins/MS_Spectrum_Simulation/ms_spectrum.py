# -*- coding: utf-8 -*-
import sys
import math
from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QWidget, QMessageBox, QPushButton, QFileDialog
from PyQt6.QtGui import QPainter, QPen, QBrush, QColor, QFont, QPalette, QLinearGradient, QGradient
from PyQt6.QtCore import Qt, QRectF, QPointF

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    Chem = None

PLUGIN_NAME = "MS Spectrum Simulation"

class MSSpectrumDialog(QDialog):
    def __init__(self, mol, parent=None):
        super().__init__(parent)
        self.setWindowTitle("MS Spectrum Simulation")
        self.resize(400, 300) # Smaller and resizable
        self.mol = mol
        self.peaks = self._calculate_peaks()
        
        # Clean white look
        self.setStyleSheet("""
            QDialog {
                background-color: #ffffff;
                color: #000000;
            }
            QLabel {
                color: #333333;
                font-size: 14px;
                font-family: 'Segoe UI', sans-serif;
                padding: 0px;
            }
        """)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 5, 10, 10) # Very tight margins
        layout.setSpacing(5) # Minimal spacing
        
        # Top Info Panel: Vertical Layout
        info_layout = QVBoxLayout()
        info_layout.setContentsMargins(20, 0, 20, 0) 
        info_layout.setSpacing(2) # Tight vertical spacing
        if self.mol:
            formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
            mw = Descriptors.MolWt(self.mol)
            exact_mw = Descriptors.ExactMolWt(self.mol)
            
            l1 = QLabel(f"<b>Formula:</b> {formula}")
            l2 = QLabel(f"<b>Avg Mass:</b> {mw:.6f}")
            l3 = QLabel(f"<b>Exact Mass:</b> {exact_mw:.6f}")
            
            info_layout.addWidget(l1)
            info_layout.addWidget(l2)
            info_layout.addWidget(l3)
        else:
            info_layout.addWidget(QLabel("No molecule loaded."))
        
        layout.addLayout(info_layout)
        
        # Container for the plot
        plot_container = QWidget()
        plot_container.setStyleSheet("background-color: #ffffff; border: 1px solid #cccccc; border-radius: 4px;")
        plot_layout = QVBoxLayout(plot_container)
        plot_layout.setContentsMargins(1, 1, 1, 1)
        
        self.plot_widget = HistogramWidget(self.peaks)
        plot_layout.addWidget(self.plot_widget)
        
        layout.addWidget(plot_container, 1)

        # Export Button
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        self.export_btn = QPushButton("Export to Image")
        self.export_btn.clicked.connect(self.export_image)
        btn_layout.addWidget(self.export_btn)
        layout.addLayout(btn_layout)

    def export_image(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Save Spectrum Image", "spectrum.png", "Images (*.png *.jpg)")
        if filename:
            pixmap = self.plot_widget.grab()
            if pixmap.save(filename):
                QMessageBox.information(self, "Success", f"Image saved to {filename}")
            else:
                QMessageBox.critical(self, "Error", "Failed to save image.")

    def _calculate_peaks(self):
        if not self.mol:
            return []

        # 1. Count elements
        atoms = self.mol.GetAtoms()
        element_counts = {}
        for atom in atoms:
            sym = atom.GetSymbol()
            element_counts[sym] = element_counts.get(sym, 0) + 1
            
        # Add implicit hydrogens
        num_implicit_h = 0
        for atom in atoms:
            num_implicit_h += atom.GetTotalNumHs()
        if num_implicit_h > 0:
             element_counts['H'] = element_counts.get('H', 0) + num_implicit_h

        # 2. Dynamic Isotope Data Fetching from RDKit
        pt = Chem.GetPeriodicTable()
        
        # Initialize Distribution [ (mass, probability), ... ]
        current_dist = [(0.0, 1.0)]
        
        for sym, count in element_counts.items():
            atomic_num = pt.GetAtomicNumber(sym)
            atom_iso_dist = []
            
            # Find isotopes by scanning around the most common isotope mass number
            # Valid range is usually within +/- a few neutrons
            center_mass = pt.GetMostCommonIsotope(atomic_num)
            # Scan a generous range to catch minor isotopes (e.g. S-36 is +4 from S-32)
            scan_range = range(max(1, center_mass - 5), center_mass + 10)
            
            for m in scan_range:
                try:
                    abundance = pt.GetAbundanceForIsotope(atomic_num, m)
                    if abundance > 0.00001: # Cutoff for very rare isotopes to save perf
                        exact_mass = pt.GetMassForIsotope(atomic_num, m)
                        atom_iso_dist.append((exact_mass, abundance))
                except RuntimeError:
                    pass # Invalid isotope check
            
            # Normalize (RDKit abundances usually sum to 1.0, but good to ensure)
            total_p = sum(p for m, p in atom_iso_dist)
            if total_p == 0: continue
            atom_iso_dist = [(m, p / total_p) for m, p in atom_iso_dist]
            
            # 3. Convolve 'count' times
            for _ in range(count):
                new_dist = []
                temp_peaks = []
                for m1, p1 in current_dist:
                    if p1 < 1e-7: continue # Prune
                    for m2, p2 in atom_iso_dist:
                        temp_peaks.append((m1 + m2, p1 * p2))
                
                temp_peaks.sort(key=lambda x: x[0])
                
                # Merge close peaks
                if not temp_peaks: continue
                
                merged = []
                curr_m, curr_p = temp_peaks[0]
                for m, p in temp_peaks[1:]:
                    if m - curr_m < 0.005: 
                        # Weighted average merge
                        total = curr_p + p
                        curr_m = (curr_m * curr_p + m * p) / total
                        curr_p = total
                    else:
                        merged.append((curr_m, curr_p))
                        curr_m, curr_p = m, p
                merged.append((curr_m, curr_p))
                
                # Prune small peaks to keep list manageable
                max_p = max(p for m, p in merged) if merged else 1
                current_dist = [x for x in merged if x[1] > max_p * 0.0001]

        # 4. Final Formatting
        current_dist.sort(key=lambda x: x[0])
        
        if not current_dist:
            return []
            
        max_intensity = max(p for m, p in current_dist)
        final_peaks = [(m, (p / max_intensity) * 100.0) for m, p in current_dist]
        
        return [p for p in final_peaks if p[1] > 0.05]

class HistogramWidget(QWidget):
    def __init__(self, peaks):
        super().__init__()
        self.peaks = peaks # list of (mass, intensity)
        self.setBackgroundRole(QPalette.ColorRole.Base)
        self.setAutoFillBackground(False)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        w = self.width()
        h = self.height()
        
        # 1. Background White
        painter.fillRect(0, 0, w, h, QColor("#ffffff"))
        
        if not self.peaks:
            painter.setPen(QColor("#999999"))
            font = painter.font()
            font.setPointSize(12)
            painter.setFont(font)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No Data")
            return

        # Margins - make plot wider
        ml, mt, mr, mb = 60, 40, 40, 50
        plot_w = w - ml - mr
        plot_h = h - mt - mb
        
        # Calculate Ranges with padding
        masses = [p[0] for p in self.peaks]
        min_mass = min(masses)
        max_mass = max(masses)
        if max_mass == min_mass:
            min_mass -= 5.0
            max_mass += 5.0
        else:
            # Widen the view: Add significant padding
            mass_span = max_mass - min_mass
            padding = max(mass_span * 0.1, 2.0) # Reduced padding to 10%
            min_mass -= padding / 2
            max_mass += padding / 2
            
        mass_range = max_mass - min_mass
        max_intensity = 110.0 

        # 2. Draw Grid
        painter.setPen(QPen(QColor("#eeeeee"), 1, Qt.PenStyle.SolidLine))
        # Horizontal grid lines
        for i in range(1, 6):
            y = (h - mb) - (i * 0.2 * plot_h)
            painter.drawLine(ml, int(y), w - mr, int(y))

        # 3. Draw Axes
        axis_pen = QPen(QColor("#000000"), 2)
        painter.setPen(axis_pen)
        # Y axis
        painter.drawLine(ml, mt, ml, h - mb)
        # X axis
        painter.drawLine(ml, h - mb, w - mr, h - mb)
        
        # 4. Draw Peaks
        font = painter.font()
        font.setPixelSize(11)
        painter.setFont(font)

        for mass, intensity in self.peaks:
            # Map coordinates
            x_ratio = (mass - min_mass) / mass_range if mass_range > 0 else 0.5
            y_ratio = intensity / max_intensity
            
            x_pos = ml + x_ratio * plot_w
            y_pos = (h - mb) - (y_ratio * plot_h)
            base_y = h - mb

            # Draw Stick (Blue)
            stick_pen = QPen(QColor("#007bff"), 3) 
            stick_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(stick_pen)
            painter.drawLine(QPointF(x_pos, base_y), QPointF(x_pos, y_pos))
            
            # Draw Labels
            # Mass label (X axis) - More precision
            painter.setPen(QPen(QColor("#000000")))
            label_rect = QRectF(x_pos - 40, h - mb + 5, 80, 20)
            painter.drawText(label_rect, Qt.AlignmentFlag.AlignCenter, f"{mass:.4f}")
            
            # Intensity label (on top of peak)
            painter.setPen(QPen(QColor("#007bff")))
            int_rect = QRectF(x_pos - 30, y_pos - 20, 60, 15)
            painter.drawText(int_rect, Qt.AlignmentFlag.AlignCenter, f"{int(intensity)}%")

        # Axis Labels
        painter.setPen(QColor("#000000"))
        font.setBold(True)
        painter.setFont(font)
        painter.drawText(QRectF(0, h - 30, w, 20), Qt.AlignmentFlag.AlignCenter, "m/z")
        
        save_state = painter.save()
        painter.translate(20, h / 2)
        painter.rotate(-90)
        painter.drawText(QRectF(-100, -100, 200, 20), Qt.AlignmentFlag.AlignCenter, "Relative Intensity (%)")
        painter.restore()

def run(main_window):
    if not main_window.current_mol:
        QMessageBox.warning(main_window, "MS Plugin", "Please load a molecule first.")
        return

    if Chem is None:
        QMessageBox.critical(main_window, "MS Plugin", "RDKit is not available.")
        return

    dialog = MSSpectrumDialog(main_window.current_mol, main_window)
    dialog.exec()
