# -*- coding: utf-8 -*-
import sys
import math
from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QWidget, QMessageBox, QPushButton, QFileDialog, QCheckBox, QDoubleSpinBox
from PyQt6.QtGui import QPainter, QPen, QBrush, QColor, QFont, QPalette, QLinearGradient, QGradient
from PyQt6.QtCore import Qt, QRectF, QPointF

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    Chem = None

__version__="2025.12.18"
__author__="HiroYokoyama"

PLUGIN_NAME = "MS Spectrum Simulation Neo"

class MSSpectrumDialog(QDialog):
    def __init__(self, mol, parent=None):
        super().__init__(parent)
        self.setWindowTitle("MS Spectrum Simulation Neo")
        self.resize(500, 550) 
        self.mol = mol
        
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
            }
            QGroupBox {
                border: 1px solid #dddddd;
                border-radius: 4px;
                margin-top: 20px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top center;
                padding: 0 5px;
                color: #555555;
            }
        """)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        
        # --- Settings Panel ---
        from PyQt6.QtWidgets import QGroupBox, QFormLayout, QLineEdit, QComboBox, QSpinBox
        
        settings_group = QGroupBox("Configuration")
        settings_layout = QFormLayout(settings_group)
        
        # 1. Formula Input
        self.formula_input = QLineEdit()
        if self.mol:
            self.formula_input.setText(Chem.rdMolDescriptors.CalcMolFormula(self.mol))
        settings_layout.addRow("Formula:", self.formula_input)
        
        # 2. Charge (Signed)
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.setValue(1) # Default +1
        # Prevent 0
        self.charge_spin.valueChanged.connect(self.on_charge_changed)
        settings_layout.addRow("Charge (z):", self.charge_spin)
        
        # 3. Adduct Species
        self.adduct_combo = QComboBox()
        self.update_adduct_options() # Populate initial
        settings_layout.addRow("Adduct:", self.adduct_combo)

        layout.addWidget(settings_group)

        # Gaussian Options
        self.gauss_check = QCheckBox("Gaussian Broadening")
        self.gauss_check.setChecked(False)
        self.gauss_check.stateChanged.connect(self.recalc_peaks)
        
        self.width_spin = QDoubleSpinBox()
        self.width_spin.setRange(0.001, 5.0)
        self.width_spin.setSingleStep(0.01)
        self.width_spin.setValue(0.04)
        self.width_spin.setSuffix(" Da")
        self.width_spin.setToolTip("Peak width (Sigma)")
        self.width_spin.valueChanged.connect(self.recalc_peaks)
        
        # Layout for Gaussian
        gauss_layout = QHBoxLayout()
        gauss_layout.addWidget(self.gauss_check)
        gauss_layout.addWidget(QLabel("Width:"))
        gauss_layout.addWidget(self.width_spin)
        gauss_layout.addStretch()
        layout.addLayout(gauss_layout)
        
        # --- Info Labels ---
        info_layout = QHBoxLayout()
        self.lbl_mw = QLabel("Neutral Avg Mass: -")
        self.lbl_em = QLabel("Neutral Exact Mass: -")
        self.lbl_ion = QLabel("Monoisotopic m/z: -")
        info_layout.addWidget(self.lbl_mw)
        info_layout.addWidget(self.lbl_em)
        info_layout.addWidget(self.lbl_ion)
        layout.addLayout(info_layout)
        
        # --- Plot ---
        plot_container = QWidget()
        plot_container.setStyleSheet("background-color: #ffffff; border: 1px solid #cccccc; border-radius: 4px;")
        plot_layout = QVBoxLayout(plot_container)
        plot_layout.setContentsMargins(1, 1, 1, 1)
        
        self.plot_widget = HistogramWidget([]) # Init empty
        plot_layout.addWidget(self.plot_widget)
        layout.addWidget(plot_container, 1)

        # Export Button
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        self.export_btn = QPushButton("Export to Image")
        self.export_btn.clicked.connect(self.export_image)
        btn_layout.addWidget(self.export_btn)

        self.btn_export_csv = QPushButton("Export CSV")
        self.btn_export_csv.clicked.connect(self.export_csv)
        btn_layout.addWidget(self.btn_export_csv)

        layout.addLayout(btn_layout)

        # Signals
        # Signals
        self.formula_input.textChanged.connect(lambda: self.recalc_peaks(reset=True))
        self.adduct_combo.currentIndexChanged.connect(lambda: self.recalc_peaks(reset=True))
        self.charge_spin.valueChanged.connect(lambda: self.recalc_peaks(reset=True))
        
        # Gaussian changes should NOT reset view (keep zoom)
        self.gauss_check.stateChanged.connect(lambda: self.recalc_peaks(reset=False))
        self.width_spin.valueChanged.connect(lambda: self.recalc_peaks(reset=False))
        
        # Initial Calc
        self.recalc_peaks(reset=True)

    def export_image(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Save Spectrum Image", "spectrum.png", "Images (*.png *.jpg)")
        if filename:
            pixmap = self.plot_widget.grab()
            if pixmap.save(filename):
                QMessageBox.information(self, "Success", f"Image saved to {filename}")
            else:
                QMessageBox.critical(self, "Error", "Failed to save image.")

    def on_charge_changed(self, val):
        # Allow 0
        
        # Always update options to reflect new charge in strings
        # Always update options to reflect new charge in strings
        self.update_adduct_options()
        # Recalc is triggered by signals if values change, but charge spin change
        # manually connected? No, above we connected charge_spin.valueChanged.
        # But wait, update_adduct_options clears combo?
        # Let's rely on signals.
        # self.recalc_peaks() -> Will be called by signal.

    def to_superscript(self, txt):
        mapping = {
            '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
            '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
            '+': '⁺', '-': '⁻'
        }
        return "".join(mapping.get(c, c) for c in str(txt))

    def to_subscript(self, txt):
        mapping = {
            '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
            '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'
        }
        return "".join(mapping.get(c, c) for c in str(txt))

    def update_adduct_options(self):
        # Save current choice
        current_idx = self.adduct_combo.currentIndex()
        if current_idx < 0: current_idx = 0
        
        self.adduct_combo.blockSignals(True)
        self.adduct_combo.clear()
        
        charge_val = self.charge_spin.value()
        z = abs(charge_val)
        
        if charge_val == 0:
            items = [
                "M (Neutral) [M]",
                "None"
            ]
            self.adduct_combo.addItems(items)
            self.adduct_combo.blockSignals(False)
            return

        sign = "⁺" if charge_val > 0 else "⁻"
        
        # Format Charge Superscript (e.g. ²⁺ or ⁺)
        z_super = self.to_superscript(z) if z > 1 else ""
        charge_str = f"{z_super}{sign}"
        
        # Format Count (coefficient) (e.g. 2 or "")
        n_str = str(z) if z > 1 else ""
        
        if charge_val > 0: # Positive
            items = [
                f"M [M]{charge_str}",
                f"H (Proton) [M+{n_str}H]{charge_str}", 
                f"Na (Sodium) [M+{n_str}Na]{charge_str}", 
                f"K (Potassium) [M+{n_str}K]{charge_str}", 
                f"NH4 (Ammonium) [M+{n_str}NH{self.to_subscript(4)}]{charge_str}", 
                f"CH3CN+H (Acetonitrile) [M+{n_str}H+{n_str}CH{self.to_subscript(3)}CN]{charge_str}"
            ]
        else: # Negative
            items = [
                f"M [M]{charge_str}",
                f"H (Deprotonation) [M-{n_str}H]{charge_str}", 
                f"Cl (Chloride) [M+{n_str}Cl]{charge_str}", 
                f"HCOO (Formate) [M+{n_str}HCOO]{charge_str}", 
                f"CH3COO (Acetate) [M+{n_str}CH{self.to_subscript(3)}COO]{charge_str}"
            ]
        
        self.adduct_combo.addItems(items)
        
        # Restore index if safe, else 0
        if current_idx < self.adduct_combo.count():
            self.adduct_combo.setCurrentIndex(current_idx)
        else:
             self.adduct_combo.setCurrentIndex(0)
             
        self.adduct_combo.blockSignals(False)

    def parse_formula_str(self, formula):
        import re
        # Tokenize: Elements (e.g., "Na", "C"), Numbers, Parentheses
        tokens = re.findall(r"([A-Z][a-z]*|\d+|\(|\))", formula)
        
        stack = [{}]
        
        i = 0
        while i < len(tokens):
            token = tokens[i]
            
            if token == '(':
                stack.append({})
                i += 1
            elif token == ')':
                # Check for multiplier after parenthesis
                multiplier = 1
                if i + 1 < len(tokens) and tokens[i+1].isdigit():
                    multiplier = int(tokens[i+1])
                    i += 1 # Consume number
                
                # Pop and merge
                if len(stack) > 1:
                    top = stack.pop()
                    for el, count in top.items():
                        stack[-1][el] = stack[-1].get(el, 0) + count * multiplier
                i += 1
            elif token.isdigit():
                # Should have been handled, but ignore if standalone
                i += 1
            elif token[0].isalpha(): # Element
                element = token
                count = 1
                if i + 1 < len(tokens) and tokens[i+1].isdigit():
                    count = int(tokens[i+1])
                    i += 1 # Consume number
                
                stack[-1][element] = stack[-1].get(element, 0) + count
                i += 1
                
        # Merge any remaining stack items
        while len(stack) > 1:
             top = stack.pop()
             for el, count in top.items():
                 stack[-1][el] = stack[-1].get(el, 0) + count
                 
        return stack[0]

    def get_adduct_delta(self, species_idx, mode, charge):
        # Map index to species logic
        # Positive: 0:M, 1:H, 2:Na, 3:K, 4:NH4, 5:CH3CN+H
        # Negative: 0:M, 1:H(-), 2:Cl, 3:HCOO, 4:CH3COO
        
        delta = {}
        
        # NOTE: Indices shifted because M is now 0
        if mode == "Positive":
            if species_idx == 0: # M
                 pass
            elif species_idx == 1: # H
                 delta = {'H': 1 * charge}
            elif species_idx == 2: # Na
                 delta = {'Na': 1 * charge}
            elif species_idx == 3: # K
                 delta = {'K': 1 * charge}
            elif species_idx == 4: # NH4
                 delta = {'N': 1 * charge, 'H': 4 * charge}
            elif species_idx == 5: # CH3CN+H
                 delta = {'C': 2 * charge, 'H': 4 * charge, 'N': 1 * charge}
            
        else: # Negative
            if species_idx == 0: # M
                 pass
            elif species_idx == 1: # H (Deprotonation -H)
                 delta = {'H': -1 * charge}
            elif species_idx == 2: # Cl
                 delta = {'Cl': 1 * charge}
            elif species_idx == 3: # HCOO
                 delta = {'C': 1 * charge, 'H': 1 * charge, 'O': 2 * charge}
            elif species_idx == 4: # CH3COO
                 delta = {'C': 2 * charge, 'H': 3 * charge, 'O': 2 * charge}

        return delta

    def export_csv(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Save Spectrum CSV", "spectrum.csv", "CSV Files (*.csv)")
        if not filename: return
        
        try:
            with open(filename, 'w') as f:
                f.write("m/z,Intensity\n")
                data = self.plot_widget.peaks
                for m, i in data:
                    f.write(f"{m:.5f},{i:.5f}\n")
            QMessageBox.information(self, "Success", f"Saved to {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save CSV: {e}")

    def recalc_peaks(self, reset=True):
        peaks, final_mw, neutral_exact_mw, ion_mz = self._calculate_peaks()
        
        # Apply Gaussian if checked
        if self.gauss_check.isChecked() and peaks:
            sigma = self.width_spin.value()
            display_peaks = self.apply_gaussian_broadening(peaks, sigma)
            self.plot_widget.draw_mode = "profile"
        else:
            display_peaks = peaks
            self.plot_widget.draw_mode = "stick"

        self.peaks = display_peaks
        self.plot_widget.peaks = display_peaks
        self.plot_widget.stick_peaks = peaks # Save theoretical peaks for labels
        
        if reset:
             self.plot_widget.reset_view()
        else:
             self.plot_widget.update()
             
        # Prepare Info Text for Graph
        formula_raw = self.formula_input.text().strip()
        
        # Format formula (C6H6 -> C₆H₆)
        sub_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        formula_formatted = formula_raw.translate(sub_map)
        
        full_adduct_str = self.adduct_combo.currentText()
        import re
        match = re.search(r"\[.*?\][⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻˙]*", full_adduct_str)
        if match:
             adduct_display = match.group(0)
        else:
             adduct_display = full_adduct_str
        
        adduct_display = adduct_display.replace("M", formula_formatted)

        if self.charge_spin.value() == 0:
             adduct_display = "Neutral"
             self.plot_widget.info_text = f"{formula_formatted}\n{adduct_display}"
        else:
             self.plot_widget.info_text = f"{adduct_display}"

        self.plot_widget.update()
        
        # Update Info
        try:
            self.lbl_ion.setText(f"<b>Monoisotopic m/z:</b> {ion_mz:.4f}")
            self.lbl_mw.setText(f"<b>Neutral Avg Mass:</b> {final_mw:.4f}")
            self.lbl_em.setText(f"<b>Neutral Exact Mass:</b> {neutral_exact_mw:.4f}")
        except:
            pass

    def apply_gaussian_broadening(self, peaks, sigma):
        import math
        if not peaks: return []
        
        # Dynamic margin to avoid clipping wide peaks
        margin = max(2.0, 6.0 * sigma)
        min_mz = min(p[0] for p in peaks) - margin
        max_mz = max(p[0] for p in peaks) + margin
        
        # Grid - Finer resolution to minimize peak shift artifacts
        import numpy as np
        # 10 points per sigma is better for appearance, 20 is good for peak top approx
        step = min(0.005, sigma / 10.0) 
        x_vals = np.arange(min_mz, max_mz, step)
        y_vals = np.zeros_like(x_vals)
        
        # Sum Gaussians
        for center, intensity in peaks:
            indices = np.where((x_vals >= center - 5*sigma) & (x_vals <= center + 5*sigma))
            if len(indices[0]) == 0: continue
            
            xs = x_vals[indices]
            ys = intensity * np.exp( - (xs - center)**2 / (2 * sigma**2) )
            y_vals[indices] += ys
            
        if np.max(y_vals) > 0:
            y_vals = (y_vals / np.max(y_vals)) * 100.0
            
        return list(zip(x_vals, y_vals))

    def _calculate_peaks(self):
        formula_str = self.formula_input.text().strip()
        if not formula_str:
            return [], 0.0, 0.0, 0.0

        # 1. Base composition
        base_counts = self.parse_formula_str(formula_str)
        if not base_counts:
             return [], 0.0, 0.0, 0.0

        # 2. Adduct adjustments
        charge_val = self.charge_spin.value() # Signed
        charge = abs(charge_val)
        mode = "Positive" if charge_val > 0 else "Negative"
        species_idx = self.adduct_combo.currentIndex()
        
        delta = self.get_adduct_delta(species_idx, mode, charge)
        
        # Merge counts
        final_counts = base_counts.copy()
        for el, count in delta.items():
            final_counts[el] = final_counts.get(el, 0) + count
            
        # Remove elements with <= 0 count
        final_counts = {k: v for k, v in final_counts.items() if v > 0}
        
        if not final_counts:
            return [], 0.0, 0.0, 0.0

        # 3. Dynamic Isotope Calc
        pt = Chem.GetPeriodicTable()
        current_dist = [(0.0, 1.0)]
        electron_mass = 0.00054858
        
        total_mw = 0.0 # Neutral average mass

        for sym, count in final_counts.items():
            # Special handling for Deuterium
            if sym == "D":
                 mass_d = 2.0141017781
                 total_mw += mass_d * count
                 atom_iso_dist = [(mass_d, 1.0)]
            else:
                try:
                    base_mass = pt.GetAtomicWeight(sym) # Average weight
                    total_mw += base_mass * count
                    
                    atomic_num = pt.GetAtomicNumber(sym)
                    atom_iso_dist = []
                    center_mass = pt.GetMostCommonIsotope(atomic_num)
                    # Scan a range
                    for m in range(max(1, center_mass - 5), center_mass + 10):
                        try:
                            abundance = pt.GetAbundanceForIsotope(atomic_num, m)
                            if abundance > 0.00001:
                                exact_mass = pt.GetMassForIsotope(atomic_num, m)
                                atom_iso_dist.append((exact_mass, abundance))
                        except RuntimeError: pass
                    
                    # Normalize
                    total_p = sum(p for m, p in atom_iso_dist)
                    if total_p == 0: continue
                    atom_iso_dist = [(m, p / total_p) for m, p in atom_iso_dist]
                except:
                    # Ignore unknown elements
                    continue

            # Convolve 'count' times
            for _ in range(count):
                new_peaks = {}
                for m1, p1 in current_dist:
                    if p1 < 1e-5: continue
                    for m2, p2 in atom_iso_dist:
                        m = m1 + m2
                        p = p1 * p2
                        m_bin = round(m, 4) 
                        new_peaks[m_bin] = new_peaks.get(m_bin, 0) + p
                
                sorted_p = sorted(new_peaks.items())
                
                merged = []
                if sorted_p:
                    curr_m, curr_p = sorted_p[0]
                    for m, p in sorted_p[1:]:
                        if m - curr_m < 0.005:
                            total = curr_p + p
                            curr_m = (curr_m * curr_p + m * p) / total
                            curr_p = total
                        else:
                            merged.append((curr_m, curr_p))
                            curr_m, curr_p = m, p
                    merged.append((curr_m, curr_p))
                
                max_p = max(p for m, p in merged) if merged else 1
                current_dist = [x for x in merged if x[1] > max_p * 0.001]
                if len(current_dist) > 200:
                        current_dist.sort(key=lambda x: x[1], reverse=True)
                        current_dist = current_dist[:200]

        # Calculate EXACT NEUTRAL MASS (Always needed)
        exact_mass_sum = 0.0
        for sym, count in final_counts.items():
            if sym == "D":
                 exact_mass_sum += 2.0141017781 * count
            else:
                 try:
                     anum = pt.GetAtomicNumber(sym)
                     m_iso = pt.GetMostCommonIsotope(anum)
                     mass_iso = pt.GetMassForIsotope(anum, m_iso)
                     exact_mass_sum += mass_iso * count
                 except: pass

        neutral_exact_mass = exact_mass_sum

        # 4. Adjust for Charge (m/z)
        if charge_val == 0:
            # Neutral Mode: Hide Spectrum
            # Return empty peaks, but valid neutral masses
            return [], total_mw, neutral_exact_mass, 0.0

        # Apply electron mass correction
        total_ion_mass_dist = []
        
        e_params = 0
        if mode == "Positive":
             e_params = -1 * charge * electron_mass
        else:
             e_params = +1 * charge * electron_mass
             
        for m, p in current_dist:
            ion_mass = m + e_params
            if ion_mass <= 0: continue
            mz = ion_mass / charge
            total_ion_mass_dist.append((mz, p))
            
        if not total_ion_mass_dist:
            return [], total_mw, neutral_exact_mass, 0.0

        max_intensity = max(p for m, p in total_ion_mass_dist)
        final_peaks = [(m, (p / max_intensity) * 100.0) for m, p in total_ion_mass_dist]
        final_peaks.sort(key=lambda x: x[0])
        
        if mode == "Positive":
             exact_mass_sum += (-1 * charge * electron_mass)
        else:
             exact_mass_sum += (+1 * charge * electron_mass)
             
        exact_mz = exact_mass_sum / charge if charge != 0 else 0

        return [p for p in final_peaks if p[1] > 0.05], total_mw, neutral_exact_mass, exact_mz


class HistogramWidget(QWidget):
    def __init__(self, peaks):
        super().__init__()
        self.peaks = peaks # list of (mass, intensity) - used for drawing (stick or profile)
        self.stick_peaks = [] # Always holds the theoretical stick peaks for labeling
        self.info_text = ""
        self.draw_mode = "stick" # "stick" or "profile"
        self.setBackgroundRole(QPalette.ColorRole.Base)
        self.setAutoFillBackground(False)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent)
        
        # View State
        self.view_min = None
        self.view_max = None
        self.last_mouse_x = None

    def reset_view(self):
        # Auto-scale based on stick peaks (stable) or actual peaks
        if self.stick_peaks:
            masses = [p[0] for p in self.stick_peaks]
        elif self.peaks:
            masses = [p[0] for p in self.peaks]
        else:
            masses = []

        if masses:
            min_m, max_m = min(masses), max(masses)
            if max_m == min_m:
                 min_m -= 5.0
                 max_m += 5.0
            else:
                 padding = 5.0
                 min_m -= padding
                 max_m += padding
            self.view_min = min_m
            self.view_max = max_m
        else:
            self.view_min = 0.0
            self.view_max = 100.0
        
        self.update()

    def wheelEvent(self, event):
        if self.view_min is None or self.view_max is None: return
        
        # Current Range
        current_range = self.view_max - self.view_min
        if current_range <= 0: return

        # Zoom Factor
        angle = event.angleDelta().y()
        factor = 0.9 if angle > 0 else 1.1
        
        # Mouse Position Ratio (0.0 to 1.0) relative to plot area
        w = self.width()
        ml, mr = 60, 40
        plot_w = w - ml - mr
        if plot_w <= 0: return
        
        mouse_x = event.position().x()
        ratio = (mouse_x - ml) / plot_w
        # Clamp ratio
        ratio = max(0.0, min(1.0, ratio))
        
        # Calculate new range
        new_range = current_range * factor
        
        # Limit zoom (e.g. max range 2000, min range 1.0)
        if new_range > 5000: new_range = 5000
        if new_range < 1.0: new_range = 1.0
        
        # Adjust min/max keeping ratio point fixed
        change = new_range - current_range
        
        self.view_min -= change * ratio
        self.view_max += change * (1.0 - ratio)
        
        self.update()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.last_mouse_x = event.position().x()

    def mouseMoveEvent(self, event):
        if self.last_mouse_x is not None and self.view_min is not None:
             dx = event.position().x() - self.last_mouse_x
             
             w = self.width()
             ml, mr = 60, 40
             plot_w = w - ml - mr
             
             if plot_w > 0:
                 current_range = self.view_max - self.view_min
                 # Pixels to Mass Units
                 mass_shift = (dx / plot_w) * current_range
                 
                 self.view_min -= mass_shift
                 self.view_max -= mass_shift
                 self.update()
             
             self.last_mouse_x = event.position().x()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.last_mouse_x = None

    def mouseDoubleClickEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.reset_view()



    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        w = self.width()
        h = self.height()
        
        # 1. Background White
        painter.fillRect(0, 0, w, h, QColor("#ffffff"))
        
        # Margins (Increased Top Margin for Labels)
        ml, mt, mr, mb = 60, 60, 40, 50


        
        if not self.peaks:
            painter.setPen(QColor("#999999"))
            font = painter.font()
            font.setPointSize(12)
            font.setBold(False)
            painter.setFont(font)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No Data")
            return

        plot_w = w - ml - mr
        plot_h = h - mt - mb
        
        if self.view_min is None or self.view_max is None:
             self.reset_view()

        min_mass = self.view_min
        max_mass = self.view_max
            
        mass_range = max_mass - min_mass
        max_intensity = 110.0 

        # 2. Draw Grid & Y-Axis Labels
        painter.setPen(QPen(QColor("#eeeeee"), 1, Qt.PenStyle.SolidLine))
        font = painter.font()
        font.setPixelSize(10)
        font.setBold(False)
        painter.setFont(font)
        
        # Grid lines and Y Labels
        for i in range(0, 6): # 0 to 5 (0, 20, 40, 60, 80, 100)
            val = i * 20
            ratio = val / max_intensity
            y = (h - mb) - (ratio * plot_h)
            
            # Grid Line
            if i > 0: # Don't draw over X axis
                painter.setPen(QPen(QColor("#eeeeee"), 1, Qt.PenStyle.SolidLine))
                painter.drawLine(ml, int(y), w - mr, int(y))
            
            # Label
            painter.setPen(QColor("#000000"))
            painter.drawText(QRectF(0, y - 10, ml - 5, 20), Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, str(val))

        # 3. Draw Axes
        axis_pen = QPen(QColor("#000000"), 2)
        painter.setPen(axis_pen)
        # Y axis
        painter.drawLine(ml, mt, ml, h - mb)
        # X axis
        painter.drawLine(ml, h - mb, w - mr, h - mb)
        
        # X-Axis Ticks
        n_ticks = 5
        tick_step = mass_range / n_ticks
        import math
        mag = 10**math.floor(math.log10(tick_step)) if tick_step > 0 else 1
        nice_step = round(tick_step / mag) * mag
        if nice_step == 0: nice_step = tick_step
        
        start_tick = math.ceil(min_mass / nice_step) * nice_step
        
        painter.setPen(QColor("#000000"))
        curr_tick = start_tick
        while curr_tick <= max_mass:
            x_ratio = (curr_tick - min_mass) / mass_range if mass_range > 0 else 0.5
            x = ml + x_ratio * plot_w
            
            if ml <= x <= w - mr:
                painter.drawLine(int(x), h - mb, int(x), h - mb + 5)
                painter.drawText(QRectF(x - 30, h - mb + 5, 60, 20), Qt.AlignmentFlag.AlignCenter, f"{curr_tick:.1f}")
            
            curr_tick += nice_step
        
        # 4. Draw Peaks


        font.setPixelSize(11)
        painter.setFont(font)
        
        base_y = h - mb

        if self.draw_mode == "profile":
            # 1. Draw Profile Curve (CLIPPED)
            painter.save()
            painter.setClipRect(ml, mt, plot_w, plot_h)
            
            path_points = []
            first_mass = self.peaks[0][0]
            first_x_ratio = (first_mass - min_mass) / mass_range if mass_range > 0 else 0.5
            first_x = ml + first_x_ratio * plot_w
            path_points.append(QPointF(first_x, base_y))
            
            # Use self.peaks (profile data) for Curve
            for mass, intensity in self.peaks:
                x_ratio = (mass - min_mass) / mass_range if mass_range > 0 else 0.5
                y_ratio = intensity / max_intensity
                
                x_pos = ml + x_ratio * plot_w
                y_pos = (h - mb) - (y_ratio * plot_h)
                path_points.append(QPointF(x_pos, y_pos))
            
            last_mass = self.peaks[-1][0]
            last_x_ratio = (last_mass - min_mass) / mass_range if mass_range > 0 else 0.5
            last_x = ml + last_x_ratio * plot_w
            path_points.append(QPointF(last_x, base_y))
            
            profile_pen = QPen(QColor("#007bff"), 2)
            painter.setPen(profile_pen)
            painter.drawPolyline(path_points)
            
            painter.restore() # Unclip for labels

            # 2. Peak Labels (Theoretical Stick Peaks) - UNCLIPPED
            painter.setPen(QPen(QColor("#007bff")))

            if self.stick_peaks:
                for s_mass, s_int in self.stick_peaks:
                    # Filter very low intensity theoretical peaks
                    if s_int < 0.05: continue 
                    
                    if s_mass < min_mass or s_mass > max_mass: continue

                    x_ratio = (s_mass - min_mass) / mass_range if mass_range > 0 else 0.5
                    y_ratio = s_int / max_intensity
                    
                    x_pos = ml + x_ratio * plot_w
                    y_pos = (h - mb) - (y_ratio * plot_h)
                    
                    if not (ml <= x_pos <= w - mr): continue
                    
                    # Mass Label - Higher Position over stick peak
                    # NOTE: We draw label at STICK height, not PROFILE height
                    label_rect = QRectF(x_pos - 40, y_pos - 25, 80, 20)
                    painter.setPen(QColor("#000000"))
                    painter.drawText(label_rect, Qt.AlignmentFlag.AlignCenter, f"{s_mass:.4f}")

        else:
            # STICK MODE
            
            # 1. Draw Sticks (CLIPPED)
            painter.save()
            painter.setClipRect(ml, mt, plot_w, plot_h)
            
            stick_pen = QPen(QColor("#007bff"), 3) 
            stick_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(stick_pen)
            
            for mass, intensity in self.peaks:
                x_ratio = (mass - min_mass) / mass_range if mass_range > 0 else 0.5
                y_ratio = intensity / max_intensity
                
                x_pos = ml + x_ratio * plot_w
                y_pos = (h - mb) - (y_ratio * plot_h)

                painter.drawLine(QPointF(x_pos, base_y), QPointF(x_pos, y_pos))
                
            painter.restore() # Unclip for labels

            # 2. Draw Labels (UNCLIPPED)
            for mass, intensity in self.peaks:
                x_ratio = (mass - min_mass) / mass_range if mass_range > 0 else 0.5
                y_ratio = intensity / max_intensity
                
                x_pos = ml + x_ratio * plot_w
                y_pos = (h - mb) - (y_ratio * plot_h)
                
                if x_pos < ml or x_pos > w - mr: continue # Simple visibility check

                painter.setPen(QPen(QColor("#000000")))
                
                # Mass Label (Line 1)
                label_rect = QRectF(x_pos - 40, y_pos - 35, 80, 15) 
                painter.drawText(label_rect, Qt.AlignmentFlag.AlignCenter, f"{mass:.4f}")
                
                # Intensity Label (Line 2)
                painter.setPen(QPen(QColor("#007bff")))
                int_rect = QRectF(x_pos - 40, y_pos - 20, 80, 15)
                painter.drawText(int_rect, Qt.AlignmentFlag.AlignCenter, f"{int(intensity)}%")
        
        # Axis Labels
        painter.setPen(QColor("#000000"))
        font.setBold(True)
        painter.setFont(font)
        painter.drawText(QRectF(0, h - 30, w, 20), Qt.AlignmentFlag.AlignCenter, "m/z")

        # Draw Info Text (Top Right) - Drawn Last to be on top of Grid
        if self.info_text:
             painter.setPen(QColor("#333333"))
             font = painter.font()
             font.setPointSize(14)
             font.setBold(True)
             painter.setFont(font)
             
             rect = QRectF(w - 250 - mr, mt, 250, 60)
             # Optional: White background for text to strictly hide grid? 
             # User just said "change drawing order", implying text on top is enough.
             # If grid is black and text is black, it might collide.
             # But usually text on top is standard.
             # Let's add a slight semi-transparent white bg if needed? 
             # No, user said "Revert", so let's go back to original simple text but drawn last.
             painter.drawText(rect, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignTop, self.info_text)
        
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
