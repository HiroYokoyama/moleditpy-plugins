# -*- coding: utf-8 -*-
import sys
import re
import logging
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QFileDialog,
    QCheckBox,
    QDoubleSpinBox,
    QComboBox,
    QGroupBox,
    QFormLayout,
    QLineEdit,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QAbstractItemView,
    QApplication,
)
from PyQt6.QtGui import (
    QPainter,
    QKeySequence,
    QColor,
    QFont,
    QPageSize,
    QTextDocument,
    QImage,
    QPageLayout,
)
from PyQt6.QtCore import (
    Qt,
    QRectF,
    QTimer,
    QByteArray,
    QBuffer,
    QIODevice,
    QSizeF,
    QMarginsF,
)
from PyQt6.QtPrintSupport import QPrinter

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, rdMolDescriptors
except ImportError:
    Chem = None
    Draw = None
    rdMolDescriptors = None

PLUGIN_ID = "elemental_analysis"
PLUGIN_VERSION = "2026.09.26"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_SUPPORTED_PYTHON_VERSION = ">=3.9, <3.15"
PLUGIN_AUTHOR = "HiroYokoyama"

PLUGIN_NAME = "Elemental Analysis"
PLUGIN_DESCRIPTION = "Calculates the theoretical elemental composition (mass %) of the currently loaded molecule, with optional solvate, and formats the 'Anal. Calcd for' line for publications."
PLUGIN_CATEGORY = "Analysis"
PLUGIN_TAGS = ["Analysis"]
PLUGIN_DEPENDENCIES = ["rdkit", "PyQt6"]

# RDKit's periodic table has no entry for the hydrogen isotope symbols.
HYDROGEN_ISOTOPE_MASSES = {"D": 2.0141017781, "T": 3.0160492779}

#: (label, formula) of the solvents most often co-crystallised or retained
#: in a sample submitted for combustion analysis.
SOLVATES = [
    ("H2O (Water)", "H2O"),
    ("CH2Cl2 (Dichloromethane)", "CH2Cl2"),
    ("CHCl3 (Chloroform)", "CHCl3"),
    ("CH3OH (Methanol)", "CH4O"),
    ("C2H5OH (Ethanol)", "C2H6O"),
    ("CH3CN (Acetonitrile)", "C2H3N"),
    ("Et2O (Diethyl ether)", "C4H10O"),
    ("EtOAc (Ethyl acetate)", "C4H8O2"),
    ("Acetone", "C3H6O"),
    ("THF (Tetrahydrofuran)", "C4H8O"),
    ("Hexane", "C6H14"),
    ("Pentane", "C5H12"),
    ("Benzene", "C6H6"),
    ("Toluene", "C7H8"),
    ("DMSO", "C2H6OS"),
    ("DMF", "C3H7NO"),
    ("HCl", "HCl"),
]

#: Elements quoted in an "Anal. Calcd" line when present (combustion analysis).
REPORT_ELEMENTS = ("C", "H", "N", "S")


def parse_formula(formula):
    """Element counts of a molecular formula, or None if it is not one.

    Handles parentheses with multipliers, ignores charge signs, and accepts
    solvate notation joined by '·' or '*' with an optional (fractional)
    leading coefficient, e.g. "C6H6·0.5H2O". Counts may therefore be floats.
    """
    formula = formula.replace(" ", "")
    if not formula:
        return {}

    total = {}
    for part in re.split(r"[·•*]", formula):
        if not part:
            return None
        m = re.match(r"^(\d+(?:\.\d+)?|\.\d+)?(.*)$", part)
        coeff = float(m.group(1)) if m.group(1) else 1
        counts = _parse_simple_formula(m.group(2))
        if not counts:
            return None
        for el, count in counts.items():
            total[el] = total.get(el, 0) + count * coeff

    return {el: _tidy_number(c) for el, c in total.items() if c > 0}


def _parse_simple_formula(formula):
    tokens = re.findall(r"([A-Z][a-z]*|\d+|\(|\)|[\+\-])", formula)
    if "".join(tokens) != formula:
        return None  # Invalid characters found

    stack = [{}]
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token == "(":
            stack.append({})
            i += 1
        elif token == ")":
            multiplier = 1
            if i + 1 < len(tokens) and tokens[i + 1].isdigit():
                multiplier = int(tokens[i + 1])
                i += 1
            if len(stack) > 1:
                top = stack.pop()
                for el, count in top.items():
                    stack[-1][el] = stack[-1].get(el, 0) + count * multiplier
            i += 1
        elif token[0].isalpha():
            count = 1
            if i + 1 < len(tokens) and tokens[i + 1].isdigit():
                count = int(tokens[i + 1])
                i += 1
            stack[-1][token] = stack[-1].get(token, 0) + count
            i += 1
        else:
            # Charge signs and stray numbers carry no composition.
            i += 1

    while len(stack) > 1:
        top = stack.pop()
        for el, count in top.items():
            stack[-1][el] = stack[-1].get(el, 0) + count

    return stack[0]


def _tidy_number(value):
    """Return an int when the value is integral, else the float."""
    if abs(value - round(value)) < 1e-9:
        return int(round(value))
    return value


def format_count(value):
    value = _tidy_number(value)
    if isinstance(value, int):
        return str(value)
    return f"{value:.3f}".rstrip("0").rstrip(".")


def hill_order(symbols):
    """Hill system: C, H first then alphabetical; purely alphabetical without C."""
    symbols = list(symbols)
    if "C" in symbols:
        head = ["C"] + (["H"] if "H" in symbols else [])
        return head + sorted(s for s in symbols if s not in ("C", "H"))
    return sorted(symbols)


def format_formula(counts):
    """Hill-ordered formula string; a count of 1 is omitted."""
    parts = []
    for sym in hill_order(counts):
        n = counts[sym]
        parts.append(sym if n == 1 else f"{sym}{format_count(n)}")
    return "".join(parts)


def atomic_weight(pt, sym):
    if sym in HYDROGEN_ISOTOPE_MASSES:
        return HYDROGEN_ISOTOPE_MASSES[sym]
    return pt.GetAtomicWeight(sym)


def monoisotopic_mass(pt, counts):
    """Mass built from each element's most abundant isotope."""
    total = 0.0
    for sym, count in counts.items():
        if sym in HYDROGEN_ISOTOPE_MASSES:
            total += HYDROGEN_ISOTOPE_MASSES[sym] * count
        else:
            anum = pt.GetAtomicNumber(sym)
            total += pt.GetMassForIsotope(anum, pt.GetMostCommonIsotope(anum)) * count
    return total


def calc_composition(pt, counts):
    """(rows, molecular_weight) for a composition.

    rows are Hill-ordered (symbol, count, mass contribution, mass %).
    Raises for an element symbol RDKit does not know.
    """
    contrib = {sym: atomic_weight(pt, sym) * n for sym, n in counts.items()}
    mw = sum(contrib.values())
    if mw <= 0:
        return [], 0.0
    rows = [
        (sym, counts[sym], contrib[sym], contrib[sym] / mw * 100.0)
        for sym in hill_order(counts)
    ]
    return rows, mw


def analysis_line(formula_label, rows, elements=REPORT_ELEMENTS):
    """Publication-style line, e.g. 'Anal. Calcd for C6H6: C, 92.26; H, 7.74.'"""
    by_sym = {r[0]: r[3] for r in rows}
    items = [f"{sym}, {by_sym[sym]:.2f}" for sym in elements if sym in by_sym]
    if not items:
        return ""
    return f"Anal. Calcd for {formula_label}: " + "; ".join(items) + "."


def rows_to_tsv(rows, header=True):
    """Tab-separated composition table; pastes as cells into Excel/Word."""
    lines = ["Element\tCount\tMass %"] if header else []
    lines += [f"{sym}\t{format_count(count)}\t{pct:.2f}" for sym, count, _c, pct in rows]
    return "\n".join(lines)


class ElementalAnalysisDialog(QDialog):
    def __init__(self, context):
        super().__init__(parent=context.get_main_window())
        self.setWindowTitle("Elemental Analysis")
        self.resize(500, 560)
        self.context = context
        self.mol = self.context.current_molecule
        self.rows = []
        self.total_formula = ""

        # Setup timer for auto-update
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.check_update)

        # Clean white look - forced light mode to prevent dark mode unreadability
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
                color: #000000;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top center;
                padding: 0 5px;
                color: #555555;
            }
            QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
                background-color: #ffffff;
                color: #000000;
            }
            QComboBox QAbstractItemView {
                background-color: #ffffff;
                color: #000000;
            }
            QTableWidget {
                background-color: #ffffff;
                color: #000000;
                gridline-color: #dddddd;
            }
            QHeaderView::section {
                background-color: #f2f2f2;
                color: #000000;
                border: 1px solid #dddddd;
                padding: 3px;
            }
            QCheckBox {
                color: #000000;
            }
            QPushButton {
                background-color: #f0f0f0;
                color: #000000;
                border: 1px solid #cccccc;
                padding: 5px;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #e0e0e0;
            }
        """)

        # Register window for V3 lifecycle management
        self.context.register_window("main_panel", self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)

        # --- Settings Panel ---
        settings_group = QGroupBox("Configuration")
        settings_layout = QFormLayout(settings_group)

        self.formula_input = QLineEdit()
        self.formula_input.setToolTip(
            "Molecular formula. Solvates may be appended, e.g. C6H6·0.5H2O"
        )
        if self.mol and Chem is not None:
            try:
                self.formula_input.setText(
                    str(rdMolDescriptors.CalcMolFormula(self.mol))
                )
            except (RuntimeError, AttributeError, TypeError, ValueError) as _e:
                logging.warning("[elemental_analysis.py] silenced: %s", _e)
        settings_layout.addRow("Formula:", self.formula_input)

        self.sync_check = QCheckBox("Sync with Main Window")
        self.sync_check.stateChanged.connect(self.toggle_sync)

        self.use_2d_check = QCheckBox("Use 2D molecule")
        self.use_2d_check.setChecked(True)
        self.use_2d_check.toggled.connect(lambda: self.check_update())

        sync_layout = QHBoxLayout()
        sync_layout.addWidget(self.sync_check)
        sync_layout.addWidget(self.use_2d_check)
        sync_layout.addStretch()

        if self.context.get_main_window():
            self.sync_check.setChecked(True)
        else:
            self.sync_check.setChecked(False)
            self.sync_check.setEnabled(False)  # Disable in standalone/mock
            self.use_2d_check.setEnabled(False)
            self.use_2d_check.setChecked(False)

        settings_layout.addRow("Options:", sync_layout)

        # Solvate: n x solvent
        self.solvate_spin = QDoubleSpinBox()
        self.solvate_spin.setDecimals(2)
        self.solvate_spin.setRange(0.0, 20.0)
        self.solvate_spin.setSingleStep(0.25)
        self.solvate_spin.setValue(0.0)
        self.solvate_spin.setToolTip("Equivalents of solvent per formula unit")

        self.solvent_combo = QComboBox()
        for label, formula in SOLVATES:
            self.solvent_combo.addItem(label, formula)

        solvate_layout = QHBoxLayout()
        solvate_layout.addWidget(self.solvate_spin)
        solvate_layout.addWidget(QLabel("×"))
        solvate_layout.addWidget(self.solvent_combo, 1)
        settings_layout.addRow("Solvate:", solvate_layout)

        layout.addWidget(settings_group)

        # --- Info Labels ---
        info_group = QGroupBox("Mass Information")
        info_layout = QVBoxLayout(info_group)
        self.lbl_formula = QLabel("Formula: -")
        self.lbl_mw = QLabel("Molecular Weight: -")
        self.lbl_em = QLabel("Exact Mass: -")
        for lbl in (self.lbl_formula, self.lbl_mw, self.lbl_em):
            lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        info_layout.addWidget(self.lbl_formula)
        info_layout.addWidget(self.lbl_mw)
        info_layout.addWidget(self.lbl_em)
        layout.addWidget(info_group)

        # --- Composition Group ---
        comp_group = QGroupBox("Composition")
        comp_layout = QVBoxLayout(comp_group)

        self.table = CopyableTable(0, 3)
        self.table.setHorizontalHeaderLabels(["Element", "Count", "Mass %"])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        comp_layout.addWidget(self.table)

        # Publication line + copy
        line_layout = QHBoxLayout()
        self.lbl_report = QLabel("")
        self.lbl_report.setWordWrap(True)
        self.lbl_report.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
        )
        line_layout.addWidget(self.lbl_report, 1)
        self.btn_copy = QPushButton("Copy")
        self.btn_copy.setToolTip("Copy the 'Anal. Calcd for' line")
        self.btn_copy.clicked.connect(self.copy_report_line)
        line_layout.addWidget(self.btn_copy)
        self.btn_copy_table = QPushButton("Copy Table")
        self.btn_copy_table.setToolTip("Copy the composition table (tab-separated)")
        self.btn_copy_table.clicked.connect(self.copy_table)
        line_layout.addWidget(self.btn_copy_table)
        comp_layout.addLayout(line_layout)

        layout.addWidget(comp_group)
        layout.addStretch(1)

        # Export Buttons
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        self.btn_export_csv = QPushButton("Export CSV")
        self.btn_export_csv.clicked.connect(self.export_csv)
        btn_layout.addWidget(self.btn_export_csv)

        self.btn_report = QPushButton("PDF Report")
        self.btn_report.clicked.connect(self.create_report)
        btn_layout.addWidget(self.btn_report)

        layout.addLayout(btn_layout)

        # Signals
        self.formula_input.textChanged.connect(lambda: self.recalc())
        self.formula_input.textEdited.connect(lambda: self.sync_check.setChecked(False))
        self.solvate_spin.valueChanged.connect(lambda: self.recalc())
        self.solvent_combo.currentIndexChanged.connect(lambda: self.recalc())

        # Initial Sync & Calc
        self.check_update()
        self.recalc()

    def _stop_sync_timer(self):
        if getattr(self, "timer", None) is not None and self.timer.isActive():
            self.timer.stop()

    def done(self, result: int):
        # Esc reaches reject() -> done() without a QCloseEvent; stop polling here too.
        self._stop_sync_timer()
        super().done(result)

    def closeEvent(self, event):
        self._stop_sync_timer()
        event.accept()

    def toggle_sync(self, state):
        if getattr(self, "timer", None) is None:
            return

        if state == 2:  # Checked
            self.timer.start(500)
            self.check_update()
        else:
            self.timer.stop()

    def check_update(self):
        if Chem is None:
            return

        try:
            if self.use_2d_check.isChecked():
                mw = self.context.get_main_window()
                if hasattr(mw, "state_manager") and hasattr(
                    mw.state_manager.data, "to_rdkit_mol"
                ):
                    new_mol = mw.state_manager.data.to_rdkit_mol()
                else:
                    new_mol = self.context.current_molecule
            else:
                new_mol = self.context.current_molecule

            if not new_mol:
                return

            self.mol = new_mol

            mol_to_calc = (
                Chem.AddHs(self.mol) if self.use_2d_check.isChecked() else self.mol
            )
            current_formula = rdMolDescriptors.CalcMolFormula(mol_to_calc)

            if self.formula_input.text() != current_formula:
                self.formula_input.setText(current_formula)
        except Exception as _e:
            # Fail silently to avoid spamming errors in timer
            logging.warning("[elemental_analysis.py] silenced: %s", _e)

    def solvate_suffix(self):
        """'·0.5H2O'-style suffix for the current solvate, or '' for none."""
        n = self.solvate_spin.value()
        if n <= 0:
            return ""
        coeff = "" if abs(n - 1.0) < 1e-9 else format_count(n)
        return f"·{coeff}{self.solvent_combo.currentData()}"

    def recalc(self):
        base_text = self.formula_input.text().strip()
        base = parse_formula(base_text) if base_text else {}
        self.rows = []
        self.total_formula = ""
        mw = exact = 0.0

        if base and Chem is not None:
            suffix = self.solvate_suffix()
            counts = parse_formula(format_formula(base) + suffix)
            try:
                pt = Chem.GetPeriodicTable()
                self.rows, mw = calc_composition(pt, counts)
                exact = monoisotopic_mass(pt, base)
                self.total_formula = format_formula(base) + suffix
            except Exception:
                # Unrecognized element -> invalid formula
                self.rows, mw, exact = [], 0.0, 0.0

        self._fill_table()

        if self.rows:
            self.lbl_formula.setText(f"<b>Formula:</b> {self.total_formula}")
            self.lbl_mw.setText(f"<b>Molecular Weight:</b> {mw:.4f}")
            self.lbl_em.setText(f"<b>Exact Mass (M):</b> {exact:.4f}")
            self.lbl_report.setText(analysis_line(self.total_formula, self.rows))
        else:
            self.lbl_formula.setText("<b>Formula:</b> -")
            self.lbl_mw.setText("<b>Molecular Weight:</b> -")
            self.lbl_em.setText("<b>Exact Mass (M):</b> -")
            self.lbl_report.setText("")

    def _fill_table(self):
        self.table.setRowCount(len(self.rows))
        for r, (sym, count, _contrib, pct) in enumerate(self.rows):
            for c, text in enumerate((sym, format_count(count), f"{pct:.2f}")):
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.table.setItem(r, c, item)
        # Fit the rows (up to 7) so a typical formula needs no scrolling.
        visible = max(1, min(len(self.rows), 7))
        self.table.setFixedHeight(
            self.table.horizontalHeader().height()
            + visible * self.table.verticalHeader().defaultSectionSize()
            + 2 * self.table.frameWidth()
        )

    def copy_report_line(self):
        text = self.lbl_report.text()
        if text:
            QApplication.clipboard().setText(text)

    def copy_table(self):
        if self.rows:
            QApplication.clipboard().setText(rows_to_tsv(self.rows))

    def export_csv(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Composition CSV", "composition.csv", "CSV Files (*.csv)"
        )
        if not filename:
            return

        try:
            with open(filename, "w", encoding="utf-8") as f:
                f.write(f"# Formula,{self.total_formula}\n")
                f.write("Element,Count,Mass contribution,Mass %\n")
                for sym, count, contrib, pct in self.rows:
                    f.write(f"{sym},{format_count(count)},{contrib:.5f},{pct:.4f}\n")
            QMessageBox.information(self, "Success", f"Saved to {filename}")
        except (OSError, RuntimeError, AttributeError, TypeError, ValueError) as e:
            QMessageBox.critical(self, "Error", f"Failed to save CSV: {e}")

    def _molecule_image_b64(self):
        """(base64 PNG, width, height) of the structure, or ('', 0, 0)."""
        if not self.mol:
            return "", 0, 0

        # 1. The user's own 2D view
        scene = self.context.scene
        if scene:
            try:
                molecule_bounds = QRectF()
                for item in scene.items():
                    if type(item).__name__ in ("AtomItem", "BondItem") and item.isVisible():
                        molecule_bounds = molecule_bounds.united(item.sceneBoundingRect())
                if molecule_bounds.isEmpty():
                    molecule_bounds = scene.itemsBoundingRect()

                if not molecule_bounds.isEmpty():
                    source_rect = molecule_bounds.adjusted(-20, -20, 20, 20)
                    w, h = int(source_rect.width()), int(source_rect.height())
                    if w > 0 and h > 0:
                        img = QImage(w, h, QImage.Format.Format_ARGB32)
                        img.fill(Qt.GlobalColor.white)
                        painter = QPainter()
                        if not img.isNull() and painter.begin(img):
                            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
                            scene.render(painter, QRectF(0, 0, w, h), source_rect)
                            painter.end()
                            ba = QByteArray()
                            buf = QBuffer(ba)
                            buf.open(QIODevice.OpenModeFlag.WriteOnly)
                            img.save(buf, "PNG")
                            return ba.toBase64().data().decode("utf-8"), w, h
            except Exception as e:
                logging.warning("Scene Capture Error: %s", e)

        # 2. Fallback to RDKit Draw
        if Draw:
            try:
                import base64
                from io import BytesIO

                img = Draw.MolToImage(self.mol, size=(300, 300))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                return base64.b64encode(buffered.getvalue()).decode("utf-8"), 300, 300
            except (ImportError, RuntimeError, AttributeError, TypeError, ValueError) as e:
                logging.warning("Mol Image Error: %s", e)

        return "", 0, 0

    def build_report_html(self, mol_b64, disp_w, disp_h):
        mw_val = self.lbl_mw.text().split("</b>", 1)[-1].strip()
        em_val = self.lbl_em.text().split("</b>", 1)[-1].strip()
        solvate = self.solvate_suffix().lstrip("·") or "None"

        cell = '<td style="border-bottom: 1px solid #ddd;" align="center">{}</td>'
        comp_rows = "".join(
            ('<tr style="background-color: #f2f2f2;">' if i % 2 == 0 else "<tr>")
            + cell.format(sym)
            + cell.format(format_count(count))
            + cell.format(f"{pct:.2f}")
            + "</tr>"
            for i, (sym, count, _contrib, pct) in enumerate(self.rows)
        )

        structure = (
            f'<img src="data:image/png;base64,{mol_b64}" width="{disp_w}" height="{disp_h}">'
            if (mol_b64 and disp_w > 0)
            else '<p style="color:#888;">No Structure</p>'
        )

        return f"""
        <html>
        <body style="font-family: sans-serif;">
            <h1 style="color: #333; border-bottom: 2px solid #007bff; padding-bottom: 10px; text-align: center;">
                Elemental Analysis Report
            </h1>
            <br>
            <table width="100%" cellpadding="0" cellspacing="0" border="0">
                <tr>
                    <td width="50%" valign="top" style="padding-bottom: 5px; padding-right: 20px;">
                        <h3 style="margin: 0; color: #555; text-align: center;">Parameters</h3>
                    </td>
                    <td width="50%" valign="top" style="padding-bottom: 5px;">
                        <h3 style="margin: 0; color: #555; text-align: center;">Structure</h3>
                    </td>
                </tr>
                <tr>
                    <td width="50%" valign="middle" style="padding-right: 20px;">
                        <table width="100%" cellpadding="5" cellspacing="0" style="border: 1px solid #ddd; border-collapse: collapse;">
                            <tr style="background-color: #f2f2f2;">
                                <td width="140" style="border-bottom: 1px solid #ddd; font-weight: bold;">Formula</td>
                                <td style="border-bottom: 1px solid #ddd;">{self.total_formula}</td>
                            </tr>
                            <tr>
                                <td style="border-bottom: 1px solid #ddd; font-weight: bold;">Solvate</td>
                                <td style="border-bottom: 1px solid #ddd;">{solvate}</td>
                            </tr>
                            <tr style="background-color: #f2f2f2;">
                                <td style="border-bottom: 1px solid #ddd; font-weight: bold;">Molecular Weight</td>
                                <td style="border-bottom: 1px solid #ddd;">{mw_val}</td>
                            </tr>
                            <tr>
                                <td style="border-bottom: 1px solid #ddd; font-weight: bold;">Exact Mass (M)</td>
                                <td style="border-bottom: 1px solid #ddd;">{em_val}</td>
                            </tr>
                        </table>
                    </td>
                    <td width="50%" valign="middle" align="center">
                        <table width="270" height="200" cellpadding="0" cellspacing="0" border="0">
                            <tr><td align="center" valign="middle" width="270" height="200">{structure}</td></tr>
                        </table>
                    </td>
                </tr>
            </table>

            <br>
            <h2>Composition</h2>
            <table width="60%" cellpadding="5" cellspacing="0" style="border: 1px solid #ddd; border-collapse: collapse;">
                <tr style="background-color: #e6e6e6;">
                    <th style="border-bottom: 1px solid #ddd;">Element</th>
                    <th style="border-bottom: 1px solid #ddd;">Count</th>
                    <th style="border-bottom: 1px solid #ddd;">Mass %</th>
                </tr>
                {comp_rows}
            </table>
            <p>{self.lbl_report.text()}</p>
        </body>
        </html>
        """

    def create_report(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save PDF Report", "elemental_analysis_report.pdf", "PDF Files (*.pdf)"
        )
        if not filename:
            return

        mol_b64, img_w, img_h = self._molecule_image_b64()
        disp_w, disp_h = 0, 0
        if mol_b64 and img_w > 0 and img_h > 0:
            ratio = min(260 / img_w, 190 / img_h)
            disp_w, disp_h = int(img_w * ratio), int(img_h * ratio)

        html = self.build_report_html(mol_b64, disp_w, disp_h)

        try:
            printer = QPrinter()
            printer.setOutputFormat(QPrinter.OutputFormat.PdfFormat)
            printer.setOutputFileName(filename)
            printer.setResolution(300)

            layout = QPageLayout(
                QPageSize(QPageSize.PageSizeId.A4),
                QPageLayout.Orientation.Portrait,
                QMarginsF(30, 30, 30, 30),
            )
            layout.setUnits(QPageLayout.Unit.Millimeter)
            printer.setPageLayout(layout)

            doc = QTextDocument()
            doc.setHtml(html)

            paint_rect = layout.paintRectPixels(printer.resolution())
            # Map logical screen pixels (96 DPI) to printer resolution
            scale_factor = printer.resolution() / 96.0
            doc.setPageSize(
                QSizeF(
                    paint_rect.width() / scale_factor,
                    paint_rect.height() / scale_factor,
                )
            )

            painter = QPainter(printer)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.save()
            painter.scale(scale_factor, scale_factor)
            doc.drawContents(painter)
            painter.restore()

            footer_text = (
                "Generated by MoleditPy Elemental Analysis Plugin (Ver. "
                + PLUGIN_VERSION
                + ")"
            )
            painter.setFont(QFont("Arial", 9))
            painter.setPen(QColor("#888888"))
            footer_rect = QRectF(0, paint_rect.height() - 100, paint_rect.width(), 100)
            painter.drawText(
                footer_rect,
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignBottom,
                footer_text,
            )
            painter.end()

            QMessageBox.information(
                self, "Report Saved", f"Report saved to:\n{filename}"
            )
        except Exception as e:
            QMessageBox.critical(
                self, "Report Logic Error", f"Failed to generate PDF content: {e}"
            )


class CopyableTable(QTableWidget):
    """QTableWidget whose Ctrl+C copies the selected cells as TSV."""

    def keyPressEvent(self, event):
        if event.matches(QKeySequence.StandardKey.Copy):
            self.copy_selection()
            return
        super().keyPressEvent(event)

    def copy_selection(self):
        ranges = self.selectedRanges()
        if not ranges:
            return
        rows = sorted({r for rg in ranges for r in range(rg.topRow(), rg.bottomRow() + 1)})
        cols = sorted(
            {c for rg in ranges for c in range(rg.leftColumn(), rg.rightColumn() + 1)}
        )
        lines = []
        for r in rows:
            cells = []
            for c in cols:
                item = self.item(r, c)
                cells.append(item.text() if item is not None and item.isSelected() else "")
            lines.append("\t".join(cells))
        QApplication.clipboard().setText("\n".join(lines))


def initialize(context):
    """
    Initialize the Elemental Analysis plugin.
    """

    def toggle_window():
        if Chem is None:
            context.show_status_message("RDKit is not available.", 5000)
            return

        win = context.get_window("main_panel")
        if win:
            win.show()
            win.raise_()
            win.activateWindow()
            win.check_update()
            if win.sync_check.isChecked() and not win.timer.isActive():
                win.timer.start(500)
            return

        new_win = ElementalAnalysisDialog(context)
        new_win.show()

    context.add_analysis_tool("Elemental Analysis", toggle_window)


class StandaloneContext:
    """Minimal stand-in for PluginContext so the dialog runs without MoleditPy."""

    current_molecule = None
    scene = None

    def __init__(self):
        self._windows = {}

    def get_main_window(self):
        return None

    def register_window(self, key, win):
        self._windows[key] = win

    def get_window(self, key):
        return self._windows.get(key)


if __name__ == "__main__":
    # Standalone: python elemental_analysis.py [formula]
    app = QApplication(sys.argv)
    dialog = ElementalAnalysisDialog(StandaloneContext())
    dialog.formula_input.setText(sys.argv[1] if len(sys.argv) > 1 else "C8H10N4O2")
    dialog.show()
    sys.exit(app.exec())
