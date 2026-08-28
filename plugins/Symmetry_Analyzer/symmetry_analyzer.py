# --- Plugin Metadata ---
PLUGIN_NAME = "Symmetry Analyzer"
PLUGIN_VERSION = "2026.08.29"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Analyzes molecular symmetry (point group) and symmetrizes structures. Refactored for MoleditPy V3.0 API."
PLUGIN_CONTEXT = None

import logging
import numpy as np
from PyQt6.QtWidgets import (
    QVBoxLayout,
    QGridLayout,
    QLabel,
    QDoubleSpinBox,
    QPushButton,
    QListWidget,
    QTextEdit,
    QGroupBox,
    QMessageBox,
    QSplitter,
    QDialog,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal

# --- RDKit Imports ---
try:
    from rdkit.Geometry import Point3D
except ImportError:
    pass  # always present in a MoleditPy environment

# --- Pymatgen Imports ---
try:
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer

    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False


class SymmetryAnalysisWorker(QThread):
    """Background worker for symmetry analysis to prevent UI freezing"""

    # Deliberately not called `finished`: QThread already has a signal by that
    # name, and shadowing it hands anyone connecting to thread completion this
    # scan's two-argument payload instead.
    analysis_finished = pyqtSignal(dict, bool)  # processed_data, found_any

    def __init__(self, mol_pmg, min_tol, max_tol, step=0.005):
        super().__init__()
        self.mol_pmg = mol_pmg
        self.min_tol = min_tol
        self.max_tol = max_tol
        self.step = step
        self._abort = False

    def abort(self):
        """Ask the scan to stop. quit() cannot: run() has no event loop."""
        self._abort = True

    def run(self):
        # A near-symmetric geometry is only reported as its true group over a
        # window a few hundredths of an Angstrom wide -- sometimes thinner than
        # one step -- so the grid is a starting point, not the answer: the
        # refinement pass below bisects whatever falls between two steps.
        step = max(float(self.step), 1e-4)
        tolerances = [
            float(t)
            for t in np.round(np.arange(self.min_tol, self.max_tol + 1e-9, step), 6)
            if t >= self.min_tol - 1e-9
        ]

        group_data = {}
        found_any = False
        seen = {}  # tolerance -> symbol, for the tolerances pymatgen resolved

        def classify(tol_val):
            """Point group at one tolerance, or None if pymatgen refuses it."""
            try:
                analyzer = PointGroupAnalyzer(self.mol_pmg, tolerance=tol_val)
                sym = analyzer.sch_symbol
            except ValueError as exc:
                # Routine: pymatgen's axis search raises ValueError ("min()
                # arg is an empty sequence") at tolerances that fit no axis
                # for this geometry. Neighbouring tolerances still resolve
                # the group -- ammonia fails 30 of 40 steps and is still
                # correctly reported as C3v.
                logging.debug(
                    "Symmetry: tolerance %.4f rejected by pymatgen: %s", tol_val, exc
                )
                return None
            except Exception:
                # Anything else is not expected from a point-group search.
                # Keep scanning the remaining tolerances -- one bad step
                # should not cost the whole result -- but log the traceback
                # instead of hiding it, which is how this stayed invisible.
                logging.exception(
                    "Symmetry: unexpected failure at tolerance %.4f", tol_val
                )
                return None

            if sym not in group_data:
                group_data[sym] = {"analyzer": analyzer, "tols": [tol_val]}
            else:
                group_data[sym]["tols"].append(tol_val)
            seen[tol_val] = sym
            return sym

        try:
            for tol in tolerances:
                if self._abort:
                    break
                if classify(float(tol)) is not None:
                    found_any = True

            # A group can occupy a window narrower than one grid step and be
            # stepped over entirely, which is how a C2 structure came back as
            # D2 with no C2 row at all. Bisect every transition instead of
            # trusting the grid to land inside the window.
            budget = 60
            edges = sorted(seen)
            work = [(a, b, 0) for a, b in zip(edges, edges[1:]) if seen[a] != seen[b]]
            while work and budget > 0 and not self._abort:
                lo, hi, depth = work.pop(0)
                if depth >= 6 or hi - lo <= 1e-4:
                    continue
                mid = round((lo + hi) / 2.0, 6)
                if mid <= lo or mid >= hi:
                    continue
                budget -= 1
                sym = classify(mid)
                if sym is None:
                    continue
                if sym != seen[lo]:
                    work.append((lo, mid, depth + 1))
                if sym != seen[hi]:
                    work.append((mid, hi, depth + 1))
        except Exception:
            # The dialog only unblocks when `analysis_finished` is emitted, so
            # an unexpected error must still report whatever the scan found
            # rather than leaving the UI waiting forever.
            logging.exception("Symmetry: scan aborted early; reporting partial results")

        # A group can drop out and come back (hydrazine is C2, then D2, then C2
        # again), so a single min-max span would hide the group in between.
        scan_order = [(seen[t], t) for t in sorted(seen)]
        for sym, data in group_data.items():
            data["tols"].sort()
            bands = []
            for i, (found_sym, tol_val) in enumerate(scan_order):
                if found_sym != sym:
                    continue
                if bands and scan_order[i - 1][0] == sym:
                    bands[-1][1] = tol_val
                else:
                    bands.append([tol_val, tol_val])
            data["bands"] = [tuple(b) for b in bands]

        self.analysis_finished.emit(group_data, found_any)


class SymmetryAnalysisPlugin(QDialog):
    """
    MoleditPy Plugin: Molecular Symmetry Analyzer & Symmetrizer

    Features:
    1. Point group determination
    2. Listing of the symmetry operations (rotations, mirrors, ...)
    3. Tolerance control
    4. Symmetrization of the structure by the projection-operator method
    """

    def __init__(self, context):
        """
        :param context: MoleditPy PluginContext
        """
        # Set parent to main window for stability
        super().__init__(parent=context.get_main_window())
        self.setWindowFlags(Qt.WindowType.Window)  # Independent window

        self.context = context
        self.analyzer = None  # pymatgen PointGroupAnalyzer instance
        self.symmetry_ops = []  # Detected symmetry operations
        self.worker = None  # QThread instance
        self._scanned_fingerprint = None  # geometry the operations belong to

        # Lighter selection colour (Light Sky Blue)
        self.setStyleSheet("""
            QListWidget::item:selected {
                background-color: #87CEFA;
                color: black;
            }
            QListWidget::item {
                padding: 1px;
                margin: 0px;
            }
        """)

        self.init_ui()

        # Namespaced window registration for V3 lifecycle management
        self.context.register_window("main_panel", self)
        self.context.register_document_reset_handler(self._invalidate_analysis)

        # Initial view update via context
        self.context.refresh_3d_view()

    def init_ui(self):
        main_layout = QVBoxLayout(self)

        # --- 1. Top Settings (Analyze & Symmetrize only) ---
        settings_group = QGroupBox("Actions")
        settings_layout = QGridLayout()

        # Min Tolerance Input
        min_tol_label = QLabel("Min Tol (Å):")
        self.min_tol_spin = QDoubleSpinBox()
        self.min_tol_spin.setRange(0.001, 10.0)
        self.min_tol_spin.setSingleStep(0.005)
        self.min_tol_spin.setDecimals(3)
        self.min_tol_spin.setValue(0.01)  # Default
        self.min_tol_spin.setToolTip(
            "Tightest tolerance to test. A near-symmetric geometry is reported "
            "as the higher point group at any tolerance above its distortion."
        )

        settings_layout.addWidget(min_tol_label, 0, 0)
        settings_layout.addWidget(self.min_tol_spin, 0, 1)

        # Max Tolerance Input
        tol_label = QLabel("Max Tol (Å):")
        self.max_tol_spin = QDoubleSpinBox()
        self.max_tol_spin.setRange(0.01, 10.0)
        self.max_tol_spin.setSingleStep(0.05)
        self.max_tol_spin.setDecimals(2)
        self.max_tol_spin.setValue(0.5)  # Default
        self.max_tol_spin.setToolTip(
            "Loosest tolerance to test. Beyond ~0.7 Å distinct atoms merge "
            "and the reported groups are artefacts."
        )

        settings_layout.addWidget(tol_label, 1, 0)
        settings_layout.addWidget(self.max_tol_spin, 1, 1)

        # Scan Step Input
        step_label = QLabel("Step (Å):")
        self.step_spin = QDoubleSpinBox()
        self.step_spin.setRange(0.001, 1.0)
        self.step_spin.setSingleStep(0.005)
        self.step_spin.setDecimals(3)
        self.step_spin.setValue(0.005)  # Default
        self.step_spin.setToolTip(
            "Spacing of the scanned tolerances. Groups that fall between two "
            "steps are still found: every transition is bisected."
        )

        settings_layout.addWidget(step_label, 2, 0)
        settings_layout.addWidget(self.step_spin, 2, 1)

        # Analyze Button
        self.calc_btn = QPushButton("Analyze (Scan)")
        self.calc_btn.setToolTip("Scan tolerances to find likely point groups.")
        self.calc_btn.clicked.connect(self.on_analyze_clicked)

        # Symmetrize Button
        self.sym_btn = QPushButton("Symmetrize Detected")
        self.sym_btn.setToolTip(
            "Symmetrize structure to match the selected point group."
        )
        self.sym_btn.clicked.connect(self.symmetrize_structure)
        self.sym_btn.setEnabled(False)

        settings_layout.addWidget(self.calc_btn, 3, 0)
        settings_layout.addWidget(self.sym_btn, 3, 1)

        settings_group.setLayout(settings_layout)
        main_layout.addWidget(settings_group)

        # --- 2. Results Area (Splitter) ---
        splitter = QSplitter(Qt.Orientation.Vertical)

        # A. Likely Groups List
        groups_box = QGroupBox("1. Likely Point Groups (Select one)")
        groups_layout = QVBoxLayout()

        self.groups_list = QListWidget()
        self.groups_list.setSelectionMode(QListWidget.SelectionMode.SingleSelection)
        self.groups_list.itemSelectionChanged.connect(self.on_group_selected)
        groups_layout.addWidget(self.groups_list)
        groups_box.setLayout(groups_layout)
        splitter.addWidget(groups_box)

        # B. Operations List
        ops_box = QGroupBox("2. Symmetry Operations")
        ops_layout = QVBoxLayout()

        # Selected Group Display
        self.selected_group_label = QLabel("Point Group: -")
        self.selected_group_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.selected_group_label.setStyleSheet(
            "QLabel { font-size: 12pt; color: #2c3e50; margin: 2px; }"
        )
        ops_layout.addWidget(self.selected_group_label)

        self.ops_list = QListWidget()
        self.ops_list.setSelectionMode(
            QListWidget.SelectionMode.ExtendedSelection
        )  # multi-select with Ctrl/Shift
        self.ops_list.setAlternatingRowColors(True)
        self.ops_list.itemSelectionChanged.connect(self.on_op_selection_changed)
        ops_layout.addWidget(self.ops_list)
        ops_box.setLayout(ops_layout)
        splitter.addWidget(ops_box)

        # C. Matrix Details
        details_box = QGroupBox("3. Operation Details")
        details_layout = QVBoxLayout()

        self.op_details = QTextEdit()
        self.op_details.setReadOnly(True)
        self.op_details.setPlaceholderText(
            "Select an operation above to view matrix details."
        )
        details_layout.addWidget(self.op_details)
        details_box.setLayout(details_layout)
        splitter.addWidget(details_box)

        # Set initial sizes for splitter (optional)
        splitter.setSizes([150, 200, 150])

        main_layout.addWidget(splitter)

        # Check Dependency
        if not HAS_PYMATGEN:
            self.calc_btn.setEnabled(False)
            QMessageBox.critical(
                self,
                "Dependency Error",
                "This plugin requires 'pymatgen'.\nPlease install it via: pip install pymatgen",
            )

        # Data storage
        self.group_data = {}  # Symbol -> {'analyzer': obj, 'tols': [float]}

    def get_pymatgen_molecule(self):
        """Convert the MoleditPy (RDKit) molecule to pymatgen form."""
        rd_mol = self.context.current_molecule
        if rd_mol is None or rd_mol.GetNumAtoms() == 0:
            return None

        try:
            conf = rd_mol.GetConformer()
        except ValueError:
            return None  # no conformer

        species = [atom.GetAtomicNum() for atom in rd_mol.GetAtoms()]
        coords = [list(conf.GetAtomPosition(i)) for i in range(rd_mol.GetNumAtoms())]

        return Molecule(species, coords)

    def _molecule_fingerprint(self):
        """Identity of the current molecule: its elements and coordinates."""
        rd_mol = self.context.current_molecule
        if rd_mol is None:
            return None

        try:
            conf = rd_mol.GetConformer()
        except ValueError:
            return None

        return (
            tuple(atom.GetAtomicNum() for atom in rd_mol.GetAtoms()),
            tuple(
                round(float(v), 6)
                for i in range(rd_mol.GetNumAtoms())
                for v in list(conf.GetAtomPosition(i))
            ),
        )

    def on_analyze_clicked(self):
        """The Analyze button doubles as the Stop button while a scan runs."""
        if self.worker is not None and self.worker.isRunning():
            self.stop_analysis()
        else:
            self.analyze_symmetry()

    def stop_analysis(self):
        """Ask the running scan to stop; it still reports what it found."""
        if self.worker is None or not self.worker.isRunning():
            return
        self.calc_btn.setText("Stopping...")
        self.calc_btn.setEnabled(False)
        self.worker.abort()

    def analyze_symmetry(self):
        """Scan the point group over a range of tolerances and list the results."""
        mol_pmg = self.get_pymatgen_molecule()
        if mol_pmg is None:
            QMessageBox.warning(self, "Error", "No molecule to analyze.")
            return

        self.groups_list.clear()
        self.ops_list.clear()
        self.selected_group_label.setText("Point Group: -")  # Reset to placeholder
        self.op_details.clear()
        self.sym_btn.setEnabled(False)
        self.group_data = {}

        # The button becomes the way to stop the scan it started.
        self.calc_btn.setText("Stop")
        self.calc_btn.setToolTip("Stop the scan and report what it found so far.")

        # The spin box floor keeps this above zero; a zero tolerance degenerates
        # to "exact coordinates only" and always reports C1.
        min_tol = self.min_tol_spin.value()
        max_tol = max(self.max_tol_spin.value(), min_tol)
        step = self.step_spin.value()

        self._scanned_fingerprint = self._molecule_fingerprint()

        # Start Worker Thread
        self.worker = SymmetryAnalysisWorker(mol_pmg, min_tol, max_tol, step)
        self.worker.analysis_finished.connect(self.on_analysis_finished)
        self.worker.start()

    def on_analysis_finished(self, group_data, found_any):
        """Handle the finished scan."""
        self.calc_btn.setEnabled(True)
        self.calc_btn.setText("Analyze (Scan)")
        self.calc_btn.setToolTip("Scan tolerances to find likely point groups.")
        self.group_data = group_data

        if not found_any:
            self.groups_list.addItem("No point groups found.")
            return

        # Strictest first: the group found at the smallest tolerance is the
        # one the geometry actually matches, and a wider range breaks ties.
        sorted_keys = sorted(
            self.group_data.keys(),
            key=lambda k: (
                min(self.group_data[k]["tols"]),
                -max(self.group_data[k]["tols"]),
            ),
        )

        for sym in sorted_keys:
            data = self.group_data[sym]
            bands = data.get("bands") or [(min(data["tols"]), max(data["tols"]))]

            # One plain row per group: "Td (Tol: 0.100 - 2.000)"
            # 4 decimals: the refinement pass locates a boundary to ~1e-4, and
            # at 3 the two bands either side of it print the same edge.
            span = ", ".join(f"{lo:.4f} - {hi:.4f}" for lo, hi in bands)
            item_text = f"{sym}  (Tol: {span} Å)"
            self.groups_list.addItem(item_text)

        # Auto-select the first row, but prefer anything over C1
        if self.groups_list.count() > 0:
            target_row = 0
            for i, sym in enumerate(sorted_keys):
                if sym != "C1":
                    target_row = i
                    break
            self.groups_list.setCurrentRow(target_row)

        # QMessageBox.information(self, "Done",
        #    f"Found {len(self.group_data)} potential point groups.\n"
        #    "Sorted by strictness (smaller tolerance first).")

    def _get_op_sort_key(self, op):
        """
        Sort key for the symmetry operations.
        Order: Identity -> Rotation(Cn) -> Reflection(sigma) -> Inversion(i) -> Improper(Sn)
        """
        m = op.rotation_matrix
        det = np.linalg.det(m)
        trace = np.trace(m)
        tol = 1e-2

        # 1. Identity (Priority: 0)
        if np.allclose(m, np.eye(3), atol=tol):
            return (0, 0)

        # 2. Proper Rotation (Priority: 1)
        if np.isclose(det, 1.0, atol=tol):
            val = np.clip((trace - 1) / 2.0, -1.0, 1.0)
            angle = np.degrees(np.arccos(val))
            order = round(360.0 / angle) if angle > 1.0 else 1
            return (1, -order)

        # 3. Reflection (Priority: 2)
        if np.isclose(trace, 1.0, atol=tol):
            return (2, 0)

        # 4. Inversion (Priority: 3)
        if np.isclose(trace, -3.0, atol=tol):
            return (3, 0)

        # 5. Improper Rotation (Priority: 4)
        val = np.clip((trace + 1) / 2.0, -1.0, 1.0)
        angle = np.degrees(np.arccos(val))
        order = round(360.0 / angle) if angle > 1.0 else 1

        # S2 (Inversion) check for safety -> Priority 3
        if order == 2:
            return (3, 0)

        return (4, -order)

    def on_group_selected(self):
        """Show the operations of the selected group."""
        item = self.groups_list.currentItem()
        if not item:
            return

        text = item.text()
        sym = text.split()[0]
        if sym not in self.group_data:
            # "No point groups found." is a selectable row; it is not a group.
            return

        # Display nicely formatted symbol
        s = self._format_symmetry_symbol(sym)
        if s.startswith("<html>"):
            # Keep HTML structure valid
            inner = s.replace("<html>", "").replace("</html>", "")
            self.selected_group_label.setText(f"<html>Point Group: {inner}</html>")
        else:
            self.selected_group_label.setText(f"Point Group: {s}")

        self.analyzer = self.group_data[sym]["analyzer"]
        ops = self.analyzer.get_symmetry_operations()
        # Identity -> Rotation -> Reflection -> Inversion -> Improper
        ops.sort(key=self._get_op_sort_key)
        self.symmetry_ops = ops

        self.update_ops_list()
        self.sym_btn.setEnabled(True)
        self.op_details.clear()

    def _format_symmetry_symbol(self, sym):
        """Render a Schoenflies symbol as HTML (italic letter + subscript)."""
        import re
        # C2v -> C, 2v
        # D3h -> D, 3h
        # Td  -> T, d
        # C*v -> C, *v -> C, ∞v

        match = re.match(r"^([A-Z]+)(.*)$", sym)
        if match:
            main = match.group(1)
            sub = match.group(2)

            # infinity axis
            sub = sub.replace("*", "∞")

            # Font tuning lives in the style sheet; this is the markup only.
            return f"<html><i>{main}</i><sub>{sub}</sub></html>"
        return sym

    def update_ops_list(self):
        """Refresh the operations list widget."""
        self.ops_list.clear()
        self.op_details.clear()

        for i, op in enumerate(self.symmetry_ops):
            label = self._get_op_label(op, i)
            self.ops_list.addItem(label)

    def _get_op_label(self, op, i):
        """
        Label an operation from its rotation matrix (Unicode notation).
        The classification order matches the sort key.
        Order: Identity -> Rotation -> Reflection -> Inversion -> Improper
        """
        m = op.rotation_matrix
        det = np.linalg.det(m)
        trace = np.trace(m)
        tol = 1e-2

        # Subscript digits
        def to_sub(n):
            return str(n).translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))

        # 1. Identity
        if np.allclose(m, np.eye(3), atol=tol):
            return f"#{i + 1}: E (Identity)"

        # 2. Proper rotation
        # Det = +1
        if np.isclose(det, 1.0, atol=tol):
            val = (trace - 1) / 2.0
            val = np.clip(val, -1.0, 1.0)
            angle = np.degrees(np.arccos(val))

            # an angle of ~0 is the identity
            if angle < 1.0:
                return f"#{i + 1}: E (Identity)"

            order = round(360.0 / angle)
            return f"#{i + 1}: C{to_sub(order)} (Rotation)"

        else:
            # Improper operations (rotoreflection, mirror, inversion)
            # Det = -1

            # 3. Reflection (trace = 1)
            if np.isclose(trace, 1.0, atol=tol):
                return f"#{i + 1}: σ (Reflection)"

            # 4. Inversion (trace = -3)
            if np.isclose(trace, -3.0, atol=tol):
                return f"#{i + 1}: i (Inversion)"

            # 5. Improper rotation
            val = (trace + 1) / 2.0
            val = np.clip(val, -1.0, 1.0)
            angle = np.degrees(np.arccos(val))

            if angle < 1.0:
                return f"#{i + 1}: σ (Reflection)"

            order = round(360.0 / angle)

            # S2 is the inversion
            if order == 2:
                return f"#{i + 1}: i (Inversion)"

            return f"#{i + 1}: S{to_sub(order)} (Improper Rotation)"

    def on_op_selection_changed(self):
        """Handle a change of the operation selection (multi-select)."""
        items = self.ops_list.selectedItems()

        ops_to_show = []
        for item in items:
            row = self.ops_list.row(item)
            if 0 <= row < len(self.symmetry_ops):
                ops_to_show.append(self.symmetry_ops[row])

        # Update the 3D visualisation
        self.visualize_ops(ops_to_show)

        # Update the detail text
        if len(ops_to_show) == 0:
            self.op_details.clear()
        elif len(ops_to_show) == 1:
            # a single selection gets its full details
            self._display_single_op_details(ops_to_show[0])
        else:
            self.op_details.setText(
                f"{len(ops_to_show)} operations selected.\n"
                "See 3D view for visualization."
            )

    def _display_single_op_details(self, op):
        """Show the details of one operation."""
        mat_str = np.array2string(op.rotation_matrix, precision=3, suppress_small=True)
        trans_str = np.array2string(
            op.translation_vector, precision=3, suppress_small=True
        )

        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            xyz_str = op.as_xyz_str()

        text = (
            f"--- Rotation Matrix ---\n{mat_str}\n\n"
            f"--- Translation Vector ---\n{trans_str}\n\n"
            f"Type: {xyz_str}"
        )
        self.op_details.setText(text)

    def visualize_ops(self, ops_list):
        """Draw the given symmetry operations in the 3D view."""
        if not self.context.plotter:
            return

        import numpy as np

        try:
            import pyvista as pv
        except ImportError:
            return

        plotter = self.context.plotter

        # Clear the previous drawing
        if getattr(self, "vis_actors", None) is None:
            self.vis_actors = []

        for actor in self.vis_actors:
            plotter.remove_actor(actor)
        self.vis_actors = []

        if not ops_list:
            plotter.render()
            return

        # Centre of mass, shared by every element drawn
        rd_mol = self.context.current_molecule

        if rd_mol:
            conf = rd_mol.GetConformer()
            positions = [conf.GetAtomPosition(i) for i in range(rd_mol.GetNumAtoms())]
            coords = np.array([[p.x, p.y, p.z] for p in positions])
            # Mass-weighted: pymatgen derives the operations in the
            # centre-of-mass frame, so the elements have to be drawn through
            # that point. This was the unweighted centroid.
            masses = np.array(
                [atom.GetMass() for atom in rd_mol.GetAtoms()], dtype=float
            )
            total_mass = masses.sum()
            if total_mass > 0:
                com = (coords * masses[:, None]).sum(axis=0) / total_mass
            else:
                com = np.mean(coords, axis=0)
            # molecule size (largest radius)
            dists = np.linalg.norm(coords - com, axis=1)
            mol_radius = np.max(dists) if len(dists) > 0 else 2.0
        else:
            com = np.array([0.0, 0.0, 0.0])
            mol_radius = 2.0

        # leave a little margin
        scale = mol_radius

        # Draw each operation
        for op in ops_list:
            self._add_op_visualization(plotter, op, com, pv, scale)

        plotter.render()

    def _add_op_visualization(self, plotter, op, com, pv, scale=2.0):
        """Add one operation to the scene (internal helper)."""
        m = op.rotation_matrix
        det = np.linalg.det(m)
        trace = np.trace(m)

        # 1. Identity
        if np.allclose(m, np.eye(3)):
            return

        # 2. Inversion
        if np.isclose(trace, -3.0, atol=1e-2):
            r = 0.25
            sphere = pv.Sphere(radius=r, center=com)
            actor = plotter.add_mesh(
                sphere, color="orange", opacity=0.8, name="sym_inversion"
            )
            self.vis_actors.append(actor)
            return

        # 3. Rotation / Reflection
        eigvals, eigvecs = np.linalg.eig(m)
        axis_idx = np.where(np.isclose(eigvals, 1.0))[0]
        normal_idx = np.where(np.isclose(eigvals, -1.0))[0]

        # A. Proper Rotation
        if np.isclose(det, 1.0) and len(axis_idx) > 0:
            axis = np.real(eigvecs[:, axis_idx[0]])
            if np.linalg.norm(axis) < 1e-3:
                return
            axis = axis / np.linalg.norm(axis)

            length = scale * 1.5
            start = com - axis * length
            end = com + axis * length
            line = pv.Line(start, end)
            actor_line = plotter.add_mesh(line, color="cyan", line_width=4)
            self.vis_actors.append(actor_line)

        # B. Mirror Reflection (sigma = S1: eigenvalues 1, 1, -1 -> trace 1)
        elif len(normal_idx) > 0 and np.isclose(trace, 1.0):
            normal = np.real(eigvecs[:, normal_idx[0]])
            if np.linalg.norm(normal) < 1e-3:
                return
            normal = normal / np.linalg.norm(normal)

            size = scale * 1.5
            disk = pv.Disc(center=com, inner=0, outer=size, normal=normal, c_res=30)
            actor_plane = plotter.add_mesh(disk, color="magenta", opacity=0.3)
            self.vis_actors.append(actor_plane)

        # C. Improper rotation S_n (n > 2). det = -1 with trace != 1, so it is
        # neither a plane nor an inversion and used to fall through both
        # branches undrawn -- methane lost all 6 of its S4 operations.
        # S_n maps its own axis v to -v, so the axis is the -1 eigenvector.
        elif np.isclose(det, -1.0) and len(normal_idx) > 0:
            axis = np.real(eigvecs[:, normal_idx[0]])
            if np.linalg.norm(axis) < 1e-3:
                return
            axis = axis / np.linalg.norm(axis)

            length = scale * 1.5
            line = pv.Line(com - axis * length, com + axis * length)
            actor_line = plotter.add_mesh(
                line, color="yellow", line_width=3, style="wireframe"
            )
            self.vis_actors.append(actor_line)

    def symmetrize_structure(self):
        if self.analyzer is None:
            return
        mol_pmg = self.get_pymatgen_molecule()
        if mol_pmg is None:
            return

        # The operations belong to the geometry that was scanned; applying them
        # to anything else silently distorts it -- water's C2v moved benzene's
        # atoms 0.15 A and still reported success.
        if self._molecule_fingerprint() != self._scanned_fingerprint:
            QMessageBox.warning(
                self,
                "Analysis Out of Date",
                "The molecule changed since the last scan.\n"
                "Run 'Analyze (Scan)' again before symmetrizing.",
            )
            self._invalidate_analysis()
            return

        original_coords = mol_pmg.cart_coords
        # pymatgen derives the symmetry operations in the centre-of-mass frame,
        # so the coordinates fed through them must be centred the same way.
        # This was the unweighted centroid, which for water sits 0.33 A away.
        center_of_mass = np.array(mol_pmg.center_of_mass)
        centered_coords = original_coords - center_of_mass

        species = mol_pmg.species
        ops = self.analyzer.get_symmetry_operations()
        if not ops:
            return

        try:
            from scipy.optimize import linear_sum_assignment

            HAS_SCIPY = True
        except ImportError:
            HAS_SCIPY = False

        if not HAS_SCIPY:
            logging.warning(
                "Symmetry: scipy is missing; falling back to greedy atom mapping."
            )

        new_coords = np.zeros_like(centered_coords)
        n_ops = len(ops)
        n_atoms = len(centered_coords)
        mapping_error = False

        for op in ops:
            # rotated_coords[j] is where operation R sends atom j
            rotated_coords = np.array([op.operate(p) for p in centered_coords])

            # --- Cost matrix ---
            cost_matrix = np.zeros((n_atoms, n_atoms))
            # Broadcasting: (N, 1, 3) - (1, N, 3) -> (N, N, 3) -> norm -> (N, N)
            # distance between target position i and rotated candidate j
            diff = centered_coords[:, np.newaxis, :] - rotated_coords[np.newaxis, :, :]
            dist_mat = np.linalg.norm(diff, axis=2)

            # elements must match
            sp_array = np.array([s.symbol for s in species])
            species_mask = sp_array[:, np.newaxis] != sp_array[np.newaxis, :]
            cost_matrix = dist_mat
            cost_matrix[species_mask] = 1e9

            # --- Mapping ---
            if HAS_SCIPY:
                row_ind, col_ind = linear_sum_assignment(cost_matrix)
            else:
                # Greedy Fallback
                row_ind = np.arange(n_atoms)
                col_ind = np.full(n_atoms, -1, dtype=int)
                used_j = set()
                for i in range(n_atoms):
                    sorted_j = np.argsort(cost_matrix[i])
                    for j in sorted_j:
                        if j not in used_j and cost_matrix[i, j] < 1.0:
                            col_ind[i] = j
                            used_j.add(j)
                            break
                    if col_ind[i] == -1:
                        mapping_error = True
                        col_ind[i] = i

            # --- Accumulate ---
            for i, j in zip(row_ind, col_ind):
                if cost_matrix[i, j] > 1.5:
                    mapping_error = True

                # BUG FIX:
                # atom i (original frame) is mapped to atom j (rotated frame).
                # Meaning: rotated_coords[j] is the position that corresponds to atom i.
                # We want to average these positions.
                # OLD (Incorrect): new_coords[i] += op.inverse.operate(centered_coords[j])
                # NEW (Correct):
                new_coords[i] += rotated_coords[j]

        new_coords /= n_ops

        if mapping_error:
            # An unmapped atom averages in a position from the wrong site, so
            # the result is not merely imprecise -- it can be a scrambled
            # geometry. Ask before overwriting the user's structure with it.
            reply = QMessageBox.question(
                self,
                "Incomplete Symmetry Mapping",
                "Some atoms could not be matched to a symmetry-equivalent "
                "partner.\nThe structure might be too distorted, or 'scipy' "
                "is missing.\n\nApplying this result may scramble the "
                "geometry. Apply anyway?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if reply != QMessageBox.StandardButton.Yes:
                return

        final_coords = new_coords + center_of_mass
        self.update_rdkit_coords(final_coords)
        self._scanned_fingerprint = self._molecule_fingerprint()

        QMessageBox.information(
            self,
            "Symmetrized",
            f"Structure symmetrized to {self.analyzer.sch_symbol}.\n"
            f"(Averaged over {len(ops)} operations)",
        )

    def update_rdkit_coords(self, new_coords):
        """Write the computed coordinates back to RDKit and refresh the view."""
        rd_mol = self.context.current_molecule
        if rd_mol is None:
            return

        # Save the undo state (Manual Section 4)
        self.context.push_undo_checkpoint()

        conf = rd_mol.GetConformer()
        for i in range(rd_mol.GetNumAtoms()):
            # cast to float: the RDKit C++ API can reject numpy.float64
            x, y, z = map(float, new_coords[i])
            conf.SetAtomPosition(i, Point3D(x, y, z))

        # Full redraw: the conformer changed
        self.context.current_mol = rd_mol
        self.context.refresh_3d_view()

    def _invalidate_analysis(self):
        """Discard the scan results and the 3D overlay; reset the UI."""
        # 1. Remove the 3D visualisation
        if self.context.plotter and getattr(self, "vis_actors", None) is not None:
            plotter = self.context.plotter
            for actor in self.vis_actors:
                plotter.remove_actor(actor)
            self.vis_actors = []
            plotter.render()

        # 2. Reset the UI and the internal data
        self.groups_list.clear()
        self.ops_list.clear()
        self.selected_group_label.setText("Point Group: -")
        self.op_details.clear()
        self.sym_btn.setEnabled(False)
        self.group_data = {}
        self.symmetry_ops = []
        self.analyzer = None
        self._scanned_fingerprint = None

    def _cleanup(self):
        """Everything that must happen when the dialog goes away."""
        if self.worker and self.worker.isRunning():
            self.worker.abort()
            self.worker.wait(5000)
            if self.worker.isRunning():
                self.worker.terminate()  # last resort

        self._invalidate_analysis()

    def closeEvent(self, event):
        self._cleanup()
        # QDialog.closeEvent routes back through reject(), so accepting here
        # rather than calling super() keeps the two paths from recursing.
        event.accept()

    def reject(self):
        # Esc closes a QDialog through reject(), which never delivers a close
        # event: the 3D overlay and the scan thread used to survive it.
        self._cleanup()
        super().reject()


def initialize(context):
    """
    Initialize the Symmetry Analyzer plugin.
    """
    global PLUGIN_CONTEXT
    PLUGIN_CONTEXT = context

    def toggle_window():
        win = context.get_window("main_panel")
        if win:
            win.show()
            win.raise_()
            win.activateWindow()
            return

        new_win = SymmetryAnalysisPlugin(context)
        new_win.setWindowTitle("Symmetry Analyzer")
        new_win.resize(400, 600)
        new_win.show()

    context.add_menu_action("3D Edit/Symmetrize...", toggle_window)
