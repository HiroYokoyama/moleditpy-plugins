"""
xTB Optimizer Plugin for MoleditPy
===================================
Runs semiempirical GFN2-xTB (or GFN1-xTB) geometry optimization using the
tblite package (via its ASE calculator interface) + ASE LBFGS optimizer.

The optimization runs on a background QThread so the UI stays responsive.
A live log, per-step energy table, and a Cancel button are provided in the
dialog.

Dependencies (install in MoleditPy's environment):
    mamba install -c conda-forge tblite-python ase
  or
    pip install tblite ase
"""

from __future__ import annotations

import logging

from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)
from rdkit import Chem

# ---------------------------------------------------------------------------
# Plugin metadata
# ---------------------------------------------------------------------------

PLUGIN_NAME = "xTB Optimizer"
PLUGIN_VERSION = "2026.07.31"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Geometry optimization using semiempirical xTB methods (GFN2-xTB, GFN1-xTB) "
    "via the tblite package. Runs on a background thread. On Windows, install "
    "tblite from conda-forge (mamba install -c conda-forge tblite-python ase) — "
    "PyPI ships no Windows wheels."
)
PLUGIN_TAGS = ["Optimization"]
PLUGIN_DEPENDENCIES = ["tblite", "ase"]

logger = logging.getLogger(__name__)


def derive_charge_and_multiplicity(mol) -> tuple[int, int]:
    """Read total charge and spin multiplicity off an RDKit molecule.

    Returns ``(charge, multiplicity)``.  The multiplicity is derived from the
    explicit radical electrons, then reconciled with the electron count:
    a molecule with an odd number of electrons cannot be a singlet, and
    handing xTB an impossible combination makes it fail outright.
    """
    charge = Chem.GetFormalCharge(mol)
    unpaired = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())

    n_electrons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms()) - charge
    # Unpaired electrons must have the same parity as the electron count.
    if (n_electrons - unpaired) % 2 != 0:
        unpaired = unpaired + 1 if unpaired == 0 else unpaired - 1

    return charge, unpaired + 1


# ---------------------------------------------------------------------------
# Worker thread
# ---------------------------------------------------------------------------


class XtbWorker(QThread):
    """
    Background worker that performs the xTB geometry optimization.

    Signals
    -------
    log_message(str)
        A plain-text line to append to the dialog log.
    step_update(int, float, float)
        Emitted after each optimizer step: (step_number, energy_eV, fmax_eV_per_A).
    finished(bool, object)
        Emitted when the run completes.  ``True`` + numpy array of new positions
        on success; ``False`` + error message string on failure.
    """

    log_message = pyqtSignal(str)
    step_update = pyqtSignal(int, float, float)
    finished = pyqtSignal(bool, object)

    def __init__(
        self,
        numbers: list[int],
        positions: list[list[float]],
        method: str,
        fmax: float,
        max_steps: int,
        charge: int = 0,
        multiplicity: int = 1,
        parent=None,
    ):
        super().__init__(parent)
        self._numbers = numbers
        self._positions = positions  # Å
        self._method = method
        self._fmax = fmax
        self._max_steps = max_steps
        self._charge = charge
        self._multiplicity = multiplicity
        self._cancelled = False

    def cancel(self):
        """Request cancellation — the optimizer loop checks this flag each step."""
        self._cancelled = True

    def run(self):
        """Entry point executed on the worker thread."""
        try:
            import numpy as np
            from ase import Atoms
            from ase.optimize import LBFGS
            from tblite.ase import TBLite

            self.log_message.emit(
                f"Starting {self._method} optimization  "
                f"(fmax={self._fmax} eV/Å, max_steps={self._max_steps}, "
                f"charge={self._charge:+d}, multiplicity={self._multiplicity})"
            )

            atoms = Atoms(
                numbers=self._numbers,
                positions=self._positions,
            )
            # charge/multiplicity must be passed explicitly.  Left at their
            # None defaults, tblite falls back to atoms.get_initial_charges()
            # and get_initial_magnetic_moments(), which are all-zero for an
            # Atoms built from numbers+positions -- so every ion was optimised
            # as neutral and every radical as closed-shell.
            atoms.calc = TBLite(
                method=self._method,
                charge=self._charge,
                multiplicity=self._multiplicity,
                verbosity=0,
            )

            step_counter = [0]
            cancelled_ref = [False]

            class _StepCallback:
                """ASE observer called after every optimizer step."""

                def __init__(cb_self):
                    pass

                def __call__(cb_self):
                    if self._cancelled:
                        cancelled_ref[0] = True
                        raise RuntimeError("Optimization cancelled by user.")

                    step_counter[0] += 1
                    try:
                        energy = float(atoms.get_potential_energy())
                        forces = atoms.get_forces()
                        fmax_val = float(np.sqrt((forces**2).sum(axis=1).max()))
                    except Exception:  # noqa: BLE001 - a bad step callback must not abort the run
                        energy = float("nan")
                        fmax_val = float("nan")

                    self.step_update.emit(step_counter[0], energy, fmax_val)
                    self.log_message.emit(
                        f"  Step {step_counter[0]:>4d}  "
                        f"E = {energy:>14.6f} eV  "
                        f"Fmax = {fmax_val:.4f} eV/Å"
                    )

            opt = LBFGS(atoms, logfile=None)
            opt.attach(_StepCallback())

            converged = opt.run(fmax=self._fmax, steps=self._max_steps)

            if cancelled_ref[0]:
                self.log_message.emit("Optimization cancelled.")
                self.finished.emit(False, "Cancelled by user.")
                return

            new_positions = atoms.get_positions().tolist()
            final_energy = float(atoms.get_potential_energy())

            if converged:
                self.log_message.emit(
                    f"\nConverged in {step_counter[0]} steps. "
                    f"Final E = {final_energy:.6f} eV"
                )
            else:
                self.log_message.emit(
                    f"\nReached max steps ({self._max_steps}) without convergence. "
                    f"Final E = {final_energy:.6f} eV"
                )

            self.finished.emit(True, new_positions)

        except RuntimeError as exc:
            if "cancelled" in str(exc).lower():
                self.log_message.emit("Optimization cancelled.")
                self.finished.emit(False, "Cancelled by user.")
            else:
                logger.exception("XtbWorker: runtime error during optimization")
                self.log_message.emit(f"\nError: {exc}")
                self.finished.emit(False, str(exc))
        except ImportError as exc:
            msg = (
                f"Missing dependency: {exc}\n"
                "Install tblite and ase:\n"
                "  mamba install -c conda-forge tblite-python ase"
            )
            self.log_message.emit(msg)
            self.finished.emit(False, msg)
        except Exception as exc:
            logger.exception("XtbWorker: unexpected error")
            self.log_message.emit(f"\nUnexpected error: {exc}")
            self.finished.emit(False, str(exc))


# ---------------------------------------------------------------------------
# Dialog
# ---------------------------------------------------------------------------


class XtbOptimizerDialog(QDialog):
    """
    Main dialog for the xTB optimizer plugin.

    Provides settings (method, convergence threshold, max steps), a live log,
    a per-step energy table, progress indicator, and Run / Cancel / Close buttons.
    """

    def __init__(self, context, parent=None):
        super().__init__(parent)
        self.context = context
        self.setWindowTitle("xTB Geometry Optimizer")
        self.setMinimumSize(520, 600)
        self.resize(600, 680)

        self._worker: XtbWorker | None = None
        self._step_count = 0

        self._build_ui()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setSpacing(8)

        # ── Settings group ────────────────────────────────────────────
        grp_settings = QGroupBox("Settings")
        form = QFormLayout(grp_settings)

        self.combo_method = QComboBox()
        self.combo_method.addItems(["GFN2-xTB", "GFN1-xTB"])
        form.addRow("Method:", self.combo_method)

        self.chk_auto_charge = QCheckBox("Read charge and spin from the molecule")
        self.chk_auto_charge.setChecked(True)
        self.chk_auto_charge.setToolTip(
            "Uses the RDKit formal charge and radical electrons. Uncheck to set "
            "them by hand, e.g. for a triplet or a charge RDKit cannot perceive."
        )
        self.chk_auto_charge.toggled.connect(self._on_auto_charge_toggled)
        form.addRow("", self.chk_auto_charge)

        self.spin_charge = QSpinBox()
        self.spin_charge.setRange(-20, 20)
        self.spin_charge.setValue(0)
        self.spin_charge.setEnabled(False)
        self.spin_charge.setToolTip("Total molecular charge handed to xTB.")
        form.addRow("Charge:", self.spin_charge)

        self.spin_multiplicity = QSpinBox()
        self.spin_multiplicity.setRange(1, 21)
        self.spin_multiplicity.setValue(1)
        self.spin_multiplicity.setEnabled(False)
        self.spin_multiplicity.setToolTip(
            "Spin multiplicity 2S+1. 1 = closed shell, 2 = doublet radical, "
            "3 = triplet."
        )
        form.addRow("Spin multiplicity:", self.spin_multiplicity)

        self.spin_fmax = QDoubleSpinBox()
        self.spin_fmax.setRange(1e-6, 10.0)
        self.spin_fmax.setDecimals(4)
        self.spin_fmax.setSingleStep(0.01)
        self.spin_fmax.setValue(0.05)
        self.spin_fmax.setToolTip("Convergence threshold (eV/Å). Smaller = tighter.")
        form.addRow("Force threshold (eV/Å):", self.spin_fmax)

        self.spin_maxsteps = QSpinBox()
        self.spin_maxsteps.setRange(1, 10000)
        self.spin_maxsteps.setValue(500)
        self.spin_maxsteps.setToolTip("Maximum number of LBFGS steps.")
        form.addRow("Max steps:", self.spin_maxsteps)

        root.addWidget(grp_settings)

        # ── Step table ────────────────────────────────────────────────
        grp_table = QGroupBox("Optimization Progress")
        tbl_layout = QVBoxLayout(grp_table)
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["Step", "Energy (eV)", "Fmax (eV/Å)"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        tbl_layout.addWidget(self.table)
        root.addWidget(grp_table)

        # ── Log area ──────────────────────────────────────────────────
        grp_log = QGroupBox("Log")
        log_layout = QVBoxLayout(grp_log)
        self.log_text = QPlainTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumBlockCount(2000)
        self.log_text.setPlaceholderText("Optimization log will appear here…")
        log_layout.addWidget(self.log_text)
        root.addWidget(grp_log)

        # ── Progress bar ──────────────────────────────────────────────
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)   # indeterminate
        self.progress_bar.setVisible(False)
        root.addWidget(self.progress_bar)

        # ── Status label ──────────────────────────────────────────────
        self.lbl_status = QLabel("Ready.")
        self.lbl_status.setAlignment(Qt.AlignmentFlag.AlignLeft)
        root.addWidget(self.lbl_status)

        # ── Buttons ───────────────────────────────────────────────────
        btn_row = QHBoxLayout()

        self.btn_run = QPushButton("▶  Run")
        self.btn_run.setDefault(True)
        self.btn_run.clicked.connect(self._on_run)

        self.btn_cancel = QPushButton("✕  Cancel")
        self.btn_cancel.setEnabled(False)
        self.btn_cancel.clicked.connect(self._on_cancel)

        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.accept)

        btn_row.addWidget(self.btn_run)
        btn_row.addWidget(self.btn_cancel)
        btn_row.addStretch()
        btn_row.addWidget(self.btn_close)
        root.addLayout(btn_row)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _on_auto_charge_toggled(self, checked: bool):
        self.spin_charge.setEnabled(not checked)
        self.spin_multiplicity.setEnabled(not checked)

    def _set_running(self, running: bool):
        self.btn_run.setEnabled(not running)
        self.btn_cancel.setEnabled(running)
        self.combo_method.setEnabled(not running)
        self.spin_fmax.setEnabled(not running)
        self.spin_maxsteps.setEnabled(not running)
        self.chk_auto_charge.setEnabled(not running)
        manual = not self.chk_auto_charge.isChecked()
        self.spin_charge.setEnabled(not running and manual)
        self.spin_multiplicity.setEnabled(not running and manual)
        self.progress_bar.setVisible(running)

    def _append_log(self, text: str):
        self.log_text.appendPlainText(text)
        # Scroll to bottom
        sb = self.log_text.verticalScrollBar()
        sb.setValue(sb.maximum())

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _on_run(self):
        mol = self.context.current_mol
        if not mol:
            QMessageBox.warning(self, PLUGIN_NAME, "No molecule loaded.")
            return
        if mol.GetNumConformers() == 0:
            QMessageBox.warning(
                self, PLUGIN_NAME,
                "The molecule has no 3D coordinates.\n"
                "Use '3D Edit ▸ Generate 3D Coordinates' first."
            )
            return

        # Gather atom data from the RDKit molecule
        conf = mol.GetConformer()
        numbers = []
        positions = []
        try:
            from rdkit.Chem import GetPeriodicTable
            pt = GetPeriodicTable()

            # Implicit hydrogens have no coordinates, so they would simply be
            # absent from the structure handed to xTB -- ethanol would be
            # optimised as a bare C-C-O fragment, and it would "succeed".
            # Inside the try: reading a corrupt molecule must still report
            # "Failed to read molecule" rather than raise out of the slot.
            implicit_h = sum(atom.GetNumImplicitHs() for atom in mol.GetAtoms())
            if implicit_h:
                QMessageBox.warning(
                    self, PLUGIN_NAME,
                    f"The molecule has {implicit_h} implicit hydrogen(s), which "
                    "have no 3D coordinates and would be left out of the "
                    "calculation.\n\nAdd explicit hydrogens before optimizing."
                )
                return

            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                if symbol == "*":
                    QMessageBox.warning(
                        self, PLUGIN_NAME,
                        "Molecule contains dummy atoms (*). Remove them before optimizing."
                    )
                    return
                numbers.append(pt.GetAtomicNumber(symbol))
                pos = conf.GetAtomPosition(atom.GetIdx())
                positions.append([pos.x, pos.y, pos.z])
        except Exception as exc:
            QMessageBox.critical(self, PLUGIN_NAME, f"Failed to read molecule: {exc}")
            return

        # Clear previous results
        self.table.setRowCount(0)
        self.log_text.clear()
        self._step_count = 0

        if self.chk_auto_charge.isChecked():
            try:
                charge, multiplicity = derive_charge_and_multiplicity(mol)
            except Exception as exc:
                logger.warning("Could not derive charge/multiplicity: %s", exc)
                charge, multiplicity = 0, 1
            # Show what was actually used, so a wrong perception is visible
            # rather than silently baked into the result.
            self.spin_charge.setValue(charge)
            self.spin_multiplicity.setValue(multiplicity)
        else:
            charge = self.spin_charge.value()
            multiplicity = self.spin_multiplicity.value()

        method = self.combo_method.currentText()
        fmax = self.spin_fmax.value()
        max_steps = self.spin_maxsteps.value()

        self._set_running(True)
        self.lbl_status.setText(f"Running {method}…")
        QApplication.processEvents()

        # Remember what we actually optimised.  The optimisation is async, so
        # the user can load a different molecule before it finishes; the
        # results must not then be written into whatever happens to be loaded.
        self._run_mol = mol
        self._run_natoms = mol.GetNumAtoms()

        self._worker = XtbWorker(
            numbers=numbers,
            positions=positions,
            method=method,
            fmax=fmax,
            max_steps=max_steps,
            charge=charge,
            multiplicity=multiplicity,
            parent=self,
        )
        self._worker.log_message.connect(self._append_log)
        self._worker.step_update.connect(self._on_step_update)
        self._worker.finished.connect(self._on_finished)
        self._worker.start()

    def _on_cancel(self):
        if self._worker and self._worker.isRunning():
            self._worker.cancel()
            self.lbl_status.setText("Cancelling…")
            self.btn_cancel.setEnabled(False)

    def _on_step_update(self, step: int, energy: float, fmax: float):
        self._step_count = step
        row = self.table.rowCount()
        self.table.insertRow(row)
        self.table.setItem(row, 0, QTableWidgetItem(str(step)))
        self.table.setItem(row, 1, QTableWidgetItem(f"{energy:.6f}"))
        self.table.setItem(row, 2, QTableWidgetItem(f"{fmax:.4f}"))
        self.table.scrollToBottom()
        self.lbl_status.setText(
            f"Step {step}  |  E = {energy:.6f} eV  |  Fmax = {fmax:.4f} eV/Å"
        )

    def _on_finished(self, success: bool, payload):
        self._set_running(False)
        self._worker = None

        if not success:
            self.lbl_status.setText(f"Stopped: {payload}")
            return

        # Apply optimized coordinates back to the original RDKit molecule
        new_positions: list[list[float]] = payload
        mol = self.context.current_mol
        if mol is None:
            self.lbl_status.setText("Error: molecule was unloaded during optimization.")
            return

        # Conformer.SetAtomPosition accepts an out-of-range index without
        # raising, so a mismatch here would silently overwrite a different
        # molecule's coordinates and then be committed by the undo checkpoint.
        if mol is not getattr(self, "_run_mol", None) or (
            mol.GetNumAtoms() != len(new_positions)
        ):
            msg = (
                "The molecule changed while the optimization was running, so "
                "the optimized coordinates were discarded."
            )
            self.lbl_status.setText(msg)
            self._append_log("\n" + msg)
            QMessageBox.warning(self, PLUGIN_NAME, msg)
            return

        try:
            conf = mol.GetConformer()
            for i, (x, y, z) in enumerate(new_positions):
                conf.SetAtomPosition(i, (x, y, z))
        except (RuntimeError, AttributeError, TypeError, ValueError) as exc:
            logger.exception("XtbOptimizer: failed to apply optimized coordinates")
            self.lbl_status.setText(f"Error applying coordinates: {exc}")
            return

        # Push updated molecule to the main window
        self.context.current_mol = mol
        self.context.refresh_3d_view()
        self.context.push_undo_checkpoint()

        self.lbl_status.setText(
            f"Optimization complete in {self._step_count} steps. "
            "Coordinates updated. ✓"
        )
        self.context.show_status_message(
            f"{PLUGIN_NAME}: optimization complete ({self._step_count} steps).", 5000
        )
        logger.info(
            "XtbOptimizer: optimization complete (%d steps, method=%s)",
            self._step_count,
            self.combo_method.currentText(),
        )

    # ------------------------------------------------------------------
    # Dialog lifecycle
    # ------------------------------------------------------------------

    def closeEvent(self, event):
        if self._worker and self._worker.isRunning():
            reply = QMessageBox.question(
                self, PLUGIN_NAME,
                "Optimization is still running. Cancel and close?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if reply == QMessageBox.StandardButton.Yes:
                self._on_cancel()
                self._worker.wait(3000)   # give the thread 3 s to exit cleanly
                event.accept()
            else:
                event.ignore()
                return
        self.context.register_window("main_panel", None)
        event.accept()

    def accept(self):
        if self._worker and self._worker.isRunning():
            self.closeEvent(
                type("FakeEvent", (), {"accept": lambda s: None, "ignore": lambda s: None})()
            )
            return
        self.context.register_window("main_panel", None)
        super().accept()


# ---------------------------------------------------------------------------
# Plugin entry points
# ---------------------------------------------------------------------------


def run_plugin(context):
    """Open (or raise) the xTB optimizer dialog."""
    mw = context.get_main_window()

    if not context.current_mol:
        QMessageBox.warning(mw, PLUGIN_NAME, "No molecule loaded.")
        return

    win = context.get_window("main_panel")
    if win:
        win.show()
        win.raise_()
        win.activateWindow()
        return

    dialog = XtbOptimizerDialog(context, parent=mw)
    context.register_window("main_panel", dialog)
    dialog.show()


_launch_fn = None


def initialize(context):
    """Register the plugin in the 3D Edit menu."""
    global _launch_fn
    _launch_fn = lambda: run_plugin(context)
    context.add_menu_action("3D Edit/xTB Optimizer…", _launch_fn)


def run(mw):
    """Legacy entry point (Plugins menu auto-registration)."""
    if hasattr(mw, "host"):
        mw = mw.host
    if _launch_fn:
        _launch_fn()
