from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QMessageBox,
    QLabel,
    QComboBox,
    QSpinBox,
    QDoubleSpinBox,
)
from PyQt6.QtCore import QTimer
from rdkit.Chem import AllChem
from rdkit import Geometry
import logging
import math

PLUGIN_NAME = "Step Optimizer"
PLUGIN_VERSION = "2026.07.12"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Interactive step-by-step force-field optimization with live 3D animation."
)


class StepOptimizerDialog(QDialog):
    def __init__(self, context, parent=None):
        super().__init__(parent)
        self.context = context
        self.main_window = context.get_main_window()
        self.setWindowTitle("Step Optimizer")
        self.resize(360, 260)

        # The molecule targeted by the currently active run (None between runs).
        self.target_mol = None
        # Original coordinates captured at the most recent Start, for revert.
        self.original_coords = []

        self.ff = None
        self.step_count = 0
        self.running = False
        self.steps_since_checkpoint = 0

        self.timer = QTimer(self)
        self.timer.setInterval(50)
        self.timer.timeout.connect(self._tick)

        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # Force Field Selection
        hbox_ff = QHBoxLayout()
        hbox_ff.addWidget(QLabel("Force Field:"))
        self.combo_ff = QComboBox()
        self.combo_ff.addItems(["MMFF94", "MMFF94s", "UFF"])
        hbox_ff.addWidget(self.combo_ff)
        hbox_ff.addStretch()
        layout.addLayout(hbox_ff)

        # Set default based on main window setting
        default_method = self.context.get_setting("optimization_method", "MMFF_RDKIT")
        if default_method:
            default_method = default_method.upper()
            if "UFF" in default_method:
                self.combo_ff.setCurrentText("UFF")
            else:
                self.combo_ff.setCurrentText("MMFF94")

        # Steps per frame
        hbox_steps = QHBoxLayout()
        hbox_steps.addWidget(QLabel("Steps per frame:"))
        self.spin_steps = QSpinBox()
        self.spin_steps.setRange(1, 200)
        self.spin_steps.setValue(1)
        hbox_steps.addWidget(self.spin_steps)
        hbox_steps.addStretch()
        layout.addLayout(hbox_steps)

        # Max move per cycle (clamp on per-atom displacement per tick)
        hbox_move = QHBoxLayout()
        hbox_move.addWidget(QLabel("Max move/cycle (Å):"))
        self.spin_max_move = QDoubleSpinBox()
        self.spin_max_move.setRange(0.00, 10.00)
        self.spin_max_move.setDecimals(2)
        self.spin_max_move.setSingleStep(0.05)
        self.spin_max_move.setValue(0.50)
        self.spin_max_move.setSpecialValueText("Off")
        hbox_move.addWidget(self.spin_max_move)
        hbox_move.addStretch()
        layout.addLayout(hbox_move)

        # Status labels
        self.lbl_step = QLabel("Step: 0")
        layout.addWidget(self.lbl_step)

        self.lbl_energy = QLabel("Energy: --")
        layout.addWidget(self.lbl_energy)

        self.lbl_state = QLabel("State: Idle")
        layout.addWidget(self.lbl_state)

        # Buttons
        btn_layout = QHBoxLayout()
        self.btn_start = QPushButton("Start")
        self.btn_start.clicked.connect(self._start)

        self.btn_stop = QPushButton("Stop")
        self.btn_stop.clicked.connect(self._stop)
        self.btn_stop.setEnabled(False)

        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.accept)

        btn_layout.addWidget(self.btn_start)
        btn_layout.addWidget(self.btn_stop)
        btn_layout.addWidget(self.btn_close)
        layout.addLayout(btn_layout)

    def _start(self):
        mol = self.context.current_mol
        if not mol:
            QMessageBox.warning(self, PLUGIN_NAME, "No molecule loaded.")
            return
        if mol.GetNumConformers() == 0:
            QMessageBox.warning(self, PLUGIN_NAME, "Molecule has no 3D coordinates.")
            return

        try:
            conf = mol.GetConformer()
            self.original_coords = [
                conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())
            ]

            selected_ff = self.combo_ff.currentText()
            if selected_ff in ("MMFF94", "MMFF94s"):
                props = AllChem.MMFFGetMoleculeProperties(
                    mol, mmffVariant=selected_ff
                )
                if not props:
                    QMessageBox.warning(
                        self,
                        PLUGIN_NAME,
                        "MMFF parameters unavailable for this molecule.",
                    )
                    return
                ff = AllChem.MMFFGetMoleculeForceField(mol, props)
            else:
                ff = AllChem.UFFGetMoleculeForceField(mol)

            if not ff:
                QMessageBox.warning(
                    self, PLUGIN_NAME, "Failed to build force field for this molecule."
                )
                return
            ff.Initialize()

            self.target_mol = mol
            self.ff = ff
            self.step_count = 0
            self.steps_since_checkpoint = 0

            self.btn_start.setEnabled(False)
            self.btn_stop.setEnabled(True)
            self.combo_ff.setEnabled(False)
            self.lbl_state.setText("State: Running")

            self.running = True
            self.timer.start()
        except Exception:
            logging.exception("Step Optimizer: failed to start optimization")
            QMessageBox.critical(
                self, PLUGIN_NAME, "Failed to start optimization. See log for details."
            )

    def _tick(self):
        try:
            if self.context.current_mol is not self.target_mol:
                self._pause_timer()
                self.lbl_state.setText("State: Molecule changed — run stopped")
                self._reenable_controls()
                return

            n_steps = self.spin_steps.value()
            max_move = self.spin_max_move.value()

            conf = self.target_mol.GetConformer()
            snapshot = None
            if max_move > 0:
                snapshot = [
                    (lambda p: (p.x, p.y, p.z))(conf.GetAtomPosition(i))
                    for i in range(self.target_mol.GetNumAtoms())
                ]

            res = self.ff.Minimize(maxIts=n_steps)

            if snapshot is not None:
                for i, (ox, oy, oz) in enumerate(snapshot):
                    p = conf.GetAtomPosition(i)
                    dx, dy, dz = p.x - ox, p.y - oy, p.z - oz
                    dist = math.sqrt(dx * dx + dy * dy + dz * dz)
                    if dist > max_move:
                        scale = max_move / dist
                        conf.SetAtomPosition(
                            i,
                            Geometry.Point3D(
                                ox + dx * scale, oy + dy * scale, oz + dz * scale
                            ),
                        )

            self.step_count += n_steps
            self.steps_since_checkpoint += n_steps

            self.lbl_step.setText(f"Step: {self.step_count}")
            energy = self.ff.CalcEnergy()
            self.lbl_energy.setText(f"Energy: {energy:.4f} kcal/mol")
            self.context.refresh_3d_view()

            # `res` (Minimize's return code) is intentionally not surfaced in
            # the UI: the run continues regardless of convergence so the user
            # can keep pulling atoms and watch the structure re-relax. The
            # state label stays "State: Running" for the whole run.
        except Exception:
            logging.exception("Step Optimizer: error during optimization tick")
            self._pause_timer()
            self.lbl_state.setText("State: Error — see log")
            self._reenable_controls()

    def _stop(self):
        self._pause_timer()
        self.lbl_state.setText("State: Stopped")
        if self.target_mol:
            self.context.push_undo_checkpoint()
            self.steps_since_checkpoint = 0
        self._reenable_controls()

    def _pause_timer(self):
        self.timer.stop()
        self.running = False

    def _reenable_controls(self):
        self.btn_start.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.combo_ff.setEnabled(True)

    def accept(self):
        self._pause_timer()
        if self.target_mol and self.steps_since_checkpoint > 0:
            self.context.push_undo_checkpoint()
            self.steps_since_checkpoint = 0
        super().accept()
        self.context.register_window("main_panel", None)

    def reject(self):
        self._pause_timer()
        if self.target_mol and self.original_coords:
            conf = self.target_mol.GetConformer()
            for i, pos in enumerate(self.original_coords):
                conf.SetAtomPosition(i, pos)
            self.context.refresh_3d_view()
        super().reject()
        self.context.register_window("main_panel", None)

    def closeEvent(self, event):
        # X button: behave like the Close button (accept), not cancel
        self.accept()
        event.ignore()


def run_plugin(context):
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

    dialog = StepOptimizerDialog(context, parent=mw)
    context.register_window("main_panel", dialog)
    dialog.show()


_launch_fn = None


def initialize(context):
    """Register the plugin in the 3D Edit menu."""
    global _launch_fn
    _launch_fn = lambda: run_plugin(context)
    context.add_menu_action("3D Edit/Step Optimizer...", _launch_fn)


def run(mw):
    if hasattr(mw, "host"):
        mw = mw.host
    if _launch_fn:
        _launch_fn()
