import numpy as np
from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QComboBox,
    QCheckBox,
    QLabel,
    QHeaderView,
    QMessageBox,
)
from PyQt6.QtCore import Qt, QTimer, QObject, QEvent
from PyQt6.QtGui import QShortcut, QKeySequence
from rdkit import Chem
from rdkit.Geometry import Point3D
import pyvista as pv
import logging
import inspect
from functools import partial


PLUGIN_NAME = "Bond Editor"
PLUGIN_VERSION = "2026.09.25"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "A table-based bond editor: add/delete bonds, change bond order, and set "
    "bond lengths by moving one side of the bond. In the 3D view pick a mode to "
    "select a bond by clicking it, or create a bond by clicking two atoms."
)
PLUGIN_CONTEXT = None

BOND_TYPE_LABELS = ["Single", "Double", "Triple", "Aromatic"]
INTERACTIVE_BOND_TYPE_LABELS = ["Single", "Double", "Triple"]
INTERACTIVE_MODE_COLOR = "#d9f7e5"
INTERACTIVE_MODE_CHECKED_COLOR = "#8ce99a"
INTERACTIVE_MODE_PRESSED_COLOR = "#69db7c"
SELECTED_BOND_COLOR = "#ff9f1c"


def bond_type_from_label(label):
    """Map a UI label to an RDKit bond type (defaults to SINGLE)."""
    return {
        "Single": Chem.BondType.SINGLE,
        "Double": Chem.BondType.DOUBLE,
        "Triple": Chem.BondType.TRIPLE,
        "Aromatic": Chem.BondType.AROMATIC,
    }.get(label, Chem.BondType.SINGLE)


def label_from_bond_type(bond_type):
    """Map an RDKit bond type to its UI label (unknown types show as Single)."""
    name = str(bond_type)
    return {
        "SINGLE": "Single",
        "DOUBLE": "Double",
        "TRIPLE": "Triple",
        "AROMATIC": "Aromatic",
    }.get(name.rsplit(".", 1)[-1], "Single")


def _has_aromatic_bonds(mol):
    """Check whether the molecule contains any aromatic bonds or atoms."""
    if not mol:
        return False
    return any(
        b.GetBondType() == Chem.BondType.AROMATIC
        or (hasattr(b, "GetIsAromatic") and b.GetIsAromatic())
        for b in mol.GetBonds()
    )


def sanitize_or_clear_aromaticity(rw):
    """Sanitize, retrying once with aromatic flags dropped if the first pass fails.

    An edit that breaks a ring leaves atoms still flagged aromatic that are no
    longer in one. UpdatePropertyCache does not clear those flags, so the
    molecule commits looking valid and only blows up later, in MolToMolBlock or
    MMFF atom typing, far from the edit that caused it.
    """
    try:
        Chem.SanitizeMol(rw)
        return
    except (RuntimeError, AttributeError, ValueError) as _e:
        logging.warning(
            "[%s] sanitize failed, clearing aromaticity: %s", PLUGIN_NAME, _e
        )
    for atom in rw.GetAtoms():
        atom.SetIsAromatic(False)
    for bond in rw.GetBonds():
        if bond.GetBondType() == Chem.BondType.AROMATIC:
            bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    try:
        Chem.SanitizeMol(rw)
    except (RuntimeError, AttributeError, ValueError) as _e:
        logging.warning("[%s] sanitize still failing: %s", PLUGIN_NAME, _e)
        try:
            rw.UpdatePropertyCache(strict=False)
            Chem.GetSSSR(rw)
        except (RuntimeError, AttributeError, ValueError) as _e2:
            logging.warning("[%s] property cache fallback: %s", PLUGIN_NAME, _e2)


class _ClickFilter(QObject):
    """Dispatch click/drag gestures and suppress camera input when editing.

    The press callback decides whether the gesture started on an editable 3D
    object. Only those gestures are consumed; empty-space drags continue to
    the VTK interactor so the camera can still rotate normally.
    """

    def __init__(self, callback, parent=None, drag_callback=None, press_callback=None):
        super().__init__(parent)
        self._callback = callback
        self._drag_callback = drag_callback
        self._press_callback = press_callback
        self._press_pos = None
        self._press_button = None
        self._consume_gesture = False

    def eventFilter(self, obj, event):
        t = event.type()
        if t == QEvent.Type.MouseButtonPress and event.button() in (
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.RightButton,
        ):
            self._press_pos = event.position().toPoint()
            self._press_button = event.button()
            self._consume_gesture = (
                bool(
                    self._press_callback(
                        self._press_pos.x(),
                        self._press_pos.y(),
                        obj,
                        event.modifiers(),
                        self._press_button,
                    )
                )
                if self._press_callback is not None
                else False
            )
            return self._consume_gesture
        if t == QEvent.Type.MouseMove and self._press_pos is not None:
            return self._consume_gesture
        if t == QEvent.Type.MouseButtonRelease and self._press_pos is not None:
            rel = event.position().toPoint()
            dx = rel.x() - self._press_pos.x()
            dy = rel.y() - self._press_pos.y()
            if dx * dx + dy * dy <= 25:
                args = (rel.x(), rel.y(), obj, event.modifiers(), self._press_button)
                try:
                    inspect.signature(self._callback).bind(*args)
                except (TypeError, ValueError):
                    self._callback(*args[:4])
                else:
                    self._callback(*args)
            elif (
                self._drag_callback is not None
                and self._press_button == Qt.MouseButton.LeftButton
            ):
                self._drag_callback(
                    rel.x(), rel.y(), obj, event.modifiers(), self._press_button
                )
            consume = self._consume_gesture
            self._press_pos = None
            self._press_button = None
            self._consume_gesture = False
            return consume
        return False


class BondEditorWindow(QWidget):
    """Table-based editor for bonds: order, existence, and length."""

    COL_IDX, COL_A1, COL_A2, COL_TYPE, COL_LEN = range(5)

    def __init__(self, context):
        super().__init__(parent=context.get_main_window())
        self.setWindowFlags(Qt.WindowType.Window)
        self.context = context
        self.setWindowTitle("Bond Editor")
        self.resize(600, 420)
        self._click_filter = None
        self._interactive_mode = True
        self._drag_start_atom = None
        self._first_pick_idx = None
        self._picked_atoms = {}
        self.init_ui()

        self.context.register_window("main_panel", self)
        self.last_seen_signature = None
        self.load_molecule()

        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.check_molecule_update)
        self.update_timer.start(500)

        self._enable_plotter_picking()

    def init_ui(self):
        layout = QVBoxLayout(self)

        self.table = QTableWidget()
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(
            ["#", "Atom 1", "Atom 2", "Type", "Length (Å)"]
        )
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.itemChanged.connect(self.on_item_changed)
        self.table.itemSelectionChanged.connect(self.highlight_selected_bonds)
        layout.addWidget(self.table)

        add_layout = QHBoxLayout()
        self.interactive_btn = QPushButton("Interactive mode")
        self.interactive_btn.setCheckable(True)
        self.interactive_btn.setStyleSheet(
            f"QPushButton {{ background-color: {INTERACTIVE_MODE_COLOR}; }} "
            f"QPushButton:checked {{ background-color: {INTERACTIVE_MODE_CHECKED_COLOR}; }} "
            f"QPushButton:pressed {{ background-color: {INTERACTIVE_MODE_PRESSED_COLOR}; }}"
        )
        self.interactive_btn.toggled.connect(self._toggle_interactive_mode)
        add_layout.addWidget(self.interactive_btn)
        self.auto_kekulize_cb = QCheckBox("Auto-kekulize")
        self.auto_kekulize_cb.setChecked(True)
        self.auto_kekulize_cb.setToolTip(
            "If aromatic bonds exist, automatically kekulize them before editing "
            "and re-perceive aromaticity afterwards."
        )
        add_layout.addWidget(self.auto_kekulize_cb)
        add_layout.addWidget(QLabel("3D click:"))
        self.click_mode_combo = QComboBox()
        self.click_mode_combo.addItems(["Select bond", "Create bond"])
        self.click_mode_combo.setToolTip(
            "Select bond: click a bond to select its row.\n"
            "Create bond: click two atoms to create a bond between them."
        )
        self.click_mode_combo.currentTextChanged.connect(self._on_click_mode_changed)
        add_layout.addWidget(self.click_mode_combo)
        add_layout.addWidget(QLabel("New bond type:"))
        self.add_type_combo = QComboBox()
        self.add_type_combo.addItems(BOND_TYPE_LABELS)
        self.add_type_combo.setToolTip(
            "Bond type used when creating a bond by clicking"
        )
        self.add_type_combo.currentTextChanged.connect(
            lambda _t: self._update_mode_overlay()
        )
        add_layout.addWidget(self.add_type_combo)
        add_layout.addStretch()
        layout.addLayout(add_layout)

        btn_layout = QHBoxLayout()
        self.adjust_h_btn = QPushButton("Adjust H")
        self.adjust_h_btn.clicked.connect(self.adjust_hydrogens)
        btn_layout.addWidget(self.adjust_h_btn)
        self.estimate_btn = QPushButton("Estimate Bonds")
        self.estimate_btn.setToolTip(
            "Estimate bonds and bond orders from 3D coordinates using RDKit"
        )
        self.estimate_btn.clicked.connect(self.estimate_bonds)
        btn_layout.addWidget(self.estimate_btn)
        self.kekulize_btn = QPushButton("Kekulize")
        self.kekulize_btn.setCheckable(True)
        self.kekulize_btn.setToolTip(
            "Toggle between Kekulized (alternating single/double) and aromatic bonds"
        )
        self.kekulize_btn.toggled.connect(self.toggle_kekulize)
        btn_layout.addWidget(self.kekulize_btn)
        self.delete_btn = QPushButton("Delete Selected Bonds")
        self.delete_btn.clicked.connect(self.delete_selected_bonds)
        btn_layout.addWidget(self.delete_btn)
        self.unselect_btn = QPushButton("Unselect")
        self.unselect_btn.setToolTip(
            "Clear the bond selection and any in-progress bond-creation pick"
        )
        self.unselect_btn.clicked.connect(self.unselect_all)
        btn_layout.addWidget(self.unselect_btn)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        hint = QLabel(
            "Edit Length to move the Atom-2 side of the bond (non-ring bonds only)."
        )
        hint.setStyleSheet("color: gray;")
        layout.addWidget(hint)

        self._del_shortcut = QShortcut(QKeySequence(Qt.Key.Key_Delete), self.table)
        self._del_shortcut.setContext(Qt.ShortcutContext.WidgetShortcut)
        self._del_shortcut.activated.connect(self.delete_selected_bonds)

        self.interactive_btn.setChecked(True)

    # ------------------------------------------------------------------
    # 3D picking (click an atom to fill the Atom 1 / Atom 2 fields)
    # ------------------------------------------------------------------

    def _enable_plotter_picking(self):
        try:
            plotter = self.context.plotter
            if plotter is None:
                return
            interactor = getattr(plotter, "interactor", None)
            if interactor is None:
                return
            self._click_filter = _ClickFilter(
                self._on_plotter_click,
                parent=self,
                drag_callback=self._on_plotter_drag,
                press_callback=self._on_plotter_press,
            )
            interactor.installEventFilter(self._click_filter)
        except (RuntimeError, AttributeError, KeyError, ValueError) as _e:
            logging.warning("[bond_editor.py:_enable_plotter_picking] silenced: %s", _e)

    def _disable_plotter_picking(self):
        try:
            plotter = self.context.plotter
            interactor = getattr(plotter, "interactor", None) if plotter else None
            if interactor and self._click_filter:
                interactor.removeEventFilter(self._click_filter)
        except (RuntimeError, AttributeError, KeyError, ValueError) as _e:
            logging.warning(
                "[bond_editor.py:_disable_plotter_picking] silenced: %s", _e
            )
        self._click_filter = None

    def _on_plotter_press(self, x, y, widget, modifiers, button):
        """Claim an interactive gesture before VTK can start camera rotation."""
        # Every gesture starts clean: a start atom left over from an earlier
        # click would turn the next camera drag into a bond.
        self._drag_start_atom = None
        if not self._interactive_mode:
            return False
        picked = self._interactive_pick(x, y, widget)
        if picked is None:
            return False
        _mol, _pos, atom, pair = picked
        if atom is not None:
            if button == Qt.MouseButton.LeftButton:
                self._drag_start_atom = atom
            return True
        return pair is not None

    def _on_plotter_click(
        self, x, y, widget, modifiers, button=Qt.MouseButton.LeftButton
    ):
        try:
            if self._interactive_mode:
                self._interactive_click(x, y, widget, button)
                return
            import vtk

            mw = self.context.get_main_window()
            v3d = getattr(mw, "view_3d_manager", None) if mw else None
            plotter = self.context.plotter
            if not v3d or not plotter:
                return

            # Qt reports click positions in logical pixels, but VTK's picker works
            # in physical device pixels. On HiDPI/Retina displays (notably macOS,
            # where devicePixelRatio == 2) the two differ, so scale by the ratio —
            # otherwise the pick lands at half the coordinates (toward the bottom-
            # left) and you have to click up and to the right of the target to hit
            # it. On Windows/Linux the ratio is 1.0, so this is a no-op there.
            ratio = widget.devicePixelRatioF()
            px = x * ratio
            vtk_y = (widget.height() - y) * ratio

            picker = vtk.vtkCellPicker()
            picker.SetTolerance(0.005)
            picker.Pick(px, vtk_y, 0, plotter.renderer)
            picked_actor = picker.GetActor()

            mol = self.context.current_mol
            if not mol or not mol.GetNumConformers():
                return

            mode = self.click_mode_combo.currentText()
            if picked_actor is None:
                # Clicked empty space: reset whatever is in progress.
                if mode == "Create bond":
                    self._cancel_bond_pick()
                else:
                    self.table.clearSelection()
                return

            pick_pos = picker.GetPickPosition()
            if mode == "Create bond":
                # Click two atoms to create the bond between them; clicking
                # anything that is not an atom cancels the pick.
                atom_actor = getattr(v3d, "atom_actor", None)
                if atom_actor is not None and picked_actor is not atom_actor:
                    self._cancel_bond_pick()
                    return
                idx = self._nearest_atom_to_point(mol, pick_pos)
                if idx is not None:
                    self._create_bond_pick(idx)
                return

            # Default (Select bond): resolve the click to the nearest bond axis,
            # whether the user clicked the bond cylinder or one of its atoms.
            pair = self._nearest_bond_to_point(mol, pick_pos)
            if pair is not None:
                self._select_bond_row_by_pair(pair)
        except Exception as _e:
            logging.warning("[bond_editor.py:_on_plotter_click] silenced: %s", _e)

    def _toggle_interactive_mode(self, enabled):
        self._interactive_mode = bool(enabled)
        self._drag_start_atom = None
        self.click_mode_combo.setEnabled(not enabled)
        # A half-finished "Create bond" pick has no meaning in interactive mode.
        self._first_pick_idx = None
        self._picked_atoms = {}
        self._update_mode_overlay()
        self._update_picked_atom_labels()
        self.context.show_status_message(
            "Interactive mode enabled." if enabled else "Interactive mode disabled."
        )

    def _interactive_pick(self, x, y, widget):
        """Resolve a viewport position to what an interactive gesture acts on.

        Returns (mol, pos, atom, pair) where at most one of *atom* (an atom
        index) and *pair* (a bonded atom pair) is set, or None when there is no
        3D molecule. Callers must act on this decision rather than re-derive
        one from *pos*: an atom's surface lies within bond-picking range of
        its own bonds, so a second lookup would edit a bond on an atom click.
        """
        import vtk

        mol = self.context.current_mol
        plotter = self.context.plotter
        if not mol or not plotter or not mol.GetNumConformers():
            return None
        picker = vtk.vtkCellPicker()
        ratio = widget.devicePixelRatioF()
        picker.SetTolerance(0.005)
        picker.Pick(x * ratio, (widget.height() - y) * ratio, 0, plotter.renderer)
        picked_actor = picker.GetActor()
        pos = picker.GetPickPosition()

        if self._pick_is_atom_actor(picked_actor):
            atom = self._screen_atom_index(x, y, widget, mol)
            if atom is None:
                atom = self._nearest_atom_to_point(mol, pos)
            return mol, pos, atom, None

        if picked_actor is not None:
            # A bond cylinder (or other scene geometry): nearest bond axis.
            return mol, pos, None, self._nearest_bond_to_point(mol, pos)

        # Nothing hit. VTK then reports the click projected onto the camera's
        # focal plane, which is close enough to catch a near miss on a thin
        # bond or small atom, with a tighter tolerance than a real hit.
        nearest_bond = self._nearest_bond_to_point(mol, pos, max_dist=0.4)
        nearest_atom = self._nearest_atom_to_point(mol, pos)
        atom_dist = (
            self._point_distance(mol, nearest_atom, pos)
            if nearest_atom is not None
            else float("inf")
        )
        bond_dist = (
            self._segment_distance(mol, nearest_bond, pos)
            if nearest_bond is not None
            else float("inf")
        )
        if nearest_bond is not None and bond_dist < atom_dist:
            return mol, pos, None, nearest_bond
        atom = self._screen_atom_index(x, y, widget, mol)
        if atom is None and atom_dist <= 0.5:
            atom = nearest_atom
        return mol, pos, atom, None

    @staticmethod
    def _point_distance(mol, idx, pos):
        p = mol.GetConformer().GetAtomPosition(idx)
        return float(
            np.linalg.norm(np.array([p.x - pos[0], p.y - pos[1], p.z - pos[2]]))
        )

    @staticmethod
    def _segment_distance(mol, pair, pos):
        """Distance from *pos* to the segment between the atoms of *pair*."""
        conf = mol.GetConformer()
        a = np.array(conf.GetAtomPosition(pair[0]))
        c = np.array(conf.GetAtomPosition(pair[1]))
        q = np.array(pos, dtype=float)
        ab = c - a
        denom = float(ab @ ab)
        t = 0.0 if denom < 1e-12 else float((q - a) @ ab / denom)
        return float(np.linalg.norm(q - (a + max(0.0, min(1.0, t)) * ab)))

    def _pick_is_atom_actor(self, picked_actor):
        """Return whether a VTK pick landed on the host's atom geometry."""
        if picked_actor is None:
            return False
        try:
            main_window = self.context.get_main_window()
            view_3d = getattr(main_window, "view_3d_manager", None)
            atom_actor = getattr(view_3d, "atom_actor", None)
            return atom_actor is not None and picked_actor is atom_actor
        except (AttributeError, RuntimeError, TypeError):
            return False

    def _screen_atom_index(self, x, y, widget, mol):
        """Use the host screen-space atom picker when available."""
        try:
            from moleditpy.ui.atom_picking import pick_atom_index_from_screen

            main_window = self.context.get_main_window()
            view_3d = getattr(main_window, "view_3d_manager", None)
            if view_3d is None:
                return None
            ratio = widget.devicePixelRatioF()
            return pick_atom_index_from_screen(
                view_3d, (int(x * ratio), int((widget.height() - y) * ratio)), mol
            )
        except (ImportError, AttributeError, RuntimeError, TypeError, ValueError):
            return None

    def _prepare_rw_for_edit(self, mol):
        rw = Chem.RWMol(mol)
        if (
            getattr(self, "auto_kekulize_cb", None)
            and self.auto_kekulize_cb.isChecked()
            and _has_aromatic_bonds(mol)
        ):
            try:
                Chem.Kekulize(rw, clearAromaticFlags=True)
            except (RuntimeError, ValueError) as e:
                logging.warning("[bond_editor] Auto-kekulize failed: %s", e)
        return rw

    def _interactive_click(self, x, y, widget, button):
        picked = self._interactive_pick(x, y, widget)
        if picked is None:
            return
        mol, _pos, _atom, pair = picked
        # Clicking an atom does nothing by itself; dragging from it makes a bond.
        if pair is None:
            return
        if button == Qt.MouseButton.RightButton:
            rw = self._prepare_rw_for_edit(mol)
            rw.RemoveBond(*pair)
            self._commit(rw, f"Deleted bond {pair[0]}-{pair[1]}.")
        else:
            self._cycle_interactive_bond_type(pair)

    def _on_plotter_drag(self, x, y, widget, modifiers, button):
        if not self._interactive_mode or self._drag_start_atom is None:
            return
        picked = self._interactive_pick(x, y, widget)
        start = self._drag_start_atom
        self._drag_start_atom = None
        if picked is not None:
            end = picked[2]
            if end is not None and end != start:
                self.add_bond(start, end)

    def _cycle_interactive_bond_type(self, pair):
        """Advance a bond through single -> double -> triple -> single.

        The current order is read after auto-kekulization: read from the
        aromatic molecule, every ring bond looks "aromatic" and always goes to
        single, which is a no-op on the Kekule single bonds.
        """
        try:
            rw = self._prepare_rw_for_edit(self.context.current_mol)
            bond = rw.GetBondBetweenAtoms(*pair)
            if bond is None:
                return
            current = label_from_bond_type(bond.GetBondType())
            cycle = INTERACTIVE_BOND_TYPE_LABELS
            index = cycle.index(current) if current in cycle else -1
            label = cycle[(index + 1) % len(cycle)]
            new_type = bond_type_from_label(label)
            bond.SetBondType(new_type)
            bond.SetIsAromatic(False)
            self._commit(rw, f"Bond {pair[0]}-{pair[1]} set to {label.lower()}.")
            self._verify_bond_type(pair, new_type, label)
        except (RuntimeError, AttributeError, ValueError) as exc:
            QMessageBox.critical(self, "Error", f"Failed to change bond type: {exc}")

    def _excess_hydrogen_indices(self, mol):
        """Indices of terminal H atoms whose removal brings over-valent heavy
        atoms back to their allowed valence, removing hydrogens with the largest
        atom ID first."""
        pt = Chem.GetPeriodicTable()
        to_remove = set()
        for atom in mol.GetAtoms():
            num = atom.GetAtomicNum()
            if num <= 1:
                continue
            try:
                allowed = pt.GetDefaultValence(num)
            except (RuntimeError, ValueError):
                continue
            if allowed <= 0:
                continue
            if num in (7, 8, 15, 16):
                allowed += atom.GetFormalCharge()
            valence = round(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()))
            excess = valence - allowed
            if excess <= 0:
                continue
            h_neighbors = sorted(
                [
                    n.GetIdx()
                    for n in atom.GetNeighbors()
                    if n.GetAtomicNum() == 1 and n.GetDegree() == 1
                ],
                reverse=True,
            )
            to_remove.update(h_neighbors[: min(excess, len(h_neighbors))])
        return sorted(to_remove)

    def adjust_hydrogens(self):
        mol = self.context.current_mol
        if not mol or not mol.GetNumConformers():
            self.context.show_status_message("No 3D molecule to adjust hydrogens on.")
            return

        try:
            rw = Chem.RWMol(mol)
            rw.UpdatePropertyCache(strict=False)

            removed = self._excess_hydrogen_indices(rw)
            for idx in sorted(removed, reverse=True):
                rw.RemoveAtom(idx)

            sanitize_or_clear_aromaticity(rw)

            new_mol = Chem.AddHs(rw, addCoords=True)
            added = new_mol.GetNumAtoms() - (mol.GetNumAtoms() - len(removed))

            if added <= 0 and not removed:
                self.context.show_status_message("Hydrogens are already explicit.")
                return

            self._set_molecule(new_mol)
            self.context.show_status_message("Hydrogens adjusted.")
        except (RuntimeError, ValueError, AttributeError) as e:
            logging.exception("[bond_editor] Failed to adjust hydrogens: %s", e)
            QMessageBox.critical(self, "Error", f"Failed to adjust hydrogens: {str(e)}")

    def estimate_from_coordinates(self):
        mol = self.context.current_molecule
        if not mol or not mol.GetNumAtoms() or not mol.GetNumConformers():
            self.context.show_status_message("No 3D molecule to estimate bonds for.")
            return

        charge = 0
        if mol.HasProp("_xyz_charge"):
            try:
                charge = int(mol.GetProp("_xyz_charge"))
            except (ValueError, TypeError):
                charge = 0
        else:
            try:
                charge = int(Chem.GetFormalCharge(mol))
            except (RuntimeError, ValueError):
                charge = 0

        applied = False
        candidate = Chem.RWMol(mol)
        try:
            from rdkit.Chem import rdDetermineBonds

            rdDetermineBonds.DetermineBonds(candidate, charge=charge)
            applied = True
        except (
            ImportError,
            RuntimeError,
            AttributeError,
            TypeError,
            ValueError,
        ) as exc:
            logging.warning("[bond_editor] rdDetermineBonds failed: %s", exc)

        if not applied:
            mw = self.context.get_main_window()
            if hasattr(mw, "io_manager") and hasattr(
                mw.io_manager, "estimate_bonds_from_distances"
            ):
                try:
                    candidate = Chem.RWMol(mol)
                    for bond in list(candidate.GetBonds()):
                        candidate.RemoveBond(
                            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                        )
                    mw.io_manager.estimate_bonds_from_distances(candidate)
                    applied = True
                except (RuntimeError, ValueError, AttributeError, TypeError) as exc:
                    logging.warning(
                        "[bond_editor] estimate_bonds_from_distances failed: %s", exc
                    )

        if not applied:
            self.context.show_status_message(
                "Failed to estimate bonds from coordinates."
            )
            return

        self._commit(candidate, "Estimated bonds from coordinates.")

    def estimate_bonds(self):
        """Estimate bonds and bond orders from 3D coordinates."""
        return self.estimate_from_coordinates()

    def toggle_kekulize(self, checked):
        """Toggle between Kekulized (alternating single/double) and aromatic bonds."""
        if checked:
            self.kekulize()
        else:
            self.aromatize()

    def kekulize(self):
        """Convert aromatic bonds into explicit alternating single and double bonds."""
        mol = self.context.current_molecule
        if not mol or not mol.GetNumAtoms():
            self.context.show_status_message("No molecule loaded.")
            self._sync_kekulize_btn(False)
            return
        rw = Chem.RWMol(mol)
        try:
            Chem.Kekulize(rw, clearAromaticFlags=True)
            self._commit(rw, "Kekulized aromatic bonds.", sanitize=False)
            self._sync_kekulize_btn(True)
        except (RuntimeError, ValueError, AttributeError) as e:
            logging.warning("[bond_editor] Kekulize failed: %s", e)
            self._sync_kekulize_btn(False)
            QMessageBox.warning(self, "Kekulize", f"Failed to kekulize molecule: {e}")

    def aromatize(self):
        """Convert alternating single and double bonds in rings back to aromatic bonds."""
        mol = self.context.current_molecule
        if not mol or not mol.GetNumAtoms():
            self.context.show_status_message("No molecule loaded.")
            self._sync_kekulize_btn(False)
            return
        rw = Chem.RWMol(mol)
        try:
            Chem.SanitizeMol(rw)
            self._commit(rw, "Aromatized bonds.", sanitize=True)
            self._sync_kekulize_btn(False)
        except (RuntimeError, ValueError, AttributeError) as e:
            logging.warning("[bond_editor] Aromatize failed: %s", e)
            self._sync_kekulize_btn(True)
            QMessageBox.warning(self, "Aromatize", f"Failed to aromatize molecule: {e}")

    def _sync_kekulize_btn(self, is_kekulized):
        if not hasattr(self, "kekulize_btn"):
            return
        self.kekulize_btn.blockSignals(True)
        self.kekulize_btn.setChecked(is_kekulized)
        self.kekulize_btn.setText("Aromatize" if is_kekulized else "Kekulize")
        self.kekulize_btn.blockSignals(False)

    def _on_click_mode_changed(self, mode):
        """Reset the two-click pick state when the 3D click mode changes."""
        self._first_pick_idx = None
        self._picked_atoms = {}
        msgs = {
            "Select bond": "Click a bond in the 3D view to select it.",
            "Create bond": "Click two atoms to create a bond between them.",
        }
        self.context.show_status_message(msgs.get(mode, ""))
        self._update_mode_overlay()
        self._update_picked_atom_labels()

    def _update_picked_atom_labels(self):
        """Label only the atom(s) picked by 3D clicks (Atom 1 / Atom 2) in the
        viewer; the label is removed when nothing is picked."""
        plotter = self.context.plotter
        if not plotter:
            return
        try:
            mol = self.context.current_mol
            points = []
            labels = []
            if mol and mol.GetNumConformers():
                conf = mol.GetConformer()
                n = mol.GetNumAtoms()
                for slot in ("Atom 1", "Atom 2"):
                    idx = self._picked_atoms.get(slot)
                    if idx is not None and 0 <= idx < n:
                        pos = conf.GetAtomPosition(idx)
                        points.append([pos.x, pos.y, pos.z])
                        labels.append(f"{slot}: {self._atom_label(mol, idx)}")
            if not points:
                plotter.remove_actor("bond_editor_atom_labels")
                plotter.render()
                return
            plotter.add_point_labels(
                points,
                labels,
                name="bond_editor_atom_labels",
                font_size=14,
                text_color="orange",
                shape=None,
                show_points=False,
                always_visible=True,
                pickable=False,
                reset_camera=False,
            )
            plotter.render()
        except Exception as _e:
            logging.warning(
                "[bond_editor.py:_update_picked_atom_labels] silenced: %s", _e
            )

    def _update_mode_overlay(self):
        """Show the 3D click mode and pick progress as a text overlay in the viewer."""
        plotter = self.context.plotter
        if not plotter:
            return
        mode = self.click_mode_combo.currentText()
        # Interactive mode ignores the click-mode combo, so its prompt would lie.
        if mode == "Create bond" and not self._interactive_mode:
            kind = self.add_type_combo.currentText().lower()
            if self._first_pick_idx is not None:
                text = (
                    f"Create {kind} bond: atom {self._first_pick_idx} picked - "
                    "click the second atom (same atom cancels)"
                )
            else:
                text = f"Create {kind} bond: click the first atom"
        else:
            text = None
        try:
            if text is None:
                plotter.remove_actor("bond_editor_mode_label")
            else:
                plotter.add_text(
                    text,
                    name="bond_editor_mode_label",
                    position="upper_left",
                    font_size=10,
                    color="orange",
                )
            plotter.render()
        except (RuntimeError, AttributeError, KeyError, ValueError) as _e:
            logging.warning("[bond_editor.py:_update_mode_overlay] silenced: %s", _e)

    def _nearest_atom_to_point(self, mol, pick_pos):
        """Return the index of the atom nearest the 3D *pick_pos*, or None."""
        conf = mol.GetConformer()
        best_idx = None
        best = float("inf")
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            dx = pos.x - pick_pos[0]
            dy = pos.y - pick_pos[1]
            dz = pos.z - pick_pos[2]
            dist = dx * dx + dy * dy + dz * dz
            if dist < best:
                best = dist
                best_idx = idx
        return best_idx

    def _cancel_bond_pick(self):
        """Drop the in-progress bond-creation pick (no-op if nothing picked)."""
        if self._first_pick_idx is None:
            return
        self._first_pick_idx = None
        self._picked_atoms = {}
        self.context.show_status_message("Bond creation pick cleared.")
        self._update_mode_overlay()
        self._update_picked_atom_labels()

    def unselect_all(self):
        """Clear the table's bond selection and any in-progress creation pick."""
        self.table.clearSelection()
        self._cancel_bond_pick()

    def _create_bond_pick(self, atom_idx):
        """Create-bond mode: first click picks the atom, second click creates
        the bond with the chosen type; clicking the same atom again cancels."""
        if self._first_pick_idx is None:
            self._first_pick_idx = atom_idx
            self._picked_atoms = {"Atom 1": atom_idx}
            self.context.show_status_message(
                f"Atom {atom_idx} picked. Click the second atom to create the bond."
            )
        elif atom_idx == self._first_pick_idx:
            self._first_pick_idx = None
            self._picked_atoms = {}
            self.context.show_status_message("Bond creation cancelled.")
        else:
            a1, self._first_pick_idx = self._first_pick_idx, None
            self._picked_atoms = {}
            self.add_bond(a1, atom_idx)
        self._update_mode_overlay()
        self._update_picked_atom_labels()

    def _nearest_bond_to_point(self, mol, pick_pos, max_dist=0.8):
        """Return the (begin, end) atom pair of the bond whose axis is closest to
        the 3D *pick_pos*, or None if no bond is within *max_dist* angstrom."""
        conf = mol.GetConformer()
        q = np.array([pick_pos[0], pick_pos[1], pick_pos[2]])
        best_pair = None
        best = max_dist * max_dist
        for bond in mol.GetBonds():
            b = bond.GetBeginAtomIdx()
            e = bond.GetEndAtomIdx()
            p1 = conf.GetAtomPosition(b)
            p2 = conf.GetAtomPosition(e)
            a = np.array([p1.x, p1.y, p1.z])
            c = np.array([p2.x, p2.y, p2.z])
            ab = c - a
            denom = float(ab @ ab)
            t = 0.0 if denom < 1e-12 else float((q - a) @ ab / denom)
            t = max(0.0, min(1.0, t))
            proj = a + t * ab
            dist = float((q - proj) @ (q - proj))
            if dist < best:
                best = dist
                best_pair = (b, e)
        return best_pair

    def _select_bond_row_by_pair(self, pair):
        """Select the single table row matching the given (begin, end) atom pair."""
        want = set(pair)
        target = None
        self.table.blockSignals(True)
        self.table.clearSelection()
        for row in range(self.table.rowCount()):
            rp = self._row_bond_atoms(row)
            if rp and set(rp) == want:
                for col in range(self.table.columnCount()):
                    item = self.table.item(row, col)
                    if item:
                        item.setSelected(True)
                target = row
                break
        self.table.blockSignals(False)
        if target is not None:
            self.table.scrollTo(self.table.model().index(target, 0))
            self.context.show_status_message(f"Selected bond {pair[0]}-{pair[1]}.")
        self.highlight_selected_bonds()

    # ------------------------------------------------------------------
    # molecule <-> table sync
    # ------------------------------------------------------------------

    def get_mol_signature(self, mol):
        if not mol:
            return None
        try:
            sig = [id(mol), mol.GetNumAtoms(), mol.GetNumBonds()]
            bond_sig = tuple(
                (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType()))
                for b in mol.GetBonds()
            )
            sig.append(hash(bond_sig))
            if mol.GetNumAtoms() > 0 and mol.GetNumConformers():
                pos_array = mol.GetConformer().GetPositions()
                sig.append(hash(np.round(pos_array, 4).tobytes()))
            return tuple(sig)
        except (RuntimeError, AttributeError, IndexError, TypeError, ValueError):
            return None

    def check_molecule_update(self):
        try:
            current_sig = self.get_mol_signature(self.context.current_molecule)
            if current_sig != self.last_seen_signature:
                self.load_molecule()
                self._update_picked_atom_labels()
        except Exception as _e:
            logging.warning("[bond_editor.py:check_molecule_update] silenced: %s", _e)

    def _atom_label(self, mol, idx):
        atom = mol.GetAtomWithIdx(idx)
        symbol = atom.GetSymbol()
        if atom.HasProp("custom_symbol"):
            symbol = atom.GetProp("custom_symbol")
        return f"{idx} ({symbol})"

    def load_molecule(self):
        self.table.blockSignals(True)
        self.table.setRowCount(0)

        mol = self.context.current_molecule
        self.last_seen_signature = self.get_mol_signature(mol)

        if not mol or not mol.GetNumAtoms():
            self.table.blockSignals(False)
            return

        conf = mol.GetConformer() if mol.GetNumConformers() else None
        for bond in mol.GetBonds():
            row = self.table.rowCount()
            self.table.insertRow(row)
            b = bond.GetBeginAtomIdx()
            e = bond.GetEndAtomIdx()

            for col, text in (
                (self.COL_IDX, str(bond.GetIdx())),
                (self.COL_A1, self._atom_label(mol, b)),
                (self.COL_A2, self._atom_label(mol, e)),
            ):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row, col, item)

            combo = QComboBox()
            combo.addItems(BOND_TYPE_LABELS)
            combo.setCurrentText(label_from_bond_type(bond.GetBondType()))
            combo.currentTextChanged.connect(partial(self.on_type_changed, row))
            self.table.setCellWidget(row, self.COL_TYPE, combo)

            if conf is not None:
                p1 = conf.GetAtomPosition(b)
                p2 = conf.GetAtomPosition(e)
                length = (
                    (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2
                ) ** 0.5
                self.table.setItem(row, self.COL_LEN, QTableWidgetItem(f"{length:.4f}"))
            else:
                item = QTableWidgetItem("n/a")
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row, self.COL_LEN, item)

        if _has_aromatic_bonds(mol):
            self._sync_kekulize_btn(False)

        self.table.blockSignals(False)

    def _row_bond_atoms(self, row):
        """(begin, end) atom indices for a table row, parsed from the labels."""
        try:
            a1 = int(self.table.item(row, self.COL_A1).text().split()[0])
            a2 = int(self.table.item(row, self.COL_A2).text().split()[0])
            return (a1, a2)
        except (AttributeError, ValueError, IndexError):
            return None

    # ------------------------------------------------------------------
    # edit operations (each commits immediately with an undo checkpoint)
    # ------------------------------------------------------------------

    def _commit(self, rw, message, sanitize=None):
        if sanitize is None:
            sanitize = not (
                hasattr(self, "kekulize_btn") and self.kekulize_btn.isChecked()
            )
        if sanitize:
            sanitize_or_clear_aromaticity(rw)
        else:
            try:
                Chem.SanitizeMol(
                    rw,
                    sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
                )
            except (RuntimeError, ValueError):
                rw.UpdatePropertyCache(strict=False)
        self._set_molecule(rw.GetMol())
        self.context.show_status_message(message)

    def _set_molecule(self, new_mol):
        """Hand an edited molecule to the host and record an undo step.

        The host's current_mol setter (current_molecule is an alias) stores
        the molecule and redraws the 3D view. Setting both names and then
        calling refresh_3d_view as well redrew the whole scene three times per
        edit; the reset_3d_camera fallback also threw away the user's view.
        """
        self.context.current_molecule = new_mol
        self.context.push_undo_checkpoint()
        self.last_seen_signature = self.get_mol_signature(new_mol)
        self.load_molecule()

    def add_bond(self, a1, a2):
        mol = self.context.current_molecule
        if not mol or not mol.GetNumAtoms():
            self.context.show_status_message("No molecule loaded.")
            return
        if a1 == a2:
            self.context.show_status_message("Cannot bond an atom to itself.")
            return
        if mol.GetBondBetweenAtoms(a1, a2) is not None:
            self.context.show_status_message(
                f"Bond {a1}-{a2} already exists — edit its type in the table."
            )
            return
        try:
            rw = (
                self._prepare_rw_for_edit(mol)
                if hasattr(self, "_prepare_rw_for_edit")
                else Chem.RWMol(mol)
            )
            rw.AddBond(a1, a2, bond_type_from_label(self.add_type_combo.currentText()))
            self._commit(
                rw, f"Added {self.add_type_combo.currentText().lower()} bond {a1}-{a2}."
            )
        except (RuntimeError, AttributeError, ValueError) as e:
            QMessageBox.critical(self, "Error", f"Failed to add bond: {str(e)}")

    def delete_selected_bonds(self):
        mol = self.context.current_molecule
        rows = sorted(set(index.row() for index in self.table.selectedIndexes()))
        if not mol or not rows:
            self.context.show_status_message("No bonds selected to delete.")
            return
        pairs = [self._row_bond_atoms(row) for row in rows]
        pairs = [p for p in pairs if p]
        if not pairs:
            return
        try:
            rw = (
                self._prepare_rw_for_edit(mol)
                if hasattr(self, "_prepare_rw_for_edit")
                else Chem.RWMol(mol)
            )
            for a1, a2 in pairs:
                rw.RemoveBond(a1, a2)
            self._commit(rw, f"Deleted {len(pairs)} bond(s).")
        except (RuntimeError, AttributeError, ValueError) as e:
            QMessageBox.critical(self, "Error", f"Failed to delete bonds: {str(e)}")

    def _verify_bond_type(self, pair, requested, label):
        """Warn when sanitization restored the previous bond type.

        SanitizeMol re-perceives aromaticity, so setting a benzene or pyridine
        ring bond to single/double silently comes back as aromatic — the edit
        did nothing while the status bar reported success.
        """
        mol = self.context.current_molecule
        if not mol:
            return True
        try:
            bond = mol.GetBondBetweenAtoms(int(pair[0]), int(pair[1]))
        except (IndexError, TypeError, ValueError) as _e:
            logging.warning("[bond_editor.py:_verify_bond_type] %s", _e)
            return True
        if bond is None or bond.GetBondType() == requested:
            return True
        actual = bond.GetBondType()
        # Only a nameable bond type is evidence of a mismatch. A stand-in
        # object compares unequal to anything and would raise a modal dialog
        # that blocks forever in a headless run.
        if str(actual).rsplit(".", 1)[-1] not in (
            "SINGLE",
            "DOUBLE",
            "TRIPLE",
            "AROMATIC",
        ):
            return True
        extra = ""
        if actual == Chem.BondType.AROMATIC:
            extra = (
                "\n\nThe bond belongs to an aromatic ring. Break the ring's "
                "aromaticity first (for example by changing an atom or adding "
                "a hydrogen), then set the bond type."
            )
        QMessageBox.warning(
            self,
            PLUGIN_NAME,
            f"Bond {pair[0]}-{pair[1]} could not be set to {label.lower()} — "
            f"it is {str(actual).split('.')[-1].lower()}.{extra}",
        )
        return False

    def on_type_changed(self, row, label):
        mol = self.context.current_molecule
        pair = self._row_bond_atoms(row)
        if not mol or not pair:
            return
        try:
            rw = Chem.RWMol(mol)
            bond = rw.GetBondBetweenAtoms(pair[0], pair[1])
            if bond is None:
                return
            new_type = bond_type_from_label(label)
            bond.SetBondType(new_type)
            aromatic = new_type == Chem.BondType.AROMATIC
            bond.SetIsAromatic(aromatic)
            if aromatic:
                bond.GetBeginAtom().SetIsAromatic(True)
                bond.GetEndAtom().SetIsAromatic(True)
            self._commit(rw, f"Bond {pair[0]}-{pair[1]} set to {label.lower()}.")
            self._verify_bond_type(pair, new_type, label)
        except (RuntimeError, AttributeError, IndexError, ValueError) as e:
            QMessageBox.critical(self, "Error", f"Failed to change bond type: {str(e)}")

    # ------------------------------------------------------------------
    # bond length editing
    # ------------------------------------------------------------------

    def _moving_side(self, mol, begin_idx, end_idx):
        """Atom indices reachable from end_idx without crossing the
        begin-end bond, or None if the bond is part of a ring."""
        visited = {end_idx}
        queue = [end_idx]
        while queue:
            cur = queue.pop()
            for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni = nb.GetIdx()
                if cur == end_idx and ni == begin_idx:
                    continue
                if ni == begin_idx:
                    return None
                if ni not in visited:
                    visited.add(ni)
                    queue.append(ni)
        return visited

    def on_item_changed(self, item):
        if item.column() != self.COL_LEN:
            return
        try:
            target = float(item.text())
        except ValueError:
            self.load_molecule()
            return
        if target <= 0:
            self.context.show_status_message("Bond length must be positive.")
            self.load_molecule()
            return
        self.set_bond_length(item.row(), target)

    def set_bond_length(self, row, target):
        """Set the bond length by translating the Atom-2-side fragment
        along the bond axis. Ring bonds are refused."""
        mol = self.context.current_molecule
        pair = self._row_bond_atoms(row)
        if not mol or not pair or not mol.GetNumConformers():
            return
        begin_idx, end_idx = pair

        side = self._moving_side(mol, begin_idx, end_idx)
        if side is None:
            self.context.show_status_message(
                f"Bond {begin_idx}-{end_idx} is in a ring — length not editable."
            )
            self.load_molecule()
            return

        conf = mol.GetConformer()
        p1 = conf.GetAtomPosition(begin_idx)
        p2 = conf.GetAtomPosition(end_idx)
        vec = np.array([p2.x - p1.x, p2.y - p1.y, p2.z - p1.z])
        dist = float(np.linalg.norm(vec))
        if dist < 1e-6:
            self.context.show_status_message(
                "Atoms are coincident — cannot set bond length."
            )
            self.load_molecule()
            return

        delta = vec * (target / dist - 1.0)
        try:
            rw = Chem.RWMol(mol)
            new_conf = rw.GetConformer()
            for idx in side:
                p = new_conf.GetAtomPosition(idx)
                new_conf.SetAtomPosition(
                    idx, Point3D(p.x + delta[0], p.y + delta[1], p.z + delta[2])
                )
            self._commit(
                rw, f"Bond {begin_idx}-{end_idx} length set to {target:.4f} Å."
            )
        except (RuntimeError, AttributeError, IndexError, ValueError) as e:
            QMessageBox.critical(self, "Error", f"Failed to set bond length: {str(e)}")

    # ------------------------------------------------------------------
    # 3D highlight of selected bonds
    # ------------------------------------------------------------------

    def highlight_selected_bonds(self):
        plotter = self.context.plotter
        if not plotter:
            return
        try:
            cam = plotter.camera_position
        except (AttributeError, RuntimeError, TypeError):
            cam = None

        rows = sorted(set(index.row() for index in self.table.selectedIndexes()))
        mol = self.context.current_molecule
        segments = []
        if mol and mol.GetNumConformers():
            conf = mol.GetConformer()
            for row in rows:
                pair = self._row_bond_atoms(row)
                if not pair:
                    continue
                p1 = conf.GetAtomPosition(pair[0])
                p2 = conf.GetAtomPosition(pair[1])
                segments.append(([p1.x, p1.y, p1.z], [p2.x, p2.y, p2.z]))

        if not segments:
            plotter.remove_actor("bond_editor_selection")
        else:
            points = []
            lines = []
            for i, (a, b) in enumerate(segments):
                points.extend([a, b])
                lines.extend([2, 2 * i, 2 * i + 1])
            poly = pv.PolyData(np.array(points), lines=np.array(lines))
            tube = poly.tube(radius=0.15)
            plotter.add_mesh(
                tube,
                name="bond_editor_selection",
                color=SELECTED_BOND_COLOR,
                opacity=0.6,
                pickable=False,
                reset_camera=False,
            )

        if cam is not None:
            try:
                plotter.camera_position = cam
            except (AttributeError, RuntimeError, TypeError):
                pass
        plotter.render()

    def closeEvent(self, event):
        try:
            if self.update_timer.isActive():
                self.update_timer.stop()
        except Exception as _e:
            logging.warning("[bond_editor.py:closeEvent] stop timer silenced: %s", _e)
        self._disable_plotter_picking()
        plotter = self.context.plotter
        if plotter:
            plotter.remove_actor("bond_editor_selection")
            plotter.remove_actor("bond_editor_mode_label")
            plotter.remove_actor("bond_editor_atom_labels")
            plotter.render()
        super().closeEvent(event)
        # Unregister so the next open builds a fresh window that reinstalls the
        # plotter click filter; reusing the closed one left 3D bonds unclickable.
        self.context.register_window("main_panel", None)


def initialize(context):
    """MoleditPy Plugin Entry Point (V4)"""
    global PLUGIN_CONTEXT
    PLUGIN_CONTEXT = context

    def show_editor():
        win = context.get_window("main_panel")
        if win:
            win.show()
            win.raise_()
            win.activateWindow()
            win.load_molecule()
            return
        win = BondEditorWindow(context)
        win.show()

    context.add_menu_action("3D Edit/Bond Editor...", show_editor)

    def on_document_reset():
        win = context.get_window("main_panel")
        if win:
            win.load_molecule()

    context.register_document_reset_handler(on_document_reset)
