import numpy as np
from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QComboBox,
    QSpinBox,
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
from functools import partial


PLUGIN_NAME = "Bond Editor"
PLUGIN_VERSION = "2026.07.16"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "A table-based bond editor: add/delete bonds, change bond order, and set "
    "bond lengths by moving one side of the bond. Click atoms in the 3D view "
    "to pick bond endpoints, or click a bond directly to select it."
)
PLUGIN_CONTEXT = None

BOND_TYPE_LABELS = ["Single", "Double", "Triple", "Aromatic"]


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


class _ClickFilter(QObject):
    """Qt event filter: detects non-drag left clicks on the 3D plotter widget."""

    def __init__(self, callback, parent=None):
        super().__init__(parent)
        self._callback = callback
        self._press_pos = None

    def eventFilter(self, obj, event):
        t = event.type()
        if t == QEvent.Type.MouseButtonPress:
            if event.button() == Qt.MouseButton.LeftButton:
                self._press_pos = event.position().toPoint()
        elif t == QEvent.Type.MouseButtonRelease:
            if (
                event.button() == Qt.MouseButton.LeftButton
                and self._press_pos is not None
            ):
                rel = event.position().toPoint()
                dx = rel.x() - self._press_pos.x()
                dy = rel.y() - self._press_pos.y()
                if dx * dx + dy * dy <= 25:  # ≤5 px → click, not drag
                    self._callback(rel.x(), rel.y(), obj, event.modifiers())
                self._press_pos = None
        return False  # never consume — camera interaction still works


class BondEditorWindow(QWidget):
    """Table-based editor for bonds: order, existence, and length."""

    COL_IDX, COL_A1, COL_A2, COL_TYPE, COL_LEN = range(5)

    def __init__(self, context):
        super().__init__(parent=context.get_main_window())
        self.setWindowFlags(Qt.WindowType.Window)
        self.context = context
        self.setWindowTitle("Bond Editor")
        self.resize(560, 420)
        self._click_filter = None
        self._pick_second = False
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
        add_layout.addWidget(QLabel("Atom 1:"))
        self.spin_a1 = QSpinBox()
        self.spin_a1.setRange(0, 0)
        add_layout.addWidget(self.spin_a1)
        add_layout.addWidget(QLabel("Atom 2:"))
        self.spin_a2 = QSpinBox()
        self.spin_a2.setRange(0, 0)
        add_layout.addWidget(self.spin_a2)
        self.add_type_combo = QComboBox()
        self.add_type_combo.addItems(BOND_TYPE_LABELS)
        add_layout.addWidget(self.add_type_combo)
        self.add_btn = QPushButton("Add Bond")
        self.add_btn.setToolTip(
            "Add a bond between Atom 1 and Atom 2 "
            "(click atoms in the 3D view to fill the fields)"
        )
        self.add_btn.clicked.connect(self.add_bond)
        add_layout.addWidget(self.add_btn)
        add_layout.addStretch()
        self.endpoint_mode_cb = QCheckBox("Endpoint pick mode")
        self.endpoint_mode_cb.setToolTip(
            "When on, clicking atoms in the 3D view fills Atom 1 / Atom 2 "
            "(for Add Bond). When off, clicking selects the nearest bond."
        )
        self.endpoint_mode_cb.toggled.connect(self._on_endpoint_mode_toggled)
        add_layout.addWidget(self.endpoint_mode_cb)
        layout.addLayout(add_layout)

        btn_layout = QHBoxLayout()
        self.delete_btn = QPushButton("Delete Selected Bonds")
        self.delete_btn.clicked.connect(self.delete_selected_bonds)
        btn_layout.addWidget(self.delete_btn)
        self.unselect_btn = QPushButton("Unselect")
        self.unselect_btn.clicked.connect(self.table.clearSelection)
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
            self._click_filter = _ClickFilter(self._on_plotter_click, parent=self)
            interactor.installEventFilter(self._click_filter)
        except Exception as _e:
            logging.warning("[bond_editor.py:_enable_plotter_picking] silenced: %s", _e)

    def _disable_plotter_picking(self):
        try:
            plotter = self.context.plotter
            interactor = getattr(plotter, "interactor", None) if plotter else None
            if interactor and self._click_filter:
                interactor.removeEventFilter(self._click_filter)
        except Exception as _e:
            logging.warning("[bond_editor.py:_disable_plotter_picking] silenced: %s", _e)
        self._click_filter = None

    def _on_plotter_click(self, x, y, widget, modifiers):
        try:
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
            if picker.GetActor() is None:
                return

            mol = self.context.current_mol
            if not mol or not mol.GetNumConformers():
                return

            pick_pos = picker.GetPickPosition()
            if self.endpoint_mode_cb.isChecked():
                # Endpoint mode: click atoms to fill the Atom 1 / Atom 2 pickers.
                idx = self._nearest_atom_to_point(mol, pick_pos)
                if idx is not None:
                    self._assign_picked_atom(idx)
                return

            # Default: resolve the click to the nearest bond axis, whether the
            # user clicked the bond cylinder or one of its atoms.
            pair = self._nearest_bond_to_point(mol, pick_pos)
            if pair is not None:
                self._select_bond_row_by_pair(pair)
        except Exception as _e:
            logging.warning("[bond_editor.py:_on_plotter_click] silenced: %s", _e)

    def _on_endpoint_mode_toggled(self, checked):
        """Reset the Atom 1/Atom 2 alternation when entering endpoint mode."""
        self._pick_second = False
        if checked:
            self.context.show_status_message(
                "Endpoint mode: click two atoms to set Atom 1 and Atom 2."
            )
        else:
            self.context.show_status_message("Bond selection mode: click a bond.")

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

    def _assign_picked_atom(self, atom_idx):
        """Alternate clicked atoms into the Atom 1 / Atom 2 pickers."""
        spin = self.spin_a2 if self._pick_second else self.spin_a1
        spin.setValue(atom_idx)
        which = "Atom 2" if self._pick_second else "Atom 1"
        self._pick_second = not self._pick_second
        self.context.show_status_message(f"{which} set to atom {atom_idx}.")

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
        except Exception:
            return None

    def check_molecule_update(self):
        try:
            current_sig = self.get_mol_signature(self.context.current_molecule)
            if current_sig != self.last_seen_signature:
                self.load_molecule()
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
            self.spin_a1.setRange(0, 0)
            self.spin_a2.setRange(0, 0)
            self.table.blockSignals(False)
            return

        n = mol.GetNumAtoms()
        self.spin_a1.setRange(0, n - 1)
        self.spin_a2.setRange(0, n - 1)

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
                self.table.setItem(
                    row, self.COL_LEN, QTableWidgetItem(f"{length:.4f}")
                )
            else:
                item = QTableWidgetItem("n/a")
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row, self.COL_LEN, item)

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

    def _commit(self, rw, message):
        try:
            Chem.SanitizeMol(rw)
        except Exception:
            rw.UpdatePropertyCache(strict=False)
            Chem.GetSSSR(rw)
        self.context.current_molecule = rw.GetMol()
        self.context.push_undo_checkpoint()
        self.last_seen_signature = self.get_mol_signature(
            self.context.current_molecule
        )
        refresh = getattr(self.context, "refresh_3d_view", None)
        if callable(refresh):
            refresh()
        else:
            self.context.reset_3d_camera()
        self.load_molecule()
        self.context.show_status_message(message)

    def add_bond(self):
        mol = self.context.current_molecule
        if not mol or not mol.GetNumAtoms():
            self.context.show_status_message("No molecule loaded.")
            return
        a1 = self.spin_a1.value()
        a2 = self.spin_a2.value()
        if a1 == a2:
            self.context.show_status_message("Cannot bond an atom to itself.")
            return
        if mol.GetBondBetweenAtoms(a1, a2) is not None:
            self.context.show_status_message(
                f"Bond {a1}-{a2} already exists — edit its type in the table."
            )
            return
        try:
            rw = Chem.RWMol(mol)
            rw.AddBond(a1, a2, bond_type_from_label(self.add_type_combo.currentText()))
            self._commit(rw, f"Added {self.add_type_combo.currentText().lower()} bond {a1}-{a2}.")
        except Exception as e:
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
            rw = Chem.RWMol(mol)
            for a1, a2 in pairs:
                rw.RemoveBond(a1, a2)
            self._commit(rw, f"Deleted {len(pairs)} bond(s).")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to delete bonds: {str(e)}")

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
        except Exception as e:
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
        except Exception as e:
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
                color="yellow",
                opacity=0.6,
                pickable=False,
            )

        if cam is not None:
            try:
                plotter.camera_position = cam
            except (AttributeError, RuntimeError, TypeError):
                pass
        plotter.render()

    def closeEvent(self, event):
        self._disable_plotter_picking()
        plotter = self.context.plotter
        if plotter:
            plotter.remove_actor("bond_editor_selection")
            plotter.render()
        super().closeEvent(event)


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
