import numpy as np
from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QSpinBox,
    QLabel,
    QHeaderView,
    QMessageBox,
)
from PyQt6.QtCore import Qt, QTimer, QObject, QEvent
from rdkit import Chem
import pyvista as pv
import logging
from functools import partial


PLUGIN_NAME = "Charge Editor"
PLUGIN_VERSION = "2026.07.16"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "A table-based editor for per-atom formal charge and radical electrons. "
    "Click atoms in the 3D view to select them, and read off the total charge "
    "and spin multiplicity for DFT input preparation."
)
PLUGIN_CONTEXT = None

CHARGE_MIN, CHARGE_MAX = -4, 4
RADICAL_MIN, RADICAL_MAX = 0, 4


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


class ChargeEditorWindow(QWidget):
    """Table-based editor for per-atom formal charge and radical electrons."""

    COL_IDX, COL_SYMBOL, COL_CHARGE, COL_RADICAL, COL_HS = range(5)

    def __init__(self, context):
        super().__init__(parent=context.get_main_window())
        self.setWindowFlags(Qt.WindowType.Window)
        self.context = context
        self.setWindowTitle("Charge Editor")
        self.resize(560, 440)
        self._click_filter = None
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
            ["#", "Atom", "Formal Charge", "Radical e⁻", "H count"]
        )
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.itemSelectionChanged.connect(self.highlight_selected_atoms)
        layout.addWidget(self.table)

        self.summary_label = QLabel("Total charge: 0    Multiplicity: 1")
        self.summary_label.setStyleSheet("font-weight: bold;")
        layout.addWidget(self.summary_label)

        btn_layout = QHBoxLayout()
        self.clear_charges_btn = QPushButton("Clear All Charges")
        self.clear_charges_btn.clicked.connect(self.clear_all_charges)
        btn_layout.addWidget(self.clear_charges_btn)
        self.clear_radicals_btn = QPushButton("Clear All Radicals")
        self.clear_radicals_btn.clicked.connect(self.clear_all_radicals)
        btn_layout.addWidget(self.clear_radicals_btn)
        self.unselect_btn = QPushButton("Unselect")
        self.unselect_btn.clicked.connect(self.table.clearSelection)
        btn_layout.addWidget(self.unselect_btn)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        hint = QLabel(
            "Edit a row's Formal Charge / Radical e⁻ directly, or click an "
            "atom in the 3D view to select its row."
        )
        hint.setStyleSheet("color: gray;")
        layout.addWidget(hint)

    # ------------------------------------------------------------------
    # 3D picking (click an atom to select its row)
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
            logging.warning(
                "[charge_editor.py:_enable_plotter_picking] silenced: %s", _e
            )

    def _disable_plotter_picking(self):
        try:
            plotter = self.context.plotter
            interactor = getattr(plotter, "interactor", None) if plotter else None
            if interactor and self._click_filter:
                interactor.removeEventFilter(self._click_filter)
        except Exception as _e:
            logging.warning(
                "[charge_editor.py:_disable_plotter_picking] silenced: %s", _e
            )
        self._click_filter = None

    def _on_plotter_click(self, x, y, widget, modifiers):
        try:
            import vtk

            plotter = self.context.plotter
            if not plotter:
                return

            # Qt reports click positions in logical pixels, but VTK's picker works
            # in physical device pixels. On HiDPI/Retina displays (notably macOS,
            # where devicePixelRatio == 2) the two differ, so scale by the ratio —
            # otherwise the pick lands at half the coordinates (toward the bottom-
            # left) and you have to click up and to the right of the atom to hit
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

            idx = self._nearest_atom_to_point(mol, picker.GetPickPosition())
            if idx is None:
                return

            ctrl_held = bool(modifiers & Qt.KeyboardModifier.ControlModifier)
            self._select_atom_row(idx, ctrl_held)
        except Exception as _e:
            logging.warning(
                "[charge_editor.py:_on_plotter_click] silenced: %s", _e
            )

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

    def _select_atom_row(self, atom_idx, ctrl_held):
        """Select the table row for *atom_idx* (Ctrl adds to the selection)."""
        target = None
        self.table.blockSignals(True)
        if not ctrl_held:
            self.table.clearSelection()
        for row in range(self.table.rowCount()):
            if self._row_atom_index(row) == atom_idx:
                already = self.table.item(row, self.COL_IDX)
                select = not (ctrl_held and already and already.isSelected())
                for col in (self.COL_IDX, self.COL_SYMBOL, self.COL_HS):
                    item = self.table.item(row, col)
                    if item:
                        item.setSelected(select)
                target = row
                break
        self.table.blockSignals(False)
        if target is not None:
            self.table.scrollTo(self.table.model().index(target, 0))
        self.highlight_selected_atoms()

    # ------------------------------------------------------------------
    # molecule <-> table sync
    # ------------------------------------------------------------------

    def get_mol_signature(self, mol):
        if not mol:
            return None
        try:
            sig = [id(mol), mol.GetNumAtoms()]
            atom_sig = tuple(
                (a.GetSymbol(), a.GetFormalCharge(), a.GetNumRadicalElectrons())
                for a in mol.GetAtoms()
            )
            sig.append(hash(atom_sig))
            return tuple(sig)
        except Exception:
            return None

    def check_molecule_update(self):
        try:
            current_sig = self.get_mol_signature(self.context.current_molecule)
            if current_sig != self.last_seen_signature:
                self.load_molecule()
        except Exception as _e:
            logging.warning(
                "[charge_editor.py:check_molecule_update] silenced: %s", _e
            )

    def _atom_symbol(self, atom):
        if atom.HasProp("custom_symbol"):
            return atom.GetProp("custom_symbol")
        return atom.GetSymbol()

    def _row_atom_index(self, row):
        """The molecule atom index shown in a table row, or None."""
        try:
            return int(self.table.item(row, self.COL_IDX).text())
        except (AttributeError, ValueError):
            return None

    def load_molecule(self):
        self.table.blockSignals(True)
        self.table.setRowCount(0)

        mol = self.context.current_molecule
        self.last_seen_signature = self.get_mol_signature(mol)

        if not mol or not mol.GetNumAtoms():
            self.table.blockSignals(False)
            self.update_summary()
            return

        for atom in mol.GetAtoms():
            row = self.table.rowCount()
            self.table.insertRow(row)
            idx = atom.GetIdx()

            for col, text in (
                (self.COL_IDX, str(idx)),
                (self.COL_SYMBOL, self._atom_symbol(atom)),
            ):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row, col, item)

            charge_spin = QSpinBox()
            charge_spin.setRange(CHARGE_MIN, CHARGE_MAX)
            charge_spin.setValue(int(atom.GetFormalCharge()))
            charge_spin.valueChanged.connect(partial(self.on_charge_changed, idx))
            self.table.setCellWidget(row, self.COL_CHARGE, charge_spin)

            radical_spin = QSpinBox()
            radical_spin.setRange(RADICAL_MIN, RADICAL_MAX)
            radical_spin.setValue(int(atom.GetNumRadicalElectrons()))
            radical_spin.valueChanged.connect(partial(self.on_radical_changed, idx))
            self.table.setCellWidget(row, self.COL_RADICAL, radical_spin)

            try:
                h_count = atom.GetTotalNumHs()
            except Exception:
                h_count = 0
            h_item = QTableWidgetItem(str(h_count))
            h_item.setFlags(h_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, self.COL_HS, h_item)

        self.table.blockSignals(False)
        self.update_summary()

    def update_summary(self):
        mol = self.context.current_molecule
        if not mol or not mol.GetNumAtoms():
            self.summary_label.setText("Total charge: 0    Multiplicity: 1")
            return
        total_charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        unpaired = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        multiplicity = unpaired + 1
        self.summary_label.setText(
            f"Total charge: {total_charge:+d}    "
            f"Multiplicity: {multiplicity} ({unpaired} unpaired e⁻)"
        )

    # ------------------------------------------------------------------
    # edit operations (each commits immediately with an undo checkpoint)
    # ------------------------------------------------------------------

    def _commit(self, rw, message):
        try:
            Chem.SanitizeMol(rw)
        except Exception:
            try:
                rw.UpdatePropertyCache(strict=False)
                Chem.GetSSSR(rw)
            except Exception as _e:
                logging.warning(
                    "[charge_editor.py:_commit] sanitize fallback: %s", _e
                )
        self.context.current_molecule = rw.GetMol()
        self.context.push_undo_checkpoint()
        self.last_seen_signature = self.get_mol_signature(
            self.context.current_molecule
        )
        refresh = getattr(self.context, "refresh_3d_view", None)
        if callable(refresh):
            refresh()
        for name in ("refresh_2d_scene", "refresh_ui"):
            fn = getattr(self.context, name, None)
            if callable(fn):
                try:
                    fn()
                except Exception as _e:
                    logging.warning(
                        "[charge_editor.py:_commit] %s: %s", name, _e
                    )
        self.load_molecule()
        self.context.show_status_message(message)

    def on_charge_changed(self, atom_idx, value):
        mol = self.context.current_molecule
        if not mol:
            return
        try:
            rw = Chem.RWMol(mol)
            atom = rw.GetAtomWithIdx(int(atom_idx))
            atom.SetFormalCharge(int(value))
            atom.SetNoImplicit(True)
            self._commit(rw, f"Atom {atom_idx} formal charge set to {value:+d}.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to set charge: {str(e)}")

    def on_radical_changed(self, atom_idx, value):
        mol = self.context.current_molecule
        if not mol:
            return
        try:
            rw = Chem.RWMol(mol)
            atom = rw.GetAtomWithIdx(int(atom_idx))
            atom.SetNumRadicalElectrons(int(value))
            atom.SetNoImplicit(True)
            self._commit(rw, f"Atom {atom_idx} radical electrons set to {value}.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to set radicals: {str(e)}")

    def clear_all_charges(self):
        mol = self.context.current_molecule
        if not mol or not mol.GetNumAtoms():
            self.context.show_status_message("No molecule loaded.")
            return
        try:
            rw = Chem.RWMol(mol)
            for atom in rw.GetAtoms():
                atom.SetFormalCharge(0)
            self._commit(rw, "Cleared all formal charges.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to clear charges: {str(e)}")

    def clear_all_radicals(self):
        mol = self.context.current_molecule
        if not mol or not mol.GetNumAtoms():
            self.context.show_status_message("No molecule loaded.")
            return
        try:
            rw = Chem.RWMol(mol)
            for atom in rw.GetAtoms():
                atom.SetNumRadicalElectrons(0)
            self._commit(rw, "Cleared all radical electrons.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to clear radicals: {str(e)}")

    # ------------------------------------------------------------------
    # 3D highlight of selected atoms
    # ------------------------------------------------------------------

    @staticmethod
    def _highlight_radius(mol, atom_idx):
        """Halo radius for an atom: 1.2x the app's 0.3-scaled VDW display radius."""
        try:
            atomic_num = mol.GetAtomWithIdx(atom_idx).GetAtomicNum()
            if atomic_num > 0:
                pt = Chem.GetPeriodicTable()
                return pt.GetRvdw(atomic_num) * 1.2 * 0.3
        except Exception as _e:
            logging.warning(
                "[charge_editor.py:_highlight_radius] silenced: %s", _e
            )
        return 1.5 * 0.3  # ghost/unknown atoms keep the old fixed size

    def highlight_selected_atoms(self):
        plotter = self.context.plotter
        if not plotter:
            return
        try:
            cam = plotter.camera_position
        except (AttributeError, RuntimeError, TypeError):
            cam = None

        rows = sorted(set(index.row() for index in self.table.selectedIndexes()))
        mol = self.context.current_molecule
        centers = []
        radii = []
        if mol and mol.GetNumConformers():
            conf = mol.GetConformer()
            for row in rows:
                idx = self._row_atom_index(row)
                if idx is None:
                    continue
                p = conf.GetAtomPosition(idx)
                centers.append([p.x, p.y, p.z])
                radii.append(self._highlight_radius(mol, idx))

        if not centers:
            plotter.remove_actor("charge_editor_selection")
        else:
            cloud = pv.PolyData(np.array(centers))
            cloud["radii"] = radii
            glyph = cloud.glyph(
                geom=pv.Sphere(radius=1.0), scale="radii", orient=False
            )
            plotter.add_mesh(
                glyph,
                name="charge_editor_selection",
                color="yellow",
                opacity=0.4,
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
            plotter.remove_actor("charge_editor_selection")
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
        win = ChargeEditorWindow(context)
        win.show()

    context.add_menu_action("3D Edit/Charge Editor...", show_editor)

    def on_document_reset():
        win = context.get_window("main_panel")
        if win:
            win.load_molecule()

    context.register_document_reset_handler(on_document_reset)
