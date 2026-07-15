import numpy as np
from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QHBoxLayout,
    QMessageBox,
    QHeaderView,
    QFileDialog,
    QCheckBox,
    QDoubleSpinBox,
    QLabel,
)
from PyQt6.QtCore import (
    Qt,
    QTimer,
    QObject,
    QEvent,
    QItemSelection,
    QItemSelectionModel,
)
from PyQt6.QtGui import QShortcut, QKeySequence
from rdkit import Chem
from rdkit.Geometry import Point3D
import pyvista as pv
import logging


PLUGIN_NAME = "XYZ Editor"
PLUGIN_VERSION = "2026.07.15"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "A table-based editor for atom coordinates and symbols, supporting ghost atoms. Refactored for V3 API."
PLUGIN_CONTEXT = None


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


class XYZEditorWindow(QWidget):
    """
    Namespaced window for editing XYZ coordinates and symbols.
    Refactored for MoleditPy V3.0 API.
    """

    def __init__(self, context):
        super().__init__(parent=context.get_main_window())
        self.setWindowFlags(Qt.WindowType.Window)
        self.context = context
        self.setWindowTitle("XYZ Editor")
        self.resize(600, 400)
        self._click_filter = None
        self.init_ui()

        # Register window for V3 lifecycle management
        self.context.register_window("main_panel", self)
        self.last_seen_signature = None
        self.load_molecule()

        # Auto-update mechanism
        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.check_molecule_update)
        self.update_timer.start(500)  # Check every 500ms

        # Install click-detection filter on 3D plotter
        self._enable_plotter_picking()

    def init_ui(self):
        layout = QVBoxLayout(self)

        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(["Index", "Symbol", "X", "Y", "Z"])
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self.table.itemChanged.connect(self.on_item_changed)
        self.table.itemSelectionChanged.connect(self.highlight_selected_atoms)

        # UI Styling to prevent overlap issues
        self.table.setStyleSheet("""
            QTableWidget {
                background-color: white;
                color: black;
                gridline-color: #d0d0d0;
            }
            QHeaderView::section {
                background-color: #f0f0f0;
                padding: 4px;
                border: 1px solid #d0d0d0;
            }
            QLineEdit {
                background-color: white;
                color: black;
                selection-background-color: #0078d7;
            }
        """)
        layout.addWidget(self.table)

        # Buttons
        btn_layout = QHBoxLayout()

        self.apply_btn = QPushButton("Apply to View")
        self.apply_btn.clicked.connect(self.apply_changes)
        btn_layout.addWidget(self.apply_btn)

        self.save_btn = QPushButton("Save as XYZ...")
        self.save_btn.clicked.connect(self.save_as_xyz)
        btn_layout.addWidget(self.save_btn)

        layout.addLayout(btn_layout)

        # Add/Remove Buttons
        edit_layout = QHBoxLayout()
        self.add_btn = QPushButton("Add Atom")
        self.add_btn.clicked.connect(self.add_atom_row)
        edit_layout.addWidget(self.add_btn)

        self.remove_btn = QPushButton("Remove Selected")
        self.remove_btn.setToolTip("Remove rows from the table only (press Apply to commit)")
        self.remove_btn.clicked.connect(self.remove_selected_rows)
        edit_layout.addWidget(self.remove_btn)

        self.delete_btn = QPushButton("Delete Atoms")
        self.delete_btn.setToolTip("Delete selected atoms from the molecule immediately")
        self.delete_btn.clicked.connect(self.delete_selected_atoms)
        edit_layout.addWidget(self.delete_btn)

        self.copy_btn = QPushButton("Copy to Clipboard")
        self.copy_btn.clicked.connect(self.copy_to_clipboard)
        edit_layout.addWidget(self.copy_btn)

        self.unselect_btn = QPushButton("Unselect")
        self.unselect_btn.clicked.connect(self.unselect_all)
        edit_layout.addWidget(self.unselect_btn)

        layout.addLayout(edit_layout)

        select_layout = QHBoxLayout()
        self.whole_mol_cb = QCheckBox("Select whole molecule on 3D click")
        self.whole_mol_cb.setToolTip(
            "Clicking an atom in the 3D view selects its entire connected fragment. "
            "Hold Ctrl to add/remove fragments or atoms from the selection."
        )
        select_layout.addWidget(self.whole_mol_cb)
        select_layout.addStretch()
        layout.addLayout(select_layout)

        dup_layout = QHBoxLayout()
        self.duplicate_btn = QPushButton("Duplicate")
        self.duplicate_btn.setToolTip(
            "Duplicate the selected atoms (whole molecule if nothing is selected), "
            "offset by the values on the right"
        )
        self.duplicate_btn.clicked.connect(self.duplicate_atoms)
        dup_layout.addWidget(self.duplicate_btn)
        dup_layout.addWidget(QLabel("Offset:"))
        self.dup_offset = []
        for axis in ("X", "Y", "Z"):
            spin = QDoubleSpinBox()
            spin.setRange(-1000.0, 1000.0)
            spin.setDecimals(3)
            spin.setSingleStep(0.5)
            spin.setValue(1.0)
            spin.setPrefix(f"{axis}: ")
            dup_layout.addWidget(spin)
            self.dup_offset.append(spin)
        dup_layout.addStretch()
        layout.addLayout(dup_layout)

        self._del_shortcut = QShortcut(QKeySequence(Qt.Key.Key_Delete), self.table)
        self._del_shortcut.setContext(Qt.ShortcutContext.WidgetShortcut)
        self._del_shortcut.activated.connect(self.remove_selected_rows)

    def _enable_plotter_picking(self):
        """Install Qt event filter on the 3D plotter interactor widget for atom click detection."""
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
            logging.warning("[xyz_editor.py:_enable_plotter_picking] silenced: %s", _e)

    def _disable_plotter_picking(self):
        """Remove the event filter from the 3D plotter interactor widget."""
        try:
            plotter = self.context.plotter
            interactor = getattr(plotter, "interactor", None) if plotter else None
            if interactor and self._click_filter:
                interactor.removeEventFilter(self._click_filter)
        except Exception as _e:
            logging.warning("[xyz_editor.py:_disable_plotter_picking] silenced: %s", _e)
        self._click_filter = None

    def _on_plotter_click(self, x, y, widget, modifiers):
        """Called when a non-drag left click is detected on the 3D plotter interactor widget."""
        try:
            import vtk

            mw = self.context.get_main_window()
            v3d = getattr(mw, "view_3d_manager", None) if mw else None
            plotter = self.context.plotter
            if not v3d or not plotter:
                return

            # Convert Qt (top-left origin) → VTK (bottom-left origin)
            # widget here is the interactor (QVTKRenderWindowInteractor)
            vtk_y = widget.height() - y

            picker = vtk.vtkCellPicker()
            picker.SetTolerance(0.005)
            picker.Pick(x, vtk_y, 0, plotter.renderer)
            picked_actor = picker.GetActor()

            atom_actor = getattr(v3d, "atom_actor", None)
            if picked_actor is None or picked_actor is not atom_actor:
                return

            # Find closest atom to pick position
            pick_pos = picker.GetPickPosition()
            mol = self.context.current_mol
            if not mol or not mol.GetNumConformers():
                return

            conf = mol.GetConformer()
            best_idx = -1
            best_dist = float("inf")
            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                pos = conf.GetAtomPosition(idx)
                dx = pos.x - pick_pos[0]
                dy = pos.y - pick_pos[1]
                dz = pos.z - pick_pos[2]
                dist = dx * dx + dy * dy + dz * dz
                if dist < best_dist:
                    best_dist = dist
                    best_idx = idx

            if best_idx < 0:
                return

            atom_indices = {best_idx}
            if self.whole_mol_cb.isChecked():
                atom_indices = self._fragment_atom_indices(mol, best_idx)

            row_map = self._atom_index_to_row_map()
            rows = {row_map[i] for i in atom_indices if i in row_map}
            if not rows:
                return
            target_row = row_map.get(best_idx, min(rows))

            ctrl_held = bool(modifiers & Qt.KeyboardModifier.ControlModifier)
            self._select_rows(rows, ctrl_held, target_row)
            self.table.scrollTo(self.table.model().index(target_row, 0))
            self.highlight_selected_atoms()
        except Exception as _e:
            logging.warning("[xyz_editor.py:223] silenced: %s", _e)

    def _fragment_atom_indices(self, mol, atom_idx):
        """Return the set of atom indices in the connected fragment containing atom_idx."""
        try:
            frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
            for frag in frags:
                if atom_idx in frag:
                    return set(frag)
        except Exception as _e:
            logging.warning("[xyz_editor.py:_fragment_atom_indices] silenced: %s", _e)
        return {atom_idx}

    def _atom_index_to_row_map(self):
        """Map molecule atom indices (column 0) to table row numbers."""
        row_map = {}
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)
            if item and item.text() not in ("", "+"):
                try:
                    row_map[int(item.text())] = row
                except ValueError as _e:
                    logging.warning("[xyz_editor.py:_atom_index_to_row_map] silenced: %s", _e)
        return row_map

    def _select_rows(self, rows, ctrl_held, anchor_row):
        """Select the given table rows. Without Ctrl the selection is replaced;
        with Ctrl the rows are added, or removed if the anchor row was already selected."""
        model = self.table.model()
        sm = self.table.selectionModel()
        last_col = self.table.columnCount() - 1
        selection = QItemSelection()
        for r in sorted(rows):
            selection.select(model.index(r, 0), model.index(r, last_col))

        if ctrl_held:
            anchor_selected = anchor_row in {
                i.row() for i in self.table.selectedIndexes()
            }
            flag = (
                QItemSelectionModel.SelectionFlag.Deselect
                if anchor_selected
                else QItemSelectionModel.SelectionFlag.Select
            )
        else:
            flag = QItemSelectionModel.SelectionFlag.ClearAndSelect

        self.table.blockSignals(True)
        sm.select(selection, flag)
        self.table.blockSignals(False)

    def closeEvent(self, event):
        self._disable_plotter_picking()
        plotter = self.context.plotter
        if plotter:
            plotter.remove_actor("xyz_selection")
            plotter.render()
        super().closeEvent(event)

    def get_mol_signature(self, mol):
        if not mol:
            return None
        try:
            # Create a lightweight signature to detect changes
            # 1. Object Identity (primary check for new files)
            # 2. Number of Atoms/Bonds (check for structure changes)
            # 3. First atom position (check for movement/conformer updates)

            sig = [id(mol), mol.GetNumAtoms(), mol.GetNumBonds()]

            if mol.GetNumAtoms() > 0:
                conf = mol.GetConformer()
                # Use a robust signature: hash of all coordinates (rounded to 4 decimals)
                # to detect any change in any atom without cancellation issues.
                pos_array = conf.GetPositions()
                coord_hash = hash(np.round(pos_array, 4).tobytes())
                sig.append(coord_hash)

            return tuple(sig)
        except Exception:
            return None

    def check_molecule_update(self):
        try:
            current_mol = self.context.current_molecule
            current_sig = self.get_mol_signature(current_mol)

            if current_sig != self.last_seen_signature:
                self.load_molecule()
        except Exception as _e:
            logging.warning("[xyz_editor.py:264] silenced: %s", _e)

    def load_molecule(self):
        self.table.blockSignals(True)
        self.table.setRowCount(0)

        mol = self.context.current_molecule
        self.last_seen_signature = self.get_mol_signature(mol)

        if not mol:
            self.table.blockSignals(False)
            return

        conf = mol.GetConformer()
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)

            # Atomic Symbol handling
            symbol = atom.GetSymbol()

            # Check for custom label property first (our new standard)
            if atom.HasProp("custom_symbol"):
                symbol = atom.GetProp("custom_symbol")
            # Fallback to dummyLabel if custom_symbol missing (legacy/other plugins)
            elif atom.HasProp("dummyLabel"):
                symbol = atom.GetProp("dummyLabel")

            # Note: We no longer map '*' to 'X' here to allow asterisks as requested.
            # RDKit dummy atoms default to '*' symbol.

            self._add_row(idx, symbol, pos.x, pos.y, pos.z)

        self.table.blockSignals(False)

    def _add_row(self, idx, symbol, x, y, z):
        row = self.table.rowCount()
        self.table.insertRow(row)

        # Index (Read-only for existing, Placeholder for new)
        idx_str = str(idx) if idx is not None else "+"
        item_idx = QTableWidgetItem(idx_str)
        item_idx.setFlags(item_idx.flags() & ~Qt.ItemFlag.ItemIsEditable)
        self.table.setItem(row, 0, item_idx)

        # Symbol
        self.table.setItem(row, 1, QTableWidgetItem(symbol))

        # Coordinates
        self.table.setItem(row, 2, QTableWidgetItem(f"{x:.5f}"))
        self.table.setItem(row, 3, QTableWidgetItem(f"{y:.5f}"))
        self.table.setItem(row, 4, QTableWidgetItem(f"{z:.5f}"))

    def add_atom_row(self):
        # Add a new row with default values
        self._add_row(None, "C", 0.0, 0.0, 0.0)
        self.table.scrollToBottom()

    def remove_selected_rows(self):
        rows = sorted(
            set(index.row() for index in self.table.selectedIndexes()), reverse=True
        )
        for row in rows:
            self.table.removeRow(row)

    def delete_selected_atoms(self):
        """Remove the selected rows and immediately apply the result to the molecule."""
        rows = set(index.row() for index in self.table.selectedIndexes())
        if not rows:
            self.context.show_status_message("No atoms selected to delete.")
            return
        self.remove_selected_rows()
        self.apply_changes()

    def _selected_atom_indices(self):
        """Atom indices (column 0) of the selected rows; unapplied '+' rows are skipped."""
        indices = set()
        for row in set(index.row() for index in self.table.selectedIndexes()):
            item = self.table.item(row, 0)
            if item and item.text() not in ("", "+"):
                try:
                    indices.add(int(item.text()))
                except ValueError as _e:
                    logging.warning(
                        "[xyz_editor.py:_selected_atom_indices] silenced: %s", _e
                    )
        return indices

    def duplicate_atoms(self):
        """Duplicate the selected atoms (or the whole molecule if nothing is
        selected) with the configured coordinate offset, preserving bonds
        between duplicated atoms. The new copy is left selected."""
        mol = self.context.current_molecule
        if not mol or not mol.GetNumAtoms() or not mol.GetNumConformers():
            self.context.show_status_message("No molecule to duplicate.")
            return

        had_selection = bool(self.table.selectedIndexes())
        sel = self._selected_atom_indices()
        if had_selection and not sel:
            self.context.show_status_message(
                "Selected rows are not applied yet — press Apply first."
            )
            return
        if not sel:
            sel = set(range(mol.GetNumAtoms()))
        dx, dy, dz = (spin.value() for spin in self.dup_offset)

        try:
            rw = Chem.RWMol(mol)
            src_conf = mol.GetConformer()
            mapping = {}
            for idx in sorted(sel):
                src = mol.GetAtomWithIdx(idx)
                atom = Chem.Atom(src.GetAtomicNum())
                atom.SetFormalCharge(src.GetFormalCharge())
                atom.SetNoImplicit(src.GetNoImplicit())
                atom.SetIsAromatic(src.GetIsAromatic())
                atom.SetNumRadicalElectrons(src.GetNumRadicalElectrons())
                if src.HasProp("custom_symbol"):
                    atom.SetProp("custom_symbol", src.GetProp("custom_symbol"))
                mapping[idx] = rw.AddAtom(atom)

            conf = rw.GetConformer()
            for idx, new_idx in mapping.items():
                pos = src_conf.GetAtomPosition(idx)
                conf.SetAtomPosition(
                    new_idx, Point3D(pos.x + dx, pos.y + dy, pos.z + dz)
                )

            for bond in mol.GetBonds():
                b = bond.GetBeginAtomIdx()
                e = bond.GetEndAtomIdx()
                if b in mapping and e in mapping:
                    rw.AddBond(mapping[b], mapping[e], bond.GetBondType())

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

            row_map = self._atom_index_to_row_map()
            new_rows = {row_map[i] for i in mapping.values() if i in row_map}
            if new_rows:
                self._select_rows(new_rows, False, min(new_rows))
                self.highlight_selected_atoms()

            self.context.show_status_message(f"Duplicated {len(mapping)} atom(s).")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to duplicate: {str(e)}")

    def unselect_all(self):
        self.table.clearSelection()

    def copy_to_clipboard(self):
        lines = self._generate_xyz_content()
        clipboard_text = "\n".join(lines)
        from PyQt6.QtGui import QGuiApplication

        clipboard = QGuiApplication.clipboard()
        clipboard.setText(clipboard_text)
        self.context.show_status_message("XYZ data copied to clipboard.")

    def save_as_xyz(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save XYZ File", "", "XYZ Files (*.xyz);;All Files (*)"
        )
        if not file_path:
            return

        if not file_path.lower().endswith(".xyz"):
            file_path += ".xyz"

        try:
            lines = self._generate_xyz_content()
            with open(file_path, "w", encoding="utf-8") as f:
                f.write("\n".join(lines))
            self.context.show_status_message(f"XYZ saved to {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save XYZ: {str(e)}")

    def _generate_xyz_content(self):
        lines = []
        atom_count = self.table.rowCount()
        lines.append(str(atom_count))
        lines.append("Generated by MoleditPy XYZ Editor")

        for row in range(atom_count):
            symbol = self.table.item(row, 1).text().strip()
            x = self.table.item(row, 2).text()
            y = self.table.item(row, 3).text()
            z = self.table.item(row, 4).text()
            lines.append(f"{symbol:<4} {x:>12} {y:>12} {z:>12}")
        return lines

    def highlight_selected_atoms(self):
        # Identify selected rows
        rows = set(index.row() for index in self.table.selectedIndexes())

        if not rows:
            plotter = self.context.plotter
            if plotter:
                try:
                    cam = plotter.camera_position
                except (AttributeError, RuntimeError, TypeError):
                    cam = None
                plotter.remove_actor("xyz_selection")
                if cam is not None:
                    try:
                        plotter.camera_position = cam
                    except (AttributeError, RuntimeError, TypeError):
                        pass
                plotter.render()
            return

        points = []
        radii = []
        pt = Chem.GetPeriodicTable()

        for row in rows:
            # Try to get coordinates from table (handles new/edited atoms)
            try:
                symbol = self.table.item(row, 1).text().strip()
                x = float(self.table.item(row, 2).text())
                y = float(self.table.item(row, 3).text())
                z = float(self.table.item(row, 4).text())
                points.append([x, y, z])

                # Determine radius based on element
                try:
                    atomic_num = pt.GetAtomicNumber(symbol.capitalize())
                    # Scaling factor 1.2 * 0.3 as requested (Highlight Halo)
                    radius = pt.GetRvdw(atomic_num) * 1.2 * 0.3
                except RuntimeError:
                    # Ghost/Unknown -> Default size
                    radius = 1.5 * 0.3

                radii.append(radius)

            except ValueError:
                continue

        if not points:
            plotter = self.context.plotter
            if plotter:
                try:
                    cam = plotter.camera_position
                except (AttributeError, RuntimeError, TypeError):
                    cam = None
                plotter.remove_actor("xyz_selection")
                if cam is not None:
                    try:
                        plotter.camera_position = cam
                    except (AttributeError, RuntimeError, TypeError):
                        pass
                plotter.render()
            return

        # Create polydata for points
        poly = pv.PolyData(points)
        poly["radii"] = radii  # Add scalar array

        # Add spheres at these points, scaled by radius
        # orient=False added to suppress PyVista UserWarning (No vector-like data to use for orient)
        spheres = poly.glyph(geom=pv.Sphere(radius=1.0), scale="radii", orient=False)

        plotter = self.context.plotter
        if plotter:
            try:
                cam = plotter.camera_position
            except (AttributeError, RuntimeError, TypeError):
                cam = None
            plotter.add_mesh(
                spheres,
                name="xyz_selection",
                color="yellow",
                opacity=0.5,
                pickable=False,
            )
            if cam is not None:
                try:
                    plotter.camera_position = cam
                except (AttributeError, RuntimeError, TypeError):
                    pass
            plotter.render()

    def on_item_changed(self, item):
        # Update highlight if the row is selected to show live movement
        if item.row() in set(index.row() for index in self.table.selectedIndexes()):
            self.highlight_selected_atoms()

    def apply_changes(self):
        # self.context.push_undo_checkpoint() # MOVED TO END
        mol = self.context.current_molecule
        # Create new editable molecule from scratch or copy
        if mol:
            Chem.RWMol(mol)
        else:
            Chem.RWMol()

        # We need to rebuild the molecule based on the table content
        # Because rows might have been deleted or added/reordered

        new_rw_mol = Chem.RWMol()
        # Store coords to set AFTER adding all atoms (requires sizing conformer)
        atom_coords = []

        # Map old_idx -> new_idx to reconstruct bonds if possible
        old_to_new_map = {}

        table_rows = self.table.rowCount()

        try:
            for row in range(table_rows):
                # Parse inputs
                symbol = self.table.item(row, 1).text().strip()
                try:
                    x = float(self.table.item(row, 2).text())
                    y = float(self.table.item(row, 3).text())
                    z = float(self.table.item(row, 4).text())
                except ValueError:
                    # Don't show error while typing, just abort update
                    return

                # Create Atom
                try:
                    pt = Chem.GetPeriodicTable()
                    # 1. Selection logic: Strict Case matching for Element Identification
                    # We only treat it as a real element if the casing is perfect (C, Ag, He).
                    at_num = -1
                    try:
                        potential_num = pt.GetAtomicNumber(symbol)
                        if pt.GetElementSymbol(potential_num) == symbol:
                            at_num = potential_num
                    except Exception as _e:
                        logging.warning("[xyz_editor.py:484] silenced: %s", _e)

                    if at_num > 0:
                        atom = Chem.Atom(at_num)
                        # No custom_symbol needed if it matches canonical perfectly
                    else:
                        # 2. Try Prefix Match (Strict Case) for things like "Ag*"
                        found_atomic_num = 0
                        for i in range(len(symbol), 0, -1):
                            prefix = symbol[:i]
                            try:
                                p_num = pt.GetAtomicNumber(prefix)
                                if pt.GetElementSymbol(p_num) == prefix:
                                    # NEW LOGIC: If we found a prefix match, ensure the next character
                                    # is NOT an alphabet (which would imply it's part of a longer
                                    # unidentified symbol like 'Bq' or 'Calpha')
                                    if i < len(symbol) and symbol[i].isalpha():
                                        continue

                                    found_atomic_num = p_num
                                    break
                            except Exception:
                                continue

                        # Create atom (defaults to dummy 0 if no prefix found)
                        atom = Chem.Atom(found_atomic_num)
                        # Store the EXACT string (c, ag, Ag*, etc.)
                        atom.SetProp("custom_symbol", symbol)
                except Exception as e:
                    logging.warning("Error creating atom for row %s: %s", row, e)
                    # Fallback to Carbon if something goes wrong
                    atom = Chem.Atom(6)

                new_idx = new_rw_mol.AddAtom(atom)
                atom_coords.append(Point3D(x, y, z))

                # Track Index Mapping
                idx_text = self.table.item(row, 0).text()
                if idx_text != "+":
                    old_idx = int(idx_text)
                    old_to_new_map[old_idx] = new_idx

            # Now create conformer with correct size
            conf = Chem.Conformer(new_rw_mol.GetNumAtoms())
            for idx, pt in enumerate(atom_coords):
                conf.SetAtomPosition(idx, pt)

            # Re-add bonds if they exist between surviving atoms
            if mol:
                for bond in mol.GetBonds():
                    b = bond.GetBeginAtomIdx()
                    e = bond.GetEndAtomIdx()
                    if b in old_to_new_map and e in old_to_new_map:
                        new_rw_mol.AddBond(
                            old_to_new_map[b], old_to_new_map[e], bond.GetBondType()
                        )

            # Add Conformer
            new_rw_mol.AddConformer(conf)

            # Commit changes
            # self.mw.edit_actions_manager.push_undo_state()  # MOVED TO START
            # Update properties and ring info to avoid RDKit errors
            try:
                Chem.SanitizeMol(new_rw_mol)
            except Exception:
                new_rw_mol.UpdatePropertyCache(strict=False)
                Chem.GetSSSR(new_rw_mol)
            self.context.current_molecule = new_rw_mol.GetMol()
            self.context.push_undo_checkpoint()
            self.last_seen_signature = self.get_mol_signature(
                self.context.current_molecule
            )

            # Refresh visualization
            self.context.reset_3d_camera()

            # Reload table to get clean indices and ensure properties stuck
            self.load_molecule()

            self.context.show_status_message("XYZ changes applied.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply changes: {str(e)}")


def initialize(context):
    """MoleditPy Plugin Entry Point (V3.0)"""
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

        win = XYZEditorWindow(context)
        win.show()

    context.add_menu_action("3D Edit/XYZ Editor...", show_editor)

    # Persistence Handling
    def save_plugin_state():
        mol = context.current_molecule
        if not mol:
            return {}
        labels = {
            a.GetIdx(): a.GetProp("custom_symbol")
            for a in mol.GetAtoms()
            if a.HasProp("custom_symbol")
        }
        return {"custom_labels": labels}

    def load_plugin_state(data):
        labels = data.get("custom_labels", {})
        if not isinstance(labels, dict):
            labels = {}
        mol = context.current_molecule
        if mol:
            for idx, lbl in labels.items():
                try:
                    mol.GetAtomWithIdx(int(idx)).SetProp("custom_symbol", lbl)
                except Exception as _e:
                    logging.warning("[xyz_editor.py:597] silenced: %s", _e)
            context.current_molecule = mol

    def on_document_reset():
        win = context.get_window("main_panel")
        if win:
            win.load_molecule()

    context.register_save_handler(save_plugin_state)
    context.register_load_handler(load_plugin_state)
    context.register_document_reset_handler(on_document_reset)


def run(mw):
    if hasattr(mw, "host"):
        mw = mw.host
    context = PLUGIN_CONTEXT
    if not context:
        return

    win = context.get_window("main_panel")
    if win is None:
        win = XYZEditorWindow(context)
    win.show()
    win.raise_()
    win.activateWindow()
