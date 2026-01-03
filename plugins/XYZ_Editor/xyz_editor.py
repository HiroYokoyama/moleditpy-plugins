import sys
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem, 
    QPushButton, QHBoxLayout, QMessageBox, QHeaderView
)
from PyQt6.QtCore import Qt, QTimer
from rdkit import Chem
from rdkit.Geometry import Point3D
import pyvista as pv

PLUGIN_NAME = "XYZ Editor"
PLUGIN_VERSION = "2026.01.04"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "A table-based editor for atom coordinates and symbols, supporting ghost atoms."

class XYZEditorWindow(QWidget):
    def __init__(self, context):
        self.mw = context.get_main_window()
        super().__init__(self.mw) # Parent to main window to stay on top of it
        self.setWindowFlags(Qt.WindowType.Window) # Ensure it's a separate window, not embedded
        self.context = context
        self.setWindowTitle("XYZ Editor")
        self.resize(600, 400)
        self.init_ui()
        self.last_seen_signature = None
        self.load_molecule()
        
        # Auto-update mechanism
        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.check_molecule_update)
        self.update_timer.start(500) # Check every 500ms

        # Connect to selection changes in the main window if possible
        # For now, we'll just handle our own selection

    def init_ui(self):
        layout = QVBoxLayout(self)

        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(["Index", "Symbol", "X", "Y", "Z"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
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

        layout.addLayout(btn_layout)

        # Add/Remove Buttons
        edit_layout = QHBoxLayout()
        self.add_btn = QPushButton("Add Atom")
        self.add_btn.clicked.connect(self.add_atom_row)
        edit_layout.addWidget(self.add_btn)

        self.remove_btn = QPushButton("Remove Selected")
        self.remove_btn.clicked.connect(self.remove_selected_rows)
        edit_layout.addWidget(self.remove_btn)
        
        self.copy_btn = QPushButton("Copy to Clipboard")
        self.copy_btn.clicked.connect(self.copy_to_clipboard)
        edit_layout.addWidget(self.copy_btn)
        
        layout.addLayout(edit_layout)

    def closeEvent(self, event):
        # Cleanup highlights when window is closed
        if hasattr(self.mw, 'plotter'):
            self.mw.plotter.remove_actor("xyz_selection")
            self.mw.plotter.render()
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
                pos = conf.GetAtomPosition(0) # distinct check
                sig.append(f"{pos.x:.4f},{pos.y:.4f},{pos.z:.4f}")
            
            return tuple(sig)
        except:
            return None

    def check_molecule_update(self):
        try:
            current_mol = self.mw.current_mol
            current_sig = self.get_mol_signature(current_mol)
            
            if current_sig != self.last_seen_signature:
                self.load_molecule()
        except Exception:
            pass
            
    def load_molecule(self):
        self.table.blockSignals(True)
        self.table.setRowCount(0)

        mol = self.mw.current_mol
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
            if atom.GetAtomicNum() == 0:
                # Check for custom label property first (our new standard)
                if atom.HasProp("custom_symbol"):
                    symbol = atom.GetProp("custom_symbol")
                # Fallback to dummyLabel if custom_symbol missing (legacy/other plugins)
                elif atom.HasProp("dummyLabel"):
                    val = atom.GetProp("dummyLabel")
                    # If it's the safe proxy "X", ignore it? No, might be legacy "X".
                    # We trust it if custom_symbol is not set.
                    symbol = val
                elif symbol == '*': 
                    symbol = 'X' 

            self._add_row(idx, symbol, pos.x, pos.y, pos.z)

        self.table.blockSignals(False)

    def _add_row(self, idx, symbol, x, y, z):
        row = self.table.rowCount()
        self.table.insertRow(row)

        # Index (Read-only for existing, Placeholder for new)
        idx_str = str(idx + 1) if idx is not None else "+"
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
        rows = sorted(set(index.row() for index in self.table.selectedIndexes()), reverse=True)
        for row in rows:
            self.table.removeRow(row)

    def copy_to_clipboard(self):
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
            
        clipboard_text = "\n".join(lines)
        from PyQt6.QtGui import QGuiApplication
        clipboard = QGuiApplication.clipboard()
        clipboard.setText(clipboard_text)
        self.mw.statusBar().showMessage("XYZ data copied to clipboard.", 3000)

    def highlight_selected_atoms(self):
        # Identify selected rows
        rows = set(index.row() for index in self.table.selectedIndexes())
        
        if not rows:
            # Clear highlights
            self.mw.plotter.remove_actor("xyz_selection")
            self.mw.plotter.render()
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
            self.mw.plotter.remove_actor("xyz_selection")
            self.mw.plotter.render()
            return
            
        # Create polydata for points
        poly = pv.PolyData(points)
        poly["radii"] = radii # Add scalar array
        
        # Add spheres at these points, scaled by radius
        spheres = poly.glyph(geom=pv.Sphere(radius=1.0), scale="radii")
        
        self.mw.plotter.add_mesh(
            spheres, 
            name="xyz_selection", 
            color="yellow", 
            opacity=0.5,
            pickable=False
        )
        self.mw.plotter.render()

    def on_item_changed(self, item):
        # Update highlight if the row is selected to show live movement
        if item.row() in set(index.row() for index in self.table.selectedIndexes()):
            self.highlight_selected_atoms()

    def apply_changes(self):
        mol = self.mw.current_mol
        # Create new editable molecule from scratch or copy
        if mol:
            rw_mol = Chem.RWMol(mol)
        else:
            rw_mol = Chem.RWMol()

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
                    # Check if it is a known element
                    pt = Chem.GetPeriodicTable()
                    # RDKit uses title case for symbols (e.g. "He", not "HE")
                    # but GetAtomicNumber is smart enough usually.
                    atomic_num = pt.GetAtomicNumber(symbol.capitalize())
                    atom = Chem.Atom(atomic_num)
                except RuntimeError:
                    # Unknown symbol -> Treat as Ghost/Dummy Atom
                    atom = Chem.Atom(0)
                    atom.SetProp("custom_symbol", symbol) # Store true label for export/editor
                    # Set dummyLabel to something safe for 3D viewer (avoid "Bq" -> crash)
                    # We use "*" which has AtomicNum 0, avoiding the "Element not found" crash.
                    #atom.SetProp("dummyLabel", "*")

                new_idx = new_rw_mol.AddAtom(atom)
                atom_coords.append(Point3D(x, y, z))
                
                # Track Index Mapping
                idx_text = self.table.item(row, 0).text()
                if idx_text != "+":
                    old_idx = int(idx_text) - 1
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
                        new_rw_mol.AddBond(old_to_new_map[b], old_to_new_map[e], bond.GetBondType())

            # Add Conformer
            new_rw_mol.AddConformer(conf)
            
            # Commit changes
            self.mw.push_undo_state()
            # Update properties and ring info to avoid RDKit errors
            try:
                Chem.SanitizeMol(new_rw_mol)
            except:
                new_rw_mol.UpdatePropertyCache(strict=False)
                Chem.GetSSSR(new_rw_mol)
            self.context.current_molecule = new_rw_mol.GetMol() 
            self.last_seen_signature = self.get_mol_signature(self.context.current_molecule)
            
            # Refresh visualization
            self.mw.plotter.reset_camera() 
            
            # Reload table to get clean indices and ensure properties stuck
            self.load_molecule()
            
            self.mw.statusBar().showMessage("XYZ changes applied.", 3000)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply changes: {str(e)}")


def initialize(context):
    # Store the window reference to prevent garbage collection
    mw = context.get_main_window()
    
    def show_editor():
        if not hasattr(mw, 'xyz_editor_window') or mw.xyz_editor_window is None:
            mw.xyz_editor_window = XYZEditorWindow(context)
        
        mw.xyz_editor_window.show()
        mw.xyz_editor_window.raise_()
        mw.xyz_editor_window.activateWindow()
        mw.xyz_editor_window.load_molecule() # Refresh on show

    context.add_menu_action("3D Edit/XYZ Editor...", show_editor)

    # Persistence Handling
    def save_plugin_state():
        # Only save name mapping (custom labels)
        mol = mw.current_mol
        if not mol:
            return {}
            
        custom_labels = {}
        for atom in mol.GetAtoms():
            if atom.HasProp("custom_symbol"):
                custom_labels[atom.GetIdx()] = atom.GetProp("custom_symbol")
            elif atom.HasProp("dummyLabel"):
                # If only dummyLabel exists (legacy), save it.
                # But check if it's our safe proxy "X" which we might not want as the "name"?
                # Actually, if custom_symbol is missing, dummyLabel IS the name.
                custom_labels[atom.GetIdx()] = atom.GetProp("dummyLabel")
            
        return {"custom_labels": custom_labels}

    def load_plugin_state(data):
        # Restore custom labels (name mapping)
        custom_labels = data.get("custom_labels")
        if not custom_labels:
            return

        mol = mw.current_mol
        if mol:
            for idx_str, label in custom_labels.items():
                try:
                    idx = int(idx_str)
                    if idx < mol.GetNumAtoms():
                        atom = mol.GetAtomWithIdx(idx)
                        atom.SetProp("custom_symbol", label)
                        # Also set safe proxy for 3D viewer
                        atom.SetProp("dummyLabel", "*")
                except:
                    pass

        # If the editor is open, refresh it
        if hasattr(mw, 'xyz_editor_window') and mw.xyz_editor_window and mw.xyz_editor_window.isVisible():
            mw.xyz_editor_window.load_molecule()

    context.register_save_handler(save_plugin_state)
    context.register_load_handler(load_plugin_state)

def run(mw):
    # Legacy support
    class LegacyContext:
        def __init__(self, mw):
            self.mw = mw
            self.current_molecule = mw.current_mol
        def get_main_window(self):
            return self.mw
        @property
        def current_molecule(self):
            return self.mw.current_mol
        @current_molecule.setter
        def current_molecule(self, mol):
            self.mw.current_mol = mol
            self.mw.draw_molecule_3d(mol)

    context = LegacyContext(mw)

    if not hasattr(mw, 'xyz_editor_window') or mw.xyz_editor_window is None:
        mw.xyz_editor_window = XYZEditorWindow(context)
    
    mw.xyz_editor_window.show()
    mw.xyz_editor_window.raise_()
    mw.xyz_editor_window.activateWindow()
    mw.xyz_editor_window.load_molecule()
