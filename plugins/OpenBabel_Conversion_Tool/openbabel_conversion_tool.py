import os
import sys
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QListWidget, QDialogButtonBox, 
    QLabel, QMessageBox, QFileDialog, QAbstractItemView
)
from PyQt6.QtCore import Qt, QTimer
from rdkit import Chem
from rdkit.Chem import AllChem

# Try importing openbabel
OBABEL_AVAILABLE = False
try:
    from openbabel import pybel
    OBABEL_AVAILABLE = True
except ImportError:
    try:
        # Fallback for some environments
        import pybel
        OBABEL_AVAILABLE = True
    except ImportError:
        OBABEL_AVAILABLE = False

PLUGIN_NAME = "OpenBabel Conversion Tool"
PLUGIN_VERSION = "2026.01.10"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Import various chemical file formats using OpenBabel with multi-molecule support."
PLUGIN_DEPENDENCIES = ["openbabel"]

def initialize(context):
    """
    Initialize the OpenBabel Conversion Tool plugin.
    """
    if not OBABEL_AVAILABLE:
        print(f"[{PLUGIN_NAME}] OpenBabel (pybel) not found. Plugin disabled.")
        return

    # Native formats to potentially exclude from automatic handling if we don't want to override
    # However, since we offer a dialog for SDF, we might want to handle it.
    # The user manual implies we can register file openers. 
    # If we register .sdf, we might override the native one.
    # Given the user request for an SDF dialog, we WILL register .sdf
    
    # We will prioritize formats that OpenBabel supports.
    # We'll also add a generic "Import via OpenBabel" menu.
    
    # helper for opening files
    def open_file_wrapper(path):
        open_file_with_openbabel(path, context)

    # 1. Register File Openers for all Pybel formats
    # pybel.informats is a dict {ext: description}
    # We exclude common native formats that we DON'T want to override 
    # (e.g. maybe .mol is fine native? But .sdf needs dialog)
    # Let's override .sdf to provide the dialog.
    # We will exclude .mol and .xyz assuming native is good enough for single structures?
    # Actually, to be safe and "all-encompassing", let's register generic handlers 
    # or specific ones.
    # Note: registering too many might clutter if not careful, but moleditpy handles it by extension.
    
    IGNORED_EXTS = {
        'mol', # Native single mol loader is robust
        'xyz', # Native loader has chatty dialog we might not want to replace or maybe we do? Native is fine.
        # 'sdf' -> We WANT to handle this for the dialog.
    }
    
    supported_exts = []
    for ext in pybel.informats:
        if ext not in IGNORED_EXTS:
            supported_exts.append(f".{ext}")
            # Register individual file opener
            # Note: context.register_file_opener replaces existing if any? 
            # The manual doesn't say, but usually last registered wins or it's a list.
            # We'll assume we can register.
            context.register_file_opener(f".{ext}", open_file_wrapper)

    # 2. Register Drop Handler
    def drop_handler(path):
        _, ext = os.path.splitext(path)
        ext = ext.lower().lstrip('.')
        if ext in pybel.informats and ext not in IGNORED_EXTS:
            open_file_with_openbabel(path, context)
            return True # Handled
        # Special case for SDF if we decided to handle it
        if ext == 'sdf':
             open_file_with_openbabel(path, context)
             return True
        return False

    context.register_drop_handler(drop_handler, priority=10)

    # 3. Register Export Action
    def export_wrapper():
        export_with_openbabel(context)
    
    context.add_export_action("Export via OpenBabel...", export_wrapper)


def open_file_with_openbabel(file_path, context):
    """
    Main logic to read file with Pybel, handle multiple molecules, convert to RDKit, and load.
    """
    mw = context.get_main_window()
    
    try:
        _, ext = os.path.splitext(file_path)
        file_fmt = ext.lower().lstrip('.')
        
        # Determine format for Pybel
        if not file_fmt:
            # try to guess or ask? 
            # Pybel readfile needs format.
            QMessageBox.warning(mw, "Error", "Cannot determine file format (no extension).")
            return

        # Handle non-ASCII paths (Windows/OpenBabel issue workaround)
        # If the path contains non-ascii characters, OpenBabel might fail to open it.
        # We assume safe-ascii workaround by copying to a temp file.
        temp_file_path = None
        path_to_open = file_path
        
        needs_temp = False
        try:
             file_path.encode('ascii')
        except UnicodeEncodeError:
             needs_temp = True
             
        if needs_temp:
            import tempfile
            import shutil
            try:
                # Create a temp file with the same extension
                # handling suffix correctly for Pybel format guessing (though we pass format explicitly)
                fd, temp_file_path = tempfile.mkstemp(suffix=f".{ext.lstrip('.')}")
                os.close(fd)
                shutil.copy2(file_path, temp_file_path)
                path_to_open = temp_file_path
            except Exception as e:
                print(f"[{PLUGIN_NAME}] Failed to create temp file for non-ascii path: {e}")
                # Fallback to original path
                path_to_open = file_path
                temp_file_path = None

        try:
            # Attempt to read file
            # pybel.readfile returns a generator
            mols = list(pybel.readfile(file_fmt, path_to_open))
        finally:
            # Cleanup temp file
            if temp_file_path and os.path.exists(temp_file_path):
                try:
                    os.remove(temp_file_path)
                except Exception:
                   pass
        
        if not mols:
            QMessageBox.warning(mw, "Error", "No molecules found in file.")
            return

        selected_mol = None

        if len(mols) == 1:
            selected_mol = mols[0]
        else:
            # Multiple molecules -> Show Dialog
            dialog = MoleculeSelectionDialog(mols, mw)
            if dialog.exec() == QDialog.DialogCode.Accepted:
                selected_index = dialog.get_selected_index()
                if selected_index is not None and 0 <= selected_index < len(mols):
                    selected_mol = mols[selected_index]
                else:
                    return # User selected nothing valid
            else:
                return # User cancelled

        if not selected_mol:
            return

        # Convert to RDKit
        # We use the MolBlock format (MDL MOL) as the bridge
        try:
            mol_block = selected_mol.write("mol")
        except Exception as e:
            QMessageBox.critical(mw, "Conversion Error", f"Failed to convert OpenBabel molecule to MolBlock:\n{e}")
            return

        rd_mol = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=True)
        
        if rd_mol:
            # Post-processing similar to native loader
            if rd_mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(rd_mol)
            
            # Assign Stereochemistry
            AllChem.AssignStereochemistry(rd_mol, cleanIt=True, force=True)
            if rd_mol.GetNumConformers() > 0:
                AllChem.WedgeMolBonds(rd_mol, rd_mol.GetConformer())
            
            # Load into Main Window
            context.current_molecule = rd_mol
            
            # Push to undo stack
            mw.push_undo_state()
            
            # Update View
            # We need to manually reset camera or ensure view is updated
            # context.current_molecule setter triggers redraw, but maybe not camera reset
            if hasattr(mw, 'plotter'):
                mw.plotter.reset_camera()
            
            # Also update 2D view if possible or fit to view
            # Using QTimer to allow UI to settle
            QTimer.singleShot(100, lambda: getattr(mw, 'fit_to_view', lambda: None)())

            # Switch to 3D only mode
            # We access the internal method on MainWindow which proxies to the UI Manager
            # We prioritize the internal name `_enter_3d_viewer_ui_mode` as verified in source.
            try:
                if hasattr(mw, '_enter_3d_viewer_ui_mode'):
                     mw._enter_3d_viewer_ui_mode()
                elif hasattr(mw, 'enter_3d_viewer_ui_mode'):
                     mw.enter_3d_viewer_ui_mode()
            except Exception as e:
                print(f"[{PLUGIN_NAME}] Failed to switch to 3D mode: {e}")

            mw.statusBar().showMessage(f"Loaded {file_path} via OpenBabel")

            # Set current file path so "Save" works and title updates
            if hasattr(mw, 'current_file_path'):
                 mw.current_file_path = file_path
            if hasattr(mw, 'has_unsaved_changes'):
                 mw.has_unsaved_changes = False
            if hasattr(mw, 'update_window_title'):
                 mw.update_window_title()

        else:
            QMessageBox.critical(mw, "Error", "Failed to create RDKit molecule from converted block.")

    except Exception as e:
        QMessageBox.critical(mw, "Import Error", f"An error occurred:\n{e}")
        import traceback
        traceback.print_exc()


def export_with_openbabel(context):
    """
    Export the current molecule using OpenBabel.
    """
    mw = context.get_main_window()
    if not OBABEL_AVAILABLE:
        QMessageBox.warning(mw, "Error", "OpenBabel is not available.")
        return
        
    mol = context.current_molecule
    if not mol:
        QMessageBox.warning(mw, "Export Error", "No molecule loaded to export.")
        return

    # Generate filter string from pybel.outformats
    # Format: "Description (*.ext);;..."
    # We sort them for easier reading
    formats = []
    for ext, desc in pybel.outformats.items():
        formats.append(f"{desc} (*.{ext})")
    formats.sort()
    
    # Add All Files option
    filter_str = "All Files (*);;" + ";;".join(formats)

    # Use QFileDialog to get path
    path, chosen_filter = QFileDialog.getSaveFileName(mw, "Export via OpenBabel", "", filter_str)
    if not path:
        return

    try:
        # Determine format from extension
        _, ext = os.path.splitext(path)
        fmt = ext.lower().lstrip('.')
        
        if not fmt:
             QMessageBox.warning(mw, "Error", "Please specify a file extension (e.g., .gjf, .xyz, .pdb).")
             return
             
        if fmt not in pybel.outformats:
             QMessageBox.warning(mw, "Error", f"Format '{fmt}' is not supported by OpenBabel.")
             return

        # RDKit -> MolBlock
        try:
            mol_block = Chem.MolToMolBlock(mol)
        except Exception as e:
            QMessageBox.critical(mw, "Conversion Error", f"Failed to convert RDKit molecule to MolBlock:\n{e}")
            return
        
        # MolBlock -> Pybel
        try:
            pybel_mol = pybel.readstring("mol", mol_block)
        except Exception as e:
            QMessageBox.critical(mw, "Conversion Error", f"Failed to read MolBlock with OpenBabel:\n{e}")
            return
        
        # Save
        pybel_mol.write(fmt, path, overwrite=True)
        
        mw.statusBar().showMessage(f"Exported to {path}")
        
    except Exception as e:
        QMessageBox.critical(mw, "Export Error", f"Failed to export:\n{e}")
        import traceback
        traceback.print_exc()


class MoleculeSelectionDialog(QDialog):
    def __init__(self, mols, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Molecule")
        self.resize(400, 300)
        
        layout = QVBoxLayout(self)
        
        lbl = QLabel(f"Found {len(mols)} molecules. Please select one:", self)
        layout.addWidget(lbl)
        
        self.list_widget = QListWidget(self)
        self.list_widget.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        
        for i, m in enumerate(mols):
            title = m.title.strip()
            name = title if title else f"Molecule {i+1}"
            formula = m.formula
            self.list_widget.addItem(f"{i+1}: {name} ({formula})")
        
        layout.addWidget(self.list_widget)
        
        # Select first by default
        self.list_widget.setCurrentRow(0)
        
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        
    def get_selected_index(self):
        rows = self.list_widget.selectedIndexes()
        if rows:
            return rows[0].row()
        return None
