import os
import sys
import importlib
import json
import tempfile
import shutil
from types import SimpleNamespace
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QListWidget,
    QDialogButtonBox,
    QLabel,
    QMessageBox,
    QFileDialog,
    QAbstractItemView,
    QProgressDialog,
)
from PyQt6.QtCore import QTimer, QProcess, QEventLoop, Qt
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

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
PLUGIN_VERSION = "2026.10.13"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Import various chemical file formats using OpenBabel with multi-molecule support."
)
PLUGIN_DEPENDENCIES = ["openbabel-wheel"]

# MoleditPy 4.11.0 lists multi-molecule SDF files in its own selector.
HOST_SDF_SELECTOR_VERSION = (4, 11, 0)



# Native parsers run in a disposable process: a Python thread cannot stop a
# hung OpenBabel C++ call, and may still hold the interpreter lock.
MAX_FILE_BYTES = 32 * 1024 * 1024
MAX_MOLECULES = 500
MAX_ATOMS = 10000
CONVERSION_TIMEOUT_MS = 30000
_CONVERSION_SCRIPT = r"""
import sys, json
try:
    from openbabel import pybel
except ImportError:
    import pybel
mode, fmt, source, output = sys.argv[1:5]
if mode == 'import':
    records = []
    size = 0
    for mol in pybel.readfile(fmt, source):
        if len(records) >= 500:
            raise ValueError('File contains more than 500 molecules; split it first.')
        if mol.OBMol.NumAtoms() > 10000:
            raise ValueError('Molecule exceeds 10000 atoms.')
        block = mol.write('mol')
        size += len(block.encode('utf-8'))
        if size > 32 * 1024 * 1024:
            raise ValueError('Converted structures exceed 32 MiB; split the file first.')
        records.append(dict(title=mol.title[:256], formula=mol.formula[:256], block=block))
    with open(output, 'w', encoding='utf-8') as handle:
        json.dump(records, handle)
else:
    with open(source, encoding='utf-8') as handle:
        mol = pybel.readstring('mol', handle.read())
    mol.write(fmt, output, overwrite=True)
"""


def _run_conversion(script, arguments, parent):
    process = QProcess(parent)
    loop = QEventLoop()
    progress = QProgressDialog('Converting with OpenBabel...', 'Cancel', 0, 0, parent)
    progress.setWindowModality(Qt.WindowModality.WindowModal)
    progress.setMinimumDuration(0)
    timer = QTimer()
    timer.setSingleShot(True)
    failure = []
    diagnostics = bytearray()

    def stop(message):
        failure.append(message)
        process.kill()
        loop.quit()

    def drain():
        data = bytes(process.readAllStandardError())
        diagnostics.extend(data[:max(0, 4096 - len(diagnostics))])
        process.readAllStandardOutput()

    process.readyReadStandardError.connect(drain)
    process.readyReadStandardOutput.connect(drain)
    process.finished.connect(loop.quit)
    process.errorOccurred.connect(loop.quit)
    progress.canceled.connect(lambda: stop('Conversion cancelled.'))
    timer.timeout.connect(lambda: stop('Conversion timed out after 30 seconds.'))
    try:
        process.start(sys.executable, ['-c', script, *arguments])
        timer.start(CONVERSION_TIMEOUT_MS)
        progress.show()
        loop.exec()
        drain()
        if failure:
            raise ValueError(failure[0])
        if process.exitStatus() != QProcess.ExitStatus.NormalExit or process.exitCode() != 0 or process.error() == QProcess.ProcessError.FailedToStart:
            raise ValueError('OpenBabel conversion failed: ' + (diagnostics.decode('utf-8', errors='replace') or process.errorString()))
    finally:
        timer.stop()
        if process.state() != QProcess.ProcessState.NotRunning:
            process.kill()
            process.waitForFinished(1000)
        progress.close()
        progress.deleteLater()
        process.deleteLater()


def _read_molecules(file_path, file_fmt, parent):
    if file_fmt not in pybel.informats:
        raise ValueError(f"Format '{file_fmt}' is not supported by OpenBabel.")
    if os.path.getsize(file_path) > MAX_FILE_BYTES:
        raise ValueError('File exceeds 32 MiB; split it into smaller files first.')
    with tempfile.TemporaryDirectory(prefix='moleditpy-obabel-') as directory:
        source = os.path.join(directory, 'input.' + file_fmt)
        output = os.path.join(directory, 'result.json')
        # A bounded copy handles Unicode source paths and files growing since stat.
        with open(file_path, 'rb') as src, open(source, 'wb') as dst:
            data = src.read(MAX_FILE_BYTES + 1)
            if len(data) > MAX_FILE_BYTES:
                raise ValueError('File exceeds 32 MiB.')
            dst.write(data)
        _run_conversion(_CONVERSION_SCRIPT, ['import', file_fmt, source, output], parent)
        with open(output, 'rb') as handle:
            raw = handle.read(MAX_FILE_BYTES * 2 + 1)
        if len(raw) > MAX_FILE_BYTES * 2:
            raise ValueError('Converted data is too large.')
        records = json.loads(raw)
        return [SimpleNamespace(title=r['title'], formula=r['formula'],
                                write=lambda fmt, block=r['block']: block)
                for r in records]


def _host_version():
    """Return the host app version as an int tuple, or None if unknown."""
    try:
        VERSION = importlib.import_module('moleditpy.utils.constants').VERSION
    except ImportError:
        try:
            VERSION = importlib.import_module('moleditpy_linux.utils.constants').VERSION
        except ImportError:
            return None
    try:
        return tuple(int(part) for part in str(VERSION).split(".")[:3])
    except ValueError:
        return None


def host_handles_sdf():
    """True when the host app has its own multi-molecule SDF selector."""
    version = _host_version()
    return version is not None and version >= HOST_SDF_SELECTOR_VERSION


def initialize(context):
    """
    Initialize the OpenBabel Conversion Tool plugin.
    """
    if not OBABEL_AVAILABLE:
        print(f"[{PLUGIN_NAME}] OpenBabel (pybel) not found. Plugin disabled.")
        return

    def open_file_wrapper(path):
        open_file_with_openbabel(path, context)

    # 1. Register file openers for every Pybel input format except the ones
    # the host loads natively. SDF is left to the host from 4.11.0 on, which
    # has its own multi-molecule selector; older hosts only read the first
    # record, so there the plugin keeps handling it.
    IGNORED_EXTS = {"mol", "xyz"}
    if host_handles_sdf():
        IGNORED_EXTS.add("sdf")

    for ext in pybel.informats:
        if ext not in IGNORED_EXTS:
            context.register_file_opener(f".{ext}", open_file_wrapper, priority=-1)

    # 2. Register Drop Handler
    def drop_handler(path):
        _, ext = os.path.splitext(path)
        ext = ext.lower().lstrip(".")
        if ext in pybel.informats and ext not in IGNORED_EXTS:
            open_file_with_openbabel(path, context)
            return True  # Handled
        # Hosts before 4.11.0 cannot pick a structure from a multi-molecule
        # SDF, so claim it there even if Pybel does not list the format.
        if ext == "sdf" and "sdf" not in IGNORED_EXTS:
            open_file_with_openbabel(path, context)
            return True
        return False

    context.register_drop_handler(drop_handler, priority=-1)

    # 3. Register Export Action
    def export_wrapper():
        export_with_openbabel(context)

    context.add_export_action("Export via OpenBabel...", export_wrapper)


def open_file_with_openbabel(file_path, context):
    """
    Main logic to read file with Pybel, handle multiple molecules, convert to RDKit, and load.
    """
    mw = context.get_main_window()
    if getattr(context, '_openbabel_import_active', False) is True:
        context.show_status_message('An OpenBabel import is already running.', 3000)
        return
    context._openbabel_import_active = True
    original_molecule = context.current_molecule

    try:
        _, ext = os.path.splitext(file_path)
        file_fmt = ext.lower().lstrip(".")

        # Determine format for Pybel
        if not file_fmt:
            # try to guess or ask?
            # Pybel readfile needs format.
            QMessageBox.warning(
                mw, "Error", "Cannot determine file format (no extension)."
            )
            return

        mols = _read_molecules(file_path, file_fmt, mw)

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
                    return  # User selected nothing valid
            else:
                return  # User cancelled

        if not selected_mol:
            return

        if context.current_molecule is not original_molecule:
            raise ValueError("Molecule changed during conversion; import cancelled.")

        # Convert to RDKit
        # We use the MolBlock format (MDL MOL) as the bridge
        try:
            mol_block = selected_mol.write("mol")
        except Exception as e:
            QMessageBox.critical(
                mw,
                "Conversion Error",
                f"Failed to convert OpenBabel molecule to MolBlock:\n{e}",
            )
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
            context.push_undo_checkpoint()

            # Update View
            # We need to manually reset camera or ensure view is updated
            # context.current_molecule setter triggers redraw, but maybe not camera reset
            context.reset_3d_camera()

            # Also update 2D view if possible or fit to view
            # Using QTimer to allow UI to settle
            QTimer.singleShot(100, context.fit_2d_view)

            # Switch to 3D only mode
            # We access the internal method on MainWindow which proxies to the UI Manager
            # We prioritize the internal name `enter_3d_viewer_ui_mode` as verified in source.
            try:
                if hasattr(context, "enter_3d_viewer_mode"):
                    context.enter_3d_viewer_mode()
                elif hasattr(mw, "ui_manager") and hasattr(
                    mw.ui_manager, "enter_3d_viewer_ui_mode"
                ):
                    mw.ui_manager.enter_3d_viewer_ui_mode()
            except Exception as e:
                logging.warning("[%s] Failed to switch to 3D mode: %s", PLUGIN_NAME, e)

            context.show_status_message(f"Loaded {file_path} via OpenBabel", 3000)

            # Set current file path so "Save" works and title updates
            if hasattr(mw, "init_manager") and hasattr(
                mw.init_manager, "current_file_path"
            ):
                mw.init_manager.current_file_path = file_path
            if hasattr(mw, "state_manager") and hasattr(
                mw.state_manager, "has_unsaved_changes"
            ):
                mw.state_manager.has_unsaved_changes = False
            context.refresh_ui()

        else:
            QMessageBox.critical(
                mw, "Error", "Failed to create RDKit molecule from converted block."
            )

    except Exception as e:
        logging.exception("Import error: %s", e)
        QMessageBox.critical(mw, "Import Error", f"An error occurred:\n{e}")
    finally:
        context._openbabel_import_active = False


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
    path, chosen_filter = QFileDialog.getSaveFileName(
        mw, "Export via OpenBabel", "", filter_str
    )
    if not path:
        return

    try:
        # Determine format from extension
        _, ext = os.path.splitext(path)
        fmt = ext.lower().lstrip(".")

        if not fmt:
            QMessageBox.warning(
                mw, "Error", "Please specify a file extension (e.g., .gjf, .xyz, .pdb)."
            )
            return

        if fmt not in pybel.outformats:
            QMessageBox.warning(
                mw, "Error", f"Format '{fmt}' is not supported by OpenBabel."
            )
            return

        if mol.GetNumAtoms() > MAX_ATOMS:
            raise ValueError("Molecule exceeds 10000 atoms.")

        # RDKit -> MolBlock
        try:
            mol_block = Chem.MolToMolBlock(mol)
        except Exception as e:  # noqa: BLE001 - any conversion failure is reported to the user
            QMessageBox.critical(
                mw,
                "Conversion Error",
                f"Failed to convert RDKit molecule to MolBlock:\n{e}",
            )
            return

        with tempfile.TemporaryDirectory(prefix='moleditpy-obabel-') as directory:
            source = os.path.join(directory, 'input.mol')
            output = os.path.join(directory, 'output.' + fmt)
            with open(source, 'w', encoding='utf-8') as handle:
                handle.write(mol_block)
            _run_conversion(_CONVERSION_SCRIPT, ['export', fmt, source, output], mw)
            if os.path.getsize(output) > MAX_FILE_BYTES:
                raise ValueError('Exported data exceeds 32 MiB.')
            # Replace only after a successful conversion, also for Unicode paths.
            fd, staged = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), suffix='.tmp')
            try:
                with os.fdopen(fd, 'wb') as handle, open(output, 'rb') as result:
                    shutil.copyfileobj(result, handle)
                os.replace(staged, path)
            finally:
                if os.path.exists(staged):
                    os.remove(staged)

        context.show_status_message(f"Exported to {path}", 3000)

    except Exception as e:
        logging.exception("Export error: %s", e)
        QMessageBox.critical(mw, "Export Error", f"Failed to export:\n{e}")


class MoleculeSelectionDialog(QDialog):
    def __init__(self, mols, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Molecule")
        self.resize(400, 300)

        layout = QVBoxLayout(self)

        lbl = QLabel(f"Found {len(mols)} molecules. Please select one:", self)
        layout.addWidget(lbl)

        self.list_widget = QListWidget(self)
        self.list_widget.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection
        )

        for i, m in enumerate(mols):
            title = m.title.strip()
            name = title if title else f"Molecule {i + 1}"
            formula = m.formula
            self.list_widget.addItem(f"{i + 1}: {name} ({formula})")

        layout.addWidget(self.list_widget)

        # Select first by default
        self.list_widget.setCurrentRow(0)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def get_selected_index(self):
        rows = self.list_widget.selectedIndexes()
        if rows:
            return rows[0].row()
        return None
