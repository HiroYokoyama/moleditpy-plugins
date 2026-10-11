import os
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (
    QApplication,
    QProgressDialog,
    QFileDialog,
    QDockWidget,
    QWidget,
    QVBoxLayout,
    QLabel,
    QHBoxLayout,
    QPushButton,
    QMessageBox,
    QDoubleSpinBox,
    QComboBox,
    QCheckBox,
    QDialog,
    QLineEdit,
    QFormLayout,
    QDialogButtonBox,
)
from PyQt6.QtCore import Qt, QTimer
import logging

# --- Dependency Management ---
try:
    from rdkit import Chem
    from rdkit import Geometry
    from rdkit.Chem import rdDetermineBonds
except ImportError:
    Chem = None
    Geometry = None
    rdDetermineBonds = None

# Periodic Table for Radii
try:
    from moleditpy.utils.constants import pt
except ImportError:
    try:
        from modules.constants import pt
    except ImportError:
        try:
            from rdkit import Chem

            pt = Chem.GetPeriodicTable()
        except (ImportError, RuntimeError, AttributeError, ValueError):
            pt = None

__author__ = "HiroYokoyama"
PLUGIN_AUTHOR = __author__
PLUGIN_NAME = "Mapped Cube Viewer"
PLUGIN_VERSION = "2026.10.13"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_DESCRIPTION = "Visualizes a property (e.g. ESP) mapped onto an isosurface (e.g. electron density) from two Cube files."

# --- Core Logic: Robust Parser from cube_viewer.py ---


def parse_cube_data(filename):
    """Parse Cube records using the PySCF Calculator's strict parser logic.

    Kept self-contained because each viewer is installed independently.
    """
    def number(token):
        return float(token.replace("D", "E").replace("d", "e"))
    if not filename or not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    if os.path.getsize(filename) > 256 * 1024 * 1024:
        raise ValueError("Cube file exceeds the 256 MiB loading limit.")

    with open(filename, "r", encoding="utf-8", errors="strict") as f:
        lines = f.readlines()

    if len(lines) < 6:
        raise ValueError("File too short to be a Cube file.")

    # --- Header Parsing ---
    try:
        # Line 3: Natoms, Origin
        tokens = lines[2].split()
        if len(tokens) < 4:
            raise ValueError("Invalid Origin line")
        n_atoms_raw = int(tokens[0])
        n_atoms = abs(n_atoms_raw)
        origin_raw = np.array([number(tokens[1]), number(tokens[2]), number(tokens[3])])

        def parse_vec(line):
            t = line.split()
            if len(t) < 4:
                raise ValueError("Invalid Axis line")
            return int(t[0]), np.array([number(t[1]), number(t[2]), number(t[3])])

        nx, x_vec_raw = parse_vec(lines[3])
        ny, y_vec_raw = parse_vec(lines[4])
        nz, z_vec_raw = parse_vec(lines[5])

        is_angstrom_header = nx < 0 or ny < 0 or nz < 0
        nx, ny, nz = abs(nx), abs(ny), abs(nz)
        if not all(np.all(np.isfinite(v)) for v in (origin_raw, x_vec_raw, y_vec_raw, z_vec_raw)):
            raise ValueError("Cube header contains nonfinite coordinates.")
        if np.linalg.matrix_rank(np.array([x_vec_raw, y_vec_raw, z_vec_raw])) != 3:
            raise ValueError("Cube axis vectors must span three dimensions.")
        if min(nx, ny, nz) == 0 or nx * ny * nz > 2_000_000:
            raise ValueError("Cube grid must contain 1 to 2,000,000 points.")

    except (ValueError, IndexError) as e:
        raise ValueError(f"Header parsing failed: {e}") from e

    # --- Atoms Parsing ---
    atoms = []
    current_line = 6

    for _ in range(n_atoms):
        if current_line >= len(lines):
            raise ValueError("Incomplete cube atom block")
        line = lines[current_line].split()
        current_line += 1

        try:
            if len(line) < 5:
                # Potentially empty line or malformed
                raise ValueError("Malformed cube atom line")
            atomic_num = int(line[0])
            if not 0 <= atomic_num <= 118:
                raise ValueError("Invalid cube atomic number")
            x, y, z = number(line[2]), number(line[3]), number(line[4])
            if not np.all(np.isfinite([number(line[1]), x, y, z])):
                raise ValueError("Nonfinite cube atom coordinates")
            atoms.append((atomic_num, np.array([x, y, z])))
        except (ValueError, IndexError) as exc:
            raise ValueError("Malformed cube atom line") from exc

    # A negative atom count means the cube holds several data sets. The
    # DSET_IDS block (count + that many ids, possibly wrapped) sits *after*
    # the atom lines -- skipping a line before them consumed the first atom.
    n_datasets = int(tokens[4]) if len(tokens) > 4 else 1
    if n_atoms_raw < 0:
        if n_datasets != 1:
            raise ValueError("Negative NATOMS requires NVAL absent or equal to one.")
        try:
            parts = lines[current_line].split()
            n_datasets = int(parts[0])
            if n_datasets < 1:
                raise ValueError("Cube dataset count must be positive")
            identifiers = [int(value) for value in parts[1:]]
            current_line += 1
            while len(identifiers) < n_datasets and current_line < len(lines):
                identifiers.extend(int(value) for value in lines[current_line].split())
                current_line += 1
            if len(identifiers) != n_datasets:
                raise ValueError("Incomplete cube dataset identifiers")
        except (ValueError, IndexError) as exc:
            raise ValueError("Malformed cube dataset identifiers") from exc
    if n_datasets < 1 or nx * ny * nz * n_datasets > 2_000_000:
        raise ValueError("Cube dataset size exceeds the loading limit.")

    # --- Volumetric Data Parsing ---
    # Find start of data
    while current_line < len(lines):
        line_content = lines[current_line].strip()
        parts = line_content.split()
        if not parts:
            current_line += 1
            continue

        break

    if current_line >= len(lines):
        # Allow header-only validation if explicitly requested?
        # But for 'data', we need data.
        # Fallback for empty data
        raise ValueError("Cube contains no volumetric data")
    else:
        # Validate every token; whitespace wrapping is independent of grid shape.
        tokens = " ".join(lines[current_line:]).split()
        try:
            data_values = np.fromiter((number(token) for token in tokens), dtype=float)
        except ValueError as exc:
            raise ValueError("Malformed cube volumetric data") from exc

    n_points = nx * ny * nz
    expected_size = n_points * n_datasets
    actual_size = len(data_values)

    if actual_size != expected_size:
        raise ValueError(f"Cube data size mismatch: expected {expected_size}, found {actual_size}")
    if not np.all(np.isfinite(data_values)):
        raise ValueError("Cube contains nonfinite volumetric values")

    if n_datasets > 1:
        # Values are interleaved point by point; the first n_points of the
        # mixed stream belong to no single orbital.
        data_values = data_values[0 : n_points * n_datasets : n_datasets]

    return {
        "atoms": atoms,
        "origin": origin_raw,
        "x_vec": x_vec_raw,
        "y_vec": y_vec_raw,
        "z_vec": z_vec_raw,
        "dims": (nx, ny, nz),
        "data_flat": data_values,
        "is_angstrom_header": is_angstrom_header,
    }


def build_grid_from_meta(meta):
    """
    Reconstructs the PyVista grid.
    Correctly handles standard Cube (Z-fast) mapping.
    """
    nx, ny, nz = meta["dims"]
    origin = meta["origin"].copy()
    x_vec = meta["x_vec"].copy()
    y_vec = meta["y_vec"].copy()
    z_vec = meta["z_vec"].copy()
    atoms = []

    BOHR_TO_ANGSTROM = 0.529177210903
    convert_to_ang = True
    if meta["is_angstrom_header"]:
        convert_to_ang = False

    if convert_to_ang:
        origin *= BOHR_TO_ANGSTROM
        x_vec *= BOHR_TO_ANGSTROM
        y_vec *= BOHR_TO_ANGSTROM
        z_vec *= BOHR_TO_ANGSTROM

    for anum, pos in meta["atoms"]:
        p = pos.copy()
        if convert_to_ang:
            p *= BOHR_TO_ANGSTROM
        atoms.append((anum, p))

    # Grid Points Generation
    x_range = np.arange(nx)
    y_range = np.arange(ny)
    z_range = np.arange(nz)

    gx, gy, gz = np.meshgrid(x_range, y_range, z_range, indexing="ij")

    # Flatten using Fortran order (X-fast) for VTK Structure
    # !!! IMPORTANT FIX: Must match Flatten order of data !!!
    gx_f = gx.flatten(order="F")
    gy_f = gy.flatten(order="F")
    gz_f = gz.flatten(order="F")

    points = (
        origin + np.outer(gx_f, x_vec) + np.outer(gy_f, y_vec) + np.outer(gz_f, z_vec)
    )

    grid = pv.StructuredGrid()
    grid.points = points
    grid.dimensions = [nx, ny, nz]

    # Data Mapping
    raw_data = meta["data_flat"]

    # Standard Cube: Z-Fast (C-order reshape)
    # Then flatten F-order to match the X-fast points we just generated
    vol_3d = raw_data.reshape((nx, ny, nz), order="C")

    # Key name "property_values" for Mapped Viewer
    grid.point_data["property_values"] = vol_3d.flatten(order="F")

    return {"atoms": atoms}, grid


def read_cube(filename):
    meta = parse_cube_data(filename)
    return build_grid_from_meta(meta)


# --- Dialog & Widget (Identical Logic, just ensured imports) ---


class MappedCubeSetupDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Mapped Cube Setup")
        self.surface_file = None
        self.property_file = None

        layout = QVBoxLayout()
        self.setLayout(layout)

        form = QFormLayout()

        # File 1: Surface
        self.le_surf = QLineEdit()
        btn_surf = QPushButton("Browse")
        btn_surf.clicked.connect(self.browse_surf)
        h1 = QHBoxLayout()
        h1.addWidget(self.le_surf)
        h1.addWidget(btn_surf)
        form.addRow("Surface (Density):", h1)

        # File 2: Property
        self.le_prop = QLineEdit()
        btn_prop = QPushButton("Browse")
        btn_prop.clicked.connect(self.browse_prop)
        h2 = QHBoxLayout()
        h2.addWidget(self.le_prop)
        h2.addWidget(btn_prop)
        form.addRow("Property (Color):", h2)

        layout.addLayout(form)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def browse_surf(self):
        f, _ = QFileDialog.getOpenFileName(
            self, "Select Surface Cube", "", "Cube (*.cube *.cub)"
        )
        if f:
            self.le_surf.setText(f)

    def browse_prop(self):
        f, _ = QFileDialog.getOpenFileName(
            self, "Select Property Cube", "", "Cube (*.cube *.cub)"
        )
        if f:
            self.le_prop.setText(f)

    def accept(self):
        self.surface_file = self.le_surf.text()
        self.property_file = self.le_prop.text()
        if not os.path.exists(self.surface_file):
            QMessageBox.warning(self, "Error", "Surface file not found.")
            return
        if not os.path.exists(self.property_file):
            QMessageBox.warning(self, "Error", "Property file not found.")
            return
        super().accept()


class MappedWidget(QWidget):
    def __init__(self, context, dock, grid_surf, grid_prop):
        super().__init__(context.get_main_window())
        self.context = context
        # [DIRECT ACCESS] to main window for legacy plotter/view access
        self.mw = context.get_main_window()
        self.dock = dock
        self.grid_surf = grid_surf
        self.grid_prop = grid_prop
        self.actor = None
        self._bound_mol = None
        self._structure_watch_timer = QTimer(self)
        self._structure_watch_timer.timeout.connect(self._check_structure_still_bound)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.warning_label = QLabel()
        self.warning_label.setWordWrap(True)
        self.layout.addWidget(self.warning_label)

        # Controls
        self.layout.addWidget(QLabel("<b>Surface Isovalue:</b>"))
        self.iso_spin = QDoubleSpinBox()
        self.iso_spin.setRange(-1000, 1000)
        self.iso_spin.setDecimals(5)
        self.iso_spin.setSingleStep(0.001)

        # Default ISO
        vals = self.grid_surf.point_data.get("property_values", np.array([0]))
        # Robust default calculation
        if len(vals) > 0:
            if vals.max() > 0.1:
                default_iso = 0.00200
            else:
                default_iso = float(np.mean(vals))
        else:
            default_iso = 0.002

        self.iso_spin.setValue(default_iso)
        self.iso_spin.setKeyboardTracking(False)
        self.iso_spin.valueChanged.connect(lambda: self.update_mesh(auto_fit=False))
        self.layout.addWidget(self.iso_spin)

        # Color Range
        self.layout.addWidget(QLabel("<b>Property Color Range:</b>"))

        cbox = QHBoxLayout()
        cbox.addWidget(QLabel("Style:"))
        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems(
            [
                "jet_r",
                "jet",
                "rainbow_r",
                "rainbow",
                "bwr_r",
                "bwr",
                "seismic_r",
                "seismic",
                "coolwarm_r",
                "coolwarm",
                "viridis_r",
                "viridis",
                "plasma_r",
                "plasma",
                "magma_r",
                "magma",
                "inferno_r",
                "inferno",
            ]
        )
        self.cmap_combo.setCurrentText("jet_r")
        self.cmap_combo.currentTextChanged.connect(
            lambda: self.update_mesh(auto_fit=False)
        )
        cbox.addWidget(self.cmap_combo)
        self.layout.addLayout(cbox)

        rbox = QHBoxLayout()

        pvals = self.grid_prop.point_data.get("property_values", np.array([-0.1, 0.1]))
        if len(pvals) > 0:
            vmin, vmax = pvals.min(), pvals.max()
        else:
            vmin, vmax = -0.1, 0.1

        self.min_spin = QDoubleSpinBox()
        self.min_spin.setRange(-1e20, 1e20)
        self.min_spin.setDecimals(6)
        self.min_spin.setSingleStep(0.01)
        self.min_spin.setKeyboardTracking(False)
        self.min_spin.setValue(vmin)
        self.min_spin.valueChanged.connect(lambda: self.update_mesh(auto_fit=False))

        self.max_spin = QDoubleSpinBox()
        self.max_spin.setRange(-1e20, 1e20)
        self.max_spin.setDecimals(6)
        self.max_spin.setSingleStep(0.01)
        self.max_spin.setKeyboardTracking(False)
        self.max_spin.setValue(vmax)
        self.max_spin.valueChanged.connect(lambda: self.update_mesh(auto_fit=False))

        rbox.addWidget(QLabel("Min:"))
        rbox.addWidget(self.min_spin)
        rbox.addWidget(QLabel("Max:"))
        rbox.addWidget(self.max_spin)
        self.layout.addLayout(rbox)

        btn_fit = QPushButton("Fit Range to Surface")
        btn_fit.clicked.connect(lambda: self.update_mesh(auto_fit=True))
        self.layout.addWidget(btn_fit)

        # Opacity
        self.layout.addWidget(QLabel("<b>Opacity:</b>"))
        self.opacity_spin = QDoubleSpinBox()
        self.opacity_spin.setRange(0, 1)
        self.opacity_spin.setSingleStep(0.1)
        self.opacity_spin.setValue(0.4)
        self.opacity_spin.valueChanged.connect(lambda: self.update_mesh(auto_fit=False))
        self.layout.addWidget(self.opacity_spin)

        # Exports
        ebox = QHBoxLayout()
        btn_view = QPushButton("Export View")
        btn_view.clicked.connect(self.export_view)
        ebox.addWidget(btn_view)

        btn_cbar = QPushButton("Export Color Bar")
        btn_cbar.clicked.connect(self.export_colorbar)
        ebox.addWidget(btn_cbar)
        self.layout.addLayout(ebox)

        # Options
        self.check_transparent = QCheckBox("Transparent Background")
        self.check_transparent.setChecked(True)
        self.layout.addWidget(self.check_transparent)

        self.layout.addStretch()
        btn = QPushButton("Close")
        btn.clicked.connect(self.close_plugin)
        self.layout.addWidget(btn)

        self.update_mesh(auto_fit=True)

    def export_view(self):
        f, _ = QFileDialog.getSaveFileName(
            self, "Save View", "mapped_view.png", "PNG Image (*.png)"
        )
        if f:
            try:
                self.context.plotter.screenshot(
                    f, transparent_background=self.check_transparent.isChecked()
                )
                QMessageBox.information(self, "Success", f"Saved View to {f}")
            except (OSError, RuntimeError, AttributeError, KeyError, ValueError) as e:
                QMessageBox.critical(self, "Error", f"Save failed: {e}")

    def export_colorbar(self):
        f, _ = QFileDialog.getSaveFileName(
            self, "Save Color Bar", "mapped_colorbar.png", "PNG Image (*.png)"
        )
        if f:
            pl = None
            try:
                vmin = self.min_spin.value()
                vmax = self.max_spin.value()
                if vmin >= vmax:
                    raise ValueError("The color minimum must be less than the maximum.")
                pl = pv.Plotter(off_screen=True, window_size=(2000, 400))
                if self.check_transparent.isChecked():
                    pl.set_background(None)
                else:
                    pl.set_background("white")

                cmap = self.cmap_combo.currentText()

                mesh = pv.Box()
                mesh.point_data["scalars"] = np.linspace(vmin, vmax, mesh.n_points)

                pl.add_mesh(
                    mesh,
                    scalars="scalars",
                    cmap=cmap,
                    clim=[vmin, vmax],
                    show_scalar_bar=True,
                    opacity=0.0,
                    scalar_bar_args={
                        "title": "",
                        "vertical": False,
                        "n_labels": 5,
                        "italic": False,
                        "bold": False,
                        "title_font_size": 1,
                        "label_font_size": 50,
                        "color": "black",
                        "height": 0.4,
                        "width": 0.8,
                        "position_x": 0.1,
                        "position_y": 0.3,
                    },
                )

                if hasattr(pl, "scalar_bar"):
                    try:
                        sb = pl.scalar_bar
                        sb.GetLabelTextProperty().SetJustificationToCentered()
                    except Exception as _e:
                        logging.warning("[mapped_cube_viewer.py:430] silenced: %s", _e)

                pl.screenshot(
                    f, transparent_background=self.check_transparent.isChecked()
                )
                QMessageBox.information(self, "Success", f"Saved Color Bar to {f}")
            except Exception as e:
                logging.exception("Save failed: %s", e)
                QMessageBox.critical(self, "Error", f"Save failed: {e}")
            finally:
                if pl is not None:
                    pl.close()

    def update_mesh(self, auto_fit=False):
        try:
            # As in PySCF's MappedVisualizer, remove the old actor before
            # contouring: an empty or failed update must leave no stale surface.
            self._remove_surface()
            self.warning_label.clear()
            iso_val = self.iso_spin.value()
            opacity = self.opacity_spin.value()

            iso = self.grid_surf.contour([iso_val], scalars="property_values")
            if iso.n_points == 0:
                self.warning_label.setText("No surface at this isovalue.")
                self.context.plotter.render()
                return

            mapped = iso.sample(self.grid_prop)
            valid = mapped.point_data.get("vtkValidPointMask")
            if valid is None or not np.all(np.asarray(valid) == 1):
                raise ValueError("The property grid does not cover the whole surface. Choose matching Cube files or a larger property grid.")
            mvals = np.asarray(mapped.point_data["property_values"])
            if not np.all(np.isfinite(mvals)):
                raise ValueError("Mapped surface contains nonfinite property values.")

            clim = [self.min_spin.value(), self.max_spin.value()]
            if clim[0] >= clim[1] and not auto_fit:
                raise ValueError("The color minimum must be less than the maximum.")
            if auto_fit:
                mvals = mapped.point_data.get("property_values", None)
                if mvals is not None and len(mvals) > 0:
                    vmin, vmax = mvals.min(), mvals.max()
                    if vmax - vmin < 1e-9:
                        vmax += 0.001
                    clim = [vmin, vmax]

                    self.min_spin.blockSignals(True)
                    self.min_spin.setValue(vmin)
                    self.min_spin.blockSignals(False)
                    self.max_spin.blockSignals(True)
                    self.max_spin.setValue(vmax)
                    self.max_spin.blockSignals(False)

            self.actor = self.context.plotter.add_mesh(
                mapped,
                scalars="property_values",
                cmap=self.cmap_combo.currentText(),
                clim=clim,
                smooth_shading=True,
                opacity=opacity,
                name="mapped_mesh_dock",
                scalar_bar_args={"title": ""},
            )

            self.context.plotter.render()

        except Exception as e:
            logging.warning("Update Error: %s", e)
            self.warning_label.setText(str(e))
            self.context.plotter.render()

    def _remove_surface(self):
        if self.actor is not None:
            self.context.plotter.remove_actor(self.actor)
            self.actor = None

    def bind_structure(self, mol):
        self._bound_mol = mol
        self._structure_watch_timer.start(1500)

    def _check_structure_still_bound(self):
        if self._bound_mol is not None and self.context.current_molecule is not self._bound_mol:
            self._structure_watch_timer.stop()
            self.close_plugin()

    def closeEvent(self, event):
        self._structure_watch_timer.stop()
        try:
            self._remove_surface()
        except (RuntimeError, AttributeError, ValueError) as exc:
            logging.warning("Mapped surface cleanup failed: %s", exc)
        super().closeEvent(event)

    def close_plugin(self):
        self._structure_watch_timer.stop()
        try:
            self._remove_surface()
            self.context.plotter.render()
        except (RuntimeError, AttributeError, KeyError, ValueError) as _e:
            logging.warning("[mapped_cube_viewer.py:485] silenced: %s", _e)
        # Restore 2D editing UI when leaving the 3D viewer mode
        try:
            owns_structure = self._bound_mol is None or self.context.current_molecule is self._bound_mol
            if owns_structure and hasattr(self.mw, "ui_manager") and hasattr(
                self.mw.ui_manager, "restore_ui_for_editing"
            ):
                self.mw.ui_manager.restore_ui_for_editing()
        except Exception as _e:
            logging.warning("[mapped_cube_viewer.py:490] silenced: %s", _e)
        if self.dock:
            self.mw.removeDockWidget(self.dock)
            self.dock.deleteLater()
            self.dock = None
        self.close()
        self.deleteLater()


# --- Entry Point ---
def run_plugin(context):
    mw = context.get_main_window()
    if Chem is None:
        QMessageBox.critical(mw, "Mapped Cube Viewer", "RDKit is required for this plugin.")
        return
    dlg = MappedCubeSetupDialog(mw)
    if dlg.exec() != QDialog.DialogCode.Accepted:
        return

    s_file = dlg.surface_file
    p_file = dlg.property_file

    if not s_file or not p_file:
        return

    progress = None
    try:
        progress = QProgressDialog("Reading Cube Files...", "Cancel", 0, 2, mw)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.show()

        # Using the robust read_cube logic
        meta1, grid_surf = read_cube(s_file)
        progress.setValue(1)
        QApplication.processEvents()
        if progress.wasCanceled():
            return

        meta2, grid_prop = read_cube(p_file)
        if progress.wasCanceled():
            return
        progress.setValue(2)
        QApplication.processEvents()
        progress.close()

        # Visualize Atoms
        atoms = meta1["atoms"]
        if Chem:
            s = f"{len(atoms)}\n\n"
            for n, p in atoms:
                sym = pt.GetElementSymbol(n) if pt else "C"
                s += f"{sym} {p[0]} {p[1]} {p[2]}\n"
            mol = Chem.MolFromXYZBlock(s)
            try:
                rdDetermineBonds.DetermineConnectivity(mol)
            except Exception as _e:  # noqa: BLE001 - bond perception is optional here
                # Connectivity is cosmetic; failing it must not stop the viewer
                # from opening, so this stays broad on purpose.
                logging.warning("[mapped_cube_viewer.py:538] silenced: %s", _e)

        # Close existing dock via context
        old_dock = context.get_window("main_panel")
        if old_dock:
            try:
                widget = old_dock.widget()
                if hasattr(widget, "close_plugin"):
                    widget.close_plugin()
                else:
                    old_dock.close()
            except (RuntimeError, AttributeError) as _e:
                logging.warning("[mapped_cube_viewer.py:559] silenced: %s", _e)

        # The old viewer may clear molecular state when it closes. Install the
        # new structure only after that cleanup completes.
        context.enter_3d_viewer_mode()
        context.current_molecule = mol

        dock = QDockWidget("Mapped Viewer", mw)
        dock.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea)
        widget = MappedWidget(context, dock, grid_surf, grid_prop)
        widget.bind_structure(mol)
        dock.setWidget(widget)

        # [DIRECT ACCESS] for dock injection
        mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)

        # Register for management
        context.register_window("main_panel", dock)

    except Exception as e:
        logging.exception("Mapped Cube Viewer Error: %s", e)
        context.show_status_message(f"Mapped Cube Viewer Error: {e}")
    finally:
        if progress is not None:
            progress.close()


def run(mw):
    if not hasattr(mw, "plugin_manager"):
        return

    from moleditpy.plugins.plugin_interface import PluginContext

    context = PluginContext(mw.plugin_manager, PLUGIN_NAME)
    if not context:
        context = PluginContext(mw.plugin_manager, PLUGIN_NAME)

    run_plugin(context)
