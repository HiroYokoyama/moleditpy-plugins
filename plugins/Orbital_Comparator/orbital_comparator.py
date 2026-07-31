"""Orbital Comparator — side-by-side comparison of cube files.

The molecular counterpart of this plugin compares structures; this one
compares volumetric data. Load any .cube files (molecular orbitals, densities,
ESP maps) and show up to four of them at once over the current structure, each
with its own lobe colours, contour level, opacity and style.

Cubes written by the ORCA Result Analyzer record the grid and margin they were
built with on their second comment line; those are read back and shown, so two
orbitals are never compared at settings that make them incomparable.
"""

import logging
import os

from PyQt6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
from PyQt6.QtCore import Qt

try:
    import numpy as np
except ImportError:  # pragma: no cover - host always ships numpy
    np = None

try:
    import pyvista as pv
except ImportError:  # pragma: no cover - host always ships pyvista
    pv = None


PLUGIN_NAME = "Orbital Comparator"
PLUGIN_VERSION = "2026.07.31"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Load and compare up to four .cube files at once (orbitals, densities, "
    "ESP), each with its own colours, isovalue, opacity and style."
)
PLUGIN_CONTEXT = None

BOHR_TO_ANG = 0.529177249
SLOT_COUNT = 4
STYLES = ["Surface", "Wireframe", "Points"]

#: The Cube File Viewers register .cube/.cub openers at the default 0, so this
#: has to be negative to stay under them: a cube named on the command line, or
#: opened from the File menu, belongs to the viewer whenever one is installed.
#: This is the fallback for when none is.
FILE_OPENER_PRIORITY = -10

#: Distinct lobe colours so four orbitals stay tellable apart.
DEFAULT_COLORS = [
    ("#ff0000", "#0000ff"),
    ("#ff8c00", "#008b8b"),
    ("#9400d3", "#7fbf00"),
    ("#ff1493", "#00bfff"),
]


# ---------------------------------------------------------------------------
# Cube reading
# ---------------------------------------------------------------------------


def read_generation_settings(filepath):
    """Grid settings recorded in a cube's second comment line.

    The ORCA Result Analyzer writes "... v3.13.2 | grid=40 | margin=4.00".
    Anything absent comes back as None so the caller says "unknown" rather
    than inventing a value.
    """
    import re

    info = {"version": None, "grid": None, "margin": None}
    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
            fh.readline()
            stamp = fh.readline()
    except OSError as e:
        logging.warning("Could not read cube header %s: %s", filepath, e)
        return info

    m = re.search(r"Analyzer v(\S+)", stamp)
    if m:
        info["version"] = m.group(1)
    m = re.search(r"grid=(\d+)", stamp)
    if m:
        info["grid"] = int(m.group(1))
    m = re.search(r"margin=([0-9.]+)", stamp)
    if m:
        try:
            info["margin"] = float(m.group(1))
        except ValueError:
            logging.warning("Unparsable margin in %s", filepath)
    return info


def describe_settings(info):
    """One-line summary of what a cube was generated with."""
    grid = info.get("grid")
    margin = info.get("margin")
    grid_s = f"{grid} pts" if grid is not None else "unknown grid"
    margin_s = f"margin {margin:.2f} Bohr" if margin is not None else "unknown margin"
    version = info.get("version")
    return f"{grid_s}, {margin_s}" + (f" — v{version}" if version else "")


def parse_cube(filename):
    """Read a Gaussian cube into dims/origin/vectors/data."""
    with open(filename, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    tokens = lines[2].split()
    n_atoms_raw = int(tokens[0])
    n_atoms = abs(n_atoms_raw)
    origin = np.array([float(x) for x in tokens[1:4]])

    nx, *x_vec = lines[3].split()
    nx = int(nx)
    x_vec = np.array([float(x) for x in x_vec])
    ny, *y_vec = lines[4].split()
    ny = int(ny)
    y_vec = np.array([float(x) for x in y_vec])
    nz, *z_vec = lines[5].split()
    nz = int(nz)
    z_vec = np.array([float(x) for x in z_vec])

    start_line = 6 + n_atoms
    n_datasets = 1
    if n_atoms_raw < 0:
        # A negative atom count flags an MO cube: a DSET_IDS block (a count
        # plus that many ids, possibly wrapped) sits before the data.
        n_ids = int(lines[start_line].split()[0])
        n_datasets = max(1, n_ids)
        consumed = len(lines[start_line].split()) - 1
        start_line += 1
        while consumed < n_ids:
            if start_line >= len(lines):
                raise ValueError("cube ends inside its DSET_IDS block")
            fields = len(lines[start_line].split())
            start_line += 1
            if fields == 0:
                # A blank line consumes no ids; without this the loop spins
                # forever and the whole app stops responding.
                continue
            consumed += fields

    data = np.array(" ".join(lines[start_line:]).split(), dtype=float)

    # Several data sets are interleaved point by point; take the first.
    n_points = abs(nx) * abs(ny) * abs(nz)
    if n_datasets > 1 and len(data) >= n_points * n_datasets:
        data = data[0 : n_points * n_datasets : n_datasets]

    return {
        "dims": (abs(nx), abs(ny), abs(nz)),
        "origin": origin,
        "vectors": (x_vec, y_vec, z_vec),
        "data": data,
        "is_angstrom": nx < 0,
    }


def build_grid(meta):
    """A pyvista StructuredGrid from parse_cube output."""
    nx, ny, nz = meta["dims"]
    xv, yv, zv = meta["vectors"]

    grid = pv.StructuredGrid()
    gx, gy, gz = np.meshgrid(
        np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij"
    )
    gx = gx.flatten(order="F")
    gy = gy.flatten(order="F")
    gz = gz.flatten(order="F")

    # A negative voxel count means the grid is already in Angstrom; treating
    # those as Bohr inflates the isosurface by 1.89x.
    scale = 1.0 if meta.get("is_angstrom") else BOHR_TO_ANG
    points = meta["origin"] + np.outer(gx, xv) + np.outer(gy, yv) + np.outer(gz, zv)
    points = points * scale

    grid.points = points
    grid.dimensions = [nx, ny, nz]
    # Cube data is Z-fastest (C order); grid points are X-fastest (F order).
    grid.point_data["values"] = meta["data"].reshape((nx, ny, nz), order="C").flatten(
        order="F"
    )
    return grid


def contrast_text(hex_c):
    """'black' or 'white', whichever stays readable on *hex_c*."""
    h = str(hex_c).lstrip("#")
    if len(h) != 6:
        return "black"
    try:
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    except ValueError:
        return "black"
    return "black" if (r * 299 + g * 587 + b * 114) / 1000 > 128 else "white"


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------


class CubeSlot:
    """One loaded cube: its widgets, its grid and its actor namespace."""

    def __init__(self, index, owner):
        self.index = index
        self.owner = owner
        self.prefix = f"orb_cmp{index}"
        self.path = None
        self.grid = None
        self.settings = {}

        self.box = QGroupBox(f"Cube {index + 1}: (empty)")
        outer = QVBoxLayout(self.box)

        row1 = QHBoxLayout()
        self.check_on = QCheckBox("Show")
        row1.addWidget(self.check_on)
        self.lbl_info = QLabel("no file")
        row1.addWidget(self.lbl_info, 1)
        outer.addLayout(row1)

        row2 = QHBoxLayout()
        self.btn_p = QPushButton("Pos (+)")
        self.btn_n = QPushButton("Neg (-)")
        self.btn_p.clicked.connect(lambda: owner.pick_color(self, "p"))
        self.btn_n.clicked.connect(lambda: owner.pick_color(self, "n"))
        row2.addWidget(self.btn_p)
        row2.addWidget(self.btn_n)
        outer.addLayout(row2)

        row3 = QHBoxLayout()
        row3.addWidget(QLabel("Iso:"))
        self.spin_iso = QDoubleSpinBox()
        self.spin_iso.setRange(0.0001, 10.0)
        self.spin_iso.setSingleStep(0.005)
        self.spin_iso.setDecimals(4)
        self.spin_iso.setValue(0.02)
        row3.addWidget(self.spin_iso)
        row3.addWidget(QLabel("Opacity:"))
        self.spin_opacity = QDoubleSpinBox()
        self.spin_opacity.setRange(0.0, 1.0)
        self.spin_opacity.setSingleStep(0.1)
        self.spin_opacity.setValue(0.5)
        row3.addWidget(self.spin_opacity)
        outer.addLayout(row3)

        row4 = QHBoxLayout()
        row4.addWidget(QLabel("Style:"))
        self.combo_style = QComboBox()
        self.combo_style.addItems(STYLES)
        row4.addWidget(self.combo_style)
        self.check_smooth = QCheckBox("Smooth")
        self.check_smooth.setChecked(True)
        row4.addWidget(self.check_smooth)
        self.btn_load = QPushButton("Load…")
        self.btn_load.clicked.connect(lambda: owner.load_into(self))
        row4.addWidget(self.btn_load)
        self.btn_clear = QPushButton("Clear")
        self.btn_clear.clicked.connect(lambda: owner.clear_slot(self))
        row4.addWidget(self.btn_clear)
        outer.addLayout(row4)

        self.set_color("p", DEFAULT_COLORS[index][0])
        self.set_color("n", DEFAULT_COLORS[index][1])

    def set_color(self, which, hex_c):
        btn = self.btn_p if which == "p" else self.btn_n
        btn.setStyleSheet(
            f"background-color: {hex_c}; color: {contrast_text(hex_c)}; "
            "font-weight: bold;"
        )

    def color(self, which):
        style = (self.btn_p if which == "p" else self.btn_n).styleSheet()
        if "background-color:" in style:
            return style.split("background-color:")[1].split(";")[0].strip()
        return "#ff0000" if which == "p" else "#0000ff"

    def is_on(self):
        return self.check_on.isChecked() and self.grid is not None

    def describe(self):
        if not self.path:
            return "no file"
        return f"{os.path.basename(self.path)} — {describe_settings(self.settings)}"


class OrbitalComparator(QWidget):
    def __init__(self, context):
        mw = context.get_main_window()
        super().__init__(mw)
        self.context = context
        self.mw = mw
        self.setWindowTitle("Orbital Comparator")
        self.setWindowFlags(Qt.WindowType.Window)
        self.setAcceptDrops(True)
        self.resize(820, 560)

        self.slots = []
        self._ready = False
        self._suspend = 0
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)
        hint = QLabel(
            "Drop .cube files here to overlay them on the current structure — "
            "several at once is fine. Colour, isovalue, opacity and style "
            "apply instantly."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet(
            "border: 1px dashed palette(mid); border-radius: 4px; padding: 8px;"
        )
        layout.addWidget(hint)

        for i in range(SLOT_COUNT):
            slot = CubeSlot(i, self)
            self.slots.append(slot)
            layout.addWidget(slot.box)

        btns = QHBoxLayout()
        # Secondary to dropping files on the window, which is the usual route.
        btn_load_many = QPushButton("Open Files…")
        btn_load_many.setToolTip(
            "Fill the empty slots from a multi-file selection. Dropping the "
            "files onto this window does the same thing."
        )
        btn_load_many.clicked.connect(self.load_many)
        btns.addWidget(btn_load_many)

        btn_sync = QPushButton("Sync Iso from Cube 1")
        btn_sync.setToolTip(
            "Give every slot Cube 1's contour level, so the shapes compare."
        )
        btn_sync.clicked.connect(self.sync_isovalue)
        btns.addWidget(btn_sync)

        btn_clear = QPushButton("Clear All")
        btn_clear.clicked.connect(self.clear_all)
        btns.addWidget(btn_clear)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btns.addWidget(btn_close)
        layout.addLayout(btns)

        self.lbl_status = QLabel("No cubes loaded.")
        layout.addWidget(self.lbl_status)

        # Wired last: everything above set widget values itself.
        for slot in self.slots:
            slot.check_on.toggled.connect(self.on_live_change)
            slot.spin_iso.valueChanged.connect(self.on_live_change)
            slot.spin_opacity.valueChanged.connect(self.on_live_change)
            slot.combo_style.currentTextChanged.connect(self.on_live_change)
            slot.check_smooth.toggled.connect(self.on_live_change)
        self._ready = True

    # -- loading -----------------------------------------------------------

    def load_into(self, slot, path=None):
        """Put one cube file into *slot*."""
        if path is None:
            path, _ = QFileDialog.getOpenFileName(
                self, "Open Cube File", "", "Cube Files (*.cube *.cub);;All Files (*)"
            )
        if not path:
            return False
        return self._read_into(slot, path)

    def _read_into(self, slot, path):
        if np is None or pv is None:
            QMessageBox.warning(
                self,
                "Missing dependency",
                "numpy and pyvista are required to read cube files.",
            )
            return False
        try:
            meta = parse_cube(path)
            slot.grid = build_grid(meta)
        except (OSError, ValueError, IndexError, KeyError) as e:
            logging.warning("Could not read cube %s: %s", path, e)
            QMessageBox.warning(
                self, "Unreadable cube", f"Could not read:\n{path}\n\n{e}"
            )
            slot.grid = None
            slot.path = None
            return False

        slot.path = path
        slot.settings = read_generation_settings(path)
        slot.lbl_info.setText(slot.describe())
        slot.box.setTitle(f"Cube {slot.index + 1}: {os.path.basename(path)}")
        self._suspend += 1
        try:
            slot.check_on.setChecked(True)
        finally:
            self._suspend -= 1
        self.render_all()
        return True

    def load_many(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Open Cube Files", "", "Cube Files (*.cube *.cub);;All Files (*)"
        )
        if paths:
            self.load_paths(paths)

    def load_paths(self, paths):
        """Fill empty slots first, then overwrite from the top."""
        empty = [s for s in self.slots if s.grid is None]
        targets = empty + [s for s in self.slots if s.grid is not None]
        loaded = 0
        self._suspend += 1
        try:
            for path, slot in zip(paths, targets):
                if self._read_into(slot, path):
                    loaded += 1
        finally:
            self._suspend -= 1
        self.render_all()
        if len(paths) > SLOT_COUNT:
            self.lbl_status.setText(
                f"Loaded {loaded}; only {SLOT_COUNT} cubes can be shown at once."
            )
        return loaded

    def clear_slot(self, slot):
        slot.grid = None
        slot.path = None
        slot.settings = {}
        slot.lbl_info.setText("no file")
        slot.box.setTitle(f"Cube {slot.index + 1}: (empty)")
        self._suspend += 1
        try:
            slot.check_on.setChecked(False)
        finally:
            self._suspend -= 1
        self._remove_actors(slot.prefix)
        self.render_all()

    def clear_all(self):
        self._suspend += 1
        try:
            for slot in self.slots:
                self.clear_slot(slot)
        finally:
            self._suspend -= 1
        self.render_all()

    # -- interaction -------------------------------------------------------

    def on_live_change(self, *_args):
        if not self._ready:
            return
        self.render_all()

    def pick_color(self, slot, which):
        from PyQt6.QtGui import QColor

        col = QColorDialog.getColor(QColor(slot.color(which)), self, "Select Color")
        if col.isValid():
            slot.set_color(which, col.name())
            self.render_all()

    def sync_isovalue(self):
        if not self.slots:
            return
        value = self.slots[0].spin_iso.value()
        self._suspend += 1
        try:
            for slot in self.slots[1:]:
                slot.spin_iso.setValue(value)
        finally:
            self._suspend -= 1
        self.render_all()

    # -- rendering ---------------------------------------------------------

    def render_all(self):
        if self._suspend or not self._ready:
            return
        plotter = getattr(self.mw, "plotter", None)
        if plotter is None:
            return

        shown = 0
        for slot in self.slots:
            self._remove_actors(slot.prefix)
            if not slot.is_on():
                continue
            try:
                iso = slot.spin_iso.value()
                pos = slot.grid.contour([iso], scalars="values")
                neg = slot.grid.contour([-iso], scalars="values")
                style = slot.combo_style.currentText().lower()
                for mesh, color, suffix in (
                    (pos, slot.color("p"), "_p"),
                    (neg, slot.color("n"), "_n"),
                ):
                    if mesh.n_points > 0:
                        plotter.add_mesh(
                            mesh,
                            color=color,
                            opacity=slot.spin_opacity.value(),
                            name=f"{slot.prefix}{suffix}",
                            style=style,
                            point_size=5,
                            smooth_shading=slot.check_smooth.isChecked(),
                        )
                shown += 1
            except (AttributeError, RuntimeError, ValueError) as e:
                logging.warning("Could not contour %s: %s", slot.path, e)

        try:
            plotter.render()
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)

        loaded = sum(1 for s in self.slots if s.grid is not None)
        self.lbl_status.setText(
            f"{loaded} cube(s) loaded, {shown} shown."
            if loaded
            else "No cubes loaded."
        )

    def _remove_actors(self, prefix):
        plotter = getattr(self.mw, "plotter", None)
        if plotter is None:
            return
        for suffix in ("_p", "_n"):
            try:
                plotter.remove_actor(f"{prefix}{suffix}")
            except (AttributeError, RuntimeError, KeyError) as _e:
                logging.warning("silenced: %s", _e)

    # -- drag and drop -----------------------------------------------------

    def dragEnterEvent(self, event):
        mime = event.mimeData()
        if mime and mime.hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        mime = event.mimeData()
        if not mime or not mime.hasUrls():
            return
        paths = [
            u.toLocalFile()
            for u in mime.urls()
            if u.toLocalFile().lower().endswith((".cube", ".cub"))
        ]
        if paths:
            self.load_paths(paths)
            event.acceptProposedAction()

    # -- teardown ----------------------------------------------------------

    def closeEvent(self, event):
        for slot in self.slots:
            self._remove_actors(slot.prefix)
        plotter = getattr(self.mw, "plotter", None)
        if plotter is not None:
            try:
                plotter.render()
            except (AttributeError, RuntimeError) as _e:
                logging.warning("silenced: %s", _e)
        event.accept()


def initialize(context):
    """V3 entry point. The menu action is added automatically by run()."""
    global PLUGIN_CONTEXT
    PLUGIN_CONTEXT = context

    def on_reset():
        mw = context.get_main_window()
        win = getattr(mw, "orbital_comparator_window", None)
        if win is not None:
            win.close()

    context.register_document_reset_handler(on_reset)

    # No register_drop_handler here on purpose: a drop on the MAIN window
    # belongs to the Cube File Viewers, and a second claimant would only make
    # which plugin answers depend on load order. Cubes reach this plugin by
    # being dropped onto its own window, which is unambiguous.

    def open_cube(file_path):
        """Last-resort opener: show the comparison window with this cube in it.

        Only reached when no Cube File Viewer is installed to take it first.
        """
        mw = context.get_main_window()
        win = getattr(mw, "orbital_comparator_window", None)
        if win is None:
            win = OrbitalComparator(context)
            mw.orbital_comparator_window = win
        win.show()
        win.raise_()
        win.load_paths([file_path])

    if hasattr(context, "register_file_opener"):
        for ext in (".cube", ".cub"):
            context.register_file_opener(
                ext, open_cube, priority=FILE_OPENER_PRIORITY
            )


def run(mw):
    if hasattr(mw, "host"):
        mw = mw.host
    context = PLUGIN_CONTEXT
    if not context:
        return

    win = getattr(mw, "orbital_comparator_window", None)
    if win is None:
        win = OrbitalComparator(context)
        mw.orbital_comparator_window = win

    if win.isVisible():
        win.close()
    else:
        win.show()
        win.raise_()
