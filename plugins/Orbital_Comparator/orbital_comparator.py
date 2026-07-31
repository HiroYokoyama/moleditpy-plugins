"""Orbital Comparator — side-by-side comparison of cube files.

The molecular counterpart of this plugin compares structures; this one
compares volumetric data. Load any .cube files (molecular orbitals, densities,
ESP maps) and show up to four of them at once over the current structure, each
with its own lobe colours, contour level, opacity and style.

Each slot reports its grid shape, atom count and units, read from the cube
format itself — the comment lines are free text and belong to whichever
program wrote them, so nothing there is interpreted.

Cubes also carry the geometry they were computed on. A spin box at the top
picks which loaded cube supplies the structure shown underneath the
isosurfaces, which matters when the files come from different jobs.
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
    QMessageBox,
    QPushButton,
    QSpinBox,
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


#: Enough of the periodic table to name the atoms a cube carries.
ELEMENTS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te "
    "I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir "
    "Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No "
    "Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


def element_symbol(z):
    """Symbol for atomic number *z*, or "X" if it is out of range."""
    try:
        z = int(z)
    except (TypeError, ValueError):
        return "X"
    return ELEMENTS[z] if 0 <= z < len(ELEMENTS) else "X"


def describe_cube(meta):
    """One line of generic cube facts: grid shape and atom count.

    Deliberately reads only what the cube format itself defines. The comment
    lines are free text — the ORCA Result Analyzer happens to record its
    version and grid settings there, but that is one program's private
    convention and not something this plugin should claim to understand.
    """
    nx, ny, nz = meta["dims"]
    n_atoms = len(meta.get("atoms") or ())
    units = "Å" if meta.get("is_angstrom") else "Bohr"
    return f"{nx}×{ny}×{nz} grid, {n_atoms} atoms, {units}"


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

    atoms = []
    for line in lines[6 : 6 + n_atoms]:
        fields = line.split()
        if len(fields) < 5:
            continue
        atoms.append(
            (
                element_symbol(float(fields[0])),
                (float(fields[2]), float(fields[3]), float(fields[4])),
            )
        )

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
        "atoms": atoms,
        "is_angstrom": nx < 0,
    }


def atoms_to_xyz(atoms, is_angstrom, source_name="cube"):
    """An XYZ block from a cube's atom records.

    Cube coordinates are Bohr unless the file flags Angstrom, and the host
    expects Angstrom — skipping the conversion shrinks the molecule to 53% and
    it no longer lines up with its own isosurface.
    """
    scale = 1.0 if is_angstrom else BOHR_TO_ANG
    lines = [str(len(atoms)), source_name]
    for symbol, (x, y, z) in atoms:
        lines.append(
            f"{symbol} {x * scale:.6f} {y * scale:.6f} {z * scale:.6f}"
        )
    return "\n".join(lines)


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
        self.meta = None

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
        if not self.meta:
            return os.path.basename(self.path)
        return f"{os.path.basename(self.path)} — {describe_cube(self.meta)}"


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
        self._loading_structure = False
        self._structure_loaded = False
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

        # Structure picker, above the slots: cubes from different jobs carry
        # different geometries, so which one the lobes sit on is a choice.
        struct_row = QHBoxLayout()
        struct_row.addWidget(QLabel("Structure from:"))
        self.spin_structure = QSpinBox()
        self.spin_structure.setRange(1, SLOT_COUNT)
        self.spin_structure.setPrefix("Cube ")
        self.spin_structure.setToolTip(
            "Which loaded cube's atoms to show under the isosurfaces."
        )
        struct_row.addWidget(self.spin_structure)

        self.btn_load_structure = QPushButton("Load Structure")
        self.btn_load_structure.setToolTip(
            "Replace the structure in the 3D view with this cube's atoms."
        )
        self.btn_load_structure.clicked.connect(self.load_structure)
        struct_row.addWidget(self.btn_load_structure)

        self.lbl_structure = QLabel("")
        struct_row.addWidget(self.lbl_structure, 1)
        layout.addLayout(struct_row)

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
        self.spin_structure.valueChanged.connect(
            lambda *_: self.refresh_structure_label()
        )
        self._ready = True
        self.refresh_structure_label()

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
        slot.meta = meta
        slot.lbl_info.setText(slot.describe())
        slot.box.setTitle(f"Cube {slot.index + 1}: {os.path.basename(path)}")
        self._suspend += 1
        try:
            slot.check_on.setChecked(True)
        finally:
            self._suspend -= 1
        self.render_all()
        self.refresh_structure_label()

        # The first cube to arrive brings the geometry with it; showing its
        # lobes over whatever was already in the viewer would be misleading.
        if not self._structure_loaded and meta.get("atoms"):
            self._suspend += 1
            try:
                self.spin_structure.setValue(slot.index + 1)
            finally:
                self._suspend -= 1
            self.refresh_structure_label()
            self.load_structure()
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
        slot.meta = None
        slot.lbl_info.setText("no file")
        slot.box.setTitle(f"Cube {slot.index + 1}: (empty)")
        self._suspend += 1
        try:
            slot.check_on.setChecked(False)
        finally:
            self._suspend -= 1
        self._remove_actors(slot.prefix)
        self.render_all()
        self.refresh_structure_label()

    def clear_all(self):
        self._suspend += 1
        try:
            for slot in self.slots:
                self.clear_slot(slot)
        finally:
            self._suspend -= 1
        self._structure_loaded = False
        self.render_all()

    # -- interaction -------------------------------------------------------

    def structure_slot(self):
        """The slot the spin box points at, 1-based in the UI."""
        index = self.spin_structure.value() - 1
        if 0 <= index < len(self.slots):
            return self.slots[index]
        return None

    def refresh_structure_label(self):
        """Say what the picked slot would load, or why it cannot."""
        slot = self.structure_slot()
        if slot is None or slot.meta is None:
            self.lbl_structure.setText("(no cube loaded in that slot)")
            self.btn_load_structure.setEnabled(False)
            return
        atoms = slot.meta.get("atoms") or []
        if not atoms:
            self.lbl_structure.setText("(that cube carries no atoms)")
            self.btn_load_structure.setEnabled(False)
            return
        self.lbl_structure.setText(f"{len(atoms)} atoms")
        self.btn_load_structure.setEnabled(True)

    def load_structure(self):
        """Show the picked cube's atoms in the host's 3D view."""
        slot = self.structure_slot()
        if slot is None or slot.meta is None:
            return
        atoms = slot.meta.get("atoms") or []
        if not atoms:
            return

        name = os.path.basename(slot.path or "cube")
        xyz = atoms_to_xyz(atoms, slot.meta.get("is_angstrom", False), source_name=name)
        try:
            # Held across the call: show_xyz_data clears the document, which
            # fires the reset handlers, and ours must not close this window
            # for a reset it caused itself.
            self._loading_structure = True
            self.context.show_xyz_data(xyz, source_name=name)
        except (AttributeError, RuntimeError, ValueError) as e:
            logging.warning("Could not load the cube structure: %s", e)
            QMessageBox.warning(
                self, "Could not load structure", f"The host refused it:\n{e}"
            )
            return
        finally:
            self._loading_structure = False

        self._structure_loaded = True
        self.enter_3d_viewer_mode()
        # The host redraws the scene, which drops the isosurfaces with it.
        self.render_all()
        self.reset_camera()

    def reset_camera(self):
        """Frame the new structure.

        The camera still points at whatever was on screen before, which for a
        differently sized or positioned molecule can leave the view empty.
        """
        plotter = getattr(self.mw, "plotter", None)
        if plotter is None:
            return
        try:
            plotter.reset_camera()
            plotter.render()
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)

    def enter_3d_viewer_mode(self):
        """Collapse the 2D panel and disable the editing tools.

        Same handshake the Cube File Viewer uses: prefer the context call,
        fall back to the ui_manager method on hosts that predate it. Purely
        cosmetic, so a host offering neither is not an error.
        """
        enter = getattr(self.context, "enter_3d_viewer_mode", None)
        if enter is not None:
            try:
                enter()
                return
            except (AttributeError, RuntimeError, TypeError) as _e:
                logging.warning("silenced: %s", _e)

        ui = getattr(self.mw, "ui_manager", None)
        legacy = getattr(ui, "enter_3d_viewer_ui_mode", None) if ui else None
        if legacy is not None:
            try:
                legacy()
            except (AttributeError, RuntimeError, TypeError) as _e:
                logging.warning("silenced: %s", _e)

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
        if win is None:
            return
        # Loading a structure goes through show_xyz_data, which clears the
        # document and fires these handlers -- so pressing Load Structure used
        # to close this window instantly. A reset we caused ourselves is not a
        # new document.
        # `is True` rather than truthiness: the flag is only ever a bool, and
        # anything else answering this attribute must not disarm a real reset.
        if getattr(win, "_loading_structure", False) is True:
            return
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
        win.enter_3d_viewer_mode()
