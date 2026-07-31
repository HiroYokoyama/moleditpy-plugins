"""
Stateful headless-Qt stubs for the Orbital Comparator tests.

The repo-wide conftest mocks Qt with MagicMocks, which is fine for pure
helpers but cannot support this plugin's tests: QWidget must be subclassable
for OrbitalComparator to be a real type, and the dialog reads back what it
wrote (``spin.setValue(0.05)`` then ``spin.value()``), which a MagicMock
forgets. So this installs a small set of widgets that keep real state, loads
one module against them, and restores ``sys.modules`` exactly as it was — no
crippled stub is allowed to outlive the load and corrupt other test modules.

Signals record their connections but do NOT fire on programmatic setters:
widgets are wired mid-``init_ui`` and auto-emitting would run slots against a
half-built window. Tests call ``signal.emit(...)`` when they want the wiring
exercised, which is also how the real binding behaves for user interaction.
"""

import importlib.util
import sys
import types
from unittest.mock import MagicMock

_SENTINEL = object()


class _Signal:
    def __init__(self):
        self._slots = []

    def connect(self, slot):
        self._slots.append(slot)

    def disconnect(self, slot=None):
        if slot is None:
            self._slots.clear()
        elif slot in self._slots:
            self._slots.remove(slot)

    def emit(self, *a):
        for slot in list(self._slots):
            slot(*a)


class _PermissiveMeta(type):
    def __getattr__(cls, name):
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        return MagicMock()


class _Base(metaclass=_PermissiveMeta):
    """Unknown attributes stay permissive; the round-tripped ones are real."""

    def __getattr__(self, attr):
        if attr.startswith("_"):
            raise AttributeError(attr)
        return MagicMock()

    def setStyleSheet(self, sheet):
        self.__dict__["_stylesheet"] = "" if sheet is None else str(sheet)

    def styleSheet(self):
        return self.__dict__.get("_stylesheet", "")

    def setToolTip(self, text):
        self.__dict__["_tooltip"] = "" if text is None else str(text)

    def toolTip(self):
        return self.__dict__.get("_tooltip", "")

    def setEnabled(self, value):
        self.__dict__["_enabled"] = bool(value)

    def isEnabled(self):
        return self.__dict__.get("_enabled", True)

    def setVisible(self, value):
        self.__dict__["_visible"] = bool(value)

    def isVisible(self):
        return self.__dict__.get("_visible", False)

    def show(self):
        self.setVisible(True)

    def hide(self):
        self.setVisible(False)

    def close(self):
        self.setVisible(False)
        return True


class _Widget(_Base):
    def __init__(self, *a, **k):
        pass


class _Label(_Base):
    def __init__(self, text="", *a, **k):
        self._text = text if isinstance(text, str) else ""

    def setText(self, t):
        self._text = "" if t is None else str(t)

    def text(self):
        return self._text


class _GroupBox(_Base):
    def __init__(self, title="", *a, **k):
        self._title = title if isinstance(title, str) else ""

    def setTitle(self, t):
        self._title = "" if t is None else str(t)

    def title(self):
        return self._title


class _Button(_Base):
    def __init__(self, *a, **k):
        self._text = a[0] if a and isinstance(a[0], str) else ""
        self.clicked = _Signal()

    def setText(self, t):
        self._text = "" if t is None else str(t)

    def text(self):
        return self._text


class _CheckBox(_Base):
    def __init__(self, *a, **k):
        self._checked = False
        self.toggled = _Signal()
        self.stateChanged = _Signal()

    def setChecked(self, v):
        self._checked = bool(v)

    def isChecked(self):
        return self._checked


class _SpinBox(_Base):
    def __init__(self, *a, **k):
        self._value = 0.0
        self._min, self._max = -1e18, 1e18
        self.valueChanged = _Signal()

    def setRange(self, lo, hi):
        self._min, self._max = lo, hi
        self._value = min(max(self._value, lo), hi)

    def setValue(self, v):
        self._value = min(max(v, self._min), self._max)

    def value(self):
        return self._value


class _ComboBox(_Base):
    def __init__(self, *a, **k):
        self._items = []
        self._idx = -1
        self.currentTextChanged = _Signal()
        self.currentIndexChanged = _Signal()

    def addItem(self, text, data=None):
        self._items.append(text)
        if self._idx < 0:
            self._idx = 0

    def addItems(self, texts):
        for t in texts:
            self.addItem(t)

    def count(self):
        return len(self._items)

    def itemText(self, i):
        return self._items[i] if 0 <= i < len(self._items) else ""

    def setCurrentIndex(self, i):
        self._idx = i

    def currentIndex(self):
        return self._idx

    def setCurrentText(self, text):
        if text in self._items:
            self._idx = self._items.index(text)

    def currentText(self):
        return self.itemText(self._idx)

    def findText(self, text):
        return self._items.index(text) if text in self._items else -1


class _LayoutItem:
    def __init__(self, widget=None, layout=None):
        self._w, self._l = widget, layout

    def widget(self):
        return self._w

    def layout(self):
        return self._l


class _Layout(_Base):
    def __init__(self, *a, **k):
        self._items = []

    def addWidget(self, w, *a, **k):
        self._items.append(_LayoutItem(widget=w))

    def addLayout(self, lay, *a, **k):
        self._items.append(_LayoutItem(layout=lay))

    def addStretch(self, *a, **k):
        self._items.append(_LayoutItem())

    def count(self):
        return len(self._items)

    def itemAt(self, i):
        return self._items[i] if 0 <= i < len(self._items) else None

    def takeAt(self, i):
        return self._items.pop(i) if 0 <= i < len(self._items) else None


class QColorStub:
    _NAMED = {"red": (255, 0, 0), "blue": (0, 0, 255), "black": (0, 0, 0)}

    def __init__(self, *a):
        r = g = b = 0
        if len(a) == 1 and isinstance(a[0], str):
            s = a[0].strip().lower()
            if s.startswith("#") and len(s) == 7:
                r, g, b = int(s[1:3], 16), int(s[3:5], 16), int(s[5:7], 16)
            else:
                r, g, b = self._NAMED.get(s, (0, 0, 0))
        elif len(a) >= 3:
            r, g, b = int(a[0]), int(a[1]), int(a[2])
        self._r, self._g, self._b = r, g, b

    def red(self):
        return self._r

    def green(self):
        return self._g

    def blue(self):
        return self._b

    def name(self):
        return f"#{self._r:02x}{self._g:02x}{self._b:02x}"

    def isValid(self):
        return True


_WIDGETS = {
    "QWidget": _Widget,
    "QDialog": _Widget,
    "QMainWindow": _Widget,
    "QVBoxLayout": _Layout,
    "QHBoxLayout": _Layout,
    "QGridLayout": _Layout,
    "QFormLayout": _Layout,
    "QLabel": _Label,
    "QGroupBox": _GroupBox,
    "QPushButton": _Button,
    "QToolButton": _Button,
    "QCheckBox": _CheckBox,
    "QDoubleSpinBox": _SpinBox,
    "QSpinBox": _SpinBox,
    "QComboBox": _ComboBox,
}

_SWAP_KEYS = (
    "PyQt6",
    "PyQt6.QtWidgets",
    "PyQt6.QtCore",
    "PyQt6.QtGui",
    "pyvista",
)


class _StubMesh:
    """A contour result: just enough for the "did it draw?" checks."""

    def __init__(self, n_points=10):
        self.n_points = n_points


class _StubStructuredGrid:
    """Stands in for pyvista's grid.

    CI installs numpy but not pyvista, so without this the plugin's dependency
    guard (``np is None or pv is None``) refuses every load and the whole
    loading suite fails there while passing locally. Building the grid is the
    only thing pyvista is needed for; the point/dimension bookkeeping this
    keeps is what the tests actually assert on.
    """

    def __init__(self):
        self.points = None
        self.dimensions = None
        self.point_data = {}
        self.contour_levels = []

    def contour(self, levels, scalars=None):
        self.contour_levels.append(levels[0])
        return _StubMesh()


def _pyvista_stub():
    mod = types.ModuleType("pyvista")
    mod.StructuredGrid = _StubStructuredGrid
    mod.PolyData = MagicMock()
    return mod


def _module(name, extra):
    mod = MagicMock()
    mod.__name__ = name
    for k, v in extra.items():
        setattr(mod, k, v)
    return mod


def load_with_stateful_qt(path, modname="orbital_comparator_isolated"):
    """Load *path* against the stateful widgets, then restore sys.modules."""
    saved = {k: sys.modules.get(k, _SENTINEL) for k in _SWAP_KEYS}

    pkg = types.ModuleType("PyQt6")
    pkg.__path__ = []
    widgets = _module("PyQt6.QtWidgets", _WIDGETS)
    core = _module("PyQt6.QtCore", {})
    gui = _module("PyQt6.QtGui", {"QColor": QColorStub})
    pkg.QtWidgets, pkg.QtCore, pkg.QtGui = widgets, core, gui
    sys.modules.update(
        {
            "PyQt6": pkg,
            "PyQt6.QtWidgets": widgets,
            "PyQt6.QtCore": core,
            "PyQt6.QtGui": gui,
            # Real pyvista would drag in VTK; the plugin only needs a grid
            # object, and CI has no pyvista at all.
            "pyvista": _pyvista_stub(),
        }
    )
    try:
        spec = importlib.util.spec_from_file_location(modname, str(path))
        mod = importlib.util.module_from_spec(spec)
        sys.modules[modname] = mod
        spec.loader.exec_module(mod)
        return mod
    finally:
        for k in _SWAP_KEYS:
            if saved[k] is _SENTINEL:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = saved[k]
