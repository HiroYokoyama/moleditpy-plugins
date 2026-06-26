"""
Shared test infrastructure for headless GUI tests.

Key difference from tests/conftest.py: PyQt6 is NOT blocked, so real Qt widgets
are created.  Chemistry/scientific libraries (rdkit, numpy, pyvista, …) are still
replaced with MagicMock so no chemistry stack needs to be installed.

Environment requirement (set by CI, or locally before running):
    QT_QPA_PLATFORM=offscreen
"""

from __future__ import annotations

import contextlib
import importlib.abc
import importlib.machinery
import importlib.util
import sys
from pathlib import Path
from typing import Generator
from unittest.mock import MagicMock

import pytest

# On Windows, PyQt6's Qt DLLs must be discoverable before any PyQt6 sub-module
# is imported.  When PySide6 is also installed the wrong Qt runtime can win the
# DLL search; os.add_dll_directory() pins the correct one.
import os
import sys as _sys

if _sys.platform == "win32":
    try:
        import PyQt6 as _pyqt6_pkg  # noqa: PLC0415

        _qt6_bin = os.path.join(os.path.dirname(_pyqt6_pkg.__file__), "Qt6", "bin")
        if os.path.isdir(_qt6_bin):
            os.add_dll_directory(_qt6_bin)
    except Exception:
        pass  # best-effort; let the actual import surface the real error


ROOT = Path(__file__).resolve().parents[1]
PLUGINS_DIR = ROOT / "plugins"

# ---------------------------------------------------------------------------
# Same as BLOCKED_TOPS in tests/conftest.py, but PyQt6 is intentionally absent
# so plugin classes that inherit from Qt bases become real widgets, not mocks.
# ---------------------------------------------------------------------------

BLOCKED_CHEMISTRY: frozenset[str] = frozenset(
    {
        "rdkit",
        "pyvista",
        "pyvistaqt",
        "numpy",
        "scipy",
        "PIL",
        "google",
        "openai",
        "anthropic",
        "truststore",
        "pymatgen",
        "openbabel",
        "pyscf",
        "ase",
        "sip",
        "modules",
        "moleditpy_linux",
        "pybel",
        "matplotlib",
        "markdown",
        "vtk",
        "vtkmodules",
        "cryptography",
        "requests",
        "moleditpy",
    }
)


# ---------------------------------------------------------------------------
# MetaPathFinder / Loader — same mechanism as tests/conftest.py
# ---------------------------------------------------------------------------


class _ChemLoader(importlib.abc.Loader):
    def create_module(self, spec: importlib.machinery.ModuleSpec) -> MagicMock:
        m = MagicMock()
        m.__name__ = spec.name
        m.__spec__ = spec
        m.__path__ = []
        m.__package__ = spec.name.split(".")[0]
        return m  # type: ignore[return-value]

    def exec_module(self, module: object) -> None:
        pass


class _ChemFinder(importlib.abc.MetaPathFinder):
    _loader = _ChemLoader()

    def find_spec(
        self,
        fullname: str,
        path: object,
        target: object = None,
    ) -> importlib.machinery.ModuleSpec | None:
        if fullname.split(".")[0] in BLOCKED_CHEMISTRY:
            return importlib.machinery.ModuleSpec(fullname, self._loader)
        return None


@contextlib.contextmanager
def mock_chemistry_imports() -> Generator[None, None, None]:
    """
    Block chemistry/scientific imports (replacing them with MagicMock) while
    letting PyQt6 through so real Qt widgets are constructed.
    """
    removed = {
        k: sys.modules.pop(k)
        for k in list(sys.modules)
        if k.split(".")[0] in BLOCKED_CHEMISTRY
    }
    finder = _ChemFinder()
    sys.meta_path.insert(0, finder)
    try:
        yield
    finally:
        sys.meta_path.remove(finder)
        sys.modules.update(removed)
        for k in list(sys.modules):
            if k.split(".")[0] in BLOCKED_CHEMISTRY and k not in removed:
                del sys.modules[k]


def load_plugin_for_gui(path: Path) -> object:
    """
    Load *path* as an isolated module with chemistry mocked.
    Must be called inside a ``mock_chemistry_imports()`` block.
    """
    mod_name = f"_gui_{path.stem}_{abs(hash(path))}"
    spec = importlib.util.spec_from_file_location(mod_name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


# ---------------------------------------------------------------------------
# Session-scoped QApplication
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def qapp():
    """
    One QApplication for the entire test session.

    pytest-qt provides its own qapp fixture, but we avoid that dependency here.
    QT_QPA_PLATFORM=offscreen must be set before this fixture is used (done by
    CI via the workflow env block, or locally by the developer).
    """
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(sys.argv[:1])
    yield app


# ---------------------------------------------------------------------------
# Registry helper (mirrors tests/conftest.py::visible_py_plugins)
# ---------------------------------------------------------------------------

import ast
import json

REGISTRY_PATH = ROOT / "REGISTRY" / "plugins.json"


def visible_py_plugins(entry_point: str | None = None) -> list[tuple[str, Path]]:
    """
    Return (name, path) for every visible single-file .py plugin in the registry.
    Mirrors tests/conftest.py so tests_gui/ tests don't need to import from tests/.
    """
    reg = json.loads(REGISTRY_PATH.read_text(encoding="utf-8-sig"))
    result = []
    for entry in reg:
        if not entry.get("visible", False):
            continue
        url = entry.get("downloadUrl", "")
        if "/plugins/" not in url:
            continue
        rel = url.split("/plugins/")[-1]
        path = PLUGINS_DIR / rel
        if path.suffix != ".py" or not path.exists():
            continue
        if entry_point is not None:
            source = path.read_text(encoding="utf-8", errors="ignore")
            try:
                tree = ast.parse(source, filename=str(path))
            except SyntaxError:
                continue
            defined = {
                node.name
                for node in ast.walk(tree)
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            }
            if entry_point not in defined:
                continue
        result.append((entry["name"], path))
    return result


# ---------------------------------------------------------------------------
# Main-app integration helpers
# ---------------------------------------------------------------------------

_APP_DIR = ROOT.parent / "python_molecular_editor"
_APP_SRC = _APP_DIR / "moleditpy" / "src"


def _ensure_app() -> Path:
    """
    Return the main-app src path, cloning the repo first if it is absent.

    The clone is --depth 1 (read-only history) placed at the sibling path
    ../python_molecular_editor, mirroring the test-api CI job layout.
    """
    if not _APP_SRC.exists():
        import subprocess

        subprocess.run(
            [
                "git",
                "clone",
                "--depth",
                "1",
                "https://github.com/HiroYokoyama/python_molecular_editor.git",
                str(_APP_DIR),
            ],
            check=True,
            text=True,
        )
    return _APP_SRC


@pytest.fixture(scope="session")
def app_src_path() -> Path:
    """Session fixture: path to the main-app Python src root (auto-clones if absent)."""
    return _ensure_app()


@pytest.fixture(scope="session")
def real_plugin_context_class(app_src_path: Path):
    """
    Return the real PluginContext class imported from the main app.

    We load plugin_interface.py *directly by file path* rather than via
    ``from moleditpy.plugins.plugin_interface import PluginContext``.
    The latter triggers moleditpy/__init__.py, which imports rdkit — a
    library not present in the CI environment.  Loading by spec bypasses
    the package __init__ entirely; plugin_interface.py itself only uses
    stdlib (logging, typing) so this is always safe.
    """
    pi_path = app_src_path / "moleditpy" / "plugins" / "plugin_interface.py"
    spec = importlib.util.spec_from_file_location("_moleditpy_plugin_interface", pi_path)
    assert spec is not None and spec.loader is not None, f"Cannot find {pi_path}"
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod.PluginContext


_STUB_DEFAULT = object()  # sentinel — distinguishes "not given" from explicit None


def make_stub_manager(*, main_window: Any = _STUB_DEFAULT) -> MagicMock:
    """
    Minimal stub that satisfies PluginContext's delegation calls.

    By default, get_main_window() returns a MagicMock (matching make_context() in
    tests/conftest.py) so plugins that access ``mw.something`` in initialize()
    get a MagicMock chain instead of crashing on NoneType.

    Pass ``main_window=None`` explicitly when you need to exercise the
    "no main window" fallback paths in PluginContext property tests.
    """
    mgr = MagicMock()
    mgr.get_main_window.return_value = (
        MagicMock() if main_window is _STUB_DEFAULT else main_window
    )
    return mgr
