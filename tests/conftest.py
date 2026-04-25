"""
Shared test infrastructure for moleditpy-plugins.

Provides:
- ``mock_optional_imports()`` — context manager that intercepts all heavy
  optional dependencies (PyQt6, rdkit, moleditpy, etc.) via a MetaPathFinder
  and replaces them with MagicMock modules so plugin tests run headlessly.
- ``load_plugin(path)`` — load a plugin .py file with deps already mocked.
- ``make_context()`` — build a stub PluginContext (MagicMock with a non-None
  main window).
- ``visible_py_plugins()`` — iterate visible single-file .py plugins from the
  registry, optionally filtered to those with a given entry-point function.
"""

from __future__ import annotations

import ast
import contextlib
import importlib.abc
import importlib.machinery
import importlib.util
import json
import sys
from pathlib import Path
from typing import Generator
from unittest.mock import MagicMock

ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "REGISTRY" / "plugins.json"
PLUGINS_DIR = ROOT / "plugins"

# ---------------------------------------------------------------------------
# Packages whose imports are replaced with MagicMock
# ---------------------------------------------------------------------------

BLOCKED_TOPS: frozenset[str] = frozenset(
    {
        "PyQt6",
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
# MetaPathFinder / Loader
# ---------------------------------------------------------------------------


class _MagicLoader(importlib.abc.Loader):
    def create_module(self, spec: importlib.machinery.ModuleSpec) -> MagicMock:
        m = MagicMock()
        m.__name__ = spec.name
        m.__spec__ = spec
        m.__path__ = []
        m.__package__ = spec.name.split(".")[0]
        return m  # type: ignore[return-value]

    def exec_module(self, module: object) -> None:
        pass


class _MagicFinder(importlib.abc.MetaPathFinder):
    _loader = _MagicLoader()

    def find_spec(
        self,
        fullname: str,
        path: object,
        target: object = None,
    ) -> importlib.machinery.ModuleSpec | None:
        if fullname.split(".")[0] in BLOCKED_TOPS:
            return importlib.machinery.ModuleSpec(fullname, self._loader)
        return None


@contextlib.contextmanager
def mock_optional_imports() -> Generator[None, None, None]:
    """
    Context manager that stubs all optional/heavy imports with MagicMock.

    Any module whose top-level package is in ``BLOCKED_TOPS`` is intercepted
    by a custom ``MetaPathFinder`` and replaced with a ``MagicMock`` for the
    duration of the block.  Previously-loaded real copies are restored on exit.
    """
    removed = {
        k: sys.modules.pop(k)
        for k in list(sys.modules)
        if k.split(".")[0] in BLOCKED_TOPS
    }
    finder = _MagicFinder()
    sys.meta_path.insert(0, finder)
    try:
        yield
    finally:
        sys.meta_path.remove(finder)
        sys.modules.update(removed)
        for k in list(sys.modules):
            if k.split(".")[0] in BLOCKED_TOPS and k not in removed:
                del sys.modules[k]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def load_plugin(path: Path) -> object:
    """
    Load *path* as an isolated module.  Must be called inside
    a ``mock_optional_imports()`` block.
    """
    mod_name = f"_smoke_{path.stem}_{abs(hash(path))}"
    spec = importlib.util.spec_from_file_location(mod_name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def make_context() -> MagicMock:
    """
    Return a stub PluginContext with a non-None main window.

    Plugins may assume ``get_main_window()`` is non-None at initialize() time
    (the host app is fully started before plugins load).
    """
    ctx = MagicMock()
    ctx.get_main_window.return_value = MagicMock()
    return ctx


def visible_py_plugins(
    entry_point: str | None = None,
) -> list[tuple[str, Path]]:
    """
    Return ``(name, path)`` for every visible single-file ``.py`` plugin in
    the registry.

    Args:
        entry_point: If given, only return plugins that define a top-level
            function with this name (e.g. ``"initialize"``, ``"run"``).
            Pass ``None`` to return all visible single-file plugins.
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
