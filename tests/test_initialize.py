"""
Smoke tests: call initialize(context) on every visible plugin that exposes it.

Verifies that initialize() completes without raising an exception when given
a realistic stub context (non-None main window, all heavy deps mocked out).

Rationale: AST-based tests (test_imports.py) catch syntax errors and missing
stdlib imports, but cannot catch runtime errors that only surface when
initialize() actually executes — e.g. missing attribute access, bad argument
count, or unguarded None dereferences.

All heavy optional dependencies (PyQt6, rdkit, pyvista, etc.) are intercepted
by a custom MetaPathFinder so these tests run without any GUI or chemistry
libraries installed.
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
from unittest.mock import MagicMock

import pytest

ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "REGISTRY" / "plugins.json"
PLUGINS_DIR = ROOT / "plugins"

# ---------------------------------------------------------------------------
# Optional dependency list — any import whose top-level package is in this
# set is intercepted and replaced with a MagicMock module.
# ---------------------------------------------------------------------------

_BLOCKED_TOPS: frozenset[str] = frozenset(
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
# MetaPathFinder that returns MagicMock modules for every blocked package
# ---------------------------------------------------------------------------


class _MagicLoader(importlib.abc.Loader):
    """Loader that creates a MagicMock as the module object."""

    def create_module(self, spec: importlib.machinery.ModuleSpec) -> MagicMock:
        m = MagicMock()
        m.__name__ = spec.name
        m.__spec__ = spec
        m.__path__ = []  # mark as package so sub-imports work
        m.__package__ = spec.name.split(".")[0]
        return m  # type: ignore[return-value]

    def exec_module(self, module: object) -> None:  # noqa: D102
        pass  # nothing to execute — the MagicMock is already fully functional


class _MagicFinder(importlib.abc.MetaPathFinder):
    """Intercept imports of blocked top-level packages and return mocks."""

    _loader = _MagicLoader()

    def find_spec(
        self,
        fullname: str,
        path: object,
        target: object = None,
    ) -> importlib.machinery.ModuleSpec | None:
        if fullname.split(".")[0] in _BLOCKED_TOPS:
            return importlib.machinery.ModuleSpec(fullname, self._loader)
        return None


@contextlib.contextmanager
def _mock_optional_imports():  # type: ignore[return]
    """Context manager: install magic finder and clean blocked caches on exit."""
    # Temporarily remove already-loaded blocked modules (e.g. the real moleditpy)
    removed = {
        k: sys.modules.pop(k)
        for k in list(sys.modules)
        if k.split(".")[0] in _BLOCKED_TOPS
    }
    finder = _MagicFinder()
    sys.meta_path.insert(0, finder)
    try:
        yield
    finally:
        sys.meta_path.remove(finder)
        # Restore real modules that were present before the context
        sys.modules.update(removed)
        # Drop any mock modules that were loaded during the test
        for k in list(sys.modules):
            if k.split(".")[0] in _BLOCKED_TOPS and k not in removed:
                del sys.modules[k]


# ---------------------------------------------------------------------------
# Parametrize helpers — mirrors the pattern in test_imports.py
# ---------------------------------------------------------------------------


def _visible_plugins_with_initialize() -> list[tuple[str, Path]]:
    """
    Return (name, path) for every visible single-file .py plugin in the
    registry that defines an ``initialize()`` function.
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
        source = path.read_text(encoding="utf-8", errors="ignore")
        try:
            tree = ast.parse(source, filename=str(path))
        except SyntaxError:
            continue
        has_initialize = any(
            isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name == "initialize"
            for node in ast.walk(tree)
        )
        if has_initialize:
            result.append((entry["name"], path))
    return result


_PLUGINS = _visible_plugins_with_initialize()


def _load_plugin(path: Path) -> object:
    """
    Load *path* as an isolated module with all optional deps mocked.

    Must be called inside a ``_mock_optional_imports()`` context.
    """
    # Use a unique name per path to avoid module cache collisions
    mod_name = f"_smoke_{path.stem}_{abs(hash(path))}"
    spec = importlib.util.spec_from_file_location(mod_name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def _make_context() -> MagicMock:
    """
    Return a stub PluginContext with a non-None main window.

    Plugins are allowed to assume main_window is not None at initialize()
    time (the app is fully started when plugins are loaded), so we return
    a MagicMock rather than None.
    """
    ctx = MagicMock()
    ctx.get_main_window.return_value = MagicMock()
    return ctx


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,path", _PLUGINS, ids=[n for n, _ in _PLUGINS]
)
def test_initialize_does_not_raise(name: str, path: Path) -> None:
    """``initialize(context)`` completes without raising any exception."""
    with _mock_optional_imports():
        mod = _load_plugin(path)
        ctx = _make_context()
        mod.initialize(ctx)  # type: ignore[attr-defined]
