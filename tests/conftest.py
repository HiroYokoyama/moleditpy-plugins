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
import textwrap
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


# ---------------------------------------------------------------------------
# Real-numpy passthrough (for export-script generators that do real vector
# math) and shared fake-rdkit structures
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def mocks_with_real_numpy() -> Generator[None, None, None]:
    """
    Like mock_optional_imports(), but numpy resolves to the real package.

    The export generators do ``import numpy as np`` at call time; placing the
    real modules in sys.modules bypasses the mocking MetaPathFinder for numpy
    only, so vector math is real while rdkit/PyQt6 stay mocked.

    ALL already-imported ``numpy.*`` submodules must be restored, not just the
    top-level package: numpy resolves internals like ``numpy._core._methods``
    through sys.modules at call time, and if the MetaPathFinder answers that
    import with a MagicMock the real numpy package is corrupted for the rest
    of the process (e.g. ``np.allclose(0, 1)`` starts returning True).
    """
    real_mods = {
        k: v
        for k, v in sys.modules.items()
        if k == "numpy" or k.startswith("numpy.")
    }
    with mock_optional_imports():
        sys.modules.update(real_mods)
        try:
            yield
        finally:
            for k in real_mods:
                sys.modules.pop(k, None)


class P3:
    """3-component point supporting .x/.y/.z and the sequence protocol."""

    def __init__(self, x, y, z):
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __len__(self):
        return 3

    def __getitem__(self, i):
        return (self.x, self.y, self.z)[i]

    def __iter__(self):
        return iter((self.x, self.y, self.z))


class FakeAtom:
    def __init__(self, idx, symbol):
        self.idx = idx
        self.symbol = symbol
        self.bonds = []

    def GetIdx(self):
        return self.idx

    def GetSymbol(self):
        return self.symbol

    def GetDegree(self):
        return len(self.bonds)

    def GetBonds(self):
        return list(self.bonds)

    def GetNeighbors(self):
        out = []
        for b in self.bonds:
            out.append(b.end_atom if b.begin_atom is self else b.begin_atom)
        return out


class FakeBond:
    def __init__(self, begin_atom, end_atom, bond_type):
        self.begin_atom = begin_atom
        self.end_atom = end_atom
        self.bond_type = bond_type
        begin_atom.bonds.append(self)
        end_atom.bonds.append(self)

    def GetBeginAtomIdx(self):
        return self.begin_atom.idx

    def GetEndAtomIdx(self):
        return self.end_atom.idx

    def GetBondType(self):
        return self.bond_type


class FakeConf:
    def __init__(self, coords):
        self.coords = {i: P3(*c) for i, c in enumerate(coords)}

    def GetAtomPosition(self, idx):
        return self.coords[idx]


class FakeMol:
    def __init__(self, symbols, coords, bonds=()):
        # bonds: iterable of (begin_idx, end_idx, bond_type)
        self.atoms = [FakeAtom(i, s) for i, s in enumerate(symbols)]
        self.conf = FakeConf(coords)
        self.bonds = [FakeBond(self.atoms[b], self.atoms[e], t) for b, e, t in bonds]

    def GetNumAtoms(self):
        return len(self.atoms)

    def GetNumBonds(self):
        return len(self.bonds)

    def GetNumConformers(self):
        return 1

    def GetConformer(self):
        return self.conf

    def GetAtoms(self):
        return list(self.atoms)

    def GetBonds(self):
        return list(self.bonds)

    def GetAtomWithIdx(self, idx):
        return self.atoms[idx]


def extract_function(
    path: Path,
    class_name: str | None,
    fn_name: str,
    extra_globals: dict | None = None,
):
    """
    Extract a function or method from source via AST and exec it standalone.

    Needed for methods of Qt-derived classes: the mocked Qt bases make the
    class object a MagicMock, so methods can't be reached via the class.

    ``class_name=None`` extracts a module-level function instead of a method.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    scope = tree
    if class_name is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == class_name:
                scope = node
                break
        else:
            raise AssertionError(f"class {class_name} not found in {path}")
    for node in scope.body:
        if isinstance(node, ast.FunctionDef) and node.name == fn_name:
            segment = ast.get_source_segment(source, node)
            ns = {"MagicMock": MagicMock}
            ns.update(extra_globals or {})
            exec(textwrap.dedent(segment), ns)  # noqa: S102
            return ns[fn_name]
    raise AssertionError(f"function {fn_name} not found in {path}")
