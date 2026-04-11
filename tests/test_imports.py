"""
Tests that all visible plugin .py files are importable (no syntax errors,
no broken top-level imports that would prevent the plugin from loading).

Heavy optional dependencies (PyQt6, pyvista, rdkit) are skipped gracefully
using importlib — we only fail on SyntaxError or missing stdlib imports.
"""

import ast
import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "REGISTRY" / "plugins.json"
PLUGINS_DIR = ROOT / "plugins"


def _visible_py_plugins():
    """Return (name, path) for every visible single-file .py plugin."""
    reg = json.loads(REGISTRY_PATH.read_text(encoding="utf-8-sig"))
    result = []
    for p in reg:
        if not p.get("visible", False):
            continue
        url = p.get("downloadUrl", "")
        if "/plugins/" not in url:
            continue
        rel = url.split("/plugins/")[-1]
        path = PLUGINS_DIR / rel
        if path.suffix == ".py" and path.exists():
            result.append((p["name"], path))
    return result


@pytest.mark.parametrize("name,path", _visible_py_plugins(), ids=[n for n, _ in _visible_py_plugins()])
def test_plugin_syntax(name, path):
    """Plugin file parses without SyntaxError."""
    source = path.read_text(encoding="utf-8", errors="ignore")
    try:
        ast.parse(source, filename=str(path))
    except SyntaxError as e:
        pytest.fail(f"{name}: SyntaxError at line {e.lineno}: {e.msg}")


@pytest.mark.parametrize("name,path", _visible_py_plugins(), ids=[n for n, _ in _visible_py_plugins()])
def test_plugin_has_run_entrypoint(name, path):
    """Plugin file defines at least one callable entry point (run/autorun/initialize)."""
    source = path.read_text(encoding="utf-8", errors="ignore")
    tree = ast.parse(source, filename=str(path))
    entrypoints = {"run", "autorun", "initialize"}
    defined = {
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    found = entrypoints & defined
    assert found, f"{name}: no entry point found (expected one of {entrypoints})"


@pytest.mark.parametrize("name,path", _visible_py_plugins(), ids=[n for n, _ in _visible_py_plugins()])
def test_plugin_no_stdlib_import_errors(name, path):
    """Plugin does not import missing stdlib modules at top level."""
    source = path.read_text(encoding="utf-8", errors="ignore")
    tree = ast.parse(source, filename=str(path))

    # Collect top-level (non-try-wrapped) imports
    stdlib_failures = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Import, ast.ImportFrom)):
            continue
        # Only check stdlib — skip known optional third-party deps
        optional = {"PyQt6", "rdkit", "pyvista", "numpy", "scipy", "PIL",
                    "google", "openai", "anthropic", "truststore", "pymatgen",
                    "openbabel", "pyscf", "ase", "moleditpy", "sip", "modules",
                    "moleditpy_linux", "pybel"}
        if isinstance(node, ast.ImportFrom):
            top = (node.module or "").split(".")[0]
        else:
            top = node.names[0].name.split(".")[0]
        if top in optional:
            continue
        spec = importlib.util.find_spec(top) if top else None
        if spec is None and top:
            stdlib_failures.append(top)

    assert not stdlib_failures, f"{name}: missing stdlib modules: {stdlib_failures}"
