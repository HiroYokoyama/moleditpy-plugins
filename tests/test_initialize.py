"""
Smoke tests: call initialize(context) on every visible plugin that exposes it.

Verifies that initialize() completes without raising an exception when given
a realistic stub context (non-None main window, all heavy deps mocked out).

Rationale: AST-based tests (test_imports.py) catch syntax errors and missing
stdlib imports, but cannot catch runtime errors that only surface when
initialize() actually executes — e.g. missing attribute access, bad argument
count, or unguarded None dereferences.

All heavy optional dependencies (PyQt6, rdkit, moleditpy, etc.) are intercepted
by a custom MetaPathFinder so these tests run without any GUI or chemistry
libraries installed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin, make_context, mock_optional_imports, visible_py_plugins

_PLUGINS = visible_py_plugins(entry_point="initialize")


@pytest.mark.parametrize(
    "name,path", _PLUGINS, ids=[n for n, _ in _PLUGINS]
)
def test_initialize_does_not_raise(name: str, path: Path) -> None:
    """``initialize(context)`` completes without raising any exception."""
    with mock_optional_imports():
        mod = load_plugin(path)
        ctx = make_context()
        mod.initialize(ctx)  # type: ignore[attr-defined]
