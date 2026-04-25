"""
Smoke tests: call run(main_window) / autorun(main_window) on every visible
plugin that uses the legacy entry-point API (no initialize()).

These plugins pre-date the PluginContext API and receive the raw MainWindow
object directly.  The test verifies that calling the entry point with a
MagicMock main window does not raise an exception.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin, make_context, mock_optional_imports, visible_py_plugins

# Collect plugins that have run() or autorun() but NOT initialize().
# visible_py_plugins(entry_point=X) returns plugins that define X; we want
# run/autorun plugins that are not also initialize() plugins.

_INIT_NAMES = {name for name, _ in visible_py_plugins(entry_point="initialize")}

_RUN_PLUGINS: list[tuple[str, Path, str]] = []
for _name, _path in visible_py_plugins():
    if _name in _INIT_NAMES:
        continue
    import ast as _ast

    _src = _path.read_text(encoding="utf-8", errors="ignore")
    try:
        _tree = _ast.parse(_src)
    except SyntaxError:
        continue
    _fns = {
        node.name
        for node in _ast.walk(_tree)
        if isinstance(node, (_ast.FunctionDef, _ast.AsyncFunctionDef))
    }
    if "autorun" in _fns:
        _RUN_PLUGINS.append((_name, _path, "autorun"))
    elif "run" in _fns:
        _RUN_PLUGINS.append((_name, _path, "run"))


@pytest.mark.parametrize(
    "name,path,entry",
    _RUN_PLUGINS,
    ids=[n for n, _, _ in _RUN_PLUGINS],
)
def test_run_does_not_raise(name: str, path: Path, entry: str) -> None:
    """``run(main_window)`` / ``autorun(main_window)`` completes without raising."""
    with mock_optional_imports():
        mod = load_plugin(path)
        ctx = make_context()
        main_window = ctx.get_main_window()
        fn = getattr(mod, entry)
        fn(main_window)  # type: ignore[attr-defined]
