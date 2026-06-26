"""
Integration tests: plugin initialize() with the real PluginContext class.

Auto-clones ../python_molecular_editor if absent (--depth 1, read-only) so
these tests always run in CI without any extra workflow steps.

How this differs from existing tests (tests/test_initialize.py and
tests/test_menu_registration.py)
--------------------------------------------------------------------
Both of those test suites pass a pure ``MagicMock`` as the context.  MagicMock
accepts any attribute access and any call signature without raising, so a
plugin that calls a non-existent PluginContext method, passes a non-string path,
or passes a non-callable as a callback would silently pass those tests.

Here the *real* ``PluginContext`` class is used.  Only the manager it wraps is
a MagicMock (we cannot spin up a full MainWindow headlessly).  This means:

- A plugin calling ``ctx.nonexistent_method()`` → AttributeError (real bug).
- A plugin passing a non-string to ``add_menu_action(path, ...)`` → caught.
- A plugin passing ``None`` instead of a callable → caught.

The tests in section 1 exercise all visible ``initialize()`` plugins against
the real API surface.  Sections 2–4 are unit tests of PluginContext itself that
verify the delegation contracts documented in plugin_interface.py — things that
tests/test_api.py (static AST analysis) cannot verify at runtime.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from conftest import (
    load_plugin_for_gui,
    make_stub_manager,
    mock_chemistry_imports,
    visible_py_plugins,
)

# ---------------------------------------------------------------------------
# Helper: extract (path, callback) pairs from add_menu_action calls
# ---------------------------------------------------------------------------

_MGR_MENU_METHOD = "register_menu_action"

# Exempt plugins that intentionally register nothing in initialize()
# (mirrored from test_menu_registration.py — must stay in sync)
_EXEMPT: frozenset[str] = frozenset(
    {
        "Dark Mode Theme",
        "All-Trans Optimizer",
        "Complex Molecule Untangler",
        "PubChem Name Resolver",
        "Vector Viewer",
        "Python Console",
    }
)

# Plugins whose initialize() passes the main window as a Qt parent (QTimer, QWidget, …)
# and therefore require a real QObject — not a MagicMock — to construct cleanly when
# Qt is real.  These pass fine in tests/ (Qt is mocked there) but cannot be exercised
# in the GUI test suite without a full live MainWindow.
_EXEMPT_NEEDS_REAL_MW: frozenset[str] = frozenset(
    {
        "Structural Updater",   # QTimer(self.mw) in __init__
        "Advanced Rendering",   # QWidget(parent=mw) in __init__
    }
)

# Combined exclusion set for all three parametrized integration test groups
_SKIP_INIT = _EXEMPT | _EXEMPT_NEEDS_REAL_MW

_INIT_PLUGINS = visible_py_plugins(entry_point="initialize")


def _menu_calls(mgr: Any) -> list[tuple[str, Any]]:
    """
    Return [(path, callback), …] from mgr.register_menu_action call records.

    PluginContext.add_menu_action(path, cb) calls:
        mgr.register_menu_action(plugin_name, path, cb, text, icon, shortcut)
    so args[1] = path, args[2] = cb.
    """
    results = []
    for c in mgr.register_menu_action.call_args_list:
        args = c.args
        if len(args) >= 3:
            results.append((args[1], args[2]))
    return results


# ===========================================================================
# 1. initialize() with the real PluginContext — catches API surface misuse
# ===========================================================================


@pytest.mark.parametrize(
    "name,path",
    [(n, p) for n, p in _INIT_PLUGINS if n not in _SKIP_INIT],
    ids=[n for n, _ in _INIT_PLUGINS if n not in _SKIP_INIT],
)
def test_initialize_real_context_no_exception(
    name: str,
    path: Path,
    real_plugin_context_class: type,
    qapp: Any,
) -> None:
    """
    initialize(ctx) must not raise when ctx is a *real* PluginContext instance.

    This complements test_initialize.py (which uses MagicMock) by catching a
    different class of bug: calling a PluginContext method that does not exist,
    passing the wrong number of arguments, or passing a non-callable where a
    callable is required.  MagicMock would silently swallow all of those.
    """
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, name)

    with mock_chemistry_imports():
        mod = load_plugin_for_gui(path)
        mod.initialize(ctx)  # type: ignore[attr-defined]


@pytest.mark.parametrize(
    "name,path",
    [(n, p) for n, p in _INIT_PLUGINS if n not in _SKIP_INIT],
    ids=[n for n, _ in _INIT_PLUGINS if n not in _SKIP_INIT],
)
def test_menu_action_path_is_nonempty_string(
    name: str,
    path: Path,
    real_plugin_context_class: type,
    qapp: Any,
) -> None:
    """
    Each path passed to ctx.add_menu_action() must be a non-empty string.

    Passing None, an integer, or "" causes the host app menu-builder to crash
    silently.  A MagicMock context cannot catch this; the real PluginContext
    records the exact values delegated to the manager, which we inspect here.

    Plugins that register nothing via add_menu_action() (e.g. file-opener-only
    plugins) are skipped — they are covered by test_menu_registration.py.
    """
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, name)

    with mock_chemistry_imports():
        mod = load_plugin_for_gui(path)
        mod.initialize(ctx)  # type: ignore[attr-defined]

    calls = _menu_calls(mgr)
    if not calls:
        pytest.skip(f"{name!r} uses a non-menu registration path; skipping path check")

    for menu_path, _ in calls:
        assert isinstance(menu_path, str) and menu_path, (
            f"Plugin {name!r}: add_menu_action() called with invalid path {menu_path!r} "
            f"(must be a non-empty str)"
        )


@pytest.mark.parametrize(
    "name,path",
    [(n, p) for n, p in _INIT_PLUGINS if n not in _SKIP_INIT],
    ids=[n for n, _ in _INIT_PLUGINS if n not in _SKIP_INIT],
)
def test_menu_action_callback_is_callable(
    name: str,
    path: Path,
    real_plugin_context_class: type,
    qapp: Any,
) -> None:
    """
    Each callback passed to ctx.add_menu_action() must be callable.

    Passing a string or None instead of a function is caught here because the
    real PluginContext faithfully records the value passed to the manager.
    A MagicMock context records the call too, but existing tests don't inspect
    the argument types.

    Plugins with no add_menu_action() calls are skipped.
    """
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, name)

    with mock_chemistry_imports():
        mod = load_plugin_for_gui(path)
        mod.initialize(ctx)  # type: ignore[attr-defined]

    calls = _menu_calls(mgr)
    if not calls:
        pytest.skip(f"{name!r} uses a non-menu registration path; skipping callback check")

    for menu_path, cb in calls:
        assert callable(cb), (
            f"Plugin {name!r}: callback registered for {menu_path!r} is not callable: {cb!r}"
        )


# ===========================================================================
# 2. PluginContext property contracts — runtime verification
# ===========================================================================
# These tests exercise the real PluginContext class itself (not plugins).
# They are not present in any existing test file because tests/ uses only
# MagicMock contexts and tests/test_api.py is static AST analysis.


def test_ctx_current_mol_returns_none_without_main_window(
    real_plugin_context_class: type,
) -> None:
    """ctx.current_mol returns None when get_main_window() returns None."""
    ctx = real_plugin_context_class(make_stub_manager(main_window=None), "t")
    assert ctx.current_mol is None


def test_ctx_current_molecule_alias_equals_current_mol(
    real_plugin_context_class: type,
) -> None:
    """ctx.current_molecule (compat alias) returns the same value as current_mol."""
    ctx = real_plugin_context_class(make_stub_manager(main_window=None), "t")
    assert ctx.current_molecule is ctx.current_mol


def test_ctx_plotter_returns_none_without_main_window(
    real_plugin_context_class: type,
) -> None:
    """ctx.plotter returns None when get_main_window() returns None."""
    ctx = real_plugin_context_class(make_stub_manager(main_window=None), "t")
    assert ctx.plotter is None


def test_ctx_get_main_window_delegates(real_plugin_context_class: type) -> None:
    """ctx.get_main_window() must delegate to manager.get_main_window()."""
    sentinel = object()
    mgr = make_stub_manager(main_window=sentinel)
    ctx = real_plugin_context_class(mgr, "t")
    assert ctx.get_main_window() is sentinel


# ===========================================================================
# 3. PluginContext delegation contracts — add_* / register_* methods
# ===========================================================================


def test_add_menu_action_delegates_with_plugin_name(
    real_plugin_context_class: type,
) -> None:
    """add_menu_action(path, cb) must forward (plugin_name, path, cb, …) to manager."""
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, "MyPlugin")
    cb = lambda: None  # noqa: E731
    ctx.add_menu_action("File/Test...", cb)
    args = mgr.register_menu_action.call_args.args
    assert args[0] == "MyPlugin"
    assert args[1] == "File/Test..."
    assert args[2] is cb


def test_add_plugin_menu_prepends_plugin_prefix(
    real_plugin_context_class: type,
) -> None:
    """add_plugin_menu('Sub/Item', cb) must register under 'Plugin/Sub/Item'."""
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, "MyPlugin")
    cb = lambda: None  # noqa: E731
    ctx.add_plugin_menu("Sub/Item", cb)
    args = mgr.register_menu_action.call_args.args
    assert args[1] == "Plugin/Sub/Item"


def test_register_window_namespaces_by_plugin(
    real_plugin_context_class: type,
) -> None:
    """register_window('panel', w) must call manager.register_window(plugin_name, 'panel', w)."""
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, "MyPlugin")
    win = object()
    ctx.register_window("panel", win)
    mgr.register_window.assert_called_once_with("MyPlugin", "panel", win)


def test_show_status_message_delegates_args(
    real_plugin_context_class: type,
) -> None:
    """show_status_message(msg, timeout) must pass both args to manager unchanged."""
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, "t")
    ctx.show_status_message("hello", 1500)
    mgr.show_status_message.assert_called_once_with("hello", 1500)


def test_add_export_action_delegates(real_plugin_context_class: type) -> None:
    """add_export_action(label, cb) must call manager.register_export_action(name, label, cb)."""
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, "Exporter")
    cb = lambda: None  # noqa: E731
    ctx.add_export_action("Export as FOO...", cb)
    mgr.register_export_action.assert_called_once_with("Exporter", "Export as FOO...", cb)


def test_register_file_opener_delegates(real_plugin_context_class: type) -> None:
    """register_file_opener(ext, cb, priority) must forward all args to manager."""
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, "Opener")
    cb = lambda p: None  # noqa: E731
    ctx.register_file_opener(".foo", cb, priority=5)
    mgr.register_file_opener.assert_called_once_with("Opener", ".foo", cb, 5)


# ===========================================================================
# 4. Backward-compat alias: register_menu_action (old 3-arg style)
# ===========================================================================


def test_register_menu_action_compat_new_style(
    real_plugin_context_class: type,
) -> None:
    """register_menu_action(path, callback) — new 2-arg style — must work."""
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, "P")
    cb = lambda: None  # noqa: E731
    ctx.register_menu_action("A/B", cb)
    args = mgr.register_menu_action.call_args.args
    assert args[1] == "A/B"
    assert args[2] is cb


def test_register_menu_action_compat_old_style(
    real_plugin_context_class: type,
) -> None:
    """register_menu_action(path, text, callback) — old 3-arg style — must work."""
    mgr = make_stub_manager()
    ctx = real_plugin_context_class(mgr, "P")
    cb = lambda: None  # noqa: E731
    ctx.register_menu_action("A/B", "Label", cb)
    args = mgr.register_menu_action.call_args.args
    assert args[1] == "A/B"
    assert args[2] is cb
