"""
Verify that every visible plugin's initialize() registers at least one
callable with PluginContext.

For each plugin that exposes initialize(), we call it with a stub context
(MagicMock that mimics PluginContext) and assert that at least one of the
standard registration methods was called.  This catches the silent bug where
a plugin sets up internal state and dialogs but forgets to wire the menu
entry with context.add_menu_action() or equivalent.

Exempt list:
- "Dark Mode Theme": uses autorun() to apply a stylesheet at load time;
  no menu entry is expected, which is intentional.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin, make_context, mock_optional_imports, visible_py_plugins

# All PluginContext V4 registration methods we consider valid proof that a
# plugin wired itself into the host application.
_REGISTRATION_ATTRS: frozenset[str] = frozenset({
    "add_menu_action",
    "add_export_action",
    "add_analysis_tool",
    "add_plugin_menu",
    "register_file_opener",
    "register_3d_style",
    "register_save_handler",
    "register_load_handler",
    "register_document_reset_handler",
    "add_toolbar_action",
})

# Plugins that intentionally perform no PluginContext menu registration in
# initialize() because they work via an alternative mechanism.
#
# run()-auto-registered plugins: the host app automatically adds these to the
# Plugin menu when it finds a run(main_window) function — calling add_menu_action()
# inside initialize() would create a duplicate entry.
#
#   - "Dark Mode Theme":              autorun() applies a stylesheet at load time;
#                                     no menu entry is expected.
#   - "Plugin Installer":             uses run(main_window) for its UI; initialize()
#                                     only schedules background checks.
#   - "All-Trans Optimizer":          run(mw) auto-registered; initialize() stores
#                                     the launch function, nothing more.
#   - "Complex Molecule Untangler":   same as above.
#   - "PubChem Name Resolver":        same as above.
#   - "Vector Viewer":                same as above; initialize() calls
#                                     show_status_message() but no menu registration.
#   - "Python Console":               run(mw) auto-registered; initialize() only stores context.
_EXEMPT: frozenset[str] = frozenset({
    "Dark Mode Theme",
    "Plugin Installer",
    "All-Trans Optimizer",
    "Complex Molecule Untangler",
    "PubChem Name Resolver",
    "Vector Viewer",
    "Python Console",
})

_PLUGINS = visible_py_plugins(entry_point="initialize")


@pytest.mark.parametrize(
    "name,path",
    _PLUGINS,
    ids=[n for n, _ in _PLUGINS],
)
def test_initialize_registers_something(name: str, path: Path) -> None:
    """
    initialize(context) must call at least one PluginContext registration
    method so the plugin is actually reachable by the user.
    """
    if name in _EXEMPT:
        pytest.skip(f"{name!r} intentionally uses autorun() — no context registration")

    with mock_optional_imports():
        mod = load_plugin(path)
        ctx = make_context()
        mod.initialize(ctx)  # type: ignore[attr-defined]

    called = {attr for attr in _REGISTRATION_ATTRS if getattr(ctx, attr).called}
    assert called, (
        f"Plugin {name!r} ({path.name}): initialize() did not call any of "
        f"{sorted(_REGISTRATION_ATTRS)}.  "
        f"Add context.add_menu_action(...) or another registration call."
    )
