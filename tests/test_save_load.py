"""
Round-trip tests for plugin save/load handlers.

For every visible plugin that:
  1. exposes ``initialize(context)``, AND
  2. calls ``context.register_save_handler(fn)`` during initialize(),

this test captures the registered handler and exercises it:

  save_data = save_handler()          # must not raise
  load_handler(save_data or {})       # must not raise

This validates the persistence code path — data serialisation and
deserialisation logic that runs when the user saves/loads a project file.
The ``save_handler`` may legitimately return ``None`` (meaning "nothing to
save right now"); in that case ``load_handler`` is called with an empty dict,
which every well-behaved load handler must tolerate gracefully.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin, make_context, mock_optional_imports, visible_py_plugins

# Only parametrize plugins that have initialize() — others don't use the
# PluginContext save/load API.
_INIT_PLUGINS = visible_py_plugins(entry_point="initialize")


def _collect() -> list[tuple[str, Path]]:
    """
    Return plugins that actually register a save handler at initialize() time.

    We do a dry-run initialize() here (at collection time) to check which
    plugins call register_save_handler.  Plugins that don't register one are
    excluded so we don't generate vacuous test cases.
    """
    result = []
    for name, path in _INIT_PLUGINS:
        try:
            with mock_optional_imports():
                mod = load_plugin(path)
                ctx = make_context()
                mod.initialize(ctx)  # type: ignore[attr-defined]
                if ctx.register_save_handler.called:
                    result.append((name, path))
        except Exception:
            # If initialize() itself fails here, test_initialize.py will
            # already catch it; skip silently during collection.
            pass
    return result


_SAVE_LOAD_PLUGINS = _collect()


@pytest.mark.parametrize(
    "name,path", _SAVE_LOAD_PLUGINS, ids=[n for n, _ in _SAVE_LOAD_PLUGINS]
)
def test_save_handler_does_not_raise(name: str, path: Path) -> None:
    """``save_handler()`` completes without raising."""
    with mock_optional_imports():
        mod = load_plugin(path)
        ctx = make_context()
        mod.initialize(ctx)  # type: ignore[attr-defined]
        save_fn = ctx.register_save_handler.call_args[0][0]
        save_fn()  # must not raise


@pytest.mark.parametrize(
    "name,path", _SAVE_LOAD_PLUGINS, ids=[n for n, _ in _SAVE_LOAD_PLUGINS]
)
def test_load_handler_does_not_raise(name: str, path: Path) -> None:
    """``load_handler(save_data or {})`` completes without raising."""
    with mock_optional_imports():
        mod = load_plugin(path)
        ctx = make_context()
        mod.initialize(ctx)  # type: ignore[attr-defined]
        save_fn = ctx.register_save_handler.call_args[0][0]
        save_data = save_fn()

        if ctx.register_load_handler.called:
            load_fn = ctx.register_load_handler.call_args[0][0]
            # Use save_data if it's a dict; fall back to {} if None
            load_fn(save_data if isinstance(save_data, dict) else {})
