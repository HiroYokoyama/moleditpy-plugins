"""
Tests for the Dummy Atom Mode plugin (initialize -> add_toolbar_action, no add_menu_action).
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
DUMMY_PATH = PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py"


class TestDummyAtomMode:
    def test_initialize_calls_add_toolbar_action(self):
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_toolbar_action.assert_called_once()

    def test_initialize_toolbar_text_is_dummy_atom_star(self):
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            kwargs = ctx.add_toolbar_action.call_args[1]
            assert kwargs.get("text") == "Dummy Atom *"

    def test_initialize_does_not_call_add_menu_action(self):
        """Dummy Atom Mode registers a toolbar button, NOT a menu item."""
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_not_called()

    def test_initialize_tooltip_mentions_dummy(self):
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            kwargs = ctx.add_toolbar_action.call_args[1]
            assert "dummy" in kwargs.get("tooltip", "").lower()
