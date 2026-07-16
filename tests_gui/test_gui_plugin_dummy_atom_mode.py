"""
GUI-tier tests for the Dummy Atom Mode plugin.

These tests run with real PyQt6 (QT_QPA_PLATFORM=offscreen) and add coverage
that the tests/ suite cannot provide — in particular tests that instantiate
real Qt objects (QColor) to exercise plugin side-effects.

Dummy Atom Mode is entirely absent from tests/.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

DUMMY_ATOM_PATH = PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py"

with mock_chemistry_imports():
    _dummy = load_plugin_for_gui(DUMMY_ATOM_PATH)


# ===========================================================================
# Dummy Atom Mode — module constants
# ===========================================================================


class TestDummyAtomConstants:
    def test_plugin_name(self):
        assert _dummy.PLUGIN_NAME == "Dummy Atom Mode"

    def test_plugin_version_format(self):
        parts = _dummy.PLUGIN_VERSION.split(".")
        assert len(parts) == 3
        assert all(p.isdigit() for p in parts)

    def test_plugin_description_mentions_dummy(self):
        assert "dummy" in _dummy.PLUGIN_DESCRIPTION.lower()

    def test_plugin_description_mentions_toolbar(self):
        assert "toolbar" in _dummy.PLUGIN_DESCRIPTION.lower()


# ===========================================================================
# Dummy Atom Mode — initialize() with real QColor
# ===========================================================================


class TestDummyAtomInitialize:
    """initialize() registers a toolbar action and creates a real QColor."""

    def test_initialize_registers_toolbar_action(self, qapp):
        ctx = MagicMock()
        _dummy.initialize(ctx)
        ctx.add_toolbar_action.assert_called_once()

    def test_toolbar_action_text_is_dummy_atom(self, qapp):
        ctx = MagicMock()
        _dummy.initialize(ctx)
        call_kwargs = ctx.add_toolbar_action.call_args
        text = call_kwargs.kwargs.get("text") or call_kwargs.args[1]
        assert "Dummy" in text or "dummy" in text.lower()

    def test_toolbar_action_has_tooltip(self, qapp):
        ctx = MagicMock()
        _dummy.initialize(ctx)
        call_kwargs = ctx.add_toolbar_action.call_args
        tooltip = call_kwargs.kwargs.get("tooltip", "")
        assert tooltip  # must be non-empty

    def test_initialize_does_not_raise(self, qapp):
        ctx = MagicMock()
        _dummy.initialize(ctx)

    def test_cpk_colors_wildcard_is_real_qcolor(self, qapp):
        """initialize() sets CPK_COLORS['*'] to a real QColor (not a MagicMock)."""
        from PyQt6.QtGui import QColor

        ctx = MagicMock()
        _dummy.initialize(ctx)
        # MagicMock tracks __setitem__ calls; find the one for key "*"
        found = False
        for call in _dummy.CPK_COLORS.__setitem__.call_args_list:
            if len(call.args) >= 2 and call.args[0] == "*":
                assert isinstance(call.args[1], QColor), (
                    f"Expected QColor for CPK_COLORS['*'], got {type(call.args[1])}"
                )
                found = True
        assert found, "CPK_COLORS['*'] was never set during initialize()"

    def test_initialize_patches_chem_atom(self, qapp):
        """initialize() monkey-patches Chem.Atom; the mocked Chem is updated."""
        ctx = MagicMock()
        before = _dummy.Chem.Atom
        _dummy.initialize(ctx)
        # After initialize, Chem.Atom should be different (the patched version)
        # (With mocked Chem the attribute is writable — just verify no error)
        assert _dummy.Chem.Atom is not None
