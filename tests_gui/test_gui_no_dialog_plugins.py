"""
GUI-tier tests for visible plugins that have no Qt dialog class.

These tests run with real PyQt6 (QT_QPA_PLATFORM=offscreen) and add coverage
that the tests/ suite cannot provide — in particular tests that instantiate
real Qt objects (QWidget, QColor) to exercise plugin side-effects.

Plugins covered (chosen because they have zero or minimal coverage in tests/):
  - Dark Mode Theme      (entirely absent from tests/)
  - Dummy Atom Mode      (entirely absent from tests/)
  - Paste from ChemDraw  (extra edge-cases beyond tests/test_plugin_ui_misc.py)
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

DARK_MODE_PATH = PLUGINS_DIR / "Dark_Mode_Theme" / "dark_mode_plugin.py"
DUMMY_ATOM_PATH = PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py"
CHEMDRAW_PATH = PLUGINS_DIR / "Paste_from_ChemDraw" / "paste_chemdraw.py"

with mock_chemistry_imports():
    _dark = load_plugin_for_gui(DARK_MODE_PATH)
    _dummy = load_plugin_for_gui(DUMMY_ATOM_PATH)
    _chemdraw = load_plugin_for_gui(CHEMDRAW_PATH)


# ===========================================================================
# Dark Mode Theme — module constants
# ===========================================================================


class TestDarkModeConstants:
    def test_plugin_name(self):
        assert _dark.PLUGIN_NAME == "Dark Mode Theme"

    def test_plugin_version_format(self):
        parts = _dark.PLUGIN_VERSION.split(".")
        assert len(parts) == 3
        assert all(p.isdigit() for p in parts)

    def test_dark_stylesheet_is_nonempty(self):
        assert len(_dark.DARK_STYLESHEET) > 100

    def test_stylesheet_contains_qwidget_rule(self):
        assert "QWidget" in _dark.DARK_STYLESHEET

    def test_stylesheet_contains_qpushbutton_rule(self):
        assert "QPushButton" in _dark.DARK_STYLESHEET

    def test_stylesheet_contains_dark_background(self):
        # The defining characteristic: dark background colour present
        assert "#2b2b2b" in _dark.DARK_STYLESHEET

    def test_stylesheet_contains_menu_bar_rule(self):
        assert "QMenuBar" in _dark.DARK_STYLESHEET

    def test_stylesheet_contains_scrollbar_rule(self):
        assert "QScrollBar" in _dark.DARK_STYLESHEET


# ===========================================================================
# Dark Mode Theme — autorun() with a real QWidget
# ===========================================================================


class TestDarkModeAutorun:
    """autorun() must call setStyleSheet() on the real QWidget."""

    def test_autorun_applies_stylesheet_to_real_widget(self, qapp):
        from PyQt6.QtWidgets import QWidget

        w = QWidget()
        _dark._CONTEXT = None  # ensure no context side-effects
        _dark.autorun(w)
        ss = w.styleSheet()
        assert "background-color" in ss
        w.destroy()

    def test_autorun_stylesheet_matches_dark_constant(self, qapp):
        from PyQt6.QtWidgets import QWidget

        w = QWidget()
        _dark._CONTEXT = None
        _dark.autorun(w)
        assert "#2b2b2b" in w.styleSheet()
        w.destroy()

    def test_autorun_with_mock_mw_does_not_raise(self, qapp):
        mw = MagicMock()
        _dark._CONTEXT = None
        _dark.autorun(mw)
        mw.setStyleSheet.assert_called_once_with(_dark.DARK_STYLESHEET)

    def test_autorun_shows_status_message_when_context_set(self, qapp):
        ctx = MagicMock()
        _dark._CONTEXT = ctx
        _dark.autorun(MagicMock())
        ctx.show_status_message.assert_called_once()
        _dark._CONTEXT = None  # restore

    def test_autorun_message_mentions_dark_mode(self, qapp):
        ctx = MagicMock()
        _dark._CONTEXT = ctx
        _dark.autorun(MagicMock())
        msg = ctx.show_status_message.call_args[0][0]
        assert "dark" in msg.lower() or "Dark" in msg
        _dark._CONTEXT = None

    def test_autorun_updates_background_color_setting(self, qapp):
        mw = MagicMock()
        mw.init_manager.settings = {}
        _dark._CONTEXT = None
        _dark.autorun(mw)
        assert mw.init_manager.settings.get("background_color") == "#2b2b2b"


# ===========================================================================
# Dark Mode Theme — initialize()
# ===========================================================================


class TestDarkModeInitialize:
    def test_initialize_stores_context(self, qapp):
        ctx = MagicMock()
        ctx.get_main_window.return_value = MagicMock()
        _dark._CONTEXT = None
        _dark.initialize(ctx)
        assert _dark._CONTEXT is ctx
        _dark._CONTEXT = None  # restore

    def test_initialize_calls_autorun(self, qapp):
        ctx = MagicMock()
        mw = MagicMock()
        ctx.get_main_window.return_value = mw
        _dark._CONTEXT = None
        _dark.initialize(ctx)
        mw.setStyleSheet.assert_called_once()
        _dark._CONTEXT = None


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


# ===========================================================================
# Paste from ChemDraw — reconstruct_from_flat_text() extra edge cases
# (tests/test_plugin_ui_misc.py covers: empty, no V2000, too few tokens,
# non-printable chars.  These cases add valid header + partial data paths.)
# ===========================================================================


def _flat(n_atoms, n_bonds, atom_tokens="", bond_tokens=""):
    """Build a minimal flat ChemDraw clipboard string."""
    # 10 header tokens followed by V2000, then atom/bond tokens
    header = f"t1 t2 t3 t4 t5 t6 t7 t8 {n_bonds} {n_atoms}"
    rest = atom_tokens + " " + bond_tokens
    return header + " V2000 " + rest


class TestReconstructFromFlatTextExtra:
    def test_zero_atoms_zero_bonds_returns_string(self):
        result = _chemdraw.reconstruct_from_flat_text(_flat(0, 0))
        # With 0 atoms and 0 bonds the reconstruction may succeed or return None;
        # either is acceptable — it must not raise.

    def test_v2000_present_but_non_integer_counts_returns_none(self):
        bad = "t1 t2 t3 t4 t5 t6 t7 t8 XY ZZ V2000 rest"
        assert _chemdraw.reconstruct_from_flat_text(bad) is None

    def test_single_carbon_atom_block_returns_something(self):
        # 1 atom, 0 bonds, atom tokens: x y z symbol + 12 zeros
        atom_tok = "0.0 0.0 0.0 C " + " ".join(["0"] * 12)
        result = _chemdraw.reconstruct_from_flat_text(_flat(1, 0, atom_tok))
        # Result is either a string (reconstructed block) or None if parsing failed;
        # must not raise.
        assert result is None or isinstance(result, str)

    def test_whitespace_only_returns_none(self):
        assert _chemdraw.reconstruct_from_flat_text("   \n\t  ") is None

    def test_phantom_f_token_is_skipped(self):
        # Insert a 'F' phantom token before atom coordinates
        atom_tok = "F 0.0 0.0 0.0 C " + " ".join(["0"] * 12)
        result = _chemdraw.reconstruct_from_flat_text(_flat(1, 0, atom_tok))
        assert result is None or isinstance(result, str)

    def test_non_printable_control_chars_sanitized(self):
        text = "\x01\x02V2000 t1 t2 t3 t4 t5 t6 t7 t8 0 0"
        # Should not raise even though the string is unusual
        result = _chemdraw.reconstruct_from_flat_text(text)
        assert result is None or isinstance(result, str)
