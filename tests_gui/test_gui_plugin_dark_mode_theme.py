"""
GUI-tier tests for the Dark Mode Theme plugin.

These tests run with real PyQt6 (QT_QPA_PLATFORM=offscreen) and add coverage
that the tests/ suite cannot provide — in particular tests that instantiate
real Qt objects (QWidget) to exercise plugin side-effects.

Dark Mode Theme is entirely absent from tests/.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

DARK_MODE_PATH = PLUGINS_DIR / "Dark_Mode_Theme" / "dark_mode_plugin.py"

with mock_chemistry_imports():
    _dark = load_plugin_for_gui(DARK_MODE_PATH)


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
