"""
Tests for the Dark Mode Theme plugin: stylesheet content, autorun settings mutation.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
DARK_PATH = PLUGINS_DIR / "Dark_Mode_Theme" / "dark_mode_plugin.py"


class TestDarkModeStylesheet:
    @pytest.fixture(scope="class")
    def mod(self):
        with mock_optional_imports():
            return load_plugin(DARK_PATH)

    @pytest.mark.parametrize("selector", [
        "QWidget", "QMainWindow", "QMenuBar", "QMenu", "QTabBar::tab",
        "QToolBar", "QPushButton", "QLineEdit", "QCheckBox", "QSlider",
        "QScrollBar:vertical", "QHeaderView::section", "QGroupBox",
        "QStatusBar", "QSplitter::handle", "QTableView",
    ])
    def test_stylesheet_covers_widget(self, mod, selector):
        assert selector in mod.DARK_STYLESHEET

    def test_stylesheet_uses_dark_background(self, mod):
        assert "#2b2b2b" in mod.DARK_STYLESHEET

    def test_stylesheet_uses_accent_color(self, mod):
        assert "#3a6ea5" in mod.DARK_STYLESHEET

    def test_stylesheet_has_balanced_braces(self, mod):
        qss = mod.DARK_STYLESHEET
        assert qss.count("{") == qss.count("}")
        assert qss.count("{") > 20


class TestDarkModeAutorun:
    def _fresh(self):
        with mock_optional_imports():
            return load_plugin(DARK_PATH)

    def test_applies_stylesheet_to_main_window(self):
        mod = self._fresh()
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        mw.setStyleSheet.assert_called_once_with(mod.DARK_STYLESHEET)

    def test_sets_dark_3d_background_in_settings(self):
        mod = self._fresh()
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        assert mw.init_manager.settings["background_color"] == "#2b2b2b"

    def test_sets_white_icon_foreground(self):
        mod = self._fresh()
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        assert mw.init_manager.settings["icon_foreground"] == "#FFFFFF"

    def test_triggers_apply_3d_settings(self):
        mod = self._fresh()
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        mw.view_3d_manager.apply_3d_settings.assert_called_once_with()

    def test_settings_error_does_not_raise(self):
        mod = self._fresh()
        mw = MagicMock()
        broken = MagicMock()
        broken.__setitem__.side_effect = RuntimeError("boom")
        mw.init_manager.settings = broken
        mod.autorun(mw)  # must be swallowed
        mw.setStyleSheet.assert_called_once()

    def test_status_message_via_context_after_initialize(self):
        mod = self._fresh()
        ctx = make_context()
        mod.initialize(ctx)
        ctx.get_main_window.return_value.setStyleSheet.assert_called_once()
        ctx.show_status_message.assert_called_once()
        msg = ctx.show_status_message.call_args[0][0]
        assert "Dark Mode" in msg

    def test_status_message_falls_back_to_statusbar(self):
        mod = self._fresh()
        assert mod._CONTEXT is None
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        mw.statusBar.return_value.showMessage.assert_called_once()
