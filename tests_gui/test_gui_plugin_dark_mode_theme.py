"""
GUI-tier tests for the Dark Mode Theme plugin.

These tests run with real PySide6 (QT_QPA_PLATFORM=offscreen) and add coverage
that the tests/ suite cannot provide — in particular tests that instantiate
real Qt objects (QWidget) to exercise plugin side-effects.

Regression coverage for settings preservation and in-memory icon tinting.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import MappingProxyType, SimpleNamespace
from unittest.mock import MagicMock

import pytest
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
        from PySide6.QtWidgets import QWidget

        w = QWidget()
        _dark._CONTEXT = None  # ensure no context side-effects
        _dark.autorun(w)
        ss = w.styleSheet()
        assert "background-color" in ss
        w.destroy()

    def test_autorun_stylesheet_matches_dark_constant(self, qapp):
        from PySide6.QtWidgets import QWidget

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

    def test_autorun_preserves_background_color_setting(self, qapp):
        mw = MagicMock()
        mw.init_manager.settings = {"background_color": "#000000"}
        _dark._CONTEXT = None
        _dark.autorun(mw)
        assert mw.init_manager.settings == {"background_color": "#000000"}


@pytest.fixture
def themed_window(qapp):
    from PySide6.QtWidgets import QMainWindow

    window = QMainWindow()
    window.init_manager = SimpleNamespace(
        settings=MappingProxyType(
            {"background_color": "#000000", "icon_foreground": "#123456"}
        ),
        mode_actions={},
        save_settings=MagicMock(),
    )
    _dark._CONTEXT = None
    yield window
    qapp.processEvents()
    window.close()
    window.deleteLater()
    qapp.processEvents()
    _dark._CONTEXT = None


def _icon(color="#123456"):
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor, QIcon, QPainter, QPixmap

    icon = QIcon()
    for size in (32, 64):
        pixmap = QPixmap(size, size)
        pixmap.fill(Qt.GlobalColor.transparent)
        painter = QPainter(pixmap)
        painter.fillRect(4, 4, 12, 12, QColor(color))
        painter.end()
        icon.addPixmap(pixmap)
    return icon


class TestRuntimeIconTint:
    def test_tints_icons_created_after_initialize_without_changing_saved_settings(
        self, qapp, themed_window, tmp_path
    ):
        from PySide6.QtGui import QAction, QIcon

        window = themed_window
        settings_before = dict(window.init_manager.settings)
        settings_path = tmp_path / "settings.json"
        settings_path.write_text(json.dumps(settings_before))
        # The host has not created its tool actions when plugins initialize.
        _dark.autorun(window)
        originals = []
        for mode in _dark._TOOL_ICON_MODES:
            action = QAction(window)
            action.setIcon(_icon())
            window.init_manager.mode_actions[mode] = action
            originals.append(action.icon().cacheKey())
        custom = QAction(window)
        custom.setIcon(_icon("#FF0000"))
        custom_key = custom.icon().cacheKey()
        window.init_manager.mode_actions["plugin_custom"] = custom
        window.init_manager.mode_actions["atom_C"] = QAction("C", window)
        assert (
            window.init_manager.mode_actions["chain"]
            .icon()
            .pixmap(32, 32)
            .toImage()
            .pixelColor(8, 8)
            .name()
            == "#123456"
        )

        qapp.processEvents()

        assert len(originals) == 16
        for mode in _dark._TOOL_ICON_MODES:
            action = window.init_manager.mode_actions[mode]
            assert action.icon().cacheKey() not in originals
            for size in (32, 64):
                image = action.icon().pixmap(size, size).toImage()
                assert image.pixelColor(8, 8).name() == "#ffffff"
                assert image.pixelColor(8, 8).alpha() == 255
                assert image.pixelColor(0, 0).alpha() == 0
            checked = (
                action.icon()
                .pixmap(32, 32, QIcon.Mode.Normal, QIcon.State.On)
                .toImage()
            )
            assert checked.pixelColor(8, 8).name() == "#ffffff"
            disabled = action.icon().pixmap(32, 32, QIcon.Mode.Disabled).toImage()
            assert disabled.pixelColor(8, 8).alpha() > 0
        assert custom.icon().cacheKey() == custom_key
        assert window.init_manager.mode_actions["atom_C"].icon().isNull()
        assert dict(window.init_manager.settings) == settings_before
        assert json.loads(settings_path.read_text()) == settings_before
        window.init_manager.save_settings.assert_not_called()

    def test_reapplying_is_safe_and_does_not_add_icon_setting(
        self, qapp, themed_window
    ):
        from PySide6.QtGui import QAction

        window = themed_window
        window.init_manager.settings = MappingProxyType({"background_color": "#F0F0F0"})
        action = QAction(window)
        action.setIcon(_icon())
        window.init_manager.mode_actions["chain"] = action
        for _ in range(2):
            _dark.autorun(window)
            qapp.processEvents()
            assert (
                action.icon().pixmap(32, 32).toImage().pixelColor(8, 8).name()
                == "#ffffff"
            )
            assert dict(window.init_manager.settings) == {"background_color": "#F0F0F0"}

    def test_null_and_unavailable_icons_are_left_alone(self, themed_window):
        from PySide6.QtGui import QAction

        empty = QAction(themed_window)
        themed_window.init_manager.mode_actions = {"chain": empty}
        _dark._tint_tool_icons(themed_window)
        assert empty.icon().isNull()
        themed_window.init_manager.mode_actions = None
        _dark._tint_tool_icons(themed_window)

    def test_invalid_action_does_not_prevent_other_icons(self, themed_window, caplog):
        from PySide6.QtGui import QAction

        broken = MagicMock()
        broken.icon.side_effect = RuntimeError("deleted action")
        valid = QAction(themed_window)
        valid.setIcon(_icon())
        themed_window.init_manager.mode_actions = {"bond_1_0": broken, "chain": valid}
        _dark._tint_tool_icons(themed_window)
        assert "Could not tint tool icon bond_1_0" in caplog.text
        assert (
            valid.icon().pixmap(32, 32).toImage().pixelColor(8, 8).name() == "#ffffff"
        )

    def test_closed_window_cancels_pending_tint(self, qapp, monkeypatch):
        from PySide6.QtCore import QCoreApplication, QEvent
        from PySide6.QtWidgets import QMainWindow

        window = QMainWindow()
        callback = MagicMock()
        monkeypatch.setattr(_dark, "_tint_tool_icons", callback)
        _dark._CONTEXT = None
        _dark.autorun(window)
        window.deleteLater()
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        qapp.processEvents()
        callback.assert_not_called()


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
