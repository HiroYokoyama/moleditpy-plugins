"""
dark_mode_plugin.py

A plugin that automatically applies a dark mode stylesheet to MoleditPy upon loading.

2026.10.30: Preserve saved settings and tint toolbar icons only in memory.
"""

import logging
from collections.abc import Mapping
from typing import Any

from PySide6.QtCore import QTimer
from PySide6.QtGui import QColor, QIcon, QPainter

PLUGIN_VERSION = "2026.10.30"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=5.0.0, <6.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "A plugin that automatically applies a dark mode stylesheet to MoleditPy upon loading."

PLUGIN_NAME = "Dark Mode Theme"
_LOGGER = logging.getLogger(__name__)

_CONTEXT = None

_TOOL_ICON_MODES = (
    "bond_1_0",
    "bond_2_0",
    "bond_3_0",
    "bond_1_1",
    "bond_1_2",
    "rs_label",
    "bond_2_5",
    "template_benzene",
    "chain",
) + tuple(f"template_{n}" for n in range(3, 10))


def initialize(context: Any) -> None:
    global _CONTEXT
    _CONTEXT = context
    autorun(context.get_main_window())


# Dark Mode QSS
DARK_STYLESHEET = """
/* General Widget Defaults */
QWidget {
    background-color: #2b2b2b;
    color: #e0e0e0;
    selection-background-color: #3a6ea5;
    selection-color: #ffffff;
}

/* Main Window & Dialogs */
QMainWindow, QDialog {
    background-color: #2b2b2b;
}

/* Menu Bar */
QMenuBar {
    background-color: #323232;
    color: #e0e0e0;
}
QMenuBar::item {
    background-color: transparent;
    padding: 6px 10px;
}
QMenuBar::item:selected {
    background-color: #3a6ea5;
    color: #ffffff;
}
QMenuBar::item:pressed {
    background-color: #444444;
}

/* Menus */
QMenu {
    background-color: #323232;
    border: 1px solid #555555;
    padding: 4px;
}
QMenu::item {
    padding: 6px 24px 6px 12px;
}
QMenu::item:selected {
    background-color: #3a6ea5;
    color: #ffffff;
}
QMenu::separator {
    height: 1px;
    background-color: #555555;
    margin: 4px 0px;
}

/* ToolTabs & TabWidgets */
QTabWidget::pane {
    border: 1px solid #444444;
    background-color: #2b2b2b;
}
QTabBar::tab {
    background-color: #323232;
    color: #b0b0b0;
    padding: 6px 12px;
    border: 1px solid #444444;
    border-bottom: none;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
    margin-right: 2px;
}
QTabBar::tab:selected {
    background-color: #2b2b2b;
    color: #ffffff;
    border-bottom: 2px solid #3a6ea5; /* Highlight selected tab */
}
QTabBar::tab:hover {
    background-color: #3e3e3e;
}

/* ToolBar */
QToolBar {
    background-color: #323232;
    border-bottom: 1px solid #444444;
    spacing: 4px; 
    padding: 4px;
}
QToolButton {
    background-color: transparent;
    color: #e0e0e0;
    border-radius: 4px;
    padding: 4px;
}
QToolButton:hover {
    background-color: #3a6ea5;
}
QToolButton:pressed {
    background-color: #254e7a;
}
QToolButton:checked {
    background-color: #3a6ea5;
    border: 1px solid #555555;
}

/* Push Buttons */
QPushButton {
    background-color: #3c3c3c;
    border: 1px solid #555555;
    border-radius: 4px;
    color: #e0e0e0;
    padding: 6px 12px;
    min-width: 60px;
}
QPushButton:hover {
    background-color: #4a4a4a;
    border-color: #666666;
}
QPushButton:pressed {
    background-color: #2a2a2a;
}
QPushButton:disabled {
    background-color: #333333;
    color: #666666;
    border-color: #444444;
}

/* Input Fields (Text, Number, Combo) */
QLineEdit, QTextEdit, QPlainTextEdit, QAbstractSpinBox, QComboBox {
    background-color: #1e1e1e;
    color: #ffffff;
    border: 1px solid #444444;
    border-radius: 4px;
    padding: 4px;
    selection-background-color: #3a6ea5;
}
QLineEdit:focus, QTextEdit:focus, QPlainTextEdit:focus, QAbstractSpinBox:focus, QComboBox:focus {
    border: 1px solid #3a6ea5;
}
QComboBox::drop-down {
    border: none;
    width: 20px;
}
QComboBox QAbstractItemView {
    background-color: #1e1e1e;
    color: #ffffff;
    border: 1px solid #444444;
    selection-background-color: #3a6ea5;
}

/* CheckBox & RadioButton */
QCheckBox, QRadioButton {
    spacing: 6px;
    color: #e0e0e0;
}
QCheckBox::indicator, QRadioButton::indicator {
    width: 16px;
    height: 16px;
    background-color: #1e1e1e;
    border: 1px solid #555555;
    border-radius: 3px;
}
QCheckBox::indicator:hover, QRadioButton::indicator:hover {
    border-color: #3a6ea5;
}
QCheckBox::indicator:checked {
    background-color: #3a6ea5;
    border-color: #3a6ea5;
    image: url(none); /* Usually would use an icon, but color is enough for now */
}
/* Ensure radio button is round */
QRadioButton::indicator {
    border-radius: 8px;
}
QRadioButton::indicator:checked {
    background-color: #3a6ea5; 
}

/* Sliders */
QSlider::groove:horizontal {
    border: 1px solid #444444;
    height: 6px;
    background: #1e1e1e;
    margin: 2px 0;
    border-radius: 3px;
}
QSlider::handle:horizontal {
    background: #555555;
    border: 1px solid #666666;
    width: 14px;
    height: 14px;
    margin: -5px 0;
    border-radius: 7px;
}
QSlider::handle:horizontal:hover {
    background: #3a6ea5;
    border-color: #3a6ea5;
}

/* ScrollBars */
QScrollBar:vertical {
    border: none;
    background-color: #2b2b2b;
    width: 12px;
    margin: 0px 0px 0px 0px;
}
QScrollBar::handle:vertical {
    background-color: #444444;
    min-height: 20px;
    border-radius: 6px;
    margin: 2px;
}
QScrollBar::handle:vertical:hover {
    background-color: #555555;
}
QScrollBar:horizontal {
    border: none;
    background-color: #2b2b2b;
    height: 12px;
    margin: 0px 0px 0px 0px;
}
QScrollBar::handle:horizontal {
    background-color: #444444;
    min-width: 20px;
    border-radius: 6px;
    margin: 2px;
}
QScrollBar::handle:horizontal:hover {
    background-color: #555555;
}

/* Header View (for Tables) */
QHeaderView::section {
    background-color: #323232;
    color: #e0e0e0;
    padding: 4px;
    border: 1px solid #444444;
}

/* GroupBox */
QGroupBox {
    border: 1px solid #444444;
    border-radius: 4px;
    margin-top: 20px;
    padding-top: 10px;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top center;
    padding: 0 5px;
    background-color: #2b2b2b; /* Match window bg to hide border behind text */
    color: #3a6ea5; /* Accent color for titles */
    font-weight: bold;
}

/* StatusBar */
QStatusBar {
    background-color: #323232;
    border-top: 1px solid #444444;
    color: #a0a0a0;
}
QStatusBar QLabel {
    color: #a0a0a0;
}

/* Splitter */
QSplitter::handle {
    background-color: #444444;
}
QSplitter::handle:hover {
    background-color: #3a6ea5;
}

/* Table View */
QTableView {
    gridline-color: #444444;
    background-color: #1e1e1e;
    alternate-background-color: #252525;
}
"""


def _tint_tool_icons(main_window: Any) -> None:
    """Tint built-in tool icons without changing settings or third-party icons."""
    manager = getattr(main_window, "init_manager", None)
    actions = getattr(manager, "mode_actions", None)
    if not isinstance(actions, Mapping):
        return

    for mode in _TOOL_ICON_MODES:
        action = actions.get(mode)
        if action is None:
            continue
        try:
            original = action.icon()
            if original.isNull():
                continue
            tinted = QIcon()
            for size in original.availableSizes():
                pixmap = original.pixmap(size)
                if pixmap.isNull():
                    continue
                painter = QPainter(pixmap)
                try:
                    painter.setCompositionMode(
                        QPainter.CompositionMode.CompositionMode_SourceIn
                    )
                    painter.fillRect(pixmap.rect(), QColor("#FFFFFF"))
                finally:
                    painter.end()
                tinted.addPixmap(pixmap)
            if not tinted.isNull():
                action.setIcon(tinted)
        except (AttributeError, RuntimeError, TypeError):
            _LOGGER.warning(
                "[%s] Could not tint tool icon %s", PLUGIN_NAME, mode, exc_info=True
            )


def autorun(main_window: Any) -> None:
    """
    Executed automatically when the plugin is loaded (and 'autorun' is detected).
    """
    _LOGGER.info("[%s] Applying Dark Mode Stylesheet...", PLUGIN_NAME)

    # 1. Apply QSS to the main window (and thus all children)
    # Note: Applying to main_window is safer than QApplication for plugins to avoid
    # affecting other parts if this were embedded in a larger system, but for a
    # value-add plugin for this specific app, it works well.
    main_window.setStyleSheet(DARK_STYLESHEET)

    # Plugins initialize before the host creates its standard tool actions.
    # The window context cancels the callback if it closes before delivery.
    try:
        QTimer.singleShot(0, main_window, lambda: _tint_tool_icons(main_window))
    except (RuntimeError, TypeError):
        _LOGGER.warning(
            "[%s] Could not schedule tool icon tinting", PLUGIN_NAME, exc_info=True
        )

    # 3. Inform user
    if _CONTEXT:
        _CONTEXT.show_status_message(f"{PLUGIN_NAME} Active: Dark Mode Applied", 5000)
    elif hasattr(main_window, "statusBar"):
        main_window.statusBar().showMessage(
            f"{PLUGIN_NAME} Active: Dark Mode Applied", 5000
        )
