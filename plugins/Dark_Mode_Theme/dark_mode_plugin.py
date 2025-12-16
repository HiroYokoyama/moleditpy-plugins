"""
dark_mode_plugin.py

A plugin that automatically applies a dark mode stylesheet to MoleditPy upon loading.
"""

from PyQt6.QtWidgets import QApplication, QMainWindow
from PyQt6.QtGui import QColor

__version__="2025.12.16"
__author__="HiroYokoyama"

PLUGIN_NAME = "Dark Mode Theme"

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

def autorun(main_window):
    """
    Executed automatically when the plugin is loaded (and 'autorun' is detected).
    """
    print(f"[{PLUGIN_NAME}] Applying Dark Mode Stylesheet...")
    
    # 1. Apply QSS to the main window (and thus all children)
    # Note: Applying to main_window is safer than QApplication for plugins to avoid
    # affecting other parts if this were embedded in a larger system, but for a 
    # value-add plugin for this specific app, it works well.
    main_window.setStyleSheet(DARK_STYLESHEET)
    
    # 2. Update 3D Background Color to match dark theme
    # The default light background (#919191 or #FFFFFF) might be too bright.
    try:
        new_bg_color = "#2b2b2b" # Dark grey
        
        # Check if settings exist and update
        if hasattr(main_window, 'settings'):
            main_window.settings['background_color'] = new_bg_color
            
            # Since we modify the settings directly, we should try to trigger an update
            # The 'apply_3d_settings' method in main_window usually handles this.
            if hasattr(main_window, 'apply_3d_settings') and hasattr(main_window, 'plotter'):
                main_window.apply_3d_settings()
            
            # Also update icon foregrounds if possible (custom logic in main_window)
            # Setting 'icon_foreground' to white makes sure icons are visible on dark bg
            main_window.settings['icon_foreground'] = '#FFFFFF'
            
            # Helper: Force update icons if main_window has the method (it's in init_ui usually)
            # Re-calling init_menu_bar might be too destructive, but we can try to refresh toolbars
            # if there is a specific method. For now, rely on restart or dynamic icon color logic
            # if it observes settings changes.
            
    except Exception as e:
        print(f"[{PLUGIN_NAME}] Warning: Could not auto-set 3D background color: {e}")

    # 3. Inform user
    if hasattr(main_window, 'statusBar'):
        main_window.statusBar().showMessage(f"{PLUGIN_NAME} Active: Dark Mode Applied", 5000)
