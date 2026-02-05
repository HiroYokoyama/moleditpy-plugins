
import os
from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QSpinBox, QPushButton, QFileDialog, QMessageBox, QApplication
import pyvista as pv

PLUGIN_NAME = "High Resolution Imager"
PLUGIN_VERSION = "2026.02.05"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Take high-resolution screenshots of the 3D view with custom resolution."

def initialize(context):
    #context.add_menu_action("View/High Resolution Screenshot", lambda: take_screenshot(context))
    context.add_export_action("High Resolution Screenshot...", lambda: take_screenshot(context))

def take_screenshot(context):
    mw = context.get_main_window()
    
    # Dialog to ask for resolution
    dialog = QDialog(mw)
    dialog.setWindowTitle("High Resolution Screenshot")
    layout = QVBoxLayout()
    
    layout.addWidget(QLabel("Note: Resolution settings apply to raster images only (PNG, JPG, BMP, TIFF)"))
    
    # Width
    row1 = QHBoxLayout()
    row1.addWidget(QLabel("Width:"))
    width_spin = QSpinBox()
    width_spin.setRange(100, 20000) # Support very high res
    width_spin.setValue(1920) # Default 2K
    row1.addWidget(width_spin)
    layout.addLayout(row1)
    
    # Height
    row2 = QHBoxLayout()
    row2.addWidget(QLabel("Height:"))
    height_spin = QSpinBox()
    height_spin.setRange(100, 20000)
    height_spin.setValue(1080) # Default 2K
    row2.addWidget(height_spin)
    layout.addLayout(row2)
    
    # Buttons for Presets
    preset_layout = QHBoxLayout()
    preset_layout.addWidget(QLabel("Presets:"))
    btn_2k = QPushButton("2K (1080p)")
    btn_4k = QPushButton("4K (2160p)")
    btn_8k = QPushButton("8K (4320p)")
    preset_layout.addWidget(btn_2k)
    preset_layout.addWidget(btn_4k)
    preset_layout.addWidget(btn_8k)
    layout.addLayout(preset_layout)

    def set_res(w, h):
        width_spin.setValue(w)
        height_spin.setValue(h)

    btn_2k.clicked.connect(lambda: set_res(1920, 1080))
    btn_4k.clicked.connect(lambda: set_res(3840, 2160))
    btn_8k.clicked.connect(lambda: set_res(7680, 4320))

    # Transparent Background
    from PyQt6.QtWidgets import QCheckBox, QColorDialog
    from PyQt6.QtGui import QColor
    
    # Background Options using Radio Buttons for clarity
    from PyQt6.QtWidgets import QRadioButton, QButtonGroup
    
    bg_group = QButtonGroup(dialog)
    
    # 1. Transparent
    rb_trans = QRadioButton("Transparent Background (Raster only)")
    rb_trans.setChecked(True)
    bg_group.addButton(rb_trans)
    layout.addWidget(rb_trans)

    # 2. Current Background
    rb_current = QRadioButton("Current Background")
    bg_group.addButton(rb_current)
    layout.addWidget(rb_current)
    
    # 3. Custom Color
    custom_layout = QHBoxLayout()
    rb_custom = QRadioButton("Custom Color")
    bg_group.addButton(rb_custom)
    
    bg_btn = QPushButton("Select Color")
    bg_btn.setEnabled(False)
    bg_color_label = QLabel("   ")
    bg_color_label.setFixedSize(20, 20)
    bg_color_label.setStyleSheet("border: 1px solid black;")
    
    custom_layout.addWidget(rb_custom)
    custom_layout.addWidget(bg_btn)
    custom_layout.addWidget(bg_color_label)
    custom_layout.addStretch()
    layout.addLayout(custom_layout)

    # Store selected string color
    selected_bg_color = ["white"] 
    bg_color_label.setStyleSheet(f"background-color: white; border: 1px solid black;")

    def on_bg_mode_changed(button):
        bg_btn.setEnabled(button == rb_custom)

    bg_group.buttonClicked.connect(on_bg_mode_changed)

    def pick_bg_color():
        c = QColorDialog.getColor(QColor(selected_bg_color[0]), dialog, "Select Background Color")
        if c.isValid():
            selected_bg_color[0] = c.name()
            bg_color_label.setStyleSheet(f"background-color: {c.name()}; border: 1px solid black;")

    bg_btn.clicked.connect(pick_bg_color)

    # Buttons
    btn_layout = QHBoxLayout()
    save_btn = QPushButton("Save Screenshot")
    cancel_btn = QPushButton("Cancel")
    btn_layout.addWidget(save_btn)
    btn_layout.addWidget(cancel_btn)
    layout.addLayout(btn_layout)
    
    dialog.setLayout(layout)
    
    # Connect signals
    save_btn.clicked.connect(dialog.accept)
    cancel_btn.clicked.connect(dialog.reject)
    
    if dialog.exec() == QDialog.DialogCode.Accepted:
        width = width_spin.value()
        height = height_spin.value()
        
        # File Dialog
        filename, filter_selected = QFileDialog.getSaveFileName(mw, "Save Screenshot", "", "PNG Images (*.png);;JPEG Images (*.jpg *.jpeg);;BMP Images (*.bmp);;TIFF Images (*.tif *.tiff);;SVG Images (*.svg);;All Files (*)")
        if filename:
            try:
                # Handle SVG (vector) separately
                # Handle Background Override
                original_bg = None
                
                # Logic for background color determined by mode
                target_bg = None
                
                if rb_custom.isChecked():
                     target_bg = selected_bg_color[0]
                elif rb_trans.isChecked():
                     # For Raster, we handle via 'transparent_background' param.
                     # For Vector, we MUST set a solid color because GL2PS transparency is unreliable.
                     # We will default to WHITE for "Transparent" vector export as a safe fallback.
                     is_vector = filename.lower().endswith(('.svg', '.pdf', '.eps', '.ps'))
                     if is_vector:
                         target_bg = "white"
                         # Optionally warn user?
                         # QMessageBox.warning(mw, "Vector Transparency", "True transparency is not supported for vector formats.\nBackground has been set to White instead.")
                
                if target_bg:
                    try:
                        original_bg = mw.settings.get('background_color', '#4f4f4f')
                        mw.plotter.set_background(target_bg)
                    except:
                        original_bg = None
                
                try:
                    # Handle SVG (vector) separately
                    if filename.lower().endswith('.svg') or filename.lower().endswith('.pdf') or filename.lower().endswith('.eps') or filename.lower().endswith('.ps'):
                        mw.plotter.save_graphic(filename)
                    else:
                        # PyVista Screenshot with custom window_size
                        param_dict = {'window_size': (width, height)}
                        if rb_trans.isChecked():
                            param_dict['transparent_background'] = True
                        
                        mw.plotter.screenshot(filename, **param_dict)
                finally:
                    # Restore background
                    if original_bg:
                        mw.plotter.set_background(original_bg)

                QMessageBox.information(mw, "Success", f"Screenshot saved to:\n{filename}")
                
            except Exception as e:
                QMessageBox.critical(mw, "Error", f"Failed to take screenshot:\n{str(e)}\n\nNote: Extremely high resolutions might fail depending on graphics card memory.")
