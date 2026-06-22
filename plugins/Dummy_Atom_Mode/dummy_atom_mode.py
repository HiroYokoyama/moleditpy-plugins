# -*- coding: utf-8 -*-
import logging
from PyQt6.QtGui import QColor
from rdkit import Chem
from moleditpy.core.molecular_data import MolecularData
from moleditpy.utils.constants import CPK_COLORS, CPK_COLORS_PV
from moleditpy.ui.ui_manager import UIManager

# Metadata
PLUGIN_NAME = "Dummy Atom Mode"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Adds a dummy atom (*) mode to the 2D editor and registers a toolbar action to place dummy atoms."
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"

def initialize(context):
    # 1. Register the CPK color for dummy atoms (*)
    # Gray color makes the asterisk look distinct from standard element items
    gray_color = QColor("#808080")
    CPK_COLORS["*"] = gray_color
    CPK_COLORS_PV["*"] = [0.5, 0.5, 0.5]
    
    # 2. Monkeypatch Chem.Atom inside molecular_data.py to handle "*"
    # Chain of responsibility: we wrap the existing Chem.Atom method, preserving any other patches.
    # We also keep a reference to the original RDKit C++ constructor to prevent recursion.
    previous_Atom = Chem.Atom
    if not hasattr(Chem, "_original_Atom"):
        Chem._original_Atom = Chem.Atom
    
    def patched_Atom(symbol_or_atomic_num):
        if symbol_or_atomic_num == "*":
            return Chem._original_Atom(0)  # Atomic number 0 is the dummy atom (symbol '*') in RDKit
        return previous_Atom(symbol_or_atomic_num)
        
    Chem.Atom = patched_Atom

    # 3. Hook into UIManager.set_mode using Chain of Responsibility pattern
    # This synchronizes the checked state of the toolbar button with the scene mode.
    previous_set_mode = UIManager.set_mode
    
    def patched_set_mode(self, mode_str):
        if previous_set_mode is not None:
            previous_set_mode(self, mode_str)
        # Check if the current mode matches the dummy atom mode
        is_dummy = (mode_str == "atom_*")
        mw = self.host
        if hasattr(mw, "init_manager") and hasattr(mw.init_manager, "plugin_toolbar"):
            for action in mw.init_manager.plugin_toolbar.actions():
                if action.text() == "Dummy Atom Mode":
                    action.setChecked(is_dummy)
                    
    UIManager.set_mode = patched_set_mode

    # 4. Define the callback to toggle dummy atom mode
    def toggle_dummy_atom_mode():
        mw = context.get_main_window()
        current_mode = getattr(mw.init_manager.scene, "mode", None)
        if current_mode == "atom_*":
            # If already active, toggle back to selection mode
            mw.ui_manager.set_mode_and_update_toolbar("select")
        else:
            mw.ui_manager.set_mode_and_update_toolbar("atom_*")

    # 5. Add the action to the plugin toolbar
    context.add_toolbar_action(
        callback=toggle_dummy_atom_mode,
        text="Dummy Atom Mode",
        tooltip="Enter dummy atom mode to place *",
    )
    
    # 6. Make the newly created QAction checkable once the UI is asynchronously populated
    from PyQt6.QtCore import QTimer
    
    def make_action_checkable():
        mw = context.get_main_window()
        if hasattr(mw, "init_manager") and hasattr(mw.init_manager, "plugin_toolbar"):
            for action in mw.init_manager.plugin_toolbar.actions():
                if action.text() == "Dummy Atom Mode":
                    action.setCheckable(True)
                    current_mode = getattr(mw.init_manager.scene, "mode", None)
                    action.setChecked(current_mode == "atom_*")
                    
    QTimer.singleShot(100, make_action_checkable)

