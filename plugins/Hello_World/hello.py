# -*- coding: utf-8 -*-
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "Hello World"

def run(main_window):
    """
    Sample plugin that displays a message box.
    """
    mol = getattr(main_window, 'current_mol', None)
    atom_count = mol.GetNumAtoms() if mol else 0
    
    msg = f"Hello from the plugin!\nCurrent molecule has {atom_count} atoms."
    QMessageBox.information(main_window, "Plugin Message", msg)
