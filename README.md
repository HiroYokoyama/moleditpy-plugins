#  Official MoleditPy Plugin Collection

This directory contains the official sample plugins for **MoleditPy**. These examples are designed to demonstrate how to extend the application's functionality, interact with the API, and manipulate molecular data.

**Contribute Your Plugin**
We believe in the power of community\! If you have created a useful plugin, we would love to include it as an official part of this collection. Please feel free to submit a Pull Request with your plugin to help us expand what MoleditPy can do.

## Available Plugins

### 1. Hello World (`hello.py`)
A simple demonstration plugin.
-   **Features**:
    -   Shows a "Hello World" message box.
    -   Prints molecule information to the console if a molecule is loaded.

### 2. MS Spectrum Simulation (`ms_spectrum.py`)
Simulates the Mass Spectrum for the currently loaded molecule using RDKit descriptors.
-   **Features**:
    -   Displays Formula, Average Mass, and Exact Mass.
    -   Interactive histogram of isotopic distribution.
    -   **Export to Image**: Save the spectrum plot as a PNG/JPG file.
-   **Dependencies**: Requires `rdkit` installed in your Python environment.



## Installation

To install a plugin:

1.  Locate your MoleditPy plugins directory:
    -   **Windows**: `C:\Users\<YourUser>\.moleditpy\plugins`
    -   **macOS/Linux**: `~/.moleditpy/plugins`
    *(If the directory does not exist, create it manually or run MoleditPy once to generate it.)*

2.  Copy the desired `.py` file (e.g., `ms_spectrum.py`) from this repository into that `plugins` folder.

3.  Restart MoleditPy. The new plugin will appear in the "Plugins" menu.

## Development

You can create your own plugins by writing a Python script with a `run(main_window)` function.

```python
# example_plugin.py
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "My Custom Plugin"

def run(main_window):
    # Access the current molecule via main_window.current_mol
    if main_window.current_mol:
        QMessageBox.information(main_window, PLUGIN_NAME, "Molecule loaded!")
    else:
        QMessageBox.warning(main_window, PLUGIN_NAME, "No molecule.")
```


