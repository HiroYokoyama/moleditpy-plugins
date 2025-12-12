#  Official MoleditPy Plugin Collection

Repo: [https://github.com/HiroYokoyama/moleditpy-plugins/](https://github.com/HiroYokoyama/moleditpy-plugins/)

This directory contains the official plugins for **MoleditPy**. 

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




### 3. Gaussian Input Generator (`gaussian_input_generator.py`)
A comprehensive tool to generate input files for Gaussian quantum chemistry calculations.
-   **Features**:
    -   Configuration of Link 0 commands (Memory, Processors, Checkpoint).
    -   Customizable Route Section, Title, Charge, and Multiplicity.
    -   Support for additional input sections (e.g., ModRedundant, Basis Sets).
    -   Saves files with `.gjf` or `.com` extensions.

### 4. ORCA Input Generator (`orca_input_generator.py`)
Generates input files for the ORCA quantum chemistry package.
-   **Features**:
    -   Setup for Parallelization (`%pal`) and Memory per core (`%maxcore`).
    -   Simple keyword input (starts with `!`).
    -   Support for Advanced Blocks (e.g., `%scf`, `%basis`) before coordinates.
    -   Saves as `.inp` files.

### 5. Animated XYZ Player (`animated_xyz.py`)
A player for viewing multi-frame XYZ files, such as molecular dynamics trajectories.
-   **Features**:
    -   Loads concatenated XYZ files.
    -   Playback controls: Play/Pause, Next/Previous frame, and Slider navigation.
    -   Adjustable playback speed (FPS).
    -   Integrates with the main 3D viewer.

### 6. Cube File Viewer (`cube_viewer.py`)
Visualizes Gaussian Cube files (.cube) containing volumetric data (e.g., orbitals, densities).
-   **Features**:
    -   Renders Isosurfaces for positive and negative values.
    -   Interactive controls for Isovalue, Color, and Opacity.
    -   Option to use complementary colors for negative lobes.
    -   Automatically generates molecular structure from the cube file headers.
-   **Dependencies**: Requires `rdkit`, `pyvista`, and `numpy`.


## Installation

To install a plugin:

1.  Locate your MoleditPy plugins directory:
    -   **Windows**: `C:\Users\<YourUser>\.moleditpy\plugins`
    -   **macOS/Linux**: `~/.moleditpy/plugins`  
    *(If the directory does not exist, create it manually or run "Open Plugin Directory" in MoleditPy application once to generate it.)*

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







