#  MoleditPy Plugin Collection

Repo: [https://github.com/HiroYokoyama/moleditpy-plugins/](https://github.com/HiroYokoyama/moleditpy-plugins/)  
Explorer: [https://hiroyokoyama.github.io/moleditpy-plugins/explorer/](https://hiroyokoyama.github.io/moleditpy-plugins/explorer/)


This directory contains the official plugins for **MoleditPy**. 

**Contribute Your Plugin**
We believe in the power of community\!  
If you have created a useful plugin, we would love to include it as an official part of this collection. Please feel free to submit a Pull Request with your plugin to help us expand what MoleditPy can do.

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


### 7. Animated XYZ Giffer (animated_xyz_giffer.py)
A player for viewing multi-frame XYZ files and recording GIF animations using PIL.  

- **Features:**
  - Loads concatenated XYZ files.
Playback controls: Play/Pause, Next/Previous frame, and Slider navigation.
  - Adjustable playback speed (FPS).
  - Exports playback as GIF animation.
  - Integrates with the main 3D viewer.
-   **Dependencies**:
Pillow
  (Pillow is required for GIF generation.)

### 8. PubChem Name Resolver (`pubchem_ressolver.py`)
Resolves chemical names and identifiers to structures using the PubChem PUG REST API.
- **Features**:
    - **Search**: Query PubChem by Name or SMILES.
    - **2D Loading**: Loads the resolved structure directly into the 2D editor.
    - **SMILES Fetching**: robustly retrieves Isomeric and Canonical SMILES.
    - **Table View**: Displays search results with Formula and Name.
- **Dependencies**: `requests`, `rdkit`.

### 9. Conformational Search (`conf_search.py`)
Performs conformational sampling and energy minimization.
- **Features**:
    - **Generation**: Uses RDKit's ETKDGv3 algorithm to generate 3D conformers.
    - **Optimization**: Optimizes geometries using **MMFF94** or **UFF** force fields.
    - **Interactive Preview**: Click on the result table to view conformers in the 3D viewer.
    - **Energy Ranking**: Sorts conformers by calculated energy (kcal/mol).
- **Dependencies**: `rdkit`.

### 10. Complex Molecule Untangler (`complex_molecule_untangler.py`)
A Monte Carlo-based tool to resolve steric clashes in complex or roughly drawn molecules.
- **Features**:
    - **Untangling**: Randomly rotates acyclic single bonds to reduce internal energy.
    - **Force Field Selection**: Supports **MMFF94** and **UFF** for energy evaluation.
    - **Real-time Progress**: Visualizes the untangling process.
- **Dependencies**: `rdkit`.

### 11. All-Trans Optimizer (`all-trans_optimizer.py`)
A geometric tool for straightening alkyl chains.
- **Features**:
    - **All-Trans**: Systematically sets dihedrals of acyclic carbon chains to 180 degrees.
    - **Structure Cleanup**: Useful for standardizing long alkyl chains.
- **Dependencies**: `rdkit`.

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

# Optional: Run automatically on startup/reload
def autorun(main_window):
    print(f"{PLUGIN_NAME} initialized locally.")
```

### Advanced Features

1.  **Subdirectories**: You can organize your plugins into subfolders within the `plugins` directory. MoleditPy will automatically create nested menus corresponding to the folder structure.
2.  **Autorun**: Define a function named `autorun(main_window)` in your plugin. This function will be executed immediately when the plugin is loaded (on application startup or when "Reload Plugins" is clicked). This is useful for plugins that need to register their own menu items, toolbars, or event listeners without waiting for the user to click the plugin name in the menu.

### Security Warning

> [!WARNING]
> **Plugins run arbitrary Python code with the same privileges as the MoleditPy application.**
> 
> *   Only install plugins from sources you trust.
> *   Be especially cautious with plugins that use `autorun`, as they execute code immediately upon loading without specific user action.
> *   Review the plugin code (`.py` file) if you are unsure about its functionality.












