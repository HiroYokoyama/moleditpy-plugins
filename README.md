#  MoleditPy Plugin Collection

Repo: [https://github.com/HiroYokoyama/moleditpy-plugins/](https://github.com/HiroYokoyama/moleditpy-plugins/)  
Explorer: [https://hiroyokoyama.github.io/moleditpy-plugins/explorer/](https://hiroyokoyama.github.io/moleditpy-plugins/explorer/)


This directory contains the official plugins for **MoleditPy**. 

**Contribute Your Plugin**
We believe in the power of community\!  
If you have created a useful plugin, we would love to include it as an official part of this collection. Please feel free to submit a Pull Request with your plugin to help us expand what MoleditPy can do.

## New in v2.2: Plugin Manager

Manage your extensions easily with the new **Plugin Manager**:
- **Install/Uninstall**: Drag-and-drop `.py` files to install, or click "Remove" to uninstall.
- **Explore**: View plugin details, authors, and versions.
- **Online Library**: Direct link to the plugin explorer.
- **Access**: Open via **Plugin > Plugin Manager...**

## Featured Plugins

### 1. Gaussian Input Generator Neo (`gaussian_input_generator_neo.py`)
High-functionality setting dialog for Gaussian input creation.
- **Features**: Supports Link 0, Route, Title, Charge/Mult, and appended data.
- **Dependencies**: `rdkit`, `PyQt6`

### 2. ORCA Input Generator Neo (`orca_input_generator_neo.py`)
High-functionality setting dialog for ORCA input creation.
- **Features**: Supports Link 0, Route, Block Builder, and Validations.
- **Dependencies**: `rdkit`, `PyQt6`

### 3. Gaussian Freq Analyzer (`gaussian_fchk_freq_analyzer.py`)
Analyzes vibrational frequencies from Gaussian FCHK files.
- **Features**: View IR spectrum, animate normal modes, and export GIF animations.
- **Dependencies**: `rdkit`, `PyQt6`, `numpy`, `Pillow`

### 4. ORCA Freq Analyzer (`orca_out_freq_analyzer.py`)
Analyzes vibrational frequencies from ORCA output files.
- **Features**: View IR spectrum, animate normal modes, and export GIF animations.
- **Dependencies**: `rdkit`, `PyQt6`, `numpy`, `Pillow`

### 5. MS Spectrum Simulation Neo (`ms_spectrum_neo.py`)
Simulates the Mass Spectrum for the currently loaded molecule using RDKit descriptors.
- **Features**: Includes Gaussian broadening and interactive zoom/pan.
- **Dependencies**: `rdkit`, `PyQt6`

### 6. Mapped Cube Viewer (`mapped_cube_viewer.py`)
Visualizes a property (e.g. ESP) mapped onto an isosurface (e.g. electron density) from two Cube files.
- **Dependencies**: `numpy`, `pyvista`, `PyQt6`, `rdkit`

### 7. Cube File Viewer (`cube_viewer.py`)
Visualizes Gaussian Cube files (.cube) containing volumetric data (e.g., orbitals, densities).
- **Features**: Renders isosurfaces with interactive controls for isovalue and color.
- **Dependencies**: `rdkit`, `pyvista`, `numpy`

### 8. PubChem Name Resolver (`pubchem_ressolver.py`)
Resolves chemical names and identifiers to structures using the PubChem PUG REST API.
- **Features**: Search by Name or SMILES, load directly into 2D editor.
- **Dependencies**: `requests`, `rdkit`

### 9. Animated XYZ Giffer (`animated_xyz_giffer.py`)
A player for viewing multi-frame XYZ files and recording GIF animations using PIL.
- **Dependencies**: `Pillow`

### 10. Chat with Molecule Neo (Gemini) (`chat_with_molecule_neo.py`)

**The AI Copilot for MoleditPy (Powered by Google Gemini)**

Revolutionize your workflow with an AI agent that controls the editor directly. This plugin utilizes **Function Calling** to bridge natural language with chemical operations—loading structures, calculating descriptors, and generating expert-level ORCA inputs instantly.

* **Text-to-Structure:** Just say "Load Cubane" to generate 3D models.
* **Context-Aware:** Analyzes the active molecule in real-time.
* **Safe Execution:** Interactive [Accept]/[Reject] workflow for all AI actions.

[**Project Detail**](https://hiroyokoyama.github.io/moleditpy-plugins/plugins/Chat_with_Molecule_Neo/)

### Gallery

<p align="center">
  <img src="img/ms-spectrum.png" width="100%" alt="MS Spectrum Simulation">
</p>
<p align="center">
  <img src="img/gaussian-input.png" width="100%" alt="Gaussian Input Generator">
</p>
<p align="center">
  <img src="img/freq-analysis.png" width="100%" alt="Frequency Analyzer">
</p>
<p align="center">
  <img src="img/cube-viewer.png" width="100%" alt="Cube File Viewer">
</p>
<p align="center">
  <img src="img/esp-map.png" width="100%" alt="Cube File Viewer">
</p>

## Installation

To install a plugin:

1.  Locate your MoleditPy plugins directory:
    -   **Windows**: `C:\Users\<YourUser>\.moleditpy\plugins`
    -   **macOS/Linux**: `~/.moleditpy/plugins`  
    *(If the directory does not exist, create it manually or run "Open Plugin Directory" in MoleditPy application once to generate it.)*

2.  Copy the desired `.py` file (e.g., `ms_spectrum.py`) from this repository into that `plugins` folder.

3.  Restart MoleditPy. The new plugin will appear in the "Plugins" menu.

## Development

MoleditPy features a robust plugin architecture allowing deep integration via the **PluginContext** API.

### Quick Start

Create a `.py` file with an `initialize(context)` function:

```python
PLUGIN_NAME = "My New Plugin"
PLUGIN_VERSION = "1.0"
PLUGIN_AUTHOR = "Your Name"

def initialize(context):
    # Register a menu action
    context.add_menu_action("My Plugin/Say Hello", lambda: print("Hello!"))
```

### Capabilities

The new API allows plugins to:
*   Add **Menu** and **Toolbar** items.
*   Register **Bond/Atom Color Overrides**.
*   Handle **File Drops** and custom **File Imports**.
*   Add custom **Export Options** and **Optimization Methods**.
*   Integrate **Analysis Tools** and persist data in **Project Files**.

For full documentation and examples, please refer to the [PLUGIN_DEVELOPMENT_MANUAL.md](PLUGIN_DEVELOPMENT_MANUAL.md).
















