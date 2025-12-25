# MoleditPy Plugin Development Manual

Welcome to the **MoleditPy Plugin Development Manual**! MoleditPy features a robust plugin system that allows you to extend the application with new tools, menus, visualizations, and file formats.

## 1. Introduction

MoleditPy plugins are Python scripts (`.py`) placed in your user plugin directory:
- **Windows**: `C:\Users\<YourName>\.moleditpy\plugins\`
- **Linux/macOS**: `~/.moleditpy/plugins/`

The application automatically discovers and loads these plugins on startup.

## 2. Plugin Structure

A modern MoleditPy plugin consists of **Metadata** and an **Initialization Hook**.

### Metadata Variables
Define these at the top level of your script:
```python
PLUGIN_NAME = "My Super Plugin"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "Jane Doe"
PLUGIN_DESCRIPTION = "Adds super powers to MoleditPy."
```

### The `initialize` Function
This is the entry point. It receives a `context` object (of type `PluginContext`) that you use to register your extensions.

```python
def initialize(context):
    # Register your hooks here
    context.add_menu_action("My Plugin/Action", my_callback)
```

---

## 3. The PluginContext API

The `context` object (Type: `PluginContext`) passed to your `initialize` function provides the following methods:

### 3.1 UI & Menus

| Method | Description |
| :--- | :--- |
| `add_menu_action(path, callback, text=None, icon=None, shortcut=None)` | Add an item to the main menu. `path` can be "MyPlugin/Action" or standard paths like "File/MyAction". |
| `add_toolbar_action(callback, text, icon=None, tooltip=None)` | Add a button to the dedicated **Plugin Toolbar**. Buttons auto-hide if no plugins use them. |
| `add_analysis_tool(label, callback)` | Add an entry to the top-level **Analysis** menu. |
| `add_export_action(label, callback)` | Add an entry to the **Export** menu button (and File > Export). |

### 3.2 File Handling & Project State

| Method | Description |
| :--- | :--- |
| `register_file_opener(extension, callback)` | Handle opening files with a specific extension (e.g., `.xyz`). |
| `register_drop_handler(callback, priority=0)` | Handle files dropped onto the main window. `callback(path)` should return `True` if handled. |
| `register_save_handler(callback)` | Register a function that returns a `dict` of data to be saved in the project file (`.pmeprj`). |
| `register_load_handler(callback)` | Register a function to receive the saved `dict` when specific project data is loaded. |

### 3.3 Computation & 3D

| Method | Description |
| :--- | :--- |
| `register_3d_style(style_name, callback)` | **New in v2.2** - Register a custom 3D visualization mode. Appears in the "3D Style" menu. |
| `register_3d_context_menu(label, callback)` | Add an action to the right-click context menu in the 3D view. |
| `register_optimization_method(name, callback)` | Register a new method for the **Compute > Optimize Geometry** menu. |
| `get_3d_controller()` | Returns a `Plugin3DController` instance for manipulating the 3D view. |
| `get_main_window()` | Returns the raw `MainWindow` instance (Use with caution). |
| `add_panel_button(text, callback, panel="right")` | Add a button to the bottom control panel ("left" or "right"). |
| `register_save_handler(callback)` | Register a function calculating data to save in `.pmeprj`. |
| `register_load_handler(callback)` | Register a function to restore state from loaded `.pmeprj` data. |

### 3.4 3D Controller API (`get_3d_controller()`)

| Method | Description |
| :--- | :--- |
| `set_atom_color(atom_index, color_hex)` | Override the color of a specific atom (e.g., `"#FF0000"`). Set to `None` to reset. |
| `set_bond_color(bond_index, color_hex)` | Override the color of a specific bond. Set to `None` to reset. |

### 3.5 Custom 3D Style API (`register_3d_style`)

Register a callback that fully handles drawing the molecule in the 3D view.

```python
def draw_my_style(main_window, mol):
    # Use mw.plotter (PyVista) to draw
    main_window.plotter.add_mesh(...)

def initialize(context):
    context.register_3d_style("My Custom Style", draw_my_style)
```

### 3.5 Power User Access (`get_main_window()`)

While the `PluginContext` API covers common use cases, you are not limited by it.

**You can change everything via `mw`**:
Code obtained via simple `context.get_main_window()` has **unrestricted access** to the entire application instance (`MainWindow`).
*   Modify any Qt widget directly.
*   Monkey-patch internal methods.
*   Access the full RDKit molecule or PyVista plotter.

```python
mw = context.get_main_window()
mw.setWindowTitle("My Custom Title") # Change the window title
mw.menuBar().clear() # Remove all menus (if you really want to!)
```


## 4. Legacy Support (Pre-2.2)

Older plugins used a different structure which is still supported but less capable.

### The `run` Function
If you define `run(main_window)`, your plugin will automatically appear in the "Plugins" menu. When clicked, `run` is executed with the raw `MainWindow` object.

```python
# Legacy Style
def run(main_window):
    mol = main_window.current_mol
    # ... direct access to implementation details ...
```

### The `autorun` Function
Legacy plugins often used `autorun(main_window)` to execute code immediately on startup. Use `initialize(context)` for this purpose in modern plugins.

```python
def autorun(main_window):
    print("I run immediately!")
```

**Note**: You can mix both! Define `initialize(context)` for startup hooks and `run(main_window)` to show a manual action in the main "Plugins" menu.

---

## 5. Full Example

```python
import os
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "Complete Example"
PLUGIN_VERSION = "2.0"
PLUGIN_AUTHOR = "MoleditPy Team"

def initialize(context):
    print("Initializing Example Plugin...")

    # 1. Add a menu item
    def show_info():
        # You can access the main window via context if needed
        mw = context.get_main_window()
        QMessageBox.information(mw, "Info", "Plugin is active!")
    
    context.add_menu_action("Example Plugin/Show Info", show_info)

    # 2. Add an analysis tool
    context.add_analysis_tool("Count Atoms", lambda: print("Counting..."))

    # 3. Handle a custom file type
    context.register_file_opener(".ex", lambda path: print(f"Opening {path}"))

# Optional: Add a 'Run' entry to the Plugins menu
def run(mw):
    QMessageBox.information(mw, "Manual Run", "You clicked me in the Plugins menu!")
```

## 6. Cookbook / Examples

Here are common example references, updated for the modern API.
Note: For advanced logic involving the 2D scene or 3D plotter, you will often use `context.get_main_window()` (`mw`) to access `mw.scene`, `mw.data`, or `mw.plotter`.

### Example 1: Analysis Tool (Molecular Weight)
Add a tool to the **Analysis** menu that calculates properties.

```python
from rdkit.Chem import Descriptors
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "Calculate Stats"

def initialize(context):
    def run_calc():
        # Access the main application to get the current molecule
        mw = context.get_main_window()
        mol = mw.current_mol
        
        if not mol:
            QMessageBox.warning(mw, "Error", "No molecule loaded!")
            return

        weight = Descriptors.MolWt(mol)
        num_atoms = mol.GetNumAtoms()
        
        QMessageBox.information(mw, "Stats", f"Molecular Weight: {weight:.2f}\nAtoms: {num_atoms}")

    context.add_analysis_tool("Show Molecular Weight", run_calc)
```

### Example 2: Custom 3D Visualization (Add Sphere)
Use the `plotter` object (PyVista) to add custom 3D geometries.

```python
import pyvista as pv

PLUGIN_NAME = "Add 3D Sphere"

def initialize(context):
    context.add_menu_action("Visuals/Add Sphere", lambda: add_sphere(context))

def add_sphere(context):
    mw = context.get_main_window()
    plotter = mw.plotter
    
    # Create and add a sphere
    sphere = pv.Sphere(radius=2.0)
    plotter.add_mesh(sphere, color="red", opacity=0.3, name="custom_sphere")
    plotter.render()
```

### Example 3: Custom Dock Panel
Add a permanent side panel UI using standard Qt widgets.

```python
from PyQt6.QtWidgets import QDockWidget, QLabel, QVBoxLayout, QWidget, QPushButton
from PyQt6.QtCore import Qt

PLUGIN_NAME = "My Panel"

def initialize(context):
    mw = context.get_main_window()
    
    # Prvenet duplicates on plugin reload
    if hasattr(mw, 'my_custom_dock'): return

    dock = QDockWidget("My Tools", mw)
    content = QWidget()
    layout = QVBoxLayout(content)
    layout.addWidget(QLabel("Hello from Plugin!"))
    layout.addWidget(QPushButton("Click Me"))
    dock.setWidget(content)
    
    mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
    
    # Store reference so we don't recreate it
    mw.my_custom_dock = dock
```

### Example 4: Custom File Importer
Register a handler for a custom file format (e.g. `.simple`).

```python
PLUGIN_NAME = "Simple Importer"

def initialize(context):
    context.register_file_opener(".simple", lambda path: load_simple(path, context))

def load_simple(path, context):
    mw = context.get_main_window()
    print(f"Parsing {path}...")
    
    # Example: Clear and add a dummy atom
    mw.clear_all()
    mw.data.add_atom("C", (0, 0))
    mw.scene.reinitialize_items()
```

### Example 5: Custom 3D Style (Native Registration)
Register a new visualization mode that appears in the 3D Style menu.

```python
import pyvista as pv

PLUGIN_NAME = "My Style Plugin"

def draw_custom_style(mw, mol):
    # 1. Base drawing (optional helper to draw standard atoms)
    # mw.main_window_view_3d.draw_standard_3d_style(mol, style_override='ball_and_stick')
    
    # 2. Add custom visualization
    mw.plotter.add_text("Custom Style Active", position='upper_left')
    mw.plotter.add_mesh(pv.Sphere(radius=5), color='blue', opacity=0.2)

def initialize(context):
    context.register_3d_style("My Blue Sphere", draw_custom_style)
```
