# MoleditPy Plugin Development Manual

Welcome to the **MoleditPy Plugin Development Manual**! MoleditPy is designed to be extensible, allowing you to add new features, analysis tools, and custom visualizations with ease using Python.

## 1. Introduction

MoleditPy searches for plugins in your user directory:
- **Windows**: `C:\Users\<YourName>\.moleditpy\plugins\`
- **Linux/macOS**: `~/.moleditpy/plugins/`

Any `.py` file placed in this directory (or its subdirectories) is automatically discovered when MoleditPy starts.

## 2. Plugin Structure

A minimum viable plugin requires just a few lines of code.

### Required Attributes

- **`PLUGIN_NAME`**: A string variable defining how your plugin appears in the menu.
- **`run(main_window)`**: A function that gets called when the user clicks your plugin in the menu.

### Optional Attributes

- **`autorun(main_window)`**: A function called automatically when MoleditPy starts up. Useful for registering background listeners or adding persistent UI elements.

### Basic Template

```python
# my_first_plugin.py
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "My First Plugin"

def run(main_window):
    QMessageBox.information(main_window, "Hello", "Hello from My First Plugin!")

def autorun(main_window):
    print("My First Plugin loaded!")
```

## 3. API Reference

When `run()` or `autorun()` is called, you receive the `main_window` instance. This object gives you full control over the application.

### Key `main_window` Attributes

| Attribute | Type | Description |
| :--- | :--- | :--- |
| `main_window.current_mol` | `rdkit.Chem.Mol` | The currently loaded RDKit molecule object. Can be `None`. |
| `main_window.data` | `MolecularData` | Wrapper for molecular data management. |
| `main_window.plotter` | `pyvista.QtInteractor` | The 3D viewer widget (based on PyVista). Use this to add custom 3D actors. |
| `main_window.view_2d` | `QGraphicsView` | The 2D editor view. |
| `main_window.scene` | `MoleculeScene` | The 2D scene containing atom/bond items. |

### Expanded Method Reference

#### 1. Data Manipulation (`main_window.data`)
Access via `main_window.data`.
- **`add_atom(symbol, pos, charge=0, radical=0)`**: Adds an atom at the given 2D position `(x, y)`. Returns an internal ID.
- **`add_bond(id1, id2, order=1, stereo=0)`**: Adds a bond between two atoms by ID.
- **`remove_atom(atom_id)`**: Removes an atom and its connected bonds.
- **`to_rdkit_mol()`**: Converts the current 2D sketch to an RDKit Mol object (2D coordinates).

#### 2. 2D Scene Control (`main_window.scene`)
Access via `main_window.scene`.
- **`create_atom(symbol, pos)`**: UI-aware method to add an atom item to the scene.
- **`create_bond(atom1_item, atom2_item, order)`**: UI-aware method to add a bond item.
- **`delete_items(items_list)`**: Safely deletes a list of `AtomItem` or `BondItem` objects.
- **`select_all()`**: Selects all items in the 2D scene.
- **`clear_all()`**: Clears the entire scene.

#### 3. 3D View Control
Methods directly on `main_window` (delegated to `MainWindowView3d`):
- **`draw_molecule_3d(mol)`**: Updates the 3D view with the provided RDKit molecule.
- **`set_3d_style(style_name)`**: Changes 3D rendering style (`'ball_and_stick'`, `'cpk'`, `'wireframe'`, `'stick'`).
- **`toggle_atom_info_display(mode)`**: Toggles labels on atoms (`'id'`, `'symbol'`, `'charge'`, etc.).
- **`export_3d_png()`**: Opens a dialog to save the current 3D view as an image.

#### 4. Computation & Optimization
Methods delegated to `MainWindowCompute`:
- **`trigger_conversion()`**: Converts 2D -> 3D.
- **`optimize_3d_structure()`**: Optimizes the geometry using force fields.
- **`set_optimization_method(method)`**: Sets the force field (`'MMFF_RDKIT'`, `'GAFF'`, etc.).
- **`check_chemistry_problems_fallback()`**: Runs basic valence checks on the molecule.

#### 5. Edit Actions
Methods delegated to `MainWindowEditActions`:
- **`copy_selection()`**: Copies selected 2D items to the system clipboard.
- **`paste_from_clipboard()`**: Pastes 2D structure from the clipboard.
- **`clean_up_2d_structure()`**: Auto-arranges the 2D layout.
- **`add_hydrogen_atoms()`**: Adds explicit H atoms to the 2D drawing.
- **`remove_hydrogen_atoms()`**: Removes explicit H atoms from the 2D drawing.

#### 6. Dialogs & UI Helpers
Methods delegated to `MainWindowDialogManager`:
- **`show_about_dialog()`**: Shows the About dialog.
- **`open_periodic_table_dialog()`**: Opens the periodic table chooser.
- **`open_settings_dialog()`**: Opens the main settings window.
- **`open_analysis_window()`**: Opens the analysis tools window.

#### 7. Project I/O
Methods delegated to `MainWindowProjectIo`:
- **`save_project()`**: Trigger the standard "Save" flow.
- **`save_as_json()`**: Return the full project state as a JSON-compatible dict.
- **`load_json_data(data)`**: Restore project state from a dict.

#### 8. App State & Undo
Methods delegated to `MainWindowAppState`:
- **`undo()`**: Perform undo.
- **`redo()`**: Perform redo.
- **`push_undo_state()`**: Manually push the current state to the undo stack (useful before making programmatic changes).

#### 9. Helper Constants
You can import useful constants from `modules.constants`.
```python
from modules.constants import CPK_COLORS, VDW_RADII, pt

# CPK_COLORS: Dictionary of symbol -> QColor
# VDW_RADII: Dictionary of symbol -> float (VDW radius)
# pt: RDKit PeriodicTable instance
```

#### 10. Advanced 3D Interaction
The `main_window.plotter` object is a `pyvistaqt.QtInteractor`. You can use standard PyVista methods to manipulate the scene.
- **`plotter.add_mesh(mesh, ...)`**: Add PyVista meshes (PolyData, UnstructuredGrid).
- **`plotter.add_actor(actor)`**: Add VTK actors directly.
- **`plotter.camera`**: Access the camera object.
- **`plotter.clear()`**: Clear the scene (be careful as this removes the molecule too!).

#### 11. Available Libraries
Your plugin runs in the same environment as MoleditPy. You can freely import:
- **`rdkit`, `rdkit.Chem`**: Full RDKit functionality.
- **`pyvista`**: For complex 3D geometry generation.
- **`PyQt6`**: For creating custom UI elements (widgets, dialogs, menus).

## 4. Examples

### Example 1: Analyze Molecule

This plugin calculates the molecular weight of the current molecule and shows it in a dialog.

```python
from PyQt6.QtWidgets import QMessageBox
from rdkit.Chem import Descriptors

PLUGIN_NAME = "Calculate Stats"

def run(main_window):
    mol = main_window.current_mol
    if not mol:
        QMessageBox.warning(main_window, "Error", "No molecule loaded!")
        return

    mw = Descriptors.MolWt(mol)
    num_atoms = mol.GetNumAtoms()
    
    QMessageBox.information(
        main_window, 
        "Molecule Stats", 
        f"Molecular Weight: {mw:.2f} g/mol\nAtom Count: {num_atoms}"
    )
```

### Example 2: Custom 3D Object

This plugin adds a transparent sphere to the 3D scene using PyVista.

```python
import pyvista as pv

PLUGIN_NAME = "Add 3D Sphere"

def run(main_window):
    plotter = main_window.plotter
    
    # Create a sphere
    sphere = pv.Sphere(radius=2.0)
    
    # Add it to the plotter
    plotter.add_mesh(sphere, color="red", opacity=0.3, name="custom_sphere")
    
    # Refresh view
    plotter.render()
```

### Example 3: Adding a Custom Dock Panel

This plugin adds a permanent side panel with a custom button.

```python
from PyQt6.QtWidgets import QDockWidget, QPushButton, QVBoxLayout, QWidget, QLabel
from PyQt6.QtCore import Qt

PLUGIN_NAME = "Custom Panel"

def autorun(main_window):
    # check if already added to avoid duplicates on reload (if applicable)
    if hasattr(main_window, 'my_custom_dock'):
        return

    dock = QDockWidget("My Tools", main_window)
    main_window.my_custom_dock = dock # Keep reference
    
    content = QWidget()
    layout = QVBoxLayout(content)
    
    lbl = QLabel("Custom Tool Panel")
    btn = QPushButton("Click Me")
    btn.clicked.connect(lambda: print("Button Clicked!"))
    
    layout.addWidget(lbl)
    layout.addWidget(btn)
    layout.addStretch()
    
    dock.setWidget(content)
    main_window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)

def run(main_window):
    # If users click the menu item, just show the dock if it was hidden
    if hasattr(main_window, 'my_custom_dock'):
        main_window.my_custom_dock.show()
```

### Example 4: Custom File Loader

This example demonstrates how to parse a custom text format and build a molecule programmatically using `main_window.data`.

```python
# Format: "Symbol x y" per line
# Example:
# C 0.0 0.0
# O 1.5 0.0

from PyQt6.QtWidgets import QFileDialog, QMessageBox

PLUGIN_NAME = "Import Simple XY"

def run(main_window):
    file_path, _ = QFileDialog.getOpenFileName(
        main_window, "Open Simple XY", "", "Text Files (*.txt)"
    )
    
    if not file_path:
        return

    try:
        # 1. Clear existing 2D scene
        main_window.clear_all()
        
        # 2. Parse file and build data
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    symbol = parts[0]
                    x = float(parts[1]) * 50 # Scale up for pixels
                    y = float(parts[2]) * 50
                    
                    # Add atom to data model
                    main_window.data.add_atom(symbol, (x, y))

        # 3. Refresh the 2D scene from data
        main_window.scene.reinitialize_items()
        
        # 4. Push to undo stack
        main_window.push_undo_state()
        
        main_window.statusBar().showMessage(f"Imported {file_path}")

    except Exception as e:
        QMessageBox.critical(main_window, "Import Error", str(e))
```

## 5. Development Tips

1.  **Dependencies**: MoleditPy uses `PyQt6` for UI, `rdkit` for chemistry, and `pyvista` for 3D rendering. You can import these directly in your plugins.
2.  **Console Output**: Use `print()` for debugging. Output appears in the terminal where you launched MoleditPy.
3.  **Error Handling**: Wrap your code in `try...except` blocks to prevent crashing the main application.
4.  **Reloading**: Currently, you must restart MoleditPy to reload plugin code changes.

Happy Coding!
