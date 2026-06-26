# Python Console

An embedded interactive Python console for MoleditPy.

## Overview

Opens a live Python interpreter inside MoleditPy. Run arbitrary scripts against the current molecule, call PluginContext API methods, or explore the application state interactively.

## Available Variables

| Variable | Type | Description |
|---|---|---|
| `mol` | `rdkit.Chem.Mol` or `None` | Current molecule (refreshed on every execution) |
| `mw` | `MainWindow` | Host application main window |
| `context` | `PluginContext` | Plugin API — preferred over `mw` for stable access |
| `Chem` | module | `rdkit.Chem` |
| `help` | `GUIHelp` | Safe help utility (non-freezing replacement for built-in `help`) |

## Usage

- **Enter** — execute the current input
- **Shift+Enter** — insert a newline (for multi-line blocks)
- **Up / Down** — navigate command history

## Examples

```python
# Atom count
mol.GetNumAtoms()

# Show a status message via PluginContext
context.show_status_message("Hello from console!", 3000)

# Get current molecule via context
m = context.current_molecule
print(m.GetNumAtoms() if m else "No molecule loaded")

# Set a new molecule
from rdkit.Chem import AllChem
m2 = Chem.MolFromSmiles("CCO")
AllChem.EmbedMolecule(m2)
context.current_molecule = m2
```

## Redrawing the 3D View

Use `context.draw_molecule_3d()` to force a redraw after modifying a molecule:

```python
# Modify the molecule in-place, then redraw
from rdkit.Chem import AllChem
mol = context.current_molecule
AllChem.MMFFOptimizeMolecule(mol)
context.draw_molecule_3d(mol)

# Or assign via the setter — sets mol AND redraws in one step
context.current_molecule = mol
```

Assigning to `context.current_molecule` is equivalent to writing the molecule and calling `draw_molecule_3d()` — use it when you want to replace the active molecule. Call `context.draw_molecule_3d(mol)` explicitly when you have mutated the existing molecule in-place and only need to refresh the view.

## Notes

- `mol` is a snapshot taken at execution time; use `context.current_molecule` inside loops or callbacks for the live value.
- `mw` exposes internal app state — prefer `context` methods where an equivalent exists.
- Output and errors appear in the output panel above the input field.
