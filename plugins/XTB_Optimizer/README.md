# xTB Optimizer

A MoleditPy plugin for semiempirical geometry optimization using the **xTB** method family (GFN2-xTB, GFN1-xTB) via the [`tblite`](https://github.com/tblite/tblite) package and the [ASE](https://wiki.fysik.dtu.dk/ase/) LBFGS optimizer.

Optimization runs on a **background thread** — the UI stays fully responsive during calculation.

---

## Features

- GFN2-xTB and GFN1-xTB methods
- Live per-step energy and force table
- Scrolling log output
- Cancel button (stops the optimizer mid-run cleanly)
- Undo-compatible (coordinates committed via `push_undo_checkpoint`)

## Installation

### 1. Install dependencies

> **Windows: use conda-forge, not pip.**
> `tblite` has no pre-built pip wheel for Windows — `pip install tblite` will
> attempt to compile from source and fail unless LAPACK is present.

```bash
mamba install -c conda-forge tblite-python ase
# or
conda install -c conda-forge tblite-python ase
```

### 2. Install the plugin

Copy the `XTB_Optimizer/` folder (or just `xtb_optimizer.py`) to your MoleditPy plugin directory:

- **Windows**: `C:\Users\<YourName>\.moleditpy\plugins\`
- **Linux / macOS**: `~/.moleditpy/plugins/`

Then restart MoleditPy or use **Plugins › Reload All Plugins**.

## Usage

1. Load a molecule and generate 3D coordinates.
2. Open **3D Edit › xTB Optimizer…**
3. Choose method, force threshold, and max steps.
4. Click **▶ Run**. Monitor progress in the table and log.
5. Close the dialog when done — optimized coordinates are already applied.

## Settings

| Setting | Default | Description |
|---|---|---|
| Method | GFN2-xTB | Hamiltonian: GFN2-xTB (more accurate) or GFN1-xTB (faster) |
| Force threshold | 0.05 eV/Å | LBFGS convergence criterion |
| Max steps | 500 | Maximum optimizer iterations |

## Notes

- Dummy atoms (`*`) are not supported and will be rejected before the run starts.
- Molecules without 3D coordinates must be converted first (**3D Edit › Generate 3D Coordinates**).
- The plugin uses the ASE **LBFGS** optimizer internally.
