# Gaussian Input Generator Neo

A [MoleditPy](https://github.com/HiroYokoyama/moleditpy) plugin for generating Gaussian quantum chemistry input files (`.gjf` / `.com`) with a full GUI, real-time preview, route builder, and preset management.

**Version:** 2026.04.13
**Author:** HiroYokoyama
**Category:** Export

---

## Features

### Route Builder
Interactive tabbed dialog for constructing the Gaussian route (`#P ...`) line:

| Tab | Options |
|-----|---------|
| **Method / Basis** | Output level (`#P` / `#` / `#T`), method type (DFT / HF), method name, basis set |
| **Job Type** | Opt+Freq, Opt, Freq, SP, Scan, IRC, Stable, Volume; convergence and frequency options |
| **Solvation / Dispersion** | PCM, CPCM, SMD, IEFPCM; 12 common solvents; GD3BJ dispersion correction |
| **Properties** | Population analysis (NBO, Hirshfeld, MK, Reg), TD-DFT, integration grid, symmetry, WFN output |

### Methods & Basis Sets
- **DFT:** WB97XD, B3LYP, M062X, PBE0, and more
- **HF / Semi-empirical / MP2** supported via manual route editing
- **Basis sets:** 6-31G(d), cc-pVTZ, def2-TZVP, Gen/GenECP (custom)

### Link 0 Resources
- `%mem` — memory with unit selector (GB / MB / MW)
- `%nprocshared` — CPU cores (1–128)
- `%chk` — checkpoint file name (auto-generated or custom)

### Charge / Multiplicity
- Spin boxes with live validation and visual feedback
- Auto-initialised from the loaded molecule

### Templates (post-coordinate input)
Insert boilerplate for common additional-input sections:

| Template | Use case |
|----------|----------|
| ModRedundant | Scan / freeze coordinates |
| Gen / GenECP | Custom basis sets and ECPs |
| ECP | Effective core potentials |
| NBO (`$NBO`) | NBO analysis keywords |
| Link1 | Multi-step jobs |
| Connectivity | `Geom=Connect` adjacency list |

### Preset Management
Save, load, and delete named presets for frequently used method/basis/job combinations.

### Real-time Preview
Editable preview pane — shows the full input file as you change settings. Manual edits are preserved until the next auto-refresh.

### Ghost Atom Support (`Bq` / `H:`)
Atoms carrying the `custom_symbol` property (set by **NICS Placer** or **XYZ Editor**) are rendered with that label in the coordinate block — no extra configuration needed:

```
Bq     0.000000     0.000000     1.000000
H:     0.000000     0.000000    -1.000000
```

---

## Workflow

1. Load a molecule with 3D coordinates in MoleditPy.
2. Open **Export → Gaussian Input…**
3. Configure Link 0 resources, route keywords, charge/multiplicity, and title.
4. Optionally open the **Route Builder** for guided keyword selection.
5. Add any post-coordinate input via the template selector.
6. Click **Save** to write the `.gjf` file.

---

## Integration with NICS Placer

The **NICS Placer** plugin places dummy atoms with `custom_symbol = "Bq"` or `"H:"`. This plugin reads those labels directly, so NICS input files are generated automatically without any manual editing.

---

## Requirements

- MoleditPy ≥ v3 (V3 plugin API)
- PyQt6
- RDKit
