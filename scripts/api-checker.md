# MoleditPy Static API Checker

This directory contains the tools used to maintain API compatibility between the MoleditPy core application and its extensive plugin ecosystem. These tools use **Static Analysis** (via the Python `ast` module) to detect broken attribute accesses without actually executing the code.

## Core Scripts

### `plugin_api_checker.py` (The Engine)
The fundamental scanning logic. It implements a two-pass visitor pattern to map the application's API and then verify plugin usage against that map.

### `check_api.py`
A repository-specific variant of the checker script. It contains the same core engine but is pre-configured for the `moleditpy-plugins` repository with:
- Pre-defined paths for the main app and plugins.
- Support for registry-based scanning (via `REGISTRY/plugins.json`).
- A customized **Allowlist** for confirmed-valid dynamic runtime attributes.

---

## Technical Architecture

### Phase 1: API Surface Extraction
Before checking any plugins, the tool builds a "truth map" of the application API:

1.  **MainWindow Analysis**: It parses `main_window.py` to index all methods, class attributes, and signals.
2.  **Manager Discovery**: It identifies attributes on the `MainWindow` that are assigned instances of known Manager classes (e.g., `View3DManager`, `StateManager`).
3.  **Cross-File Parsing**: For every discovered Manager, it locates the source file, parses the class definition, and indexes:
    *   Public methods and `@property` getters.
    *   Attributes explicitly assigned to `self` in `__init__` and setup methods.
4.  **Host-Delegation Support**: It scans every `.py` file in the core app for the `self.host.xxx = ...` pattern. This handles cases where a Manager dynamically adds an attribute to the `MainWindow` (Host).

### Phase 2: Plugin Verification
The tool then scans plugin files in two passes:

1.  **Pass 1 (Alias Tracking)**: It identifies which variables hold references to the application. It recognizes:
    *   Standard names: `mw`, `main_window`, `context`, `ctx`.
    *   Method calls: `context.get_main_window()`.
    *   Propagated aliases: `self.mw = mw`.
2.  **Pass 2 (Access Checking)**: It examines every attribute access (e.g., `mw.view_3d_manager.current_mol`).
    *   **Direct Access**: `mw.some_method` must exist on `MainWindow`.
    *   **Chained Access**: `mw.manager.attr` must verify that `manager` exists on `MainWindow` and `attr` exists on the corresponding Manager class.

---

## Handling False Positives

Static analysis cannot see everything that happens at runtime. The tool uses several heuristics to reduce noise:

-   **Inherited Methods**: It has a built-in list (`_QT_INHERITED`) of standard `QWidget`/`QMainWindow` methods (like `setVisible`, `show`, `resize`) that are always permitted.
-   **Try-Block Heuristic**: Accesses inside a `try: ... except:` block are flagged with a `[try]` tag and are often ignored in CI, as they imply the developer has provided a fallback.
-   **The Manager Allowlist** (`--default-allowlist`): Suppresses confirmed-valid dynamic attributes on manager objects — those assigned via `self.host.init_manager.X = ...` patterns invisible to AST (e.g. `init_manager.scene`, `state_manager.data`). Safe to enable in all scans.
-   **The MW Allowlist** (`--mw-allowlist`): Suppresses known legacy compat bridge attrs accessed directly on the MainWindow (e.g. `mw.host`, `mw.view3d`). **Off by default** — enabling it will hide real V3 migration bugs where plugins still call methods that have moved to managers.

---

## Usage Guide

### Basic Scan (strictest — recommended for CI)
Check a specific plugin against the local app repository. Direct `mw.X` accesses are all reported.
```powershell
python scripts/check_api.py --app ../python_molecular_editor --plugin plugins/MyPlugin --default-allowlist
```

### Suppress manager false positives only (recommended for routine use)
```powershell
python scripts/check_api.py --app ../python_molecular_editor --plugin plugins/MyPlugin --default-allowlist
```

### Also suppress legacy mw.X compat bridges (opt-in)
Use this when you know the `mw.host` / `mw.view3d` patterns are intentional and don't want noise from them:
```powershell
python scripts/check_api.py --app ../python_molecular_editor --plugin plugins/MyPlugin --default-allowlist --mw-allowlist
```

### Suppress try-block issues
```powershell
python scripts/check_api.py --app ../python_molecular_editor --plugin plugins/MyPlugin --default-allowlist --skip-try
```

### Registry Scan
Check all active plugins listed in the registry:
```powershell
python scripts/check_api.py --app ../python_molecular_editor --registry --default-allowlist
```

### Intra-App Consistency Checks
While designed for plugins, the tool can also be used to verify that the core application's own internal logic correctly calls the refactored manager API:
```powershell
python scripts/check_api.py --app ../python_molecular_editor --plugin ../python_molecular_editor/moleditpy/src --default-allowlist
```

### Show the detected API surface
If you want to see exactly what the tool thinks the "valid API" is:
```powershell
python scripts/check_api.py --app ../python_molecular_editor --show-api
```

---

## Allowlist Design

There are two separate allowlists with different risk profiles:

| Flag | Allowlist | What it suppresses | Risk |
|---|---|---|---|
| `--default-allowlist` | `_MANAGER_ALLOWLIST` | Manager attrs set via `self.host.manager.X = ...` (AST-invisible) | **Low** — confirmed false positives |
| `--mw-allowlist` | `_MW_ALLOWLIST` | Direct `mw.X` legacy compat bridge attrs | **High** — can hide real V3 migration bugs |

**Rule of thumb**: always pass `--default-allowlist`. Only pass `--mw-allowlist` when auditing a plugin that you know intentionally uses legacy compat patterns.

---

## Maintenance & Troubleshooting

### Scanned Attribute "Not Found"
If a valid attribute is being flagged as missing, follow these steps:

1.  **Check Visibility**: Ensure the attribute is assigned to `self` in the Manager's `__init__` or a method.
    *   *Bad*: `self.data = ...` (assigned in a hidden helper).
    *   *Better*: Add `self.data: Optional[MoleculeData] = None` to `__init__`.
2.  **Check for Host-Delegation**: If the attribute is set via `self.host.attr = ...`, ensure the script is scanning the file where that assignment happens.
3.  **Update the Manager Allowlist**: If the attribute is genuinely runtime-injected and confirmed valid, add it to `_MANAGER_ALLOWLIST` in `check_api.py`. Do **not** put it in `_MW_ALLOWLIST` unless it is a direct `mw.X` legacy bridge.

### Error Codes
- `[UNKNOWN_MW_ATTR]`: The attribute was called directly on the `MainWindow`.
- `[UNKNOWN_MANAGER_ATTR]`: The attribute was called on a sub-manager (e.g., `mw.io_manager.xxx`).
- `[UNKNOWN_CONTEXT_ATTR]`: The attribute was called on the `PluginContext` object.
