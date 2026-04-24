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
-   **The Allowlist**: Confirmed-valid dynamic attributes (like `splitter`, `edit_3d_action`) are maintained in the `_DEFAULT_ALLOWLIST` in `check_api.py`.

---

## Usage Guide

### Basic Scan
Check a specific plugin against the local app repository:
```powershell
python scripts/check_api.py --app ../python_molecular_editor --plugin plugins/MyPlugin --default-allowlist
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

## Maintenance & Troubleshooting

### Scanned Attribute "Not Found"
If a valid attribute is being flagged as missing, follow these steps:

1.  **Check Visibility**: Ensure the attribute is assigned to `self` in the Manager's `__init__` or a method.
    *   *Bad*: `self.data = ...` (assigned in a hidden helper).
    *   *Better*: Add `self.data: Optional[MoleculeData] = None` to `__init__`.
2.  **Check for Host-Delegation**: If the attribute is set via `self.host.attr = ...`, ensure the script is scanning the file where that assignment happens.
3.  **Update the Allowlist**: If the attribute is injected by external tools (like a `.ui` file or `QtDesigner`), add it to the `_DEFAULT_ALLOWLIST` in `check_api.py`.

### Error Codes
- `[UNKNOWN_MW_ATTR]`: The attribute was called directly on the `MainWindow`.
- `[UNKNOWN_MANAGER_ATTR]`: The attribute was called on a sub-manager (e.g., `mw.io_manager.xxx`).
- `[UNKNOWN_CONTEXT_ATTR]`: The attribute was called on the `PluginContext` object.
