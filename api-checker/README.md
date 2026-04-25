# MoleditPy Static API Checker

This directory contains the tools used to maintain API compatibility between the MoleditPy core application and its extensive plugin ecosystem. These tools use **Static Analysis** (via the Python `ast` module) to detect broken attribute accesses without actually executing the code.

---

## Core Scripts

### `plugin_api_checker.py` (The Engine)
The core scanning logic. It is a **standalone, portable tool** that requires explicit paths for both the application and the plugin. Use this for general-purpose auditing or when working on a single plugin outside of this repository.

### `check_api.py` (The Repo Runner)
A **repository-specific wrapper** pre-configured for the `moleditpy-plugins` environment. Its primary advantages are:
- **Registry Support**: Can scan all plugins listed in `REGISTRY/plugins.json` in one pass.
- **Pre-configured Paths**: Defaults to the standard relative paths for this workspace layout.
- **CI Integration**: Used by the automated test suite to ensure the entire collection remains compatible.

---

## Technical Architecture

### Phase 1: API Surface Extraction
Before checking any plugins, the tool builds a truth map of the application API:

1. **MainWindow Analysis**: Parses `main_window.py` to index all methods, class attributes, and signals.
2. **Manager Discovery**: Identifies attributes on `MainWindow` that are assigned instances of known Manager classes (e.g. `View3DManager`, `StateManager`).
3. **Cross-File Parsing**: For every discovered Manager, locates the source file and indexes public methods, `@property` getters, and attributes assigned to `self` in `__init__` and setup methods.
4. **Host-Delegation Support**: Scans every `.py` file in the core app for the `self.host.xxx = ...` pattern, handling cases where a Manager dynamically adds an attribute to `MainWindow`.

### Phase 2: Plugin Verification
Plugin files are scanned in two passes:

1. **Pass 1 (Alias Tracking)**: Identifies which variables hold references to the application — standard names (`mw`, `main_window`), method calls (`context.get_main_window()`), and propagated aliases (`self.mw = mw`).
2. **Pass 2 (Access Checking)**: Examines every attribute access.
    - **Direct access**: `mw.some_method` must exist on `MainWindow`.
    - **Chained access**: `mw.manager.attr` must verify both `manager` on `MainWindow` and `attr` on the Manager class.

---

## Allowlist Design

There are three allowlist mechanisms with different scopes and risk profiles:

| Source | What it suppresses | Activated by | Risk |
|---|---|---|---|
| `_MANAGER_ALLOWLIST` | Manager attrs set via `self.host.manager.X = ...` (AST-invisible) | `--default-allowlist` | **Low** — confirmed false positives |
| `_MW_ALLOWLIST` | Direct `mw.X` legacy compat bridge attrs (`mw.host`, `mw.view3d`, …) | `--mw-allowlist` | **High** — can hide real V3 migration bugs |
| `.moleditpy-api-allowlist` | Per-repo site-specific attrs (e.g. `hasattr`-guarded V2 compat patterns) | Auto-discovered | **Low** — opt-in per plugin repo |

**Rule of thumb**: always pass `--default-allowlist`. Only pass `--mw-allowlist` when you know the plugin intentionally uses V2 compat patterns. Use `.moleditpy-api-allowlist` in the plugin repo for per-repo suppressions.

### `.moleditpy-api-allowlist` — per-repo site allowlist

Place this JSON file in the plugin repo root. Both `plugin_api_checker.py` and `check_api.py` auto-discover it by walking up from the `--plugin` path — no flag required.

Each value accepts a **list** (attrs only) or an **object** (`attr → reason`). The reason strings are ignored by the checker and serve as inline comments explaining why each attr is safe to skip:

```json
{
  "mw": {
    "select_all": "hasattr-guarded V2 compat — interaction.py:702",
    "last_open_path": "hasattr-guarded V2 compat — mode_manager.py:817"
  },
  "manager": {
    "state_manager": {
      "data": "runtime-injected by molecule loader, not in __init__"
    }
  }
}
```

List form also accepted (no comments):

```json
{ "mw": ["attr1", "attr2"] }
```

Use this for attrs that are accessed via `hasattr()` guards in the plugin — the scanner cannot see through `hasattr` at the static level, so these are safe false positives specific to that repo.

Example: `moleditpy_reaction_sketcher_plugin/.moleditpy-api-allowlist` suppresses 6 V2 compat attrs that are all guarded by `hasattr(mw, ...)` in `interaction.py`, `mode_manager.py`, and `patcher.py`.

---

## Usage Guide

### For External Developers (Recommended)
Use `plugin_api_checker.py` for standalone plugin development. It is portable and requires explicit paths.

**Routine scan** (suppress manager false positives):
```powershell
python api-checker/plugin_api_checker.py --app path/to/python_molecular_editor --plugin path/to/MyPlugin --default-allowlist
```

**Full audit** (also suppress legacy `mw.X` compat bridges):
```powershell
python api-checker/plugin_api_checker.py --app path/to/python_molecular_editor --plugin path/to/MyPlugin --default-allowlist --mw-allowlist
```

---

### For `moleditpy-plugins` Contributors (Internal)
The primary way to verify the repository is to scan all **visible** plugins listed in the registry.

**Full Registry Scan** (Recommended for maintenance):
This command parses `REGISTRY/plugins.json` and scans only the plugins where `"visible": true`.
```powershell
python api-checker/check_api.py --app ../python_molecular_editor --registry --default-allowlist --mw-allowlist
```

**Individual Plugin Scan**:
Use this while developing a specific plugin before adding it to the registry.
```powershell
python api-checker/check_api.py --app ../python_molecular_editor --plugin plugins/MyPlugin --default-allowlist
```

**Main app source/test scan**:
Verify that the core application's own internal logic correctly uses the V3 manager API.
```powershell
# Using the generic engine for intra-app audit:
python api-checker/plugin_api_checker.py --app ../python_molecular_editor --plugin ../python_molecular_editor/moleditpy/src --default-allowlist

# Scan the main app's test suite:
python api-checker/plugin_api_checker.py --app ../python_molecular_editor --plugin ../python_molecular_editor/tests --default-allowlist
```

### Show the detected API surface
```powershell
python api-checker/plugin_api_checker.py --app path/to/python_molecular_editor --show-api
```

---

## Handling False Positives

Static analysis cannot see everything that happens at runtime. The tool uses several mechanisms to reduce noise:

- **Inherited Methods**: A built-in list (`_QT_INHERITED`) of standard `QWidget`/`QMainWindow` methods (`setVisible`, `show`, `resize`, …) that are always permitted.
- **Try-Block Heuristic**: Accesses inside a `try: ... except:` block are flagged with a `[try]` tag. These often indicate the developer has provided a fallback and are safe to suppress with `--skip-try`.
- **Manager Allowlist** (`--default-allowlist`): Suppresses manager attrs set via `self.host.manager.X = ...` patterns invisible to AST (e.g. `init_manager.scene`, `state_manager.data`). Safe to enable in all scans.
- **MW Allowlist** (`--mw-allowlist`): Suppresses direct `mw.X` legacy compat bridge attrs. Off by default — enabling it can hide real V3 migration bugs.

### `mw.host` and `mw.view3d` — why they are false positives

The scanner flags these as `UNKNOWN_MW_ATTR` because neither exists as a declared attribute on `MainWindow` in V3. However, all plugin usages are guarded by `hasattr`:

```python
def run(mw):
    if hasattr(mw, 'host'):      # V2 compat: unwrap host wrapper object
        mw = mw.host
```

```python
if hasattr(mw, 'view3d'):        # V2 compat: optional direct 3D view shortcut
    mw.view3d.draw_standard_3d_style(mol, style_override='ball_and_stick')
```

- `mw.host` — V2-era pattern where plugins were passed a wrapper object holding the real `MainWindow` as `.host`. In V3 the wrapper is gone; `hasattr` returns `False` and the unwrap is skipped.
- `mw.view3d` — V2 shortcut to the 3D view widget, moved to `mw.view_3d_manager` in V3. The `hasattr` guard makes the old call a no-op, not a crash.

Suppress with `--mw-allowlist` for routine scans. Leave visible when auditing a plugin that should be fully migrated to V3.

### Scanning Test Code

Scanning a test directory produces false positives from patterns the AST scanner cannot see:

**Pattern 1 — `cls.X = lambda` monkey-patching** (e.g. `gui/conftest.py`)

Test conftest files patch `MainWindow` with proxy lambdas for the duration of the test session:

```python
cls.halt_all_calculations    = lambda self: self.compute_manager.halt_conversion()
cls.toggle_atom_info_display = lambda self, m: self.view_3d_manager.toggle_atom_info_display(m)
cls.close_all_3d_edit_dialogs = lambda self: self.edit_3d_manager.close_all_3d_edit_dialogs()
```

When those methods are later called on `main_window`, the scanner flags them as `UNKNOWN_MW_ATTR` — it saw the usage but not the `cls.X = lambda` definition, because the scanner only tracks `def X(self)` style definitions.

**Pattern 2 — `MainWindow` subclassing** (e.g. `unit/test_app_state_persistence.py`)

Unit tests sometimes define a thin subclass that forwards manager methods:

```python
class TestableMainWindow(MainWindow):
    def create_json_data(self): return super().create_json_data()
    def clear_2d_editor(self, push_to_undo=True): ...
```

The scanner sees `mw.create_json_data()` and flags it, because those methods are not declared on `MainWindow` itself.

**Verdict**: These are false positives. Do not add them to `_MW_ALLOWLIST` and do not fix the test code. Ignore `UNKNOWN_MW_ATTR` hits when scanning test directories, or pass `--mw-allowlist` scoped to test runs only.

---

## Maintenance & Troubleshooting

### Attribute flagged as "Not Found" but it is valid

1. **Check visibility**: Ensure the attribute is assigned to `self` in the Manager's `__init__` or a named setup method — not inside a helper called indirectly.
2. **Check host-delegation**: If set via `self.host.attr = ...`, ensure the scanner is pointed at the file containing that assignment.
3. **Update the Manager Allowlist**: If the attribute is genuinely runtime-injected and confirmed valid, add it to `_MANAGER_ALLOWLIST` in `check_api.py`. Do **not** use `_MW_ALLOWLIST` unless it is a direct `mw.X` legacy bridge.

### Error Codes

| Code | Meaning |
|---|---|
| `[UNKNOWN_MW_ATTR]` | Attribute called directly on `MainWindow` — not found in its API surface. |
| `[UNKNOWN_MANAGER_ATTR]` | Attribute called on a sub-manager (e.g. `mw.io_manager.xxx`) — not found on that Manager class. |
| `[UNKNOWN_CONTEXT_ATTR]` | Attribute called on the `PluginContext` object — not found in the plugin interface. |
| `[try]` prefix | The access is inside a `try: ... except:` block — likely has a fallback; lower priority. |
