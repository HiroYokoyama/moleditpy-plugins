# Plugin Installer

![coverage](https://img.shields.io/badge/coverage-%3E90%25-brightgreen)

A management plugin for MoleditPy that installs, updates, and removes plugins directly from the official registry.

## Features

- **Non-blocking registry fetch**: Background `QThread` fetches the PyPI app version and remote registry JSON concurrently — the UI stays responsive while loading.
- **Alphabetical plugin table**: All plugins sorted A–Z; status badge and action button on each row.
- **Update count label**: Shows "N plugin(s) have updates available" (amber) or "All plugins are up to date" (green) above the Update All button.
- **Chunked download with progress**: Downloads in 8 KB chunks with a live progress dialog; cancellable mid-transfer.
- **Batch update**: "Update All" iterates every "Update Available" row in one pass with a single batch progress dialog.
- **AST version read**: Installed plugin version is read directly from the `.py` or `__init__.py` file on disk via `ast.parse` — no import, always reflects the actual file on disk.
- **Session install tracking**: Plugins installed during the session are tracked in `_pending_installs` so the table shows the correct status immediately without waiting for `discover_plugins`.
- **Reload only when needed**: Plugin reload (`discover_plugins` + menu rebuild) runs exactly once, via the `finished` signal, only if at least one plugin was installed or updated that session.
- **Application Update Detection**: Detects newer versions of MoleditPy on PyPI and shows a "Copy upgrade command" button.
- **Dependency Constraint Validation**: Parses and verifies PEP-508 style dependency constraints (`numpy>=1.20`, `rdkit~=2022.03.1`, etc.).
- **Python version compatibility check**: The registry's `supported_python_version` spec (e.g. `>=3.9, <3.15`) is checked against the running interpreter — incompatible plugins are marked in the table (tooltip shows the required range) and a warning dialog asks for confirmation before install/update. Entries without the field are treated as compatible.
- **OS compatibility check**: The registry's `supported_os` list (`Windows` / `macOS` / `Linux` / `WSL`) is checked against the running machine, the same way as the Python version — incompatible plugins are marked and confirmed before install. Entries without the field are treated as compatible. WSL reports as Linux, so a plain `Linux` entry satisfies it.
- **Dependency install hints**: Highlights missing vs installed dependencies in the details dialog with pip-safe quoted commands.
- **Optional dependencies**: `optional_dependencies` entries get their own details section — missing ones read "Not installed" (amber) with a separate copy-install-command button. They never trigger the missing-dependency prompt, so an install is never gated on a package the plugin only needs for extras.
- **State Preservation**: Backs up and restores plugin `settings.json` when overwriting an existing installation.
- **Startup check**: On first launch asks whether to enable automatic update checks; thereafter runs a silent background check at startup if opted in.

## Entry Point

This plugin uses the V4 API:

```python
def initialize(context: PluginContext) -> None:
    context.add_menu_action("Plugin/Plugin Installer...", _open_installer)
```

The menu item opens `PluginInstallerWindow` as a modal dialog.

## Network Calls

Exactly two network requests per check cycle (both handled by `_FetchWorker`):

| Request | URL | Purpose |
|---|---|---|
| PyPI metadata | `https://pypi.org/pypi/{package}/json` | Detect app updates |
| Registry JSON | `https://hiroyokoyama.github.io/moleditpy-plugins/REGISTRY/plugins.json` | Plugin list and versions |

Plugin downloads are a third request, one per plugin, made by `_download_chunked`.

## Key Internal Components

| Symbol | Description |
|---|---|
| `_FetchWorker(QThread)` | Background thread: fetches PyPI version + registry JSON, emits `done(pypi_ver, remote_entries)` |
| `_read_plugin_version_ast(filepath)` | Reads `PLUGIN_VERSION` from a `.py` or `__init__.py` via AST |
| `PluginInstallerWindow` | Main dialog (QDialog); opened with `exec()` |
| `_pending_installs` | Dict of plugins installed this session, keyed by name |
| `_needs_plugin_reload` | Set `True` after any successful install; cleared in `_on_finished` |
| `_on_finished(result)` | Connected to `finished` signal — runs `discover_plugins` + menu rebuild once on close |
| `_update_status_label(count)` | Updates the status label above Update All |
| `PluginDetailsDialog._add_dependency_section(...)` | Renders one dependency block (required or optional) and returns its missing entries |
| `_download_chunked(url, dest, on_progress)` | 8 KB chunked download; `on_progress(received, total) -> bool` callback |

## Module-Level Helpers (unit-testable)

- `_read_plugin_version_ast(filepath)` — AST version extraction
- `parse_dependency(dep_str)` — parses PEP-508 string into name and specifier
- `check_dependency_satisfied(dep_str)` — checks installed distributions against requirements
- `sanitize_and_quote_dependency(dep_str)` — formats pip command segments safely
- `is_app_version_compatible(app_version, specifier)` — evaluates version specifier against running app
- `get_python_version()` — running interpreter version as `X.Y.Z` (checked with `is_app_version_compatible` against `supported_python_version`)
- `get_current_os()` — this machine as one of the `supported_os` tokens
- `is_os_compatible(current_os, supported_os)` — checks it against a registry entry's list; a missing or malformed list means compatible
- `format_supported_os(supported_os)` — the list as displayed in the table and warning dialog

## Plugin Metadata Constants

The installer reads these constants from each plugin file (via the registry or AST):

| Constant | Purpose |
|---|---|
| `PLUGIN_NAME` | Unique display and index name |
| `PLUGIN_VERSION` | Current version (date-based `YYYY.MM.DD`) |
| `PLUGIN_AUTHOR` | Author name |
| `PLUGIN_DESCRIPTION` | One-line description |
| `PLUGIN_SUPPORTED_MOLEDITPY_VERSION` | PEP-440 specifier checked against the running app |
| `PLUGIN_SUPPORTED_PYTHON_VERSION` | Optional PEP-440-style specifier checked against the running Python (registry field `supported_python_version`; visible plugins default to `>=3.9, <3.15`) |
| `PLUGIN_SUPPORTED_OS` | Optional list of `Windows` / `macOS` / `Linux` / `WSL` (registry field `supported_os`). Declare it only for a plugin that actually needs an OS-restricted backend; omitting it means every OS. |
| `PLUGIN_TAGS` | List of category strings |
| `PLUGIN_DEPENDENCIES` | List of PEP-508 requirement strings the plugin requires |
| `PLUGIN_OPTIONAL_DEPENDENCIES` | Optional list of PEP-508 requirement strings for extra features (registry field `optional_dependencies`); shown in its own details section and never warned about on install |
