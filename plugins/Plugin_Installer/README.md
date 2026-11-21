# Plugin Installer

![coverage](https://img.shields.io/badge/coverage-%3E90%25-brightgreen)

A management plugin for MoleditPy that installs, updates, and removes plugins directly from the official registry.

## Features

- **Non-blocking registry fetch**: Background `QThread` fetches the PyPI app version and remote registry JSON concurrently — the UI stays responsive while loading.
- **Alphabetical plugin table**: All plugins sorted A–Z; status badge and action button on each row.
- **Hidden installed plugins**: Registry entries with `visible: false` are shown as `Stale` when already installed, with migration guidance in the name/status tooltips and details dialog. Disabled entries show `Stale (Disabled)`. Hidden entries remain excluded from new installs and Update All; manual reinstallation and deletion remain available.
- **Update count label**: Shows "N plugin(s) have updates available" (amber) or "All plugins are up to date" (green) above the Update All button.
- **Chunked download with progress**: Downloads in 8 KB chunks with a live progress dialog; cancellable mid-transfer.
- **Batch update**: "Update All" iterates every "Update Available" row in one pass with a single batch progress dialog.
- **AST version read**: Installed plugin version is read directly from the `.py` or `__init__.py` file on disk via `ast.parse` — no import, always reflects the actual file on disk.
- **Session install tracking**: Plugins installed during the session are tracked in `_pending_installs` so the table shows the correct status immediately without waiting for `discover_plugins`.
- **Reload only when needed**: Plugin reload (`discover_plugins` + menu rebuild) runs exactly once, via the `finished` signal, only if at least one plugin was installed or updated that session.
- **Application Update Detection**: Detects newer versions of MoleditPy on PyPI and shows a "Copy upgrade command" button.
- **Dependency Constraint Validation**: Parses and verifies PEP-508 style dependency constraints (`numpy>=1.20`, `rdkit~=2022.03.1`, etc.).
- **Python version compatibility check**: The registry's `supported_python_version` spec (e.g. `>=3.10, <3.15`) is checked against the running interpreter — incompatible plugins are marked in the table (tooltip shows the required range) and a warning dialog asks for confirmation before install/update. Entries without the field are treated as compatible.
- **OS compatibility check**: The registry's `supported_os` list (`Windows` / `macOS` / `Linux` / `WSL`) is checked against the running machine, the same way as the Python version — incompatible plugins are marked and confirmed before install. Entries without the field are treated as compatible. WSL reports as Linux, so a plain `Linux` entry satisfies it.
- **Dependency install hints**: Highlights missing vs installed dependencies in the details dialog with pip-safe quoted commands.
- **Optional dependencies**: `optional_dependencies` entries get their own details section — missing ones read "Not installed" (amber) with a separate copy-install-command button. They never trigger the missing-dependency prompt, so an install is never gated on a package the plugin only needs for extras.
- **State Preservation**: Backs up and restores plugin `settings.json` when overwriting an existing installation.
- **Startup check**: On first launch asks whether to enable automatic update checks; thereafter runs a silent background check at startup if opted in.
- **Command line**: With MoleditPy 5.3.2 or later, `moleditpy -r plugin_installer` lists, checks, updates and installs plugins without opening the editor (see below).

## Command Line (MoleditPy 5.3.2+)

MoleditPy 5.3.2 runs single-file plugins that have a top-level `if __name__ == "__main__":` block, so the installer works from a terminal. It uses the same registry, compatibility rules, SHA-256 check and install steps as the dialog.

```sh
moleditpy -r plugin_installer list              # installed plugins and their status
moleditpy -r plugin_installer list --all        # also plugins that are not installed
moleditpy -r plugin_installer check             # MoleditPy and plugin updates
moleditpy -r plugin_installer update --all      # update every out-of-date plugin
moleditpy -r plugin_installer update xyz_editor # update named plugins
moleditpy -r plugin_installer install pyscf_calculator
```

| Command | What it does |
|---|---|
| `list [--all] [--json]` | Name, registry id, installed and latest version, and status — the same statuses as the dialog (*Update Available*, *Up to date*, *Stale*, *Disabled*, …). |
| `check [--json] [--exit-code]` | Reports a newer MoleditPy on PyPI, plugins with updates, and stale plugins with their replacement. `--exit-code` exits with `1` when anything can be updated. |
| `update PLUGIN... \| --all` | Updates installed plugins. `--all` takes every plugin marked *Update Available*, disabled ones included and stale ones left out. A named plugin that is up to date is skipped unless `--reinstall` is given. |
| `install PLUGIN...` | Installs plugins that are not installed yet, and lists missing dependencies with a `pip install` command. |

- A plugin is named by its registry id (the `Id` column of `list --all`) or by its name, in any case; quote names with spaces.
- `update` and `install` list what they will do and ask `Proceed? [y/N]`. `-y` / `--yes` answers yes; without a terminal they refuse to run unless `--yes` is given. `--dry-run` only lists.
- A plugin that requires another MoleditPy, Python or OS is skipped with the reason; `--force` installs it anyway. A registry entry without a SHA-256 checksum is refused unless `--allow-unverified` is given. A checksum mismatch always stops that plugin.
- Exit status: `0` success or nothing to do, `1` a plugin failed, was skipped or the confirmation was declined, `2` a usage error, `3` the registry could not be fetched or MoleditPy is older than 5.3.2.
- An editor that is already running keeps the plugin code it loaded; restart MoleditPy after updating.
- An installer older than 2026.11.22 has no command line; update it once from the dialog (or with `moleditpy --install-plugin`) first.

The command line is `cli_main(argv)` at the end of `plugin_installer.py`. It finds installed plugins through `moleditpy.plugins.plugin_runner.installed_plugins` and installs through the host's `PluginManager.install_plugin`, and shares the module-level helpers below with the dialog.

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

Plugin downloads are a third request, one per plugin, made by `_download_file` (wrapped by the dialog's `_download_chunked`). The command line makes the same requests, without the thread.

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
| `_download_chunked(url, dest, on_progress)` | Dialog wrapper around `_download_file`; shows a message box on `DownloadError` |
| `cli_main(argv=None, out=None, err=None)` | Command-line entry point (`list`, `check`, `update`, `install`); returns the exit status |

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
- `compare_versions(v1, v2)` / `get_app_version()` / `get_package_name()` — version comparison and the running app's version and PyPI name
- `_fetch_pypi_version(package)` / `_fetch_registry(url=None)` — the two network fetches, used by `_FetchWorker` and the command line
- `_download_file(url, dest, on_progress)` — 8 KB chunked download; `on_progress(received, total) -> bool` can cancel; raises `DownloadError`
- `_resolve_download_url(url)` / `_download_filename(final_url, target_file)` — the URL and file name a registry entry downloads to
- `_install_download(...)` — installs a downloaded `.py` or ZIP, overwriting a folder plugin while keeping its `settings.json`
- `_classify_plugin(name, remote_info, local_info, app_ver, registry)` — one table row (status, compatibility, stale migration target); `_row_color(row)` gives its colour

## Plugin Metadata Constants

Migration guidance is curated registry metadata: a hidden entry's optional
`migration_to` field contains the replacement plugin's ID. The installer resolves
the ID to its current display name. Without a migration target, only the Stale
notice is shown. Plugins missing from the registry remain `Not in Registry`.
No migration is performed automatically.

Maintainers can record a target without modifying a retired plugin's source:

```bash
python scripts/update_registry_migrations.py --set orca_xyz2inp_gui=orca_input_generator_pro
```

The installer reads these constants from each plugin file (via the registry or AST):

| Constant | Purpose |
|---|---|
| `PLUGIN_NAME` | Unique display and index name |
| `PLUGIN_VERSION` | Current version (date-based `YYYY.MM.DD`) |
| `PLUGIN_AUTHOR` | Author name |
| `PLUGIN_DESCRIPTION` | One-line description |
| `PLUGIN_SUPPORTED_MOLEDITPY_VERSION` | PEP-440 specifier checked against the running app |
| `PLUGIN_SUPPORTED_PYTHON_VERSION` | Optional PEP-440-style specifier checked against the running Python (registry field `supported_python_version`; visible plugins default to `>=3.10, <3.15`) |
| `PLUGIN_SUPPORTED_OS` | Optional list of `Windows` / `macOS` / `Linux` / `WSL` (registry field `supported_os`). Declare it only for a plugin that actually needs an OS-restricted backend; omitting it means every OS. |
| `PLUGIN_TAGS` | List of category strings |
| `PLUGIN_DEPENDENCIES` | List of PEP-508 requirement strings the plugin requires |
| `PLUGIN_OPTIONAL_DEPENDENCIES` | Optional list of PEP-508 requirement strings for extra features (registry field `optional_dependencies`); shown in its own details section and never warned about on install |
