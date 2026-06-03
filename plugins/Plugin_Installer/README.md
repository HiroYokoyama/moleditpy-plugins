# Plugin Installer

A built-in utility plugin for MoleditPy to manage, install, and update plugins directly from the registry.

## Features

- **Registry Syncing**: Fetches list of available plugins from the remote registry automatically.
- **Application Update Detection**: Detects newer versions of MoleditPy on PyPI and allows copying the upgrade command to the clipboard.
- **Dependency Constraint Validation**: Parses and verifies PEP-508 style dependency constraints (e.g. `numpy>=1.20`, `rdkit>=2022.03,<2023.0`, `rdkit~=2022.03.1`).
- **Installed vs. Missing Indicators**: Highlights which dependencies are installed or missing inside the details dialog.
- **Safe Command Formatting**: Quotes dependency version requirements properly (e.g., `pip install "numpy>=1.20"`) to prevent shell redirection errors and rejects unsafe execution patterns.
- **Warning-Only Version Constraints**: Prompts warnings for incompatible application versions but does not block the user from proceeding.
- **State Preservation**: Wipes old package folders while backing up and restoring plugin `settings.json` configurations seamlessly.

## Metadata Constants Used

The installer looks up these constants defined at the top level of plugin files:
- `PLUGIN_NAME`: The unique indexing name.
- `PLUGIN_VERSION`: The current version of the plugin.
- `PLUGIN_AUTHOR`: The plugin author.
- `PLUGIN_DESCRIPTION`: Description of the plugin.
- `PLUGIN_SUPPORTED_MOLEDITPY_VERSION`: Version specifier matching the current MoleditPy app version.
- `PLUGIN_TAGS`: A list of categories/tags.
- `PLUGIN_DEPENDENCIES`: A list of package requirement strings.

## Internal Structure & Testing

Helper checking functions are located at the module-level to facilitate headless unit tests:
- `parse_dependency(dep_str)`: Parses a PEP-508 string into name and specifier, rejecting colon formats.
- `check_dependency_satisfied(dep_str)`: Looks up installed distribution versions and matches against requirements.
- `sanitize_and_quote_dependency(dep_str)`: Securely formats terminal command segments.
