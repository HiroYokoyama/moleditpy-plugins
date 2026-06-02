---
name: Request Plugin Registration / Update
about: Submit a request to register a new remote plugin or update an existing one in the MoleditPy registry.
title: "[Registration Request] <Plugin Name>"
labels: [plugin-registration]
assignees: []

---

Please provide the details below to request the registration of a new plugin or an update to an existing one. Maintainers will run the automated registration workflow to update `REGISTRY/plugins.json`.

### 1. Plugin Source URL
Provide the direct link to your plugin's GitHub Release asset (must end in `.py` or `.zip`).
**URL:** 

### 2. Target Action
* [ ] Register a New Plugin
* [ ] Update an Existing Plugin

### 3. Plugin ID
The unique lowercase identifier (e.g., `my_awesome_plugin`).
*For new plugins, if left blank, this will automatically be derived from the filename.*
**Plugin ID:** 

### 4. Expected SHA-256 Hash
To verify integrity and security, please calculate and paste the SHA-256 hash of your release file. 
*Note: This is strictly required for security validation.*
**SHA-256:** 

### 5. Tags (Optional)
Comma-separated categories for the Plugin Manager (e.g., `Visualization, Analysis, Utility, File IO`).
**Tags:** 

### 6. Dependencies (Optional)
Comma-separated list of required Python packages (e.g., `numpy, rdkit, PyQt6`).
**Dependencies:** 

### 7. Visibility
* [ ] Yes, make this plugin visible in the registry manager.
* [ ] No, keep it hidden for now.

### 8. Supported MoleditPy Version
Provide the supported version of MoleditPy (e.g., `3.*`).
*Note: This is strictly required if the plugin is visible.*
**Supported MoleditPy Version:** 

---

### Verification Instructions for Developers

1. Ensure your plugin defines the required metadata constants at the top of your python file (`__init__.py` or single `.py` file):
   * `PLUGIN_NAME`
   * `PLUGIN_VERSION` (must match the release tag version, e.g. `1.2.0`)
   * `PLUGIN_AUTHOR`
   * `PLUGIN_DESCRIPTION`
2. Calculate the SHA-256 hash of the asset:
   * **Windows (PowerShell)**: `Get-FileHash ./plugin_file.py -Algorithm SHA256`
   * **Linux / macOS (Bash)**: `shasum -a 256 plugin_file.py`
