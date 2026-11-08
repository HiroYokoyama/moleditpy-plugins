---
name: Request Plugin Registration / Update
about: Submit a request to register a new remote plugin or update an existing one in the MoleditPy registry.
title: "[Registration Request] <Plugin Name>"
labels: [plugin-registration]
assignees: []

---

Please provide the details below to request the registration of a new plugin or an update to an existing one. Maintainers will review the submission and run the automated registration workflow to update `REGISTRY/plugins.json`.

### 1. Target Action
* [ ] **Register a New Plugin**
* [ ] **Update an Existing Plugin (New Release / Version Bump)**

*Note: Name, description, tags, dependencies, and supported versions are read directly from your release code on every version bump — to change any of them, update the `PLUGIN_*` constants and publish a new release.*

### 2. Plugin Information
* **Plugin Name:** 
* **Plugin ID:** 
  *(Unique lowercase identifier, e.g. `my_awesome_plugin`. For new plugins, make it the same as the file name without `.py` or package name; if left blank, this will automatically be derived from the asset filename.)*

### 3. GitHub Repository
Link to the public GitHub repository hosting the plugin source code:
* **Repository URL:** `https://github.com/<owner>/<repo>`
  *(The repository owner must match your GitHub username and the `PLUGIN_AUTHOR` metadata constant.)*

### 4. Plugin Description & Overview
Provide a brief description of the plugin and the scientific workflow or problem it solves:
* **Description:** 
* **Key Features / Workflow:** 
  - 

### 5. Plugin Source URL (Release Asset)
Provide the direct link to your plugin's GitHub Release asset (must end in `.py` or `.zip`):
* **Release Asset URL:** 

### 6. Expected SHA-256 Hash
To verify integrity and security, calculate and paste the SHA-256 hash of the release file:
* **SHA-256:** 
  *(Strictly required for security validation.)*

### 7. Tags (Optional)
Comma-separated categories for the Plugin Manager (e.g., `Visualization, Analysis, Utility, File, Calculation`). Only needed if your code does not declare `PLUGIN_TAGS` in code:
* **Tags:** 

### 8. Dependencies (Optional)
Comma-separated list of required Python packages (e.g., `numpy, rdkit, PySide6`):
* **Dependencies:** 

### 9. Visibility
* [ ] Yes, make this plugin visible in the registry manager.
* [ ] No, keep it hidden for now.

### 10. Supported Versions
* **Supported MoleditPy Version:** *(e.g., `>=5.0.0, <6.0.0`; required if visible)*
* **Supported Python Version (optional):** *(e.g., `>=3.10, <3.15`; defaults to `>=3.10, <3.15` if omitted)*

---

### Verification Instructions for Developers

1. Ensure your plugin defines the required metadata constants at the top of your python file (`__init__.py` or single `.py` file):
   * `PLUGIN_NAME`
   * `PLUGIN_VERSION` (must match the release tag version, e.g. `1.2.0`)
   * `PLUGIN_AUTHOR` (must match your GitHub username/repository owner)
   * `PLUGIN_DESCRIPTION`
2. Calculate the SHA-256 hash of the release asset:
   * **Windows (PowerShell)**: `Get-FileHash ./plugin_file.py -Algorithm SHA256`
   * **Linux / macOS (Bash)**: `shasum -a 256 plugin_file.py`
