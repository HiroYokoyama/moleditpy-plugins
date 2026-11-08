---
name: Request Plugin Registration / Update
about: Submit a request to register a new remote plugin or update an existing one in the MoleditPy registry.
title: "[Registration Request] <Plugin Name>"
labels: [plugin-registration]
assignees: []

---

Please provide the details below to request the registration of a new plugin or an update to an existing one. Maintainers will review the submission and run the automated registration workflow to update `REGISTRY/plugins.json`.

> [!NOTE]
> **Declare metadata constants in your plugin code (`PLUGIN_*`)**
> The registry reads metadata directly from your release code (`PLUGIN_*` constants) on registration and on every version release. Please define all required and applicable constants directly in your plugin code so that future updates seamlessly inherit your metadata.

### 1. Target Action
* [ ] **Register a New Plugin**
* [ ] **Update an Existing Plugin (New Release / Version Bump)**

### 2. Plugin Information
* **Plugin Name:** 
* **Plugin ID:** 
  *(Unique lowercase identifier, e.g. `my_awesome_plugin`. For new plugins, make it the same as the file name without `.py` or package name; if left blank, this will automatically be derived from the asset filename.)*

### 3. GitHub Repository
Link to the public GitHub repository hosting the plugin source code:
* **Repository URL:** `https://github.com/<owner>/<repo>`
  *(The repository owner must match your GitHub username and the `PLUGIN_AUTHOR` metadata constant.)*

### 4. Plugin Description & Overview
Provide a brief description of the plugin and the scientific workflow or problem it solves (should match `PLUGIN_DESCRIPTION` in code):
* **Description:** 
* **Key Features / Workflow:** 

### 5. Plugin Source URL (Release Asset)
Provide the direct link to your plugin's GitHub Release asset (must end in `.py` or `.zip`):
* **Release Asset URL:** 

### 6. Expected SHA-256 Hash
To verify integrity and security, calculate and paste the SHA-256 hash of the release file:
* **SHA-256:** 
  *(Strictly required for security validation.)*

### 7. Tags (Optional)
Comma-separated categories for the Plugin Manager (e.g., `Visualization, Analysis, Utility, File, Calculation`). If used, please declare in your code as `PLUGIN_TAGS`:
* **Tags:** 

### 8. Dependencies (Optional)
Required and optional Python packages. If external packages are needed, please declare them in your code as `PLUGIN_DEPENDENCIES` (and `PLUGIN_OPTIONAL_DEPENDENCIES` if applicable):
* **Dependencies:** 
* **Optional Dependencies:** 

### 9. Visibility
* [ ] Yes, make this plugin visible in the registry manager.
* [ ] No, keep it hidden for now.

### 10. Supported Versions (Required)
Please declare both `PLUGIN_SUPPORTED_MOLEDITPY_VERSION` and `PLUGIN_SUPPORTED_PYTHON_VERSION` in your plugin code:
* **Supported MoleditPy Version:** *(Required, e.g., `>=5.0.0, <6.0.0`)*
* **Supported Python Version:** *(Required, e.g., `>=3.10, <3.15`)*

---

### Verification Instructions for Developers

1. **Declare metadata constants in your plugin code** (`__init__.py` or single `.py` file):
   * **Required constants:**
     * `PLUGIN_NAME = "..."`
     * `PLUGIN_VERSION = "..."` (must match the release tag version, e.g. `1.2.0`)
     * `PLUGIN_AUTHOR = "..."` (must match your GitHub username / repository owner)
     * `PLUGIN_DESCRIPTION = "..."` (brief summary of functionality)
     * `PLUGIN_SUPPORTED_MOLEDITPY_VERSION = "..."` (e.g. `>=5.0.0, <6.0.0`)
     * `PLUGIN_SUPPORTED_PYTHON_VERSION = "..."` (e.g. `>=3.10, <3.15`)
   * **Optional constants (declare if applicable):**
     * `PLUGIN_LICENSE = "..."` (SPDX identifier, e.g. `LGPL-3.0-or-later`, `MIT`)
     * `PLUGIN_TAGS = [...]` (e.g. `["Visualization", "Analysis"]`)
     * `PLUGIN_DEPENDENCIES = [...]` (e.g. `["numpy>=1.20", "rdkit>=2022.03"]`)
     * `PLUGIN_OPTIONAL_DEPENDENCIES = [...]` (optional, for extra-feature packages)
     * `PLUGIN_SUPPORTED_OS = [...]` (optional, e.g. `["Windows", "macOS", "Linux", "WSL"]` if OS-restricted)
2. **Calculate the SHA-256 hash** of the release asset:
   * **Windows (PowerShell)**: `Get-FileHash ./plugin_file.py -Algorithm SHA256`
   * **Linux / macOS (Bash)**: `shasum -a 256 plugin_file.py`
