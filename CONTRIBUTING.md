# Contributing to MoleditPy Plugins

Thank you for your interest in contributing to the MoleditPy plugin collection! This repository serves as a hub for official and community-developed plugins.

## Plugin Development Requirements

All plugins must follow the standard MoleditPy plugin structure. For detailed API documentation, please refer to the [Plugin Development Manual](https://github.com/HiroYokoyama/python_molecular_editor/blob/main/docs/PLUGIN_DEVELOPMENT_MANUAL_V3.md).

### Basic Structure
A plugin should define its metadata constants and an `initialize` function:

```python
PLUGIN_NAME = "My Plugin"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "YourGitHubUsername"                  # Must match your GitHub username
PLUGIN_DESCRIPTION = "Brief description of what the plugin does."
PLUGIN_TAGS = ["Visualization", "Utility"]             # Optional: parsed automatically by script
PLUGIN_DEPENDENCIES = ["numpy", "rdkit"]                 # Optional: parsed automatically by script

def initialize(context):
    # Register hooks, menus, and tools using the context API
    context.add_menu_action("My Plugin/Action", my_callback)
```

## Contribution Workflow

Contributions are handled differently depending on whether the plugin's source code is hosted inside this repository (intra-repo) or externally (inter-repo).

### 1. External (Inter-Repo) Plugins
*Plugins hosted in external repositories, typically distributed via GitHub Releases.*

All new registrations and updates to external plugins **must be requested by opening a GitHub Issue**. Direct Pull Requests modifying the registry file (`REGISTRY/plugins.json`) for external URLs will be closed.

#### Standard Flow for External Developers:
1.  **Host your plugin**: Upload your plugin file (`.py`) or package (`.zip`) to a public repository. **The plugin must be hosted in your own GitHub repository** (the repository owner must match your GitHub username and the `PLUGIN_AUTHOR` metadata constant). We recommend using GitHub Releases to host stable tag versions (e.g. `v1.2.0` or `1.2.0`).
2.  **Ensure Metadata Consistency**: Your plugin code **must** define the required metadata constants (`PLUGIN_NAME`, `PLUGIN_VERSION`, etc.) at the top of your python file. The version constant must match the release tag version.
3.  **Calculate SHA-256**: Calculate the SHA-256 hash of your release file.
4.  **Open a Registration Issue**: Open a new issue on GitHub using the **Request Plugin Registration / Update** template.
5.  **Review**: A repository maintainer will trigger the automated workflow with your inputs, which downloads, strictly validates (version tag alignment and SHA-256 match), and commits the update to `plugins.json`.

### 2. Internal (Intra-Repo) Plugins
*Plugins whose source code lives directly inside the `plugins/` directory of this repository.*

> [!IMPORTANT]
> **No New Internal Plugins by Third Parties**: Third-party developers are **not allowed to add new plugins** to the `plugins/` directory. Any new third-party plugins must be hosted in your own GitHub repository and registered as **External (Inter-Repo) Plugins** (see the flow above).
> 
> However, Pull Requests from any contributor to update, bugfix, or improve **existing** internal plugins in the `plugins/` folder are highly welcome.

#### Standard Flow for Updating Internal Plugins:
1.  **Clone / Fork** the repository: `git clone https://github.com/HiroYokoyama/moleditpy-plugins.git`
2.  **Create a branch**: `git checkout -b update/plugin-name`
3.  **Update Source**: Modify the existing plugin file or folder located in the `plugins/` directory.
4.  **Update Registry**: Run `python scripts/update_intra_repo_metadata.py` to automatically update the version, SHA-256, and timestamps in `REGISTRY/plugins.json` based on your source code.
5.  **Submit PR**: Commit and push your changes to your branch and open a Pull Request for review.

## Registering in `plugins.json`

The `REGISTRY/plugins.json` file is the registry that the application uses to discover and install plugins. 

### Field Descriptions
- `id`: A unique string identifier (e.g., `my_plugin_id`).
- `visible`: Whether the plugin is shown in the registry interface (`true` or `false`).
- `supported_moleditpy_version`: The supported version of MoleditPy (e.g., `3.*`). Required for visible plugins.
- `name`: **Critical Identifier**. This must be unique across all plugins. It is used as the primary identifier by the application for indexing, updates, and UI display.
- `version`: Version string (e.g., `2026.02.18` or `1.0.0`).
- `author`: Your GitHub username (must match your release repository owner).
- `authorUrl`: URL to your GitHub profile.
- `projectUrl`: URL to the plugin's source code repository.
- `description`: A clear, concise description.
- `tags`: Array of categories (e.g., `["Analysis", "Visualization"]`).
- `dependencies`: Array of required Python packages.
- `downloadUrl`: 
    - For external: `https://github.com/.../my_plugin.zip`
    - For internal: `../plugins/My_Plugin/my_plugin.py`
- `sha256`: The SHA256 hash of the downloaded file.
- `lastUpdated`: The date of the last update (YYYY-MM-DD).
- `firstAppeared`: The date the plugin was first added to the registry (YYYY-MM-DD).

### Example JSON Block
When adding or updating a plugin, use the following structure:

```json
{
  "id": "my_awesome_plugin",
  "visible": true,
  "supported_moleditpy_version": "3.*",
  "name": "My Awesome Plugin",
  "version": "1.0.0",
  "author": "Your Name",
  "authorUrl": "https://github.com/YourUsername",
  "projectUrl": "https://github.com/YourUsername/my-awesome-plugin",
  "description": "A brief description of what your plugin does.",
  "tags": ["Utility", "Visualization"],
  "dependencies": ["rdkit", "numpy"],
  "downloadUrl": "https://github.com/YourUsername/my-awesome-plugin/releases/download/v1.0.0/my_plugin.py",
  "sha256": "YOUR_SHA256_HASH_HERE",
  "lastUpdated": "2026-02-18",
  "firstAppeared": "2026-02-18"
}
```

> [!IMPORTANT]
> **Security Verification**: The **Plugin Installer** plugin automatically checks the `sha256` hash after downloading a plugin. If the hash does not match, the installation will be aborted to protect users from tampered files.

### Calculating SHA256
On Windows (PowerShell):
```powershell
Get-FileHash ./my_plugin.py -Algorithm SHA256
```
On Linux/macOS:
```bash
shasum -a 256 my_plugin.py
```

## Best Practices
- **Isolation**: Ensure your plugin does not interfere with core application functions unless intended.
- **Error Handling**: Wrap UI callbacks in `try...except` blocks to prevent application crashes.
- **Dependencies**: Keep dependencies to a minimum. If you use external libraries, list them in the `dependencies` field in `plugins.json`.

## Security and Verification Policy

To ensure the safety and integrity of the MoleditPy ecosystem, the following policies are strictly enforced:

1.  **Human Verification Required**: Only plugins that have been manually reviewed and verified by a human maintainer are acceptable for inclusion in the official registry (`plugins.json`).
2.  **No Autonomous AI Submissions**: Direct pushes or automated registry updates by autonomous AI agents (e.g., **OpenClaw**) are strictly prohibited.
3.  **Mandatory Registration Request**: All contributions must be submitted via a **GitHub Issue** using the registration template. Direct Pull Requests to edit the registry file bypass our security validations and will be closed.
4.  **Registry Access**: The application's **Plugin Installer** will only provide access to plugins that have passed this verification process.
5.  **Repository Ownership**: All registered remote plugins must be hosted in the contributor's own GitHub repository. The repository owner name parsed from the release URL must match the `PLUGIN_AUTHOR` constant declared in the plugin's source code.
6.  **Maintainer Plugin Deprecation**: If a new, more powerful external plugin is contributed by the community, the maintainer (`HiroYokoyama`) will deprecate and hide their own corresponding official plugin (by setting its `"visible"` flag to `false` in the registry) in favor of the superior community alternative.


