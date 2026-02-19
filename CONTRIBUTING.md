# Contributing to MoleditPy Plugins

Thank you for your interest in contributing to the MoleditPy plugin collection! This repository serves as a hub for official and community-developed plugins.

## Plugin Development Requirements

All plugins must follow the standard MoleditPy plugin structure. For detailed API documentation, please refer to the [Plugin Development Manual](PLUGIN_DEVELOPMENT_MANUAL.md).

### Basic Structure
A plugin should define at least the following metadata and an `initialize` function:

```python
PLUGIN_NAME = "My Plugin"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "Your Name"
PLUGIN_DESCRIPTION = "Brief description of what the plugin does."

def initialize(context):
    # Register hooks, menus, and tools using the context API
    context.add_menu_action("My Plugin/Action", my_callback)
```

## Contribution Workflow

All contributions—whether they are new plugin registrations or updates to existing ones—must be submitted via **Pull Request**.

### Standard Git Flow
1.  **Clone** the repository: `git clone https://github.com/HiroYokoyama/moleditpy-plugins.git`
2.  **Pull** the latest changes: `git pull origin main`
3.  **Create a branch** for your changes: `git checkout -b feature/my-new-plugin`
4.  **Implement** changes and **Commit**: `git commit -am "Add my new plugin"`
5.  **Push** and create a **Pull Request** on GitHub.

### 1. New Plugins
**Rule: New plugins must be hosted in external repositories.** 
To add a new plugin to the collection:
1.  **Host your plugin**: Upload your plugin file (`.py`) or package (`.zip`) to a public repository (e.g., GitHub). We recommend using GitHub Releases to host stable versions.
2.  **Register the plugin**: Add a new entry to [explorer/plugins.json](explorer/plugins.json) via a Pull Request.
3.  **Download URL**: The `downloadUrl` in `plugins.json` should point to the direct download link (e.g., a GitHub Release asset or a raw file URL).

### 2. Updates to Existing Plugins
**Rule: Updates to existing plugins already in this repository are welcome via Pull Request.**
If you are updating a plugin that is already located in the `plugins/` directory:
1.  **Update the source**: Replace the `.py` file or folder content in the `plugins/` directory.
2.  **Update Metadata**: Increment the version number and update the metadata in the script.
3.  **Update plugins.json**: Update the corresponding entry in `explorer/plugins.json` with the new version, `lastUpdated` timestamp, and a fresh `sha256` hash.
4.  **Submit PR**: Push your changes to a branch and open a Pull Request.

## Registering in `plugins.json`

The `explorer/plugins.json` file is the registry that the application uses to discover and install plugins. 

### Field Descriptions
- `id`: A unique string identifier (e.g., `my_plugin_id`).
- `name`: **Critical Identifier**. This must be unique across all plugins. It is used as the primary identifier by the application for indexing, updates, and UI display.
- `version`: Version string (e.g., `2026.02.18` or `1.0.0`).
- `author`: Your name.
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
3.  **Mandatory Pull Requests**: All contributions, whether human-authored or AI-assisted, must be submitted via a **Pull Request** for human review. Autonomous agent activity bypasses our security model and will be rejected or reverted.
4.  **Registry Access**: The application's **Plugin Installer** will only provide access to plugins that have passed this verification process.

