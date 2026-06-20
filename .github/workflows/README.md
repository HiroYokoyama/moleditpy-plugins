# GitHub Actions Workflows Documentation

This directory contains the automated workflows for the MoleditPy Plugins Registry.

---

## 1. Register or Update Remote Plugin (`register-remote-plugin.yml`)

A manual workflow (`workflow_dispatch`) used to automatically register new third-party/remote plugins or update existing remote plugins in the registry (`REGISTRY/plugins.json`).

### Workflow Inputs

| Input Field | Required? | Description |
| :--- | :--- | :--- |
| `release_url` | **Yes** | The direct download URL of the GitHub Release asset (must be a `.py` file or a `.zip` file). Format: `https://github.com/{owner}/{repo}/releases/download/{tag}/{filename}`. |
| `plugin_id` | No | The unique ID for the plugin. For new plugins, if omitted, it will be automatically derived from the release file name (stem). |
| `tags` | No | A comma-separated list of tags (only used when registering a new plugin, e.g., `Analysis, Visualization`). |
| `dependencies` | No | A comma-separated list of required Python packages (only used when registering a new plugin, e.g., `numpy, rdkit`). |
| `visible` | **Yes** | Visibility flag in the registry (`true` or `false`). Defaults to `true`. |
| `expected_sha256` | *Conditional* | The expected SHA-256 hash. **Mandatory** for security verification if the repository owner is not `HiroYokoyama`. |
| `date` | No | Override registration/update date (`YYYY-MM-DD`). If omitted or empty, automatically falls back to the current system date. |
| `supported_version` | No | Supported MoleditPy version (e.g., `3.*`). If omitted, falls back to `PLUGIN_SUPPORTED_MOLEDITPY_VERSION` in the plugin code or the existing registry entry value. |
| `dry_run` | **Yes** | If set to `true`, the workflow performs all downloads and verification checks but **does not** commit or push changes back to the registry. |
---

## Registry Metadata Mappings

When registering or updating a plugin, the entry in `REGISTRY/plugins.json` is generated and mapped as follows:

| Field in `plugins.json` | Source | How it is Derived / Formatted |
| :--- | :--- | :--- |
| `id` | **Input / Derived** | Used directly if `plugin_id` is supplied as a workflow input. If left blank, it is derived from the release file name (stem) converted to lowercase with dashes replaced by underscores. |
| `visible` | **Input** | Directly from the `visible` selection input in the workflow (defaults to `true`). |
| `supported_moleditpy_version` | **Input / Code Constant / Registry** | Prioritizes: 1. `supported_version` input from the workflow/CLI (if provided), 2. `PLUGIN_SUPPORTED_MOLEDITPY_VERSION` defined at the top of the downloaded python file, 3. The existing registry value (when updating). Mandatory for visible plugins. |
| `name` | **Code Constant** | Extracted from `PLUGIN_NAME` defined at the top of the downloaded `.py` or `__init__.py` file. |
| `version` | **Code Constant** | Extracted from `PLUGIN_VERSION` in the code. Normalised to remove leading `v/V`. Checked for tag consistency. |
| `author` | **Code Constant** | Extracted from `PLUGIN_AUTHOR` in the code. Must match the GitHub owner of the repository. |
| `authorUrl` | **Derived** | Generated automatically as `https://github.com/{owner}` where `{owner}` is parsed from the `release_url`. |
| `projectUrl` | **Derived** | Generated automatically as `https://github.com/{owner}/{repo}` where `{owner}/{repo}` is parsed from the `release_url`. |
| `description` | **Code Constant** | Extracted from `PLUGIN_DESCRIPTION` defined in the downloaded code. |
| `tags` | **Code Constant / Input** | Extracted from `PLUGIN_TAGS` list/string in the code. If missing in code, falls back to the `tags` input list from the workflow. |
| `dependencies` | **Code Constant / Input** | Extracted from `PLUGIN_DEPENDENCIES` list/string in the code. If missing in code, falls back to the `dependencies` input list from the workflow. |
| `downloadUrl` | **Input** | Set to the exact provided `release_url`. |
| `sha256` | **Computed** | Calculated as the SHA-256 hash of the downloaded asset file. (Must match `expected_sha256` for external plugins). |
| `lastUpdated` | **System Date / Input** | Set to the provided `date` input (if valid `YYYY-MM-DD`), otherwise defaults to the current system date. |
| `firstAppeared` | **System Date / Input** | Set to the provided `date` input (if valid `YYYY-MM-DD`), otherwise defaults to the current system date (only populated when registering a new plugin). |

---

### Security Policy: Third-Party Repository Verification

To prevent unauthorized updates or malicious code injection from external repositories, the registration script enforces a strict validation policy:

* **Owner = `HiroYokoyama`**: Skip manual SHA-256 verification (it will still check code version consistency and verify the version tag).
* **Owner != `HiroYokoyama`**: A manual `expected_sha256` input **must** be provided. The workflow will download the file, compute its actual SHA-256 hash, and compare it with your manual input. If they do not match, or if the input is missing, the workflow will fail and abort.
  > [!CAUTION]
  > Always download and inspect external plugin code locally, compute its SHA-256, and supply it to the workflow to confirm you trust the release.

---

### Dry Run Mode & Version Warning Bypass

* When testing or verifying a release URL, enable `dry_run: true`.
* If a dry run is executed with a version tag that is **already registered** in `plugins.json`, the validation succeeds but logs a `Warning: [Dry Run] New version 'X.X.X' is equal to existing version...` rather than exiting with an error. This allows maintainers to test integration and integrity checks without version bumps.
* In normal execution (when `dry_run` is `false`), identical or downgraded versions will raise a strict error and abort the update.

---

## 2. Auto-Register Remote Plugin (`auto-register-remote-plugin.yml`)

An automated workflow triggered by a `repository_dispatch` event (`plugin_release`) sent from external plugin repositories when they publish a new GitHub Release. This is the receiving end of the cross-repo release pipeline.

### Trigger

```
Event type : repository_dispatch
types      : [plugin_release]
```

External plugin repos send this event from their own `release.yml` after a successful release (see [Section 4](#4-release-workflow-for-external-plugin-repos-releaseyml) for the sending side).

### Payload

| Field | Description |
| :--- | :--- |
| `client_payload.repo` | Full `owner/repo` of the releasing plugin (e.g. `HiroYokoyama/moleditpy_cif_viewer`) |
| `client_payload.tag` | Release tag without `v` prefix (e.g. `0.11.0`) |

### What It Does

1. Uses `gh api` to look up the GitHub Release for the given `repo` + `tag` and resolves the `.zip` asset download URL automatically.
2. Calls `scripts/register_remote_plugin.py` with that URL — which downloads the zip, extracts all metadata constants (`PLUGIN_VERSION`, `PLUGIN_SUPPORTED_MOLEDITPY_VERSION`, etc.) via AST, validates version consistency, and updates `REGISTRY/plugins.json`.
3. Runs `validate_json.py` and `tests/test_registry.py` to verify the updated registry.
4. Commits and pushes the registry change to `main` if any fields changed.

### No Manual Input Needed

`supported_moleditpy_version` is **read directly from the plugin source** inside the zip — no extra input is required from the caller. The priority order inside `register_remote_plugin.py` is:
1. `--supported-version` CLI flag (not used by this workflow)
2. `PLUGIN_SUPPORTED_MOLEDITPY_VERSION` constant in the downloaded code ← this path
3. Existing registry value (fallback)

### Permissions Required

This workflow uses `GITHUB_TOKEN` (already scoped to `contents: write`) to resolve release assets and push registry changes. No extra secrets needed on the moleditpy-plugins side.

---

## 3. Plugin Tests (`test-plugins.yml`)

An automated CI workflow triggered on every `push` and `pull_request` to `main`. It runs:
1. Registry JSON validation (`scripts/validate_json.py`) to check formatting and prevent duplicate IDs or SHA-256 hashes.
2. Full automated pytest execution on `tests/` (including validation of local files, source metadata conventions, and remote registration script unit tests).
3. API compatibility checks against the main `python_molecular_editor` repository.

---

## 4. Release Workflow for External Plugin Repos (`release.yml`)

Each external plugin repository (e.g. `moleditpy_cif_viewer`, `moleditpy_reaction_sketcher_plugin`) contains its own `.github/workflows/release.yml`. This is the **sending side** of the cross-repo pipeline.

### One-Time Setup: `REGISTRY_PAT` Secret

The release workflow needs permission to dispatch events to this registry repo. Set this up once per external plugin repo:

**Step 1 — Create a Personal Access Token**

1. Go to GitHub → **Settings** → **Developer settings** → **Personal access tokens** → **Fine-grained tokens** (recommended) or **Tokens (classic)**.
2. **Fine-grained token** (recommended):
   - Resource owner: `HiroYokoyama`
   - Repository access: **Only select repositories** → `moleditpy-plugins`
   - Permissions → Repository permissions → **Contents**: `Read and write`
3. **Classic token** (alternative):
   - Scope: `repo` (full control of private repositories)
4. Copy the generated token value.

**Step 2 — Add the Secret to the Plugin Repo**

1. Open the external plugin repository on GitHub (e.g. `HiroYokoyama/moleditpy_cif_viewer`).
2. Go to **Settings** → **Secrets and variables** → **Actions** → **New repository secret**.
3. Name: `REGISTRY_PAT`
4. Value: paste the token from Step 1.
5. Click **Add secret**.

Repeat Step 2 for each external plugin repo. The token itself only needs to be created once.

### How to Publish a Release

1. Update `PLUGIN_VERSION` in the plugin source file and commit + push.
2. Go to the plugin repo on GitHub → **Actions** → **Release** → **Run workflow**.
3. Enter the version number **without a leading `v`** (e.g. `0.11.0`).
4. Click **Run workflow**.

The workflow will:
- Verify the entered version matches `PLUGIN_VERSION` in the source file (fails if mismatched).
- Build a `.zip` of the plugin package.
- Create a GitHub Release with that tag and attach the zip.
- Dispatch a `plugin_release` event to `moleditpy-plugins`, which triggers `auto-register-remote-plugin.yml` to update the registry automatically.

> [!NOTE]
> The registry update step has `continue-on-error: true` and only fires if `REGISTRY_PAT` is set. A missing secret will not fail the release itself — the zip and GitHub Release are still created. The registry update simply won't happen automatically; you can always trigger `register-remote-plugin.yml` manually as a fallback.

### What Gets Packaged

Each repo's `release.yml` zips its plugin package directory (e.g. `cif_viewer/`) excluding `__pycache__` and `.git`. The zip is named `<prefix>_<version>.zip` (e.g. `cif_viewer_0.11.0.zip`).

### Required Plugin Source Constants

For the release and registry update to succeed, the plugin source must define:

| Constant | Required | Purpose |
| :--- | :--- | :--- |
| `PLUGIN_VERSION` | **Yes** | Must match the version entered at release time |
| `PLUGIN_NAME` | **Yes** | Used as the registry display name |
| `PLUGIN_AUTHOR` | **Yes** | Must match the GitHub repo owner (`HiroYokoyama`) |
| `PLUGIN_DESCRIPTION` | **Yes** | Registry description |
| `PLUGIN_SUPPORTED_MOLEDITPY_VERSION` | **Yes** (visible plugins) | Auto-synced to `supported_moleditpy_version` in registry |
| `PLUGIN_TAGS` | No | Registry tags (list or comma-separated string) |
| `PLUGIN_DEPENDENCIES` | No | Required pip packages |
