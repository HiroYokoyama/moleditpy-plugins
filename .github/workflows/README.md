# GitHub Actions Workflows Documentation

This directory contains the automated workflows for the MoleditPy Plugins Registry.

---

## 1. Register or Update Remote Plugin (`register-remote-plugin.yml`)

A manual workflow (`workflow_dispatch`) used to automatically register new third-party/remote plugins or update existing remote plugins in the registry (`REGISTRY/plugins.json`).

### Workflow Inputs

| Input Field | Required? | Description |
| :--- | :--- | :--- |
| `release_url` | **Yes** | The direct download URL of the GitHub Release asset (must be a `.py` file or a `.zip` file). Format: `https://github.com/{owner}/{repo}/releases/download/{tag}/{filename}`. |
| `plugin_id` | No | The unique ID for the plugin. For new plugins, if omitted, it will be automatically derived from the repository name. |
| `tags` | No | A comma-separated list of tags (only used when registering a new plugin, e.g., `Analysis, Visualization`). |
| `dependencies` | No | A comma-separated list of required Python packages (only used when registering a new plugin, e.g., `numpy, rdkit`). |
| `visible` | **Yes** | Visibility flag in the registry (`true` or `false`). Defaults to `true`. |
| `expected_sha256` | *Conditional* | The expected SHA-256 hash. **Mandatory** for security verification if the repository owner is not `HiroYokoyama`. |
| `dry_run` | **Yes** | If set to `true`, the workflow performs all downloads and verification checks but **does not** commit or push changes back to the registry. |

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

## 2. Plugin Tests (`test-plugins.yml`)

An automated CI workflow triggered on every `push` and `pull_request` to `main`. It runs:
1. Registry JSON validation (`scripts/validate_json.py`) to check formatting and prevent duplicate IDs or SHA-256 hashes.
2. Full automated pytest execution on `tests/` (including validation of local files, source metadata conventions, and remote registration script unit tests).
3. API compatibility checks against the main `python_molecular_editor` repository.
