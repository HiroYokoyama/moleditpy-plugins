#!/usr/bin/env python3
from __future__ import annotations

"""
Update metadata only for an existing remote plugin in REGISTRY/plugins.json.

This script parses a GitHub Release asset URL, downloads the release package,
extracts metadata constants from the code (e.g. PLUGIN_TAGS, PLUGIN_DESCRIPTION,
PLUGIN_DEPENDENCIES, etc.), and updates the metadata fields of an existing entry
in REGISTRY/plugins.json without requiring a version increase or changing the
registered version, downloadUrl, or sha256.

A constant the code declares always wins; CLI flags only fill fields the code
leaves undeclared:
  --tags, --description, --name, --supported-version, --supported-python,
  --supported-os, --dependencies, --optional-dependencies, --visible
"""

import argparse
import datetime
import json
import os
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

# Import helpers from register_remote_plugin
sys.path.insert(0, str(Path(__file__).resolve().parent))
from register_remote_plugin import (
    DEFAULT_OS_LIST,
    DEFAULT_PYTHON_SPEC,
    _set_after,
    canonicalize_os_list,
    extract_metadata_from_file,
    find_existing_plugin,
    parse_github_release_url,
)


def update_metadata_only(
    release_url: str,
    plugin_id: str | None = None,
    tags: str | None = None,
    description: str | None = None,
    name: str | None = None,
    supported_version: str | None = None,
    supported_python: str | None = None,
    supported_os: str | None = None,
    dependencies: str | None = None,
    optional_dependencies: str | None = None,
    visible: str | None = None,
    sync_all_from_code: bool = True,
    dry_run: bool = False,
) -> dict:
    """Updates metadata fields of an existing plugin in plugins.json from URL or CLI flags."""
    url_info = parse_github_release_url(release_url)
    owner = url_info["owner"]
    repo = url_info["repo"]
    filename = url_info["filename"]

    # Locate registry file
    repo_root = Path(__file__).resolve().parents[1]
    registry_path = repo_root / "REGISTRY" / "plugins.json"
    if not registry_path.exists():
        raise FileNotFoundError(f"Registry file not found at '{registry_path}'")

    with open(registry_path, "r", encoding="utf-8-sig") as f:
        plugins = json.load(f)

    existing_entry = find_existing_plugin(plugins, owner, repo, filename, plugin_id)
    if not existing_entry:
        raise ValueError(
            f"No existing plugin entry found in registry for {owner}/{repo} ({filename}). "
            f"Plugin must already be registered to update its metadata."
        )

    # Download release asset to tempfile to extract metadata from source code
    print(f"Downloading asset from: {release_url}")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_asset = Path(tmpdir) / filename
        req = urllib.request.Request(
            release_url,
            headers={"User-Agent": "MoleditPy-Plugin-Metadata-Updater"}
        )
        with urllib.request.urlopen(req) as resp, open(tmp_asset, "wb") as out:
            out.write(resp.read())

        code_meta = extract_metadata_from_file(tmp_asset)

    print(f"Target Plugin ID: {existing_entry.get('id')}")
    print(f"Current Registry Name: {existing_entry.get('name')}")
    print(f"Extracted Code Metadata Keys: {list(code_meta.keys())}")

    changed_fields = []

    def _split(raw):
        return [x.strip() for x in raw.split(",") if x.strip()] if raw is not None else None

    def _pick(field, code_value, input_value):
        """A value the code declares always wins; an input only fills a gap."""
        if sync_all_from_code and code_value is not None:
            if input_value is not None and input_value != code_value:
                print(f"Note: ignoring the {field} input -- the code declares {code_value!r}.")
            return code_value, "code"
        if input_value is not None:
            return input_value, "input"
        return None, None

    def _code_list(key):
        if key not in code_meta or code_meta[key] is None:
            return None
        return [str(x).strip() for x in code_meta[key] if str(x).strip()]

    def _code_str(key):
        value = code_meta.get(key)
        return value.strip() if isinstance(value, str) and value.strip() else None

    # 1. Tags (an empty PLUGIN_TAGS is treated as undeclared)
    new_tags, src = _pick("tags", _code_list("tags") or None, _split(tags))
    if src and existing_entry.get("tags") != new_tags:
        existing_entry["tags"] = new_tags
        changed_fields.append(f"tags (from {src})")

    # 2. Description
    new_desc, src = _pick("description", _code_str("description"),
                          description.strip() if description is not None else None)
    if src and existing_entry.get("description") != new_desc:
        existing_entry["description"] = new_desc
        changed_fields.append(f"description (from {src})")

    # 3. Name
    new_name, src = _pick("name", _code_str("name"), name.strip() if name is not None else None)
    if src and existing_entry.get("name") != new_name:
        existing_entry["name"] = new_name
        changed_fields.append(f"name (from {src})")

    # 4. Supported MoleditPy version
    sup_ver, src = _pick("supported_version", _code_str("supported_moleditpy_version"),
                         supported_version.strip() if supported_version else None)
    if src and existing_entry.get("supported_moleditpy_version") != sup_ver:
        if "supported_moleditpy_version" in existing_entry:
            existing_entry["supported_moleditpy_version"] = sup_ver
        else:
            _set_after(existing_entry, "visible", "supported_moleditpy_version", sup_ver)
        changed_fields.append("supported_moleditpy_version")

    # 5. Visible (registry-only: no code constant)
    if visible is not None:
        vis_bool = visible.lower() == "true"
        if existing_entry.get("visible") != vis_bool:
            existing_entry["visible"] = vis_bool
            changed_fields.append("visible")

    # 6. Dependencies (a declared empty list means "needs nothing")
    deps, src = _pick("dependencies", _code_list("dependencies"), _split(dependencies))
    if src and existing_entry.get("dependencies") != deps:
        existing_entry["dependencies"] = deps
        changed_fields.append(f"dependencies (from {src})")

    # 7. Optional dependencies
    opt_deps, src = _pick("optional_dependencies", _code_list("optional_dependencies"),
                          _split(optional_dependencies))
    if src and existing_entry.get("optional_dependencies") != opt_deps:
        if opt_deps or "optional_dependencies" in existing_entry:
            _set_after(existing_entry, "dependencies", "optional_dependencies", opt_deps)
            changed_fields.append(f"optional_dependencies (from {src})")

    # 8. Supported Python version
    sup_py, src = _pick("supported_python", _code_str("supported_python_version"),
                        supported_python.strip() if supported_python else None)
    if src and existing_entry.get("supported_python_version") != sup_py:
        existing_entry["supported_python_version"] = sup_py
        changed_fields.append("supported_python_version")

    # 9. Supported OS
    code_os = canonicalize_os_list(code_meta["supported_os"]) if code_meta.get("supported_os") else None
    input_os = canonicalize_os_list(_split(supported_os)) if supported_os else None
    new_os, src = _pick("supported_os", code_os or None, input_os or None)
    if src and existing_entry.get("supported_os") != new_os:
        existing_entry["supported_os"] = new_os
        changed_fields.append("supported_os")

    if changed_fields:
        print(f"Updated fields for '{existing_entry.get('name')}': {', '.join(changed_fields)}")
    else:
        print("No metadata changes detected between source code / inputs and registry.")

    if dry_run:
        print("[Dry Run] Skipping writing to registry file.")
        print(json.dumps(existing_entry, indent=2, ensure_ascii=False))
    else:
        with open(registry_path, "w", encoding="utf-8", newline="\n") as f:
            json.dump(plugins, f, indent=2, ensure_ascii=False)
            f.write("\n")
        print("Saved updated metadata to REGISTRY/plugins.json")

    return {
        "plugin_id": existing_entry.get("id"),
        "plugin_name": existing_entry.get("name"),
        "version": existing_entry.get("version"),
        "changed_fields": changed_fields,
    }


def main():
    parser = argparse.ArgumentParser(description="Update metadata only for an existing remote plugin in registry.")
    parser.add_argument("release_url", help="GitHub Release file URL (.zip or .py)")
    parser.add_argument("--id", dest="plugin_id", help="Plugin ID (optional, matched by repository or filename if omitted)")
    parser.add_argument("--tags", help="Comma-separated tags, used only when the code declares no PLUGIN_TAGS.")
    parser.add_argument("--description", help="Description, used only when the code declares no PLUGIN_DESCRIPTION.")
    parser.add_argument("--name", help="Plugin name, used only when the code declares no PLUGIN_NAME.")
    parser.add_argument("--supported-version", dest="supported_version", help="Supported MoleditPy version spec")
    parser.add_argument("--supported-python", dest="supported_python", help="Supported Python version spec")
    parser.add_argument("--supported-os", dest="supported_os", help="Comma-separated list of supported OS tokens")
    parser.add_argument("--dependencies", help="Comma-separated list of dependencies")
    parser.add_argument("--optional-dependencies", dest="optional_dependencies", help="Comma-separated list of optional dependencies")
    parser.add_argument("--visible", help="Visibility (true or false)")
    parser.add_argument("--dry-run", action="store_true", help="Validate and display changes without modifying registry")

    args = parser.parse_args()

    res = update_metadata_only(
        release_url=args.release_url,
        plugin_id=args.plugin_id,
        tags=args.tags,
        description=args.description,
        name=args.name,
        supported_version=args.supported_version,
        supported_python=args.supported_python,
        supported_os=args.supported_os,
        dependencies=args.dependencies,
        optional_dependencies=args.optional_dependencies,
        visible=args.visible,
        dry_run=args.dry_run,
    )

    if "GITHUB_OUTPUT" in os.environ:
        with open(os.environ["GITHUB_OUTPUT"], "a") as f:
            f.write(f"plugin_id={res['plugin_id']}\n")
            f.write(f"plugin_name={res['plugin_name']}\n")
            f.write(f"version={res['version']}\n")
            f.write(f"changed_fields={','.join(res['changed_fields'])}\n")


if __name__ == "__main__":
    main()
