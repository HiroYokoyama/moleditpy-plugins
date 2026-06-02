#!/usr/bin/env python3
"""
Register or update a remote plugin in REGISTRY/plugins.json.

This script parses a GitHub Release URL, downloads the corresponding .py or .zip asset,
performs strict validation on the version tag and internal code version constants,
and automatically updates or registers the plugin entry in REGISTRY/plugins.json.

Usage:
    python scripts/register_remote_plugin.py <release_url> [--id <plugin_id>] [--tags <tags>] [--dependencies <dependencies>] [--visible <true|false>]
"""

import argparse
import datetime
import hashlib
import json
import os
import re
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

# Match metadata constants in python files
PLUGIN_NAME_RE = re.compile(r"^\s*PLUGIN_NAME\s*=\s*(['\"])(?P<val>.*?)\1", re.MULTILINE)
PLUGIN_VERSION_RE = re.compile(r"^\s*PLUGIN_VERSION\s*=\s*(['\"])(?P<val>.*?)\1", re.MULTILINE)
PLUGIN_AUTHOR_RE = re.compile(r"^\s*PLUGIN_AUTHOR\s*=\s*(['\"])(?P<val>.*?)\1", re.MULTILINE)
PLUGIN_DESCRIPTION_RE = re.compile(r"^\s*PLUGIN_DESCRIPTION\s*=\s*(['\"])(?P<val>.*?)\1", re.MULTILINE)

def parse_github_release_url(url: str) -> dict:
    """
    Parses GitHub release asset URL.
    Expected format: https://github.com/{owner}/{repo}/releases/download/{tag}/{filename}
    """
    pattern = r"^https://github\.com/(?P<owner>[^/]+)/(?P<repo>[^/]+)/releases/download/(?P<tag>[^/]+)/(?P<filename>[^/]+)$"
    match = re.match(pattern, url.strip())
    if not match:
        raise ValueError(
            f"Invalid GitHub Release asset URL format: '{url}'.\n"
            f"Expected: https://github.com/{{owner}}/{{repo}}/releases/download/{{tag}}/{{filename}}"
        )
    return match.groupdict()

def parse_version(v_str: str) -> tuple:
    """Parses a version string into a tuple for semantic comparison."""
    v_str = v_str.strip().lstrip('vV')
    parts = []
    for part in v_str.split('.'):
        try:
            parts.append(int(part))
        except ValueError:
            parts.append(part)
    return tuple(parts)

def extract_metadata_from_code(code_text: str) -> dict:
    """Extracts plugin metadata constants from python code."""
    meta = {}
    for key, regex in [
        ("name", PLUGIN_NAME_RE),
        ("version", PLUGIN_VERSION_RE),
        ("author", PLUGIN_AUTHOR_RE),
        ("description", PLUGIN_DESCRIPTION_RE),
    ]:
        match = regex.search(code_text)
        if match:
            meta[key] = match.group("val").strip()
    return meta

def extract_metadata_from_file(file_path: Path) -> dict:
    """Extracts metadata from a .py or .zip file."""
    suffix = file_path.suffix.lower()
    
    if suffix == ".py":
        try:
            content = file_path.read_text(encoding="utf-8", errors="ignore")
            return extract_metadata_from_code(content)
        except Exception as e:
            raise ValueError(f"Failed to read python source file metadata: {e}")
            
    elif suffix == ".zip":
        try:
            with zipfile.ZipFile(file_path, 'r') as z:
                # Find all .py files inside zip
                py_files = [name for name in z.namelist() if name.endswith(".py")]
                
                # Check for __init__.py files first
                init_files = [name for name in py_files if name.endswith("__init__.py")]
                search_order = init_files + [name for name in py_files if not name.endswith("__init__.py")]
                
                for py_file in search_order:
                    with z.open(py_file) as f:
                        content = f.read().decode("utf-8", errors="ignore")
                        meta = extract_metadata_from_code(content)
                        if "version" in meta:
                            print(f"Extracted metadata from {py_file} inside zip.")
                            return meta
        except Exception as e:
            raise ValueError(f"Failed to read metadata from zip package: {e}")
            
    raise ValueError(f"Unsupported file format '{suffix}'. Only .py and .zip are supported.")

def find_existing_plugin(plugins: list, owner: str, repo: str, filename: str, plugin_id: str = None) -> dict:
    """Searches for an existing plugin entry in registry."""
    if plugin_id:
        for p in plugins:
            if p.get("id") == plugin_id:
                return p

    target_repo = f"github.com/{owner}/{repo}".lower()
    for p in plugins:
        dl_url = p.get("downloadUrl", "")
        if dl_url.startswith("http") and target_repo in dl_url.lower():
            # Match the file name or file name stem (excluding version tags in zip names if necessary)
            existing_filename = dl_url.split("/")[-1]
            existing_stem = Path(existing_filename).stem
            new_stem = Path(filename).stem
            if existing_filename == filename or existing_stem == new_stem:
                return p
                
        project_url = p.get("projectUrl", "").lower().rstrip(".git").rstrip("/")
        target_project = f"https://github.com/{owner}/{repo}".lower().rstrip(".git").rstrip("/")
        if project_url == target_project:
            return p
            
    return None

def main():
    parser = argparse.ArgumentParser(description="Register or update a remote plugin in the registry.")
    parser.add_argument("release_url", help="GitHub Release file URL")
    parser.add_argument("--id", dest="plugin_id", help="Plugin ID (optional, derived for new plugins if omitted)")
    parser.add_argument("--tags", help="Comma-separated tags for new plugins")
    parser.add_argument("--dependencies", help="Comma-separated dependencies for new plugins")
    parser.add_argument("--visible", default="true", help="Set visibility of the plugin (default: true)")
    parser.add_argument("--dry-run", action="store_true", help="Perform checks and downloads but do not write to the registry")
    
    args = parser.parse_args()
    
    # 1. Parse and validate GitHub Release URL
    try:
        url_info = parse_github_release_url(args.release_url)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
        
    owner = url_info["owner"]
    repo = url_info["repo"]
    tag = url_info["tag"]
    filename = url_info["filename"]
    
    print(f"GitHub Release Details:")
    print(f"  Owner:    {owner}")
    print(f"  Repo:     {repo}")
    print(f"  Tag:      {tag}")
    print(f"  Filename: {filename}")
    if args.dry_run:
        print("  Running in DRY-RUN mode.")
    
    # Validate tag format (X.X.X or vX.X.X)
    tag_match = re.match(r"^v?(?P<ver>\d+\.\d+\.\d+)$", tag)
    if not tag_match:
        print(f"Error: Release tag '{tag}' does not match expected format (X.X.X or vX.X.X).", file=sys.stderr)
        sys.exit(1)
    tag_version = tag_match.group("ver")
    print(f"  Semantic Tag Version: {tag_version}")
    
    # 2. Locate registry file
    repo_root = Path(__file__).resolve().parents[1]
    registry_path = repo_root / "REGISTRY" / "plugins.json"
    if not registry_path.exists():
        print(f"Error: Registry file not found at '{registry_path}'", file=sys.stderr)
        sys.exit(1)
        
    with open(registry_path, "r", encoding="utf-8-sig") as f:
        plugins = json.load(f)
        
    # 3. Detect Mode (Update vs Add)
    existing_entry = find_existing_plugin(plugins, owner, repo, filename, args.plugin_id)
    mode = "UPDATE" if existing_entry else "ADD"
    print(f"Detected execution mode: {mode}")
    
    # If adding, validate plugin ID uniqueness
    plugin_id = args.plugin_id
    if mode == "ADD":
        if not plugin_id:
            # Derive ID from repo name (strip prefix/suffix)
            plugin_id = repo.lower().replace("moleditpy_", "").replace("_plugin", "").replace("-", "_")
            print(f"Derived plugin ID: '{plugin_id}'")
            
        # Ensure unique ID
        for p in plugins:
            if p.get("id") == plugin_id:
                print(f"Error: Plugin ID '{plugin_id}' already exists in registry.", file=sys.stderr)
                sys.exit(1)
    else:
        plugin_id = existing_entry["id"]
        print(f"Matching existing entry: ID='{plugin_id}', Name='{existing_entry['name']}'")
        
    # 4. Download Release Asset
    print(f"Downloading release asset from '{args.release_url}'...")
    try:
        with urllib.request.urlopen(args.release_url) as response:
            asset_data = response.read()
    except Exception as e:
        print(f"Error downloading asset: {e}", file=sys.stderr)
        sys.exit(1)
        
    sha256_hash = hashlib.sha256(asset_data).hexdigest()
    print(f"Downloaded successfully. SHA-256: {sha256_hash}")
    
    # Save to temp file to read metadata
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = Path(temp_dir) / filename
        temp_file_path.write_bytes(asset_data)
        
        # 5. Extract and Validate Metadata
        try:
            meta = extract_metadata_from_file(temp_file_path)
        except ValueError as e:
            print(f"Error validating downloaded asset: {e}", file=sys.stderr)
            sys.exit(1)
            
    code_version = meta.get("version")
    if not code_version:
        print("Error: Strict Validation Failed: PLUGIN_VERSION constant is missing from the plugin source code.", file=sys.stderr)
        sys.exit(1)
        
    # Verify version matches tag
    normalized_code_version = code_version.strip().lstrip('vV')
    if normalized_code_version != tag_version:
        print(
            f"Error: Strict Consistency Failed: Version in code ('{code_version}') "
            f"does not match tag version in Release URL ('{tag}')",
            file=sys.stderr
        )
        sys.exit(1)
        
    print(f"Metadata extracted successfully:")
    print(f"  Name:        {meta.get('name')}")
    print(f"  Version:     {code_version}")
    print(f"  Author:      {meta.get('author')}")
    print(f"  Description: {meta.get('description')}")
    
    # 6. Apply Registry updates
    today_str = datetime.date.today().isoformat()
    
    if mode == "UPDATE":
        # Check that version increases
        old_version = existing_entry.get("version", "0.0.0")
        try:
            if parse_version(normalized_code_version) <= parse_version(old_version):
                print(
                    f"Error: Version Check Failed: New version '{code_version}' "
                    f"is not greater than existing version '{old_version}' in registry.",
                    file=sys.stderr
                )
                sys.exit(1)
        except SystemExit:
            raise
        except Exception as e:
            print(f"Warning: Error comparing versions ({e}), defaulting current version to 0.0.0 for safety.")
            old_version = "0.0.0"
            
        existing_entry["version"] = normalized_code_version
        existing_entry["downloadUrl"] = args.release_url
        existing_entry["sha256"] = sha256_hash
        existing_entry["lastUpdated"] = today_str
        plugin_name = existing_entry["name"]
        
    else: # mode == "ADD"
        # Validate other required fields are present in the code metadata
        missing = [k for k in ["name", "author", "description"] if not meta.get(k)]
        if missing:
            print(f"Error: Strict Validation Failed: Missing required metadata constants: {', '.join([f'PLUGIN_{k.upper()}' for k in missing])}", file=sys.stderr)
            sys.exit(1)
            
        # Parse inputs
        tags_list = [t.strip() for t in args.tags.split(",")] if args.tags else []
        deps_list = [d.strip() for d in args.dependencies.split(",")] if args.dependencies else []
        visible_flag = args.visible.lower() == "true"
        
        new_entry = {
            "id": plugin_id,
            "visible": visible_flag,
            "name": meta["name"],
            "version": normalized_code_version,
            "author": meta["author"],
            "authorUrl": f"https://github.com/{owner}",
            "description": meta["description"],
            "tags": tags_list,
            "dependencies": deps_list,
            "downloadUrl": args.release_url,
            "projectUrl": f"https://github.com/{owner}/{repo}",
            "lastUpdated": today_str,
            "sha256": sha256_hash,
            "firstAppeared": today_str
        }
        plugins.append(new_entry)
        plugin_name = meta["name"]
        
    if args.dry_run:
        print("[Dry Run] Verification complete. Skipping writing modifications to registry.")
    else:
        # Write back to JSON
        with open(registry_path, "w", encoding="utf-8") as f:
            json.dump(plugins, f, indent=2, ensure_ascii=False)
            f.write("\n")
        print(f"Successfully saved changes to REGISTRY/plugins.json!")
    
    # 7. Write outputs to GitHub Actions environment if present
    if "GITHUB_OUTPUT" in os.environ:
        with open(os.environ["GITHUB_OUTPUT"], "a") as f:
            f.write(f"plugin_id={plugin_id}\n")
            f.write(f"plugin_name={plugin_name}\n")
            f.write(f"version={normalized_code_version}\n")
            f.write(f"mode={mode}\n")
            
    print("Execution completed successfully.")

if __name__ == "__main__":
    main()
