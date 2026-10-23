"""Set collection ZIP registry entries from the packages built for one tag."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path, PurePosixPath


def update_registry(registry_path, manifest_path, repository, tag):
    if not re.fullmatch(r"v\d{4}\.\d{2}\.\d{2}", tag):
        raise ValueError("Expected a vYYYY.MM.DD collection tag")
    if not re.fullmatch(r"[\w.-]+/[\w.-]+", repository):
        raise ValueError("Invalid GitHub repository")
    registry_path = Path(registry_path)
    entries = json.loads(registry_path.read_text(encoding="utf-8"))
    packages = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    by_id = {entry["id"]: entry for entry in entries}
    seen = set()
    for package in packages:
        entry = by_id[package["id"]]
        source = package["sourcePath"]
        asset = package["assetName"]
        digest = package["sha256"]
        if package["id"] in seen:
            raise ValueError("Duplicate package ID")
        seen.add(package["id"])
        if (
            not source.startswith("../plugins/")
            or ".." in PurePosixPath(source).parts[1:]
            or "/_old/" in source
            or "\\" in source
        ):
            raise ValueError("Invalid package source path")
        if not re.fullmatch(r"[\w.-]+\.zip", asset) or not re.fullmatch(
            r"[0-9a-f]{64}", digest
        ):
            raise ValueError("Invalid package asset metadata")
        if not entry.get("visible", True):
            raise ValueError("Retired plugins cannot be republished")
        metadata = package["metadata"]
        allowed = {
            "name",
            "author",
            "description",
            "tags",
            "dependencies",
            "optional_dependencies",
            "supported_os",
            "license",
            "version",
            "supported_moleditpy_version",
            "supported_python_version",
        }
        if set(metadata) - allowed:
            raise ValueError("Unexpected package metadata fields")
        current = entry.get("version", "")
        version = metadata["version"]
        if (
            re.fullmatch(r"\d{4}\.\d{2}\.\d{2}", current)
            and re.fullmatch(r"\d{4}\.\d{2}\.\d{2}", version)
            and current > version
        ):
            continue  # A manual backfill must never downgrade newer metadata.
        prefix = f"https://github.com/{repository}/releases/download/"
        url = entry.get("downloadUrl", "")
        if not url.startswith(("../plugins/", prefix)):
            raise ValueError("The registry entry is no longer a collection package")
        if url.startswith(prefix) and url[len(prefix) :].split("/")[0] > tag:
            continue
        entry.update(metadata)
        entry["sourcePath"] = source
        entry["downloadUrl"] = (
            f"https://github.com/{repository}/releases/download/{tag}/{asset}"
        )
        entry["sha256"] = digest
        version = entry["version"]
        if re.fullmatch(r"\d{4}\.\d{2}\.\d{2}", version):
            entry["lastUpdated"] = version.replace(".", "-")
    registry_path.write_text(
        json.dumps(entries, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--repository", required=True)
    parser.add_argument("--tag", required=True)
    args = parser.parse_args()
    update_registry(args.registry, args.manifest, args.repository, args.tag)
