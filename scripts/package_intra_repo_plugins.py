"""Build deterministic registry ZIPs from tracked plugin sources and notices."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

ROOT = Path(__file__).resolve().parents[1]


def build_packages(root=ROOT, output_dir=None):
    from update_intra_repo_metadata import (
        infer_display_metadata_from_target,
        infer_optional_dependencies_from_target,
        infer_supported_python_from_target,
        infer_supported_version_from_target,
        infer_version_from_target,
    )

    root = Path(root).resolve()
    output_dir = Path(output_dir).resolve() if output_dir is not None else None
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    entries = json.loads((root / "REGISTRY/plugins.json").read_text(encoding="utf-8"))
    for entry in entries:
        url = entry.get("sourcePath", entry["downloadUrl"])
        if (
            not entry.get("visible", True)
            or not url.startswith("../plugins/")
            or not url.endswith(".zip")
        ):
            continue
        target = (root / "REGISTRY" / url).resolve()
        relative = target.relative_to(root / "plugins")
        if "_old" in relative.parts:
            continue
        package = target.parent / target.stem
        tracked = (
            subprocess.check_output(
                ["git", "ls-files", "-z", "--", str(package.relative_to(root))],
                cwd=root,
            )
            .decode("utf-8")
            .split("\0")
        )
        sources = [root / name for name in tracked if name]
        if not sources or not (package / "__init__.py").is_file():
            raise ValueError(f"No tracked package sources for {target}")
        contents = {
            str(path.relative_to(package.parent)).replace("\\", "/"): path.read_bytes()
            for path in sources
        }
        for name in ("LICENSE", "COPYING", "LICENSING.md"):
            contents[f"{package.name}/{name}"] = (root / name).read_bytes()
        destination = output_dir / target.name if output_dir is not None else target
        if any(item["assetName"] == destination.name for item in manifest):
            raise ValueError(f"Duplicate release asset name: {destination.name}")
        with ZipFile(destination, "w", compression=ZIP_DEFLATED) as archive:
            for name, data in sorted(contents.items()):
                info = ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.create_system = 3
                info.external_attr = 0o100644 << 16
                archive.writestr(info, data, compress_type=ZIP_DEFLATED)
        metadata = infer_display_metadata_from_target(target)
        for key, value in {
            "version": infer_version_from_target(target),
            "supported_moleditpy_version": infer_supported_version_from_target(target),
            "supported_python_version": infer_supported_python_from_target(target),
            "optional_dependencies": infer_optional_dependencies_from_target(target),
        }.items():
            if value is not None:
                metadata[key] = value
        if "version" not in metadata:
            raise ValueError(f"No source version for {entry['id']}")
        manifest.append(
            {
                "id": entry["id"],
                "sourcePath": url,
                "assetName": destination.name,
                "sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
                "metadata": metadata,
            }
        )
        print(destination)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    packages = build_packages(output_dir=args.output_dir)
    if args.manifest:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(
            json.dumps(packages, indent=2) + "\n", encoding="utf-8"
        )
