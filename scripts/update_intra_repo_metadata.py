import hashlib
import json
import re
import sys
from pathlib import Path


PLUGIN_VERSION_RE = re.compile(r"^\s*PLUGIN_VERSION\s*=\s*(['\"])(?P<version>.+?)\1", re.MULTILINE)
DUUNDER_VERSION_RE = re.compile(r"^\s*__version__\s*=\s*(['\"])(?P<version>.+?)\1", re.MULTILINE)
GENERIC_VERSION_RE = re.compile(r"^\s*VERSION\s*=\s*(['\"])(?P<version>.+?)\1", re.MULTILINE)
DATE_VERSION_RE = re.compile(r"^\d{4}\.\d{2}\.\d{2}$")


def sha256_of_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_plugin_version(path: Path) -> str | None:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None

    for pattern in (PLUGIN_VERSION_RE, DUUNDER_VERSION_RE, GENERIC_VERSION_RE):
        match = pattern.search(text)
        if match:
            return match.group("version").strip()
    return None


def read_package_version(package_dir: Path) -> str | None:
    init_py = package_dir / "__init__.py"
    if init_py.exists():
        version = read_plugin_version(init_py)
        if version:
            return version

    # Fallback: scan package python files for a version constant.
    for py_path in sorted(package_dir.rglob("*.py")):
        version = read_plugin_version(py_path)
        if version:
            return version
    return None


def infer_version_from_target(target: Path) -> str | None:
    if target.suffix.lower() == ".py":
        return read_plugin_version(target)

    if target.suffix.lower() == ".zip":
        package_dir = target.parent / target.stem
        if package_dir.exists() and package_dir.is_dir():
            return read_package_version(package_dir)
    return None


def find_target_json_files(repo_root: Path) -> list[Path]:
    cli_paths = [Path(arg) for arg in sys.argv[1:]]
    if cli_paths:
        return [(repo_root / p).resolve() if not p.is_absolute() else p.resolve() for p in cli_paths]

    registry_target = repo_root / "REGISTRY" / "plugins.json"
    if registry_target.exists():
        return [registry_target]
    return []


def update_single_json(json_path: Path) -> tuple[int, int, int, list[str]]:
    data = json.loads(json_path.read_text(encoding="utf-8-sig"))
    if not isinstance(data, list):
        raise ValueError(f"{json_path} root must be a list")

    updated_sha = 0
    updated_ver = 0
    updated_date = 0
    missing = []

    for plugin in data:
        download_url = plugin.get("downloadUrl")
        if not isinstance(download_url, str) or not download_url.startswith("../plugins/"):
            continue

        target = (json_path.parent / download_url).resolve()
        if not target.exists() or not target.is_file():
            missing.append(download_url)
            continue

        new_sha = sha256_of_file(target)
        if plugin.get("sha256") != new_sha:
            plugin["sha256"] = new_sha
            updated_sha += 1

        new_version = infer_version_from_target(target)
        if new_version and plugin.get("version") != new_version:
            plugin["version"] = new_version
            updated_ver += 1
        if new_version and DATE_VERSION_RE.match(new_version):
            last_updated = new_version.replace(".", "-")
            if plugin.get("lastUpdated") != last_updated:
                plugin["lastUpdated"] = last_updated
                updated_date += 1

    json_path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return updated_sha, updated_ver, updated_date, missing


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    targets = find_target_json_files(repo_root)
    if not targets:
        print("No target plugins.json found (expected REGISTRY/plugins.json).")
        return 1

    total_sha = 0
    total_ver = 0
    total_date = 0
    total_missing = 0

    for json_path in targets:
        updated_sha, updated_ver, updated_date, missing = update_single_json(json_path)
        total_sha += updated_sha
        total_ver += updated_ver
        total_date += updated_date
        total_missing += len(missing)
        rel = json_path.relative_to(repo_root)
        print(f"[{rel}] Updated sha256: {updated_sha}")
        print(f"[{rel}] Updated version: {updated_ver}")
        print(f"[{rel}] Updated lastUpdated: {updated_date}")
        print(f"[{rel}] Missing local targets: {len(missing)}")
        for item in missing:
            print(f"  - {item}")

    print(f"Total updated sha256: {total_sha}")
    print(f"Total updated version: {total_ver}")
    print(f"Total updated lastUpdated: {total_date}")
    print(f"Total missing local targets: {total_missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
