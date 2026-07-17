"""
Update sha256, version, and lastUpdated fields in REGISTRY/plugins.json.

Reads PLUGIN_VERSION (or __version__) from each plugin's source file,
computes sha256 of the download target (.py or .zip), and writes back
to the registry. Run after modifying any plugin file.

Usage:
  python scripts/update_intra_repo_metadata.py
"""
import hashlib
import json
import re
import sys
from pathlib import Path


PLUGIN_VERSION_RE = re.compile(r"^\s*PLUGIN_VERSION\s*=\s*(['\"])(?P<version>.+?)\1", re.MULTILINE)
DUUNDER_VERSION_RE = re.compile(r"^\s*__version__\s*=\s*(['\"])(?P<version>.+?)\1", re.MULTILINE)
GENERIC_VERSION_RE = re.compile(r"^\s*VERSION\s*=\s*(['\"])(?P<version>.+?)\1", re.MULTILINE)
DATE_VERSION_RE = re.compile(r"^\d{4}\.\d{2}\.\d{2}$")
SUPPORTED_VER_RE = re.compile(r"^\s*PLUGIN_SUPPORTED_MOLEDITPY_VERSION\s*=\s*(['\"])(?P<ver>.+?)\1", re.MULTILINE)
SUPPORTED_PY_RE = re.compile(r"^\s*PLUGIN_SUPPORTED_PYTHON_VERSION\s*=\s*(['\"])(?P<ver>.+?)\1", re.MULTILINE)

# Applied to visible plugins whose source does not declare
# PLUGIN_SUPPORTED_PYTHON_VERSION.
DEFAULT_PYTHON_SPEC = ">=3.9, <3.15"


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


def _read_constant(path: Path, pattern: re.Pattern) -> str | None:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    match = pattern.search(text)
    return match.group("ver").strip() if match else None


def read_supported_version(path: Path) -> str | None:
    return _read_constant(path, SUPPORTED_VER_RE)


def _infer_constant_from_target(target: Path, pattern: re.Pattern) -> str | None:
    if target.suffix.lower() == ".py":
        return _read_constant(target, pattern)
    if target.suffix.lower() == ".zip":
        package_dir = target.parent / target.stem
        if package_dir.exists() and package_dir.is_dir():
            init_py = package_dir / "__init__.py"
            if init_py.exists():
                ver = _read_constant(init_py, pattern)
                if ver:
                    return ver
            for py_path in sorted(package_dir.rglob("*.py")):
                ver = _read_constant(py_path, pattern)
                if ver:
                    return ver
    return None


def infer_supported_version_from_target(target: Path) -> str | None:
    return _infer_constant_from_target(target, SUPPORTED_VER_RE)


def infer_supported_python_from_target(target: Path) -> str | None:
    return _infer_constant_from_target(target, SUPPORTED_PY_RE)


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


def update_single_json(json_path: Path) -> tuple[int, int, int, int, int, list[str]]:
    data = json.loads(json_path.read_text(encoding="utf-8-sig"))
    if not isinstance(data, list):
        raise ValueError(f"{json_path} root must be a list")

    updated_sha = 0
    updated_ver = 0
    updated_date = 0
    updated_supported = 0
    updated_supported_py = 0
    missing = []

    # Every visible entry (external ones included) gets a python spec; the
    # per-target pass below overrides it from PLUGIN_SUPPORTED_PYTHON_VERSION.
    for plugin in data:
        if plugin.get("visible", True) and not plugin.get("supported_python_version"):
            plugin["supported_python_version"] = DEFAULT_PYTHON_SPEC
            updated_supported_py += 1

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

        new_supported = infer_supported_version_from_target(target)
        if new_supported and plugin.get("supported_moleditpy_version") != new_supported:
            plugin["supported_moleditpy_version"] = new_supported
            updated_supported += 1

        new_supported_py = infer_supported_python_from_target(target)
        if new_supported_py and plugin.get("supported_python_version") != new_supported_py:
            plugin["supported_python_version"] = new_supported_py
            updated_supported_py += 1

    with json_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    return updated_sha, updated_ver, updated_date, updated_supported, updated_supported_py, missing


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    targets = find_target_json_files(repo_root)
    if not targets:
        print("No target plugins.json found (expected REGISTRY/plugins.json).")
        return 1

    total_sha = 0
    total_ver = 0
    total_date = 0
    total_supported = 0
    total_supported_py = 0
    total_missing = 0

    for json_path in targets:
        updated_sha, updated_ver, updated_date, updated_supported, updated_supported_py, missing = update_single_json(json_path)
        total_sha += updated_sha
        total_ver += updated_ver
        total_date += updated_date
        total_supported += updated_supported
        total_supported_py += updated_supported_py
        total_missing += len(missing)
        rel = json_path.relative_to(repo_root)
        print(f"[{rel}] Updated sha256: {updated_sha}")
        print(f"[{rel}] Updated version: {updated_ver}")
        print(f"[{rel}] Updated lastUpdated: {updated_date}")
        print(f"[{rel}] Updated supported_moleditpy_version: {updated_supported}")
        print(f"[{rel}] Updated supported_python_version: {updated_supported_py}")
        print(f"[{rel}] Missing local targets: {len(missing)}")
        for item in missing:
            print(f"  - {item}")

    print(f"Total updated sha256: {total_sha}")
    print(f"Total updated version: {total_ver}")
    print(f"Total updated lastUpdated: {total_date}")
    print(f"Total updated supported_moleditpy_version: {total_supported}")
    print(f"Total updated supported_python_version: {total_supported_py}")
    print(f"Total missing local targets: {total_missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
