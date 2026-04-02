import hashlib
import json
import re
from pathlib import Path


PLUGIN_VERSION_RE = re.compile(
    r"^\s*PLUGIN_VERSION\s*=\s*(['\"])(?P<version>.+?)\1",
    re.MULTILINE,
)


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
    match = PLUGIN_VERSION_RE.search(text)
    return match.group("version").strip() if match else None


def infer_version_from_target(target: Path) -> str | None:
    if target.suffix.lower() == ".py":
        return read_plugin_version(target)

    if target.suffix.lower() == ".zip":
        init_py = target.parent / target.stem / "__init__.py"
        if init_py.exists():
            return read_plugin_version(init_py)
    return None


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    json_path = repo_root / "explorer" / "plugins.json"

    data = json.loads(json_path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("plugins.json root must be a list")

    updated_sha = 0
    updated_ver = 0
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

    json_path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"Updated sha256: {updated_sha}")
    print(f"Updated version: {updated_ver}")
    print(f"Missing local targets: {len(missing)}")
    for item in missing:
        print(f"  - {item}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
