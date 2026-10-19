"""Build deterministic registry ZIPs from tracked plugin sources and notices."""
import json
from pathlib import Path
import subprocess
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

ROOT = Path(__file__).resolve().parents[1]


def build_packages(root=ROOT):
    entries = json.loads((root / "REGISTRY/plugins.json").read_text(encoding="utf-8"))
    for entry in entries:
        url = entry["downloadUrl"]
        if not entry.get("visible", True) or not url.startswith("../plugins/") or not url.endswith(".zip"):
            continue
        target = (root / "REGISTRY" / url).resolve()
        package = target.parent / target.stem
        tracked = subprocess.check_output(
            ["git", "ls-files", "-z", "--", str(package.relative_to(root))],
            cwd=root,
        ).decode("utf-8").split("\0")
        sources = [root / name for name in tracked if name]
        if not sources or not (package / "__init__.py").is_file():
            raise ValueError(f"No tracked package sources for {target}")
        contents = {str(path.relative_to(package.parent)).replace("\\", "/"): path.read_bytes() for path in sources}
        for name in ("LICENSE", "COPYING", "LICENSING.md"):
            contents[f"{package.name}/{name}"] = (root / name).read_bytes()
        with ZipFile(target, "w", compression=ZIP_DEFLATED) as archive:
            for name, data in sorted(contents.items()):
                info = ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.create_system = 3
                info.external_attr = 0o100644 << 16
                archive.writestr(info, data, compress_type=ZIP_DEFLATED)
        print(target.relative_to(root))


if __name__ == "__main__":
    build_packages()
