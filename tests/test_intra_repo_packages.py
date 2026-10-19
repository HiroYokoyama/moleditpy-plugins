"""Check what users download, independently of the unpacked source metadata."""
import ast
import json
from pathlib import Path
import subprocess
from zipfile import ZipFile

import pytest

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = [entry for entry in json.loads((ROOT / "REGISTRY/plugins.json").read_text(encoding="utf-8"))
            if entry.get("visible", True) and entry["downloadUrl"].startswith("../plugins/")
            and entry["downloadUrl"].endswith(".zip")]


@pytest.mark.parametrize("entry", PACKAGES, ids=lambda entry: entry["id"])
def test_downloaded_zip_contains_current_sources_and_metadata(entry):
    target = (ROOT / "REGISTRY" / entry["downloadUrl"]).resolve()
    package = target.parent / target.stem
    with ZipFile(target) as archive:
        files = subprocess.check_output(
            ["git", "ls-files", "-z", "--", str(package.relative_to(ROOT))], cwd=ROOT,
        ).decode("utf-8").split("\0")
        for name in filter(None, files):
            path = ROOT / name
            member = path.relative_to(package.parent).as_posix()
            assert archive.read(member) == path.read_bytes(), f"Stale packaged source: {member}"
        metadata = {}
        for node in ast.parse(archive.read(f"{package.name}/__init__.py").decode("utf-8")).body:
            if isinstance(node, ast.Assign):
                for field in node.targets:
                    if isinstance(field, ast.Name) and field.id in ("PLUGIN_VERSION", "PLUGIN_SUPPORTED_MOLEDITPY_VERSION"):
                        metadata[field.id] = ast.literal_eval(node.value)
        assert metadata["PLUGIN_VERSION"] == entry["version"]
        assert metadata["PLUGIN_SUPPORTED_MOLEDITPY_VERSION"] == entry["supported_moleditpy_version"]
        for name in ("LICENSE", "COPYING", "LICENSING.md"):
            assert archive.read(f"{package.name}/{name}") == (ROOT / name).read_bytes()
