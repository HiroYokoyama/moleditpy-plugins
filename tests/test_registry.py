"""
Tests for REGISTRY/plugins.json integrity.

Checks that all visible plugins have required fields,
sha256 hashes match the actual files, and versions are
consistent between the registry and plugin source files.
"""

import hashlib
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "REGISTRY" / "plugins.json"
PLUGINS_DIR = ROOT / "plugins"

REQUIRED_FIELDS = ["name", "version", "description", "downloadUrl", "sha256"]
REQUIRED_METADATA = ["PLUGIN_NAME", "PLUGIN_VERSION", "PLUGIN_AUTHOR", "PLUGIN_DESCRIPTION"]


def _load_registry():
    return json.loads(REGISTRY_PATH.read_text(encoding="utf-8-sig"))


def _visible_plugins():
    return [p for p in _load_registry() if p.get("visible", False)]


def _resolve_local_path(download_url: str) -> Path:
    """Convert ../plugins/Foo/foo.py style URL to absolute path."""
    if "/plugins/" in download_url:
        rel = download_url.split("/plugins/")[-1]
        return PLUGINS_DIR / rel
    return None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Parametrize over visible plugins
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def visible_plugins():
    return _visible_plugins()


@pytest.mark.parametrize("field", REQUIRED_FIELDS)
def test_registry_required_fields(field):
    """Every visible plugin has the required registry field."""
    missing = [
        p.get("name", "?")
        for p in _visible_plugins()
        if not p.get(field)
    ]
    assert not missing, f"Visible plugins missing '{field}': {missing}"


def test_registry_names_unique():
    """All visible plugin names are unique."""
    names = [p["name"] for p in _visible_plugins()]
    assert len(names) == len(set(names)), "Duplicate plugin names in registry"


def test_registry_sha256_match():
    """sha256 in registry matches actual file on disk."""
    mismatches = []
    for p in _visible_plugins():
        url = p.get("downloadUrl", "")
        local = _resolve_local_path(url)
        if local is None or not local.exists():
            continue
        expected = p.get("sha256", "")
        actual = _sha256(local)
        if actual != expected:
            mismatches.append(f"{p['name']}: expected {expected[:12]}... got {actual[:12]}...")
    assert not mismatches, "sha256 mismatches:\n" + "\n".join(mismatches)


def test_local_files_exist():
    """Every visible plugin's downloadUrl resolves to an existing file (skips external URLs)."""
    missing = []
    for p in _visible_plugins():
        url = p.get("downloadUrl", "")
        if url.startswith("http"):
            continue  # external / hosted elsewhere — not checked here
        local = _resolve_local_path(url)
        if local is None or not local.exists():
            missing.append(f"{p['name']}: {url}")
    assert not missing, "Missing local files:\n" + "\n".join(missing)


def test_plugin_source_metadata():
    """Each visible plugin's .py file has all required PLUGIN_* metadata constants."""
    incomplete = []
    for p in _visible_plugins():
        url = p.get("downloadUrl", "")
        local = _resolve_local_path(url)
        if local is None or not local.exists() or local.suffix != ".py":
            continue
        txt = local.read_text(encoding="utf-8", errors="ignore")
        miss = [k for k in REQUIRED_METADATA if not re.search(rf"^{k}\s*=", txt, re.MULTILINE)]
        if miss:
            incomplete.append(f"{p['name']}: missing {miss}")
    assert not incomplete, "Incomplete plugin metadata:\n" + "\n".join(incomplete)


def test_version_consistency():
    """PLUGIN_VERSION in source file matches version in registry."""
    mismatches = []
    pat = re.compile(r'PLUGIN_VERSION\s*=\s*[\"]([\d.]+)[\"]')
    for p in _visible_plugins():
        url = p.get("downloadUrl", "")
        local = _resolve_local_path(url)
        if local is None or not local.exists() or local.suffix != ".py":
            continue
        txt = local.read_text(encoding="utf-8", errors="ignore")
        m = pat.search(txt)
        if m and m.group(1) != p.get("version", ""):
            mismatches.append(f"{p['name']}: file={m.group(1)} registry={p.get('version')}")
    assert not mismatches, "Version mismatches:\n" + "\n".join(mismatches)
