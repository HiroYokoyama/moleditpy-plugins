"""Verify published ZIPs match tracked sources and registry publication is safe."""

import hashlib
import importlib
import json
import subprocess
import zipfile
from pathlib import Path

import pytest


@pytest.fixture
def scripts(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    return (
        importlib.import_module("package_intra_repo_plugins"),
        importlib.import_module("update_release_packages"),
    )


@pytest.fixture
def package_repo(tmp_path):
    root = tmp_path / "repo"
    package = root / "plugins" / "Example" / "example"
    package.mkdir(parents=True)
    (root / "REGISTRY").mkdir()
    (package / "__init__.py").write_text(
        'PLUGIN_NAME = "Example"\nPLUGIN_VERSION = "2026.10.23"\n'
        'PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=5.0.0, <6.0.0"\n',
        encoding="utf-8",
    )
    (package / "engine.py").write_text("VALUE = 42\n", encoding="utf-8")
    (package / "settings.json").write_text('{"private": true}', encoding="utf-8")
    for name in ("LICENSE", "COPYING", "LICENSING.md"):
        (root / name).write_text("Public license notice\n", encoding="utf-8")
    legacy = package.with_suffix(".zip")
    legacy.write_bytes(b"frozen legacy download")
    entry = {
        "id": "example",
        "visible": True,
        "version": "2026.10.20",
        "downloadUrl": "../plugins/Example/example.zip",
        "sha256": "0" * 64,
    }
    (root / "REGISTRY/plugins.json").write_text(json.dumps([entry]), encoding="utf-8")
    subprocess.run(["git", "init", "--quiet", str(root)], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "add",
            "--",
            "plugins/Example/example/__init__.py",
            "plugins/Example/example/engine.py",
        ],
        check=True,
    )
    return root


def test_zip_is_deterministic_and_excludes_untracked_private_files(
    scripts, package_repo, tmp_path
):
    packager, _ = scripts
    destination = tmp_path / "assets"
    first = packager.build_packages(package_repo, destination)
    archive = destination / "example.zip"
    raw = archive.read_bytes()
    assert packager.build_packages(package_repo, destination) == first
    assert archive.read_bytes() == raw
    assert first[0]["sha256"] == hashlib.sha256(raw).hexdigest()
    assert first[0]["metadata"]["version"] == "2026.10.23"
    with zipfile.ZipFile(archive) as zipped:
        assert set(zipped.namelist()) == {
            "example/__init__.py",
            "example/engine.py",
            "example/LICENSE",
            "example/COPYING",
            "example/LICENSING.md",
        }
        assert (
            zipped.read("example/engine.py")
            == (package_repo / "plugins/Example/example/engine.py").read_bytes()
        )
    assert (
        package_repo / "plugins/Example/example.zip"
    ).read_bytes() == b"frozen legacy download"


def test_future_build_finds_sources_after_url_migration(
    scripts, package_repo, tmp_path
):
    packager, _ = scripts
    registry = package_repo / "REGISTRY/plugins.json"
    entries = json.loads(registry.read_text())
    entries[0]["sourcePath"] = entries[0]["downloadUrl"]
    entries[0]["downloadUrl"] = (
        "https://github.com/owner/repo/releases/download/v2026.10.22/example.zip"
    )
    registry.write_text(json.dumps(entries), encoding="utf-8")
    assert (
        packager.build_packages(package_repo, tmp_path / "assets")[0]["id"] == "example"
    )


def _publish(scripts, root, tmp_path, mutate=None):
    packager, publisher = scripts
    packages = packager.build_packages(root, tmp_path / "assets")
    if mutate:
        mutate(packages[0])
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(packages), encoding="utf-8")
    registry = root / "REGISTRY/plugins.json"
    publisher.update_registry(registry, manifest, "owner/repo", "v2026.10.23")
    return json.loads(registry.read_text(encoding="utf-8"))


def test_publishing_preserves_other_plugin_registrations(
    scripts, package_repo, tmp_path
):
    registry = package_repo / "REGISTRY/plugins.json"
    entries = json.loads(registry.read_text())
    external = {
        "id": "external",
        "version": "9.8.7",
        "downloadUrl": "https://example.org/released.zip",
        "sha256": "f" * 64,
    }
    registry.write_text(json.dumps(entries + [external]), encoding="utf-8")
    result = _publish(scripts, package_repo, tmp_path)
    assert result[1] == external
    assert (
        result[0]["downloadUrl"]
        == "https://github.com/owner/repo/releases/download/v2026.10.23/example.zip"
    )
    assert result[0]["sourcePath"] == "../plugins/Example/example.zip"
    assert result[0]["version"] == "2026.10.23"
    assert result[0]["lastUpdated"] == "2026-10-23"
    assert (
        result[0]["sha256"]
        == hashlib.sha256((tmp_path / "assets/example.zip").read_bytes()).hexdigest()
    )


def test_backfill_cannot_replace_newer_published_version(
    scripts, package_repo, tmp_path
):
    registry = package_repo / "REGISTRY/plugins.json"
    entries = json.loads(registry.read_text())
    entries[0]["version"] = "2026.10.24"
    registry.write_text(json.dumps(entries), encoding="utf-8")
    assert _publish(scripts, package_repo, tmp_path) == entries


@pytest.mark.parametrize(
    "key,value",
    [
        ("assetName", "../bad.zip"),
        ("sha256", "bad"),
        ("sourcePath", "../plugins/../../private.zip"),
        ("sourcePath", "../plugins/_old/private.zip"),
    ],
)
def test_invalid_manifest_does_not_modify_registry(
    scripts, package_repo, tmp_path, key, value
):
    registry = package_repo / "REGISTRY/plugins.json"
    original = registry.read_bytes()
    with pytest.raises(ValueError):
        _publish(
            scripts,
            package_repo,
            tmp_path,
            lambda package: package.update({key: value}),
        )
    assert registry.read_bytes() == original


def test_workflow_builds_from_tag_before_publishing_zip_registry():
    workflow = (
        Path(__file__).resolve().parents[1] / ".github/workflows/release.yml"
    ).read_text(encoding="utf-8")
    assert (
        workflow.index('git checkout --detach "refs/tags/$TAG"')
        < workflow.index("Build ZIP assets from the tag")
        < workflow.index("gh release upload")
        < workflow.index("Publish ZIP URLs and checksums to main")
    )
    assert "pages/builds" in workflow
    assert "refs/heads/deploy" not in workflow
