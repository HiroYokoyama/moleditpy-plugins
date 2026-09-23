"""
Tests for scripts/update_intra_repo_metadata.py.

Covers the optional_dependencies sync for in-repo plugins: the registry is
generated from the plugin source, so a declared PLUGIN_OPTIONAL_DEPENDENCIES
has to reach the entry (in canonical field order), and a plugin that declares
nothing must leave the entry untouched.
"""

import json
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1] / "scripts"))

import update_intra_repo_metadata as sync


def _write_registry(tmp_path: Path, entry: dict) -> Path:
    registry_dir = tmp_path / "REGISTRY"
    registry_dir.mkdir()
    path = registry_dir / "plugins.json"
    path.write_text(json.dumps([entry], indent=2) + "\n", encoding="utf-8")
    return path


def _write_plugin(tmp_path: Path, body: str) -> None:
    plugin_dir = tmp_path / "plugins" / "Demo"
    plugin_dir.mkdir(parents=True)
    (plugin_dir / "demo.py").write_text(body, encoding="utf-8")


def _entry() -> dict:
    return {
        "id": "demo",
        "visible": True,
        "supported_moleditpy_version": ">=4.0.0, <5.0.0",
        "supported_python_version": ">=3.9, <3.15",
        "name": "Demo",
        "version": "2026.01.01",
        "dependencies": ["numpy"],
        "downloadUrl": "../plugins/Demo/demo.py",
        "sha256": "stale",
    }


def test_read_optional_dependencies_from_source(tmp_path):
    src = tmp_path / "demo.py"
    src.write_text(
        'PLUGIN_OPTIONAL_DEPENDENCIES = [\n    "matplotlib>=3.5",\n    "pillow",\n]\n',
        encoding="utf-8",
    )
    assert sync.infer_optional_dependencies_from_target(src) == [
        "matplotlib>=3.5",
        "pillow",
    ]


def test_undeclared_constant_reads_as_none(tmp_path):
    src = tmp_path / "demo.py"
    src.write_text('PLUGIN_DEPENDENCIES = ["numpy"]\n', encoding="utf-8")
    assert sync.infer_optional_dependencies_from_target(src) is None


def test_sync_writes_optional_dependencies_after_dependencies(tmp_path):
    registry = _write_registry(tmp_path, _entry())
    _write_plugin(
        tmp_path,
        'PLUGIN_VERSION = "2026.01.01"\nPLUGIN_OPTIONAL_DEPENDENCIES = ["matplotlib"]\n',
    )

    sync.update_single_json(registry)

    data = json.loads(registry.read_text(encoding="utf-8"))
    assert data[0]["optional_dependencies"] == ["matplotlib"]
    keys = list(data[0])
    assert keys.index("optional_dependencies") == keys.index("dependencies") + 1


def test_sync_leaves_entry_without_the_constant_alone(tmp_path):
    registry = _write_registry(tmp_path, _entry())
    _write_plugin(tmp_path, 'PLUGIN_VERSION = "2026.01.01"\n')

    sync.update_single_json(registry)

    data = json.loads(registry.read_text(encoding="utf-8"))
    assert "optional_dependencies" not in data[0]


def test_sync_clears_optional_dependencies_when_emptied(tmp_path):
    entry = _entry()
    entry["optional_dependencies"] = ["matplotlib"]
    registry = _write_registry(tmp_path, entry)
    _write_plugin(
        tmp_path,
        'PLUGIN_VERSION = "2026.01.01"\nPLUGIN_OPTIONAL_DEPENDENCIES = []\n',
    )

    sync.update_single_json(registry)

    data = json.loads(registry.read_text(encoding="utf-8"))
    assert data[0]["optional_dependencies"] == []


def test_sync_takes_name_description_and_tags_from_code(tmp_path):
    entry = _entry()
    entry["description"] = "Old wording."
    entry["tags"] = ["Utility"]
    registry = _write_registry(tmp_path, entry)
    _write_plugin(
        tmp_path,
        'PLUGIN_VERSION = "2026.01.02"\n'
        'PLUGIN_NAME = "Demo Plus"\n'
        "PLUGIN_DESCRIPTION = (\n"
        '    "New wording, "\n'
        '    "split over two lines."\n'
        ")\n"
        'PLUGIN_TAGS = ["Analysis", "Visualization"]\n',
    )

    sync.update_single_json(registry)

    data = json.loads(registry.read_text(encoding="utf-8"))
    assert data[0]["name"] == "Demo Plus"
    assert data[0]["description"] == "New wording, split over two lines."
    assert data[0]["tags"] == ["Analysis", "Visualization"]


def test_sync_keeps_registry_values_the_code_does_not_declare(tmp_path):
    entry = _entry()
    entry["description"] = "Curated wording."
    entry["tags"] = ["Utility"]
    registry = _write_registry(tmp_path, entry)
    _write_plugin(tmp_path, 'PLUGIN_VERSION = "2026.01.01"\n')

    sync.update_single_json(registry)

    data = json.loads(registry.read_text(encoding="utf-8"))
    assert data[0]["name"] == "Demo"
    assert data[0]["description"] == "Curated wording."
    assert data[0]["tags"] == ["Utility"]


def test_sync_leaves_retired_old_plugins_frozen(tmp_path):
    entry = _entry()
    entry["downloadUrl"] = "../plugins/_old/Demo/demo.py"
    entry["tags"] = ["Utility"]
    registry = _write_registry(tmp_path, entry)
    plugin_dir = tmp_path / "plugins" / "_old" / "Demo"
    plugin_dir.mkdir(parents=True)
    (plugin_dir / "demo.py").write_text(
        'PLUGIN_VERSION = "2026.01.01"\nPLUGIN_TAGS = ["Analysis"]\n',
        encoding="utf-8",
    )

    sync.update_single_json(registry)

    data = json.loads(registry.read_text(encoding="utf-8"))
    assert data[0]["tags"] == ["Utility"]


def test_sync_takes_dependencies_and_supported_os_from_code(tmp_path):
    entry = _entry()
    entry["supported_os"] = ["Windows", "macOS", "Linux", "WSL"]
    registry = _write_registry(tmp_path, entry)
    _write_plugin(
        tmp_path,
        'PLUGIN_VERSION = "2026.01.02"\n'
        'PLUGIN_DEPENDENCIES = ["numpy", "scipy"]\n'
        'PLUGIN_SUPPORTED_OS = ["macOS", "Linux", "WSL"]\n',
    )

    sync.update_single_json(registry)

    data = json.loads(registry.read_text(encoding="utf-8"))
    assert data[0]["dependencies"] == ["numpy", "scipy"]
    assert data[0]["supported_os"] == ["macOS", "Linux", "WSL"]


def test_sync_writes_declared_empty_dependencies(tmp_path):
    registry = _write_registry(tmp_path, _entry())
    _write_plugin(tmp_path, 'PLUGIN_VERSION = "2026.01.02"\nPLUGIN_DEPENDENCIES = []\n')

    sync.update_single_json(registry)

    data = json.loads(registry.read_text(encoding="utf-8"))
    assert data[0]["dependencies"] == []


def test_sync_keeps_dependencies_the_code_does_not_declare(tmp_path):
    registry = _write_registry(tmp_path, _entry())
    _write_plugin(tmp_path, 'PLUGIN_VERSION = "2026.01.02"\n')

    sync.update_single_json(registry)

    data = json.loads(registry.read_text(encoding="utf-8"))
    assert data[0]["dependencies"] == ["numpy"]
