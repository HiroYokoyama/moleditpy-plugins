"""Curated migration targets are validated atomically and survive regeneration."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/update_registry_migrations.py"
SPEC = importlib.util.spec_from_file_location("registry_migrations", SCRIPT)
migrations = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(migrations)


@pytest.fixture
def registry(tmp_path):
    path = tmp_path / "plugins.json"
    path.write_text(
        json.dumps(
            [
                {
                    "id": "old",
                    "name": "Old",
                    "visible": False,
                    "version": "1",
                    "sha256": "keep",
                },
                {"id": "old2", "name": "Old 2", "visible": False},
                {"id": "new", "name": "New", "visible": True},
                {"id": "new2", "name": "New 2", "visible": True},
            ]
        ),
        encoding="utf-8",
    )
    return path


def test_sets_multiple_targets_without_touching_existing_metadata(registry):
    before = json.loads(registry.read_text())
    assert migrations.update_migrations(registry, ["old=new", "old2=new2"]) == 2
    after = json.loads(registry.read_text())
    assert after[0].pop("migration_to") == "new"
    assert after[1].pop("migration_to") == "new2"
    assert after == before


@pytest.mark.parametrize(
    "assignment",
    [
        "missing=new",
        "old=missing",
        "old",
        "old=",
        "old=old",
        "old=old2",
        "new=new2",
    ],
)
def test_invalid_assignment_does_not_partially_write(registry, assignment):
    before = registry.read_bytes()
    with pytest.raises(ValueError):
        migrations.update_migrations(registry, ["old=new", assignment])
    assert registry.read_bytes() == before


def test_can_update_existing_target(registry):
    migrations.update_migrations(registry, ["old=new"])
    migrations.update_migrations(registry, ["old=new2"])
    assert json.loads(registry.read_text())[0]["migration_to"] == "new2"


def test_cli_documents_curated_id_assignments():
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--set" in result.stdout
    assert "old_plugin_id=new_plugin_id" in result.stdout


def test_real_registry_migrations_point_to_visible_plugins():
    path = SCRIPT.parents[1] / "REGISTRY/plugins.json"
    entries = json.loads(path.read_text(encoding="utf-8"))
    by_id = {entry["id"]: entry for entry in entries}
    for entry in entries:
        if "migration_to" not in entry:
            continue
        assert entry["visible"] is False
        assert entry["migration_to"] != entry["id"]
        assert by_id[entry["migration_to"]]["visible"] is True


def test_metadata_regeneration_preserves_migration(registry):
    sys.path.insert(0, str(SCRIPT.parent))
    import update_intra_repo_metadata as metadata

    migrations.update_migrations(registry, ["old=new"])
    metadata.update_single_json(registry)
    assert json.loads(registry.read_text())[0]["migration_to"] == "new"
