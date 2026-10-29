"""Set curated migration targets for hidden plugins without editing their sources.

Usage: python scripts/update_registry_migrations.py --set old_plugin_id=new_plugin_id
Multiple --set arguments are allowed. Targets must be visible registry entries.
"""

import argparse
import json
from pathlib import Path


def update_migrations(registry_path: Path, assignments: list[str]) -> int:
    plugins = json.loads(registry_path.read_text(encoding="utf-8"))
    by_id = {entry["id"]: entry for entry in plugins}
    reviewed = {}
    for assignment in assignments:
        plugin_id, separator, target_id = assignment.partition("=")
        if not separator or plugin_id not in by_id or target_id not in by_id:
            raise ValueError("Expected existing old_plugin_id=new_plugin_id values")
        if by_id[plugin_id].get("visible", True) is not False:
            raise ValueError("Migration sources must be hidden registry entries")
        if plugin_id == target_id or not by_id[target_id].get("visible", True):
            raise ValueError("Migration targets must be different, visible plugins")
        reviewed[plugin_id] = target_id
    for plugin_id, target_id in reviewed.items():
        by_id[plugin_id]["migration_to"] = target_id
    with registry_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(plugins, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return len(reviewed)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--set", dest="assignments", action="append", required=True)
    args = parser.parse_args()
    registry = Path(__file__).resolve().parents[1] / "REGISTRY/plugins.json"
    print("Updated migration targets:", update_migrations(registry, args.assignments))


if __name__ == "__main__":
    main()
