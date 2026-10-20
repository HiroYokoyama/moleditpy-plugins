"""Set reviewed registry-only license values without modifying plugin sources.

Usage: python scripts/update_registry_licenses.py --set plugin_id=LGPL-3.0-or-later
Multiple --set arguments are allowed. Review the actual plugin license first.
"""

import argparse
import json
from pathlib import Path

from license_metadata import normalize_license


def update_licenses(registry_path, assignments):
    plugins = json.loads(registry_path.read_text(encoding="utf-8"))
    by_id = {p["id"]: p for p in plugins}
    reviewed = {}
    for assignment in assignments:
        plugin_id, separator, value = assignment.partition("=")
        if not separator or plugin_id not in by_id:
            raise ValueError("Expected an existing plugin_id=license assignment")
        if not by_id[plugin_id].get("visible", False):
            raise ValueError("Retired/hidden entries are left unchanged")
        value = normalize_license(value)
        if value is None:
            raise ValueError("Use NOASSERTION when a license is unconfirmed")
        reviewed[plugin_id] = value
    for plugin_id, value in reviewed.items():
        by_id[plugin_id]["license"] = value
    with registry_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(plugins, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return len(reviewed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--set", dest="assignments", action="append", required=True)
    args = parser.parse_args()
    registry = Path(__file__).resolve().parents[1] / "REGISTRY/plugins.json"
    print("Updated licenses:", update_licenses(registry, args.assignments))


if __name__ == "__main__":
    main()
