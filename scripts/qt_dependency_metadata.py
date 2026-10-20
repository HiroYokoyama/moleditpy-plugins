"""Repair inherited Qt metadata for the migrated MoleditPy 5 registry."""


def normalize_inherited_qt_dependencies(entry: dict) -> bool:
    """Replace the obsolete binding only for visible, MoleditPy 5 entries."""
    spec = str(entry.get("supported_moleditpy_version", "")).replace(" ", "")
    if (
        not entry.get("visible", True)
        or spec != ">=5.0.0,<6.0.0"
        or "/_old/" in str(entry.get("downloadUrl", ""))
    ):
        return False
    dependencies = entry.get("dependencies", [])
    if "PyQt6" not in dependencies:
        return False
    entry["dependencies"] = list(
        dict.fromkeys(
            "PySide6-Essentials" if dependency == "PyQt6" else dependency
            for dependency in dependencies
        )
    )
    return True
