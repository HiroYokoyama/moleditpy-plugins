"""Legacy dependency values must not reinstall PyQt6 into a PySide6 host."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from qt_dependency_metadata import normalize_inherited_qt_dependencies


def test_migrated_plugin_inherited_dependencies_are_repaired():
    entry = {
        "visible": True,
        "supported_moleditpy_version": ">=5.0.0, <6.0.0",
        "dependencies": ["numpy", "PyQt6", "PySide6-Essentials"],
    }
    assert normalize_inherited_qt_dependencies(entry)
    assert entry["dependencies"] == ["numpy", "PySide6-Essentials"]
    assert not normalize_inherited_qt_dependencies(entry)


@pytest.mark.parametrize(
    "overrides",
    [
        {"visible": False},
        {"supported_moleditpy_version": ">=4.0.0, <5.0.0"},
        {"downloadUrl": "../plugins/_old/example.py"},
    ],
)
def test_historical_or_hidden_entries_stay_unchanged(overrides):
    entry = {
        "visible": True,
        "supported_moleditpy_version": ">=5.0.0, <6.0.0",
        "dependencies": ["PyQt6"],
        **overrides,
    }
    assert not normalize_inherited_qt_dependencies(entry)
    assert entry["dependencies"] == ["PyQt6"]
