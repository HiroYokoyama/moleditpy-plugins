"""
Tests for the VDW Radii Overlay plugin: load_settings / save_settings.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
VDW_PATH = PLUGINS_DIR / "VDW_Radii_Overlay" / "vdw_radii_overlay.py"

with mock_optional_imports():
    _vdw = load_plugin(VDW_PATH)


class TestVDWSettings:
    def test_round_trip_occupancy(self, tmp_path, monkeypatch):
        fresh = {"occupancy": 0.3, "resolution": 0.125}
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "vdw.json"))
        monkeypatch.setattr(_vdw, "_vdw_settings", fresh)
        fresh["occupancy"] = 0.75
        _vdw.save_settings()
        fresh["occupancy"] = 0.3  # reset
        _vdw.load_settings()
        assert fresh["occupancy"] == pytest.approx(0.75)

    def test_round_trip_resolution(self, tmp_path, monkeypatch):
        fresh = {"occupancy": 0.3, "resolution": 0.125}
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "vdw.json"))
        monkeypatch.setattr(_vdw, "_vdw_settings", fresh)
        fresh["resolution"] = 0.25
        _vdw.save_settings()
        fresh["resolution"] = 0.125
        _vdw.load_settings()
        assert fresh["resolution"] == pytest.approx(0.25)

    def test_load_missing_file_is_noop(self, tmp_path, monkeypatch):
        fresh = {"occupancy": 0.3, "resolution": 0.125}
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "nonexistent.json"))
        monkeypatch.setattr(_vdw, "_vdw_settings", fresh)
        _vdw.load_settings()  # must not raise
        assert fresh["occupancy"] == pytest.approx(0.3)

    def test_saved_file_is_valid_json(self, tmp_path, monkeypatch):
        fresh = {"occupancy": 0.5, "resolution": 0.2}
        path = tmp_path / "vdw.json"
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(path))
        monkeypatch.setattr(_vdw, "_vdw_settings", fresh)
        _vdw.save_settings()
        data = json.loads(path.read_text())
        assert "occupancy" in data
        assert "resolution" in data
