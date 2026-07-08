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




# ---------------------------------------------------------------------------
# settings edge cases + initialize
# ---------------------------------------------------------------------------

from unittest.mock import MagicMock

from conftest import make_context


class TestVDWSettingsEdgeCases:
    def _mod(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
        return mod

    def test_partial_file_keeps_other_default(self, tmp_path):
        mod = self._mod(tmp_path)
        Path(mod.SETTINGS_FILE).write_text(json.dumps({"occupancy": 0.7}))
        mod.load_settings()
        assert mod._vdw_settings["occupancy"] == pytest.approx(0.7)
        assert mod._vdw_settings["resolution"] == pytest.approx(0.125)

    def test_string_values_coerced_to_float(self, tmp_path):
        mod = self._mod(tmp_path)
        Path(mod.SETTINGS_FILE).write_text(
            json.dumps({"occupancy": "0.55", "resolution": "0.25"})
        )
        mod.load_settings()
        assert mod._vdw_settings["occupancy"] == pytest.approx(0.55)
        assert mod._vdw_settings["resolution"] == pytest.approx(0.25)

    def test_unknown_keys_ignored(self, tmp_path):
        mod = self._mod(tmp_path)
        Path(mod.SETTINGS_FILE).write_text(
            json.dumps({"bogus": 1, "occupancy": 0.4})
        )
        mod.load_settings()
        assert "bogus" not in mod._vdw_settings
        assert mod._vdw_settings["occupancy"] == pytest.approx(0.4)

    def test_save_writes_valid_json_with_both_keys(self, tmp_path):
        mod = self._mod(tmp_path)
        mod._vdw_settings["occupancy"] = 0.9
        mod._vdw_settings["resolution"] = 0.2
        mod.save_settings()
        on_disk = json.loads(Path(mod.SETTINGS_FILE).read_text())
        assert on_disk == {"occupancy": 0.9, "resolution": 0.2}

    def test_defaults(self, tmp_path):
        mod = self._mod(tmp_path)
        assert mod._vdw_settings["occupancy"] == pytest.approx(0.3)
        assert mod._vdw_settings["resolution"] == pytest.approx(0.125)


class TestVDWInitialize:
    def test_registers_vdw_overlay_style(self):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        ctx.register_3d_style.assert_called_once()
        name, drawer = ctx.register_3d_style.call_args[0]
        assert name == "vdw_overlay"
        assert drawer is mod.draw_vdw_overlay

    def test_run_without_context_is_noop(self):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        assert mod.PLUGIN_CONTEXT is None
        mod.run(MagicMock())
