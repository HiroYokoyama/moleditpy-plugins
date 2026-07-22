"""
Tests for the VDW Radii Overlay plugin: load_settings / save_settings.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as real_numpy  # noqa: F401  (imported so mocks_with_real_numpy can snapshot it)
import pytest

from conftest import load_plugin, mock_optional_imports, mocks_with_real_numpy

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
        assert on_disk == {
            "occupancy": 0.9,
            "resolution": 0.2,
            "base_style": "default",
        }

    def test_defaults(self, tmp_path):
        mod = self._mod(tmp_path)
        assert mod._vdw_settings["occupancy"] == pytest.approx(0.3)
        assert mod._vdw_settings["resolution"] == pytest.approx(0.125)
        assert mod._vdw_settings["base_style"] == "default"


class TestVDWBaseStyleSetting:
    def _mod(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
        return mod

    def test_round_trip_base_style_stick(self, tmp_path):
        mod = self._mod(tmp_path)
        mod._vdw_settings["base_style"] = "stick"
        mod.save_settings()
        mod._vdw_settings["base_style"] = "default"
        mod.load_settings()
        assert mod._vdw_settings["base_style"] == "stick"

    def test_invalid_base_style_in_file_ignored(self, tmp_path):
        mod = self._mod(tmp_path)
        Path(mod.SETTINGS_FILE).write_text(json.dumps({"base_style": "bogus"}))
        mod.load_settings()
        assert mod._vdw_settings["base_style"] == "default"

    def test_missing_base_style_key_keeps_default(self, tmp_path):
        mod = self._mod(tmp_path)
        Path(mod.SETTINGS_FILE).write_text(json.dumps({"occupancy": 0.6}))
        mod.load_settings()
        assert mod._vdw_settings["base_style"] == "default"

    def test_base_style_to_override_map(self, tmp_path):
        mod = self._mod(tmp_path)
        assert mod._BASE_STYLE_TO_OVERRIDE["default"] == "ball_and_stick"
        assert mod._BASE_STYLE_TO_OVERRIDE["stick"] == "stick"


# ---------------------------------------------------------------------------
# draw_vdw_overlay: base model style_override selection
# ---------------------------------------------------------------------------


class TestDrawVdwOverlayBaseStyle:
    def _make_mw_and_mol(self):
        mw = MagicMock()
        mw.view_3d_manager._plugin_color_overrides = {}
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 0  # skip the heavy overlay path
        return mw, mol

    def test_default_base_style_uses_ball_and_stick(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
        mod._vdw_settings["base_style"] = "default"
        mw, mol = self._make_mw_and_mol()
        mod.draw_vdw_overlay(mw, mol)
        mw.view_3d_manager.draw_standard_3d_style.assert_called_once_with(
            mol, style_override="ball_and_stick"
        )

    def test_stick_base_style_uses_stick_override(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
        mod._vdw_settings["base_style"] = "stick"
        mw, mol = self._make_mw_and_mol()
        mod.draw_vdw_overlay(mw, mol)
        mw.view_3d_manager.draw_standard_3d_style.assert_called_once_with(
            mol, style_override="stick"
        )


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


# ---------------------------------------------------------------------------
# draw_vdw_overlay: full render path with real numpy math
# ---------------------------------------------------------------------------


class _FakePeriodicTable:
    """Fixed VDW radii so the SDF math below is deterministic."""

    def __init__(self, radii):
        self._radii = radii

    def GetRvdw(self, atomic_num):
        return self._radii.get(atomic_num, 1.5)


def _make_atom(symbol, atomic_num):
    a = MagicMock()
    a.GetSymbol.return_value = symbol
    a.GetAtomicNum.return_value = atomic_num
    return a


def _make_mol_with_conformer(atoms_spec, positions):
    """atoms_spec: list of (symbol, atomic_num); positions: list of (x, y, z)."""
    mol = MagicMock()
    mol.GetNumAtoms.return_value = len(atoms_spec)
    mol.GetNumConformers.return_value = 1
    atoms = [_make_atom(sym, num) for sym, num in atoms_spec]
    mol.GetAtomWithIdx.side_effect = lambda i: atoms[i]

    conf = MagicMock()

    def _get_pos(i):
        x, y, z = positions[i]
        p = MagicMock()
        p.x, p.y, p.z = x, y, z
        return p

    conf.GetAtomPosition.side_effect = _get_pos
    mol.GetConformer.return_value = conf
    return mol


def _wire_fake_grid(mod, mesh_points):
    """Replace mod.pv.ImageData() with a fake grid whose .contour() returns a
    fake mesh exposing real numpy point coordinates, so the vertex-coloring
    loop runs with real math."""
    mesh = MagicMock()
    mesh.points = mesh_points
    mesh.point_data = {}
    grid = MagicMock()
    grid.point_data = {}
    grid.contour.return_value = mesh
    mod.pv.ImageData.return_value = grid
    return grid, mesh


class TestDrawVdwOverlayFullRenderPath:
    def _base_mw(self):
        mw = MagicMock()
        mw.view_3d_manager._plugin_color_overrides = {}
        mw.view_3d_manager.plotter = MagicMock()
        return mw

    def test_full_path_adds_mesh_with_default_cpk_colors(self, tmp_path):
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
            mod._vdw_settings["resolution"] = 0.5
            mod._vdw_settings["occupancy"] = 0.4
            mod.pt = _FakePeriodicTable({6: 1.7, 8: 1.52})
            mod.CPK_COLORS_PV = {"C": [0.5, 0.5, 0.5], "O": [1.0, 0.0, 0.0]}

            mol = _make_mol_with_conformer(
                [("C", 6), ("O", 8)], [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)]
            )
            mesh_points = mod.np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
            grid, mesh = _wire_fake_grid(mod, mesh_points)

            mw = self._base_mw()
            mod.draw_vdw_overlay(mw, mol)

            mw.view_3d_manager.draw_standard_3d_style.assert_called_once_with(
                mol, style_override="ball_and_stick"
            )
            grid.contour.assert_called_once()
            mw.view_3d_manager.plotter.add_mesh.assert_called_once()
            call = mw.view_3d_manager.plotter.add_mesh.call_args
            assert call.kwargs["opacity"] == pytest.approx(0.4)
            assert call.kwargs["name"] == "vdw_overlay_mesh"
            # Nearest-atom coloring assigned to the fake mesh's point_data.
            colors = mesh.point_data["AtomColors"]
            assert colors.shape == (2, 3)

    def test_full_path_hex_color_override(self, tmp_path):
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
            mod._vdw_settings["resolution"] = 0.5
            mod.pt = _FakePeriodicTable({6: 1.7})

            mol = _make_mol_with_conformer([("C", 6)], [(0.0, 0.0, 0.0)])
            mesh_points = mod.np.array([[0.0, 0.0, 0.0]])
            grid, mesh = _wire_fake_grid(mod, mesh_points)

            mw = self._base_mw()
            mw.view_3d_manager._plugin_color_overrides = {0: "#00ff00"}
            mod.draw_vdw_overlay(mw, mol)

            colors = mesh.point_data["AtomColors"]
            # Green: R=0, G=1, B=0
            assert colors[0][0] == pytest.approx(0.0)
            assert colors[0][1] == pytest.approx(1.0)
            assert colors[0][2] == pytest.approx(0.0)

    def test_full_path_legacy_255_scale_color_override(self, tmp_path):
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
            mod._vdw_settings["resolution"] = 0.5
            mod.pt = _FakePeriodicTable({6: 1.7})

            mol = _make_mol_with_conformer([("C", 6)], [(0.0, 0.0, 0.0)])
            mesh_points = mod.np.array([[0.0, 0.0, 0.0]])
            grid, mesh = _wire_fake_grid(mod, mesh_points)

            mw = self._base_mw()
            mw.view_3d_manager._plugin_color_overrides = {0: [255, 0, 0]}
            mod.draw_vdw_overlay(mw, mol)

            colors = mesh.point_data["AtomColors"]
            assert colors[0][0] == pytest.approx(1.0)
            assert colors[0][1] == pytest.approx(0.0)

    def test_full_path_legacy_custom_atom_colors_fallback(self, tmp_path):
        """When view_3d_manager has no overrides, fall back to mw.custom_atom_colors."""
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
            mod._vdw_settings["resolution"] = 0.5
            mod.pt = _FakePeriodicTable({6: 1.7})

            mol = _make_mol_with_conformer([("C", 6)], [(0.0, 0.0, 0.0)])
            mesh_points = mod.np.array([[0.0, 0.0, 0.0]])
            grid, mesh = _wire_fake_grid(mod, mesh_points)

            mw = self._base_mw()
            mw.view_3d_manager._plugin_color_overrides = {}
            mw.custom_atom_colors = {0: [0.25, 0.5, 0.75]}
            mod.draw_vdw_overlay(mw, mol)

            colors = mesh.point_data["AtomColors"]
            assert colors[0][0] == pytest.approx(0.25)

    def test_no_view_3d_manager_uses_legacy_view3d_fallback(self, tmp_path):
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")

            mol = MagicMock()
            mol.GetNumAtoms.return_value = 0

            mw = MagicMock(spec=["view3d", "custom_atom_colors"])
            mod.draw_vdw_overlay(mw, mol)

            mw.view3d.draw_standard_3d_style.assert_called_once_with(
                mol, style_override="ball_and_stick"
            )

    def test_no_conformer_skips_overlay_mesh(self, tmp_path):
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")

            mol = MagicMock()
            mol.GetNumAtoms.return_value = 2
            mol.GetNumConformers.return_value = 0

            mw = self._base_mw()
            mod.draw_vdw_overlay(mw, mol)

            mw.view_3d_manager.plotter.add_mesh.assert_not_called()

    def test_resolution_below_minimum_is_clamped(self, tmp_path):
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
            mod._vdw_settings["resolution"] = 0.0001  # below the 0.01 floor
            mod.pt = _FakePeriodicTable({6: 1.7})

            mol = _make_mol_with_conformer([("C", 6)], [(0.0, 0.0, 0.0)])
            mesh_points = mod.np.array([[0.0, 0.0, 0.0]])
            grid, mesh = _wire_fake_grid(mod, mesh_points)

            mw = self._base_mw()
            mod.draw_vdw_overlay(mw, mol)

            # spacing kwarg is set on the fake grid from the clamped value.
            assert grid.spacing == (0.01, 0.01, 0.01)

    def test_no_plotter_skips_add_mesh(self, tmp_path):
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
            mod._vdw_settings["resolution"] = 0.5
            mod.pt = _FakePeriodicTable({6: 1.7})

            mol = _make_mol_with_conformer([("C", 6)], [(0.0, 0.0, 0.0)])
            mesh_points = mod.np.array([[0.0, 0.0, 0.0]])
            _wire_fake_grid(mod, mesh_points)

            mw = self._base_mw()
            mw.view_3d_manager.plotter = None
            mod.draw_vdw_overlay(mw, mol)  # must not raise

    def test_empty_contour_mesh_skips_coloring(self, tmp_path):
        """An empty iso-surface (n_mesh_pts == 0) must not crash."""
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
            mod._vdw_settings["resolution"] = 0.5
            mod.pt = _FakePeriodicTable({6: 1.7})

            mol = _make_mol_with_conformer([("C", 6)], [(0.0, 0.0, 0.0)])
            mesh_points = mod.np.zeros((0, 3))
            grid, mesh = _wire_fake_grid(mod, mesh_points)

            mw = self._base_mw()
            mod.draw_vdw_overlay(mw, mol)

            mw.view_3d_manager.plotter.add_mesh.assert_not_called()

    def test_exception_during_render_is_caught_and_logged(self, tmp_path, caplog):
        with mocks_with_real_numpy():
            mod = load_plugin(VDW_PATH)
            mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
            mod.pt = _FakePeriodicTable({6: 1.7})

            mol = _make_mol_with_conformer([("C", 6)], [(0.0, 0.0, 0.0)])
            # pv.ImageData raises -> exercises the except branch.
            mod.pv.ImageData.side_effect = RuntimeError("boom")

            mw = self._base_mw()
            mod.draw_vdw_overlay(mw, mol)  # must not raise
            assert any("VDW Overlay Error" in r.message for r in caplog.records)


class TestOpenSettings:
    def test_creates_new_window_when_none_registered(self):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        ctx = MagicMock()
        ctx.get_window.return_value = None
        created = MagicMock()
        mod.VDWConfigWindow = MagicMock(return_value=created)
        mod.open_settings(ctx)
        created.show.assert_called_once()

    def test_reuses_existing_window(self):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        ctx = MagicMock()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        mod.open_settings(ctx)
        existing.show.assert_called_once()
        existing.raise_.assert_called_once()
        existing.activateWindow.assert_called_once()


class TestRunEntryPoint:
    def test_run_with_context_opens_settings(self):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        ctx = MagicMock()
        mod.initialize(ctx)
        ctx.get_window.return_value = None
        mod.VDWConfigWindow = MagicMock()
        mw = MagicMock(spec=["host"])
        mw.host = MagicMock(spec=[])
        mod.run(mw)
        ctx.get_window.assert_called_once()

    def test_run_unwraps_host_attribute(self):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        ctx = MagicMock()
        mod.initialize(ctx)
        ctx.get_window.return_value = MagicMock()
        mw = MagicMock()
        mw.host = MagicMock(spec=[])
        mod.run(mw)
        ctx.get_window.return_value.show.assert_called_once()
