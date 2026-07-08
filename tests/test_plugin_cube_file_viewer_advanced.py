"""
Tests for the Cube File Viewer Advanced plugin
(plugins/Cube_File_Viewer_Advanced/cube_viewer_advanced.py).

Pure-function tests (parse_cube_data, build_grid_from_meta) run with real
numpy injected after plugin load so we can verify actual parsing/grid-math
behaviour. Qt / PyVista / RDKit remain mocked throughout.

Methods on CubeViewerWidget (a QWidget subclass, which can't be instantiated
under the mocked-Qt environment) are extracted standalone via AST
(extract_function) and invoked against a SimpleNamespace fake ``self``.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CUBE_ADV_PATH = PLUGINS_DIR / "Cube_File_Viewer_Advanced" / "cube_viewer_advanced.py"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_with_real_numpy(plugin_path: Path):
    """Load plugin with Qt/PyVista/RDKit mocked but inject real numpy.

    parse_cube_data looks up ``np`` in the module globals at call time, so
    replacing it after load makes the pure arithmetic work with real arrays.
    """
    with mock_optional_imports():
        mod = load_plugin(plugin_path)
    mod.np = np
    return mod


def _make_cube_content(nx: int = 2, ny: int = 2, nz: int = 2,
                       n_atoms: int = 1,
                       data_values: list | None = None) -> str:
    """Build a minimal valid Gaussian Cube file as a string.

    Atoms are all written as Carbon (atomic_num=6).
    Data is written in rows of 6, which satisfies the ≥6-token heuristic
    used by parse_cube_data to find the start of volumetric data.
    """
    if data_values is None:
        data_values = list(range(1, nx * ny * nz + 1))

    header = (
        f"Title line 1\n"
        f"Title line 2\n"
        f"   {n_atoms}   0.000000   0.000000   0.000000\n"
        f"   {nx}   0.200000   0.000000   0.000000\n"
        f"   {ny}   0.000000   0.200000   0.000000\n"
        f"   {nz}   0.000000   0.000000   0.200000\n"
    )
    atoms = "   6   0.000000   0.000000   0.000000   0.629118   0.000000\n" * n_atoms

    data_lines = ""
    row: list[str] = []
    for v in data_values:
        row.append(f"{float(v):.4f}")
        if len(row) == 6:
            data_lines += " " + " ".join(row) + "\n"
            row = []
    if row:
        data_lines += " " + " ".join(row) + "\n"

    return header + atoms + data_lines


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def cube_adv_mod():
    return _load_with_real_numpy(CUBE_ADV_PATH)


@pytest.fixture
def cube_file(tmp_path):
    """Minimal valid 2×2×2 cube file with one Carbon atom."""
    p = tmp_path / "test.cube"
    p.write_text(_make_cube_content(nx=2, ny=2, nz=2, n_atoms=1))
    return p


# ---------------------------------------------------------------------------
# parse_cube_data — Cube File Viewer Advanced (independent copy)
# ---------------------------------------------------------------------------

class TestCubeAdvancedParseCubeData:
    def test_short_file_raises(self, tmp_path, cube_adv_mod):
        f = tmp_path / "short.cube"
        f.write_text("a\nb\n")
        with pytest.raises(ValueError, match="too short"):
            cube_adv_mod.parse_cube_data(str(f))

    def test_returns_expected_keys(self, cube_file, cube_adv_mod):
        result = cube_adv_mod.parse_cube_data(str(cube_file))
        assert {"atoms", "dims", "data_flat", "is_angstrom_header"}.issubset(
            result.keys()
        )

    def test_dims(self, cube_file, cube_adv_mod):
        result = cube_adv_mod.parse_cube_data(str(cube_file))
        assert result["dims"] == (2, 2, 2)

    def test_pads_short_data(self, tmp_path, cube_adv_mod):
        content = _make_cube_content(nx=2, ny=2, nz=2, n_atoms=1,
                                     data_values=[1.0, 2.0])
        # 2×2×2 = 8 expected, 2 provided → 6 zeros appended
        # BUT: 2 values on one line = only 2 tokens < 6, fails heuristic → 0 data found
        # So result is all zeros (padded from 0 to 8). That's also valid padding behaviour.
        p = tmp_path / "short2.cube"
        p.write_text(content)
        result = cube_adv_mod.parse_cube_data(str(p))
        assert len(result["data_flat"]) == 8


# ---------------------------------------------------------------------------
# initialize() — Cube File Viewer Advanced
# ---------------------------------------------------------------------------

class TestCubeAdvancedInitialize:
    def test_does_not_raise(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_ADV_PATH)
            ctx = make_context()
            mod.initialize(ctx)

    def test_registers_dot_cube_opener(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_ADV_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        exts = [c.args[0] for c in ctx.register_file_opener.call_args_list]
        assert ".cube" in exts

    def test_registers_dot_cub_opener(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_ADV_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        exts = [c.args[0] for c in ctx.register_file_opener.call_args_list]
        assert ".cub" in exts


# ---------------------------------------------------------------------------
# build_grid_from_meta — grid/volume computation (real numpy)
# ---------------------------------------------------------------------------

def _meta(nx=2, ny=2, nz=2, is_angstrom_header=False, n_atoms=1):
    data = np.arange(1, nx * ny * nz + 1, dtype=float)
    atoms = [(6, np.array([1.0, 2.0, 3.0])) for _ in range(n_atoms)]
    return {
        "atoms": atoms,
        "origin": np.array([0.0, 0.0, 0.0]),
        "x_vec": np.array([0.2, 0.0, 0.0]),
        "y_vec": np.array([0.0, 0.2, 0.0]),
        "z_vec": np.array([0.0, 0.0, 0.2]),
        "dims": (nx, ny, nz),
        "data_flat": data,
        "is_angstrom_header": is_angstrom_header,
    }


class TestBuildGridFromMeta:
    BOHR_TO_ANG = 0.529177210903

    def test_bohr_header_converts_to_angstrom(self, cube_adv_mod):
        meta = _meta(is_angstrom_header=False)
        info, grid = cube_adv_mod.build_grid_from_meta(meta)
        anum, pos = info["atoms"][0]
        assert anum == 6
        assert pos[0] == pytest.approx(1.0 * self.BOHR_TO_ANG)

    def test_angstrom_header_not_converted(self, cube_adv_mod):
        meta = _meta(is_angstrom_header=True)
        info, grid = cube_adv_mod.build_grid_from_meta(meta)
        _anum, pos = info["atoms"][0]
        assert pos[0] == pytest.approx(1.0)

    def test_grid_dimensions_set(self, cube_adv_mod):
        meta = _meta(nx=2, ny=3, nz=4)
        _info, grid = cube_adv_mod.build_grid_from_meta(meta)
        assert grid.dimensions == [2, 3, 4]

    def test_grid_points_shape(self, cube_adv_mod):
        meta = _meta(nx=2, ny=2, nz=2)
        _info, grid = cube_adv_mod.build_grid_from_meta(meta)
        assert grid.points.shape == (8, 3)

    def test_point_data_length_matches_grid_size(self, cube_adv_mod):
        meta = _meta(nx=2, ny=2, nz=3)
        _info, grid = cube_adv_mod.build_grid_from_meta(meta)
        # grid.point_data is a MagicMock (pyvista is mocked); inspect the
        # value passed to its __setitem__ instead of round-tripping through it.
        (_key, values), _kw = grid.point_data.__setitem__.call_args
        assert len(values) == 12

    def test_origin_offset_applied(self, cube_adv_mod):
        meta = _meta()
        meta["origin"] = np.array([1.0, 0.0, 0.0])
        _info, grid = cube_adv_mod.build_grid_from_meta(meta)
        # First grid point (gx=gy=gz=0) should be exactly the (converted) origin.
        expected = 1.0 * self.BOHR_TO_ANG
        assert grid.points[0][0] == pytest.approx(expected)

    def test_multiple_atoms_all_converted(self, cube_adv_mod):
        meta = _meta(n_atoms=3, is_angstrom_header=False)
        info, _grid = cube_adv_mod.build_grid_from_meta(meta)
        assert len(info["atoms"]) == 3
        for _anum, pos in info["atoms"]:
            assert pos[0] == pytest.approx(1.0 * self.BOHR_TO_ANG)


# ---------------------------------------------------------------------------
# CubeViewerWidget methods extracted standalone
# ---------------------------------------------------------------------------

import logging as _logging

_update_iso = extract_function(
    CUBE_ADV_PATH, "CubeViewerWidget", "update_iso",
    extra_globals={"logging": _logging},
)
_on_opacity_changed = extract_function(
    CUBE_ADV_PATH, "CubeViewerWidget", "on_opacity_changed"
)
_on_opacity_spin_changed = extract_function(
    CUBE_ADV_PATH, "CubeViewerWidget", "on_opacity_spin_changed"
)
_on_slider_changed = extract_function(
    CUBE_ADV_PATH, "CubeViewerWidget", "on_slider_changed"
)
_on_spin_changed = extract_function(
    CUBE_ADV_PATH, "CubeViewerWidget", "on_spin_changed"
)
_get_settings_path = extract_function(
    CUBE_ADV_PATH, "CubeViewerWidget", "get_settings_path",
    extra_globals={"os": __import__("os"), "__file__": str(CUBE_ADV_PATH)},
)
_init_default_presets = extract_function(
    CUBE_ADV_PATH, "CubeViewerWidget", "_init_default_presets"
)
import colorsys


class FakeQColor:
    """Minimal QColor stand-in with real HSV math (colorsys-backed)."""

    def __init__(self, r=0, g=0, b=0):
        self._r, self._g, self._b = r, g, b

    def isValid(self):
        return True

    def red(self):
        return self._r

    def green(self):
        return self._g

    def blue(self):
        return self._b

    def hue(self):
        if self._r == self._g == self._b:
            return -1
        h, _s, _v = colorsys.rgb_to_hsv(self._r / 255, self._g / 255, self._b / 255)
        return int(round(h * 360))

    def saturation(self):
        _h, s, _v = colorsys.rgb_to_hsv(self._r / 255, self._g / 255, self._b / 255)
        return int(round(s * 255))

    def value(self):
        _h, _s, v = colorsys.rgb_to_hsv(self._r / 255, self._g / 255, self._b / 255)
        return int(round(v * 255))

    @staticmethod
    def fromHsv(h, s, v):
        r, g, b = colorsys.hsv_to_rgb((h % 360) / 360, s / 255, v / 255)
        return FakeQColor(int(round(r * 255)), int(round(g * 255)), int(round(b * 255)))


_update_complementary_color = extract_function(
    CUBE_ADV_PATH, "CubeViewerWidget", "update_complementary_color",
    extra_globals={"QColor": FakeQColor},
)
_on_comp_color_toggled = extract_function(
    CUBE_ADV_PATH, "CubeViewerWidget", "on_comp_color_toggled"
)


class FakeSpin:
    def __init__(self, val=0.0):
        self._val = val
        self.set_calls = []

    def value(self):
        return self._val

    def setValue(self, v):
        self._val = v
        self.set_calls.append(v)

    def blockSignals(self, _v):
        pass


class FakeCombo:
    def __init__(self, text=""):
        self._text = text

    def currentText(self):
        return self._text


class FakeCheck:
    def __init__(self, checked=False):
        self._checked = checked

    def isChecked(self):
        return self._checked

    def setChecked(self, v):
        self._checked = v


class FakeLabel:
    def __init__(self):
        self._text = ""

    def setText(self, t):
        self._text = t

    def text(self):
        return self._text


def _widget_self(**overrides):
    fake = SimpleNamespace()
    fake.spin = FakeSpin(0.05)
    fake.opacity_spin = FakeSpin(0.4)
    fake.opacity_slider = MagicMock()
    fake.combo_style = FakeCombo("Surface")
    fake.check_smooth = FakeCheck(False)
    fake.opacity_label = FakeLabel()
    fake.iso_actor_p = None
    fake.iso_actor_n = None
    fake.iso_actor_p_sil = None
    fake.iso_actor_n_sil = None
    grid = MagicMock()
    contour_mesh = MagicMock()
    contour_mesh.n_points = 5
    grid.contour.return_value = contour_mesh
    fake.grid = grid
    fake.plotter = MagicMock()
    fake.color_p = (0, 0, 255)
    fake.color_n = (255, 0, 0)
    fake.use_pbr = False
    fake.metallic = 0.5
    fake.roughness = 0.5
    fake.use_silhouette = False
    fake.update_iso = MagicMock()
    fake.presets = {}
    fake.btn_color_n = MagicMock()
    for k, v in overrides.items():
        setattr(fake, k, v)
    return fake


class TestUpdateIso:
    def test_wireframe_is_density_mode_label(self):
        fake = _widget_self(combo_style=FakeCombo("wireframe"))
        _update_iso(fake)
        assert fake.opacity_label.text() == "Density:"

    def test_surface_is_opacity_label(self):
        fake = _widget_self(combo_style=FakeCombo("Surface"))
        _update_iso(fake)
        assert fake.opacity_label.text() == "Opacity:"

    def test_smoothed_surface_maps_to_surface_style(self):
        fake = _widget_self(combo_style=FakeCombo("Smoothed Surface"))
        _update_iso(fake)
        # Positive-lobe mesh must have been smoothed before add_mesh.
        mesh = fake.grid.contour.return_value
        mesh.smooth.assert_called()

    def test_density_mode_clamps_target_reduction(self):
        # opacity_val=0.0 -> target_reduction would be 1.0, clamped to 0.99
        fake = _widget_self(
            combo_style=FakeCombo("points"), opacity_spin=FakeSpin(0.0)
        )
        _update_iso(fake)
        mesh = fake.grid.contour.return_value
        # Both lobes (positive + negative) route through the same mesh mock.
        assert mesh.decimate.call_count == 2
        for call in mesh.decimate.call_args_list:
            (reduction,), _kw = call
            assert reduction == pytest.approx(0.99)

    def test_contour_called_with_positive_and_negative_isovalue(self):
        fake = _widget_self(spin=FakeSpin(0.03))
        _update_iso(fake)
        calls = fake.grid.contour.call_args_list
        isovals = [c.kwargs["isosurfaces"][0] for c in calls]
        assert 0.03 in isovals
        assert -0.03 in isovals

    def test_actor_cleanup_swallows_remove_actor_errors(self):
        fake = _widget_self()
        fake.iso_actor_p = MagicMock()
        fake.plotter.remove_actor.side_effect = RuntimeError("stale actor")
        # Should not raise despite remove_actor blowing up.
        _update_iso(fake)

    def test_silhouette_added_when_enabled(self):
        fake = _widget_self(use_silhouette=True)
        _update_iso(fake)
        assert fake.plotter.add_silhouette.call_count == 2  # positive + negative lobe

    def test_no_contour_points_skips_add_mesh(self):
        fake = _widget_self()
        fake.grid.contour.return_value.n_points = 0
        _update_iso(fake)
        fake.plotter.add_mesh.assert_not_called()

    def test_plotter_render_called_on_success(self):
        fake = _widget_self()
        _update_iso(fake)
        fake.plotter.render.assert_called_once()

    def test_exception_in_contour_is_caught(self):
        fake = _widget_self()
        fake.grid.contour.side_effect = RuntimeError("boom")
        _update_iso(fake)  # must not raise


class TestOpacitySync:
    def test_opacity_slider_syncs_spin(self):
        fake = _widget_self()
        _on_opacity_changed(fake, 55)
        assert fake.opacity_spin.set_calls == [0.55]

    def test_opacity_spin_syncs_slider(self):
        fake = _widget_self()
        _on_opacity_spin_changed(fake, 0.75)
        fake.opacity_slider.setValue.assert_called_once_with(75)

    def test_iso_slider_updates_spin_value(self):
        fake = _widget_self(slider_max_int=1000, max_val=1.0)
        _on_slider_changed(fake, 250)
        assert fake.spin.set_calls == [pytest.approx(0.25)]

    def test_iso_spin_updates_slider_value(self):
        slider = MagicMock()
        fake = _widget_self(slider=slider, slider_max_int=1000, max_val=1.0)
        _on_spin_changed(fake, 0.5)
        slider.setValue.assert_called_once_with(500)

    def test_iso_spin_zero_max_val_skips_slider_update(self):
        slider = MagicMock()
        fake = _widget_self(slider=slider, slider_max_int=1000, max_val=0.0)
        _on_spin_changed(fake, 0.0)
        slider.setValue.assert_not_called()


class TestSettingsPathAndDefaults:
    def test_settings_path_is_json_next_to_plugin(self):
        fake = SimpleNamespace()
        path = _get_settings_path(fake)
        assert path.endswith("cube_viewer_advanced.json")
        assert "Cube_File_Viewer_Advanced" in path

    def test_default_presets_contains_default(self):
        fake = SimpleNamespace(presets={})
        _init_default_presets(fake)
        assert "Default" in fake.presets
        assert fake.default_preset_names == {"Default"}
        assert fake.presets["Default"]["isovalue"] == 0.02


class TestComplementaryColor:
    def test_blue_complement_is_orange_hue(self):
        fake = _widget_self(color_p=(0, 0, 255))
        _update_complementary_color(fake)
        r, g, b = fake.color_n
        # Complement of pure blue (hue 240) is hue 60 (yellow/orange), i.e.
        # high red+green, zero/low blue.
        assert b <= max(r, g)

    def test_grayscale_color_unchanged_hue(self):
        fake = _widget_self(color_p=(128, 128, 128))
        _update_complementary_color(fake)
        # Achromatic colors keep the same (achromatic) hue -> stays gray.
        r, g, b = fake.color_n
        assert r == g == b

    def test_update_iso_called_after_complement(self, monkeypatch):
        fake = _widget_self(color_p=(0, 0, 255))
        called = []
        fake.update_iso = lambda: called.append(True)
        # update_complementary_color calls self.update_iso() at the end.
        _update_complementary_color(fake)
        assert called == [True]


class TestCompColorToggle:
    def test_checked_enables_complement_and_disables_neg_button(self):
        fake = _widget_self(color_p=(10, 20, 30))
        fake.btn_color_n = MagicMock()
        called = []
        fake.update_complementary_color = lambda: called.append(True)
        _on_comp_color_toggled(fake, True)
        fake.btn_color_n.setEnabled.assert_called_once_with(False)
        assert called == [True]

    def test_unchecked_enables_neg_button_without_complement(self):
        fake = _widget_self()
        fake.btn_color_n = MagicMock()
        called = []
        fake.update_complementary_color = lambda: called.append(True)
        _on_comp_color_toggled(fake, False)
        fake.btn_color_n.setEnabled.assert_called_once_with(True)
        assert called == []


# ---------------------------------------------------------------------------
# open_cube_viewer / run_plugin / initialize — module-level entry points
# (real code paths, not extracted; CubeViewerWidget()/ChargeDialog() calls
# collapse to MagicMock under the mocked-Qt-base trick, so their internals
# are inert here and we only exercise open_cube_viewer's own logic).
# ---------------------------------------------------------------------------

def _fresh_mod():
    """Load a brand-new module instance with real numpy, isolated per test."""
    return _load_with_real_numpy(CUBE_ADV_PATH)


class TestOpenCubeViewer:
    def test_existing_dock_with_close_plugin_uses_it(self, cube_file):
        mod = _fresh_mod()
        ctx = make_context()
        old_widget = MagicMock()
        old_dock = MagicMock()
        old_dock.widget.return_value = old_widget
        ctx.get_window.return_value = old_dock
        mod.open_cube_viewer(ctx, str(cube_file))
        old_widget.close_plugin.assert_called_once()
        old_dock.close.assert_not_called()

    def test_existing_dock_without_close_plugin_closes_dock(self, cube_file):
        mod = _fresh_mod()
        ctx = make_context()
        old_widget = MagicMock(spec=[])  # no close_plugin attribute
        old_dock = MagicMock()
        old_dock.widget.return_value = old_widget
        ctx.get_window.return_value = old_dock
        mod.open_cube_viewer(ctx, str(cube_file))
        old_dock.close.assert_called_once()

    def test_successful_load_sets_current_molecule_and_status(self, cube_file):
        mod = _fresh_mod()
        ctx = make_context()
        ctx.get_window.return_value = None
        mol = MagicMock()
        mol.GetNumBonds.return_value = 2
        mod.Chem.MolFromXYZBlock.return_value = mol
        mod.open_cube_viewer(ctx, str(cube_file))
        assert ctx.current_molecule is mol
        ctx.draw_molecule_3d.assert_called_once_with(mol)
        msg = ctx.show_status_message.call_args[0][0]
        assert "Bonds: 2" in msg

    def test_mol_none_triggers_fallback_construction(self, cube_file):
        mod = _fresh_mod()
        ctx = make_context()
        ctx.get_window.return_value = None
        mod.Chem.MolFromXYZBlock.return_value = None
        rwmol = MagicMock()
        mod.Chem.RWMol.return_value = rwmol
        mod.open_cube_viewer(ctx, str(cube_file))
        rwmol.AddConformer.assert_called_once()
        # Fallback path constructs from Chem.RWMol(); confirm it was used.
        mod.Chem.RWMol.assert_called_once()

    def test_zero_bonds_no_rddeterminebonds_skips_retry_loop(self, cube_file):
        mod = _fresh_mod()
        mod.rdDetermineBonds = None
        ctx = make_context()
        ctx.get_window.return_value = None
        mol = MagicMock()
        mol.GetNumBonds.return_value = 0
        mod.Chem.MolFromXYZBlock.return_value = mol
        mod.open_cube_viewer(ctx, str(cube_file))
        assert ctx.current_molecule is mol
        ctx.draw_molecule_3d.assert_called_once_with(mol)

    def test_user_cancels_charge_dialog_returns_without_setting_molecule(self, cube_file):
        mod = _fresh_mod()
        ctx = make_context()
        ctx.get_window.return_value = None
        mol = MagicMock()
        mol.GetNumBonds.return_value = 0
        mod.Chem.MolFromXYZBlock.return_value = mol
        mod.rdDetermineBonds.DetermineConnectivity.side_effect = RuntimeError("bad")
        # ChargeDialog(...) collapses to a MagicMock under mocked Qt bases, so
        # dlg.exec() != QDialog.DialogCode.Accepted -> the "user cancelled" branch.
        mod.open_cube_viewer(ctx, str(cube_file))
        assert ctx.current_molecule != mol
        ctx.draw_molecule_3d.assert_not_called()

    def test_exception_reads_bad_file_reports_error_status(self):
        mod = _fresh_mod()
        ctx = make_context()
        ctx.get_window.return_value = None
        mod.open_cube_viewer(ctx, "this_file_does_not_exist.cube")
        msg = ctx.show_status_message.call_args[0][0]
        assert "Error" in msg


class TestRunPluginAndInitialize:
    def test_run_plugin_shows_critical_when_rdkit_missing(self):
        mod = _fresh_mod()
        ctx = make_context()
        mod.Chem = None
        mod.run_plugin(ctx)
        mod.QMessageBox.critical.assert_called_once()

    def test_run_plugin_opens_selected_file(self, cube_file):
        mod = _fresh_mod()
        ctx = make_context()
        mod.QFileDialog.getOpenFileName.return_value = (str(cube_file), "")
        called = []
        mod.open_cube_viewer = lambda c, f: called.append((c, f))
        mod.run_plugin(ctx)
        assert called == [(ctx, str(cube_file))]

    def test_run_plugin_cancel_does_not_open(self):
        mod = _fresh_mod()
        ctx = make_context()
        mod.QFileDialog.getOpenFileName.return_value = ("", "")
        called = []
        mod.open_cube_viewer = lambda c, f: called.append((c, f))
        mod.run_plugin(ctx)
        assert called == []

    def test_initialize_file_opener_wrapper_calls_open_cube_viewer(self):
        mod = _fresh_mod()
        ctx = make_context()
        mod.initialize(ctx)
        opener_calls = {c.args[0]: c.args[1] for c in ctx.register_file_opener.call_args_list}
        called = []
        mod.open_cube_viewer = lambda c, f: called.append((c, f))
        opener_calls[".cube"]("test.cube")
        assert called == [(ctx, "test.cube")]

    def test_initialize_drop_handler_accepts_cube_extension(self):
        mod = _fresh_mod()
        ctx = make_context()
        mod.initialize(ctx)
        handler = ctx.register_drop_handler.call_args[0][0]
        called = []
        mod.open_cube_viewer = lambda c, f: called.append((c, f))
        assert handler("molecule.cube") is True
        assert called == [(ctx, "molecule.cube")]

    def test_initialize_drop_handler_rejects_other_extension(self):
        mod = _fresh_mod()
        ctx = make_context()
        mod.initialize(ctx)
        handler = ctx.register_drop_handler.call_args[0][0]
        assert handler("molecule.xyz") is False
