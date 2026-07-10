"""
Tests for the Mapped Cube Viewer plugin
(plugins/Mapped_Cube_Viewer/mapped_cube_viewer.py).

Pure-function tests (parse_cube_data) run with real numpy injected after
plugin load so we can verify actual parsing behaviour.
Qt / PyVista / RDKit remain mocked throughout.
"""
from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest

from conftest import extract_function, load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
MAPPED_PATH = PLUGINS_DIR / "Mapped_Cube_Viewer" / "mapped_cube_viewer.py"


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
def mapped_mod():
    return _load_with_real_numpy(MAPPED_PATH)


@pytest.fixture
def cube_file(tmp_path):
    """Minimal valid 2×2×2 cube file with one Carbon atom."""
    p = tmp_path / "test.cube"
    p.write_text(_make_cube_content(nx=2, ny=2, nz=2, n_atoms=1))
    return p


# ---------------------------------------------------------------------------
# parse_cube_data — Mapped Cube Viewer (independent copy)
# ---------------------------------------------------------------------------

class TestMappedCubeParseCubeData:
    def test_short_file_raises(self, tmp_path, mapped_mod):
        f = tmp_path / "short.cube"
        f.write_text("x\ny\nz\n")
        with pytest.raises(ValueError, match="too short"):
            mapped_mod.parse_cube_data(str(f))

    def test_returns_expected_keys(self, cube_file, mapped_mod):
        result = mapped_mod.parse_cube_data(str(cube_file))
        assert {"atoms", "dims", "data_flat", "is_angstrom_header"}.issubset(
            result.keys()
        )

    def test_dims(self, cube_file, mapped_mod):
        result = mapped_mod.parse_cube_data(str(cube_file))
        assert result["dims"] == (2, 2, 2)

    def test_data_length(self, cube_file, mapped_mod):
        result = mapped_mod.parse_cube_data(str(cube_file))
        nx, ny, nz = result["dims"]
        assert len(result["data_flat"]) == nx * ny * nz


# ---------------------------------------------------------------------------
# Mapped Cube Viewer — entry points
# ---------------------------------------------------------------------------

class TestMappedCubeEntryPoints:
    def test_has_run_function(self):
        with mock_optional_imports():
            mod = load_plugin(MAPPED_PATH)
        assert callable(getattr(mod, "run", None))

    def test_has_run_plugin_function(self):
        with mock_optional_imports():
            mod = load_plugin(MAPPED_PATH)
        assert callable(getattr(mod, "run_plugin", None))

    def test_no_initialize_function(self):
        """Mapped Cube Viewer uses legacy run() only — no initialize()."""
        with mock_optional_imports():
            mod = load_plugin(MAPPED_PATH)
        assert not hasattr(mod, "initialize")


# ---------------------------------------------------------------------------
# parse_cube_data — additional edge cases
# ---------------------------------------------------------------------------


class TestParseCubeDataEdgeCases:
    def test_negative_atom_count_skips_extra_header_line(self, tmp_path, mapped_mod):
        """A negative n_atoms in the header (MO-info cube) is preceded by an
        extra metadata line before the atom records begin."""
        content = _make_cube_content(nx=2, ny=1, nz=1, n_atoms=1)
        lines = content.splitlines(keepends=True)
        # lines[2] holds "   1   0.0 0.0 0.0\n" -> make it negative and insert
        # a non-5-token MO-info line right after the vector headers (line 6).
        lines[2] = lines[2].replace("   1", "  -1", 1)
        lines.insert(6, "    1    1\n")  # 2 tokens => triggers the skip branch
        f = tmp_path / "neg.cube"
        f.write_text("".join(lines))
        result = mapped_mod.parse_cube_data(str(f))
        assert len(result["atoms"]) == 1
        assert result["atoms"][0][0] == 6

    def test_data_trimmed_when_excess_values(self, tmp_path, mapped_mod):
        # The data-start heuristic requires >=6 tokens on the first data
        # line, so the first row alone already supplies 6 values.
        content = _make_cube_content(
            nx=2, ny=1, nz=1, n_atoms=1, data_values=[9, 9, 9, 9, 9, 9, 1, 2]
        )
        f = tmp_path / "excess.cube"
        f.write_text(content)
        result = mapped_mod.parse_cube_data(str(f))
        # expected_size = 2; trimmed from the start -> last 2 values kept
        assert list(result["data_flat"]) == [1.0, 2.0]

    def test_data_padded_when_missing_values(self, tmp_path, mapped_mod):
        # expected_size = 9 but only one 6-value row is supplied.
        content = _make_cube_content(
            nx=3, ny=3, nz=1, n_atoms=1, data_values=[1, 2, 3, 4, 5, 6]
        )
        f = tmp_path / "short_data.cube"
        f.write_text(content)
        result = mapped_mod.parse_cube_data(str(f))
        assert len(result["data_flat"]) == 9
        assert list(result["data_flat"]) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0]

    def test_malformed_atom_coords_falls_back_to_origin(self, tmp_path, mapped_mod):
        header = (
            "T1\nT2\n   1   0.0 0.0 0.0\n   1   0.2 0.0 0.0\n"
            "   1   0.0 0.2 0.0\n   1   0.0 0.0 0.2\n"
        )
        bad_atom_line = "   6   0.0   bad   bad\n"
        data = " 1.0\n"
        f = tmp_path / "bad_atom.cube"
        f.write_text(header + bad_atom_line + data)
        result = mapped_mod.parse_cube_data(str(f))
        assert result["atoms"][0][0] == 6
        assert list(result["atoms"][0][1]) == [0.0, 0.0, 0.0]

    def test_positive_dims_means_not_angstrom_header(self, cube_file, mapped_mod):
        result = mapped_mod.parse_cube_data(str(cube_file))
        assert result["is_angstrom_header"] is False

    def test_negative_dims_means_angstrom_header(self, tmp_path, mapped_mod):
        content = _make_cube_content(nx=2, ny=1, nz=1, n_atoms=1)
        lines = content.splitlines(keepends=True)
        lines[3] = lines[3].replace("   2", "  -2", 1)
        f = tmp_path / "angstrom.cube"
        f.write_text("".join(lines))
        result = mapped_mod.parse_cube_data(str(f))
        assert result["is_angstrom_header"] is True
        assert result["dims"] == (2, 1, 1)


# ---------------------------------------------------------------------------
# build_grid_from_meta / read_cube
# ---------------------------------------------------------------------------


class _FakeGrid:
    def __init__(self):
        self.points = None
        self.dimensions = None
        self.point_data = {}


def _install_fake_structured_grid(mod):
    mod.pv.StructuredGrid = lambda: _FakeGrid()


def _make_meta(is_angstrom_header):
    return {
        "atoms": [(6, np.array([1.0, 0.0, 0.0]))],
        "origin": np.array([0.0, 0.0, 0.0]),
        "x_vec": np.array([1.0, 0.0, 0.0]),
        "y_vec": np.array([0.0, 1.0, 0.0]),
        "z_vec": np.array([0.0, 0.0, 1.0]),
        "dims": (2, 2, 2),
        "data_flat": np.arange(8, dtype=float),
        "is_angstrom_header": is_angstrom_header,
    }


BOHR_TO_ANGSTROM = 0.529177210903


class TestBuildGridFromMeta:
    def test_converts_bohr_to_angstrom_by_default(self, mapped_mod):
        _install_fake_structured_grid(mapped_mod)
        meta = _make_meta(is_angstrom_header=False)
        atoms_dict, grid = mapped_mod.build_grid_from_meta(meta)
        expected_pos = np.array([1.0, 0.0, 0.0]) * BOHR_TO_ANGSTROM
        assert np.allclose(atoms_dict["atoms"][0][1], expected_pos)

    def test_leaves_angstrom_header_unconverted(self, mapped_mod):
        _install_fake_structured_grid(mapped_mod)
        meta = _make_meta(is_angstrom_header=True)
        atoms_dict, grid = mapped_mod.build_grid_from_meta(meta)
        assert np.allclose(atoms_dict["atoms"][0][1], [1.0, 0.0, 0.0])

    def test_grid_dimensions_match_meta_dims(self, mapped_mod):
        _install_fake_structured_grid(mapped_mod)
        meta = _make_meta(is_angstrom_header=True)
        _, grid = mapped_mod.build_grid_from_meta(meta)
        assert grid.dimensions == [2, 2, 2]

    def test_grid_points_span_origin_to_far_corner(self, mapped_mod):
        _install_fake_structured_grid(mapped_mod)
        meta = _make_meta(is_angstrom_header=True)  # no scaling, easy to check
        _, grid = mapped_mod.build_grid_from_meta(meta)
        assert grid.points.shape == (8, 3)
        assert np.allclose(grid.points.min(axis=0), [0.0, 0.0, 0.0])
        assert np.allclose(grid.points.max(axis=0), [1.0, 1.0, 1.0])

    def test_property_values_mapping_matches_reshape_flatten_contract(self, mapped_mod):
        """Regression guard for the documented "IMPORTANT FIX": data must be
        reshaped C-order then re-flattened F-order to align with the
        X-fast point ordering generated above."""
        _install_fake_structured_grid(mapped_mod)
        meta = _make_meta(is_angstrom_header=True)
        _, grid = mapped_mod.build_grid_from_meta(meta)
        nx, ny, nz = meta["dims"]
        expected = meta["data_flat"].reshape((nx, ny, nz), order="C").flatten(order="F")
        assert np.allclose(grid.point_data["property_values"], expected)

    def test_read_cube_combines_parse_and_build(self, tmp_path, mapped_mod):
        _install_fake_structured_grid(mapped_mod)
        f = tmp_path / "rc.cube"
        f.write_text(_make_cube_content(nx=2, ny=1, nz=1, n_atoms=1))
        atoms_dict, grid = mapped_mod.read_cube(str(f))
        assert grid.dimensions == [2, 1, 1]
        assert len(atoms_dict["atoms"]) == 1


# ---------------------------------------------------------------------------
# MappedCubeSetupDialog.accept() — file existence validation
# ---------------------------------------------------------------------------


def _dialog_accept_fn():
    globs = {
        "os": __import__("os"),
        "QMessageBox": MagicMock(),
        # accept() calls super().accept(); zero-arg super() needs a __class__
        # cell only available in a real class body, so stub the builtin.
        "super": lambda *a: SimpleNamespace(accept=lambda: None),
    }
    fn = extract_function(MAPPED_PATH, "MappedCubeSetupDialog", "accept", globs)
    return fn, globs["QMessageBox"]


class TestMappedCubeSetupDialogAccept:
    def test_rejects_when_both_files_missing(self, tmp_path):
        fn, qmb = _dialog_accept_fn()
        self_ = SimpleNamespace(
            le_surf=SimpleNamespace(text=lambda: str(tmp_path / "nope1.cube")),
            le_prop=SimpleNamespace(text=lambda: str(tmp_path / "nope2.cube")),
        )
        fn(self_)
        qmb.warning.assert_called_once()
        assert "Surface file not found" in qmb.warning.call_args[0][2]

    def test_rejects_when_property_file_missing(self, tmp_path):
        fn, qmb = _dialog_accept_fn()
        surf = tmp_path / "surf.cube"
        surf.write_text("x")
        self_ = SimpleNamespace(
            le_surf=SimpleNamespace(text=lambda: str(surf)),
            le_prop=SimpleNamespace(text=lambda: str(tmp_path / "nope.cube")),
        )
        fn(self_)
        qmb.warning.assert_called_once()
        assert "Property file not found" in qmb.warning.call_args[0][2]

    def test_accepts_when_both_files_exist(self, tmp_path):
        fn, qmb = _dialog_accept_fn()
        surf = tmp_path / "surf.cube"
        prop = tmp_path / "prop.cube"
        surf.write_text("x")
        prop.write_text("y")
        self_ = SimpleNamespace(
            le_surf=SimpleNamespace(text=lambda: str(surf)),
            le_prop=SimpleNamespace(text=lambda: str(prop)),
        )
        fn(self_)
        qmb.warning.assert_not_called()
        assert self_.surface_file == str(surf)
        assert self_.property_file == str(prop)


# ---------------------------------------------------------------------------
# MappedWidget.update_mesh — extracted, driven with stub grids
# ---------------------------------------------------------------------------


def _update_mesh_fn():
    globs = {"logging": logging}
    return extract_function(MAPPED_PATH, "MappedWidget", "update_mesh", globs)


class _FakeContourResult:
    def __init__(self, n_points, sample_result=None):
        self.n_points = n_points
        self._sample_result = sample_result

    def sample(self, other):
        return self._sample_result


class _FakeMapped:
    def __init__(self, values):
        self.point_data = {"property_values": np.array(values)}


def _widget_stub(iso_val=0.002, opacity=0.4, min_v=-0.1, max_v=0.1,
                  contour_n_points=10, mapped_values=None, actor=None):
    mapped = _FakeMapped(mapped_values if mapped_values is not None else [1.0, 2.0])
    contour_result = _FakeContourResult(contour_n_points, sample_result=mapped)
    grid_surf = SimpleNamespace(contour=lambda vals, scalars: contour_result)
    plotter = MagicMock()
    plotter.add_mesh.return_value = "ACTOR"
    return SimpleNamespace(
        iso_spin=SimpleNamespace(value=lambda: iso_val),
        opacity_spin=SimpleNamespace(value=lambda: opacity),
        grid_surf=grid_surf,
        grid_prop=SimpleNamespace(),
        min_spin=SimpleNamespace(
            value=lambda: min_v, setValue=MagicMock(), blockSignals=MagicMock()
        ),
        max_spin=SimpleNamespace(
            value=lambda: max_v, setValue=MagicMock(), blockSignals=MagicMock()
        ),
        cmap_combo=SimpleNamespace(currentText=lambda: "jet_r"),
        actor=actor,
        context=SimpleNamespace(plotter=plotter),
    ), plotter


class TestUpdateMesh:
    def test_empty_contour_is_noop(self):
        fn = _update_mesh_fn()
        stub, plotter = _widget_stub(contour_n_points=0)
        fn(stub)
        plotter.add_mesh.assert_not_called()

    def test_adds_mesh_with_manual_clim(self):
        fn = _update_mesh_fn()
        stub, plotter = _widget_stub(min_v=-5.0, max_v=5.0)
        fn(stub, auto_fit=False)
        plotter.add_mesh.assert_called_once()
        kwargs = plotter.add_mesh.call_args.kwargs
        assert kwargs["clim"] == [-5.0, 5.0]
        assert kwargs["cmap"] == "jet_r"
        assert stub.actor == "ACTOR"
        plotter.render.assert_called_once()

    def test_removes_previous_actor(self):
        fn = _update_mesh_fn()
        stub, plotter = _widget_stub(actor="OLD_ACTOR")
        fn(stub, auto_fit=False)
        plotter.remove_actor.assert_called_once_with("OLD_ACTOR")

    def test_auto_fit_recomputes_clim_from_mapped_values(self):
        fn = _update_mesh_fn()
        stub, plotter = _widget_stub(mapped_values=[3.0, -2.0, 7.0])
        fn(stub, auto_fit=True)
        kwargs = plotter.add_mesh.call_args.kwargs
        assert kwargs["clim"] == [-2.0, 7.0]
        stub.min_spin.setValue.assert_called_once_with(-2.0)
        stub.max_spin.setValue.assert_called_once_with(7.0)

    def test_auto_fit_widens_degenerate_range(self):
        fn = _update_mesh_fn()
        stub, plotter = _widget_stub(mapped_values=[4.0, 4.0, 4.0])
        fn(stub, auto_fit=True)
        kwargs = plotter.add_mesh.call_args.kwargs
        assert kwargs["clim"][0] == 4.0
        assert kwargs["clim"][1] == pytest.approx(4.001)

    def test_exception_in_contour_is_caught(self):
        fn = _update_mesh_fn()
        stub, plotter = _widget_stub()

        def _boom(vals, scalars):
            raise RuntimeError("boom")

        stub.grid_surf.contour = _boom
        fn(stub)  # must not raise
        plotter.add_mesh.assert_not_called()


# ---------------------------------------------------------------------------
# MappedWidget.close_plugin
# ---------------------------------------------------------------------------


def _close_plugin_fn():
    globs = {"logging": logging}
    return extract_function(MAPPED_PATH, "MappedWidget", "close_plugin", globs)


class TestClosePlugin:
    def test_removes_actor_closes_dock_and_self(self):
        fn = _close_plugin_fn()
        plotter = MagicMock()
        dock = MagicMock()
        mw = SimpleNamespace(
            ui_manager=SimpleNamespace(restore_ui_for_editing=MagicMock()),
            edit_actions_manager=MagicMock(),
        )
        stub = SimpleNamespace(
            context=SimpleNamespace(plotter=plotter),
            actor="ACTOR",
            dock=dock,
            mw=mw,
            close=MagicMock(),
        )
        fn(stub)
        plotter.remove_actor.assert_called_once_with("ACTOR")
        plotter.render.assert_called_once()
        mw.ui_manager.restore_ui_for_editing.assert_called_once()
        dock.close.assert_called_once()
        mw.edit_actions_manager.clear_all.assert_called_once()
        stub.close.assert_called_once()

    def test_missing_ui_manager_does_not_crash(self):
        fn = _close_plugin_fn()
        stub = SimpleNamespace(
            context=SimpleNamespace(plotter=MagicMock()),
            actor=None,
            dock=None,
            mw=SimpleNamespace(),  # no ui_manager, no edit_actions_manager
            close=MagicMock(),
        )
        fn(stub)  # must not raise
        stub.close.assert_called_once()
