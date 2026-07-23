"""
Tests for the Cube File Viewer plugin (plugins/Cube_File_Viewer/cube_viewer.py).

Pure-function tests (parse_cube_data) run with real numpy injected after
plugin load so we can verify actual parsing behaviour.
Qt / PyVista / RDKit remain mocked throughout.
"""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CUBE_PATH = PLUGINS_DIR / "Cube_File_Viewer" / "cube_viewer.py"


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
def cube_mod():
    return _load_with_real_numpy(CUBE_PATH)


@pytest.fixture
def cube_file(tmp_path):
    """Minimal valid 2×2×2 cube file with one Carbon atom."""
    p = tmp_path / "test.cube"
    p.write_text(_make_cube_content(nx=2, ny=2, nz=2, n_atoms=1))
    return p


# ---------------------------------------------------------------------------
# parse_cube_data — Cube File Viewer
# ---------------------------------------------------------------------------

class TestParseCubeData:
    def test_short_file_raises(self, tmp_path, cube_mod):
        f = tmp_path / "short.cube"
        f.write_text("line1\nline2\nline3\n")
        with pytest.raises(ValueError, match="too short"):
            cube_mod.parse_cube_data(str(f))

    def test_returns_expected_keys(self, cube_file, cube_mod):
        result = cube_mod.parse_cube_data(str(cube_file))
        required = {
            "atoms", "origin", "x_vec", "y_vec", "z_vec",
            "dims", "data_flat", "is_angstrom_header",
        }
        assert required.issubset(result.keys())

    def test_dims_parsed_correctly(self, cube_file, cube_mod):
        result = cube_mod.parse_cube_data(str(cube_file))
        assert result["dims"] == (2, 2, 2)

    def test_atom_count_and_atomic_number(self, cube_file, cube_mod):
        result = cube_mod.parse_cube_data(str(cube_file))
        assert len(result["atoms"]) == 1
        atomic_num, _pos = result["atoms"][0]
        assert atomic_num == 6  # Carbon

    def test_data_flat_length_equals_grid_size(self, cube_file, cube_mod):
        result = cube_mod.parse_cube_data(str(cube_file))
        nx, ny, nz = result["dims"]
        assert len(result["data_flat"]) == nx * ny * nz

    def test_bohr_header_flag_false_when_positive_n(self, cube_file, cube_mod):
        """Positive NX/NY/NZ → Bohr units → is_angstrom_header False."""
        result = cube_mod.parse_cube_data(str(cube_file))
        assert result["is_angstrom_header"] is False

    def test_negative_nx_sets_angstrom_header(self, tmp_path, cube_mod):
        """Negative NX in the file signals Angstrom units."""
        content = _make_cube_content(nx=2, ny=2, nz=2)
        # Flip sign of NX line
        content = content.replace(
            "   2   0.200000   0.000000   0.000000",
            "  -2   0.200000   0.000000   0.000000",
            1,
        )
        p = tmp_path / "angstrom.cube"
        p.write_text(content)
        result = cube_mod.parse_cube_data(str(p))
        assert result["is_angstrom_header"] is True
        assert result["dims"] == (2, 2, 2)  # abs() must have been applied

    def test_pads_short_data_with_zeros(self, tmp_path, cube_mod):
        """data_flat is zero-padded at the end when fewer values than grid size."""
        # 3×3×3 = 27 expected; provide only 6
        content = _make_cube_content(nx=3, ny=3, nz=3, n_atoms=1,
                                     data_values=list(range(1, 7)))
        p = tmp_path / "short_data.cube"
        p.write_text(content)
        result = cube_mod.parse_cube_data(str(p))
        assert result["dims"] == (3, 3, 3)
        assert len(result["data_flat"]) == 27
        assert result["data_flat"][6] == pytest.approx(0.0)

    def test_trims_excess_data_from_start(self, tmp_path, cube_mod):
        """Excess data values are trimmed from the beginning of the array."""
        # 2×2×2 = 8 expected; provide 10 values [1…10]
        content = _make_cube_content(nx=2, ny=2, nz=2, n_atoms=1,
                                     data_values=list(range(1, 11)))
        p = tmp_path / "excess_data.cube"
        p.write_text(content)
        result = cube_mod.parse_cube_data(str(p))
        assert len(result["data_flat"]) == 8
        # First 2 values trimmed → array starts at 3.0
        assert result["data_flat"][0] == pytest.approx(3.0)

    def test_multiple_atoms(self, tmp_path, cube_mod):
        content = _make_cube_content(nx=2, ny=2, nz=2, n_atoms=3)
        p = tmp_path / "multi_atom.cube"
        p.write_text(content)
        result = cube_mod.parse_cube_data(str(p))
        assert len(result["atoms"]) == 3

    def test_origin_vector_parsed(self, cube_file, cube_mod):
        result = cube_mod.parse_cube_data(str(cube_file))
        origin = result["origin"]
        assert len(origin) == 3
        assert origin[0] == pytest.approx(0.0)

    def test_x_vec_parsed(self, cube_file, cube_mod):
        result = cube_mod.parse_cube_data(str(cube_file))
        x_vec = result["x_vec"]
        assert len(x_vec) == 3
        assert x_vec[0] == pytest.approx(0.2)


# ---------------------------------------------------------------------------
# initialize() — Cube File Viewer
# ---------------------------------------------------------------------------

class TestCubeViewerInitialize:
    def test_does_not_raise(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            ctx = make_context()
            mod.initialize(ctx)

    def test_registers_dot_cube_opener(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        exts = [c.args[0] for c in ctx.register_file_opener.call_args_list]
        assert ".cube" in exts

    def test_registers_dot_cub_opener(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        exts = [c.args[0] for c in ctx.register_file_opener.call_args_list]
        assert ".cub" in exts

    def test_registers_drop_handler(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        assert ctx.register_drop_handler.called
        _, kwargs = ctx.register_drop_handler.call_args
        assert kwargs.get("priority") == 10

    def test_drop_handler_handles_cube_file(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            mod.open_cube_viewer = MagicMock()
            handler = ctx.register_drop_handler.call_args[0][0]
            assert handler("/some/path/mo.CUBE") is True
            mod.open_cube_viewer.assert_called_once_with(ctx, "/some/path/mo.CUBE")

    def test_drop_handler_ignores_other_extension(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            mod.open_cube_viewer = MagicMock()
            handler = ctx.register_drop_handler.call_args[0][0]
            assert handler("/some/path/mol.xyz") is False
            mod.open_cube_viewer.assert_not_called()

    def test_file_opener_wrapper_calls_open_cube_viewer(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            mod.open_cube_viewer = MagicMock()
            opener = ctx.register_file_opener.call_args_list[0].args[1]
            opener("/x/y.cube")
            mod.open_cube_viewer.assert_called_once_with(ctx, "/x/y.cube")


# ---------------------------------------------------------------------------
# RDKit truly unavailable — import-time fallback (Chem/Geometry/rdDetermineBonds
# set to None) and the plugin's user-facing guards for that state.
# ---------------------------------------------------------------------------

class TestRDKitUnavailable:
    def _load_without_rdkit(self):
        with mock_optional_imports():
            sys.modules["rdkit"] = None
            try:
                mod = load_plugin(CUBE_PATH)
            finally:
                del sys.modules["rdkit"]
            return mod

    def test_chem_is_none(self):
        mod = self._load_without_rdkit()
        assert mod.Chem is None
        assert mod.Geometry is None
        assert mod.rdDetermineBonds is None


    def test_open_cube_viewer_shows_error(self):
        mod = self._load_without_rdkit()
        ctx = make_context()
        mod.open_cube_viewer(ctx, "somefile.cube")
        mod.QMessageBox.critical.assert_called_once()

    def test_run_shows_error_and_no_file_dialog(self):
        mod = self._load_without_rdkit()
        mod.run(MagicMock())
        mod.QMessageBox.critical.assert_called_once()
        mod.QFileDialog.getOpenFileName.assert_not_called()


# ---------------------------------------------------------------------------
# run() — module-level PLUGIN_CONTEXT gate and file-dialog wiring
# ---------------------------------------------------------------------------

class TestRun:
    def test_run_returns_when_no_context(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            mod.PLUGIN_CONTEXT = None
            mod.QFileDialog.getOpenFileName = MagicMock()
            mod.run(MagicMock())
        mod.QFileDialog.getOpenFileName.assert_not_called()

    def test_run_skips_open_when_no_filename(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            mod.open_cube_viewer = MagicMock()
            mod.QFileDialog.getOpenFileName = MagicMock(return_value=("", ""))
            mod.run(MagicMock())
        mod.open_cube_viewer.assert_not_called()

    def test_run_opens_cube_when_filename_chosen(self):
        with mock_optional_imports():
            mod = load_plugin(CUBE_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            mod.open_cube_viewer = MagicMock()
            mod.QFileDialog.getOpenFileName = MagicMock(
                return_value=("/some/file.cube", "Cube Files (*.cube)")
            )
            mw = MagicMock()
            mod.run(mw)
        mod.open_cube_viewer.assert_called_once_with(ctx, "/some/file.cube")


# ---------------------------------------------------------------------------
# parse_cube_data — additional header/atom/data edge cases
# ---------------------------------------------------------------------------

class TestParseCubeDataEdgeCases:
    def test_negative_natoms_skips_non_atom_metadata_line(self, tmp_path, cube_mod):
        """n_atoms_raw<0 (MO cube) + a metadata line that isn't a 5-token
        atom line right after the vectors -> that line is skipped."""
        content = _make_cube_content(nx=2, ny=2, nz=2, n_atoms=1)
        content = content.replace(
            "   1   0.000000   0.000000   0.000000\n",
            "  -1   0.000000   0.000000   0.000000\n",
            1,
        )
        lines = content.splitlines(keepends=True)
        # Insert a non-atom metadata line ("1 5") right after the vector lines.
        lines.insert(6, "    1    5\n")
        p = tmp_path / "mo_meta.cube"
        p.write_text("".join(lines))
        result = cube_mod.parse_cube_data(str(p))
        assert len(result["atoms"]) == 1
        assert result["atoms"][0][0] == 6

    def test_negative_natoms_no_trailing_lines_hits_except(self, tmp_path, cube_mod):
        """n_atoms_raw<0 but the file ends right at the header -> the smart
        lookahead's IndexError is swallowed."""
        content = (
            "Title 1\nTitle 2\n"
            "  -1   0.0   0.0   0.0\n"
            "   1   0.2   0.0   0.0\n"
            "   1   0.0   0.2   0.0\n"
            "   1   0.0   0.0   0.2\n"
        )
        p = tmp_path / "truncated.cube"
        p.write_text(content)
        with pytest.raises(IndexError):
            cube_mod.parse_cube_data(str(p))

    def test_malformed_atom_line_falls_back_to_origin(self, tmp_path, cube_mod):
        content = _make_cube_content(nx=2, ny=2, nz=2, n_atoms=1)
        content = content.replace(
            "   6   0.000000   0.000000   0.000000   0.629118   0.000000\n",
            "   6   0.000000   bad   bad   bad\n",
            1,
        )
        p = tmp_path / "bad_atom.cube"
        p.write_text(content)
        result = cube_mod.parse_cube_data(str(p))
        atomic_num, pos = result["atoms"][0]
        assert atomic_num == 6
        assert list(pos) == [0.0, 0.0, 0.0]

    def test_skips_blank_and_short_and_nonnumeric_lines_before_data(
        self, tmp_path, cube_mod
    ):
        content = _make_cube_content(nx=2, ny=2, nz=2, n_atoms=1)
        header, atoms_and_data = content.split(
            "   6   0.000000   0.000000   0.000000   0.629118   0.000000\n", 1
        )
        content = (
            header
            + "   6   0.000000   0.000000   0.000000   0.629118   0.000000\n"
            + "\n"  # blank line
            + "   1   150\n"  # short metadata line (<6 tokens)
            + "comment word two three four five\n"  # non-numeric leading token, 6 tokens
            + atoms_and_data
        )
        p = tmp_path / "skippy.cube"
        p.write_text(content)
        result = cube_mod.parse_cube_data(str(p))
        assert len(result["data_flat"]) == 8

    def test_fromstring_exception_yields_zero_padded_array(
        self, cube_file, cube_mod, monkeypatch
    ):
        def _raiser(*_a, **_k):
            raise ValueError("boom")

        monkeypatch.setattr(cube_mod.np, "fromstring", _raiser)
        result = cube_mod.parse_cube_data(str(cube_file))
        nx, ny, nz = result["dims"]
        assert len(result["data_flat"]) == nx * ny * nz
        assert np.all(result["data_flat"] == 0.0)


# ---------------------------------------------------------------------------
# build_grid_from_meta / read_cube — grid reconstruction (real numpy, mocked
# pyvista — the plugin only assigns attributes on the StructuredGrid, so a
# MagicMock records them faithfully for inspection).
# ---------------------------------------------------------------------------

def _synthetic_meta(is_angstrom_header):
    nx = ny = nz = 2
    return {
        "atoms": [(6, np.array([1.0, 0.0, 0.0]))],
        "origin": np.array([0.0, 0.0, 0.0]),
        "x_vec": np.array([1.0, 0.0, 0.0]),
        "y_vec": np.array([0.0, 1.0, 0.0]),
        "z_vec": np.array([0.0, 0.0, 1.0]),
        "dims": (nx, ny, nz),
        "data_flat": np.arange(nx * ny * nz, dtype=float),
        "is_angstrom_header": is_angstrom_header,
    }


class TestBuildGridFromMeta:
    def test_bohr_header_converts_to_angstrom(self, cube_mod):
        meta = _synthetic_meta(is_angstrom_header=False)
        result, grid = cube_mod.build_grid_from_meta(meta)
        BOHR_TO_ANGSTROM = 0.529177210903
        assert result["atoms"][0][1][0] == pytest.approx(1.0 * BOHR_TO_ANGSTROM)
        # Max x-coordinate among the 8 points is 1 unit-cell step (converted).
        assert grid.points[:, 0].max() == pytest.approx(BOHR_TO_ANGSTROM)

    def test_angstrom_header_no_conversion(self, cube_mod):
        meta = _synthetic_meta(is_angstrom_header=True)
        result, grid = cube_mod.build_grid_from_meta(meta)
        assert result["atoms"][0][1][0] == pytest.approx(1.0)
        assert grid.points[:, 0].max() == pytest.approx(1.0)

    def test_grid_dimensions_set(self, cube_mod):
        meta = _synthetic_meta(is_angstrom_header=True)
        _, grid = cube_mod.build_grid_from_meta(meta)
        assert grid.dimensions == [2, 2, 2]

    def test_point_data_values_shape(self, cube_mod):
        meta = _synthetic_meta(is_angstrom_header=True)
        _, grid = cube_mod.build_grid_from_meta(meta)
        set_args = grid.point_data.__setitem__.call_args
        key, values = set_args.args
        assert key == "values"
        assert len(values) == 8

    def test_read_cube_end_to_end(self, cube_file, cube_mod):
        meta_result, grid = cube_mod.read_cube(str(cube_file))
        assert len(meta_result["atoms"]) == 1
        assert grid.dimensions == [2, 2, 2]
