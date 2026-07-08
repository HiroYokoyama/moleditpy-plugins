"""
Tests for the Cube File Viewer plugin (plugins/Cube_File_Viewer/cube_viewer.py).

Pure-function tests (parse_cube_data) run with real numpy injected after
plugin load so we can verify actual parsing behaviour.
Qt / PyVista / RDKit remain mocked throughout.
"""
from __future__ import annotations

from pathlib import Path

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
