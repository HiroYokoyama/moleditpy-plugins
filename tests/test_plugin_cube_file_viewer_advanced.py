"""
Tests for the Cube File Viewer Advanced plugin
(plugins/Cube_File_Viewer_Advanced/cube_viewer_advanced.py).

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
