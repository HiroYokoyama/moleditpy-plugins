"""
Tests for the Mapped Cube Viewer plugin
(plugins/Mapped_Cube_Viewer/mapped_cube_viewer.py).

Pure-function tests (parse_cube_data) run with real numpy injected after
plugin load so we can verify actual parsing behaviour.
Qt / PyVista / RDKit remain mocked throughout.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from conftest import load_plugin, mock_optional_imports

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
