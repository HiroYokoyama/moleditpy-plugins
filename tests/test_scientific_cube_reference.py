"""Independent coordinate/scalar oracle for all separately shipped Cube readers."""

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest
from conftest import load_plugin, mocks_with_real_numpy

ROOT = Path(__file__).resolve().parents[1] / "plugins"
PATHS = [
    ROOT / "Cube_File_Viewer/cube_viewer.py",
    ROOT / "Cube_File_Viewer_Advanced/cube_viewer_advanced.py",
    ROOT / "Mapped_Cube_Viewer/mapped_cube_viewer.py",
    ROOT / "Orbital_Comparator/orbital_comparator.py",
    ROOT / "Gaussian_MO_Analyzer/gaussian_fchk_mo_analyzer/vis.py",
]


class Grid:
    def __init__(self):
        self.point_data = {}


@pytest.mark.parametrize("path", PATHS, ids=lambda p: p.parent.name)
@pytest.mark.parametrize("angstrom", [False, True], ids=["bohr", "angstrom"])
@pytest.mark.parametrize("datasets", [1, 2], ids=["scalar", "interleaved-MOs"])
def test_skew_grid_matches_independent_affine_field(path, angstrom, datasets, tmp_path):
    # Unequal dimensions and skew axes distinguish every axis permutation.
    origin = np.array([0.3, -0.7, 1.1])
    vectors = np.array([[0.4, 0.1, 0.0], [-0.2, 0.5, 0.1], [0.1, -0.1, 0.6]])
    scale = 1.0 if angstrom else 0.529177210903
    dims = (2, 3, 4)
    lines = [
        "Independent affine field",
        "z fastest",
        "1 0.3 -0.7 1.1 2" if datasets == 2 else "1 0.3 -0.7 1.1",
    ]
    for count, vector in zip(dims, vectors):
        lines.append(f"{(-count if angstrom else count)} " + " ".join(map(str, vector)))
    lines.append("8 8 0.25 -0.5 0.75")
    for i in range(2):
        for j in range(3):
            for k in range(4):
                x, y, z = origin + i * vectors[0] + j * vectors[1] + k * vectors[2]
                value = 2 * x - 3 * y + 5 * z + 7
                lines.append(
                    f"{value:.15g}"
                    + (f" {10000 + value:.15g}" if datasets == 2 else "")
                )
    source = tmp_path / "reference.cube"
    source.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with mocks_with_real_numpy():
        module = load_plugin(path)
    module.np = np
    module.pv = SimpleNamespace(StructuredGrid=Grid)
    if path.name == "vis.py":
        viewer = module.CubeVisualizer(MagicMock())
        meta = viewer._parse_cube(source)
        grid = viewer._build_grid(meta)
    elif path.name == "orbital_comparator.py":
        meta = module.parse_cube(source)
        grid = module.build_grid(meta)
    else:
        meta = module.parse_cube_data(source)
        atoms, grid = module.build_grid_from_meta(meta)
        np.testing.assert_allclose(
            atoms["atoms"][0][1], [0.25 * scale, -0.5 * scale, 0.75 * scale], atol=1e-12
        )
    # Oracle associates each physical point with its analytic scalar; it never
    # duplicates the implementation's C/F reshape sequence.
    # Some legacy viewers use the 1986 Bohr constant (relative difference7e-8).
    actual_scale = 1.0 if angstrom else getattr(module, "BOHR_TO_ANG", scale)
    assert abs(actual_scale / scale - 1) < 1e-7
    points = np.asarray(grid.points) / actual_scale
    values = next(iter(grid.point_data.values()))
    np.testing.assert_allclose(
        values, points @ np.array([2.0, -3.0, 5.0]) + 7, atol=1e-12
    )
    expected_points = [
        origin + i * vectors[0] + j * vectors[1] + k * vectors[2]
        for i in range(2)
        for j in range(3)
        for k in range(4)
    ]
    np.testing.assert_allclose(
        sorted(map(tuple, np.round(points, 12))),
        sorted(map(tuple, np.round(expected_points, 12))),
        atol=1e-12,
    )
