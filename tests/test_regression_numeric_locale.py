"""Coordinates must survive a comma-decimal numeric locale (HiroYokoyama/python_molecular_editor#185).

RDKit's Chem.MolFromXYZBlock parses numbers through the C runtime locale, so
under e.g. de_DE it truncates every coordinate to its integer part and bond
perception then invents bonds. The plugins below therefore build molecules
straight from numeric coordinates; these tests run that path with real RDKit,
switching to a comma-decimal locale when the system has one.
"""

import contextlib
import locale
import sys
from pathlib import Path

import numpy as np
import pytest
from conftest import extract_function, load_plugin, mocks_with_real_numpy

Chem = pytest.importorskip("rdkit.Chem")
Geometry = pytest.importorskip("rdkit.Geometry")
rdDetermineBonds = pytest.importorskip("rdkit.Chem.rdDetermineBonds")

ROOT = Path(__file__).resolve().parents[1] / "plugins"
CUBE_PATHS = [
    ROOT / "Cube_File_Viewer/cube_viewer.py",
    ROOT / "Cube_File_Viewer_Advanced/cube_viewer_advanced.py",
    ROOT / "Mapped_Cube_Viewer/mapped_cube_viewer.py",
]

# Ethanol, fractional parts on every coordinate.
ATOMS = [
    (6, (-0.0465, 0.5591, 0.0001)),
    (6, (1.2645, -0.2122, -0.0002)),
    (8, (-1.1104, -0.3812, 0.0001)),
    (1, (-0.1017, 1.1953, 0.8885)),
    (1, (-0.1017, 1.1957, -0.8880)),
    (1, (2.1053, 0.4867, -0.0004)),
    (1, (1.3352, -0.8494, -0.8862)),
    (1, (1.3352, -0.8490, 0.8860)),
    (1, (-1.9549, 0.0952, 0.0004)),
]
XYZ = "\n".join(
    [str(len(ATOMS)), "ethanol"]
    + [
        f"{Chem.GetPeriodicTable().GetElementSymbol(z)} {x:.4f} {y:.4f} {zz:.4f}"
        for z, (x, y, zz) in ATOMS
    ]
)


@contextlib.contextmanager
def comma_numeric_locale():
    """Switch LC_NUMERIC to a comma-decimal locale when one is installed."""
    previous = locale.setlocale(locale.LC_NUMERIC)
    for name in ("de_DE.UTF-8", "de_DE.utf8", "de_DE", "de-DE", "German", "fr_FR.UTF-8"):
        try:
            locale.setlocale(locale.LC_NUMERIC, name)
        except locale.Error:
            continue
        if locale.localeconv()["decimal_point"] == ",":
            break
    try:
        yield
    finally:
        locale.setlocale(locale.LC_NUMERIC, previous)


def assert_ethanol(mol):
    assert mol.GetNumAtoms() == len(ATOMS)
    positions = mol.GetConformer().GetPositions()
    np.testing.assert_allclose(positions, [p for _, p in ATOMS], atol=1e-9)
    work = Chem.RWMol(mol)
    rdDetermineBonds.DetermineConnectivity(work)
    assert work.GetNumBonds() == 8


def _load(path):
    with mocks_with_real_numpy():
        mod = load_plugin(path)
    # The loader mocks rdkit; swap in the real modules for these paths.
    mod.Chem, mod.Geometry = Chem, Geometry
    return mod


@pytest.mark.parametrize("path", CUBE_PATHS, ids=lambda p: p.parent.name)
def test_cube_viewers_keep_fractional_coordinates(path):
    mod = _load(path)
    with comma_numeric_locale():
        assert_ethanol(mod.mol_from_cube_atoms(ATOMS))


def test_gaussian_mo_reader_keeps_fractional_coordinates():
    with mocks_with_real_numpy():
        mod = load_plugin(ROOT / "Gaussian_MO_Analyzer/gaussian_fchk_mo_analyzer/analyzer.py")
    reader = mod.FCHKReader.__new__(mod.FCHKReader)
    reader.data = {
        "Atomic numbers": np.array([z for z, _ in ATOMS]),
        "Current cartesian coordinates": (
            np.array([p for _, p in ATOMS]) / 0.529177249
        ).ravel(),
    }
    with comma_numeric_locale():
        assert_ethanol(reader.get_rdkit_mol())


def test_animated_xyz_giffer_frame_keeps_fractional_coordinates():
    # AnimatedXYZPlayer derives from a mocked QDialog, so take the method by AST.
    build = extract_function(
        ROOT / "Animated_XYZ_Giffer/animated_xyz_giffer.py",
        "AnimatedXYZPlayer",
        "_frame_atoms_mol",
        extra_globals={"Chem": Chem, "rdGeometry": Geometry},
    )
    frame = {
        "symbols": [Chem.GetPeriodicTable().GetElementSymbol(z) for z, _ in ATOMS],
        "coords": [p for _, p in ATOMS],
    }
    with comma_numeric_locale():
        assert_ethanol(build(None, frame))


def test_molecule_comparator_xyz_text_keeps_fractional_coordinates():
    mod = _load(ROOT / "Molecule_Comparator/molecule_comparator.py")
    mod.Point3D = Geometry.Point3D
    mod.rdDetermineBonds = rdDetermineBonds
    with comma_numeric_locale():
        assert_ethanol(mod.mol_from_xyz_text(XYZ))


def test_rdkit_xyz_reader_is_locale_sensitive_here():
    """Documents why the plugins avoid MolFromXYZBlock; skips where untestable."""
    with comma_numeric_locale():
        if locale.localeconv()["decimal_point"] != ",":
            pytest.skip("no comma-decimal locale installed")
        mol = Chem.MolFromXYZBlock(XYZ)
    if mol is not None and abs(mol.GetConformer().GetAtomPosition(1).x - 1.2645) < 1e-6:
        pytest.skip(f"this RDKit build parses XYZ locale-independently ({sys.platform})")
    assert mol is None or mol.GetConformer().GetAtomPosition(1).x != pytest.approx(1.2645)
