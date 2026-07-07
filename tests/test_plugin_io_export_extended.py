"""
Extended tests for FILE-IO / EXPORT / UTILITY plugins:

  - Blender Export        (generate_blender_script content, double-bond offset)
  - POV-Ray Export        (generate_povray_scene content, resolution scaling,
                           triple-bond off_dir regression)
  - OpenBabel Conversion  (export_with_openbabel success status regression)
  - Structural Updater    (finalize() restore target regression, trigger routing)
  - Encrypted Project     (export/import round-trip file format, patched routing)
  - XYZ Editor            (_generate_xyz_content, _ClickFilter, get_mol_signature)
  - Python Console        (run_code execution / output capture)
  - Animated XYZ Giffer   (frame navigation, fps, status label)
  - Dummy Atom Mode       (patched Chem.Atom routing)
  - Paste from ChemDraw   (reconstructed MolBlock content)
  - Paste XYZ             (_resolve_mol_with_charge routing)

Script-generation tests run the real generator functions with the plugin's
heavy deps mocked, but with the REAL numpy injected into sys.modules (numpy is
installed on CI). Molecules are plain-Python fakes, so no chemistry libraries
are required.
"""

from __future__ import annotations

import ast
import contextlib
import sys
import textwrap
import types
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as real_numpy
import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

BLENDER_PATH = PLUGINS_DIR / "Blender_Export" / "blender_export.py"
POVRAY_PATH = PLUGINS_DIR / "POV-Ray_Export" / "povray_export.py"
OPENBABEL_PATH = PLUGINS_DIR / "OpenBabel_Conversion_Tool" / "openbabel_conversion_tool.py"
STRUCTURAL_PATH = PLUGINS_DIR / "Structural_Updater" / "structural_updater.py"
ENCRYPTED_PATH = PLUGINS_DIR / "Encrypted_Project" / "encrypted_project.py"
XYZ_EDITOR_PATH = PLUGINS_DIR / "XYZ_Editor" / "xyz_editor.py"
CONSOLE_PATH = PLUGINS_DIR / "Python_Console" / "console.py"
GIFFER_PATH = PLUGINS_DIR / "Animated_XYZ_Giffer" / "animated_xyz_giffer.py"
DUMMY_ATOM_PATH = PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py"
CHEMDRAW_PATH = PLUGINS_DIR / "Paste_from_ChemDraw" / "paste_chemdraw.py"
PASTE_XYZ_PATH = PLUGINS_DIR / "Paste_XYZ" / "paste_xyz.py"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def mocks_with_real_numpy():
    """
    Like mock_optional_imports(), but numpy resolves to the real package.

    The export generators do ``import numpy as np`` at call time; placing the
    real modules in sys.modules bypasses the mocking MetaPathFinder for numpy
    only, so vector math is real while rdkit/PyQt6 stay mocked.

    ALL already-imported ``numpy.*`` submodules must be restored, not just the
    top-level package: numpy resolves internals like ``numpy._core._methods``
    through sys.modules at call time, and if the MetaPathFinder answers that
    import with a MagicMock the real numpy package is corrupted for the rest
    of the process (e.g. ``np.allclose(0, 1)`` starts returning True).
    """
    real_mods = {
        k: v
        for k, v in sys.modules.items()
        if k == "numpy" or k.startswith("numpy.")
    }
    with mock_optional_imports():
        sys.modules.update(real_mods)
        try:
            yield
        finally:
            for k in real_mods:
                sys.modules.pop(k, None)


class P3:
    """3-component point supporting .x/.y/.z and the sequence protocol."""

    def __init__(self, x, y, z):
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __len__(self):
        return 3

    def __getitem__(self, i):
        return (self.x, self.y, self.z)[i]

    def __iter__(self):
        return iter((self.x, self.y, self.z))


class FakeAtom:
    def __init__(self, idx, symbol):
        self.idx = idx
        self.symbol = symbol
        self.bonds = []

    def GetIdx(self):
        return self.idx

    def GetSymbol(self):
        return self.symbol

    def GetDegree(self):
        return len(self.bonds)

    def GetBonds(self):
        return list(self.bonds)

    def GetNeighbors(self):
        out = []
        for b in self.bonds:
            out.append(b.end_atom if b.begin_atom is self else b.begin_atom)
        return out


class FakeBond:
    def __init__(self, begin_atom, end_atom, bond_type):
        self.begin_atom = begin_atom
        self.end_atom = end_atom
        self.bond_type = bond_type
        begin_atom.bonds.append(self)
        end_atom.bonds.append(self)

    def GetBeginAtomIdx(self):
        return self.begin_atom.idx

    def GetEndAtomIdx(self):
        return self.end_atom.idx

    def GetBondType(self):
        return self.bond_type


class FakeConf:
    def __init__(self, coords):
        self.coords = {i: P3(*c) for i, c in enumerate(coords)}

    def GetAtomPosition(self, idx):
        return self.coords[idx]


class FakeMol:
    def __init__(self, symbols, coords, bonds=()):
        # bonds: iterable of (begin_idx, end_idx, bond_type)
        self.atoms = [FakeAtom(i, s) for i, s in enumerate(symbols)]
        self.conf = FakeConf(coords)
        self.bonds = [FakeBond(self.atoms[b], self.atoms[e], t) for b, e, t in bonds]

    def GetNumAtoms(self):
        return len(self.atoms)

    def GetNumBonds(self):
        return len(self.bonds)

    def GetNumConformers(self):
        return 1

    def GetConformer(self):
        return self.conf

    def GetAtoms(self):
        return list(self.atoms)

    def GetBonds(self):
        return list(self.bonds)

    def GetAtomWithIdx(self, idx):
        return self.atoms[idx]


def make_export_mw(style="ball_and_stick", settings=None, color_map=None, plotter=None):
    """Main-window stub with real dict settings so .get() defaults apply."""
    return SimpleNamespace(
        init_manager=SimpleNamespace(settings=settings or {}, current_file_path=None),
        view_3d_manager=SimpleNamespace(
            _3d_color_map=color_map or {},
            current_3d_style=style,
            glyph_source=None,
        ),
        plotter=plotter,
    )


def _extract_fn(path: Path, class_name: str | None, fn_name: str, extra_globals=None):
    """
    Extract a function or method from source via AST and exec it standalone.

    Needed for methods of Qt-derived classes: the mocked Qt bases make the
    class object a MagicMock, so methods can't be reached via the class.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    scope = tree
    if class_name is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == class_name:
                scope = node
                break
        else:
            raise AssertionError(f"class {class_name} not found in {path}")
    for node in scope.body:
        if isinstance(node, ast.FunctionDef) and node.name == fn_name:
            segment = ast.get_source_segment(source, node)
            ns = {"MagicMock": MagicMock}
            ns.update(extra_globals or {})
            exec(textwrap.dedent(segment), ns)  # noqa: S102
            return ns[fn_name]
    raise AssertionError(f"function {fn_name} not found in {path}")


# ---------------------------------------------------------------------------
# Blender Export — _calculate_double_bond_offset
# ---------------------------------------------------------------------------


class TestBlenderDoubleBondOffset:
    def _offset(self, mol, bond):
        with mocks_with_real_numpy():
            mod = load_plugin(BLENDER_PATH)
            return mod._calculate_double_bond_offset(mol, bond, mol.GetConformer())

    def test_zero_length_bond_returns_z(self):
        mol = FakeMol(["C", "C"], [(0, 0, 0), (0, 0, 0)], [(0, 1, "D")])
        off = self._offset(mol, mol.bonds[0])
        assert list(off) == [0, 0, 1]

    def test_result_is_unit_vector(self):
        mol = FakeMol(["C", "C"], [(0, 0, 0), (1.4, 0, 0)], [(0, 1, "D")])
        off = self._offset(mol, mol.bonds[0])
        assert abs(real_numpy.linalg.norm(off) - 1.0) < 1e-9

    def test_result_perpendicular_to_bond(self):
        mol = FakeMol(["C", "C"], [(0, 0, 0), (0, 1.4, 0)], [(0, 1, "D")])
        off = self._offset(mol, mol.bonds[0])
        assert abs(real_numpy.dot(off, [0, 1, 0])) < 1e-9

    def test_neighbor_defines_plane(self):
        # Ethylene-like fragment in the XY plane: offset must stay in-plane.
        mol = FakeMol(
            ["C", "C", "H"],
            [(0, 0, 0), (1.33, 0, 0), (-0.5, 0.9, 0)],
            [(0, 1, "D"), (0, 2, "S")],
        )
        off = self._offset(mol, mol.bonds[0])
        assert abs(off[2]) < 1e-9
        assert abs(abs(off[1]) - 1.0) < 1e-9


# ---------------------------------------------------------------------------
# Blender Export — generate_blender_script content
# ---------------------------------------------------------------------------


def _blender_script(mol, mw):
    with mocks_with_real_numpy():
        mod = load_plugin(BLENDER_PATH)
        return mod.generate_blender_script(mol, mw), mod


def _simple_mol(bond_types):
    """Linear chain C-C(-C...) with the given rdkit-mocked bond types."""
    n = len(bond_types) + 1
    symbols = ["C"] * n
    coords = [(1.5 * i, 0.0, 0.0) for i in range(n)]
    bonds = [(i, i + 1, bond_types[i]) for i in range(len(bond_types))]
    return FakeMol(symbols, coords, bonds)


class TestBlenderScriptGeneration:
    def _generate(self, style="ball_and_stick", bond_kinds=("SINGLE",), color_map=None):
        with mocks_with_real_numpy():
            mod = load_plugin(BLENDER_PATH)
            import rdkit  # mocked; stable within this context

            bt = rdkit.Chem.rdchem.BondType
            kinds = [getattr(bt, k) for k in bond_kinds]
            mol = _simple_mol(kinds)
            mw = make_export_mw(style=style, color_map=color_map)
            return mod.generate_blender_script(mol, mw), mod

    def test_header_contains_plugin_version(self):
        script, mod = self._generate()
        assert f"{mod.PLUGIN_NAME} v{mod.PLUGIN_VERSION}" in script

    def test_script_is_valid_python(self):
        script, _ = self._generate(bond_kinds=("SINGLE", "DOUBLE", "TRIPLE"))
        ast.parse(script)  # must not raise

    def test_clear_scene_flag_present(self):
        script, _ = self._generate()
        assert "CLEAR_SCENE = True" in script

    def test_atom_objects_created_for_each_atom(self):
        script, _ = self._generate(bond_kinds=("SINGLE", "SINGLE"))
        assert script.count("obj = create_atom(") == 3
        assert 'name="Atom_0_C"' in script
        assert 'name="Atom_2_C"' in script

    def test_single_bond_fallback_color(self):
        script, _ = self._generate(bond_kinds=("SINGLE",))
        assert script.count("obj = create_bond(") == 1
        assert "color=(0.6, 0.6, 0.6)," in script

    def test_double_bond_creates_four_segments(self):
        # Empty color map → each of the 2 parallel cylinders is CPK-split in 2.
        script, _ = self._generate(bond_kinds=("DOUBLE",))
        assert script.count("obj = create_bond(") == 4

    def test_triple_bond_creates_six_segments(self):
        # Center + 2 side cylinders, each split at the midpoint.
        script, _ = self._generate(bond_kinds=("TRIPLE",))
        assert script.count("obj = create_bond(") == 6

    def test_atom_color_map_override(self):
        script, _ = self._generate(color_map={"atom_0": (255, 0, 0)})
        assert "color=(1.0, 0.0, 0.0)," in script

    def test_cpk_style_has_no_bonds(self):
        script, _ = self._generate(style="cpk", bond_kinds=("SINGLE",))
        assert script.count("obj = create_bond(") == 0
        assert script.count("obj = create_atom(") == 2

    def test_wireframe_style_has_no_atoms(self):
        script, _ = self._generate(style="wireframe", bond_kinds=("SINGLE",))
        assert script.count("obj = create_atom(") == 0
        assert script.count("obj = create_bond(") == 1

    def test_fallback_camera_centered_on_molecule(self):
        script, _ = self._generate(bond_kinds=("SINGLE",))
        # Two atoms at x=0 and x=1.5 → center x = 0.75
        assert "look_at=(0.750000, 0.000000, 0.000000)," in script


# ---------------------------------------------------------------------------
# POV-Ray Export — generate_povray_scene content
# ---------------------------------------------------------------------------


class _StubQColor:
    """Minimal QColor stand-in parsing '#rrggbb'."""

    def __init__(self, hexstr):
        h = str(hexstr).lstrip("#")
        self.r = int(h[0:2], 16) / 255.0
        self.g = int(h[2:4], 16) / 255.0
        self.b = int(h[4:6], 16) / 255.0

    def redF(self):
        return self.r

    def greenF(self):
        return self.g

    def blueF(self):
        return self.b


@contextlib.contextmanager
def povray_env():
    """Mock context with real numpy and a usable QColor stub."""
    with mocks_with_real_numpy():
        mod = load_plugin(POVRAY_PATH)
        import PyQt6.QtGui as qtgui  # mocked module object

        qtgui.QColor = _StubQColor
        yield mod


class TestPovraySceneGeneration:
    def _generate(self, bond_kinds=("SINGLE",), style="ball_and_stick",
                  settings=None, plotter=None):
        with povray_env() as mod:
            import rdkit

            bt = rdkit.Chem.rdchem.BondType
            mol = _simple_mol([getattr(bt, k) for k in bond_kinds])
            mw = make_export_mw(style=style, settings=settings, plotter=plotter)
            return mod.generate_povray_scene(mol, mw), mod

    def test_returns_scene_and_scaled_resolution(self):
        (scene, w, h), _ = self._generate()
        # Fallback 1920x1080 → 1.5x smart scaling
        assert (w, h) == (2880, 1620)
        assert isinstance(scene, str) and scene

    def test_small_window_gets_2x_scaling(self):
        plotter = SimpleNamespace(window_size=(800, 600))
        (_, w, h), _ = self._generate(plotter=plotter)
        assert (w, h) == (1600, 1200)

    def test_sphere_per_atom(self):
        (scene, _, _), _ = self._generate(bond_kinds=("SINGLE", "SINGLE"))
        assert scene.count("\nsphere {") == 3

    def test_single_bond_fallback_gray_cylinder(self):
        (scene, _, _), _ = self._generate()
        assert scene.count("cylinder {") == 1
        assert "color rgb <0.600, 0.600, 0.600>" in scene

    def test_double_bond_creates_four_cylinders(self):
        (scene, _, _), _ = self._generate(bond_kinds=("DOUBLE",))
        assert scene.count("cylinder {") == 4

    def test_background_color_from_settings(self):
        (scene, _, _), _ = self._generate(settings={"background_color": "#ff0000"})
        assert "background { color rgb <1.000, 0.000, 0.000> }" in scene

    def test_end_of_scene_summary(self):
        (scene, _, _), _ = self._generate(bond_kinds=("SINGLE", "SINGLE"))
        assert "// End of scene - 3 atoms, 2 bonds" in scene

    def test_camera_from_plotter(self):
        cam = SimpleNamespace(
            GetPosition=lambda: (1.0, 2.0, 3.0),
            GetFocalPoint=lambda: (0.0, 0.0, 0.0),
            GetViewUp=lambda: (0.0, 1.0, 0.0),
            GetParallelProjection=lambda: False,
            GetViewAngle=lambda: 30.0,
        )
        plotter = SimpleNamespace(window_size=(1920, 1080), camera=cam)
        (scene, _, _), _ = self._generate(plotter=plotter)
        assert "location <1.0000, 2.0000, 3.0000>" in scene
        assert "angle 30.0" in scene

    def test_wireframe_has_no_spheres(self):
        (scene, _, _), _ = self._generate(style="wireframe")
        assert scene.count("\nsphere {") == 0
        assert scene.count("cylinder {") == 1


class TestPovrayTripleBondRegression:
    """
    Regression: the TRIPLE-bond branch used ``off_dir`` without assigning it.
    A molecule whose first multiple bond is a triple bond raised
    UnboundLocalError (or reused a stale offset from an earlier double bond).
    """

    def test_triple_bond_first_does_not_raise(self):
        (scene, _, _), _ = TestPovraySceneGeneration()._generate(
            bond_kinds=("TRIPLE",)
        )
        # Center split (2) + two offset sides split (2 each) = 6 cylinders
        assert scene.count("cylinder {") == 6

    def test_triple_bond_sides_are_offset(self):
        (scene, _, _), _ = TestPovraySceneGeneration()._generate(
            bond_kinds=("TRIPLE",)
        )
        # ball_and_stick: cyl_radius=0.1, triple_offset_factor=2.0 → ±0.2
        # Bond along x → offset along ±y (cross with z-arbitrary vector).
        assert "<0.0000, -0.2000, 0.0000>," in scene or "<0.0000, 0.2000, 0.0000>," in scene
