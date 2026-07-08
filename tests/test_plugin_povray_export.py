"""
Tests for the POV-Ray Export plugin (plugins/POV-Ray_Export/povray_export.py).

All heavy deps (PyQt6, rdkit, numpy) are stubbed via mock_optional_imports().
"""

from __future__ import annotations

import contextlib
from pathlib import Path
from types import SimpleNamespace

from conftest import FakeMol, load_plugin, make_context, mock_optional_imports, mocks_with_real_numpy

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
POVRAY_PATH = PLUGINS_DIR / "POV-Ray_Export" / "povray_export.py"

with mock_optional_imports():
    _povray = load_plugin(POVRAY_PATH)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def _simple_mol(bond_types):
    """Linear chain C-C(-C...) with the given rdkit-mocked bond types."""
    n = len(bond_types) + 1
    symbols = ["C"] * n
    coords = [(1.5 * i, 0.0, 0.0) for i in range(n)]
    bonds = [(i, i + 1, bond_types[i]) for i in range(len(bond_types))]
    return FakeMol(symbols, coords, bonds)


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


# ---------------------------------------------------------------------------
# POV-Ray Export — generate_povray_scene content
# ---------------------------------------------------------------------------


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


class TestPOVRayExportInitialize:
    def test_initialize_calls_add_export_action(self):
        with mock_optional_imports():
            mod = load_plugin(POVRAY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_export_action_label_contains_povray(self):
        with mock_optional_imports():
            mod = load_plugin(POVRAY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        label = ctx.add_export_action.call_args[0][0]
        assert "POV-Ray" in label

    def test_export_action_callback_is_callable(self):
        with mock_optional_imports():
            mod = load_plugin(POVRAY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        callback = ctx.add_export_action.call_args[0][1]
        assert callable(callback)

    def test_plugin_version_constant_present(self):
        assert hasattr(_povray, "PLUGIN_VERSION")
        assert _povray.PLUGIN_VERSION
