"""
Tests for the Blender Export plugin (plugins/Blender_Export/blender_export.py).

All heavy deps (PyQt6, rdkit, numpy) are stubbed via mock_optional_imports().
"""

from __future__ import annotations

import ast
from pathlib import Path
from types import SimpleNamespace

import numpy as real_numpy
import pytest

from conftest import FakeMol, load_plugin, make_context, mock_optional_imports, mocks_with_real_numpy

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
BLENDER_PATH = PLUGINS_DIR / "Blender_Export" / "blender_export.py"

with mock_optional_imports():
    _blender = load_plugin(BLENDER_PATH)


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


class TestBlenderExportInitialize:
    def test_initialize_calls_add_export_action(self):
        with mock_optional_imports():
            mod = load_plugin(BLENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_export_action_label_contains_blender(self):
        with mock_optional_imports():
            mod = load_plugin(BLENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        label = ctx.add_export_action.call_args[0][0]
        assert "Blender" in label

    def test_export_action_callback_is_callable(self):
        with mock_optional_imports():
            mod = load_plugin(BLENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        callback = ctx.add_export_action.call_args[0][1]
        assert callable(callback)

    def test_plugin_version_constant_present(self):
        assert hasattr(_blender, "PLUGIN_VERSION")
        assert _blender.PLUGIN_VERSION
