"""
Tests for the Paste XYZ plugin (parse_xyz_lines, initialize smoke).

All tests run headlessly — PyQt6 / rdkit / etc. are replaced with MagicMock.
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

with mock_optional_imports():
    PASTE_XYZ_MOD = load_plugin(
        PLUGINS_DIR / "Paste_XYZ" / "paste_xyz.py"
    )

# parse_xyz_lines is a module-level function in the plugin so tests import it directly.
_parse = PASTE_XYZ_MOD.parse_xyz_lines


class TestPasteXYZParser:
    def test_valid_three_atom_block(self):
        xyz = "C  0.0  0.0  0.0\nH  1.0  0.0  0.0\nH -1.0  0.0  0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 3
        assert atoms[0] == ("C", 0.0, 0.0, 0.0)
        assert atoms[1] == ("H", 1.0, 0.0, 0.0)

    def test_header_lines_skipped(self):
        xyz = "3\nComment line\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 3

    def test_short_lines_skipped(self):
        xyz = "C 0.0 0.0\nH 1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 1
        assert atoms[0][0] == "H"

    def test_non_float_coords_skipped(self):
        xyz = "C abc 0.0 0.0\nH 1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 1
        assert atoms[0][0] == "H"

    def test_non_alpha_symbol_skipped(self):
        xyz = "1 0.0 0.0 0.0\nH 1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 1
        assert atoms[0][0] == "H"

    def test_empty_input_returns_empty(self):
        assert _parse("") == []

    def test_blank_lines_skipped(self):
        xyz = "\n\nC 0.0 0.0 0.0\n\n"
        atoms = _parse(xyz)
        assert len(atoms) == 1

    def test_fractional_coords_parsed(self):
        xyz = "O -0.12345 6.789 -3.141"
        atoms = _parse(xyz)
        assert len(atoms) == 1
        assert atoms[0][0] == "O"
        assert abs(atoms[0][1] - (-0.12345)) < 1e-9

    def test_mixed_case_symbols_accepted(self):
        xyz = "Ca 0.0 0.0 0.0\nFe 1.0 0.0 0.0"
        atoms = _parse(xyz)
        assert len(atoms) == 2
        assert atoms[0][0] == "Ca"
        assert atoms[1][0] == "Fe"


class TestPasteXYZSmoke:
    def test_initialize_registers_menu_action(self):
        """initialize(ctx) must register exactly one menu action under File/."""
        with mock_optional_imports():
            ctx = make_context()
            PASTE_XYZ_MOD.initialize(ctx)
        ctx.add_menu_action.assert_called_once()
        path = ctx.add_menu_action.call_args[0][0]
        assert path.startswith("File/")

    def test_initialize_action_callable(self):
        """The callback registered by initialize must be callable without raising."""
        with mock_optional_imports():
            ctx = make_context()
            PASTE_XYZ_MOD.initialize(ctx)
        callback = ctx.add_menu_action.call_args[0][1]
        assert callable(callback)
