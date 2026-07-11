"""
Tests for the Paste XYZ plugin (parse_xyz_lines, initialize smoke).

All tests run headlessly — PyQt6 / rdkit / etc. are replaced with MagicMock.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

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


class TestPasteXYZParserMore:
    def test_scientific_notation_coords_parsed(self):
        atoms = _parse("C 1e0 -2.5E-1 0.0")
        assert atoms == [("C", 1.0, -0.25, 0.0)]

    def test_tab_separated_columns_parsed(self):
        atoms = _parse("C\t0.0\t1.0\t2.0")
        assert atoms == [("C", 0.0, 1.0, 2.0)]

    def test_extra_trailing_columns_ignored(self):
        # e.g. XYZ variants with charges/gradients appended
        atoms = _parse("O 0.0 0.0 0.0 -0.834 extra")
        assert atoms == [("O", 0.0, 0.0, 0.0)]

    def test_atom_count_header_not_misread_as_atom(self):
        # "5" alone and a comment full of words must both be skipped
        atoms = _parse("5\nenergy = -76.4 hartree scf done\nH 0 0 0")
        assert [a[0] for a in atoms] == ["H"]


PASTE_XYZ_PATH = PLUGINS_DIR / "Paste_XYZ" / "paste_xyz.py"


class TestApplyBondsFallback:
    """_apply_bonds must fall back to distance-based estimation when
    rdDetermineBonds fails, and stamp _xyz_charge either way."""

    def _get_fn(self, chem):
        return extract_function(
            PASTE_XYZ_PATH, None, "_apply_bonds",
            extra_globals={"Chem": chem, "logging": __import__("logging")},
        )

    def test_determine_bonds_failure_falls_back_to_estimator(self):
        chem = MagicMock()
        chem.RWMol.side_effect = RuntimeError("no rdkit")  # force fallback
        fn = self._get_fn(chem)

        result_mol = MagicMock()
        rwmol = MagicMock()
        rwmol.GetMol.return_value = result_mol
        mw = MagicMock()

        out = fn(rwmol, -1, mw)

        mw.io_manager.estimate_bonds_from_distances.assert_called_once_with(rwmol)
        assert out is result_mol
        result_mol.SetIntProp.assert_called_once_with("_xyz_charge", -1)

    def test_estimator_exception_is_swallowed(self):
        chem = MagicMock()
        chem.RWMol.side_effect = RuntimeError("no rdkit")
        fn = self._get_fn(chem)

        rwmol = MagicMock()
        mw = MagicMock()
        mw.io_manager.estimate_bonds_from_distances.side_effect = ValueError("bad geom")

        out = fn(rwmol, 0, mw)  # must not raise
        assert out is rwmol.GetMol.return_value

    def test_mw_without_io_manager_still_returns_mol(self):
        chem = MagicMock()
        chem.RWMol.side_effect = RuntimeError("no rdkit")
        fn = self._get_fn(chem)

        rwmol = MagicMock()
        mw = SimpleNamespace()  # no io_manager at all
        out = fn(rwmol, 2, mw)
        out.SetIntProp.assert_called_once_with("_xyz_charge", 2)


class TestResolveMolWithCharge:
    """Charge-prompt loop of _resolve_mol_with_charge."""

    def _get_fn(self, apply_bonds, prompt_charge):
        return extract_function(
            PASTE_XYZ_PATH, None, "_resolve_mol_with_charge",
            extra_globals={
                "_apply_bonds": apply_bonds,
                "_prompt_charge": prompt_charge,
                "logging": __import__("logging"),
            },
        )

    def test_default_path_uses_charge_zero_without_prompting(self):
        sentinel = object()
        apply_bonds = MagicMock(return_value=sentinel)
        prompt = MagicMock()
        fn = self._get_fn(apply_bonds, prompt)

        ctx = MagicMock()
        ctx.get_setting.return_value = False  # always_ask off
        out = fn(MagicMock(), ctx, MagicMock())

        assert out is sentinel
        apply_bonds.assert_called_once()
        assert apply_bonds.call_args[0][1] == 0
        prompt.assert_not_called()

    def test_cancelling_prompt_returns_none(self):
        apply_bonds = MagicMock(side_effect=RuntimeError("fail"))
        prompt = MagicMock(return_value=(0, False, False))  # user cancels
        fn = self._get_fn(apply_bonds, prompt)

        ctx = MagicMock()
        ctx.get_setting.return_value = True  # always ask
        assert fn(MagicMock(), ctx, MagicMock()) is None

    def test_skip_marks_molecule_and_uses_estimator(self):
        prompt = MagicMock(return_value=(0, True, True))  # skip chemistry
        fn = self._get_fn(MagicMock(), prompt)

        ctx = MagicMock()
        ctx.get_setting.return_value = True
        rwmol = MagicMock()
        mw = MagicMock()
        out = fn(rwmol, ctx, mw)

        mw.io_manager.estimate_bonds_from_distances.assert_called_once_with(rwmol)
        out.SetIntProp.assert_called_once_with("_xyz_skip_checks", 1)

    def test_retry_after_failed_charge_then_success(self):
        sentinel = object()
        apply_bonds = MagicMock(side_effect=[RuntimeError("first"), RuntimeError("bad"), sentinel])
        prompt = MagicMock(side_effect=[(1, True, False), (0, True, False)])
        fn = self._get_fn(apply_bonds, prompt)

        ctx = MagicMock()
        ctx.get_setting.return_value = False  # default path fails -> loop
        out = fn(MagicMock(), ctx, MagicMock())

        assert out is sentinel
        assert prompt.call_count == 2
        ctx.show_status_message.assert_called_once()


class _FakeDialogCode:
    Accepted = 1


class _FakeQDialog:
    DialogCode = _FakeDialogCode


class TestRunPluginCancelKeepsCanvas:
    """Regression: cancelling the charge prompt must NOT clear the canvas
    (before 2026.07.11 the canvas was wiped before charge resolution)."""

    def _run(self, resolved_mol):
        dialog = MagicMock()
        dialog.exec.return_value = _FakeDialogCode.Accepted
        dialog.get_data.return_value = "C 0 0 0"

        ctx = MagicMock()
        ctx.get_setting.return_value = False  # skip_chemistry_checks off

        fn = extract_function(
            PASTE_XYZ_PATH, None, "run_plugin",
            extra_globals={
                "Chem": MagicMock(),
                "QDialog": _FakeQDialog,
                "QMessageBox": MagicMock(),
                "PasteXYZDialog": MagicMock(return_value=dialog),
                "parse_xyz_lines": lambda text: [("C", 0.0, 0.0, 0.0)],
                "_build_rwmol": MagicMock(),
                "_resolve_mol_with_charge": MagicMock(return_value=resolved_mol),
                "logging": __import__("logging"),
            },
        )
        fn(ctx)
        return ctx

    def test_cancel_leaves_canvas_untouched(self):
        ctx = self._run(resolved_mol=None)
        ctx.clear_canvas.assert_not_called()
        ctx.push_undo_checkpoint.assert_not_called()

    def test_success_clears_canvas_and_loads(self):
        mol = MagicMock()
        ctx = self._run(resolved_mol=mol)
        ctx.clear_canvas.assert_called_once_with(push_to_undo=False)
        assert ctx.current_molecule is mol
        ctx.push_undo_checkpoint.assert_called_once()


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
