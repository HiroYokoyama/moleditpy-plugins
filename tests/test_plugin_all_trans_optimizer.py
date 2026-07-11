"""
Tests for the All-Trans Optimizer plugin.

Covers:
  1. initialize() must register at least one menu/export/plugin action
  2. No-molecule guard paths (run_plugin with mol=None)
  3. Dialog accept/reject round-trips
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
ALL_TRANS_PATH = PLUGINS_DIR / "All-Trans_Optimizer" / "all-trans_optimizer.py"


class TestAllTransOptimizer:
    def test_initialize_registers_3d_edit_menu_action(self):
        """initialize() registers exactly one menu action in the 3D Edit menu."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()
            path = ctx.add_menu_action.call_args[0][0]
            assert path == "3D Edit/All-Trans Optimizer"

    def test_initialize_menu_action_invokes_run_plugin(self):
        """The registered callback runs run_plugin(context) without raising."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            ctx.current_mol = None
            mod.initialize(ctx)
            callback = ctx.add_menu_action.call_args[0][1]
            with patch.object(mod.QMessageBox, "warning"):
                callback()

    def test_run_plugin_no_mol_does_not_raise(self):
        """run_plugin with mol=None should show a warning, not raise."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            mock_warn.assert_called_once()

    def test_run_plugin_no_mol_shows_correct_message(self):
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            args = mock_warn.call_args[0]
            assert "molecule" in args[2].lower() or "No molecule" in args[2]

    def test_module_has_no_legacy_run_entry_point(self):
        """Pure V4 plugin: run() must not exist (would duplicate the Plugin menu entry)."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            assert not hasattr(mod, "run")


class TestAllTransRunPlugin:
    def _load(self):
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
        return mod

    def _context(self, mol):
        ctx = MagicMock()
        ctx.get_main_window.return_value = MagicMock()
        ctx.current_mol = mol
        return ctx

    def test_no_conformer_warns_and_stops(self):
        mod = self._load()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 0
        ctx = self._context(mol)
        mod.QMessageBox.reset_mock()
        mod.rdMolTransforms.reset_mock()
        mod.run_plugin(ctx)
        mod.QMessageBox.warning.assert_called_once()
        mod.rdMolTransforms.SetDihedralDeg.assert_not_called()

    def test_matches_set_to_180_and_view_refreshed(self):
        mod = self._load()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        conf = MagicMock()
        mol.GetConformer.return_value = conf
        mol.GetSubstructMatches.return_value = [(0, 1, 2, 3), (4, 5, 6, 7)]
        ctx = self._context(mol)
        mod.rdMolTransforms.reset_mock()
        mod.QMessageBox.reset_mock()

        mod.run_plugin(ctx)

        calls = mod.rdMolTransforms.SetDihedralDeg.call_args_list
        assert len(calls) == 2
        assert calls[0].args == (conf, 0, 1, 2, 3, 180.0)
        assert calls[1].args == (conf, 4, 5, 6, 7, 180.0)
        ctx.refresh_3d_view.assert_called_once()
        ctx.push_undo_checkpoint.assert_called_once()
        assert "2" in ctx.show_status_message.call_args[0][0]

    def test_no_matches_informs_without_refresh(self):
        mod = self._load()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        mol.GetSubstructMatches.return_value = []
        ctx = self._context(mol)
        mod.rdMolTransforms.reset_mock()
        mod.QMessageBox.reset_mock()

        mod.run_plugin(ctx)

        mod.QMessageBox.information.assert_called_once()
        mod.rdMolTransforms.SetDihedralDeg.assert_not_called()
        ctx.refresh_3d_view.assert_not_called()
        ctx.push_undo_checkpoint.assert_not_called()

    def test_duplicate_central_bond_applied_only_once(self):
        """A branched backbone atom yields several matches sharing one central
        bond; only one dihedral must be set for that bond so later calls do not
        undo earlier ones."""
        mod = self._load()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        conf = MagicMock()
        mol.GetConformer.return_value = conf
        # Three matches, first two share the central bond (1, 2).
        mol.GetSubstructMatches.return_value = [
            (0, 1, 2, 3),
            (9, 1, 2, 8),  # same central bond (1, 2) - must be dropped
            (4, 5, 6, 7),
        ]
        ctx = self._context(mol)
        mod.rdMolTransforms.reset_mock()
        mod.QMessageBox.reset_mock()

        mod.run_plugin(ctx)

        calls = mod.rdMolTransforms.SetDihedralDeg.call_args_list
        assert len(calls) == 2
        assert calls[0].args == (conf, 0, 1, 2, 3, 180.0)
        assert calls[1].args == (conf, 4, 5, 6, 7, 180.0)
        assert "2" in ctx.show_status_message.call_args[0][0]

    def test_reversed_central_bond_treated_as_duplicate(self):
        """Central bond (2, 1) and (1, 2) are the same rotatable bond."""
        mod = self._load()
        matches = [(0, 1, 2, 3), (9, 2, 1, 8)]
        selected = mod._select_torsions(matches)
        assert selected == [(0, 1, 2, 3)]


class TestAllTransSmarts:
    def _load(self):
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
        return mod

    def test_smarts_includes_heteroatoms(self):
        """Backbone SMARTS must accept O/N/S/P/Si, not just carbon."""
        mod = self._load()
        smarts = mod._ALL_TRANS_SMARTS
        for atomic_num in ("#7", "#8", "#15", "#16", "#14"):
            assert atomic_num in smarts, f"{atomic_num} missing from backbone SMARTS"

    def test_smarts_central_bond_is_acyclic_single(self):
        mod = self._load()
        assert "-;!@" in mod._ALL_TRANS_SMARTS
