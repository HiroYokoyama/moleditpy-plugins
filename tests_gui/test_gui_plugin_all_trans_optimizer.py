"""
Headless GUI tests for the All-Trans Optimizer plugin.

Covers: _select_torsions, run_plugin guards (no dialog class).

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_all_trans_optimizer.py
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ALL_TRANS_PATH = PLUGINS_DIR / "All-Trans_Optimizer" / "all-trans_optimizer.py"

with mock_chemistry_imports():
    _all_trans = load_plugin_for_gui(ALL_TRANS_PATH)


def _ctx(setting="MMFF_RDKIT") -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.get_setting.return_value = setting
    return ctx


@pytest.fixture
def no_msgbox(monkeypatch):
    """Replace the module's QMessageBox so warnings/info don't block offscreen."""
    box = MagicMock()
    monkeypatch.setattr(_all_trans, "QMessageBox", box)
    return {_all_trans.__name__: box}


# ===========================================================================
# All-Trans Optimizer  (no dialog — torsion selection + guards)
# ===========================================================================


class TestAllTransSelectTorsions:
    def test_keeps_one_match_per_central_bond(self):
        matches = [(0, 1, 2, 3), (4, 1, 2, 5), (6, 1, 2, 7)]
        assert _all_trans._select_torsions(matches) == [(0, 1, 2, 3)]

    def test_reversed_central_bond_is_same_bond(self):
        matches = [(0, 1, 2, 3), (5, 2, 1, 6)]
        assert _all_trans._select_torsions(matches) == [(0, 1, 2, 3)]

    def test_distinct_bonds_all_kept_in_order(self):
        matches = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
        assert _all_trans._select_torsions(matches) == matches

    def test_empty_matches(self):
        assert _all_trans._select_torsions([]) == []


class TestAllTransRunPlugin:
    def test_no_molecule_warns(self, no_msgbox):
        ctx = _ctx()
        _all_trans.run_plugin(ctx)
        no_msgbox[_all_trans.__name__].warning.assert_called_once()

    def test_no_conformer_warns(self, no_msgbox):
        ctx = _ctx()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 0
        ctx.current_mol = mol
        _all_trans.run_plugin(ctx)
        msg = no_msgbox[_all_trans.__name__].warning.call_args[0][2]
        assert "3D" in msg

    def test_no_torsions_informs(self, no_msgbox):
        ctx = _ctx()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        mol.GetSubstructMatches.return_value = ()
        ctx.current_mol = mol
        _all_trans.run_plugin(ctx)
        no_msgbox[_all_trans.__name__].information.assert_called_once()
        ctx.push_undo_checkpoint.assert_not_called()

    def test_applies_one_dihedral_per_bond(self, no_msgbox, monkeypatch):
        set_dihedral = MagicMock()
        monkeypatch.setattr(
            _all_trans.rdMolTransforms, "SetDihedralDeg", set_dihedral
        )
        ctx = _ctx()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        mol.GetSubstructMatches.return_value = ((0, 1, 2, 3), (4, 1, 2, 5))
        ctx.current_mol = mol
        _all_trans.run_plugin(ctx)
        assert set_dihedral.call_count == 1
        assert set_dihedral.call_args[0][-1] == 180.0
        ctx.refresh_3d_view.assert_called_once()
        ctx.push_undo_checkpoint.assert_called_once()
        assert "1 torsions" in ctx.show_status_message.call_args[0][0]

    def test_initialize_registers_menu_action(self):
        ctx = _ctx()
        _all_trans.initialize(ctx)
        assert ctx.add_menu_action.call_args[0][0] == "3D Edit/All-Trans Optimizer"

    def test_plugin_name(self):
        assert _all_trans.PLUGIN_NAME == "All-Trans Optimizer"
