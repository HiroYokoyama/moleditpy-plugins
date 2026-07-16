"""
Headless GUI tests for the GAMESS Input Generator plugin.

Covers: GamessSetupDialog.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

GAMESS_PATH = PLUGINS_DIR / "Gamess_Input_Generator" / "gamess_input_generator.py"

with mock_chemistry_imports():
    _gamess = load_plugin_for_gui(GAMESS_PATH)


# ===========================================================================
# GamessSetupDialog  (visible plugin: "GAMESS Input Generator")
# ===========================================================================


class TestGamessSetupDialog:
    """GamessSetupDialog with mol=None."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _gamess.GamessSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "GAMESS Input Generator"

    def test_run_type_default(self, dlg):
        assert dlg.run_type.currentText() == "OPTIMIZE"

    def test_scf_type_default(self, dlg):
        assert dlg.scf_type.currentText() == "RHF"

    def test_nosym_checked_by_default(self, dlg):
        assert dlg.chk_nosym.isChecked()

    def test_basis_ngauss_default(self, dlg):
        assert dlg.basis_ngauss.value() == 6

    def test_mem_spin_default(self, dlg):
        assert dlg.mem_spin.value() == 100
