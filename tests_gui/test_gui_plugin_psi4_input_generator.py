"""
Headless GUI tests for the Psi4 Input Generator plugin.

Covers: Psi4SetupDialog.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PSI4_PATH = PLUGINS_DIR / "Psi4_Input_Generator" / "psi4_input_generator.py"

with mock_chemistry_imports():
    _psi4 = load_plugin_for_gui(PSI4_PATH)


# ===========================================================================
# Psi4SetupDialog  (visible plugin: "Psi4 Input Generator")
# ===========================================================================


class TestPsi4SetupDialog:
    """Psi4SetupDialog with mol=None."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _psi4.Psi4SetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Psi4 Input Generator"

    def test_method_default(self, dlg):
        assert dlg.method_combo.currentText() == "b3lyp"

    def test_basis_default(self, dlg):
        assert dlg.basis_combo.currentText() == "def2-svp"

    def test_reference_default(self, dlg):
        assert dlg.ref_combo.currentText() == "rks"

    def test_mem_spin_default(self, dlg):
        assert dlg.mem_spin.value() == 2

    def test_thread_spin_default(self, dlg):
        assert dlg.thread_spin.value() == 4
