"""
Headless GUI tests for the MOPAC Input Generator plugin.

Covers: MopacSetupDialog.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MOPAC_PATH = PLUGINS_DIR / "Mopac_Input_Generator" / "mopac_input_generator.py"

with mock_chemistry_imports():
    _mopac = load_plugin_for_gui(MOPAC_PATH)


# ===========================================================================
# MopacSetupDialog  (visible plugin: "MOPAC Input Generator")
# ===========================================================================


class TestMopacSetupDialog:
    """MopacSetupDialog with mol=None — exercises the no-molecule UI path."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _mopac.MopacSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "MOPAC Input Generator"

    def test_default_keywords_field(self, dlg):
        assert dlg.keywords_edit.text() == "PM7 PRECISE"

    def test_default_title_field(self, dlg):
        assert dlg.title_edit.text() == "MOPAC Calculation"

    def test_charge_spinbox_range(self, dlg):
        assert dlg.charge_spin.minimum() == -10
        assert dlg.charge_spin.maximum() == 10

    def test_mult_spinbox_minimum_is_one(self, dlg):
        assert dlg.mult_spin.minimum() == 1

    def test_template_combo_populated(self, dlg):
        assert dlg.template_combo.count() > 0

    def test_nosym_checkbox_initially_unchecked(self, dlg):
        assert not dlg.chk_nosym.isChecked()

    def test_preview_area_is_editable(self, dlg):
        assert not dlg.preview_text.isReadOnly()

    def test_preset_combo_initially_empty(self, dlg):
        assert dlg.preset_combo.count() == 0
