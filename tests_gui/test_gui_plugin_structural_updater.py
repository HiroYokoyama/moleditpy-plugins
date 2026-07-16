"""
Headless GUI tests for the Structural Updater plugin.

Covers: SettingsDialog + StructuralUpdaterPlugin.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

UPDATER_PATH = PLUGINS_DIR / "Structural_Updater" / "structural_updater.py"

with mock_chemistry_imports():
    _updater = load_plugin_for_gui(UPDATER_PATH)


# ===========================================================================
# Structural Updater — SettingsDialog + StructuralUpdaterPlugin
# ===========================================================================


class TestUpdaterSettingsDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _updater.SettingsDialog(parent=None, current_enabled=True)
        yield d
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Structural Updater Settings"

    def test_checkbox_reflects_current_enabled(self, dlg):
        assert dlg.chk_enable.isChecked()

    def test_disabled_state_propagates(self, qapp):
        d = _updater.SettingsDialog(parent=None, current_enabled=False)
        assert not d.chk_enable.isChecked()
        d.destroy()

    def test_accept_picks_up_checkbox_state(self, dlg):
        dlg.chk_enable.setChecked(False)
        dlg.accept()
        assert dlg.enabled is False

    def test_reject_keeps_previous_value(self, dlg):
        dlg.chk_enable.setChecked(False)
        dlg.reject()
        assert dlg.enabled is True


class TestStructuralUpdaterPlugin:
    @pytest.fixture
    def plugin(self, qapp):
        from PyQt6.QtWidgets import QWidget

        # QTimer needs a real QObject parent; managers stay MagicMock.
        w = QWidget()
        w.init_manager = MagicMock()
        w.init_manager.convert_button.text.return_value = "Convert 2D to 3D"
        w.compute_manager = MagicMock()
        ctx = MagicMock()
        ctx.get_main_window.return_value = w
        ctx.current_molecule = None
        saved_originals = dict(_updater._ORIGINAL_METHODS)
        _updater._ORIGINAL_METHODS.clear()
        p = _updater.StructuralUpdaterPlugin(ctx)
        yield p
        p.timer.stop()
        _updater._ORIGINAL_METHODS.clear()
        _updater._ORIGINAL_METHODS.update(saved_originals)
        w.destroy()

    def test_enabled_by_default(self, plugin):
        assert plugin.enabled is True

    def test_menu_action_registered(self, plugin):
        args = plugin.context.add_menu_action.call_args.args
        assert args[0] == "Settings/Structural Updater..."

    def test_timer_running_at_one_second(self, plugin):
        assert plugin.timer.isActive()
        assert plugin.timer.interval() == 1000

    def test_patch_stores_originals(self, plugin):
        assert "trigger_conversion" in _updater._ORIGINAL_METHODS
        assert "on_calculation_finished" in _updater._ORIGINAL_METHODS

    def test_check_state_enters_apply_mode(self, plugin):
        atom = MagicMock()
        atom.HasProp.return_value = True
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 3
        mol.GetAtoms.return_value = [atom]
        plugin.context.current_molecule = mol
        plugin.check_state()
        assert plugin.apply_mode_active is True
        plugin.mw.init_manager.convert_button.setText.assert_called_with(
            "Apply 2D Changes to 3D"
        )

    def test_check_state_exits_apply_mode_without_molecule(self, plugin):
        plugin.apply_mode_active = True
        plugin.context.current_molecule = None
        plugin.check_state()
        assert plugin.apply_mode_active is False
        plugin.mw.init_manager.convert_button.setText.assert_called_with(
            "Convert 2D to 3D"
        )

    def test_check_state_skips_during_running_calculation(self, plugin):
        plugin.mw.init_manager.convert_button.text.return_value = "Halt conversion"
        plugin.apply_mode_active = False
        atom = MagicMock()
        atom.HasProp.return_value = True
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 3
        mol.GetAtoms.return_value = [atom]
        plugin.context.current_molecule = mol
        plugin.check_state()
        assert plugin.apply_mode_active is False

    def test_trigger_falls_back_to_original_when_not_apply_mode(self, plugin):
        plugin.apply_mode_active = False
        plugin.new_trigger_conversion()
        _updater._ORIGINAL_METHODS["trigger_conversion"].assert_called_once()

    def test_trigger_respects_temp_mode_override(self, plugin):
        plugin.apply_mode_active = True
        plugin.mw._temp_conv_mode = "force_full"
        plugin.new_trigger_conversion()
        _updater._ORIGINAL_METHODS["trigger_conversion"].assert_called_once()
