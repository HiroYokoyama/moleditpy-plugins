"""
Headless GUI tests for the Settings Saver plugin.

Covers: SettingsSaverDialog.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

SETTINGS_PATH = PLUGINS_DIR / "Settings_Saver" / "settings_saver.py"

with mock_chemistry_imports():
    _settings = load_plugin_for_gui(SETTINGS_PATH)


# ===========================================================================
# SettingsSaverDialog  (Settings Saver)
# ===========================================================================


def _settings_context() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


class TestSettingsSaverDialog:
    """SettingsSaverDialog with MagicMock context and parent=None."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _settings_context()
        d = _settings.SettingsSaverDialog(context=ctx, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Settings Saver Manager"

    def test_preset_list_initially_empty(self, dlg):
        assert dlg.preset_list.count() == 0

    def test_load_button_exists(self, dlg):
        assert dlg.btn_load.text() == "Load Preset"

    def test_save_button_exists(self, dlg):
        assert dlg.btn_save.text() == "Save New..."

    def test_delete_button_exists(self, dlg):
        assert dlg.btn_delete.text() == "Delete"

    def test_embed_checkbox_unchecked_by_default(self, dlg):
        assert not dlg.chk_embed.isChecked()

    def test_global_default_button_disabled_when_embed_off(self, dlg):
        assert not dlg.btn_set_global.isEnabled()

    def test_export_button_has_menu(self, dlg):
        assert dlg.btn_export.menu() is not None

    def test_close_button_exists(self, dlg):
        assert dlg.btn_close.text() == "Close"
