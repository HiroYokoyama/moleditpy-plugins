"""
Tests for the Structural Updater plugin.

Visible: true in REGISTRY/plugins.json.
All tests run headlessly — PyQt6, rdkit, etc. are mocked.
"""
from __future__ import annotations

import json
from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
STRUCTURAL_PATH = PLUGINS_DIR / "Structural_Updater" / "structural_updater.py"

# Load module once at module scope
with mock_optional_imports():
    SU = load_plugin(STRUCTURAL_PATH)


# ---------------------------------------------------------------------------
# Structural Updater — initialize / finalize
# ---------------------------------------------------------------------------

class TestStructuralUpdaterInitialize:
    def test_initialize_does_not_raise(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)

    def test_initialize_sets_plugin_instance(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)
        assert SU._PLUGIN_INSTANCE is not None

    def test_initialize_adds_menu_action(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)
        ctx.add_menu_action.assert_called()
        path = ctx.add_menu_action.call_args[0][0]
        assert "Structural Updater" in path

    def test_finalize_does_not_raise(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)
            SU.finalize()


# ---------------------------------------------------------------------------
# Structural Updater — settings load / save
# ---------------------------------------------------------------------------

class TestStructuralUpdaterSettings:
    def test_default_enabled_is_true(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)
        assert SU._PLUGIN_INSTANCE.enabled is True

    def test_save_creates_file(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        plugin.save_settings()
        assert (tmp_path / "su.json").exists()

    def test_save_writes_enabled_true(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        plugin.enabled = True
        plugin.save_settings()
        data = json.loads((tmp_path / "su.json").read_text())
        assert data == {"enabled": True}

    def test_save_writes_enabled_false(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        plugin.enabled = False
        plugin.save_settings()
        data = json.loads((tmp_path / "su.json").read_text())
        assert data == {"enabled": False}

    def test_load_reads_enabled_false(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        (tmp_path / "su.json").write_text(json.dumps({"enabled": False}))
        plugin.load_settings()
        assert plugin.enabled is False

    def test_load_reads_enabled_true(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        (tmp_path / "su.json").write_text(json.dumps({"enabled": True}))
        plugin.load_settings()
        assert plugin.enabled is True

    def test_load_missing_file_creates_default(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        settings_file = tmp_path / "su.json"
        plugin.settings_file = str(settings_file)
        plugin.load_settings()
        assert settings_file.exists()  # save_settings() called automatically

    def test_round_trip(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        plugin.enabled = False
        plugin.save_settings()
        plugin.enabled = True  # change in memory
        plugin.load_settings()  # reload from file
        assert plugin.enabled is False


# ---------------------------------------------------------------------------
# Structural Updater — check_state
# ---------------------------------------------------------------------------

class TestStructuralUpdaterCheckState:
    def test_check_state_disabled_exits_early(self):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.enabled = False
        plugin.check_state()  # must not raise

    def test_check_state_enabled_does_not_raise(self):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.enabled = True
        # MagicMock() > 0 raises TypeError — configure GetNumAtoms to return int
        plugin.context.current_molecule.GetNumAtoms.return_value = 0
        plugin.check_state()  # must not raise
