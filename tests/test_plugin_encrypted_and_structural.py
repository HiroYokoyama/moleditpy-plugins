"""
Tests for Encrypted Project and Structural Updater plugins.

Both are visible: true in REGISTRY/plugins.json.
All tests run headlessly — PyQt6, rdkit, cryptography, etc. are mocked.
"""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
ENCRYPTED_PATH = PLUGINS_DIR / "Encrypted_Project" / "encrypted_project.py"
STRUCTURAL_PATH = PLUGINS_DIR / "Structural_Updater" / "structural_updater.py"

# Load modules once at module scope
with mock_optional_imports():
    ENC = load_plugin(ENCRYPTED_PATH)

with mock_optional_imports():
    SU = load_plugin(STRUCTURAL_PATH)


# ---------------------------------------------------------------------------
# Encrypted Project — on_drop
# ---------------------------------------------------------------------------

class TestEncryptedProjectOnDrop:
    def _make_plugin(self):
        with mock_optional_imports():
            ctx = make_context()
            return ENC.PmeencPlugin(ctx)

    def test_pmeenc_extension_returns_true(self):
        assert self._make_plugin().on_drop("file.pmeenc") is True

    def test_pmeenc_extension_case_insensitive(self):
        assert self._make_plugin().on_drop("FILE.PMEENC") is True

    def test_pmeenc_full_path_returns_true(self):
        assert self._make_plugin().on_drop("/home/user/mol.pmeenc") is True

    def test_mol_extension_returns_false(self):
        assert self._make_plugin().on_drop("file.mol") is False

    def test_xyz_extension_returns_false(self):
        assert self._make_plugin().on_drop("file.xyz") is False

    def test_empty_string_returns_false(self):
        assert self._make_plugin().on_drop("") is False

    def test_pdf_extension_returns_false(self):
        assert self._make_plugin().on_drop("report.pdf") is False


# ---------------------------------------------------------------------------
# Encrypted Project — on_document_reset
# ---------------------------------------------------------------------------

class TestEncryptedProjectDocumentReset:
    def test_clears_key_and_salt(self):
        with mock_optional_imports():
            ctx = make_context()
            plugin = ENC.PmeencPlugin(ctx)
        plugin.current_key = b"somekey"
        plugin.current_salt = b"somesalt"
        plugin.on_document_reset()
        assert plugin.current_key is None
        assert plugin.current_salt is None

    def test_reset_from_none_does_not_raise(self):
        with mock_optional_imports():
            ctx = make_context()
            plugin = ENC.PmeencPlugin(ctx)
        # current_key / current_salt already None
        plugin.on_document_reset()  # must not raise


# ---------------------------------------------------------------------------
# Encrypted Project — _patch_method / _unpatch_method
# ---------------------------------------------------------------------------

class TestEncryptedProjectPatchMethod:
    def _make_original(self):
        """Return a new plain function each call.

        Using a local def avoids two pitfalls:
        - MagicMock: getattr(mock, "_is_pmeenc_patch", False) returns a truthy
          MagicMock, so _patch_method's early-exit guard fires and does nothing.
        - Bound methods: self.some_method accessed twice creates two distinct
          wrapper objects, breaking `is` identity assertions.
        """
        def original_fn(): pass
        return original_fn

    def test_patch_replaces_method(self):
        original = self._make_original()
        with mock_optional_imports():
            ctx = make_context()
            plugin = ENC.PmeencPlugin(ctx)
        plugin.mw.save_project = original
        plugin._patch_method("save_project", lambda mw, *a, **kw: None, [])
        assert plugin.mw.save_project is not original

    def test_patched_method_has_flag(self):
        original = self._make_original()
        with mock_optional_imports():
            ctx = make_context()
            plugin = ENC.PmeencPlugin(ctx)
        plugin.mw.save_project = original
        plugin._patch_method("save_project", lambda mw, *a, **kw: None, [])
        assert getattr(plugin.mw.save_project, "_is_pmeenc_patch", False) is True

    def test_patch_stores_original_on_plugin(self):
        original = self._make_original()
        with mock_optional_imports():
            ctx = make_context()
            plugin = ENC.PmeencPlugin(ctx)
        plugin.mw.save_project = original
        plugin._patch_method("save_project", lambda mw, *a, **kw: None, [])
        assert plugin._original_save_project is original

    def test_unpatch_restores_original(self):
        original = self._make_original()
        with mock_optional_imports():
            ctx = make_context()
            plugin = ENC.PmeencPlugin(ctx)
        plugin.mw.save_project = original
        plugin._patch_method("save_project", lambda mw, *a, **kw: None, ["Save Project"])
        plugin._unpatch_method("save_project", ["Save Project"])
        assert plugin.mw.save_project is original

    def test_no_double_patch(self):
        """Second _patch_method call on already-patched method is a no-op."""
        original = self._make_original()
        with mock_optional_imports():
            ctx = make_context()
            plugin = ENC.PmeencPlugin(ctx)
        plugin.mw.save_project = original
        plugin._patch_method("save_project", lambda mw, *a, **kw: None, [])
        first_patch = plugin.mw.save_project
        plugin._patch_method("save_project", lambda mw, *a, **kw: None, [])
        assert plugin.mw.save_project is first_patch

    def test_patch_missing_method_does_not_raise(self):
        """When mw doesn't have the method, _patch_method returns silently."""
        with mock_optional_imports():
            ctx = make_context()
            plugin = ENC.PmeencPlugin(ctx)
        # Use a plain object so hasattr returns False for unknown attrs
        plugin.mw = object()
        plugin._patch_method("no_such_method", lambda: None, [])  # must not raise


# ---------------------------------------------------------------------------
# Encrypted Project — initialize() registrations
# ---------------------------------------------------------------------------

class TestEncryptedProjectInitialize:
    def test_initialize_registers_file_opener_for_pmeenc(self):
        with mock_optional_imports():
            ctx = make_context()
            ENC.initialize(ctx)
        assert ctx.register_file_opener.called
        ext = ctx.register_file_opener.call_args[0][0]
        assert ext == ".pmeenc"

    def test_initialize_registers_export_action(self):
        with mock_optional_imports():
            ctx = make_context()
            ENC.initialize(ctx)
        assert ctx.add_export_action.called

    def test_initialize_registers_document_reset_handler(self):
        with mock_optional_imports():
            ctx = make_context()
            ENC.initialize(ctx)
        assert ctx.register_document_reset_handler.called

    def test_initialize_registers_drop_handler(self):
        with mock_optional_imports():
            ctx = make_context()
            ENC.initialize(ctx)
        assert ctx.register_drop_handler.called


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
