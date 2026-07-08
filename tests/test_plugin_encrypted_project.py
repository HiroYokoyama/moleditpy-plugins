"""
Tests for the Encrypted Project plugin.

Visible: true in REGISTRY/plugins.json.
All tests run headlessly — PyQt6, rdkit, cryptography, etc. are mocked.
"""
from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
ENCRYPTED_PATH = PLUGINS_DIR / "Encrypted_Project" / "encrypted_project.py"

# Load module once at module scope
with mock_optional_imports():
    ENC = load_plugin(ENCRYPTED_PATH)


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
# derive_key, export/import flow, salt header
# ---------------------------------------------------------------------------

import base64
import json
from unittest.mock import MagicMock, patch

import pytest


def _enc_plugin(mod):
    ctx = make_context()
    plugin = mod.PmeencPlugin(ctx)
    return plugin, ctx


class TestEncryptedDeriveKey:
    def test_pbkdf2_parameters(self):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            captured = {}

            class FakeKDF:
                def __init__(self, **kwargs):
                    captured.update(kwargs)

                def derive(self, pw_bytes):
                    captured["password"] = pw_bytes
                    return b"K" * 32

            mod.PBKDF2HMAC = FakeKDF
            plugin, _ = _enc_plugin(mod)
            key = plugin.derive_key("hunter2", b"S" * 16)

            assert captured["length"] == 32
            assert captured["iterations"] == 100000
            assert captured["salt"] == b"S" * 16
            assert captured["password"] == b"hunter2"
            assert key == base64.urlsafe_b64encode(b"K" * 32)


class TestEncryptedExport:
    def _export(self, mod, tmp_path, cached=False, cancel=False, fname="proj.pmeenc"):
        plugin, ctx = _enc_plugin(mod)
        target = tmp_path / fname

        mod.QInputDialog.getText.reset_mock()
        if cancel:
            mod.QInputDialog.getText.return_value = ("", False)
        else:
            mod.QInputDialog.getText.return_value = ("pw", True)

        plugin.derive_key = MagicMock(return_value=b"DERIVED_KEY")
        plugin._apply_patches = MagicMock()
        plugin.mw.state_manager.create_json_data.return_value = {"atoms": [1]}

        fernet_obj = MagicMock()
        fernet_obj.encrypt.return_value = b"ENCRYPTED_BLOB"
        mod.Fernet = MagicMock(return_value=fernet_obj)

        if cached:
            plugin.current_key = b"CACHED_KEY"
            plugin.current_salt = b"C" * 16
            plugin.mw.init_manager.current_file_path = str(target)
        else:
            plugin.mw.init_manager.current_file_path = None

        with patch("os.urandom", return_value=b"R" * 16):
            plugin.export_encrypted(str(target))
        return plugin, ctx, target

    def test_file_layout_salt_then_ciphertext(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, _, target = self._export(mod, tmp_path)
            data = target.read_bytes()
            assert data[:16] == b"R" * 16
            assert data[16:] == b"ENCRYPTED_BLOB"

    def test_key_and_salt_cached_after_save(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, _, target = self._export(mod, tmp_path)
            assert plugin.current_key == b"DERIVED_KEY"
            assert plugin.current_salt == b"R" * 16
            assert plugin.mw.state_manager.has_unsaved_changes is False
            assert plugin.mw.init_manager.current_file_path == str(target)

    def test_patches_applied_after_successful_save(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, _, _ = self._export(mod, tmp_path)
            plugin._apply_patches.assert_called_once()

    def test_cached_key_skips_password_prompt(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, _, target = self._export(mod, tmp_path, cached=True)
            mod.QInputDialog.getText.assert_not_called()
            data = target.read_bytes()
            assert data[:16] == b"C" * 16

    def test_cancel_writes_nothing(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, ctx, target = self._export(mod, tmp_path, cancel=True)
            assert not target.exists()
            ctx.show_status_message.assert_called_with("Export cancelled.")


class TestEncryptedPatchedSave:
    def test_pmeenc_path_routes_to_encrypted_export(self):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, _ = _enc_plugin(mod)
            plugin.export_encrypted = MagicMock()
            mw_instance = MagicMock()
            mw_instance.current_file_path = "C:/data/proj.PMEENC"
            plugin._patched_save_project(mw_instance)
            plugin.export_encrypted.assert_called_once_with(
                "C:/data/proj.PMEENC"
            )

    def test_plain_path_falls_back_to_original(self):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, _ = _enc_plugin(mod)
            plugin.export_encrypted = MagicMock()
            original = MagicMock()
            plugin._original_save_project = original
            mw_instance = MagicMock()
            mw_instance.current_file_path = "C:/data/proj.pme"
            plugin._patched_save_project(mw_instance)
            original.assert_called_once()
            plugin.export_encrypted.assert_not_called()

    def test_original_typeerror_falls_back_to_no_args(self):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, _ = _enc_plugin(mod)
            calls = []

            def original(*args):
                if args:
                    raise TypeError("no args expected")
                calls.append("bare")

            plugin._original_save_project = original
            mw_instance = MagicMock()
            mw_instance.current_file_path = None
            plugin._patched_save_project(mw_instance, True)
            assert calls == ["bare"]


class TestEncryptedImport:
    def test_cancel_raises_value_error(self):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, _ = _enc_plugin(mod)
            mod.QInputDialog.getText.return_value = ("", False)
            with pytest.raises(ValueError, match="cancelled"):
                plugin.on_import("x.pmeenc")

    def test_salt_read_from_first_16_bytes(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(ENCRYPTED_PATH)
            plugin, _ = _enc_plugin(mod)
            target = tmp_path / "in.pmeenc"
            target.write_bytes(b"S" * 16 + b"CIPHER")

            mod.QInputDialog.getText.return_value = ("pw", True)
            plugin.derive_key = MagicMock(return_value=b"KEY")
            fernet_obj = MagicMock()
            fernet_obj.decrypt.return_value = json.dumps({"a": 1}).encode()
            mod.Fernet = MagicMock(return_value=fernet_obj)

            plugin._apply_patches = MagicMock()
            plugin.on_import(str(target))

            plugin.derive_key.assert_called_once_with("pw", b"S" * 16)
            fernet_obj.decrypt.assert_called_once_with(b"CIPHER")
            assert plugin.current_salt == b"S" * 16
            plugin.mw.state_manager.load_from_json_data.assert_called_once_with(
                {"a": 1}
            )
            plugin._apply_patches.assert_called_once()
