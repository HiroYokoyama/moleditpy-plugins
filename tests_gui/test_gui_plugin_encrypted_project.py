"""
Headless GUI tests for the Encrypted Project plugin.

Covers: PmeencPlugin (patch/unpatch with real QActions).
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest import mock
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ENC_PATH = PLUGINS_DIR / "Encrypted_Project" / "encrypted_project.py"

with mock_chemistry_imports():
    _enc = load_plugin_for_gui(ENC_PATH)


# ===========================================================================
# PmeencPlugin  (Encrypted Project)
# ===========================================================================


def _enc_plugin(mw=None):
    ctx = MagicMock()
    ctx.get_main_window.return_value = mw if mw is not None else MagicMock()
    plugin = _enc.PmeencPlugin(ctx)
    return plugin, ctx


class TestPmeencPlugin:
    def test_cryptography_flag_true_in_test_env(self):
        assert _enc.CRYPTOGRAPHY_AVAILABLE

    def test_initialize_registers_hooks(self, qapp):
        plugin, ctx = _enc_plugin()
        plugin.initialize()
        ctx.register_file_opener.assert_called_once()
        assert ctx.register_file_opener.call_args.args[0] == ".pmeenc"
        ctx.add_export_action.assert_called_once()
        ctx.register_document_reset_handler.assert_called_once()
        ctx.register_drop_handler.assert_called_once()

    def test_drop_accepts_pmeenc(self, qapp):
        plugin, _ = _enc_plugin()
        plugin.on_import = MagicMock()  # keep the deferred QTimer call harmless
        assert plugin.on_drop("project.pmeenc") is True

    def test_drop_rejects_other_extensions(self, qapp):
        plugin, _ = _enc_plugin()
        assert plugin.on_drop("project.pmeprj") is False

    def test_document_reset_clears_keys(self, qapp):
        plugin, ctx = _enc_plugin()
        plugin.current_key = b"key"
        plugin.current_salt = b"salt"
        ctx.current_molecule = None
        plugin.on_document_reset()
        assert plugin.current_key is None
        assert plugin.current_salt is None

    def test_document_reset_clears_colorizer_overrides(self, qapp):
        plugin, ctx = _enc_plugin()
        plugin.mw._plugin_color_overrides = {0: "#ff0000"}
        plugin.on_document_reset()
        assert plugin.mw._plugin_color_overrides == {}
        ctx.draw_molecule_3d.assert_called_once()

    @pytest.fixture
    def patched_widget(self, qapp):
        from PyQt6.QtGui import QAction
        from PyQt6.QtWidgets import QWidget

        w = QWidget()
        action = QAction("Save Project", w)
        w.addAction(action)
        calls = []
        w.save_project = lambda: calls.append("orig")
        plugin, _ = _enc_plugin(mw=w)
        yield plugin, w, action, calls
        w.destroy()

    def test_patch_method_replaces_and_reconnects(self, patched_widget):
        plugin, w, action, calls = patched_widget
        plugin._patch_method(
            "save_project",
            lambda mw, *a, **k: calls.append("patched"),
            ["Save Project"],
        )
        assert getattr(w.save_project, "_is_pmeenc_patch", False)
        action.trigger()
        assert calls == ["patched"]

    def test_patch_method_is_idempotent(self, patched_widget):
        plugin, w, action, calls = patched_widget
        patch_fn = lambda mw, *a, **k: calls.append("patched")
        plugin._patch_method("save_project", patch_fn, ["Save Project"])
        first = w.save_project
        plugin._patch_method("save_project", patch_fn, ["Save Project"])
        assert w.save_project is first

    def test_unpatch_method_restores_original(self, patched_widget):
        plugin, w, action, calls = patched_widget
        original = w.save_project
        plugin._patch_method(
            "save_project",
            lambda mw, *a, **k: calls.append("patched"),
            ["Save Project"],
        )
        plugin._unpatch_method("save_project", ["Save Project"])
        assert w.save_project is original
        action.trigger()
        assert calls == ["orig"]

    def test_patched_clear_all_forwards_to_original(self, qapp):
        plugin, _ = _enc_plugin()
        called = []
        plugin._original_clear_all = lambda *a, **k: called.append(a)
        plugin._patched_clear_all(plugin.mw, True)
        assert called == [(True,)]

    def test_patched_clear_all_typeerror_falls_back_to_bare_call(self, qapp):
        plugin, _ = _enc_plugin()
        calls = []

        def original(*args):
            if args:
                raise TypeError("no args expected")
            calls.append("bare")

        plugin._original_clear_all = original
        plugin._patched_clear_all(plugin.mw, True)
        assert calls == ["bare"]

    def test_apply_patches_calls_patch_method_twice(self, qapp):
        plugin, _ = _enc_plugin()
        plugin._patch_method = MagicMock()
        plugin._apply_patches()
        assert plugin._patch_method.call_count == 2
        names = [c.args[0] for c in plugin._patch_method.call_args_list]
        assert names == ["save_project", "clear_all"]

    def test_initialize_disabled_without_cryptography(self, qapp, monkeypatch, capsys):
        monkeypatch.setattr(_enc, "CRYPTOGRAPHY_AVAILABLE", False)
        plugin, ctx = _enc_plugin()
        plugin.initialize()
        ctx.register_file_opener.assert_not_called()
        out = capsys.readouterr().out
        assert "disabled" in out

    def test_drop_safe_import_swallows_value_error(self, qapp):
        plugin, _ = _enc_plugin()
        plugin.on_import = MagicMock(side_effect=ValueError("cancelled"))
        assert plugin.on_drop("project.pmeenc") is True
        qapp.processEvents()
        plugin.on_import.assert_called_once_with("project.pmeenc")

    def test_drop_safe_import_shows_critical_on_unexpected_error(self, qapp):
        plugin, _ = _enc_plugin()
        plugin.on_import = MagicMock(side_effect=RuntimeError("boom"))
        with mock.patch.object(_enc.QMessageBox, "critical") as crit:
            assert plugin.on_drop("project.pmeenc") is True
            qapp.processEvents()
            crit.assert_called_once()
            assert "boom" in crit.call_args.args[-1]

    def test_patch_method_disconnect_failure_is_swallowed(self, patched_widget):
        """Disconnecting an action with no existing connection raises TypeError,
        which _patch_method should silently swallow before connecting the patch."""
        plugin, w, action, calls = patched_widget
        # Disconnect the auto-connected slot set up by addAction/creation so the
        # subsequent explicit disconnect() inside _patch_method has nothing to
        # remove and raises.
        plugin._patch_method(
            "save_project",
            lambda mw, *a, **k: calls.append("patched"),
            ["Save Project"],
        )
        # A second, independent action targeting the *original* (now-replaced)
        # method forces the disconnect() call inside _patch_method to fail.
        plugin._unpatch_method("save_project", ["Save Project"])
        plugin._patch_method(
            "save_project",
            lambda mw, *a, **k: calls.append("patched2"),
            ["Save Project"],
        )
        action.trigger()
        assert calls == ["patched2"]


class TestPmeencExportImportFlow:
    """Real-Qt tests for on_export/export_encrypted/on_import driving the
    QInputDialog/QFileDialog/QMessageBox code paths for coverage."""

    def _plugin(self, qapp):
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.state_manager = MagicMock()
        mw.init_manager = MagicMock()
        mw.ui_manager = MagicMock()
        plugin, ctx = _enc_plugin(mw)
        plugin.mw.state_manager.data.atoms = [1]
        plugin.mw.state_manager.create_json_data.return_value = {"atoms": [1]}
        plugin._apply_patches = MagicMock()
        plugin.derive_key = MagicMock(return_value=b"DERIVED_KEY")
        return plugin, ctx, mw

    def test_on_export_no_data_shows_status_message(self, qapp):
        plugin, ctx, mw = self._plugin(qapp)
        mw.state_manager.data.atoms = []
        ctx.current_molecule = None
        plugin.on_export()
        ctx.show_status_message.assert_called_once_with("Error: Nothing to save.")

    def test_on_export_cancelled_file_dialog_returns_early(self, qapp, tmp_path):
        plugin, ctx, mw = self._plugin(qapp)
        mw.init_manager.current_file_path = str(tmp_path / "proj.pmeprj")
        with mock.patch.object(
            _enc.QFileDialog, "getSaveFileName", return_value=("", "")
        ):
            plugin.on_export()
        assert not (tmp_path / "proj.pmeenc").exists()

    def test_on_export_appends_extension_and_encrypts(self, qapp, tmp_path):
        plugin, ctx, mw = self._plugin(qapp)
        mw.init_manager.current_file_path = None
        target = tmp_path / "myproj"
        with mock.patch.object(
            _enc.QFileDialog, "getSaveFileName", return_value=(str(target), "")
        ), mock.patch.object(
            _enc.QInputDialog, "getText", return_value=("pw", True)
        ), mock.patch.object(
            _enc, "Fernet"
        ) as fernet_cls:
            fernet_obj = MagicMock()
            fernet_obj.encrypt.return_value = b"CIPHER"
            fernet_cls.return_value = fernet_obj
            plugin.on_export()
        out = target.with_suffix(".pmeenc")
        assert out.exists()
        assert out.read_bytes()[16:] == b"CIPHER"

    def test_export_encrypted_empty_password_warns_then_succeeds(self, qapp, tmp_path):
        plugin, ctx, mw = self._plugin(qapp)
        target = tmp_path / "proj.pmeenc"
        mw.init_manager.current_file_path = None
        answers = iter([("", True), ("pw", True)])
        with mock.patch.object(
            _enc.QInputDialog, "getText", side_effect=lambda *a, **k: next(answers)
        ), mock.patch.object(
            _enc.QMessageBox, "warning"
        ) as warn, mock.patch.object(
            _enc, "Fernet"
        ) as fernet_cls:
            fernet_obj = MagicMock()
            fernet_obj.encrypt.return_value = b"CIPHER"
            fernet_cls.return_value = fernet_obj
            plugin.export_encrypted(str(target))
        warn.assert_called_once()
        assert target.exists()

    def test_export_encrypted_exception_shows_critical(self, qapp, tmp_path):
        plugin, ctx, mw = self._plugin(qapp)
        target = tmp_path / "proj.pmeenc"
        mw.init_manager.current_file_path = None
        mw.state_manager.create_json_data.side_effect = RuntimeError("nope")
        with mock.patch.object(
            _enc.QInputDialog, "getText", return_value=("pw", True)
        ), mock.patch.object(_enc.QMessageBox, "critical") as crit:
            plugin.export_encrypted(str(target))
        crit.assert_called_once()
        assert not target.exists()

    def test_on_import_empty_password_warns_then_retries(self, qapp, tmp_path):
        plugin, ctx, mw = self._plugin(qapp)
        target = tmp_path / "in.pmeenc"
        target.write_bytes(b"S" * 16 + b"CIPHER")
        answers = iter([("", True), ("pw", True)])
        with mock.patch.object(
            _enc.QInputDialog, "getText", side_effect=lambda *a, **k: next(answers)
        ), mock.patch.object(
            _enc.QMessageBox, "warning"
        ) as warn, mock.patch.object(
            _enc, "Fernet"
        ) as fernet_cls:
            fernet_obj = MagicMock()
            fernet_obj.decrypt.return_value = json.dumps({"a": 1}).encode()
            fernet_cls.return_value = fernet_obj
            plugin.on_import(str(target))
        warn.assert_called_once()
        mw.state_manager.load_from_json_data.assert_called_once_with({"a": 1})

    def test_on_import_decrypt_failure_warns_and_retries(self, qapp, tmp_path):
        plugin, ctx, mw = self._plugin(qapp)
        target = tmp_path / "in.pmeenc"
        target.write_bytes(b"S" * 16 + b"CIPHER")
        answers = iter([("wrongpw", True), ("rightpw", True)])
        with mock.patch.object(
            _enc.QInputDialog, "getText", side_effect=lambda *a, **k: next(answers)
        ), mock.patch.object(
            _enc.QMessageBox, "warning"
        ) as warn, mock.patch.object(
            _enc, "Fernet"
        ) as fernet_cls:
            fernet_obj = MagicMock()
            fernet_obj.decrypt.side_effect = [Exception("bad"), json.dumps({"b": 2}).encode()]
            fernet_cls.return_value = fernet_obj
            plugin.on_import(str(target))
        warn.assert_called_once()
        assert "Decryption Failed" in warn.call_args.args[1]
        mw.state_manager.load_from_json_data.assert_called_once_with({"b": 2})

    def test_on_import_unexpected_error_shows_critical_and_reraises(self, qapp, tmp_path):
        plugin, ctx, mw = self._plugin(qapp)
        target = tmp_path / "in.pmeenc"
        target.write_bytes(b"S" * 16 + b"CIPHER")
        plugin.derive_key = MagicMock(side_effect=RuntimeError("kaboom"))
        with mock.patch.object(
            _enc.QInputDialog, "getText", return_value=("pw", True)
        ), mock.patch.object(_enc.QMessageBox, "critical") as crit:
            with pytest.raises(RuntimeError, match="kaboom"):
                plugin.on_import(str(target))
        crit.assert_called_once()
