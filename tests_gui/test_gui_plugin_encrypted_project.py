"""
Headless GUI tests for the Encrypted Project plugin.

Covers: PmeencPlugin (patch/unpatch with real QActions).
"""

from __future__ import annotations

from pathlib import Path
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
