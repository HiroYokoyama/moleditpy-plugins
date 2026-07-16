"""
Headless GUI tests for the Blender Export plugin.

Covers module constants, initialize(), and run() (no dialog class; export
logic is callable-only).

All tests use real PyQt6 (QT_QPA_PLATFORM=offscreen).
Chemistry libraries are mocked via mock_chemistry_imports().
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

BLENDER_PATH = PLUGINS_DIR / "Blender_Export" / "blender_export.py"

with mock_chemistry_imports():
    _blender = load_plugin_for_gui(BLENDER_PATH)


def _make_ctx() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# Blender Export — module-level and initialize()/run() tests
# (no dialog class; export logic is callable-only)
# ===========================================================================


class TestBlenderExportModule:
    """Module constants, initialize(), and run() for Blender Export."""

    def test_plugin_name(self):
        assert _blender.PLUGIN_NAME == "Blender Export"

    def test_plugin_version_is_string(self):
        assert isinstance(_blender.PLUGIN_VERSION, str)

    def test_plugin_context_initially_none(self):
        assert _blender.PLUGIN_CONTEXT is None

    def test_initialize_registers_export_action(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_blender(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        label = ctx.add_export_action.call_args.args[0]
        assert "Blender" in label

    def test_initialize_callback_is_callable(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        cb = ctx.add_export_action.call_args.args[1]
        assert callable(cb)

    def test_run_noop_without_plugin_manager(self):
        mw = MagicMock(spec=[])  # no plugin_manager attribute
        _blender.run(mw)  # must not raise

    def test_run_noop_when_context_is_none(self):
        _blender.PLUGIN_CONTEXT = None
        mw = MagicMock()  # has plugin_manager (MagicMock auto-attribute)
        _blender.run(mw)  # PLUGIN_CONTEXT is None → early return, no error


# ===========================================================================
# export_to_blender — control flow (guards / cancel / save / error)
# ===========================================================================


def _mol(n_atoms=3, n_confs=1):
    mol = MagicMock()
    mol.GetNumAtoms.return_value = n_atoms
    mol.GetNumConformers.return_value = n_confs
    return mol


class TestBlenderExportFlow:
    @pytest.fixture(autouse=True)
    def _qt(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog, QMessageBox

        self.warning = MagicMock()
        self.info = MagicMock()
        self.critical = MagicMock()
        monkeypatch.setattr(QMessageBox, "warning", self.warning)
        monkeypatch.setattr(QMessageBox, "information", self.info)
        monkeypatch.setattr(QMessageBox, "critical", self.critical)
        self._save_ret = ("", "")
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: self._save_ret
        )
        monkeypatch.setattr(
            _blender,
            "generate_blender_script",
            MagicMock(return_value="# blender script"),
        )

    def _ctx(self, mol):
        ctx = MagicMock()
        ctx.get_main_window.return_value = MagicMock()
        ctx.current_molecule = mol
        _blender.PLUGIN_CONTEXT = ctx
        return ctx

    def test_no_molecule_warns_and_skips_generate(self):
        ctx = self._ctx(None)
        _blender.export_to_blender(ctx)
        self.warning.assert_called_once()
        _blender.generate_blender_script.assert_not_called()

    def test_empty_molecule_warns(self):
        ctx = self._ctx(_mol(n_atoms=0))
        _blender.export_to_blender(ctx)
        self.warning.assert_called_once()

    def test_no_conformer_warns_about_3d(self):
        ctx = self._ctx(_mol(n_confs=0))
        _blender.export_to_blender(ctx)
        assert "3D" in self.warning.call_args.args[1]

    def test_cancel_writes_nothing(self, tmp_path):
        ctx = self._ctx(_mol())
        self._save_ret = ("", "")
        _blender.export_to_blender(ctx)
        _blender.generate_blender_script.assert_not_called()
        assert list(tmp_path.iterdir()) == []

    def test_save_appends_py_extension_and_writes(self, tmp_path):
        ctx = self._ctx(_mol())
        target = tmp_path / "scene"
        self._save_ret = (str(target), "")
        _blender.export_to_blender(ctx)
        written = tmp_path / "scene.py"
        assert written.exists()
        assert written.read_text(encoding="utf-8") == "# blender script"
        self.info.assert_called_once()
        ctx.show_status_message.assert_called()

    def test_generation_error_shows_critical(self, tmp_path):
        ctx = self._ctx(_mol())
        self._save_ret = (str(tmp_path / "x.py"), "")
        _blender.generate_blender_script.side_effect = RuntimeError("boom")
        _blender.export_to_blender(ctx)
        self.critical.assert_called_once()
