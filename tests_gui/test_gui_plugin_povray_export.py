"""
Headless GUI tests for the POV-Ray Export plugin.

Covers module constants, initialize(), and run().

All tests use real PyQt6 (QT_QPA_PLATFORM=offscreen).
Chemistry libraries are mocked via mock_chemistry_imports().
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

POVRAY_PATH = PLUGINS_DIR / "POV-Ray_Export" / "povray_export.py"

with mock_chemistry_imports():
    _povray = load_plugin_for_gui(POVRAY_PATH)


def _make_ctx() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# POV-Ray Export — module-level and initialize()/run() tests
# ===========================================================================


class TestPovRayExportModule:
    """Module constants, initialize(), and run() for POV-Ray Export."""

    def test_plugin_name(self):
        assert _povray.PLUGIN_NAME == "POV-Ray Export"

    def test_plugin_version_is_string(self):
        assert isinstance(_povray.PLUGIN_VERSION, str)

    def test_plugin_context_initially_none(self):
        assert _povray.PLUGIN_CONTEXT is None

    def test_initialize_registers_export_action(self):
        ctx = _make_ctx()
        _povray.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_povray(self):
        ctx = _make_ctx()
        _povray.initialize(ctx)
        label = ctx.add_export_action.call_args.args[0]
        assert "POV" in label

    def test_initialize_callback_is_callable(self):
        ctx = _make_ctx()
        _povray.initialize(ctx)
        cb = ctx.add_export_action.call_args.args[1]
        assert callable(cb)

    def test_run_noop_without_plugin_manager(self):
        mw = MagicMock(spec=[])
        _povray.run(mw)  # must not raise

    def test_run_noop_when_context_is_none(self):
        _povray.PLUGIN_CONTEXT = None
        mw = MagicMock()
        _povray.run(mw)  # PLUGIN_CONTEXT is None → early return


# ===========================================================================
# export_to_povray — control flow (guards / cancel / save / error)
# ===========================================================================


def _mol(n_atoms=3, n_confs=1):
    mol = MagicMock()
    mol.GetNumAtoms.return_value = n_atoms
    mol.GetNumConformers.return_value = n_confs
    return mol


class TestPovRayExportFlow:
    @pytest.fixture(autouse=True)
    def _qt(self, qapp, monkeypatch):
        # export_to_povray imports these from PyQt6.QtWidgets at call time.
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
        # Isolate the flow from the heavy geometry generator.
        monkeypatch.setattr(
            _povray, "generate_povray_scene", MagicMock(return_value=("SCENE", 800, 600))
        )

    def _ctx(self, mol):
        ctx = MagicMock()
        ctx.get_main_window.return_value = MagicMock()
        ctx.current_molecule = mol
        _povray.PLUGIN_CONTEXT = ctx
        return ctx

    def test_no_molecule_warns_and_skips_generate(self):
        ctx = self._ctx(None)
        _povray.export_to_povray(ctx)
        self.warning.assert_called_once()
        _povray.generate_povray_scene.assert_not_called()

    def test_empty_molecule_warns(self):
        ctx = self._ctx(_mol(n_atoms=0))
        _povray.export_to_povray(ctx)
        self.warning.assert_called_once()

    def test_no_conformer_warns_about_3d(self):
        ctx = self._ctx(_mol(n_confs=0))
        _povray.export_to_povray(ctx)
        assert "3D" in self.warning.call_args.args[1]

    def test_cancel_writes_nothing(self, tmp_path):
        ctx = self._ctx(_mol())
        self._save_ret = ("", "")
        _povray.export_to_povray(ctx)
        _povray.generate_povray_scene.assert_not_called()
        assert list(tmp_path.iterdir()) == []

    def test_save_appends_pov_extension_and_writes(self, tmp_path):
        ctx = self._ctx(_mol())
        target = tmp_path / "scene"
        self._save_ret = (str(target), "")
        _povray.export_to_povray(ctx)
        written = tmp_path / "scene.pov"
        assert written.exists()
        assert written.read_text(encoding="utf-8") == "SCENE"
        self.info.assert_called_once()
        ctx.show_status_message.assert_called()

    def test_generation_error_shows_critical(self, tmp_path):
        ctx = self._ctx(_mol())
        self._save_ret = (str(tmp_path / "x.pov"), "")
        _povray.generate_povray_scene.side_effect = RuntimeError("boom")
        _povray.export_to_povray(ctx)
        self.critical.assert_called_once()
