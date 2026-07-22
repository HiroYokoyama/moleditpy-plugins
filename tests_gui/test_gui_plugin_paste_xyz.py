"""
Headless GUI tests for the Paste XYZ plugin.

Covers: PasteXYZDialog.

These tests use real PyQt6 (QT_QPA_PLATFORM=offscreen) rather than mocking
all Qt classes.  Chemistry/scientific libraries (rdkit, numpy, …) are still
replaced with MagicMock via mock_chemistry_imports() so no installed chemistry
stack is required.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from PyQt6.QtWidgets import QDialog, QDialogButtonBox, QLineEdit, QPushButton, QWidget

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PASTE_XYZ_PATH = PLUGINS_DIR / "Paste_XYZ" / "paste_xyz.py"

with mock_chemistry_imports():
    _paste_xyz = load_plugin_for_gui(PASTE_XYZ_PATH)


def _real_mw(io_manager=None):
    """A real QWidget usable as a QDialog parent, with an optional io_manager
    attribute duck-typed on (MagicMock() itself is not a valid QWidget parent)."""
    mw = QWidget()
    if io_manager is not None:
        mw.io_manager = io_manager
    return mw


def _click(dialog, text):
    """Click the first QPushButton child whose label matches *text*."""
    for b in dialog.findChildren(QPushButton):
        if b.text() == text:
            b.click()
            return
    raise AssertionError(f"button {text!r} not found in dialog")


@pytest.fixture
def patch_exec(monkeypatch):
    """Replace QDialog.exec() with a queue of scripts run against the real
    dialog instance instead of blocking the (offscreen) event loop."""
    scripts = []

    def _exec(self):
        script = scripts.pop(0)
        script(self)
        return self.result()

    monkeypatch.setattr(QDialog, "exec", _exec)

    def _register(script):
        scripts.append(script)

    return _register


# ===========================================================================
# PasteXYZDialog  (visible plugin: "Paste XYZ")
# ===========================================================================


class TestPasteXYZDialog:
    """PasteXYZDialog — pure-Qt dialog; rdkit guarded by try/except."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _paste_xyz.PasteXYZDialog(parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Paste XYZ"

    def test_default_size_wide_enough(self, dlg):
        assert dlg.width() >= 400

    def test_get_data_initially_empty(self, dlg):
        assert dlg.get_data() == ""

    def test_get_data_returns_typed_text(self, dlg):
        dlg.text_edit.setPlainText("C 0.0 0.0 0.0\nH 1.0 0.0 0.0")
        result = dlg.get_data()
        assert "C 0.0 0.0 0.0" in result
        assert "H 1.0 0.0 0.0" in result

    def test_multiline_xyz_round_trip(self, dlg):
        xyz = "C 0.0 0.0 0.0\nH 1.089 0.0 0.0\nH -0.363 1.027 0.0"
        dlg.text_edit.setPlainText(xyz)
        assert dlg.get_data().strip() == xyz

    def test_load_button_label(self, dlg):
        assert dlg.load_btn.text() == "Load"

    def test_cancel_button_label(self, dlg):
        assert dlg.cancel_btn.text() == "Cancel"

    def test_text_edit_placeholder_mentions_xyz(self, dlg):
        ph = dlg.text_edit.placeholderText().lower()
        assert "xyz" in ph

    def test_clear_resets_get_data(self, dlg):
        dlg.text_edit.setPlainText("something")
        dlg.text_edit.clear()
        assert dlg.get_data() == ""


# ===========================================================================
# _build_rwmol / _apply_bonds — pure functions, real Qt not required but
# exercised here (rather than tests/) so coverage tracks the real module.
# ===========================================================================


class TestBuildRwmol:
    def test_creates_atom_and_conformer_per_entry(self, qapp):
        atoms = [("C", 0.0, 0.0, 0.0), ("H", 1.0, 0.0, 0.0)]
        rwmol = _paste_xyz._build_rwmol(atoms)
        assert rwmol is not None
        assert _paste_xyz.Chem.Atom.call_count >= 2


class TestApplyBonds:
    def test_success_path_sets_charge(self, qapp):
        rwmol = MagicMock()
        mw = MagicMock()
        result = _paste_xyz._apply_bonds(rwmol, 1, mw)
        result.SetIntProp.assert_called_once_with("_xyz_charge", 1)

    def test_fallback_when_determine_bonds_raises(self, qapp, monkeypatch):
        monkeypatch.setattr(_paste_xyz.Chem.RWMol, "side_effect", RuntimeError("boom"))
        rwmol = MagicMock()
        mw = MagicMock()
        result = _paste_xyz._apply_bonds(rwmol, -1, mw)
        mw.io_manager.estimate_bonds_from_distances.assert_called_once_with(rwmol)
        result.SetIntProp.assert_called_once_with("_xyz_charge", -1)

    def test_fallback_without_io_manager(self, qapp, monkeypatch):
        monkeypatch.setattr(_paste_xyz.Chem.RWMol, "side_effect", RuntimeError("boom"))
        rwmol = MagicMock()
        result = _paste_xyz._apply_bonds(rwmol, 2, SimpleNamespace())
        result.SetIntProp.assert_called_once_with("_xyz_charge", 2)

    def test_fallback_estimator_exception_swallowed(self, qapp, monkeypatch):
        monkeypatch.setattr(_paste_xyz.Chem.RWMol, "side_effect", RuntimeError("boom"))
        rwmol = MagicMock()
        mw = MagicMock()
        mw.io_manager.estimate_bonds_from_distances.side_effect = ValueError("bad geom")
        result = _paste_xyz._apply_bonds(rwmol, 0, mw)  # must not raise
        result.SetIntProp.assert_called_once_with("_xyz_charge", 0)


# ===========================================================================
# _prompt_charge — real QDialog, driven via patch_exec.
# ===========================================================================


class TestPromptCharge:
    def test_valid_charge_accepted(self, qapp, patch_exec):
        def script(dialog):
            dialog.findChild(QLineEdit).setText("-2")
            dialog.findChild(QDialogButtonBox).button(QDialogButtonBox.StandardButton.Ok).click()

        patch_exec(script)
        assert _paste_xyz._prompt_charge(None) == (-2, True, False)

    def test_invalid_text_keeps_dialog_open_then_valid_accepted(self, qapp, patch_exec):
        def script(dialog):
            line_edit = dialog.findChild(QLineEdit)
            box = dialog.findChild(QDialogButtonBox)
            line_edit.setText("abc")
            box.button(QDialogButtonBox.StandardButton.Ok).click()
            assert dialog.result() != QDialog.DialogCode.Accepted
            line_edit.setText("3")
            box.button(QDialogButtonBox.StandardButton.Ok).click()

        patch_exec(script)
        assert _paste_xyz._prompt_charge(None) == (3, True, False)

    def test_cancel_returns_not_ok(self, qapp, patch_exec):
        def script(dialog):
            dialog.findChild(QDialogButtonBox).button(QDialogButtonBox.StandardButton.Cancel).click()

        patch_exec(script)
        assert _paste_xyz._prompt_charge(None) == (0, False, False)

    def test_skip_button_returns_skip(self, qapp, patch_exec):
        def script(dialog):
            _click(dialog, "Skip chemistry")

        patch_exec(script)
        assert _paste_xyz._prompt_charge(None) == (0, True, True)


# ===========================================================================
# _resolve_mol_with_charge
# ===========================================================================


class TestResolveMolWithCharge:
    def test_default_path_no_prompt(self, qapp):
        ctx = MagicMock()
        ctx.get_setting.return_value = False
        rwmol = MagicMock()
        mw = MagicMock()
        result = _paste_xyz._resolve_mol_with_charge(rwmol, ctx, mw)
        result.SetIntProp.assert_called_once_with("_xyz_charge", 0)

    def test_cancel_prompt_returns_none(self, qapp, patch_exec):
        ctx = MagicMock()
        ctx.get_setting.return_value = True  # always_ask_charge

        def script(dialog):
            dialog.findChild(QDialogButtonBox).button(QDialogButtonBox.StandardButton.Cancel).click()

        patch_exec(script)
        assert _paste_xyz._resolve_mol_with_charge(MagicMock(), ctx, _real_mw()) is None

    def test_skip_uses_estimator(self, qapp, patch_exec):
        ctx = MagicMock()
        ctx.get_setting.return_value = True
        io_manager = MagicMock()
        mw = _real_mw(io_manager)
        rwmol = MagicMock()

        def script(dialog):
            _click(dialog, "Skip chemistry")

        patch_exec(script)
        result = _paste_xyz._resolve_mol_with_charge(rwmol, ctx, mw)
        io_manager.estimate_bonds_from_distances.assert_called_once_with(rwmol)
        result.SetIntProp.assert_called_once_with("_xyz_skip_checks", 1)

    def test_retry_after_failed_charge_then_success(self, qapp, patch_exec, monkeypatch):
        ctx = MagicMock()
        ctx.get_setting.return_value = True
        monkeypatch.setattr(_paste_xyz.Chem.RWMol, "side_effect", RuntimeError("bad"))
        rwmol = MagicMock()
        rwmol.GetMol.side_effect = [RuntimeError("boom"), MagicMock()]

        def script(dialog):
            dialog.findChild(QLineEdit).setText("1")
            dialog.findChild(QDialogButtonBox).button(QDialogButtonBox.StandardButton.Ok).click()

        patch_exec(script)
        patch_exec(script)
        result = _paste_xyz._resolve_mol_with_charge(rwmol, ctx, _real_mw())
        ctx.show_status_message.assert_called_once()
        assert result is not None


# ===========================================================================
# run_plugin — end-to-end through the real dialogs.
# ===========================================================================


class TestRunPlugin:
    def test_cancel_paste_dialog_leaves_canvas(self, qapp, patch_exec):
        ctx = MagicMock()
        ctx.get_main_window.return_value = _real_mw()

        def script(dialog):
            _click(dialog, "Cancel")

        patch_exec(script)
        _paste_xyz.run_plugin(ctx)
        ctx.clear_canvas.assert_not_called()

    def test_empty_text_returns(self, qapp, patch_exec):
        ctx = MagicMock()
        ctx.get_main_window.return_value = _real_mw()

        def script(dialog):
            dialog.text_edit.setPlainText("   ")
            _click(dialog, "Load")

        patch_exec(script)
        _paste_xyz.run_plugin(ctx)
        ctx.clear_canvas.assert_not_called()

    def test_invalid_xyz_shows_warning(self, qapp, patch_exec, monkeypatch):
        warn = MagicMock()
        monkeypatch.setattr(_paste_xyz.QMessageBox, "warning", warn)
        ctx = MagicMock()
        ctx.get_main_window.return_value = _real_mw()

        def script(dialog):
            dialog.text_edit.setPlainText("this has no coordinates")
            _click(dialog, "Load")

        patch_exec(script)
        _paste_xyz.run_plugin(ctx)
        warn.assert_called_once()
        ctx.clear_canvas.assert_not_called()

    def test_skip_checks_with_io_manager(self, qapp, patch_exec):
        ctx = MagicMock()
        io_manager = MagicMock()
        ctx.get_main_window.return_value = _real_mw(io_manager)
        ctx.get_setting.side_effect = lambda key, default=None: (
            True if key == "skip_chemistry_checks" else default
        )

        def script(dialog):
            dialog.text_edit.setPlainText("C 0 0 0")
            _click(dialog, "Load")

        patch_exec(script)
        _paste_xyz.run_plugin(ctx)
        ctx.clear_canvas.assert_called_once_with(push_to_undo=False)
        io_manager.estimate_bonds_from_distances.assert_called_once()

    def test_skip_checks_without_io_manager(self, qapp, patch_exec):
        ctx = MagicMock()
        ctx.get_main_window.return_value = _real_mw()  # no io_manager attr
        ctx.get_setting.side_effect = lambda key, default=None: (
            True if key == "skip_chemistry_checks" else default
        )

        def script(dialog):
            dialog.text_edit.setPlainText("C 0 0 0")
            _click(dialog, "Load")

        patch_exec(script)
        _paste_xyz.run_plugin(ctx)  # must not raise despite no io_manager
        ctx.clear_canvas.assert_called_once_with(push_to_undo=False)

    def test_charge_prompt_cancel_leaves_canvas(self, qapp, patch_exec):
        ctx = MagicMock()
        ctx.get_main_window.return_value = _real_mw()
        ctx.get_setting.side_effect = lambda key, default=None: {
            "skip_chemistry_checks": False,
            "always_ask_charge": True,
        }.get(key, default)

        def paste_script(dialog):
            dialog.text_edit.setPlainText("C 0 0 0")
            _click(dialog, "Load")

        def charge_script(dialog):
            dialog.findChild(QDialogButtonBox).button(QDialogButtonBox.StandardButton.Cancel).click()

        patch_exec(paste_script)
        patch_exec(charge_script)
        _paste_xyz.run_plugin(ctx)
        ctx.clear_canvas.assert_not_called()
        ctx.push_undo_checkpoint.assert_not_called()

    def test_full_success_with_charge_prompt(self, qapp, patch_exec):
        ctx = MagicMock()
        ctx.get_main_window.return_value = _real_mw()
        ctx.get_setting.side_effect = lambda key, default=None: {
            "skip_chemistry_checks": False,
            "always_ask_charge": True,
        }.get(key, default)

        def paste_script(dialog):
            dialog.text_edit.setPlainText("C 0 0 0\nH 1 0 0")
            _click(dialog, "Load")

        def charge_script(dialog):
            dialog.findChild(QLineEdit).setText("0")
            dialog.findChild(QDialogButtonBox).button(QDialogButtonBox.StandardButton.Ok).click()

        patch_exec(paste_script)
        patch_exec(charge_script)
        _paste_xyz.run_plugin(ctx)
        ctx.clear_canvas.assert_called_once_with(push_to_undo=False)
        ctx.push_undo_checkpoint.assert_called_once()
        ctx.check_chemistry_problems.assert_called_once()
        ctx.refresh_ui.assert_called_once()
        ctx.enter_3d_mode.assert_called_once()
        ctx.reset_3d_camera.assert_called_once()
        ctx.show_status_message.assert_called_once()

    def test_exception_during_parse_shows_critical(self, qapp, patch_exec, monkeypatch):
        monkeypatch.setattr(_paste_xyz, "parse_xyz_lines", MagicMock(side_effect=RuntimeError("boom")))
        crit = MagicMock()
        monkeypatch.setattr(_paste_xyz.QMessageBox, "critical", crit)
        ctx = MagicMock()
        ctx.get_main_window.return_value = _real_mw()

        def script(dialog):
            dialog.text_edit.setPlainText("C 0 0 0")
            _click(dialog, "Load")

        patch_exec(script)
        _paste_xyz.run_plugin(ctx)
        crit.assert_called_once()
        ctx.clear_canvas.assert_not_called()

    def test_no_rdkit_shows_critical(self, qapp, monkeypatch):
        monkeypatch.setattr(_paste_xyz, "Chem", None)
        crit = MagicMock()
        monkeypatch.setattr(_paste_xyz.QMessageBox, "critical", crit)
        ctx = MagicMock()
        _paste_xyz.run_plugin(ctx)
        crit.assert_called_once()
