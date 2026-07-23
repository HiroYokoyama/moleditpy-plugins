"""
Headless GUI tests for the Atom Colorizer plugin.

Covers: AtomColorizerWindow.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ATOM_COLORIZER_PATH = PLUGINS_DIR / "Atom_Colorizer" / "atom_colorizer.py"

with mock_chemistry_imports():
    _atom_colorizer = load_plugin_for_gui(ATOM_COLORIZER_PATH)


# ===========================================================================
# AtomColorizerWindow  (visible plugin: "Atom Colorizer")
# ===========================================================================


def _colorizer_context() -> MagicMock:
    """Minimal stub: get_main_window() returns None so no QObject parent issues."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


class TestAtomColorizerWindow:
    """AtomColorizerWindow with a no-main-window context."""

    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "Atom Colorizer"

    def test_indices_placeholder(self, win):
        assert win.le_indices.placeholderText() == "e.g. 0, 1, 5"

    def test_choose_color_button_exists(self, win):
        assert win.btn_color.text() == "Choose Color"

    def test_timer_is_active(self, win):
        assert win.sel_timer.isActive()

    def test_is_non_modal(self, win):
        assert not win.isModal()

    def test_default_color_is_red(self, win):
        assert win.current_color.name() == "#ff0000"

    def test_registered_as_main_panel(self, win):
        win.context.register_window.assert_any_call("main_panel", win)


class TestChooseColor:
    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

    def test_valid_color_updates_state_and_style(self, win, monkeypatch):
        from PyQt6.QtGui import QColor

        monkeypatch.setattr(
            _atom_colorizer.QColorDialog,
            "getColor",
            staticmethod(lambda *a, **k: QColor("#00ff00")),
        )
        win.choose_color()
        assert win.current_color.name() == "#00ff00"
        assert "#00ff00" in win.btn_color.styleSheet()

    def test_cancelled_dialog_keeps_current_color(self, win, monkeypatch):
        from PyQt6.QtGui import QColor

        monkeypatch.setattr(
            _atom_colorizer.QColorDialog,
            "getColor",
            staticmethod(lambda *a, **k: QColor()),  # invalid == cancelled
        )
        win.choose_color()
        assert win.current_color.name() == "#ff0000"


class TestApplyAndResetColor:
    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

    @pytest.fixture
    def warnings(self, monkeypatch):
        calls = []
        monkeypatch.setattr(
            _atom_colorizer.QMessageBox,
            "warning",
            staticmethod(lambda *a, **k: calls.append(a)),
        )
        return calls

    def test_empty_selection_warns_and_skips_controller(self, win, warnings):
        win.le_indices.setText("")
        win.apply_color()
        assert len(warnings) == 1
        assert "No atoms selected" in warnings[0][2]
        win.context.get_3d_controller.assert_not_called()

    def test_invalid_indices_warns(self, win, warnings):
        win.le_indices.setText("1, x, 3")
        win.apply_color()
        assert len(warnings) == 1
        assert "Invalid indices" in warnings[0][2]

    def test_no_molecule_warns(self, win, warnings):
        win.context.current_molecule = None
        win.le_indices.setText("0,1")
        win.apply_color()
        assert len(warnings) == 1
        assert "No molecule" in warnings[0][2]

    def test_apply_colors_in_range_indices_only(self, win, warnings):
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 3
        win.context.current_molecule = mol
        win.le_indices.setText("0, 2, 5")  # 5 is out of range
        win.apply_color()
        controller = win.context.get_3d_controller.return_value
        assert [c.args for c in controller.set_atom_color.call_args_list] == [
            (0, "#ff0000"),
            (2, "#ff0000"),
        ]
        win.context.refresh_3d_view.assert_called_once()
        assert warnings == []

    def test_reset_colors_clears_every_atom(self, win):
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 2
        win.context.current_molecule = mol
        win.reset_colors()
        controller = win.context.get_3d_controller.return_value
        assert [c.args for c in controller.set_atom_color.call_args_list] == [
            (0, None),
            (1, None),
        ]
        win.context.refresh_3d_view.assert_called_once()

    def test_reset_colors_no_molecule_is_noop(self, win):
        win.context.current_molecule = None
        win.reset_colors()
        win.context.get_3d_controller.assert_not_called()


class TestSelectionSyncAndClose:
    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

    def test_selection_from_context_sorted_into_field(self, win):
        win.context.get_selected_atom_indices.return_value = [2, 0, 1]
        win.get_selection_from_viewer()
        assert win.le_indices.text() == "0,1,2"

    def test_auto_update_selection_runs_when_field_unfocused(self, win):
        win.context.get_selected_atom_indices.return_value = [3, 1]
        win._auto_update_selection()
        assert win.le_indices.text() == "1,3"

    def test_selection_sync_error_is_silenced(self, win):
        win.context.get_selected_atom_indices.side_effect = RuntimeError("boom")
        win.le_indices.setText("keep")
        win.get_selection_from_viewer()  # must not raise
        assert win.le_indices.text() == "keep"

    def test_close_stops_timer_and_unregisters(self, win):
        assert win.sel_timer.isActive()
        win.close()
        assert not win.sel_timer.isActive()
        assert win.context.register_window.call_args.args == ("main_panel", None)


class _FakeEdit3D:
    """Plain object (not QWidget) standing in for edit_3d_manager."""

    def __init__(self):
        self.measurement_mode = False
        self.is_3d_edit_mode = False
        self.selected_atoms_3d = set()
        self.selected_atoms_for_measurement = []
        self.toggle_measurement_mode = MagicMock()
        self.clear_measurement_selection = MagicMock()
        self.update_3d_selection_display = MagicMock()


class _FakeMW:
    def __init__(self):
        self.edit_3d_manager = _FakeEdit3D()
        self.init_manager = MagicMock()
        self.ui_manager = MagicMock()


class TestEnterAndRestoreSelectMode:
    """Real (non-MagicMock parent) mw stand-in to exercise the mode
    save/restore bodies without constructing them as a QWidget parent."""

    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

    def test_enter_select_mode_no_main_window_is_noop(self, win):
        win.context.get_main_window.return_value = None
        win._enter_select_mode()  # must not raise

    def test_enter_select_mode_forces_measurement_when_off(self, win):
        fake_mw = _FakeMW()
        win.context.get_main_window.return_value = fake_mw
        win._forced_measurement_mode = False
        win._enter_select_mode()
        fake_mw.edit_3d_manager.toggle_measurement_mode.assert_called_once_with(True)
        fake_mw.init_manager.measurement_action.setChecked.assert_called_once_with(True)
        assert win._forced_measurement_mode is True

    def test_enter_select_mode_skips_toggle_when_already_active(self, win):
        fake_mw = _FakeMW()
        fake_mw.edit_3d_manager.measurement_mode = True
        win.context.get_main_window.return_value = fake_mw
        win._enter_select_mode()
        fake_mw.edit_3d_manager.toggle_measurement_mode.assert_not_called()
        assert win._restore_measurement_mode is True

    def test_enter_select_mode_exception_silenced(self, win):
        fake_mw = _FakeMW()
        fake_mw.edit_3d_manager.toggle_measurement_mode.side_effect = RuntimeError("x")
        win.context.get_main_window.return_value = fake_mw
        win._enter_select_mode()  # must not raise

    def test_restore_select_mode_no_main_window_is_noop(self, win):
        win.context.get_main_window.return_value = None
        win._restore_select_mode()  # must not raise

    def test_restore_select_mode_clears_selection_and_restores_modes(self, win):
        fake_mw = _FakeMW()
        fake_mw.edit_3d_manager.selected_atoms_3d = {1, 2}
        win.context.get_main_window.return_value = fake_mw
        win._restore_measurement_mode = True
        win._restore_edit_mode = True
        win._restore_select_mode()
        fake_mw.edit_3d_manager.clear_measurement_selection.assert_called_once()
        assert fake_mw.edit_3d_manager.selected_atoms_3d == set()
        fake_mw.edit_3d_manager.update_3d_selection_display.assert_called_once()
        fake_mw.init_manager.measurement_action.setChecked.assert_called_once_with(True)
        fake_mw.edit_3d_manager.toggle_measurement_mode.assert_called_once_with(True)
        fake_mw.ui_manager.toggle_3d_edit_mode.assert_called_once_with(True)

    def test_restore_select_mode_selection_clear_exception_silenced(self, win):
        fake_mw = _FakeMW()
        fake_mw.edit_3d_manager.clear_measurement_selection.side_effect = RuntimeError("x")
        win.context.get_main_window.return_value = fake_mw
        win._restore_select_mode()  # must not raise
        # restore-mode block still runs despite the selection-clear exception
        fake_mw.edit_3d_manager.toggle_measurement_mode.assert_called_once()

    def test_restore_select_mode_toggle_exception_silenced(self, win):
        fake_mw = _FakeMW()
        fake_mw.edit_3d_manager.toggle_measurement_mode.side_effect = RuntimeError("x")
        win.context.get_main_window.return_value = fake_mw
        win._restore_select_mode()  # must not raise


class TestGetSelectionFromViewer3D:
    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

    def test_merges_3d_and_measurement_picks(self, win):
        fake_mw = _FakeMW()
        fake_mw.edit_3d_manager.selected_atoms_3d = {2, 5}
        fake_mw.edit_3d_manager.selected_atoms_for_measurement = (7, 5)
        win.context.get_main_window.return_value = fake_mw
        win.context.get_selected_atom_indices.return_value = [0, 2]
        win.get_selection_from_viewer()
        assert win.le_indices.text() == "0,2,5,7"

    def test_non_collection_measurement_attribute_ignored(self, win):
        fake_mw = _FakeMW()
        fake_mw.edit_3d_manager.selected_atoms_for_measurement = "garbage"
        win.context.get_main_window.return_value = fake_mw
        win.context.get_selected_atom_indices.return_value = [4]
        win.get_selection_from_viewer()
        assert win.le_indices.text() == "4"


class TestApplyResetExceptionPaths:
    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

    @pytest.fixture
    def criticals(self, monkeypatch):
        calls = []
        monkeypatch.setattr(
            _atom_colorizer.QMessageBox,
            "critical",
            staticmethod(lambda *a, **k: calls.append(a)),
        )
        return calls

    def test_apply_color_controller_exception_shows_critical(self, win, criticals):
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 2
        win.context.current_molecule = mol
        win.context.get_3d_controller.side_effect = RuntimeError("boom")
        win.le_indices.setText("0")
        win.apply_color()
        assert len(criticals) == 1
        assert "Failed to apply color" in criticals[0][2]

    def test_reset_colors_controller_exception_shows_critical(self, win, criticals):
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 2
        win.context.current_molecule = mol
        win.context.get_3d_controller.side_effect = RuntimeError("boom")
        win.reset_colors()
        assert len(criticals) == 1
        assert "Failed to reset colors" in criticals[0][2]


class TestCloseEventTimerExceptionSilenced:
    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.destroy()

    def test_stop_timer_exception_is_silenced(self, win):
        bad_timer = MagicMock()
        bad_timer.isActive.side_effect = RuntimeError("boom")
        win.sel_timer = bad_timer
        win.close()  # must not raise
        assert win.context.register_window.call_args.args == ("main_panel", None)


class TestLaunch:
    def test_launch_creates_new_window_when_none_registered(self, qapp):
        ctx = _colorizer_context()
        ctx.get_window.return_value = None
        _atom_colorizer.launch(ctx)
        win = ctx.register_window.call_args.args[1]
        assert isinstance(win, _atom_colorizer.AtomColorizerWindow)
        win.sel_timer.stop()
        win.close()

    def test_launch_reuses_existing_window(self, qapp):
        ctx = _colorizer_context()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _atom_colorizer.launch(ctx)
        existing.show.assert_called_once()
        existing.raise_.assert_called_once()
        existing.activateWindow.assert_called_once()


class TestRunLegacyEntryPoint:
    def test_run_without_plugin_manager_is_noop(self, qapp):
        class _NoPluginManagerMW:
            pass

        _atom_colorizer.run(_NoPluginManagerMW())  # must not raise / launch

    def test_run_with_context_launches(self, qapp, monkeypatch):
        ctx = _colorizer_context()
        _atom_colorizer.PLUGIN_CONTEXT = ctx
        ctx.get_window.return_value = None
        called = []
        monkeypatch.setattr(_atom_colorizer, "launch", lambda c: called.append(c))

        class _MW:
            plugin_manager = MagicMock()

        _atom_colorizer.run(_MW())
        assert called == [ctx]
        _atom_colorizer.PLUGIN_CONTEXT = None
