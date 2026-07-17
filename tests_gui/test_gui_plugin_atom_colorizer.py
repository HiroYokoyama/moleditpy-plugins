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
