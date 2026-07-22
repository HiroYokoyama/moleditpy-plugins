"""
Headless GUI tests for the Mapped Cube Viewer plugin.

Covers: MappedCubeSetupDialog, MappedWidget, run_plugin(), run().

Chemistry libs (pyvista, numpy, rdkit, …) are mocked; real PyQt6 is used.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MAPPED_PATH = PLUGINS_DIR / "Mapped_Cube_Viewer" / "mapped_cube_viewer.py"

with mock_chemistry_imports():
    _mapped = load_plugin_for_gui(MAPPED_PATH)


# ===========================================================================
# MappedCubeSetupDialog  (Mapped Cube Viewer)
# ===========================================================================


class TestMappedCubeSetupDialog:
    """MappedCubeSetupDialog: file-picker dialog for surface + property cubes."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _mapped.MappedCubeSetupDialog(parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Mapped Cube Setup"

    def test_surface_file_initially_none(self, dlg):
        assert dlg.surface_file is None

    def test_property_file_initially_none(self, dlg):
        assert dlg.property_file is None

    def test_surface_line_edit_empty(self, dlg):
        assert dlg.le_surf.text() == ""

    def test_property_line_edit_empty(self, dlg):
        assert dlg.le_prop.text() == ""

    def test_line_edits_accept_text(self, dlg):
        dlg.le_surf.setText("/some/surf.cube")
        dlg.le_prop.setText("/some/prop.cube")
        assert dlg.le_surf.text() == "/some/surf.cube"
        assert dlg.le_prop.text() == "/some/prop.cube"


# ===========================================================================
# MappedCubeSetupDialog — browse_surf / browse_prop / accept (real Qt)
# ===========================================================================


class TestMappedCubeSetupDialogBrowse:
    @pytest.fixture
    def dlg(self, qapp):
        d = _mapped.MappedCubeSetupDialog(parent=None)
        yield d
        d.destroy()

    def test_browse_surf_sets_line_edit(self, dlg, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog

        monkeypatch.setattr(
            QFileDialog, "getOpenFileName", lambda *a, **k: ("/picked/surf.cube", "")
        )
        dlg.browse_surf()
        assert dlg.le_surf.text() == "/picked/surf.cube"

    def test_browse_surf_cancelled_leaves_line_edit_empty(self, dlg, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog

        monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: ("", ""))
        dlg.browse_surf()
        assert dlg.le_surf.text() == ""

    def test_browse_prop_sets_line_edit(self, dlg, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog

        monkeypatch.setattr(
            QFileDialog, "getOpenFileName", lambda *a, **k: ("/picked/prop.cube", "")
        )
        dlg.browse_prop()
        assert dlg.le_prop.text() == "/picked/prop.cube"

    def test_browse_prop_cancelled_leaves_line_edit_empty(self, dlg, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog

        monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: ("", ""))
        dlg.browse_prop()
        assert dlg.le_prop.text() == ""


class TestMappedCubeSetupDialogAcceptReal:
    @pytest.fixture
    def dlg(self, qapp):
        d = _mapped.MappedCubeSetupDialog(parent=None)
        yield d
        d.destroy()

    @pytest.fixture(autouse=True)
    def _no_modal_warnings(self, monkeypatch):
        from PyQt6.QtWidgets import QMessageBox

        self.warning = MagicMock()
        monkeypatch.setattr(QMessageBox, "warning", self.warning)

    def test_missing_surface_file_warns(self, dlg, tmp_path):
        dlg.le_surf.setText(str(tmp_path / "nope1.cube"))
        dlg.le_prop.setText(str(tmp_path / "nope2.cube"))
        dlg.accept()
        self.warning.assert_called_once()
        assert "Surface file not found" in self.warning.call_args.args[2]

    def test_missing_property_file_warns(self, dlg, tmp_path):
        surf = tmp_path / "surf.cube"
        surf.write_text("x")
        dlg.le_surf.setText(str(surf))
        dlg.le_prop.setText(str(tmp_path / "nope.cube"))
        dlg.accept()
        self.warning.assert_called_once()
        assert "Property file not found" in self.warning.call_args.args[2]

    def test_both_files_present_accepts(self, dlg, tmp_path):
        surf = tmp_path / "surf.cube"
        prop = tmp_path / "prop.cube"
        surf.write_text("x")
        prop.write_text("y")
        dlg.le_surf.setText(str(surf))
        dlg.le_prop.setText(str(prop))
        dlg.accept()
        self.warning.assert_not_called()
        assert dlg.surface_file == str(surf)
        assert dlg.property_file == str(prop)
        assert dlg.result() == 1  # QDialog.Accepted


# ===========================================================================
# MappedWidget — real construction / update_mesh / export / close
# ===========================================================================


class _Vals:
    """Minimal array-like stand-in (real numpy is mocked in this suite)."""

    def __init__(self, values):
        self._v = list(values)

    def __len__(self):
        return len(self._v)

    def max(self):
        return max(self._v)

    def min(self):
        return min(self._v)

    def mean(self):
        return sum(self._v) / len(self._v)


class _FakeContour:
    def __init__(self, n_points, sample_result=None):
        self.n_points = n_points
        self._sample_result = sample_result

    def sample(self, other):
        return self._sample_result


class _FakeMapped:
    def __init__(self, values):
        self.point_data = {"property_values": _Vals(values)}


def _make_mw():
    from PyQt6.QtWidgets import QWidget

    mw = QWidget()
    mw.addDockWidget = MagicMock()
    return mw


def _make_grids(surf_vals=(0.05, 0.5, 1.0), prop_vals=(-0.2, 0.3),
                 contour_n_points=5, mapped_vals=(1.0, 2.0)):
    grid_surf = MagicMock()
    grid_surf.point_data.get.return_value = _Vals(surf_vals)
    mapped = _FakeMapped(mapped_vals)
    grid_surf.contour.return_value = _FakeContour(contour_n_points, sample_result=mapped)

    grid_prop = MagicMock()
    grid_prop.point_data.get.return_value = _Vals(prop_vals)
    return grid_surf, grid_prop


class TestMappedWidgetConstruction:
    @pytest.fixture
    def widget(self, qapp):
        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        dock = MagicMock()
        grid_surf, grid_prop = _make_grids()
        w = _mapped.MappedWidget(ctx, dock, grid_surf, grid_prop)
        yield w, ctx, mw, dock
        w.destroy()

    def test_creates_without_error(self, widget):
        w, *_ = widget
        assert w is not None

    def test_iso_default_uses_high_iso_when_surf_max_above_threshold(self, widget):
        w, *_ = widget
        assert w.iso_spin.value() == pytest.approx(0.002)

    def test_min_max_spin_reflect_property_range(self, qapp):
        # contour_n_points=0 keeps the __init__-time auto_fit update_mesh() a
        # no-op, so min/max stay at the values set directly from grid_prop.
        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf, grid_prop = _make_grids(contour_n_points=0)
        w = _mapped.MappedWidget(ctx, MagicMock(), grid_surf, grid_prop)
        assert w.min_spin.value() == pytest.approx(-0.2)
        assert w.max_spin.value() == pytest.approx(0.3)
        w.destroy()

    def test_cmap_combo_default(self, widget):
        w, *_ = widget
        assert w.cmap_combo.currentText() == "jet_r"

    def test_actor_set_after_init_update_mesh(self, widget):
        w, ctx, *_ = widget
        # add_mesh() called during the auto_fit update_mesh() at end of __init__
        ctx.plotter.add_mesh.assert_called()
        assert w.actor is ctx.plotter.add_mesh.return_value

    def test_opacity_default(self, widget):
        w, *_ = widget
        assert w.opacity_spin.value() == pytest.approx(0.4)

    def test_transparent_checkbox_checked_by_default(self, widget):
        w, *_ = widget
        assert w.check_transparent.isChecked()


class TestMappedWidgetInitDefaultBranches:
    """Empty / low-max grid values exercise the else branches in __init__."""

    def test_empty_surf_values_uses_fallback_iso(self, qapp):
        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf, grid_prop = _make_grids(surf_vals=())
        w = _mapped.MappedWidget(ctx, MagicMock(), grid_surf, grid_prop)
        assert w.iso_spin.value() == pytest.approx(0.002)
        w.destroy()

    def test_low_max_surf_values_uses_mean_iso(self, qapp, monkeypatch):
        # np is mocked in this suite; MagicMock's default __float__ returns
        # 1.0, which would mask the real `float(np.mean(vals))` computation.
        monkeypatch.setattr(_mapped.np, "mean", lambda v: sum(v._v) / len(v._v))
        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf, grid_prop = _make_grids(surf_vals=(0.01, 0.02, 0.03))
        w = _mapped.MappedWidget(ctx, MagicMock(), grid_surf, grid_prop)
        assert w.iso_spin.value() == pytest.approx(0.02)
        w.destroy()

    def test_empty_prop_values_uses_fallback_range(self, qapp):
        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf, grid_prop = _make_grids(prop_vals=(), contour_n_points=0)
        w = _mapped.MappedWidget(ctx, MagicMock(), grid_surf, grid_prop)
        assert w.min_spin.value() == pytest.approx(-0.1)
        assert w.max_spin.value() == pytest.approx(0.1)
        w.destroy()


class TestMappedWidgetUpdateMeshSignals:
    @pytest.fixture
    def widget(self, qapp):
        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf, grid_prop = _make_grids()
        w = _mapped.MappedWidget(ctx, MagicMock(), grid_surf, grid_prop)
        yield w, ctx
        w.destroy()

    def test_iso_spin_change_triggers_update_mesh(self, widget):
        w, ctx = widget
        ctx.plotter.reset_mock()
        w.iso_spin.setValue(0.5)
        ctx.plotter.add_mesh.assert_called()

    def test_cmap_change_triggers_update_mesh(self, widget):
        w, ctx = widget
        ctx.plotter.reset_mock()
        w.cmap_combo.setCurrentText("viridis")
        ctx.plotter.add_mesh.assert_called()
        assert ctx.plotter.add_mesh.call_args.kwargs["cmap"] == "viridis"

    def test_min_spin_change_triggers_update_mesh(self, widget):
        w, ctx = widget
        ctx.plotter.reset_mock()
        w.min_spin.setValue(-3.0)
        ctx.plotter.add_mesh.assert_called()

    def test_max_spin_change_triggers_update_mesh(self, widget):
        w, ctx = widget
        ctx.plotter.reset_mock()
        w.max_spin.setValue(3.0)
        ctx.plotter.add_mesh.assert_called()

    def test_opacity_change_triggers_update_mesh(self, widget):
        w, ctx = widget
        ctx.plotter.reset_mock()
        w.opacity_spin.setValue(0.9)
        ctx.plotter.add_mesh.assert_called()
        assert ctx.plotter.add_mesh.call_args.kwargs["opacity"] == pytest.approx(0.9)

    def test_fit_button_triggers_auto_fit_update(self, widget):
        w, ctx = widget
        ctx.plotter.reset_mock()
        fit_btn = [
            b for b in w.findChildren(_mapped.QPushButton)
            if b.text() == "Fit Range to Surface"
        ][0]
        fit_btn.click()
        ctx.plotter.add_mesh.assert_called()

    def test_empty_contour_result_is_noop(self, qapp):
        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf, grid_prop = _make_grids(contour_n_points=0)
        w = _mapped.MappedWidget(ctx, MagicMock(), grid_surf, grid_prop)
        ctx.plotter.add_mesh.assert_not_called()
        w.destroy()

    def test_exception_during_update_mesh_is_caught(self, qapp):
        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf = MagicMock()
        grid_surf.point_data.get.return_value = _Vals([0.05, 0.5])
        grid_surf.contour.side_effect = RuntimeError("boom")
        grid_prop = MagicMock()
        grid_prop.point_data.get.return_value = _Vals([-0.1, 0.1])
        w = _mapped.MappedWidget(ctx, MagicMock(), grid_surf, grid_prop)  # must not raise
        ctx.plotter.add_mesh.assert_not_called()
        w.destroy()


class TestMappedWidgetClosePlugin:
    def _widget(self, qapp, mw_extra=None):
        mw = _make_mw()
        if mw_extra:
            mw_extra(mw)
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        dock = MagicMock()
        grid_surf, grid_prop = _make_grids()
        w = _mapped.MappedWidget(ctx, dock, grid_surf, grid_prop)
        return w, ctx, mw, dock

    def test_close_plugin_removes_actor_and_closes_dock(self, qapp):
        def _extra(mw):
            mw.ui_manager = MagicMock()
            mw.edit_actions_manager = MagicMock()

        w, ctx, mw, dock = self._widget(qapp, _extra)
        actor_before_close = w.actor
        ctx.plotter.reset_mock()
        w.close_plugin()
        ctx.plotter.remove_actor.assert_called_once_with(actor_before_close)
        ctx.plotter.render.assert_called()
        mw.ui_manager.restore_ui_for_editing.assert_called_once()
        dock.close.assert_called_once()
        mw.edit_actions_manager.clear_all.assert_called_once()

    def test_close_plugin_without_ui_manager_does_not_crash(self, qapp):
        w, ctx, mw, dock = self._widget(qapp)
        w.close_plugin()  # mw has no ui_manager / edit_actions_manager -> silenced warnings
        dock.close.assert_called_once()

    def test_close_plugin_without_dock_does_not_crash(self, qapp):
        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf, grid_prop = _make_grids()
        w = _mapped.MappedWidget(ctx, None, grid_surf, grid_prop)
        w.close_plugin()  # no dock -> skip dock.close()


class TestMappedWidgetExportView:
    @pytest.fixture
    def widget(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog, QMessageBox

        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf, grid_prop = _make_grids()
        w = _mapped.MappedWidget(ctx, MagicMock(), grid_surf, grid_prop)

        self.info = MagicMock()
        self.critical = MagicMock()
        monkeypatch.setattr(QMessageBox, "information", self.info)
        monkeypatch.setattr(QMessageBox, "critical", self.critical)
        self._save_ret = ("", "")
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: self._save_ret)

        yield w, ctx
        w.destroy()

    def test_cancel_does_nothing(self, widget):
        w, ctx = widget
        ctx.plotter.reset_mock()
        w.export_view()
        ctx.plotter.screenshot.assert_not_called()

    def test_success_shows_information(self, widget, tmp_path):
        w, ctx = widget
        f = str(tmp_path / "out.png")
        self._save_ret = (f, "")
        w.export_view()
        ctx.plotter.screenshot.assert_called_once_with(
            f, transparent_background=True
        )
        self.info.assert_called_once()

    def test_screenshot_failure_shows_critical(self, widget, tmp_path):
        w, ctx = widget
        ctx.plotter.screenshot.side_effect = RuntimeError("disk full")
        f = str(tmp_path / "out.png")
        self._save_ret = (f, "")
        w.export_view()
        self.critical.assert_called_once()

    def test_transparent_flag_follows_checkbox(self, widget, tmp_path):
        w, ctx = widget
        w.check_transparent.setChecked(False)
        f = str(tmp_path / "out.png")
        self._save_ret = (f, "")
        w.export_view()
        ctx.plotter.screenshot.assert_called_once_with(
            f, transparent_background=False
        )


class TestMappedWidgetExportColorbar:
    @pytest.fixture
    def widget(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog, QMessageBox

        mw = _make_mw()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        grid_surf, grid_prop = _make_grids()
        w = _mapped.MappedWidget(ctx, MagicMock(), grid_surf, grid_prop)

        self.info = MagicMock()
        self.critical = MagicMock()
        monkeypatch.setattr(QMessageBox, "information", self.info)
        monkeypatch.setattr(QMessageBox, "critical", self.critical)
        self._save_ret = ("", "")
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: self._save_ret)

        yield w
        w.destroy()

    def test_cancel_does_nothing(self, widget):
        widget.export_colorbar()
        self.info.assert_not_called()
        self.critical.assert_not_called()

    def test_success_shows_information(self, widget, tmp_path):
        f = str(tmp_path / "cbar.png")
        self._save_ret = (f, "")
        widget.export_colorbar()
        self.info.assert_called_once()

    def test_transparent_background_selected(self, widget, tmp_path):
        widget.check_transparent.setChecked(True)
        f = str(tmp_path / "cbar.png")
        self._save_ret = (f, "")
        widget.export_colorbar()
        self.info.assert_called_once()

    def test_non_transparent_background_selected(self, widget, tmp_path):
        widget.check_transparent.setChecked(False)
        f = str(tmp_path / "cbar.png")
        self._save_ret = (f, "")
        widget.export_colorbar()
        self.info.assert_called_once()

    def test_plotter_construction_failure_shows_critical(self, widget, tmp_path, monkeypatch):
        monkeypatch.setattr(
            _mapped.pv, "Plotter", MagicMock(side_effect=RuntimeError("no display"))
        )
        f = str(tmp_path / "cbar.png")
        self._save_ret = (f, "")
        widget.export_colorbar()
        self.critical.assert_called_once()


# ===========================================================================
# run_plugin() — full entry-point flow
# ===========================================================================


class _FakeDialog:
    """Stand-in for MappedCubeSetupDialog, fully controllable."""

    _accepted = True
    _surface_file = "s.cube"
    _property_file = "p.cube"

    def __init__(self, parent=None):
        pass

    def exec(self):
        from PyQt6.QtWidgets import QDialog

        return (
            QDialog.DialogCode.Accepted
            if self._accepted
            else QDialog.DialogCode.Rejected
        )

    @property
    def surface_file(self):
        return self._surface_file

    @property
    def property_file(self):
        return self._property_file


def _fake_meta_grid(name):
    grid = MagicMock()
    grid.point_data.get.return_value = _Vals([0.05, 0.5, 1.0])
    grid.contour.return_value = _FakeContour(1, sample_result=_FakeMapped([1.0, 2.0]))
    meta = {"atoms": [(6, (0.0, 0.0, 0.0))]}
    return meta, grid


class TestRunPlugin:
    @pytest.fixture(autouse=True)
    def _patch_dialog_and_read_cube(self, monkeypatch):
        _FakeDialog._accepted = True
        _FakeDialog._surface_file = "s.cube"
        _FakeDialog._property_file = "p.cube"
        monkeypatch.setattr(_mapped, "MappedCubeSetupDialog", _FakeDialog)

        def _fake_read_cube(f):
            return _fake_meta_grid(f)

        monkeypatch.setattr(_mapped, "read_cube", _fake_read_cube)

    def _ctx(self, mw, old_dock=None):
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.get_window.return_value = old_dock
        return ctx

    def test_dialog_cancelled_does_nothing(self, qapp):
        _FakeDialog._accepted = False
        mw = _make_mw()
        ctx = self._ctx(mw)
        _mapped.run_plugin(ctx)
        ctx.enter_3d_mode.assert_not_called()
        ctx.register_window.assert_not_called()

    def test_missing_files_returns_early(self, qapp):
        _FakeDialog._surface_file = ""
        mw = _make_mw()
        ctx = self._ctx(mw)
        _mapped.run_plugin(ctx)
        ctx.register_window.assert_not_called()

    def test_success_path_registers_dock(self, qapp):
        mw = _make_mw()
        ctx = self._ctx(mw)
        _mapped.run_plugin(ctx)
        ctx.enter_3d_mode.assert_called_once()
        mw.addDockWidget.assert_called_once()
        ctx.register_window.assert_called_once()
        assert ctx.register_window.call_args.args[0] == "main_panel"

    def test_old_dock_with_close_plugin_widget_is_closed(self, qapp):
        mw = _make_mw()
        old_widget = MagicMock()
        old_widget.close_plugin = MagicMock()
        old_dock = MagicMock()
        old_dock.widget.return_value = old_widget
        ctx = self._ctx(mw, old_dock=old_dock)
        _mapped.run_plugin(ctx)
        old_widget.close_plugin.assert_called_once()

    def test_old_dock_without_close_plugin_is_closed_directly(self, qapp):
        mw = _make_mw()
        old_widget = MagicMock(spec=[])  # no close_plugin attribute
        old_dock = MagicMock()
        old_dock.widget.return_value = old_widget
        ctx = self._ctx(mw, old_dock=old_dock)
        _mapped.run_plugin(ctx)
        old_dock.close.assert_called_once()

    def test_read_cube_error_shows_status_message(self, qapp, monkeypatch):
        def _boom(f):
            raise RuntimeError("bad cube")

        monkeypatch.setattr(_mapped, "read_cube", _boom)
        mw = _make_mw()
        ctx = self._ctx(mw)
        _mapped.run_plugin(ctx)  # must not raise
        ctx.show_status_message.assert_called_once()
        assert "Mapped Cube Viewer Error" in ctx.show_status_message.call_args.args[0]


# ===========================================================================
# run(mw) — legacy entry point
# ===========================================================================


class TestRunEntryPoint:
    def test_noop_without_plugin_manager(self, qapp):
        mw = MagicMock(spec=[])  # no plugin_manager attribute
        _mapped.run(mw)  # must not raise

    def test_delegates_to_run_plugin(self, qapp, monkeypatch):
        called = MagicMock()
        monkeypatch.setattr(_mapped, "run_plugin", called)
        mw = MagicMock()  # has plugin_manager (auto MagicMock attribute)
        with mock_chemistry_imports():
            _mapped.run(mw)
        called.assert_called_once()
