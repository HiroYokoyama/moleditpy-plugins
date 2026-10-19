"""Real VTK interpolation and Qt lifecycle checks for the Cube viewers."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from PySide6.QtWidgets import QFileDialog, QMessageBox, QWidget

np = pytest.importorskip("numpy")
pv = pytest.importorskip("pyvista")
pytest.importorskip("rdkit")

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def module():
    path = ROOT / "plugins/Mapped_Cube_Viewer/mapped_cube_viewer.py"
    spec = importlib.util.spec_from_file_location("mapped_cube_real_data", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def grids():
    indices = np.indices((3, 3, 3)).astype(float) - 1
    points = np.column_stack([axis.flatten(order="F") for axis in indices])
    surface = pv.StructuredGrid()
    surface.points = points
    surface.dimensions = (3, 3, 3)
    surface["property_values"] = np.linalg.norm(points, axis=1)
    prop = surface.copy()
    prop["property_values"] = points @ np.array([1., 2., 3.])
    return surface, prop


@pytest.fixture
def viewer(module, qapp):
    parent = QWidget()
    plotter = MagicMock()
    plotter.add_mesh.return_value = object()
    context = SimpleNamespace(get_main_window=lambda: parent, plotter=plotter, current_molecule=None)
    surface, prop = grids()
    widget = module.MappedWidget(context, None, surface, prop)
    yield widget, context
    widget.deleteLater()
    parent.deleteLater()


def test_samples_property_values_instead_of_density(viewer):
    widget, context = viewer
    widget.iso_spin.setValue(0.75)
    mapped = context.plotter.add_mesh.call_args.args[0]
    np.testing.assert_allclose(mapped["property_values"], mapped.points @ [1., 2., 3.], atol=1e-7)
    assert np.all(mapped["vtkValidPointMask"] == 1)
    assert widget.warning_label.text() == ""


def test_empty_contour_removes_previous_surface(viewer):
    widget, context = viewer
    old_actor = widget.actor
    context.plotter.reset_mock()
    widget.iso_spin.setValue(999)
    context.plotter.remove_actor.assert_called_once_with(old_actor)
    context.plotter.add_mesh.assert_not_called()
    assert widget.actor is None
    assert "No surface" in widget.warning_label.text()


def test_uncovered_property_points_are_not_colored_as_zero(viewer):
    widget, context = viewer
    widget.grid_prop = widget.grid_prop.translate([10, 0, 0])
    context.plotter.reset_mock()
    widget.update_mesh(auto_fit=True)
    context.plotter.add_mesh.assert_not_called()
    assert widget.actor is None
    assert "does not cover" in widget.warning_label.text()


def test_invalid_color_range_removes_surface_and_reports_error(viewer):
    widget, context = viewer
    widget.min_spin.setValue(10)
    assert widget.actor is None
    assert "minimum" in widget.warning_label.text()


def test_mapped_view_detaches_when_molecule_changes(viewer):
    widget, context = viewer
    mol = object()
    context.current_molecule = mol
    widget.bind_structure(mol)
    context.current_molecule = object()
    widget._check_structure_still_bound()
    assert widget.actor is None
    assert not widget._structure_watch_timer.isActive()


def test_closing_widget_removes_surface(viewer):
    widget, context = viewer
    old_actor = widget.actor
    context.plotter.reset_mock()
    widget.close()
    context.plotter.remove_actor.assert_called_once_with(old_actor)
    assert widget.actor is None


def test_export_colorbar_closes_plotter_when_save_fails(module, viewer, monkeypatch):
    widget, _ = viewer
    plotter = MagicMock()
    plotter.screenshot.side_effect = OSError("disk full")
    monkeypatch.setattr(module.pv, "Plotter", lambda **kwargs: plotter)
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *args, **kwargs: ("out.png", ""))
    critical = MagicMock()
    monkeypatch.setattr(QMessageBox, "critical", critical)
    widget.export_colorbar()
    critical.assert_called_once()
    plotter.close.assert_called_once()


def test_export_view_reports_filesystem_error(viewer, monkeypatch):
    widget, context = viewer
    context.plotter.screenshot.side_effect = OSError("disk full")
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *args, **kwargs: ("out.png", ""))
    critical = MagicMock()
    monkeypatch.setattr(QMessageBox, "critical", critical)
    widget.export_view()
    critical.assert_called_once()


@pytest.mark.parametrize("directory,filename", (("Cube_File_Viewer", "cube_viewer.py"), ("Cube_File_Viewer_Advanced", "cube_viewer_advanced.py"), ("Mapped_Cube_Viewer", "mapped_cube_viewer.py")))
def test_grid_values_and_skewed_coordinates_stay_aligned(directory, filename):
    path = ROOT / "plugins" / directory / filename
    spec = importlib.util.spec_from_file_location("cube_grid_real_data", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    dims = (2, 3, 4)
    vectors = np.array([[1., .1, 0], [.2, 2, 0], [0, .3, 3]])
    origin = np.array([.5, -.5, 1])
    data = np.arange(24)
    meta = dict(dims=dims, origin=origin, x_vec=vectors[0], y_vec=vectors[1], z_vec=vectors[2], atoms=[], data_flat=data, is_angstrom_header=True)
    _, grid = mod.build_grid_from_meta(meta)
    key = "property_values" if directory == "Mapped_Cube_Viewer" else "values"
    for z in range(4):
        for y in range(3):
            for x in range(2):
                row = x + 2 * y + 6 * z
                np.testing.assert_allclose(grid.points[row], origin + np.array([x, y, z]) @ vectors)
                assert grid[key][row] == data[x * 12 + y * 4 + z]


@pytest.mark.parametrize("directory,filename", (("Cube_File_Viewer", "cube_viewer.py"), ("Cube_File_Viewer_Advanced", "cube_viewer_advanced.py"), ("Mapped_Cube_Viewer", "mapped_cube_viewer.py")))
def test_title_bar_close_runs_plugin_cleanup(directory, filename, qapp):
    path = ROOT / "plugins" / directory / filename
    spec = importlib.util.spec_from_file_location("cube_dock_real_qt", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    dock = mod.CubeViewerDock("Cube")
    child = QWidget()
    child.close_plugin = MagicMock()
    dock.setWidget(child)
    dock.show()
    assert dock.close()
    child.close_plugin.assert_called_once_with()
    dock.deleteLater()
