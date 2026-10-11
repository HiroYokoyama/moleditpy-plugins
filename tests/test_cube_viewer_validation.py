"""Shared regressions for the independently installed Cube viewers."""

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from conftest import load_plugin, mock_optional_imports, extract_function


ROOT = Path(__file__).resolve().parents[1]
PATHS = (
    ROOT / "plugins/Cube_File_Viewer/cube_viewer.py",
    ROOT / "plugins/Cube_File_Viewer_Advanced/cube_viewer_advanced.py",
    ROOT / "plugins/Mapped_Cube_Viewer/mapped_cube_viewer.py",
)


@pytest.fixture(params=PATHS, ids=lambda path: path.parent.name)
def module(request):
    with mock_optional_imports():
        mod = load_plugin(request.param)
    mod.np = np
    return mod


def cube_text(data="1 2\n3 4\n5 6\n7 8\n", natoms=1, nval=None, ids=""):
    nval_field = "" if nval is None else f" {nval}"
    return (
        f"Cube\nComment\n{natoms} 0 0 0{nval_field}\n"
        "2 1 0 0\n2 0 1 0\n2 0 0 1\n"
        "6 6 0 0 0\n" + ids + data
    )


def parse(module, tmp_path, text):
    path = tmp_path / "example.cube"
    path.write_text(text, encoding="utf-8")
    return module.parse_cube_data(path)


@pytest.mark.parametrize("width", (1, 2, 3, 5, 6, 8))
def test_data_wrapping_does_not_change_scalars(module, tmp_path, width):
    data = "\n".join(" ".join(str(value) for value in range(start, min(9, start + width))) for start in range(1, 9, width))
    meta = parse(module, tmp_path, cube_text(data))
    np.testing.assert_array_equal(meta["data_flat"], np.arange(1, 9))


@pytest.mark.parametrize("data", ("1 2", "1 2 3 4 5 6 7 8 9", "1 2 3 4 5 6 7 bad", "1 2 3 4 5 6 7 nan", "1 2 3 4 5 6 7 inf"))
def test_invalid_scalar_payload_is_rejected(module, tmp_path, data):
    with pytest.raises(ValueError):
        parse(module, tmp_path, cube_text(data))


def test_positive_nval_is_deinterleaved(module, tmp_path):
    data = " ".join(str(value) for pair in zip(range(1, 9), range(101, 109)) for value in pair)
    meta = parse(module, tmp_path, cube_text(data, nval=2))
    np.testing.assert_array_equal(meta["data_flat"], np.arange(1, 9))


def test_negative_natoms_wrapped_ids_are_after_atoms(module, tmp_path):
    data = " ".join(str(value) for pair in zip(range(1, 9), range(101, 109)) for value in pair)
    meta = parse(module, tmp_path, cube_text(data, natoms=-1, ids="2 42\n43\n"))
    assert meta["atoms"][0][0] == 6
    np.testing.assert_array_equal(meta["data_flat"], np.arange(1, 9))


@pytest.mark.parametrize("ids", ("", "0\n", "-2 1 2\n", "2 1\n", "1 42 43\n"))
def test_invalid_dataset_ids_are_rejected(module, tmp_path, ids):
    with pytest.raises(ValueError):
        parse(module, tmp_path, cube_text(natoms=-1, ids=ids))


def test_fortran_d_exponents_are_supported(module, tmp_path):
    meta = parse(module, tmp_path, cube_text(" ".join(f"{value}.0D+00" for value in range(1, 9))).replace("6 6 0 0 0", "6 6 1D-2 0 0"))
    assert meta["atoms"][0][1][0] == pytest.approx(0.01)
    np.testing.assert_array_equal(meta["data_flat"], np.arange(1, 9))


@pytest.mark.parametrize("old,new", (("2 1 0 0", "0 1 0 0"), ("2 0 1 0", "2 1 0 0"), ("2 1 0 0", "2000001 1 0 0"), ("6 6 0 0 0", "6 6 nan 0 0")))
def test_invalid_or_excessive_grid_and_atoms_are_rejected(module, tmp_path, old, new):
    with pytest.raises(ValueError):
        parse(module, tmp_path, cube_text().replace(old, new))


@pytest.mark.parametrize("path", PATHS[:2], ids=lambda path: path.parent.name)
def test_failed_read_does_not_close_existing_viewer(path):
    with mock_optional_imports():
        mod = load_plugin(path)
    context = MagicMock()
    old_dock = context.get_window.return_value
    mod.read_cube = MagicMock(side_effect=ValueError("invalid data"))
    mod.open_cube_viewer(context, "bad.cube")
    old_dock.widget.assert_not_called()
    context.get_main_window.return_value.plotter.clear.assert_not_called()
    context.register_window.assert_not_called()


@pytest.mark.parametrize("path", PATHS[:2], ids=lambda path: path.parent.name)
def test_close_does_not_erase_a_replacement_molecule(path):
    close = extract_function(path, "CubeViewerWidget", "close_plugin", {"logging": MagicMock()})
    widget = MagicMock()
    widget._bound_mol = object()
    widget.mw.current_mol = object()
    close(widget)
    widget._detach_stale_view.assert_called_once()
    widget.plotter.clear.assert_not_called()


def test_excess_data_conversion_is_bounded(module, tmp_path, monkeypatch):
    converted = []
    real_fromiter = np.fromiter

    def track_values(values, dtype):
        def tracked():
            for value in values:
                converted.append(value)
                yield value
        return real_fromiter(tracked(), dtype=dtype)

    monkeypatch.setattr(np, "fromiter", track_values)
    with pytest.raises(ValueError, match="size mismatch"):
        parse(module, tmp_path, cube_text(" ".join(str(value) for value in range(10000))))
    assert len(converted) == 9
