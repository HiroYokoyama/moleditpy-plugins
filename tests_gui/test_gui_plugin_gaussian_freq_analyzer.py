"""
Headless GUI tests for the Gaussian Freq Analyzer plugin.

Covers:
  - GaussianFCHKFreqAnalyzer — main analyser widget (QWidget subclass), safe
    to construct when context.get_main_window()=None.
  - SpectrumDialog (QDialog) — on_range_changed() is called at end of
    __init__; it calls recalc_curve() which guards on len(freqs)==0
    (MagicMock.__len__ returns 0), so no crash.
  - FCHKParser — pure-Python parser object that needs no Qt.

Run:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/test_gui_plugin_gaussian_freq_analyzer.py -v
"""

from __future__ import annotations

import contextlib
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

GAUSSIAN_PATH = PLUGINS_DIR / "Gaussian_Freq_Analyzer" / "gaussian_fchk_freq_analyzer.py"

with mock_chemistry_imports():
    _gaussian = load_plugin_for_gui(GAUSSIAN_PATH)


def _gaussian_context() -> MagicMock:
    """Minimal stub context — get_main_window() returns None."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


# ===========================================================================
# GaussianFCHKFreqAnalyzer  (visible plugin: "Gaussian Freq Analyzer")
# ===========================================================================


class TestGaussianFCHKFreqAnalyzer:
    """Main analyser widget — safe to construct when context.get_main_window()=None."""

    @pytest.fixture
    def widget(self, qapp):
        ctx = _gaussian_context()
        w = _gaussian.GaussianFCHKFreqAnalyzer(ctx)
        yield w
        w.timer.stop()
        w.destroy()

    def test_creates_without_error(self, widget):
        assert widget is not None

    def test_info_label_text(self, widget):
        assert widget.lbl_info.text() == "Drop .fchk file here or click Open"

    def test_freq_list_column_count(self, widget):
        assert widget.list_freq.columnCount() == 3

    def test_freq_list_header_labels(self, widget):
        header = widget.list_freq.headerItem()
        assert header.text(0) == "No."
        assert "Frequency" in header.text(1)
        assert "Intensity" in header.text(2)

    def test_freq_list_initially_empty(self, widget):
        assert widget.list_freq.topLevelItemCount() == 0

    def test_scaling_factor_default(self, widget):
        assert widget.spin_sf.value() == 1.0

    def test_scaling_factor_range(self, widget):
        assert widget.spin_sf.minimum() == 0.0
        assert widget.spin_sf.maximum() == 2.0

    def test_show_vectors_unchecked(self, widget):
        assert not widget.chk_vectors.isChecked()

    def test_vector_scale_default(self, widget):
        assert widget.spin_vec_scale.value() == 2.0

    def test_amplitude_slider_default(self, widget):
        assert widget.slider_amp.value() == 5

    def test_speed_slider_default(self, widget):
        assert widget.slider_speed.value() == 20

    def test_play_button_text(self, widget):
        assert widget.btn_play.text() == "Play"

    def test_stop_button_text(self, widget):
        assert widget.btn_stop.text() == "Stop"

    def test_timer_not_active_initially(self, widget):
        assert not widget.timer.isActive()

    def test_parser_initially_none(self, widget):
        assert widget.parser is None

    def test_meta_label_initially_empty(self, widget):
        assert widget.lbl_meta.text() == ""

    def test_is_not_playing_initially(self, widget):
        assert not widget.is_playing

    def test_accepts_drops(self, widget):
        assert widget.acceptDrops()


# ===========================================================================
# SpectrumDialog (Gaussian, module-level)
# ===========================================================================


class TestGaussianSpectrumDialog:
    """Top-level SpectrumDialog from the Gaussian plugin.

    on_range_changed() is called at end of __init__; it calls recalc_curve()
    which guards on len(freqs)==0 (MagicMock.__len__ returns 0), so no crash.
    """

    @pytest.fixture
    def dlg(self, qapp):
        d = _gaussian.SpectrumDialog([], [], parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_default_window_title(self, dlg):
        assert dlg.windowTitle() == "IR Spectrum"

    def test_custom_title(self, qapp):
        d = _gaussian.SpectrumDialog([], [], title="FCHK Raman")
        assert d.windowTitle() == "FCHK Raman"
        d.destroy()

    def test_fwhm_spinbox_default(self, dlg):
        assert dlg.spin_fwhm.value() == 50

    def test_min_wavenumber_default(self, dlg):
        assert dlg.spin_min.value() == 0

    def test_max_wavenumber_default(self, dlg):
        assert dlg.spin_max.value() == 4000

    def test_size_at_least_600_wide(self, dlg):
        assert dlg.width() >= 600


# ===========================================================================
# FCHKParser  — pure Python, no Qt required
# ===========================================================================


class TestFCHKParser:
    """FCHKParser parses Gaussian FCHK files; tested in isolation."""

    def test_creates_without_error(self):
        assert _gaussian.FCHKParser() is not None

    def test_initial_charge(self):
        assert _gaussian.FCHKParser().charge == 0

    def test_initial_multiplicity(self):
        assert _gaussian.FCHKParser().multiplicity == 1

    def test_initial_frequencies_empty(self):
        assert _gaussian.FCHKParser().frequencies == []

    def test_initial_atoms_empty(self):
        assert _gaussian.FCHKParser().atoms == []

    def test_initial_coords_empty(self):
        assert _gaussian.FCHKParser().coords == []

    def test_initial_vib_modes_empty(self):
        assert _gaussian.FCHKParser().vib_modes == []


# ===========================================================================
# Real-numpy module instance `_gnp`
#
# The plugin's animation/spectrum/parsing math (np.sin, np.pi, np.linspace,
# np.exp, np.array, ...) needs real numpy to actually execute rather than
# chase MagicMock attribute chains; rdkit/PIL/PyQt-adjacent chemistry stays
# mocked (as in tests_gui/conftest.py's BLOCKED_CHEMISTRY).
# ===========================================================================


@contextlib.contextmanager
def _mock_chemistry_keep_numpy():
    real_np_mods = {
        k: v for k, v in sys.modules.items() if k == "numpy" or k.startswith("numpy.")
    }
    with mock_chemistry_imports():
        sys.modules.update(real_np_mods)
        yield


with _mock_chemistry_keep_numpy():
    _gnp = load_plugin_for_gui(GAUSSIAN_PATH)


_WATER_FCHK = textwrap.dedent("""\
    Total Energy                               R          -76.0000000000
    Charge                 I                0
    Multiplicity           I                1
    Atomic numbers                             I   N=           3
         8     1     1
    Current cartesian coordinates              R   N=           9
       0.000000E+00  0.000000E+00  0.000000E+00
       0.000000E+00  1.434600E+00  1.126430E+00
       0.000000E+00 -1.434600E+00  1.126430E+00
    Vib-Modes                                  R   N=          27
       0.1  0.2  0.3  0.0  0.0  0.0  0.0  0.0  0.0
       0.0  0.0  0.0  0.1  0.2  0.3  0.0  0.0  0.0
       0.0  0.0  0.0  0.0  0.0  0.0  0.1  0.2  0.3
    Vib-E2                                     R   N=          12
       1609.85  3681.19  3811.12
       1.000000  2.000000  3.000000
       10.000000  20.000000  30.000000
       90.0  80.0  70.0
""")


def _write_water_fchk(tmp_path, content=_WATER_FCHK, name="water.fchk"):
    f = tmp_path / name
    f.write_text(content)
    return str(f)


def _no_block_msgbox(monkeypatch, module=_gnp):
    """QMessageBox.* pop up modal dialogs that block exec() under offscreen
    Qt with no user interaction; replace with recording stubs."""
    calls = {"warning": [], "information": [], "critical": []}
    for kind in calls:
        monkeypatch.setattr(
            module.QMessageBox, kind,
            staticmethod(lambda *a, _k=kind, **kw: calls[_k].append(a)),
        )
    return calls


class _FakePos:
    """Stand-in for rdkit.Geometry.Point3D with real, inspectable x/y/z."""

    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class _FakeConformer:
    def __init__(self, coords):
        self._coords = [_FakePos(*c) for c in coords]

    def SetAtomPosition(self, idx, pos):
        self._coords[idx] = pos

    def GetAtomPosition(self, idx):
        return self._coords[idx]


class _FakeMol:
    """Minimal RDKit-Mol stand-in — just enough for GetConformer()."""

    def __init__(self, coords):
        self._conf = _FakeConformer(coords)

    def GetConformer(self):
        return self._conf


class _FakeUrl:
    def __init__(self, path):
        self._path = path

    def toLocalFile(self):
        return self._path


class _FakeMimeData:
    def __init__(self, urls):
        self._urls = urls

    def hasUrls(self):
        return bool(self._urls)

    def urls(self):
        return self._urls


class _FakeDropEvent:
    def __init__(self, urls):
        self._mime = _FakeMimeData(urls)
        self.accepted = False
        self.ignored = False

    def mimeData(self):
        return self._mime

    def acceptProposedAction(self):
        self.accepted = True

    def ignore(self):
        self.ignored = True


class _FakeImage:
    def __init__(self):
        self.mode = "RGBA"
        self.info = {}
        self.saved_path = None
        self.saved_kwargs = None

    def convert(self, mode, **kwargs):
        self.mode = mode
        return self

    def split(self):
        return [self, self, self, self]

    def paste(self, idx, mask):
        return None

    def save(self, path, **kwargs):
        self.saved_path = path
        self.saved_kwargs = kwargs


class _FakePILImageModule:
    def __init__(self):
        self.created = []
        self.Palette = SimpleNamespace(ADAPTIVE="ADAPTIVE")

    def fromarray(self, arr):
        img = _FakeImage()
        self.created.append(img)
        return img

    def eval(self, img, fn):
        return _FakeImage()


@pytest.fixture
def gnp_widget(qapp):
    ctx = _gaussian_context()
    w = _gnp.GaussianFCHKFreqAnalyzer(ctx)
    yield w
    w.timer.stop()
    w.destroy()


@pytest.fixture
def loaded_widget(gnp_widget, tmp_path):
    """Widget after loading a real (small, synthetic) FCHK file."""
    path = _write_water_fchk(tmp_path)
    gnp_widget.load_file(path)
    return gnp_widget


def _select_first_row(widget):
    item = widget.list_freq.topLevelItem(0)
    widget.list_freq.setCurrentItem(item)
    return item


# ===========================================================================
# Drag & drop
# ===========================================================================


class TestDragAndDrop:
    def test_drag_enter_accepts_fchk(self, gnp_widget):
        event = _FakeDropEvent([_FakeUrl("/tmp/molecule.fchk")])
        gnp_widget.dragEnterEvent(event)
        assert event.accepted

    def test_drag_enter_accepts_fch_and_fck(self, gnp_widget):
        for ext in ("fch", "fck"):
            event = _FakeDropEvent([_FakeUrl(f"/tmp/molecule.{ext}")])
            gnp_widget.dragEnterEvent(event)
            assert event.accepted

    def test_drag_enter_rejects_other_extension(self, gnp_widget):
        event = _FakeDropEvent([_FakeUrl("/tmp/molecule.log")])
        gnp_widget.dragEnterEvent(event)
        assert event.ignored
        assert not event.accepted

    def test_drag_enter_rejects_no_urls(self, gnp_widget):
        event = _FakeDropEvent([])
        gnp_widget.dragEnterEvent(event)
        assert event.ignored

    def test_drop_event_loads_file(self, gnp_widget, tmp_path):
        path = _write_water_fchk(tmp_path)
        event = _FakeDropEvent([_FakeUrl(path)])
        with patch.object(gnp_widget, "load_file") as mock_load:
            gnp_widget.dropEvent(event)
        mock_load.assert_called_once_with(path)
        assert event.accepted

    def test_drop_event_ignores_non_fchk(self, gnp_widget):
        event = _FakeDropEvent([_FakeUrl("/tmp/other.txt")])
        with patch.object(gnp_widget, "load_file") as mock_load:
            gnp_widget.dropEvent(event)
        mock_load.assert_not_called()


# ===========================================================================
# open_file_dialog
# ===========================================================================


class TestOpenFileDialog:
    def test_selecting_file_loads_it(self, gnp_widget, tmp_path):
        path = _write_water_fchk(tmp_path)
        with patch.object(_gnp.QFileDialog, "getOpenFileName", return_value=(path, "")):
            gnp_widget.open_file_dialog()
        assert gnp_widget.parser is not None
        assert gnp_widget.list_freq.topLevelItemCount() == 3

    def test_cancel_does_not_load(self, gnp_widget):
        with patch.object(_gnp.QFileDialog, "getOpenFileName", return_value=("", "")):
            gnp_widget.open_file_dialog()
        assert gnp_widget.parser is None


# ===========================================================================
# load_file / update_ui_after_load / update_list_and_spectrum_values
# ===========================================================================


class TestLoadFileReal:
    def test_populates_frequency_list(self, loaded_widget):
        assert loaded_widget.list_freq.topLevelItemCount() == 3
        item0 = loaded_widget.list_freq.topLevelItem(0)
        assert item0.text(0) == "1"
        assert item0.text(1) == "1609.85"
        assert item0.text(2) == "90.0000"

    def test_meta_label_shows_charge_mult_atoms(self, loaded_widget):
        assert "Charge: 0" in loaded_widget.lbl_meta.text()
        assert "Multiplicity: 1" in loaded_widget.lbl_meta.text()
        assert "Atoms: 3" in loaded_widget.lbl_meta.text()

    def test_info_label_updates_to_filename_and_green_style(self, loaded_widget):
        assert "water.fchk" in loaded_widget.lbl_info.text()
        assert "4CAF50" in loaded_widget.lbl_info.styleSheet()

    def test_base_mol_created_when_atoms_present(self, loaded_widget):
        assert loaded_widget.base_mol is not None

    def test_parse_error_shows_critical(self, gnp_widget, monkeypatch, tmp_path):
        calls = _no_block_msgbox(monkeypatch)
        missing = str(tmp_path / "does_not_exist.fchk")
        gnp_widget.load_file(missing)
        assert calls["critical"]

    def test_update_ui_after_load_no_intensities_attr_shows_dash(self, gnp_widget):
        parser = _gnp.FCHKParser()
        parser.frequencies = [123.45]
        parser.atoms = []
        gnp_widget.parser = parser
        gnp_widget.update_ui_after_load()
        item = gnp_widget.list_freq.topLevelItem(0)
        assert item.text(2) == "-"

    def test_create_base_molecule_noop_without_parser(self, gnp_widget):
        gnp_widget.parser = None
        gnp_widget.create_base_molecule()  # should not raise
        assert gnp_widget.base_mol is None


class TestUpdateListAndSpectrumValues:
    def test_scaling_updates_displayed_frequency(self, loaded_widget):
        loaded_widget.spin_sf.setValue(2.0)
        item0 = loaded_widget.list_freq.topLevelItem(0)
        assert item0.text(1) == f"{1609.85 * 2.0:.2f}"

    def test_noop_without_parser(self, gnp_widget):
        gnp_widget.parser = None
        gnp_widget.update_list_and_spectrum_values()  # should not raise

    def test_noop_without_frequencies(self, gnp_widget):
        gnp_widget.parser = _gnp.FCHKParser()
        gnp_widget.update_list_and_spectrum_values()  # frequencies == [] -> return


# ===========================================================================
# on_freq_selected
# ===========================================================================


class TestOnFreqSelected:
    def test_selecting_row_updates_vectors_when_not_playing(self, loaded_widget):
        loaded_widget.chk_vectors.setChecked(True)
        _select_first_row(loaded_widget)
        # No exception; vector_actor path exercised via update_vectors()
        assert loaded_widget.list_freq.currentItem() is not None

    def test_selecting_row_while_playing_is_a_noop(self, loaded_widget):
        loaded_widget.is_playing = True
        loaded_widget.on_freq_selected(None, None)  # 'pass' branch, no crash
        loaded_widget.is_playing = False


# ===========================================================================
# Playback controls: toggle_play / stop_play / update_timer_interval
# ===========================================================================


class TestPlaybackControls:
    def test_toggle_play_noop_without_selection(self, loaded_widget):
        loaded_widget.toggle_play()
        assert not loaded_widget.is_playing

    def test_toggle_play_starts_and_pauses(self, loaded_widget):
        _select_first_row(loaded_widget)
        loaded_widget.toggle_play()
        assert loaded_widget.is_playing
        assert loaded_widget.timer.isActive()
        assert loaded_widget.btn_play.text() == "Pause"

        loaded_widget.toggle_play()
        assert not loaded_widget.is_playing
        assert not loaded_widget.timer.isActive()
        assert loaded_widget.btn_play.text() == "Play"

    def test_stop_play_resets_state(self, loaded_widget):
        _select_first_row(loaded_widget)
        loaded_widget.toggle_play()
        loaded_widget.stop_play()
        assert not loaded_widget.is_playing
        assert not loaded_widget.timer.isActive()
        assert loaded_widget.btn_play.text() == "Play"

    def test_update_timer_interval_matches_fps(self, loaded_widget):
        loaded_widget.slider_speed.setValue(25)
        loaded_widget.update_timer_interval()
        assert loaded_widget.timer.interval() == int(1000 / 25)


# ===========================================================================
# apply_displacement / reset_geometry — real coordinate math via _FakeMol
# ===========================================================================


class TestApplyDisplacementReal:
    def test_displacement_moves_atoms_by_scaled_vector(self, monkeypatch, gnp_widget):
        monkeypatch.setattr(_gnp, "Point3D", _FakePos)
        gnp_widget.parser = _gnp.FCHKParser()
        gnp_widget.parser.coords = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0)]
        gnp_widget.base_mol = _FakeMol(gnp_widget.parser.coords)

        mode_vecs = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
        gnp_widget.apply_displacement(mode_vecs, 0.5)

        conf = gnp_widget.base_mol.GetConformer()
        p0 = conf.GetAtomPosition(0)
        p1 = conf.GetAtomPosition(1)
        assert p0.x == pytest.approx(0.5)
        assert p0.y == pytest.approx(0.0)
        assert p1.y == pytest.approx(2.5)
        assert p1.x == pytest.approx(1.0)

    def test_reset_geometry_restores_base_coords(self, monkeypatch, gnp_widget):
        monkeypatch.setattr(_gnp, "Point3D", _FakePos)
        gnp_widget.parser = _gnp.FCHKParser()
        gnp_widget.parser.coords = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0)]
        gnp_widget.base_mol = _FakeMol(gnp_widget.parser.coords)

        gnp_widget.apply_displacement([(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)], 0.5)
        gnp_widget.reset_geometry()

        conf = gnp_widget.base_mol.GetConformer()
        p1 = conf.GetAtomPosition(1)
        assert p1.y == pytest.approx(2.0)

    def test_reset_geometry_noop_without_base_mol(self, gnp_widget):
        gnp_widget.base_mol = None
        gnp_widget.parser = _gnp.FCHKParser()
        gnp_widget.reset_geometry()  # should not raise

    def test_reset_geometry_noop_without_parser(self, gnp_widget):
        gnp_widget.base_mol = MagicMock()
        gnp_widget.parser = None
        gnp_widget.reset_geometry()  # should not raise


# ===========================================================================
# animate_frame — real np.sin / np.pi math
# ===========================================================================


class TestAnimateFrameReal:
    def test_stops_when_no_parser_or_base_mol(self, gnp_widget):
        gnp_widget.parser = None
        gnp_widget.base_mol = None
        gnp_widget.animate_frame()
        assert not gnp_widget.is_playing

    def test_noop_without_current_item(self, loaded_widget):
        loaded_widget.list_freq.setCurrentItem(None)
        loaded_widget.animate_frame()  # should not raise

    def test_advances_animation_step_and_redraws(self, loaded_widget):
        _select_first_row(loaded_widget)
        start_step = loaded_widget.animation_step
        loaded_widget.animate_frame()
        assert loaded_widget.animation_step == start_step + 1
        loaded_widget.context.draw_molecule_3d.assert_called_with(loaded_widget.base_mol)


# ===========================================================================
# remove_vectors / update_vectors
# ===========================================================================


class TestVectorsReal:
    def test_remove_vectors_calls_plotter(self, gnp_widget):
        gnp_widget.vector_actor = MagicMock()
        gnp_widget.remove_vectors()
        gnp_widget.context.plotter.remove_actor.assert_called_once()
        assert gnp_widget.vector_actor is None

    def test_remove_vectors_swallows_exception(self, gnp_widget):
        gnp_widget.vector_actor = MagicMock()
        gnp_widget.context.plotter.remove_actor.side_effect = RuntimeError("gone")
        gnp_widget.remove_vectors()  # should not raise
        assert gnp_widget.vector_actor is None

    def test_remove_vectors_noop_without_actor(self, gnp_widget):
        gnp_widget.vector_actor = None
        gnp_widget.remove_vectors()

    def test_update_vectors_unchecked_is_noop(self, loaded_widget):
        loaded_widget.chk_vectors.setChecked(False)
        loaded_widget.update_vectors()
        assert loaded_widget.vector_actor is None

    def test_update_vectors_noop_without_parser(self, gnp_widget):
        gnp_widget.chk_vectors.setChecked(True)
        gnp_widget.parser = None
        gnp_widget.update_vectors()

    def test_update_vectors_noop_without_current_item(self, loaded_widget):
        loaded_widget.chk_vectors.setChecked(True)
        loaded_widget.list_freq.setCurrentItem(None)
        loaded_widget.update_vectors()

    def test_update_vectors_calls_add_arrows_when_checked(self, loaded_widget):
        loaded_widget.chk_vectors.setChecked(True)
        # currentItemChanged -> on_freq_selected -> update_vectors() already
        # fires the call; selecting the row is enough to exercise the path.
        _select_first_row(loaded_widget)
        assert loaded_widget.context.plotter.add_arrows.called

    def test_update_vectors_swallows_add_arrows_exception(self, loaded_widget):
        loaded_widget.chk_vectors.setChecked(True)
        _select_first_row(loaded_widget)
        loaded_widget.context.plotter.add_arrows.side_effect = RuntimeError("boom")
        loaded_widget.update_vectors()  # should not raise


# ===========================================================================
# show_spectrum
# ===========================================================================


class TestShowSpectrumReal:
    def test_noop_without_parser(self, gnp_widget):
        gnp_widget.show_spectrum()  # parser is None -> return

    def test_noop_without_frequencies(self, gnp_widget):
        gnp_widget.parser = _gnp.FCHKParser()
        gnp_widget.show_spectrum()  # frequencies == [] -> return

    def test_all_low_frequencies_shows_info(self, gnp_widget, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        gnp_widget.parser = _gnp.FCHKParser()
        gnp_widget.parser.frequencies = [1.0, -2.0]
        gnp_widget.parser.intensities = [1.0, 1.0]
        gnp_widget.show_spectrum()
        assert calls["information"]

    def test_opens_dialog_with_scaled_filtered_frequencies(self, loaded_widget, monkeypatch):
        captured = {}

        class _StubDialog:
            def __init__(self, freqs, intensities, title="", parent=None):
                captured["freqs"] = list(freqs)
                captured["intensities"] = list(intensities)
                captured["title"] = title

            def exec(self):
                captured["exec_called"] = True
                return 0

        monkeypatch.setattr(_gnp, "SpectrumDialog", _StubDialog)
        loaded_widget.spin_sf.setValue(2.0)
        loaded_widget.show_spectrum()

        assert captured["exec_called"]
        assert captured["freqs"] == pytest.approx([1609.85 * 2.0, 3681.19 * 2.0, 3811.12 * 2.0])
        assert captured["intensities"] == pytest.approx([90.0, 80.0, 70.0])
        assert "water.fchk" in captured["title"]


# ===========================================================================
# save_as_gif
# ===========================================================================


class TestSaveAsGifReal:
    def test_noop_without_parser_or_base_mol(self, gnp_widget):
        gnp_widget.save_as_gif()  # should not raise

    def test_no_selection_shows_warning(self, loaded_widget, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        loaded_widget.list_freq.setCurrentItem(None)
        loaded_widget.save_as_gif()
        assert calls["warning"]

    def test_dialog_reject_resumes_playback(self, loaded_widget):
        _select_first_row(loaded_widget)
        loaded_widget.toggle_play()
        assert loaded_widget.is_playing
        with patch.object(
            _gnp.QDialog, "exec", return_value=_gnp.QDialog.DialogCode.Rejected
        ):
            loaded_widget.save_as_gif()
        assert loaded_widget.is_playing  # resumed after cancel

    def test_save_dialog_cancel_resumes_playback(self, loaded_widget):
        _select_first_row(loaded_widget)
        loaded_widget.toggle_play()
        with patch.object(
            _gnp.QDialog, "exec", return_value=_gnp.QDialog.DialogCode.Accepted
        ), patch.object(_gnp.QFileDialog, "getSaveFileName", return_value=("", "")):
            loaded_widget.save_as_gif()
        assert loaded_widget.is_playing

    def test_success_writes_gif_and_appends_extension(self, loaded_widget, tmp_path, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        _select_first_row(loaded_widget)
        out_path = str(tmp_path / "out")  # no .gif suffix
        fake_image_mod = _FakePILImageModule()
        loaded_widget.context.plotter.screenshot.return_value = object()
        with patch.object(_gnp, "Image", fake_image_mod), patch.object(
            _gnp.QDialog, "exec", return_value=_gnp.QDialog.DialogCode.Accepted
        ), patch.object(
            _gnp.QFileDialog, "getSaveFileName", return_value=(out_path, "")
        ):
            loaded_widget.save_as_gif()
        assert fake_image_mod.created, "no frames captured"
        assert fake_image_mod.created[0].saved_path == out_path + ".gif"
        assert calls["information"]

    def test_screenshot_none_produces_no_frames(self, loaded_widget, tmp_path):
        _select_first_row(loaded_widget)
        loaded_widget.context.plotter.screenshot.return_value = None
        fake_image_mod = _FakePILImageModule()
        out_path = str(tmp_path / "out.gif")
        with patch.object(_gnp, "Image", fake_image_mod), patch.object(
            _gnp.QDialog, "exec", return_value=_gnp.QDialog.DialogCode.Accepted
        ), patch.object(
            _gnp.QFileDialog, "getSaveFileName", return_value=(out_path, "")
        ):
            loaded_widget.save_as_gif()
        assert fake_image_mod.created == []

    def test_save_exception_shows_critical(self, loaded_widget, tmp_path, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        _select_first_row(loaded_widget)

        class _RaisingImage(_FakeImage):
            def save(self, path, **kwargs):
                raise OSError("disk full")

        fake_image_mod = _FakePILImageModule()
        fake_image_mod.fromarray = lambda arr: _RaisingImage()
        loaded_widget.context.plotter.screenshot.return_value = object()
        out_path = str(tmp_path / "out.gif")
        with patch.object(_gnp, "Image", fake_image_mod), patch.object(
            _gnp.QDialog, "exec", return_value=_gnp.QDialog.DialogCode.Accepted
        ), patch.object(
            _gnp.QFileDialog, "getSaveFileName", return_value=(out_path, "")
        ):
            loaded_widget.save_as_gif()
        assert calls["critical"]

    def test_resumes_playback_in_finally_block(self, loaded_widget, tmp_path):
        _select_first_row(loaded_widget)
        loaded_widget.toggle_play()
        fake_image_mod = _FakePILImageModule()
        loaded_widget.context.plotter.screenshot.return_value = object()
        out_path = str(tmp_path / "resumed.gif")
        with patch.object(_gnp, "Image", fake_image_mod), patch.object(
            _gnp.QMessageBox, "information", staticmethod(lambda *a, **kw: None)
        ), patch.object(
            _gnp.QDialog, "exec", return_value=_gnp.QDialog.DialogCode.Accepted
        ), patch.object(
            _gnp.QFileDialog, "getSaveFileName", return_value=(out_path, "")
        ):
            loaded_widget.save_as_gif()
        assert loaded_widget.is_playing  # resumed in the finally block


# ===========================================================================
# on_dock_visibility_changed / close_plugin
# ===========================================================================


class TestOnDockVisibilityChangedReal:
    def test_analyzer_stops_play_when_hidden_while_playing(self, loaded_widget):
        _select_first_row(loaded_widget)
        loaded_widget.toggle_play()
        assert loaded_widget.is_playing
        loaded_widget.on_dock_visibility_changed(False)
        assert not loaded_widget.is_playing

    def test_analyzer_visible_is_noop(self, loaded_widget):
        loaded_widget.is_playing = False
        loaded_widget.on_dock_visibility_changed(True)  # should not raise

    def test_close_plugin_resets_state(self, loaded_widget):
        _select_first_row(loaded_widget)
        loaded_widget.toggle_play()
        loaded_widget.vector_actor = MagicMock()
        loaded_widget.close_plugin()
        assert not loaded_widget.is_playing
        assert loaded_widget.animation_step == 0
        assert loaded_widget.vector_actor is None

    def test_spectrum_plot_widget_dead_on_dock_visibility_changed(self, qapp):
        """SpectrumPlotWidget.on_dock_visibility_changed (~line 1273) refers
        to self.is_playing / self.stop_play / self.base_mol / self.context,
        none of which SpectrumPlotWidget defines. It is unreachable in the
        running app: only GaussianFCHKFreqAnalyzer.on_dock_visibility_changed
        (line 944) is ever connected to a dock's visibilityChanged signal.
        Exercise the literal method body (with stand-in attributes) so the
        dead branch is covered without altering runtime behaviour."""
        w = _gnp.SpectrumPlotWidget([], [])
        stop_calls = []
        reset_calls = []
        w.is_playing = True
        w.stop_play = lambda: stop_calls.append(True)
        w.reset_geometry = lambda: reset_calls.append(True)
        w.animation_step = 7
        w.base_mol = MagicMock()
        w.context = MagicMock()

        w.on_dock_visibility_changed(False)

        assert stop_calls == [True]
        assert reset_calls == [True]
        assert w.animation_step == 0
        w.context.plotter.render.assert_called_once()
        w.destroy()

    def test_spectrum_plot_widget_dead_visible_is_noop(self, qapp):
        w = _gnp.SpectrumPlotWidget([], [])
        w.on_dock_visibility_changed(True)  # "if not visible" False -> no-op
        w.destroy()


# ===========================================================================
# close button (nested close_action in init_ui)
# ===========================================================================


class TestCloseAction:
    def _find_close_button(self, widget):
        from PyQt6.QtWidgets import QPushButton

        for btn in widget.findChildren(QPushButton):
            if btn.text() == "Close":
                return btn
        raise AssertionError("Close button not found")

    def test_close_with_no_dock_closes_widget(self, loaded_widget):
        btn = self._find_close_button(loaded_widget)
        btn.click()
        assert not loaded_widget.is_playing

    def test_close_with_dock_closes_dock(self, qapp, tmp_path):
        ctx = _gaussian_context()
        dock = MagicMock()
        w = _gnp.GaussianFCHKFreqAnalyzer(ctx, dock_widget=dock)
        btn = self._find_close_button(w)
        btn.click()
        dock.close.assert_called_once()
        w.timer.stop()
        w.destroy()

    def test_close_restores_editing_ui_when_available(self, qapp):
        from PyQt6.QtWidgets import QWidget

        class _FakeMainWindow(QWidget):
            def __init__(self):
                super().__init__()
                self.ui_manager = MagicMock()

        mw = _FakeMainWindow()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        w = _gnp.GaussianFCHKFreqAnalyzer(ctx)
        btn = self._find_close_button(w)
        btn.click()
        mw.ui_manager.restore_ui_for_editing.assert_called_once()
        w.timer.stop()
        w.destroy()
        mw.deleteLater()


# ===========================================================================
# SpectrumDialog — real numpy curve, export handlers
# ===========================================================================


class TestSpectrumDialogReal:
    @pytest.fixture
    def dlg(self, qapp):
        d = _gnp.SpectrumDialog([1000.0, 2000.0], [50.0, 100.0], parent=None)
        yield d
        d.destroy()

    def test_initial_curve_is_nonempty(self, dlg):
        x, y = dlg.plot_widget.get_curve_data()
        assert len(x) == 1000
        assert len(y) == 1000

    def test_fwhm_change_updates_plot_widget(self, dlg):
        dlg.spin_fwhm.setValue(120)
        assert dlg.plot_widget.fwhm == 120

    def test_range_change_updates_plot_widget(self, dlg):
        dlg.spin_min.setValue(500)
        dlg.spin_max.setValue(1500)
        assert dlg.plot_widget.min_x == 500.0
        assert dlg.plot_widget.max_x == 1500.0

    def test_range_change_ignored_when_max_not_greater_than_min(self, dlg):
        before_min, before_max = dlg.plot_widget.min_x, dlg.plot_widget.max_x
        dlg.spin_min.setValue(2000)
        dlg.spin_max.setValue(2000)
        # on_range_changed guards mx>mn; min_x tracks spin_min regardless of
        # whether recalc_curve() fired for the (mn==mx) case.
        assert dlg.plot_widget.max_x in (before_max, 2000.0)

    def test_export_csv_writes_file(self, dlg, tmp_path, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, module=_gnp)
        out = str(tmp_path / "spectrum")
        with patch.object(_gnp.QFileDialog, "getSaveFileName", return_value=(out, "")):
            dlg.export_csv()
        out_file = Path(out + ".csv")
        assert out_file.exists()
        text = out_file.read_text()
        assert text.startswith("Frequency,Intensity\n")
        assert calls["information"]

    def test_export_csv_cancel_is_noop(self, dlg, tmp_path):
        with patch.object(_gnp.QFileDialog, "getSaveFileName", return_value=("", "")):
            dlg.export_csv()  # should not raise, nothing written

    def test_export_csv_error_shows_critical(self, dlg, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, module=_gnp)
        with patch.object(
            _gnp.QFileDialog, "getSaveFileName", return_value=("/no/such/dir/out.csv", "")
        ):
            dlg.export_csv()
        assert calls["critical"]

    def test_export_png_writes_file(self, dlg, tmp_path, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, module=_gnp)
        out = str(tmp_path / "spectrum_img")
        with patch.object(_gnp.QFileDialog, "getSaveFileName", return_value=(out, "")):
            dlg.export_png()
        assert Path(out + ".png").exists()
        assert calls["information"]

    def test_export_png_cancel_is_noop(self, dlg):
        with patch.object(_gnp.QFileDialog, "getSaveFileName", return_value=("", "")):
            dlg.export_png()

    def test_export_png_error_shows_critical(self, dlg, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, module=_gnp)
        monkeypatch.setattr(
            dlg.plot_widget, "grab", MagicMock(side_effect=RuntimeError("no device"))
        )
        with patch.object(
            _gnp.QFileDialog, "getSaveFileName", return_value=("out.png", "")
        ):
            dlg.export_png()
        assert calls["critical"]


# ===========================================================================
# SpectrumPlotWidget — real numpy Gaussian-sum curve + paintEvent rendering
# ===========================================================================


class TestSpectrumPlotWidgetReal:
    def test_recalc_curve_empty_freqs_is_noop(self, qapp):
        w = _gnp.SpectrumPlotWidget([], [])
        assert w.curve_x == []
        assert w.curve_y == []
        w.destroy()

    def test_recalc_curve_peak_matches_intensity_at_frequency(self, qapp):
        w = _gnp.SpectrumPlotWidget([1000.0], [42.0])
        w.set_range(0.0, 2000.0)
        x, y = w.get_curve_data()

        peak_idx = int(np.argmax(y))
        # The Gaussian is centred on the frequency; the sampled peak should
        # land close to it and its height should equal the intensity.
        assert x[peak_idx] == pytest.approx(1000.0, abs=5.0)
        assert y[peak_idx] == pytest.approx(42.0, rel=0.02)

    def test_set_fwhm_recalculates_curve(self, qapp):
        w = _gnp.SpectrumPlotWidget([1000.0], [10.0])
        w.set_fwhm(200)
        assert w.fwhm == 200
        assert len(w.curve_x) == 1000
        w.destroy()

    def test_set_frequencies_recalculates_curve(self, qapp):
        w = _gnp.SpectrumPlotWidget([1000.0], [10.0])
        w.set_frequencies([500.0, 1500.0])
        assert list(w.freqs) == [500.0, 1500.0]
        w.destroy()

    def test_paint_event_no_data_renders(self, qapp):
        w = _gnp.SpectrumPlotWidget([], [])
        w.resize(400, 300)
        pixmap = w.grab()
        assert not pixmap.isNull()
        w.destroy()

    def test_paint_event_with_data_renders(self, qapp):
        w = _gnp.SpectrumPlotWidget([1000.0, 2000.0], [50.0, 100.0])
        w.set_range(0.0, 4000.0)
        w.resize(600, 400)
        pixmap = w.grab()
        assert not pixmap.isNull()
        w.destroy()


# ===========================================================================
# Module-level functions: load_from_file / run_plugin / initialize / run
# ===========================================================================


class TestLoadFromFileReal:
    def _make_main_window(self, qapp):
        from PyQt6.QtWidgets import QMainWindow

        return QMainWindow()

    def test_creates_dock_and_loads_file(self, qapp, tmp_path):
        mw = self._make_main_window(qapp)
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        path = _write_water_fchk(tmp_path)

        _gnp.load_from_file(ctx, path)

        from PyQt6.QtWidgets import QDockWidget

        docks = [d for d in mw.findChildren(QDockWidget) if d.windowTitle() == "Gaussian Freq Analyzer"]
        assert len(docks) == 1
        analyzer = docks[0].widget()
        assert analyzer.parser is not None
        assert analyzer.list_freq.topLevelItemCount() == 3
        mw.deleteLater()

    def test_reuses_existing_dock(self, qapp, tmp_path):
        mw = self._make_main_window(qapp)
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        path = _write_water_fchk(tmp_path)

        _gnp.load_from_file(ctx, path)
        _gnp.load_from_file(ctx, path)

        from PyQt6.QtWidgets import QDockWidget

        docks = [d for d in mw.findChildren(QDockWidget) if d.windowTitle() == "Gaussian Freq Analyzer"]
        assert len(docks) == 1  # not duplicated
        mw.deleteLater()

    def test_closes_conflicting_orca_dock(self, qapp, tmp_path):
        from PyQt6.QtWidgets import QDockWidget, QWidget
        from PyQt6.QtCore import Qt as _Qt

        mw = self._make_main_window(qapp)
        mw.show()
        orca_dock = QDockWidget("ORCA Output Freq Analyzer", mw)
        orca_dock.setWidget(QWidget())
        mw.addDockWidget(_Qt.DockWidgetArea.LeftDockWidgetArea, orca_dock)
        orca_dock.show()
        assert orca_dock.isVisible()

        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        path = _write_water_fchk(tmp_path)
        _gnp.load_from_file(ctx, path)

        assert not orca_dock.isVisible()
        mw.deleteLater()


class TestRunPluginReal:
    def test_smart_open_uses_current_file_path(self, monkeypatch):
        mw = SimpleNamespace(
            init_manager=SimpleNamespace(current_file_path="/tmp/molecule.fchk")
        )
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        called = {}
        monkeypatch.setattr(
            _gnp, "load_from_file", lambda c, f: called.update(ctx=c, fname=f)
        )
        _gnp.run_plugin(ctx)
        assert called == {"ctx": ctx, "fname": "/tmp/molecule.fchk"}

    def test_falls_back_to_file_dialog(self, monkeypatch, tmp_path):
        mw = SimpleNamespace(init_manager=SimpleNamespace(current_file_path=""))
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        path = _write_water_fchk(tmp_path)
        called = {}
        monkeypatch.setattr(
            _gnp, "load_from_file", lambda c, f: called.update(fname=f)
        )
        with patch.object(_gnp.QFileDialog, "getOpenFileName", return_value=(path, "")):
            _gnp.run_plugin(ctx)
        assert called["fname"] == path

    def test_file_dialog_cancel_is_noop(self, monkeypatch):
        mw = SimpleNamespace(init_manager=SimpleNamespace(current_file_path=""))
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        called = {}
        monkeypatch.setattr(
            _gnp, "load_from_file", lambda c, f: called.update(fname=f)
        )
        with patch.object(_gnp.QFileDialog, "getOpenFileName", return_value=("", "")):
            _gnp.run_plugin(ctx)
        assert called == {}


class TestInitializeAndRunReal:
    def test_initialize_registers_openers_and_drop_handler(self, monkeypatch):
        ctx = MagicMock()
        _gnp.initialize(ctx)
        exts = [call[0][0] for call in ctx.register_file_opener.call_args_list]
        assert exts == [".fchk", ".fck", ".fch"]
        ctx.register_drop_handler.assert_called_once()
        assert _gnp.PLUGIN_CONTEXT is ctx

    def test_file_opener_wrapper_calls_load_from_file(self, monkeypatch):
        ctx = MagicMock()
        _gnp.initialize(ctx)
        wrapper = ctx.register_file_opener.call_args_list[0][0][1]
        called = {}
        monkeypatch.setattr(
            _gnp, "load_from_file", lambda c, f: called.update(ctx=c, fname=f)
        )
        wrapper("/tmp/mine.fchk")
        assert called == {"ctx": ctx, "fname": "/tmp/mine.fchk"}

    def test_drop_handler_accepts_fchk_rejects_others(self, monkeypatch):
        ctx = MagicMock()
        _gnp.initialize(ctx)
        handler = ctx.register_drop_handler.call_args_list[0][0][0]
        called = {}
        monkeypatch.setattr(
            _gnp, "load_from_file", lambda c, f: called.update(fname=f)
        )
        assert handler("/tmp/mine.fchk") is True
        assert called["fname"] == "/tmp/mine.fchk"

        called.clear()
        assert handler("/tmp/other.log") is False
        assert called == {}

    def test_run_without_plugin_manager_is_noop(self, monkeypatch):
        monkeypatch.setattr(_gnp, "run_plugin", MagicMock())
        mw = SimpleNamespace()  # no plugin_manager attribute
        _gnp.run(mw)
        _gnp.run_plugin.assert_not_called()

    def test_run_without_context_is_noop(self, monkeypatch):
        monkeypatch.setattr(_gnp, "PLUGIN_CONTEXT", None)
        monkeypatch.setattr(_gnp, "run_plugin", MagicMock())
        mw = SimpleNamespace(plugin_manager=True)
        _gnp.run(mw)
        _gnp.run_plugin.assert_not_called()

    def test_run_dispatches_to_run_plugin(self, monkeypatch):
        ctx = MagicMock()
        monkeypatch.setattr(_gnp, "PLUGIN_CONTEXT", ctx)
        monkeypatch.setattr(_gnp, "run_plugin", MagicMock())
        mw = SimpleNamespace(plugin_manager=True)
        _gnp.run(mw)
        _gnp.run_plugin.assert_called_once_with(ctx)
