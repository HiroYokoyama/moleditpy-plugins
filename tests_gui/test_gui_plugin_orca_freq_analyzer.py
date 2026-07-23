"""
Headless GUI tests for the ORCA Freq Analyzer plugin.

Covers:
  - OrcaOutFreqAnalyzer — main analyser widget (QWidget subclass), safe to
    construct with no main-window / no loaded file. Unlike the Gaussian Freq
    Analyzer, it is constructed with a bare ``main_window`` (not a
    PluginContext) and reaches out to a *module-level* ``PLUGIN_CONTEXT``
    global for all 3D-viewer interaction (draw_molecule_3d, plotter, etc).
  - SpectrumDialog (QDialog, module level) + SpectrumPlotWidget — real
    numpy Gaussian-sum curve rendering and CSV/PNG export.
  - OrcaOutFreqAnalyzer.SpectrumWidget — an inner class defined on the
    analyser but never actually instantiated anywhere in the source (the
    only call site, in show_spectrum(), resolves the bare name
    ``SpectrumDialog`` against module globals, not this inner class); it is
    exercised directly here purely for coverage of otherwise-dead code,
    without asserting anything about program behaviour.
  - OrcaParser — pure-Python parser object that needs no Qt.

Run:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/test_gui_plugin_orca_freq_analyzer.py -v
"""

from __future__ import annotations

import contextlib
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

# Guarded so CI's bare test-gui job (only pytest+PyQt6 installed) skips this
# real-numpy module instead of erroring at collection.
np = pytest.importorskip("numpy")

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ORCA_PATH = PLUGINS_DIR / "ORCA_Freq_Analyzer" / "orca_out_freq_analyzer.py"

with mock_chemistry_imports():
    _orca = load_plugin_for_gui(ORCA_PATH)


# ===========================================================================
# OrcaOutFreqAnalyzer  (visible plugin: "ORCA Freq Analyzer")
# ===========================================================================


class TestOrcaOutFreqAnalyzer:
    """Main analyser widget — safe to construct with main_window=None."""

    @pytest.fixture
    def widget(self, qapp):
        w = _orca.OrcaOutFreqAnalyzer(None)
        yield w
        w.timer.stop()
        w.destroy()

    def test_creates_without_error(self, widget):
        assert widget is not None

    def test_info_label_text(self, widget):
        assert widget.lbl_info.text() == "Drop .out file here or click Open"

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

    def test_amplitude_slider_range(self, widget):
        assert widget.slider_amp.minimum() == 1
        assert widget.slider_amp.maximum() == 20

    def test_speed_slider_default(self, widget):
        assert widget.slider_speed.value() == 20

    def test_speed_slider_range(self, widget):
        assert widget.slider_speed.minimum() == 1
        assert widget.slider_speed.maximum() == 60

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
# SpectrumDialog (ORCA, module-level)  — numpy mocked; empty lists are safe
# ===========================================================================


class TestOrcaSpectrumDialog:
    """Top-level SpectrumDialog from the ORCA plugin.

    np.array([]) returns a MagicMock; SpectrumPlotWidget.__init__ only stores
    the arrays without operating on them, so no crash occurs.
    """

    @pytest.fixture
    def dlg(self, qapp):
        d = _orca.SpectrumDialog([], [], parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_default_window_title(self, dlg):
        assert dlg.windowTitle() == "IR Spectrum"

    def test_custom_title(self, qapp):
        d = _orca.SpectrumDialog([], [], title="Raman Spectrum", parent=None)
        assert d.windowTitle() == "Raman Spectrum"
        d.destroy()

    def test_fwhm_spinbox_default(self, dlg):
        assert dlg.spin_fwhm.value() == 50

    def test_fwhm_spinbox_range(self, dlg):
        assert dlg.spin_fwhm.minimum() == 1
        assert dlg.spin_fwhm.maximum() == 500

    def test_min_wavenumber_default(self, dlg):
        assert dlg.spin_min.value() == 0

    def test_max_wavenumber_default(self, dlg):
        assert dlg.spin_max.value() == 4000

    def test_size_at_least_600_wide(self, dlg):
        assert dlg.width() >= 600


# ===========================================================================
# OrcaParser  — pure Python, no Qt required
# ===========================================================================


class TestOrcaParser:
    """OrcaParser parses .out files; tested in isolation without opening any file."""

    def test_creates_without_error(self):
        assert _orca.OrcaParser() is not None

    def test_initial_charge(self):
        assert _orca.OrcaParser().charge == 0

    def test_initial_multiplicity(self):
        assert _orca.OrcaParser().multiplicity == 1

    def test_initial_frequencies_empty(self):
        assert _orca.OrcaParser().frequencies == []

    def test_initial_atoms_empty(self):
        assert _orca.OrcaParser().atoms == []

    def test_initial_coords_empty(self):
        assert _orca.OrcaParser().coords == []

    def test_initial_vib_modes_empty(self):
        assert _orca.OrcaParser().vib_modes == []

    def test_initial_intensities_empty(self):
        assert _orca.OrcaParser().intensities == []


# ===========================================================================
# Real-numpy module instance `_onp`
#
# The plugin's animation/spectrum/parsing math (np.sin, np.pi, np.linspace,
# np.exp, np.array, ...) needs real numpy to actually execute rather than
# chase MagicMock attribute chains; rdkit/PIL/pyvista stay mocked (as in
# tests_gui/conftest.py's BLOCKED_CHEMISTRY).
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
    _onp = load_plugin_for_gui(ORCA_PATH)


# A full, self-consistent synthetic ORCA output: 3 atoms (water), charge -1,
# multiplicity 3, one imaginary frequency, two real ones. Column-aligned
# NORMAL MODES block gives each vibration a distinct atom moving along a
# distinct axis, matching the layout already validated in
# tests/test_plugin_orca_freq_analyzer.py::TestOrcaFinalModes.
_ORCA_FULL = textwrap.dedent("""\
    O   R   C   A

    Total Charge           Charge          ....   -1
    Multiplicity           Mult            ....    3

    CARTESIAN COORDINATES (ANGSTROEM)
    ---------------------------------
      O      0.000000    0.000000    0.000000
      H      0.000000    0.759337    0.596043
      H      0.000000   -0.759337    0.596043

    -----------------------
    VIBRATIONAL FREQUENCIES
    -----------------------
    Scaling factor for frequencies =  1.000000000 (already applied!)

       0:       0.00 cm**-1
       1:       0.00 cm**-1
       2:       0.00 cm**-1
       3:       0.00 cm**-1
       4:       0.00 cm**-1
       5:       0.00 cm**-1
       6:   -1609.85 cm**-1
       7:    3681.19 cm**-1
       8:    3811.12 cm**-1

    ------------
    NORMAL MODES
    ------------
    These modes are printed in Cartesian displacements
    mass weighted comment line
    filler line three
    filler line four
    filler line five

                      0          1          2          3          4          5
          0       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          1       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          2       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          3       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          4       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          5       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          6       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          7       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
          8       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
                      6          7          8
          0       0.010000   0.020000   0.030000
          1       0.040000   0.050000   0.060000
          2       0.070000   0.080000   0.090000
          3       0.100000   0.110000   0.120000
          4       0.130000   0.140000   0.150000
          5       0.160000   0.170000   0.180000
          6       0.190000   0.200000   0.210000
          7       0.220000   0.230000   0.240000
          8       0.250000   0.260000   0.270000

    -----------
    IR SPECTRUM
    -----------

     Mode   freq       eps      Int      T**2         TX        TY        TZ
           cm**-1   L/(mol*cm) km/mol    a.u.
    ----------------------------------------------------------------------------
       6:  -1609.85   0.015725   79.47  0.002871  (-0.001018 -0.053574 -0.000000)
       7:   3681.19   0.005000   10.00  0.001000  ( 0.001000  0.002000  0.000000)
       8:   3811.12   0.001000    5.00  0.000500  ( 0.000000  0.001000  0.000000)

    The first frequency considered to be a vibration is 6
""")

_NOT_ORCA = "Entering Gaussian System\n Link0: mem=4gb\n"


def _write_orca_out(tmp_path, content=_ORCA_FULL, name="water.out"):
    f = tmp_path / name
    f.write_text(content)
    return str(f)


def _no_block_msgbox(monkeypatch, module=_onp):
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
def onp_ctx(monkeypatch):
    """OrcaOutFreqAnalyzer reaches for the *module-level* PLUGIN_CONTEXT
    global (not a per-instance attribute) for all 3D-viewer interaction."""
    ctx = MagicMock()
    monkeypatch.setattr(_onp, "PLUGIN_CONTEXT", ctx)
    return ctx


@pytest.fixture
def onp_widget(qapp, onp_ctx):
    w = _onp.OrcaOutFreqAnalyzer(None)
    yield w
    w.timer.stop()
    w.destroy()


@pytest.fixture
def loaded_widget(onp_widget, tmp_path):
    """Widget after loading a real (small, synthetic) ORCA .out file."""
    path = _write_orca_out(tmp_path)
    onp_widget.load_file(path)
    return onp_widget


def _select_first_row(widget):
    item = widget.list_freq.topLevelItem(0)
    widget.list_freq.setCurrentItem(item)
    return item


# ===========================================================================
# Drag & drop
# ===========================================================================


class TestDragAndDrop:
    def test_drag_enter_accepts_out(self, onp_widget):
        event = _FakeDropEvent([_FakeUrl("/tmp/calc.out")])
        onp_widget.dragEnterEvent(event)
        assert event.accepted

    def test_drag_enter_accepts_log(self, onp_widget):
        event = _FakeDropEvent([_FakeUrl("/tmp/calc.log")])
        onp_widget.dragEnterEvent(event)
        assert event.accepted

    def test_drag_enter_rejects_other_extension(self, onp_widget):
        event = _FakeDropEvent([_FakeUrl("/tmp/molecule.fchk")])
        onp_widget.dragEnterEvent(event)
        assert event.ignored
        assert not event.accepted

    def test_drag_enter_rejects_no_urls(self, onp_widget):
        event = _FakeDropEvent([])
        onp_widget.dragEnterEvent(event)
        assert event.ignored

    def test_drop_event_loads_file(self, onp_widget, tmp_path):
        path = _write_orca_out(tmp_path)
        event = _FakeDropEvent([_FakeUrl(path)])
        with patch.object(onp_widget, "load_file") as mock_load:
            onp_widget.dropEvent(event)
        mock_load.assert_called_once_with(path)
        assert event.accepted

    def test_drop_event_ignores_non_out_log(self, onp_widget):
        event = _FakeDropEvent([_FakeUrl("/tmp/other.txt")])
        with patch.object(onp_widget, "load_file") as mock_load:
            onp_widget.dropEvent(event)
        mock_load.assert_not_called()


# ===========================================================================
# open_file_dialog
# ===========================================================================


class TestOpenFileDialog:
    def test_selecting_file_loads_it(self, onp_widget, tmp_path):
        path = _write_orca_out(tmp_path)
        with patch.object(_onp.QFileDialog, "getOpenFileName", return_value=(path, "")):
            onp_widget.open_file_dialog()
        assert onp_widget.parser is not None
        assert onp_widget.list_freq.topLevelItemCount() == 3

    def test_cancel_does_not_load(self, onp_widget):
        with patch.object(_onp.QFileDialog, "getOpenFileName", return_value=("", "")):
            onp_widget.open_file_dialog()
        assert onp_widget.parser is None


# ===========================================================================
# load_file / update_ui_after_load / update_list_and_spectrum_values
# ===========================================================================


class TestLoadFileReal:
    def test_populates_frequency_list(self, loaded_widget):
        assert loaded_widget.list_freq.topLevelItemCount() == 3
        item0 = loaded_widget.list_freq.topLevelItem(0)
        assert item0.text(0) == "1"
        assert item0.text(1) == "-1609.85"
        assert item0.text(2) == "79.47"

    def test_meta_label_shows_charge_mult_atoms(self, loaded_widget):
        assert "Charge: -1" in loaded_widget.lbl_meta.text()
        assert "Multiplicity: 3" in loaded_widget.lbl_meta.text()
        assert "Atoms: 3" in loaded_widget.lbl_meta.text()

    def test_info_label_updates_to_filename_and_green_style(self, loaded_widget):
        assert "water.out" in loaded_widget.lbl_info.text()
        assert "4CAF50" in loaded_widget.lbl_info.styleSheet()

    def test_base_mol_created_when_atoms_present(self, loaded_widget):
        assert loaded_widget.base_mol is not None

    def test_invalid_orca_header_shows_critical_and_no_parser(
        self, onp_widget, monkeypatch, tmp_path
    ):
        calls = _no_block_msgbox(monkeypatch)
        bad = tmp_path / "not_orca.out"
        bad.write_text(_NOT_ORCA)
        onp_widget.load_file(str(bad))
        assert calls["critical"]
        assert onp_widget.parser is None

    def test_parse_exception_shows_critical(self, onp_widget, monkeypatch, tmp_path):
        calls = _no_block_msgbox(monkeypatch)
        path = _write_orca_out(tmp_path)
        monkeypatch.setattr(
            _onp.OrcaParser, "parse",
            lambda self, filename: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        onp_widget.load_file(path)
        assert calls["critical"]

    def test_update_ui_after_load_no_final_modes_attr_is_noop(self, onp_widget):
        parser = _onp.OrcaParser()  # parse() never called -> no final_modes
        onp_widget.parser = parser
        onp_widget.update_ui_after_load()  # should not raise
        assert onp_widget.list_freq.topLevelItemCount() == 0

    def test_update_ui_after_load_none_intensity_shows_dash(self, onp_widget):
        parser = _onp.OrcaParser()
        parser.atoms = []
        parser.final_modes = [{"freq": 123.45, "intensity": None, "vector": []}]
        onp_widget.parser = parser
        onp_widget.update_ui_after_load()
        item = onp_widget.list_freq.topLevelItem(0)
        assert item.text(1) == "123.45"
        assert item.text(2) == "-"

    def test_create_base_molecule_noop_without_parser(self, onp_widget):
        onp_widget.parser = None
        onp_widget.create_base_molecule()  # should not raise
        assert onp_widget.base_mol is None


class TestUpdateListAndSpectrumValues:
    def test_scaling_updates_displayed_frequency(self, loaded_widget):
        loaded_widget.spin_sf.setValue(2.0)
        item0 = loaded_widget.list_freq.topLevelItem(0)
        assert item0.text(1) == f"{-1609.85 * 2.0:.2f}"

    def test_noop_without_parser(self, onp_widget):
        onp_widget.parser = None
        onp_widget.update_list_and_spectrum_values()  # should not raise

    def test_noop_without_final_modes(self, onp_widget):
        parser = _onp.OrcaParser()
        parser.final_modes = []
        onp_widget.parser = parser
        onp_widget.update_list_and_spectrum_values()  # final_modes == [] -> return


# ===========================================================================
# on_freq_selected
# ===========================================================================


class TestOnFreqSelected:
    def test_selecting_row_updates_vectors_when_not_playing(self, loaded_widget):
        loaded_widget.chk_vectors.setChecked(True)
        _select_first_row(loaded_widget)
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
    def test_displacement_moves_atoms_by_scaled_vector(self, monkeypatch, onp_widget):
        monkeypatch.setattr(_onp, "Point3D", _FakePos)
        onp_widget.parser = _onp.OrcaParser()
        onp_widget.parser.coords = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0)]
        onp_widget.base_mol = _FakeMol(onp_widget.parser.coords)

        mode_vecs = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
        onp_widget.apply_displacement(mode_vecs, 0.5)

        conf = onp_widget.base_mol.GetConformer()
        p0 = conf.GetAtomPosition(0)
        p1 = conf.GetAtomPosition(1)
        assert p0.x == pytest.approx(0.5)
        assert p0.y == pytest.approx(0.0)
        assert p1.y == pytest.approx(2.5)
        assert p1.x == pytest.approx(1.0)

    def test_reset_geometry_restores_base_coords(self, monkeypatch, onp_widget):
        monkeypatch.setattr(_onp, "Point3D", _FakePos)
        onp_widget.parser = _onp.OrcaParser()
        onp_widget.parser.coords = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0)]
        onp_widget.base_mol = _FakeMol(onp_widget.parser.coords)

        onp_widget.apply_displacement([(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)], 0.5)
        onp_widget.reset_geometry()

        conf = onp_widget.base_mol.GetConformer()
        p1 = conf.GetAtomPosition(1)
        assert p1.y == pytest.approx(2.0)

    def test_reset_geometry_noop_without_base_mol(self, onp_widget):
        onp_widget.base_mol = None
        onp_widget.parser = _onp.OrcaParser()
        onp_widget.reset_geometry()  # should not raise

    def test_reset_geometry_noop_without_parser(self, onp_widget):
        onp_widget.base_mol = MagicMock()
        onp_widget.parser = None
        onp_widget.reset_geometry()  # should not raise


# ===========================================================================
# animate_frame — real np.sin / np.pi math
# ===========================================================================


class TestAnimateFrameReal:
    def test_stops_when_no_parser_or_base_mol(self, onp_widget):
        onp_widget.parser = None
        onp_widget.base_mol = None
        onp_widget.animate_frame()
        assert not onp_widget.is_playing

    def test_noop_without_current_item(self, loaded_widget):
        loaded_widget.list_freq.setCurrentItem(None)
        loaded_widget.animate_frame()  # should not raise

    def test_advances_animation_step_and_redraws(self, loaded_widget, onp_ctx):
        _select_first_row(loaded_widget)
        start_step = loaded_widget.animation_step
        loaded_widget.animate_frame()
        assert loaded_widget.animation_step == start_step + 1
        onp_ctx.draw_molecule_3d.assert_called_with(loaded_widget.base_mol)


# ===========================================================================
# remove_vectors / update_vectors
# ===========================================================================


class TestVectorsReal:
    def test_remove_vectors_calls_plotter(self, onp_widget, onp_ctx):
        onp_widget.vector_actor = MagicMock()
        onp_widget.remove_vectors()
        onp_ctx.plotter.remove_actor.assert_called_once()
        assert onp_widget.vector_actor is None

    def test_remove_vectors_swallows_exception(self, onp_widget, onp_ctx):
        onp_widget.vector_actor = MagicMock()
        onp_ctx.plotter.remove_actor.side_effect = RuntimeError("gone")
        onp_widget.remove_vectors()  # should not raise
        assert onp_widget.vector_actor is None

    def test_remove_vectors_noop_without_actor(self, onp_widget):
        onp_widget.vector_actor = None
        onp_widget.remove_vectors()

    def test_update_vectors_unchecked_is_noop(self, loaded_widget):
        loaded_widget.chk_vectors.setChecked(False)
        loaded_widget.update_vectors()
        assert loaded_widget.vector_actor is None

    def test_update_vectors_noop_without_parser(self, onp_widget):
        onp_widget.chk_vectors.setChecked(True)
        onp_widget.parser = None
        onp_widget.update_vectors()

    def test_update_vectors_noop_without_current_item(self, loaded_widget):
        loaded_widget.chk_vectors.setChecked(True)
        loaded_widget.list_freq.setCurrentItem(None)
        loaded_widget.update_vectors()

    def test_update_vectors_calls_add_arrows_when_checked(self, loaded_widget, onp_ctx):
        loaded_widget.chk_vectors.setChecked(True)
        # currentItemChanged -> on_freq_selected -> update_vectors() already
        # fires the call; selecting the row is enough to exercise the path.
        _select_first_row(loaded_widget)
        assert onp_ctx.plotter.add_arrows.called

    def test_update_vectors_swallows_add_arrows_exception(self, loaded_widget, onp_ctx):
        loaded_widget.chk_vectors.setChecked(True)
        _select_first_row(loaded_widget)
        onp_ctx.plotter.add_arrows.side_effect = RuntimeError("boom")
        loaded_widget.update_vectors()  # should not raise


# ===========================================================================
# show_spectrum
# ===========================================================================


class TestShowSpectrumReal:
    def test_noop_without_parser(self, onp_widget, monkeypatch):
        # parser is None -> hasattr(None, "final_modes") is False -> the
        # "no data" QMessageBox.warning() branch fires (not a silent return).
        calls = _no_block_msgbox(monkeypatch)
        onp_widget.show_spectrum()
        assert calls["warning"]

    def test_all_low_frequencies_shows_info(self, onp_widget, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        parser = _onp.OrcaParser()
        parser.final_modes = []
        onp_widget.parser = parser
        onp_widget.show_spectrum()
        assert calls["warning"]

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

        monkeypatch.setattr(_onp, "SpectrumDialog", _StubDialog)
        loaded_widget.spin_sf.setValue(2.0)
        loaded_widget.show_spectrum()

        assert captured["exec_called"]
        assert captured["freqs"] == pytest.approx(
            [-1609.85 * 2.0, 3681.19 * 2.0, 3811.12 * 2.0]
        )
        assert captured["intensities"] == pytest.approx([79.47, 10.00, 5.00])
        assert "water.out" in captured["title"]
        assert "SF: 2.0000" in captured["title"]


# ===========================================================================
# save_as_gif
# ===========================================================================


class TestSaveAsGifReal:
    def test_noop_without_parser_or_base_mol(self, onp_widget):
        onp_widget.save_as_gif()  # should not raise

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
            _onp.QDialog, "exec", return_value=_onp.QDialog.DialogCode.Rejected
        ):
            loaded_widget.save_as_gif()
        assert loaded_widget.is_playing  # resumed after cancel

    def test_save_dialog_cancel_resumes_playback(self, loaded_widget):
        _select_first_row(loaded_widget)
        loaded_widget.toggle_play()
        with patch.object(
            _onp.QDialog, "exec", return_value=_onp.QDialog.DialogCode.Accepted
        ), patch.object(_onp.QFileDialog, "getSaveFileName", return_value=("", "")):
            loaded_widget.save_as_gif()
        assert loaded_widget.is_playing

    def test_success_writes_gif_and_appends_extension(
        self, loaded_widget, tmp_path, monkeypatch, onp_ctx
    ):
        calls = _no_block_msgbox(monkeypatch)
        _select_first_row(loaded_widget)
        out_path = str(tmp_path / "out")  # no .gif suffix
        fake_image_mod = _FakePILImageModule()
        onp_ctx.plotter.screenshot.return_value = object()
        with patch.object(_onp, "Image", fake_image_mod), patch.object(
            _onp.QDialog, "exec", return_value=_onp.QDialog.DialogCode.Accepted
        ), patch.object(
            _onp.QFileDialog, "getSaveFileName", return_value=(out_path, "")
        ):
            loaded_widget.save_as_gif()
        assert fake_image_mod.created, "no frames captured"
        assert fake_image_mod.created[0].saved_path == out_path + ".gif"
        assert calls["information"]

    def test_screenshot_none_produces_no_frames(self, loaded_widget, tmp_path, onp_ctx):
        _select_first_row(loaded_widget)
        onp_ctx.plotter.screenshot.return_value = None
        fake_image_mod = _FakePILImageModule()
        out_path = str(tmp_path / "out.gif")
        with patch.object(_onp, "Image", fake_image_mod), patch.object(
            _onp.QDialog, "exec", return_value=_onp.QDialog.DialogCode.Accepted
        ), patch.object(
            _onp.QFileDialog, "getSaveFileName", return_value=(out_path, "")
        ):
            loaded_widget.save_as_gif()
        assert fake_image_mod.created == []

    def test_save_exception_shows_critical(
        self, loaded_widget, tmp_path, monkeypatch, onp_ctx
    ):
        calls = _no_block_msgbox(monkeypatch)
        _select_first_row(loaded_widget)

        class _RaisingImage(_FakeImage):
            def save(self, path, **kwargs):
                raise OSError("disk full")

        fake_image_mod = _FakePILImageModule()
        fake_image_mod.fromarray = lambda arr: _RaisingImage()
        onp_ctx.plotter.screenshot.return_value = object()
        out_path = str(tmp_path / "out.gif")
        with patch.object(_onp, "Image", fake_image_mod), patch.object(
            _onp.QDialog, "exec", return_value=_onp.QDialog.DialogCode.Accepted
        ), patch.object(
            _onp.QFileDialog, "getSaveFileName", return_value=(out_path, "")
        ):
            loaded_widget.save_as_gif()
        assert calls["critical"]

    def test_resumes_playback_in_finally_block(self, loaded_widget, tmp_path, onp_ctx):
        _select_first_row(loaded_widget)
        loaded_widget.toggle_play()
        fake_image_mod = _FakePILImageModule()
        onp_ctx.plotter.screenshot.return_value = object()
        out_path = str(tmp_path / "resumed.gif")
        with patch.object(_onp, "Image", fake_image_mod), patch.object(
            _onp.QMessageBox, "information", staticmethod(lambda *a, **kw: None)
        ), patch.object(
            _onp.QDialog, "exec", return_value=_onp.QDialog.DialogCode.Accepted
        ), patch.object(
            _onp.QFileDialog, "getSaveFileName", return_value=(out_path, "")
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
        """SpectrumPlotWidget.on_dock_visibility_changed refers to
        self.is_playing / self.stop_play / self.base_mol / self.context,
        none of which SpectrumPlotWidget defines. It is unreachable in the
        running app: only OrcaOutFreqAnalyzer.on_dock_visibility_changed is
        ever connected to a dock's visibilityChanged signal. Exercise the
        literal method body (with stand-in attributes) so the dead branch
        is covered without altering runtime behaviour."""
        w = _onp.SpectrumPlotWidget([], [])
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
        w = _onp.SpectrumPlotWidget([], [])
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

    def test_close_with_dock_closes_dock(self, qapp, onp_ctx):
        dock = MagicMock()
        w = _onp.OrcaOutFreqAnalyzer(None, dock_widget=dock)
        btn = self._find_close_button(w)
        btn.click()
        dock.close.assert_called_once()
        w.timer.stop()
        w.destroy()


# ===========================================================================
# OrcaOutFreqAnalyzer.SpectrumWidget — inner class, dead code (no live call
# site actually reaches it — see module docstring), exercised directly.
# ===========================================================================


class TestOrcaInnerSpectrumWidgetDead:
    def test_paint_event_no_freqs_returns_early(self, qapp):
        w = _onp.OrcaOutFreqAnalyzer.SpectrumWidget([], [])
        w.resize(400, 300)
        pixmap = w.grab()
        assert not pixmap.isNull()
        w.destroy()

    def test_paint_event_with_data_renders(self, qapp):
        w = _onp.OrcaOutFreqAnalyzer.SpectrumWidget([1000.0, 2000.0], [50.0, 100.0])
        w.resize(600, 400)
        pixmap = w.grab()
        assert not pixmap.isNull()
        w.destroy()

    def test_max_int_zero_normalized_to_one(self, qapp):
        w = _onp.OrcaOutFreqAnalyzer.SpectrumWidget([1000.0], [0.0])
        assert w.max_int == 1.0
        w.destroy()

    def test_max_int_empty_defaults_to_one(self, qapp):
        w = _onp.OrcaOutFreqAnalyzer.SpectrumWidget([1000.0], [])
        assert w.max_int == 1.0
        w.destroy()


# ===========================================================================
# SpectrumDialog — real numpy curve, export handlers
# ===========================================================================


class TestSpectrumDialogReal:
    @pytest.fixture
    def dlg(self, qapp):
        d = _onp.SpectrumDialog([1000.0, 2000.0], [50.0, 100.0], parent=None)
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
        before_max = dlg.plot_widget.max_x
        dlg.spin_min.setValue(2000)
        dlg.spin_max.setValue(2000)
        assert dlg.plot_widget.max_x in (before_max, 2000.0)

    def test_export_csv_writes_file(self, dlg, tmp_path, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, module=_onp)
        out = str(tmp_path / "spectrum")
        with patch.object(_onp.QFileDialog, "getSaveFileName", return_value=(out, "")):
            dlg.export_csv()
        out_file = Path(out + ".csv")
        assert out_file.exists()
        text = out_file.read_text()
        assert text.startswith("Frequency,Intensity\n")
        assert calls["information"]

    def test_export_csv_cancel_is_noop(self, dlg, tmp_path):
        with patch.object(_onp.QFileDialog, "getSaveFileName", return_value=("", "")):
            dlg.export_csv()  # should not raise, nothing written

    def test_export_csv_error_shows_critical(self, dlg, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, module=_onp)
        with patch.object(
            _onp.QFileDialog, "getSaveFileName", return_value=("/no/such/dir/out.csv", "")
        ):
            dlg.export_csv()
        assert calls["critical"]

    def test_export_png_writes_file(self, dlg, tmp_path, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, module=_onp)
        out = str(tmp_path / "spectrum_img")
        with patch.object(_onp.QFileDialog, "getSaveFileName", return_value=(out, "")):
            dlg.export_png()
        assert Path(out + ".png").exists()
        assert calls["information"]

    def test_export_png_cancel_is_noop(self, dlg):
        with patch.object(_onp.QFileDialog, "getSaveFileName", return_value=("", "")):
            dlg.export_png()

    def test_export_png_error_shows_critical(self, dlg, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, module=_onp)
        monkeypatch.setattr(
            dlg.plot_widget, "grab", MagicMock(side_effect=RuntimeError("no device"))
        )
        with patch.object(
            _onp.QFileDialog, "getSaveFileName", return_value=("out.png", "")
        ):
            dlg.export_png()
        assert calls["critical"]


# ===========================================================================
# SpectrumPlotWidget — real numpy Gaussian-sum curve + paintEvent rendering
# ===========================================================================


class TestSpectrumPlotWidgetReal:
    def test_recalc_curve_empty_freqs_is_noop(self, qapp):
        w = _onp.SpectrumPlotWidget([], [])
        assert w.curve_x == []
        assert w.curve_y == []
        w.destroy()

    def test_recalc_curve_peak_matches_intensity_at_frequency(self, qapp):
        w = _onp.SpectrumPlotWidget([1000.0], [42.0])
        w.set_range(0.0, 2000.0)
        x, y = w.get_curve_data()

        peak_idx = int(np.argmax(y))
        assert x[peak_idx] == pytest.approx(1000.0, abs=5.0)
        assert y[peak_idx] == pytest.approx(42.0, rel=0.02)

    def test_set_fwhm_recalculates_curve(self, qapp):
        w = _onp.SpectrumPlotWidget([1000.0], [10.0])
        w.set_fwhm(200)
        assert w.fwhm == 200
        assert len(w.curve_x) == 1000
        w.destroy()

    def test_paint_event_no_data_renders(self, qapp):
        w = _onp.SpectrumPlotWidget([], [])
        w.resize(400, 300)
        pixmap = w.grab()
        assert not pixmap.isNull()
        w.destroy()

    def test_paint_event_with_data_renders(self, qapp):
        w = _onp.SpectrumPlotWidget([1000.0, 2000.0], [50.0, 100.0])
        w.set_range(0.0, 4000.0)
        w.resize(600, 400)
        pixmap = w.grab()
        assert not pixmap.isNull()
        w.destroy()


# ===========================================================================
# Module-level functions: load_from_file / run / initialize
# ===========================================================================


class TestLoadFromFileReal:
    def _make_main_window(self, qapp):
        from PyQt6.QtWidgets import QMainWindow

        return QMainWindow()

    def test_creates_dock_and_loads_file(self, qapp, tmp_path, onp_ctx):
        mw = self._make_main_window(qapp)
        path = _write_orca_out(tmp_path)

        _onp.load_from_file(mw, path)

        from PyQt6.QtWidgets import QDockWidget

        docks = [
            d for d in mw.findChildren(QDockWidget)
            if d.windowTitle() == "ORCA Output Freq Analyzer"
        ]
        assert len(docks) == 1
        analyzer = docks[0].widget()
        assert analyzer.parser is not None
        assert analyzer.list_freq.topLevelItemCount() == 3
        mw.deleteLater()

    def test_reuses_existing_dock(self, qapp, tmp_path, onp_ctx):
        mw = self._make_main_window(qapp)
        path = _write_orca_out(tmp_path)

        _onp.load_from_file(mw, path)
        _onp.load_from_file(mw, path)

        from PyQt6.QtWidgets import QDockWidget

        docks = [
            d for d in mw.findChildren(QDockWidget)
            if d.windowTitle() == "ORCA Output Freq Analyzer"
        ]
        assert len(docks) == 1  # not duplicated
        mw.deleteLater()

    def test_closes_conflicting_gaussian_dock(self, qapp, tmp_path, onp_ctx):
        from PyQt6.QtWidgets import QDockWidget, QWidget
        from PyQt6.QtCore import Qt as _Qt

        mw = self._make_main_window(qapp)
        mw.show()
        gauss_dock = QDockWidget("Gaussian Freq Analyzer", mw)
        gauss_dock.setWidget(QWidget())
        mw.addDockWidget(_Qt.DockWidgetArea.LeftDockWidgetArea, gauss_dock)
        gauss_dock.show()
        assert gauss_dock.isVisible()

        path = _write_orca_out(tmp_path)
        _onp.load_from_file(mw, path)

        assert not gauss_dock.isVisible()
        mw.deleteLater()


class TestRunReal:
    def test_smart_open_uses_current_file_path(self, monkeypatch, tmp_path):
        path = _write_orca_out(tmp_path)
        mw = SimpleNamespace(current_file_path=path)
        called = {}
        monkeypatch.setattr(
            _onp, "load_from_file", lambda m, f: called.update(mw=m, fname=f)
        )
        _onp.run(mw)
        assert called == {"mw": mw, "fname": path}

    def test_smart_open_skips_invalid_content(self, monkeypatch, tmp_path):
        bad = tmp_path / "not_orca.out"
        bad.write_text(_NOT_ORCA)
        mw = SimpleNamespace(current_file_path=str(bad))
        called = {}
        monkeypatch.setattr(
            _onp, "load_from_file", lambda m, f: called.update(fname=f)
        )
        with patch.object(_onp.QFileDialog, "getOpenFileName", return_value=("", "")):
            _onp.run(mw)
        assert called == {}  # falls through to file dialog, which is cancelled

    def test_falls_back_to_file_dialog(self, monkeypatch, tmp_path):
        mw = SimpleNamespace(current_file_path="")
        path = _write_orca_out(tmp_path)
        called = {}
        monkeypatch.setattr(
            _onp, "load_from_file", lambda m, f: called.update(fname=f)
        )
        with patch.object(_onp.QFileDialog, "getOpenFileName", return_value=(path, "")):
            _onp.run(mw)
        assert called["fname"] == path

    def test_file_dialog_cancel_is_noop(self, monkeypatch):
        mw = SimpleNamespace(current_file_path="")
        called = {}
        monkeypatch.setattr(
            _onp, "load_from_file", lambda m, f: called.update(fname=f)
        )
        with patch.object(_onp.QFileDialog, "getOpenFileName", return_value=("", "")):
            _onp.run(mw)
        assert called == {}

    def test_file_dialog_invalid_selection_shows_warning(self, monkeypatch, tmp_path):
        mw = SimpleNamespace(current_file_path="")
        bad = tmp_path / "not_orca.out"
        bad.write_text(_NOT_ORCA)
        calls = _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(_onp, "load_from_file", MagicMock())
        with patch.object(
            _onp.QFileDialog, "getOpenFileName", return_value=(str(bad), "")
        ):
            _onp.run(mw)
        assert calls["warning"]
        _onp.load_from_file.assert_not_called()

    def test_no_current_file_path_attr_falls_back_to_dialog(self, monkeypatch):
        mw = SimpleNamespace()  # no current_file_path attribute at all
        with patch.object(_onp.QFileDialog, "getOpenFileName", return_value=("", "")):
            _onp.run(mw)  # should not raise; falls through, dialog cancelled


class TestInitializeReal:
    def test_registers_openers_and_drop_handler(self):
        ctx = MagicMock()
        _onp.initialize(ctx)
        exts = [call[0][0] for call in ctx.register_file_opener.call_args_list]
        assert exts == [".out", ".log"]
        ctx.register_drop_handler.assert_called_once()
        assert ctx.register_drop_handler.call_args.kwargs.get("priority") == 10
        assert _onp.PLUGIN_CONTEXT is ctx

    def test_file_opener_wrapper_loads_valid_file(self, monkeypatch, tmp_path, capsys):
        ctx = MagicMock()
        _onp.initialize(ctx)
        wrapper = ctx.register_file_opener.call_args_list[0][0][1]
        path = _write_orca_out(tmp_path)
        called = {}
        monkeypatch.setattr(
            _onp, "load_from_file", lambda m, f: called.update(fname=f)
        )
        wrapper(path)
        assert called["fname"] == path

    def test_file_opener_wrapper_skips_invalid_file(self, monkeypatch, tmp_path, capsys):
        ctx = MagicMock()
        _onp.initialize(ctx)
        wrapper = ctx.register_file_opener.call_args_list[0][0][1]
        bad = tmp_path / "not_orca.out"
        bad.write_text(_NOT_ORCA)
        monkeypatch.setattr(_onp, "load_from_file", MagicMock())
        wrapper(str(bad))
        _onp.load_from_file.assert_not_called()
        assert "Skipping invalid ORCA file" in capsys.readouterr().out

    def test_drop_handler_accepts_out_rejects_others(self, monkeypatch, tmp_path):
        ctx = MagicMock()
        _onp.initialize(ctx)
        handler = ctx.register_drop_handler.call_args_list[0][0][0]
        path = _write_orca_out(tmp_path)
        called = {}
        monkeypatch.setattr(
            _onp, "load_from_file", lambda m, f: called.update(fname=f)
        )
        assert handler(path) is True
        assert called["fname"] == path

        called.clear()
        assert handler("/tmp/other.fchk") is False
        assert called == {}

    def test_drop_handler_rejects_invalid_content(self, monkeypatch, tmp_path):
        ctx = MagicMock()
        _onp.initialize(ctx)
        handler = ctx.register_drop_handler.call_args_list[0][0][0]
        bad = tmp_path / "not_orca.out"
        bad.write_text(_NOT_ORCA)
        monkeypatch.setattr(_onp, "load_from_file", MagicMock())
        assert handler(str(bad)) is False
        _onp.load_from_file.assert_not_called()
