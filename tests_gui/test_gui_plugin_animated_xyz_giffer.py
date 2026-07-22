"""
Headless GUI tests for the Animated XYZ Giffer plugin.

Covers: AnimatedXYZPlayer (QDialog) + parse_multi_frame_xyz.

All tests use real PyQt6 (QT_QPA_PLATFORM=offscreen).
Chemistry libraries are mocked via mock_chemistry_imports().
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

GIFFER_PATH = PLUGINS_DIR / "Animated_XYZ_Giffer" / "animated_xyz_giffer.py"

with mock_chemistry_imports():
    _giffer = load_plugin_for_gui(GIFFER_PATH)


# ===========================================================================
# AnimatedXYZPlayer  (Animated XYZ Giffer — real QDialog)
# ===========================================================================
# try_import_from_mainwindow() accesses self.mw.init_manager which fails if
# mw is None; patch it to a no-op so the dialog constructs cleanly.


def _giffer_ctx() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None  # keeps QDialog parent = None
    return ctx


class TestAnimatedXYZPlayer:
    """AnimatedXYZPlayer dialog with try_import_from_mainwindow patched."""

    @pytest.fixture
    def player(self, qapp):
        ctx = _giffer_ctx()
        with patch.object(
            _giffer.AnimatedXYZPlayer,
            "try_import_from_mainwindow",
            lambda self: None,
        ):
            w = _giffer.AnimatedXYZPlayer(context=ctx)
        yield w
        w.timer.stop()
        w.destroy()

    def test_creates_without_error(self, player):
        assert player is not None

    def test_window_title(self, player):
        assert player.windowTitle() == "Animated XYZ Player"

    def test_slider_initially_disabled(self, player):
        assert not player.slider.isEnabled()

    def test_play_button_initially_disabled(self, player):
        assert not player.btn_play.isEnabled()

    def test_prev_button_initially_disabled(self, player):
        assert not player.btn_prev.isEnabled()

    def test_next_button_initially_disabled(self, player):
        assert not player.btn_next.isEnabled()

    def test_fps_spin_default(self, player):
        assert player.spin_fps.value() == 10

    def test_fps_spin_range(self, player):
        assert player.spin_fps.minimum() == 1
        assert player.spin_fps.maximum() == 60

    def test_dynamic_bonds_checked_by_default(self, player):
        assert player.chk_dynamic_bonds.isChecked()

    def test_loop_checked_by_default(self, player):
        assert player.chk_loop.isChecked()

    def test_save_gif_button_disabled_initially(self, player):
        assert not player.btn_save_gif.isEnabled()

    def test_status_label_initial_text(self, player):
        assert "0" in player.lbl_status.text()

    def test_frames_initially_empty(self, player):
        assert player.frames == []

    def test_is_playing_initially_false(self, player):
        assert player.is_playing is False

    def test_timer_not_active_initially(self, player):
        assert not player.timer.isActive()


# ===========================================================================
# parse_multi_frame_xyz  (method on AnimatedXYZPlayer — pure file parser)
# ===========================================================================


@pytest.fixture(scope="module")
def _player_instance(qapp):
    """One AnimatedXYZPlayer used only to call parse_multi_frame_xyz."""
    ctx = _giffer_ctx()
    with patch.object(
        _giffer.AnimatedXYZPlayer,
        "try_import_from_mainwindow",
        lambda self: None,
    ):
        w = _giffer.AnimatedXYZPlayer(context=ctx)
    yield w
    w.timer.stop()
    w.destroy()


_SINGLE_FRAME_XYZ = """\
3
water molecule
O  0.000  0.000  0.117
H  0.757  0.000 -0.468
H -0.757  0.000 -0.468
"""

_TWO_FRAME_XYZ = """\
2
frame 1
C  0.0  0.0  0.0
H  1.0  0.0  0.0
2
frame 2
C  0.0  0.0  0.1
H  1.0  0.0  0.1
"""

_EMPTY_XYZ = ""


class TestAnimatedXYZParser:
    """parse_multi_frame_xyz() — pure text parsing, no chemistry or Qt needed."""

    def test_empty_file_returns_empty_list(self, tmp_path, _player_instance):
        p = tmp_path / "empty.xyz"
        p.write_text(_EMPTY_XYZ)
        assert _player_instance.parse_multi_frame_xyz(str(p)) == []

    def test_single_frame_count(self, tmp_path, _player_instance):
        p = tmp_path / "single.xyz"
        p.write_text(_SINGLE_FRAME_XYZ)
        frames = _player_instance.parse_multi_frame_xyz(str(p))
        assert len(frames) == 1

    def test_single_frame_symbols(self, tmp_path, _player_instance):
        p = tmp_path / "single.xyz"
        p.write_text(_SINGLE_FRAME_XYZ)
        frames = _player_instance.parse_multi_frame_xyz(str(p))
        assert frames[0]["symbols"] == ["O", "H", "H"]

    def test_single_frame_coords(self, tmp_path, _player_instance):
        p = tmp_path / "single.xyz"
        p.write_text(_SINGLE_FRAME_XYZ)
        frames = _player_instance.parse_multi_frame_xyz(str(p))
        ox, oy, oz = frames[0]["coords"][0]
        assert abs(ox) < 1e-6
        assert abs(oz - 0.117) < 1e-6

    def test_single_frame_comment(self, tmp_path, _player_instance):
        p = tmp_path / "single.xyz"
        p.write_text(_SINGLE_FRAME_XYZ)
        frames = _player_instance.parse_multi_frame_xyz(str(p))
        assert frames[0]["comment"] == "water molecule"

    def test_two_frames_count(self, tmp_path, _player_instance):
        p = tmp_path / "two.xyz"
        p.write_text(_TWO_FRAME_XYZ)
        frames = _player_instance.parse_multi_frame_xyz(str(p))
        assert len(frames) == 2

    def test_two_frames_second_coord_z(self, tmp_path, _player_instance):
        p = tmp_path / "two.xyz"
        p.write_text(_TWO_FRAME_XYZ)
        frames = _player_instance.parse_multi_frame_xyz(str(p))
        _, _, z = frames[1]["coords"][0]
        assert abs(z - 0.1) < 1e-6

    def test_incomplete_frame_skipped(self, tmp_path, _player_instance):
        bad = "3\nincomplete\nC 0 0 0\n"  # only 1 atom line of 3
        p = tmp_path / "bad.xyz"
        p.write_text(bad)
        frames = _player_instance.parse_multi_frame_xyz(str(p))
        assert frames == []

    def test_non_integer_atom_count_skipped(self, tmp_path, _player_instance):
        content = "NOT_A_NUMBER\n2\nframe1\nH 0 0 0\nO 0 0 1\n"
        p = tmp_path / "garbage.xyz"
        p.write_text(content)
        frames = _player_instance.parse_multi_frame_xyz(str(p))
        assert len(frames) == 1

    def test_bad_coord_line_skipped(self, tmp_path, _player_instance):
        content = "2\nbad\nH 0 0 0\nO not_a_float 0 1\n"
        p = tmp_path / "bad_coord.xyz"
        p.write_text(content)
        frames = _player_instance.parse_multi_frame_xyz(str(p))
        assert frames[0]["symbols"] == ["H"]


# ===========================================================================
# Fake chemistry primitives — real RDKit is mocked out in this test module,
# so we drive AnimatedXYZPlayer's chemistry-touching methods with small,
# deterministic stand-ins for Chem/rdGeometry/rdDetermineBonds/PIL.Image.
# This lets us exercise real Qt widget instances (buttons, sliders, dialogs)
# together with realistic (not MagicMock-random) molecule bookkeeping.
# ===========================================================================


class FakeBond:
    def __init__(self, a, b, order="SINGLE"):
        self._a = a
        self._b = b
        self._order = order

    def GetBeginAtomIdx(self):
        return self._a

    def GetEndAtomIdx(self):
        return self._b

    def GetBondType(self):
        return self._order


class FakeConformer:
    def __init__(self, n):
        self.positions = [None] * n

    def SetAtomPosition(self, idx, pt):
        self.positions[idx] = pt


class FakeAtom:
    def __init__(self, sym):
        self.sym = sym


class FakeMol:
    """Minimal stand-in for an RDKit RWMol/Mol as used by the giffer."""

    def __init__(self, n_atoms=0, bonds=None):
        self._n_atoms = n_atoms
        self.bonds = list(bonds or [])
        self._conformer = None

    def AddAtom(self, atom):
        self._n_atoms += 1
        return self._n_atoms - 1

    def AddConformer(self, conf):
        self._conformer = conf

    def GetConformer(self):
        if self._conformer is None:
            self._conformer = FakeConformer(self._n_atoms)
        return self._conformer

    def GetNumAtoms(self):
        return self._n_atoms

    def GetNumBonds(self):
        return len(self.bonds)

    def GetBonds(self):
        return list(self.bonds)

    def GetBondWithIdx(self, i):
        return self.bonds[i]

    def RemoveBond(self, a, b):
        self.bonds = [
            bd for bd in self.bonds
            if not (bd.GetBeginAtomIdx() == a and bd.GetEndAtomIdx() == b)
        ]

    def AddBond(self, a, b, order="SINGLE"):
        self.bonds.append(FakeBond(a, b, order))

    def GetMol(self):
        return self

    def copy(self):
        return FakeMol(
            self._n_atoms,
            [FakeBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondType()) for b in self.bonds],
        )


class FakeChem(SimpleNamespace):
    """Stand-in for the ``rdkit.Chem`` module used by the giffer."""

    def __init__(self, xyz_block_mol_factory=None):
        super().__init__()
        self._xyz_factory = xyz_block_mol_factory or (lambda block: FakeMol(0))
        self.RWMol = lambda: FakeMol()
        self.Atom = FakeAtom
        self.Conformer = lambda n: FakeConformer(n)
        self.Mol = lambda m: m.copy()
        self.MolFromXYZBlock = lambda block: self._xyz_factory(block)


class FakeImage:
    def __init__(self):
        self.mode = "RGB"
        self.saved_path = None
        self.saved_kwargs = None

    def convert(self, mode):
        self.mode = mode
        return self

    def split(self):
        return [self, self, self, self]

    def quantize(self, colors=255):
        return self

    def paste(self, idx, mask):
        return None

    def save(self, path, **kwargs):
        self.saved_path = path
        self.saved_kwargs = kwargs


class FakePILImageModule(SimpleNamespace):
    def __init__(self):
        super().__init__()
        self.created = []
        self.fromarray = self._fromarray
        self.eval = lambda img, fn: FakeImage()

    def _fromarray(self, arr):
        img = FakeImage()
        self.created.append(img)
        return img


def _make_multiframe_xyz(tmp_path, n_frames=3, n_atoms=2, name="traj.xyz"):
    lines = []
    for i in range(n_frames):
        lines.append(str(n_atoms))
        lines.append(f"frame {i}")
        for a in range(n_atoms):
            lines.append(f"C {float(a)} {float(i)} 0.0")
    p = tmp_path / name
    p.write_text("\n".join(lines), encoding="utf-8")
    return str(p)


@pytest.fixture
def loaded_player(qapp, tmp_path):
    """A real AnimatedXYZPlayer with a 3-frame trajectory loaded, using
    deterministic fake chemistry primitives instead of MagicMock chemistry."""
    ctx = _giffer_ctx()
    ctx.plotter = MagicMock()
    fake_chem = FakeChem()
    patches = [
        patch.object(_giffer, "Chem", fake_chem),
        patch.object(_giffer, "rdGeometry", SimpleNamespace(Point3D=lambda x, y, z: (x, y, z))),
        patch.object(_giffer, "rdDetermineBonds", SimpleNamespace(
            DetermineConnectivity=lambda m: None,
            DetermineBondOrders=lambda m: None,
        )),
        patch.object(_giffer, "HAS_PIL", True),
        patch.object(_giffer.AnimatedXYZPlayer, "try_import_from_mainwindow", lambda self: None),
    ]
    for p in patches:
        p.start()
    w = None
    try:
        w = _giffer.AnimatedXYZPlayer(context=ctx)
        w.mw = MagicMock()
        path = _make_multiframe_xyz(tmp_path)
        w.load_from_path(path)
        yield w
    finally:
        if w is not None:
            w.timer.stop()
            w._reload_timer.stop()
            w.destroy()
        for p in patches:
            p.stop()


class TestGifferLoadedPlayerBasics:
    def test_frames_loaded(self, loaded_player):
        assert len(loaded_player.frames) == 3

    def test_base_mol_created(self, loaded_player):
        assert loaded_player.base_mol is not None
        assert loaded_player.base_mol.GetNumAtoms() == 2

    def test_slider_range_matches_frame_count(self, loaded_player):
        assert loaded_player.slider.maximum() == 2

    def test_controls_enabled_after_load(self, loaded_player):
        assert loaded_player.slider.isEnabled()
        assert loaded_player.btn_prev.isEnabled()
        assert loaded_player.btn_play.isEnabled()
        assert loaded_player.btn_next.isEnabled()
        assert loaded_player.btn_save_gif.isEnabled()

    def test_status_label_reflects_first_frame(self, loaded_player):
        assert loaded_player.lbl_status.text() == "Frame: 1 / 3"

    def test_comment_label_set(self, loaded_player):
        assert loaded_player.lbl_comment.text() == "frame 0"


class TestGifferFrameNavigationReal:
    def test_next_frame_advances_and_updates_ui(self, loaded_player):
        w = loaded_player
        w.next_frame()
        assert w.current_frame_idx == 1
        assert w.lbl_status.text() == "Frame: 2 / 3"

    def test_prev_frame_from_start_wraps_with_loop(self, loaded_player):
        w = loaded_player
        assert w.chk_loop.isChecked()
        w.prev_frame()
        assert w.current_frame_idx == 2

    def test_next_frame_wraps_at_end_with_loop(self, loaded_player):
        w = loaded_player
        w.slider.setValue(2)
        w.next_frame()
        assert w.current_frame_idx == 0

    def test_next_frame_stops_at_end_without_loop(self, loaded_player):
        w = loaded_player
        w.chk_loop.setChecked(False)
        w.slider.setValue(2)
        w.next_frame()
        assert w.current_frame_idx == 2  # unchanged, no wrap

    def test_slider_change_schedules_update(self, loaded_player):
        w = loaded_player
        w.slider.setValue(1)
        assert w.current_frame_idx == 1
        assert w.lbl_status.text() == "Frame: 2 / 3"


class TestGifferPlaybackReal:
    def test_toggle_play_starts_and_stops_timer(self, loaded_player):
        w = loaded_player
        w.toggle_play()
        assert w.is_playing is True
        assert w.timer.isActive()
        assert w.btn_play.text() == "Pause"

        w.toggle_play()
        assert w.is_playing is False
        assert not w.timer.isActive()
        assert w.btn_play.text() == "Play"

    def test_on_timer_advances_frame(self, loaded_player):
        w = loaded_player
        w.on_timer()
        assert w.current_frame_idx == 1

    def test_set_fps_updates_value(self, loaded_player):
        w = loaded_player
        w.spin_fps.setValue(30)
        assert w.fps == 30

    def test_set_fps_restarts_timer_while_playing(self, loaded_player):
        w = loaded_player
        w.toggle_play()
        w.set_fps(25)
        assert w.fps == 25
        assert w.timer.isActive()

    def test_play_from_last_frame_restarts_at_zero(self, loaded_player):
        w = loaded_player
        w.slider.setValue(2)
        assert w.current_frame_idx == 2
        w.toggle_play()
        assert w.current_frame_idx == 0

    def test_pause_pushes_last_display_mol_to_context(self, loaded_player):
        w = loaded_player
        w.toggle_play()
        w.toggle_play()
        assert w.context.current_molecule is w.last_display_mol


class TestGifferDynamicBondsReal:
    def test_dynamic_bonds_default_reconstructs_topology(self, loaded_player):
        """Default (checked) path: MolFromXYZBlock returns a fresh mol per frame."""
        w = loaded_player
        assert w.chk_dynamic_bonds.isChecked()
        w.next_frame()
        assert w.last_display_mol is not None

    def test_dynamic_bonds_off_restores_original_topology(self, loaded_player):
        w = loaded_player
        w.chk_dynamic_bonds.setChecked(False)
        # Simulate a stale extra bond only on base_mol (not on original_topology)
        w.base_mol.AddBond(0, 1, "SINGLE")
        assert w.base_mol.GetNumBonds() == 1
        w.next_frame()
        assert w.base_mol.GetNumBonds() == 0  # restored from original_topology

    def test_dynamic_bonds_xyz_block_none_falls_back_to_coords(self, loaded_player, tmp_path):
        ctx = _giffer_ctx()
        ctx.plotter = MagicMock()
        fake_chem = FakeChem(xyz_block_mol_factory=lambda block: None)
        with patch.object(_giffer, "Chem", fake_chem), \
             patch.object(_giffer, "rdGeometry", SimpleNamespace(Point3D=lambda x, y, z: (x, y, z))), \
             patch.object(_giffer, "rdDetermineBonds", SimpleNamespace(
                 DetermineConnectivity=lambda m: None,
                 DetermineBondOrders=lambda m: None,
             )), \
             patch.object(_giffer.AnimatedXYZPlayer, "try_import_from_mainwindow", lambda self: None):
            w = _giffer.AnimatedXYZPlayer(context=ctx)
            w.mw = MagicMock()
            path = _make_multiframe_xyz(tmp_path, name="fallback.xyz")
            w.load_from_path(path)
            w.next_frame()
            assert w.last_display_mol is w.base_mol
            w.timer.stop()
            w._reload_timer.stop()
            w.destroy()

    def test_dynamic_bonds_exception_falls_back_to_base_mol(self, loaded_player, tmp_path):
        ctx = _giffer_ctx()
        ctx.plotter = MagicMock()
        fake_chem = FakeChem()

        def _raise(m):
            raise RuntimeError("boom")

        with patch.object(_giffer, "Chem", fake_chem), \
             patch.object(_giffer, "rdGeometry", SimpleNamespace(Point3D=lambda x, y, z: (x, y, z))), \
             patch.object(_giffer, "rdDetermineBonds", SimpleNamespace(
                 DetermineConnectivity=_raise,
                 DetermineBondOrders=lambda m: None,
             )), \
             patch.object(_giffer.AnimatedXYZPlayer, "try_import_from_mainwindow", lambda self: None):
            w = _giffer.AnimatedXYZPlayer(context=ctx)
            w.mw = MagicMock()
            path = _make_multiframe_xyz(tmp_path, name="exc.xyz")
            w.load_from_path(path)
            w.next_frame()
            assert w.last_display_mol is w.base_mol
            w.timer.stop()
            w._reload_timer.stop()
            w.destroy()


class TestGifferTryImportFromMainwindowReal:
    def test_real_import_path_loads_xyz(self, qapp, tmp_path):
        """Exercise the actual try_import_from_mainwindow body (unpatched)."""
        ctx = _giffer_ctx()
        ctx.plotter = MagicMock()
        path = _make_multiframe_xyz(tmp_path, name="auto.xyz")
        mw = MagicMock()
        mw.init_manager.current_file_path = path
        fake_chem = FakeChem()
        with patch.object(_giffer, "Chem", fake_chem), \
             patch.object(_giffer, "rdGeometry", SimpleNamespace(Point3D=lambda x, y, z: (x, y, z))), \
             patch.object(_giffer, "rdDetermineBonds", SimpleNamespace(
                 DetermineConnectivity=lambda m: None,
                 DetermineBondOrders=lambda m: None,
             )):
            ctx.get_main_window.return_value = None  # QDialog parent stays valid (None)
            with patch.object(
                _giffer.AnimatedXYZPlayer, "try_import_from_mainwindow", lambda self: None
            ):
                w = _giffer.AnimatedXYZPlayer(context=ctx)
            # try_import_from_mainwindow was skipped during __init__ (self.mw was
            # still None at that point); set self.mw now and invoke the real,
            # unpatched body directly to exercise it deterministically.
            w.mw = mw
            w.try_import_from_mainwindow()
            assert len(w.frames) == 3
            w.timer.stop()
            w._reload_timer.stop()
            w.destroy()

    def test_non_xyz_path_is_skipped(self, qapp):
        ctx = _giffer_ctx()
        mw = MagicMock()
        mw.init_manager.current_file_path = "molecule.pdb"
        with patch.object(_giffer.AnimatedXYZPlayer, "try_import_from_mainwindow", lambda self: None):
            w = _giffer.AnimatedXYZPlayer(context=ctx)
        w.mw = mw
        w.try_import_from_mainwindow()
        assert w.frames == []
        w.timer.stop()
        w._reload_timer.stop()
        w.destroy()


class TestGifferLoadFileReal:
    def test_load_file_dialog_selection_loads_trajectory(self, loaded_player, tmp_path):
        w = loaded_player
        new_path = _make_multiframe_xyz(tmp_path, n_frames=5, name="picked.xyz")
        with patch.object(_giffer.QFileDialog, "getOpenFileName", return_value=(new_path, "")):
            w.load_file()
        assert len(w.frames) == 5

    def test_load_file_dialog_cancel_noop(self, loaded_player):
        w = loaded_player
        with patch.object(_giffer.QFileDialog, "getOpenFileName", return_value=("", "")):
            w.load_file()
        assert len(w.frames) == 3  # unchanged


class TestGifferLoadFromPathErrors:
    def test_empty_frames_shows_warning(self, qapp, tmp_path):
        ctx = _giffer_ctx()
        p = tmp_path / "empty.xyz"
        p.write_text("", encoding="utf-8")
        with patch.object(_giffer.AnimatedXYZPlayer, "try_import_from_mainwindow", lambda self: None), \
             patch.object(_giffer, "QMessageBox", MagicMock()) as qmb:
            w = _giffer.AnimatedXYZPlayer(context=ctx)
            w.load_from_path(str(p))
            qmb.warning.assert_called_once()
        assert w.base_mol is None
        w.timer.stop()
        w._reload_timer.stop()
        w.destroy()

    def test_parse_exception_shows_critical(self, qapp):
        ctx = _giffer_ctx()
        with patch.object(_giffer.AnimatedXYZPlayer, "try_import_from_mainwindow", lambda self: None), \
             patch.object(_giffer, "QMessageBox", MagicMock()) as qmb:
            w = _giffer.AnimatedXYZPlayer(context=ctx)
            w.load_from_path("/does/not/exist.xyz")
            qmb.critical.assert_called_once()
        w.timer.stop()
        w._reload_timer.stop()
        w.destroy()


class TestGifferSaveAsGifReal:
    def test_no_frames_is_noop(self, qapp):
        ctx = _giffer_ctx()
        with patch.object(_giffer.AnimatedXYZPlayer, "try_import_from_mainwindow", lambda self: None):
            w = _giffer.AnimatedXYZPlayer(context=ctx)
        with patch.object(_giffer, "QDialog") as qd:
            w.save_as_gif()
            qd.assert_not_called()
        w.timer.stop()
        w._reload_timer.stop()
        w.destroy()

    def test_dialog_cancel_resumes_playback(self, loaded_player):
        w = loaded_player
        w.toggle_play()
        assert w.is_playing is True
        with patch.object(_giffer.QDialog, "exec", return_value=_giffer.QDialog.DialogCode.Rejected):
            w.save_as_gif()
        assert w.is_playing is True  # resumed after cancel

    def test_file_dialog_cancel_resumes_playback(self, loaded_player):
        w = loaded_player
        w.toggle_play()
        with patch.object(_giffer.QDialog, "exec", return_value=_giffer.QDialog.DialogCode.Accepted), \
             patch.object(_giffer.QFileDialog, "getSaveFileName", return_value=("", "")):
            w.save_as_gif()
        assert w.is_playing is True

    def test_success_writes_gif_and_appends_extension(self, loaded_player, tmp_path):
        w = loaded_player
        out_path = str(tmp_path / "out")  # no .gif suffix
        fake_image_mod = FakePILImageModule()
        w.context.plotter.screenshot.return_value = object()
        with patch.object(_giffer, "Image", fake_image_mod), \
             patch.object(_giffer, "QMessageBox", MagicMock()) as qmb, \
             patch.object(_giffer.QFileDialog, "getSaveFileName", return_value=(out_path, "")), \
             patch.object(_giffer.QDialog, "exec", return_value=_giffer.QDialog.DialogCode.Accepted):
            w.save_as_gif()
        assert fake_image_mod.created, "no frames captured"
        assert fake_image_mod.created[0].saved_path == out_path + ".gif"
        qmb.information.assert_called_once()

    def test_screenshot_none_shows_warning(self, loaded_player):
        w = loaded_player
        w.context.plotter.screenshot.return_value = None
        with patch.object(_giffer, "QMessageBox", MagicMock()) as qmb, \
             patch.object(_giffer.QFileDialog, "getSaveFileName", return_value=("out.gif", "")), \
             patch.object(_giffer.QDialog, "exec", return_value=_giffer.QDialog.DialogCode.Accepted):
            w.save_as_gif()
        qmb.warning.assert_called_once()

    def test_save_exception_shows_critical(self, loaded_player):
        w = loaded_player

        class RaisingImage(FakeImage):
            def save(self, path, **kwargs):
                raise OSError("disk full")

        fake_image_mod = FakePILImageModule()
        fake_image_mod.fromarray = lambda arr: RaisingImage()
        w.context.plotter.screenshot.return_value = object()
        with patch.object(_giffer, "Image", fake_image_mod), \
             patch.object(_giffer, "QMessageBox", MagicMock()) as qmb, \
             patch.object(_giffer.QFileDialog, "getSaveFileName", return_value=("out.gif", "")), \
             patch.object(_giffer.QDialog, "exec", return_value=_giffer.QDialog.DialogCode.Accepted):
            w.save_as_gif()
        qmb.critical.assert_called_once()

    def test_pauses_running_playback_then_resumes(self, loaded_player):
        w = loaded_player
        w.toggle_play()
        fake_image_mod = FakePILImageModule()
        w.context.plotter.screenshot.return_value = object()
        with patch.object(_giffer, "Image", fake_image_mod), \
             patch.object(_giffer, "QMessageBox", MagicMock()), \
             patch.object(_giffer.QFileDialog, "getSaveFileName", return_value=("resumed.gif", "")), \
             patch.object(_giffer.QDialog, "exec", return_value=_giffer.QDialog.DialogCode.Accepted):
            w.save_as_gif()
        assert w.is_playing is True  # resumed in the finally block


class TestGifferDocumentResetReal:
    def test_clears_state(self, loaded_player):
        w = loaded_player
        w.on_document_reset()
        assert w.frames == []
        assert w.base_mol is None
        assert w.lbl_file.text() == "No file loaded"
        assert not w.slider.isEnabled()
        assert not w.btn_play.isEnabled()

    def test_stops_playback(self, loaded_player):
        w = loaded_player
        w.toggle_play()
        assert w.is_playing is True
        w.on_document_reset()
        assert w.is_playing is False
        assert not w.timer.isActive()

    def test_reload_poll_loads_new_path_when_visible(self, loaded_player, tmp_path):
        w = loaded_player
        w.on_document_reset()
        new_path = _make_multiframe_xyz(tmp_path, n_frames=4, name="reloaded.xyz")
        w.mw.init_manager.current_file_path = new_path
        with patch.object(w, "isVisible", return_value=True):
            w._on_reload_poll()
        assert len(w.frames) == 4
        assert not w._reload_timer.isActive()

    def test_reload_poll_keeps_polling_without_path(self, loaded_player):
        w = loaded_player
        w.on_document_reset()
        w.mw.init_manager.current_file_path = None
        with patch.object(w, "isVisible", return_value=True):
            w._on_reload_poll()
        assert w._reload_attempts == 1

    def test_reload_poll_stops_when_hidden(self, loaded_player):
        w = loaded_player
        w.on_document_reset()
        with patch.object(w, "isVisible", return_value=False):
            w._on_reload_poll()
        assert not w._reload_timer.isActive()

    def test_reload_poll_gives_up_after_max_attempts(self, loaded_player):
        w = loaded_player
        w.on_document_reset()
        w._reload_attempts = 14
        w.mw.init_manager.current_file_path = None
        with patch.object(w, "isVisible", return_value=True):
            w._on_reload_poll()
        assert not w._reload_timer.isActive()


class TestGifferCloseEventReal:
    def test_close_pushes_last_display_mol(self, loaded_player):
        from PyQt6.QtGui import QCloseEvent

        w = loaded_player
        w.next_frame()
        assert w.last_display_mol is not None
        w.closeEvent(QCloseEvent())
        assert w.context.current_molecule is w.last_display_mol
        w.context.push_undo_checkpoint.assert_called_once()

    def test_close_does_not_call_draw_molecule_3d(self, loaded_player):
        """Regression: draw_molecule_3d must NOT be called on close (v1.3.2 fix)."""
        from PyQt6.QtGui import QCloseEvent

        w = loaded_player
        w.closeEvent(QCloseEvent())
        w.context.draw_molecule_3d.assert_not_called()
