"""
Headless GUI tests for export / rendering plugin dialogs.

Covers:
  - High Resolution Imager  (no dialog class — module + initialize tests)
  - Blender Export           (no dialog class — module + initialize tests)
  - POV-Ray Export           (no dialog class — module + initialize tests)
  - Animated XYZ Giffer     (AnimatedXYZPlayer QDialog + parse_multi_frame_xyz)

All tests use real PyQt6 (QT_QPA_PLATFORM=offscreen).
Chemistry libraries are mocked via mock_chemistry_imports().
"""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

HI_RES_PATH = PLUGINS_DIR / "High_Resolution_Imager" / "high_res_imager.py"
BLENDER_PATH = PLUGINS_DIR / "Blender_Export" / "blender_export.py"
POVRAY_PATH = PLUGINS_DIR / "POV-Ray_Export" / "povray_export.py"
GIFFER_PATH = PLUGINS_DIR / "Animated_XYZ_Giffer" / "animated_xyz_giffer.py"

# ---------------------------------------------------------------------------
# Load all plugins once at collection time — chemistry mocked, Qt real.
# ---------------------------------------------------------------------------

with mock_chemistry_imports():
    _hi_res = load_plugin_for_gui(HI_RES_PATH)
    _blender = load_plugin_for_gui(BLENDER_PATH)
    _povray = load_plugin_for_gui(POVRAY_PATH)
    _giffer = load_plugin_for_gui(GIFFER_PATH)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ctx() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# High Resolution Imager — module-level and initialize() tests
# (dialog is built inline in take_screenshot(); not a separate class)
# ===========================================================================


class TestHighResImagerModule:
    """Module constants and initialize() for High Resolution Imager."""

    def test_plugin_name(self):
        assert _hi_res.PLUGIN_NAME == "High Resolution Imager"

    def test_plugin_version_is_date_string(self):
        assert isinstance(_hi_res.PLUGIN_VERSION, str)
        assert len(_hi_res.PLUGIN_VERSION) == 10  # YYYY.MM.DD

    def test_plugin_context_initially_none(self):
        assert _hi_res.PLUGIN_CONTEXT is None

    def test_initialize_registers_export_action(self):
        ctx = _make_ctx()
        _hi_res.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_screenshot(self):
        ctx = _make_ctx()
        _hi_res.initialize(ctx)
        label = ctx.add_export_action.call_args.args[0]
        assert "Screenshot" in label

    def test_initialize_callback_is_callable(self):
        ctx = _make_ctx()
        _hi_res.initialize(ctx)
        cb = ctx.add_export_action.call_args.args[1]
        assert callable(cb)

    def test_initialize_stores_context(self):
        ctx = _make_ctx()
        _hi_res.PLUGIN_CONTEXT = None
        _hi_res.initialize(ctx)
        assert _hi_res.PLUGIN_CONTEXT is ctx
        _hi_res.PLUGIN_CONTEXT = None  # reset


# ===========================================================================
# Blender Export — module-level and initialize()/run() tests
# (no dialog class; export logic is callable-only)
# ===========================================================================


class TestBlenderExportModule:
    """Module constants, initialize(), and run() for Blender Export."""

    def test_plugin_name(self):
        assert _blender.PLUGIN_NAME == "Blender Export"

    def test_plugin_version_is_string(self):
        assert isinstance(_blender.PLUGIN_VERSION, str)

    def test_plugin_context_initially_none(self):
        assert _blender.PLUGIN_CONTEXT is None

    def test_initialize_registers_export_action(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_blender(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        label = ctx.add_export_action.call_args.args[0]
        assert "Blender" in label

    def test_initialize_callback_is_callable(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        cb = ctx.add_export_action.call_args.args[1]
        assert callable(cb)

    def test_run_noop_without_plugin_manager(self):
        mw = MagicMock(spec=[])  # no plugin_manager attribute
        _blender.run(mw)  # must not raise

    def test_run_noop_when_context_is_none(self):
        _blender.PLUGIN_CONTEXT = None
        mw = MagicMock()  # has plugin_manager (MagicMock auto-attribute)
        _blender.run(mw)  # PLUGIN_CONTEXT is None → early return, no error


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
