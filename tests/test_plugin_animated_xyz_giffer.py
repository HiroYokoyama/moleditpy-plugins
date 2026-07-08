"""
Tests for the Animated XYZ Giffer plugin (parse_multi_frame_xyz — pure file parser).
"""

from __future__ import annotations

import ast
import logging
import textwrap
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
GIFFER_PATH = PLUGINS_DIR / "Animated_XYZ_Giffer" / "animated_xyz_giffer.py"


def _write_xyz(tmp_path, frames, filename="test.xyz"):
    """Write a multi-frame XYZ file and return its string path."""
    lines = []
    for f in frames:
        atoms = f["atoms"]  # list of (sym, x, y, z)
        comment = f.get("comment", "frame")
        lines.append(str(len(atoms)))
        lines.append(comment)
        for sym, x, y, z in atoms:
            lines.append(f"{sym} {x} {y} {z}")
    p = tmp_path / filename
    p.write_text("\n".join(lines), encoding="utf-8")
    return str(p)


def _extract_method_as_fn(path: Path, class_name: str, method_name: str, extra_globals: dict | None = None):
    """
    Use AST to extract a class method as a standalone callable.

    Needed because Qt base classes are MagicMock instances whose metaclass call
    returns a MagicMock instead of a real type, so the class definition doesn't
    produce a usable type and object.__new__ fails.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == method_name:
                    func_src = ast.get_source_segment(source, item)
                    if func_src:
                        local_ns: dict = {}
                        globs = {"logging": logging, **(extra_globals or {})}
                        exec(textwrap.dedent(func_src), globs, local_ns)
                        return local_ns[method_name]
    return None


# Extract parse_multi_frame_xyz as a standalone function (self unused in body)
_parse_xyz_raw = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "parse_multi_frame_xyz")


def _parse_xyz(file_path: str):
    return _parse_xyz_raw(None, file_path)  # self=None: method never touches self


class TestParseMultiFrameXYZ:
    def test_single_frame_two_atoms(self, tmp_path):
        p = _write_xyz(tmp_path, [
            {"atoms": [("H", 0.0, 0.0, 0.0), ("O", 0.0, 0.0, 1.0)], "comment": "water"},
        ])
        frames = _parse_xyz(p)
        assert len(frames) == 1
        assert frames[0]["symbols"] == ["H", "O"]
        assert frames[0]["comment"] == "water"
        assert len(frames[0]["coords"]) == 2

    def test_two_frames(self, tmp_path):
        p = _write_xyz(tmp_path, [
            {"atoms": [("C", 0.0, 0.0, 0.0)], "comment": "frame1"},
            {"atoms": [("C", 1.0, 0.0, 0.0)], "comment": "frame2"},
        ])
        frames = _parse_xyz(p)
        assert len(frames) == 2
        assert frames[1]["comment"] == "frame2"
        assert frames[1]["coords"][0] == (1.0, 0.0, 0.0)

    def test_empty_file_returns_empty_list(self, tmp_path):
        p = tmp_path / "empty.xyz"
        p.write_text("", encoding="utf-8")
        assert _parse_xyz(str(p)) == []

    def test_blank_lines_between_frames(self, tmp_path):
        content = "\n\n2\nframe1\nH 0 0 0\nO 0 0 1\n\n2\nframe2\nH 1 0 0\nO 1 0 1\n"
        p = tmp_path / "blanks.xyz"
        p.write_text(content, encoding="utf-8")
        frames = _parse_xyz(str(p))
        assert len(frames) == 2

    def test_incomplete_frame_at_eof_skipped(self, tmp_path):
        content = "2\nframe1\nH 0 0 0\nO 0 0 1\n5\nincomplete\nH 0 0 0\n"
        p = tmp_path / "incomplete.xyz"
        p.write_text(content, encoding="utf-8")
        frames = _parse_xyz(str(p))
        assert len(frames) == 1  # the incomplete frame is dropped

    def test_coord_values_parsed_correctly(self, tmp_path):
        p = _write_xyz(tmp_path, [
            {"atoms": [("N", 1.5, -2.3, 0.001)], "comment": "c"},
        ])
        x, y, z = _parse_xyz(p)[0]["coords"][0]
        assert abs(x - 1.5) < 1e-9
        assert abs(y - (-2.3)) < 1e-9
        assert abs(z - 0.001) < 1e-9

    def test_bad_coord_line_skipped(self, tmp_path):
        """Atom lines with a non-float coordinate are silently skipped."""
        content = "2\nbad\nH 0 0 0\nO not_a_float 0 1\n"
        p = tmp_path / "bad.xyz"
        p.write_text(content, encoding="utf-8")
        frames = _parse_xyz(str(p))
        assert len(frames) == 1
        assert frames[0]["symbols"] == ["H"]  # O line dropped

    def test_non_integer_atom_count_line_skipped(self, tmp_path):
        """Lines that can't be parsed as an int atom count are skipped."""
        content = "GARBAGE\n2\nframe1\nH 0 0 0\nO 0 0 1\n"
        p = tmp_path / "garbage.xyz"
        p.write_text(content, encoding="utf-8")
        frames = _parse_xyz(str(p))
        assert len(frames) == 1

    def test_giffer_run_does_not_raise(self):
        """run(mw) returns early when plugin_manager is absent (MagicMock lacks it by default)."""
        with mock_optional_imports():
            mod = load_plugin(GIFFER_PATH)
            mw = MagicMock(spec=[])  # no attributes → plugin_manager absent
            mod.run(mw)




# ---------------------------------------------------------------------------
# frame stepping, play/pause, update scheduling
# ---------------------------------------------------------------------------

from unittest.mock import MagicMock


class FakeCheck:
    def __init__(self, checked=False):
        self._checked = checked

    def isChecked(self):
        return self._checked

    def setChecked(self, v):
        self._checked = v


class FakeLabel:
    def __init__(self):
        self._text = ""

    def setText(self, t):
        self._text = t

    def text(self):
        return self._text


_g_next = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "next_frame")
_g_prev = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "prev_frame")
_g_toggle = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "toggle_play")
_g_fps = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "set_fps")
_g_sched = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "schedule_update")
_g_status = _extract_method_as_fn(
    GIFFER_PATH, "AnimatedXYZPlayer", "update_status_silent"
)


def _giffer_self(n_frames=5, current=0, loop=True, playing=False, fps=10):
    fake = SimpleNamespace()
    fake.frames = [{"coords": []} for _ in range(n_frames)]
    fake.current_frame_idx = current
    fake.target_frame_idx = current
    fake.chk_loop = FakeCheck(loop)
    fake.is_playing = playing
    fake.fps = fps
    fake.schedule_update = MagicMock()
    fake.toggle_play = MagicMock()
    fake.timer = MagicMock()
    fake.btn_play = MagicMock()
    fake.context = MagicMock()
    fake.base_mol = "BASE"
    fake.last_display_mol = None
    fake.is_updating_view = False
    fake.pending_update = False
    fake.do_effective_update = MagicMock()
    fake.lbl_status = FakeLabel()
    fake.slider = MagicMock()
    return fake


class TestGifferFrameStepping:
    def test_next_advances(self):
        fake = _giffer_self(current=1)
        _g_next(fake)
        assert fake.target_frame_idx == 2
        fake.schedule_update.assert_called_once()

    def test_next_wraps_with_loop(self):
        fake = _giffer_self(current=4, loop=True)
        _g_next(fake)
        assert fake.target_frame_idx == 0

    def test_next_at_end_no_loop_stops_playback(self):
        fake = _giffer_self(current=4, loop=False, playing=True)
        _g_next(fake)
        fake.toggle_play.assert_called_once()
        fake.schedule_update.assert_not_called()

    def test_next_at_end_no_loop_not_playing_noop(self):
        fake = _giffer_self(current=4, loop=False, playing=False)
        _g_next(fake)
        fake.toggle_play.assert_not_called()
        fake.schedule_update.assert_not_called()

    def test_prev_steps_back(self):
        fake = _giffer_self(current=3)
        _g_prev(fake)
        assert fake.target_frame_idx == 2

    def test_prev_wraps_with_loop(self):
        fake = _giffer_self(current=0, loop=True)
        _g_prev(fake)
        assert fake.target_frame_idx == 4

    def test_prev_at_start_no_loop_noop(self):
        fake = _giffer_self(current=0, loop=False)
        _g_prev(fake)
        fake.schedule_update.assert_not_called()


class TestGifferPlayback:
    def test_play_starts_timer_with_fps_interval(self):
        fake = _giffer_self(playing=False, fps=20)
        _g_toggle(fake)
        assert fake.is_playing is True
        fake.timer.start.assert_called_once_with(50)
        fake.btn_play.setText.assert_called_with("Pause")

    def test_play_from_last_frame_restarts(self):
        fake = _giffer_self(current=4, playing=False)
        _g_toggle(fake)
        assert fake.target_frame_idx == 0
        fake.schedule_update.assert_called_once()

    def test_pause_stops_timer_and_pushes_molecule(self):
        fake = _giffer_self(playing=True)
        fake.last_display_mol = "FRAME_MOL"
        _g_toggle(fake)
        assert fake.is_playing is False
        fake.timer.stop.assert_called_once()
        assert fake.context.current_molecule == "FRAME_MOL"

    def test_pause_falls_back_to_base_mol(self):
        fake = _giffer_self(playing=True)
        fake.last_display_mol = None
        _g_toggle(fake)
        assert fake.context.current_molecule == "BASE"

    def test_set_fps_restarts_timer_only_when_playing(self):
        fake = _giffer_self(playing=True)
        _g_fps(fake, 25)
        assert fake.fps == 25
        fake.timer.start.assert_called_once_with(40)

        fake2 = _giffer_self(playing=False)
        _g_fps(fake2, 25)
        fake2.timer.start.assert_not_called()


class TestGifferScheduling:
    def test_schedule_runs_when_idle(self):
        fake = _giffer_self()
        _g_sched(fake)
        fake.do_effective_update.assert_called_once()
        assert fake.is_updating_view is True  # cleared by do_effective_update IRL

    def test_schedule_defers_when_busy(self):
        fake = _giffer_self()
        fake.is_updating_view = True
        _g_sched(fake)
        fake.do_effective_update.assert_not_called()
        assert fake.pending_update is True

    def test_status_label_and_slider_sync(self):
        fake = _giffer_self(n_frames=8, current=2)
        _g_status(fake)
        assert fake.lbl_status.text() == "Frame: 3 / 8"
        fake.slider.blockSignals.assert_any_call(True)
        fake.slider.setValue.assert_called_once_with(2)
        fake.slider.blockSignals.assert_any_call(False)
