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


# ---------------------------------------------------------------------------
# on_slider_changed / on_timer
# ---------------------------------------------------------------------------

_g_slider_changed = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "on_slider_changed")
_g_on_timer = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "on_timer")


class TestGifferSliderAndTimer:
    def test_slider_changed_sets_target_and_schedules(self):
        fake = _giffer_self()
        _g_slider_changed(fake, 3)
        assert fake.target_frame_idx == 3
        fake.schedule_update.assert_called_once()

    def test_on_timer_calls_next_frame(self):
        fake = SimpleNamespace()
        fake.next_frame = MagicMock()
        _g_on_timer(fake)
        fake.next_frame.assert_called_once()


# ---------------------------------------------------------------------------
# try_import_from_mainwindow / load_from_path / load_file
# ---------------------------------------------------------------------------

_g_try_import = _extract_method_as_fn(
    GIFFER_PATH, "AnimatedXYZPlayer", "try_import_from_mainwindow"
)
_g_load_from_path = _extract_method_as_fn(
    GIFFER_PATH,
    "AnimatedXYZPlayer",
    "load_from_path",
    extra_globals={
        "os": __import__("os"),
        "QMessageBox": MagicMock(),
        "HAS_PIL": True,
    },
)
_g_load_file = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "load_file")


class TestGifferTryImportFromMainwindow:
    def test_xyz_path_triggers_load(self):
        mw = MagicMock()
        mw.init_manager.current_file_path = "traj.xyz"
        fake = SimpleNamespace(mw=mw)
        fake.load_from_path = MagicMock()
        _g_try_import(fake)
        fake.load_from_path.assert_called_once_with("traj.xyz")

    def test_extxyz_path_triggers_load(self):
        mw = MagicMock()
        mw.init_manager.current_file_path = "traj.extxyz"
        fake = SimpleNamespace(mw=mw)
        fake.load_from_path = MagicMock()
        _g_try_import(fake)
        fake.load_from_path.assert_called_once_with("traj.extxyz")

    def test_non_xyz_path_skipped(self):
        mw = MagicMock()
        mw.init_manager.current_file_path = "molecule.pdb"
        fake = SimpleNamespace(mw=mw)
        fake.load_from_path = MagicMock()
        _g_try_import(fake)
        fake.load_from_path.assert_not_called()

    def test_no_current_file_path_skipped(self):
        mw = MagicMock()
        mw.init_manager.current_file_path = ""
        fake = SimpleNamespace(mw=mw)
        fake.load_from_path = MagicMock()
        _g_try_import(fake)
        fake.load_from_path.assert_not_called()


def _load_path_self(parse_result):
    fake = SimpleNamespace()
    fake.parse_multi_frame_xyz = MagicMock(return_value=parse_result)
    fake.lbl_file = FakeLabel()
    fake.setWindowTitle = MagicMock()
    fake._reload_timer = MagicMock()
    fake.slider = MagicMock()
    fake.btn_prev = MagicMock()
    fake.btn_play = MagicMock()
    fake.btn_next = MagicMock()
    fake.btn_save_gif = MagicMock()
    fake.create_base_molecule = MagicMock()
    fake.update_view = MagicMock()
    fake.update_status = MagicMock()
    return fake


class TestGifferLoadFromPath:
    def test_valid_frames_populate_state(self):
        frames = [{"symbols": ["H"], "coords": [(0, 0, 0)], "comment": "f1"}]
        fake = _load_path_self(frames)
        _g_load_from_path(fake, "/tmp/traj.xyz")
        assert fake.frames == frames
        assert fake.lbl_file.text() == "traj.xyz"
        fake.setWindowTitle.assert_called_once_with(
            "Animated XYZ Player - traj.xyz"
        )
        assert fake.current_frame_idx == 0
        assert fake.target_frame_idx == 0
        fake.slider.setRange.assert_called_once_with(0, 0)
        fake.slider.setEnabled.assert_called_once_with(True)
        fake.create_base_molecule.assert_called_once()
        fake.update_view.assert_called_once()
        fake.update_status.assert_called_once()

    def test_empty_frames_warns_and_stops(self):
        fake = _load_path_self([])
        _g_load_from_path(fake, "/tmp/empty.xyz")
        fake.create_base_molecule.assert_not_called()

    def test_parse_exception_reports_critical(self):
        fake = _load_path_self(None)
        fake.parse_multi_frame_xyz.side_effect = RuntimeError("bad file")
        _g_load_from_path(fake, "/tmp/broken.xyz")
        fake.create_base_molecule.assert_not_called()


class TestGifferLoadFile:
    def test_file_selected_delegates_to_load_from_path(self):
        qfd = MagicMock()
        qfd.getOpenFileName.return_value = ("/tmp/chosen.xyz", "")
        fn = _extract_method_as_fn(
            GIFFER_PATH, "AnimatedXYZPlayer", "load_file",
            extra_globals={"QFileDialog": qfd},
        )
        fake = SimpleNamespace()
        fake.load_from_path = MagicMock()
        fn(fake)
        fake.load_from_path.assert_called_once_with("/tmp/chosen.xyz")

    def test_cancel_does_not_load(self):
        qfd = MagicMock()
        qfd.getOpenFileName.return_value = ("", "")
        fn = _extract_method_as_fn(
            GIFFER_PATH, "AnimatedXYZPlayer", "load_file",
            extra_globals={"QFileDialog": qfd},
        )
        fake = SimpleNamespace()
        fake.load_from_path = MagicMock()
        fn(fake)
        fake.load_from_path.assert_not_called()


# ---------------------------------------------------------------------------
# create_base_molecule
# ---------------------------------------------------------------------------

def _base_mol_self(frames, mw_has_estimator=True):
    fake = SimpleNamespace()
    fake.frames = frames
    fake.mw = MagicMock()
    if not mw_has_estimator:
        fake.mw.io_manager = MagicMock(spec=[])
    fake.context = MagicMock()
    return fake


class TestGifferCreateBaseMolecule:
    def test_no_frames_is_noop(self):
        chem = MagicMock()
        fn = _extract_method_as_fn(
            GIFFER_PATH, "AnimatedXYZPlayer", "create_base_molecule",
            extra_globals={"Chem": chem, "rdGeometry": MagicMock()},
        )
        fake = _base_mol_self([])
        fn(fake)
        chem.RWMol.assert_not_called()

    def test_builds_mol_and_sets_context(self):
        chem = MagicMock()
        rdgeom = MagicMock()
        fn = _extract_method_as_fn(
            GIFFER_PATH, "AnimatedXYZPlayer", "create_base_molecule",
            extra_globals={"Chem": chem, "rdGeometry": rdgeom},
        )
        frames = [{"symbols": ["H", "O"], "coords": [(0, 0, 0), (0, 0, 1)]}]
        fake = _base_mol_self(frames)
        fn(fake)
        chem.RWMol.return_value.AddAtom.assert_called()
        fake.context.enter_3d_mode.assert_called_once()
        fake.context.reset_3d_camera.assert_called_once()
        assert fake.context.current_molecule is fake.base_mol

    def test_estimate_bonds_called_when_available(self):
        chem = MagicMock()
        fn = _extract_method_as_fn(
            GIFFER_PATH, "AnimatedXYZPlayer", "create_base_molecule",
            extra_globals={"Chem": chem, "rdGeometry": MagicMock()},
        )
        frames = [{"symbols": ["H"], "coords": [(0, 0, 0)]}]
        fake = _base_mol_self(frames, mw_has_estimator=True)
        fn(fake)
        fake.mw.io_manager.estimate_bonds_from_distances.assert_called_once()

    def test_estimate_bonds_skipped_when_unavailable(self):
        chem = MagicMock()
        fn = _extract_method_as_fn(
            GIFFER_PATH, "AnimatedXYZPlayer", "create_base_molecule",
            extra_globals={"Chem": chem, "rdGeometry": MagicMock()},
        )
        frames = [{"symbols": ["H"], "coords": [(0, 0, 0)]}]
        fake = _base_mol_self(frames, mw_has_estimator=False)
        fn(fake)  # must not raise even though estimate_bonds_from_distances is absent


# ---------------------------------------------------------------------------
# closeEvent
# ---------------------------------------------------------------------------

def _fake_super_factory():
    return SimpleNamespace(closeEvent=lambda event: None)


_g_close_event = _extract_method_as_fn(
    GIFFER_PATH, "AnimatedXYZPlayer", "closeEvent",
    extra_globals={"super": _fake_super_factory},
)


class TestGifferCloseEvent:
    def test_stops_timer(self):
        fake = SimpleNamespace()
        fake.timer = MagicMock()
        fake._reload_timer = MagicMock()
        fake.context = MagicMock()
        fake.last_display_mol = None
        fake.base_mol = None
        _g_close_event(fake, MagicMock())
        fake.timer.stop.assert_called_once()

    def test_pushes_last_display_mol_when_present(self):
        fake = SimpleNamespace()
        fake.timer = MagicMock()
        fake._reload_timer = MagicMock()
        fake.context = MagicMock()
        fake.last_display_mol = "LAST"
        fake.base_mol = "BASE"
        _g_close_event(fake, MagicMock())
        assert fake.context.current_molecule == "LAST"
        fake.context.push_undo_checkpoint.assert_called_once()

    def test_falls_back_to_base_mol(self):
        fake = SimpleNamespace()
        fake.timer = MagicMock()
        fake._reload_timer = MagicMock()
        fake.context = MagicMock()
        fake.last_display_mol = None
        fake.base_mol = "BASE"
        _g_close_event(fake, MagicMock())
        assert fake.context.current_molecule == "BASE"

    def test_exception_in_push_is_swallowed(self):
        fake = SimpleNamespace()
        fake.timer = MagicMock()
        fake._reload_timer = MagicMock()
        fake.context = MagicMock()
        fake.context.push_undo_checkpoint.side_effect = RuntimeError("boom")
        fake.last_display_mol = None
        fake.base_mol = None
        _g_close_event(fake, MagicMock())  # must not raise


# ---------------------------------------------------------------------------
# on_document_reset / initialize
# ---------------------------------------------------------------------------

_g_doc_reset = _extract_method_as_fn(
    GIFFER_PATH, "AnimatedXYZPlayer", "on_document_reset"
)
_g_reload_poll = _extract_method_as_fn(
    GIFFER_PATH, "AnimatedXYZPlayer", "_on_reload_poll"
)


def _doc_reset_self():
    fake = SimpleNamespace()
    fake.timer = MagicMock()
    fake.is_playing = True
    fake.frames = [{"coords": []}]
    fake.current_frame_idx = 3
    fake.target_frame_idx = 3
    fake.base_mol = "BASE"
    fake.original_topology = "TOPO"
    fake.last_display_mol = "LAST"
    fake.lbl_file = FakeLabel()
    fake.lbl_comment = FakeLabel()
    fake.lbl_status = FakeLabel()
    fake.setWindowTitle = MagicMock()
    fake.slider = MagicMock()
    fake.btn_prev = MagicMock()
    fake.btn_play = MagicMock()
    fake.btn_next = MagicMock()
    fake.btn_save_gif = MagicMock()
    fake.toggle_play = MagicMock()
    fake._reload_timer = MagicMock()
    fake._reload_attempts = 7
    return fake


class TestGifferDocumentReset:
    def test_clears_state_and_ui(self):
        fake = _doc_reset_self()
        _g_doc_reset(fake)
        assert fake.frames == []
        assert fake.current_frame_idx == 0
        assert fake.target_frame_idx == 0
        assert fake.base_mol is None
        assert fake.original_topology is None
        assert fake.last_display_mol is None
        assert fake.lbl_file.text() == "No file loaded"
        assert fake.lbl_comment.text() == ""
        assert fake.lbl_status.text() == "Frame: 0 / 0"
        fake.setWindowTitle.assert_called_once_with("Animated XYZ Player")
        fake.slider.setEnabled.assert_called_once_with(False)
        for btn in (fake.btn_prev, fake.btn_play, fake.btn_next, fake.btn_save_gif):
            btn.setEnabled.assert_called_once_with(False)

    def test_stops_playback_without_toggle_play(self):
        """toggle_play() would push the stale molecule into the cleared document."""
        fake = _doc_reset_self()
        _g_doc_reset(fake)
        fake.timer.stop.assert_called_once()
        assert fake.is_playing is False
        fake.btn_play.setText.assert_called_once_with("Play")
        fake.toggle_play.assert_not_called()

    def test_starts_reload_poll(self):
        """Re-import is deferred: the host sets the new path after handlers ran."""
        fake = _doc_reset_self()
        _g_doc_reset(fake)
        assert fake._reload_attempts == 0
        fake._reload_timer.start.assert_called_once_with(200)


class TestGifferReloadPoll:
    @staticmethod
    def _poll_self(visible=True, path=None, attempts=0):
        fake = SimpleNamespace()
        fake.isVisible = MagicMock(return_value=visible)
        fake.mw = MagicMock()
        fake.mw.init_manager.current_file_path = path
        fake._reload_attempts = attempts
        fake._reload_timer = MagicMock()
        fake.load_from_path = MagicMock()
        return fake

    def test_loads_new_xyz_and_stops(self):
        fake = self._poll_self(path="C:/data/new_traj.xyz")
        _g_reload_poll(fake)
        fake._reload_timer.stop.assert_called_once()
        fake.load_from_path.assert_called_once_with("C:/data/new_traj.xyz")

    def test_extxyz_also_loads(self):
        fake = self._poll_self(path="run.extxyz")
        _g_reload_poll(fake)
        fake.load_from_path.assert_called_once_with("run.extxyz")

    def test_keeps_polling_while_no_path(self):
        fake = self._poll_self(path=None, attempts=3)
        _g_reload_poll(fake)
        fake._reload_timer.stop.assert_not_called()
        fake.load_from_path.assert_not_called()
        assert fake._reload_attempts == 4

    def test_non_xyz_path_never_loaded(self):
        fake = self._poll_self(path="molecule.pdb", attempts=3)
        _g_reload_poll(fake)
        fake.load_from_path.assert_not_called()

    def test_gives_up_after_max_attempts(self):
        fake = self._poll_self(path=None, attempts=14)
        _g_reload_poll(fake)
        fake._reload_timer.stop.assert_called_once()
        fake.load_from_path.assert_not_called()

    def test_stops_when_window_hidden(self):
        fake = self._poll_self(visible=False, path="new.xyz")
        _g_reload_poll(fake)
        fake._reload_timer.stop.assert_called_once()
        fake.load_from_path.assert_not_called()


class TestGifferInitialize:
    def test_registers_document_reset_handler(self):
        with mock_optional_imports():
            mod = load_plugin(GIFFER_PATH)
            ctx = MagicMock()
            mod.initialize(ctx)
            ctx.register_document_reset_handler.assert_called_once()
            ctx.add_menu_action.assert_not_called()  # run() auto-registers the menu

    def test_handler_dispatches_to_visible_window(self):
        with mock_optional_imports():
            mod = load_plugin(GIFFER_PATH)
            ctx = MagicMock()
            mod.initialize(ctx)
            handler = ctx.register_document_reset_handler.call_args[0][0]

            win = MagicMock()
            win.isVisible.return_value = True
            ctx.get_window.return_value = win
            handler()
            win.on_document_reset.assert_called_once()
            ctx.get_window.assert_called_with("main_panel")

    def test_handler_skips_hidden_or_missing_window(self):
        with mock_optional_imports():
            mod = load_plugin(GIFFER_PATH)
            ctx = MagicMock()
            mod.initialize(ctx)
            handler = ctx.register_document_reset_handler.call_args[0][0]

            ctx.get_window.return_value = None
            handler()  # must not raise

            win = MagicMock()
            win.isVisible.return_value = False
            ctx.get_window.return_value = win
            handler()
            win.on_document_reset.assert_not_called()


# ---------------------------------------------------------------------------
# run_plugin / run entry points
# ---------------------------------------------------------------------------

class TestGifferRunPlugin:
    def test_run_plugin_shows_and_registers_window(self):
        with mock_optional_imports():
            mod = load_plugin(GIFFER_PATH)
            ctx = MagicMock()
            ctx.get_main_window.return_value = MagicMock()
            mod.run_plugin(ctx)
            ctx.register_window.assert_called_once()
            assert ctx.register_window.call_args[0][0] == "main_panel"

    def test_run_with_plugin_manager_invokes_run_plugin(self):
        with mock_optional_imports():
            mod = load_plugin(GIFFER_PATH)
            called = []
            mod.run_plugin = lambda ctx: called.append(ctx)
            mw = MagicMock()
            mod.run(mw)
            assert len(called) == 1
