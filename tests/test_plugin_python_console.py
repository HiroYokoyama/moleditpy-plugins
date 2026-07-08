"""
Tests for the Python Console plugin: ConsoleInput.append_history, GUIHelp.__repr__.

Qt-inheriting classes become MagicMock when loaded under mock_optional_imports
(their metaclass is the mocked PyQt6 type). append_history is extracted from
the plugin source via AST and compiled as a standalone function so we can
test pure logic that does not touch self's Qt attributes.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from conftest import extract_function, load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CONSOLE_PATH = PLUGINS_DIR / "Python_Console" / "console.py"

with mock_optional_imports():
    _console = load_plugin(CONSOLE_PATH)

# ConsoleInput.append_history — only touches self.history / self.history_index
_append_history = extract_function(CONSOLE_PATH, "ConsoleInput", "append_history")


class _HistoryStub:
    """Minimal stand-in for ConsoleInput — only carries history state."""

    def __init__(self):
        self.history: list = []
        self.history_index: int = 0


class TestAppendHistory:
    def test_adds_first_entry(self):
        stub = _HistoryStub()
        _append_history(stub, "print('hello')")
        assert stub.history == ["print('hello')"]
        assert stub.history_index == 1

    def test_ignores_empty_string(self):
        stub = _HistoryStub()
        _append_history(stub, "")
        assert stub.history == []
        assert stub.history_index == 0

    def test_no_consecutive_duplicate(self):
        stub = _HistoryStub()
        _append_history(stub, "x = 1")
        _append_history(stub, "x = 1")
        assert stub.history == ["x = 1"]
        assert stub.history_index == 1

    def test_non_consecutive_duplicate_is_allowed(self):
        stub = _HistoryStub()
        for cmd in ["x = 1", "y = 2", "x = 1"]:
            _append_history(stub, cmd)
        assert stub.history == ["x = 1", "y = 2", "x = 1"]

    def test_index_always_points_past_end(self):
        stub = _HistoryStub()
        for cmd in ["a", "b", "c"]:
            _append_history(stub, cmd)
        assert stub.history_index == len(stub.history)

    def test_multiple_distinct_entries(self):
        stub = _HistoryStub()
        cmds = ["import math", "math.sqrt(4)", "print(2)"]
        for cmd in cmds:
            _append_history(stub, cmd)
        assert stub.history == cmds
        assert stub.history_index == 3


class TestGUIHelp:
    def test_repr_mentions_help_object(self):
        h = _console.GUIHelp()
        assert "help(object)" in repr(h)

    def test_repr_warns_about_freezing(self):
        r = repr(_console.GUIHelp())
        lowered = r.lower()
        assert "disabled" in lowered or "prevent" in lowered

    def test_call_no_args_prints_repr(self, capsys):
        h = _console.GUIHelp()
        h()
        out = capsys.readouterr().out
        assert "help(object)" in out

    def test_call_with_object_delegates_to_pydoc_help(self):
        h = _console.GUIHelp()
        with patch.object(_console.pydoc, "help") as mock_help:
            h(str)
        mock_help.assert_called_once_with(str)

    def test_call_uses_plainpager_and_restores_original(self):
        h = _console.GUIHelp()
        sentinel_pager = object()
        _console.pydoc.pager = sentinel_pager
        try:
            with patch.object(_console.pydoc, "help"):
                h(str)
            assert _console.pydoc.pager is sentinel_pager
        finally:
            _console.pydoc.pager = sentinel_pager

    def test_call_swallows_pydoc_exceptions(self, capsys):
        h = _console.GUIHelp()
        with patch.object(_console.pydoc, "help", side_effect=RuntimeError("boom")):
            h(str)  # must not raise
        err = capsys.readouterr().err
        assert "help() error" in err
        assert "boom" in err


class TestInitializeAndRun:
    def test_initialize_sets_plugin_context(self):
        ctx = MagicMock()
        _console.initialize(ctx)
        assert _console.PLUGIN_CONTEXT is ctx

    def test_run_returns_early_without_context(self):
        _console.PLUGIN_CONTEXT = None
        mw = MagicMock()
        _console.run(mw)  # must not raise

    def test_run_unwraps_host_attribute(self):
        ctx = MagicMock()
        win = MagicMock()
        ctx.get_window.return_value = win
        _console.initialize(ctx)
        outer_mw = MagicMock()
        real_mw = MagicMock()
        outer_mw.host = real_mw
        _console.run(outer_mw)
        win.show.assert_called_once()
        win.raise_.assert_called_once()
        win.activateWindow.assert_called_once()

    def test_run_creates_dialog_when_no_existing_window(self):
        ctx = MagicMock()
        ctx.get_window.return_value = None
        _console.initialize(ctx)
        with patch.object(_console, "PythonConsoleDialog") as dialog_cls:
            dialog_cls.return_value = MagicMock()
            mw = MagicMock()
            del mw.host
            _console.run(mw)
        dialog_cls.assert_called_once_with(ctx)


class _SuperProxy:
    """Stand-in for the QPlainTextEdit base reached via super() in keyPressEvent."""

    def __init__(self):
        self.calls = []

    def keyPressEvent(self, event):
        self.calls.append(event)


def _key_press_fn(super_proxy):
    class _FakeKey:
        Key_Return = "Return"
        Key_Up = "Up"
        Key_Down = "Down"
        Key_A = "A"

    class _FakeModifier:
        ShiftModifier = 1
        NoModifier = 0

    class _FakeMoveOp:
        End = "End"

    fake_qt = SimpleNamespace(
        Key=_FakeKey, KeyboardModifier=_FakeModifier
    )
    globs = {
        "Qt": fake_qt,
        "QTextCursor": SimpleNamespace(MoveOperation=_FakeMoveOp),
        "super": lambda *a, **kw: super_proxy,
    }
    return extract_function(CONSOLE_PATH, "ConsoleInput", "keyPressEvent", globs)


class _ConsoleInputStub:
    """Minimal stand-in for ConsoleInput carrying only what keyPressEvent touches."""

    def __init__(self, history=None, history_index=0, block_number=0, block_count=1):
        self.history = history or []
        self.history_index = history_index
        self._block_number = block_number
        self._block_count = block_count
        self.set_texts = []
        self.moved = []
        self.cleared = False

    def textCursor(self):
        return SimpleNamespace(blockNumber=lambda: self._block_number)

    def blockCount(self):
        return self._block_count

    def setPlainText(self, text):
        self.set_texts.append(text)

    def moveCursor(self, op):
        self.moved.append(op)

    def clear(self):
        self.cleared = True


def _event(key, shift=False):
    return SimpleNamespace(
        key=lambda: key,
        modifiers=lambda: (1 if shift else 0),
    )


class TestConsoleInputKeyPressEvent:
    def test_enter_executes_and_does_not_insert_newline(self):
        proxy = _SuperProxy()
        fn = _key_press_fn(proxy)
        stub = _ConsoleInputStub()
        stub.execute_signal = MagicMock()
        fn(stub, _event("Return", shift=False))
        stub.execute_signal.emit.assert_called_once()
        assert proxy.calls == []  # newline not inserted

    def test_shift_enter_inserts_newline_via_super(self):
        proxy = _SuperProxy()
        fn = _key_press_fn(proxy)
        stub = _ConsoleInputStub()
        stub.execute_signal = MagicMock()
        ev = _event("Return", shift=True)
        fn(stub, ev)
        stub.execute_signal.emit.assert_not_called()
        assert proxy.calls == [ev]

    def test_up_at_top_line_navigates_history(self):
        proxy = _SuperProxy()
        fn = _key_press_fn(proxy)
        stub = _ConsoleInputStub(history=["a", "b"], history_index=2, block_number=0)
        fn(stub, _event("Up"))
        assert stub.history_index == 1
        assert stub.set_texts == ["b"]
        assert stub.moved == ["End"]

    def test_up_not_at_top_line_moves_cursor_normally(self):
        proxy = _SuperProxy()
        fn = _key_press_fn(proxy)
        stub = _ConsoleInputStub(history=["a"], history_index=1, block_number=1)
        ev = _event("Up")
        fn(stub, ev)
        assert proxy.calls == [ev]
        assert stub.set_texts == []

    def test_up_with_empty_history_moves_cursor_normally(self):
        proxy = _SuperProxy()
        fn = _key_press_fn(proxy)
        stub = _ConsoleInputStub(history=[], history_index=0, block_number=0)
        ev = _event("Up")
        fn(stub, ev)
        assert proxy.calls == [ev]

    def test_down_at_last_line_navigates_forward_in_history(self):
        proxy = _SuperProxy()
        fn = _key_press_fn(proxy)
        stub = _ConsoleInputStub(
            history=["a", "b"], history_index=0, block_number=0, block_count=1
        )
        fn(stub, _event("Down"))
        assert stub.history_index == 1
        assert stub.set_texts == ["b"]

    def test_down_past_end_of_history_clears_input(self):
        proxy = _SuperProxy()
        fn = _key_press_fn(proxy)
        stub = _ConsoleInputStub(
            history=["a", "b"], history_index=1, block_number=0, block_count=1
        )
        fn(stub, _event("Down"))
        assert stub.history_index == 2
        assert stub.cleared is True

    def test_down_not_at_last_line_moves_cursor_normally(self):
        proxy = _SuperProxy()
        fn = _key_press_fn(proxy)
        stub = _ConsoleInputStub(
            history=["a", "b"], history_index=1, block_number=0, block_count=2
        )
        ev = _event("Down")
        fn(stub, ev)
        assert proxy.calls == [ev]

    def test_other_key_falls_through_to_super(self):
        proxy = _SuperProxy()
        fn = _key_press_fn(proxy)
        stub = _ConsoleInputStub()
        ev = _event("A")
        fn(stub, ev)
        assert proxy.calls == [ev]


def _console_self(command, mol="MOL"):
    import code as code_mod

    s = SimpleNamespace()
    s.input_area = MagicMock()
    s.input_area.toPlainText.return_value = command
    s.outputs = []  # (text, color)
    s.append_output = lambda text, color=None: s.outputs.append((text, color))
    s.output_area = MagicMock()
    s.context = MagicMock()
    s._get_best_mol = lambda: mol
    s.local_scope = {}
    s.interpreter = code_mod.InteractiveInterpreter(s.local_scope)
    return s


def _run_code_fn():
    import io
    import traceback
    from contextlib import redirect_stderr, redirect_stdout

    globs = {
        "io": io,
        "traceback": traceback,
        "redirect_stdout": redirect_stdout,
        "redirect_stderr": redirect_stderr,
    }
    return extract_function(CONSOLE_PATH, "PythonConsoleDialog", "run_code", globs)


class TestConsoleRunCode:
    def test_empty_command_is_ignored(self):
        fn = _run_code_fn()
        s = _console_self("   \n  ")
        fn(s)
        s.input_area.append_history.assert_not_called()
        assert s.outputs == []

    def test_print_output_captured(self):
        fn = _run_code_fn()
        s = _console_self("print(21 * 2)")
        fn(s)
        s.output_area.append.assert_any_call("42")

    def test_expression_result_echoed_in_single_mode(self):
        fn = _run_code_fn()
        s = _console_self("1 + 1")
        fn(s)
        s.output_area.append.assert_any_call("2")

    # code.InteractiveInterpreter routes the traceback through sys.excepthook,
    # which pytest-qt's exception capture would report as a Qt-loop error.
    @pytest.mark.qt_no_exception_capture
    def test_exception_written_in_error_color(self):
        fn = _run_code_fn()
        s = _console_self("1/0")
        fn(s)
        err = [t for t, c in s.outputs if c == "#FF5252"]
        assert err and "ZeroDivisionError" in err[0]

    def test_incomplete_block_warns(self):
        fn = _run_code_fn()
        s = _console_self("def f():")
        fn(s)
        assert any("Incomplete" in t for t, _ in s.outputs)

    def test_multiline_runs_in_exec_mode(self):
        fn = _run_code_fn()
        s = _console_self("a = 6\nprint(a * 7)")
        fn(s)
        s.output_area.append.assert_any_call("42")

    def test_command_stored_in_history(self):
        fn = _run_code_fn()
        s = _console_self("x = 1")
        fn(s)
        s.input_area.append_history.assert_called_once_with("x = 1")
        s.input_area.clear.assert_called_once()

    def test_mol_and_mw_synced_into_scope(self):
        fn = _run_code_fn()
        s = _console_self("x = 1", mol="THEMOL")
        mw = object()
        s.context.get_main_window.return_value = mw
        fn(s)
        assert s.local_scope["mol"] == "THEMOL"
        assert s.local_scope["mw"] is mw

    def test_none_mol_warning_when_command_uses_mol(self):
        fn = _run_code_fn()
        s = _console_self("mol", mol=None)
        fn(s)
        assert any("'mol' is None" in t for t, _ in s.outputs)

    def test_no_mol_warning_for_unrelated_command(self):
        fn = _run_code_fn()
        s = _console_self("x = 5", mol=None)
        fn(s)
        assert not any("'mol' is None" in t for t, _ in s.outputs)

    def test_input_echoed_with_prompt_markers(self):
        fn = _run_code_fn()
        s = _console_self("a = 1\nb = 2")
        fn(s)
        texts = [t for t, _ in s.outputs]
        assert ">>> a = 1" in texts
        assert "... b = 2" in texts
