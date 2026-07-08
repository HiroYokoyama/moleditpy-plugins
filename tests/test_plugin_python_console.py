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
from unittest.mock import MagicMock

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
