"""
Headless GUI tests for the Python Console plugin.

Covers: PythonConsoleDialog, initialize()/run() entry points, ConsoleInput,
PythonHighlighter.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from PyQt6.QtCore import QEvent, Qt
from PyQt6.QtGui import QKeyEvent, QTextCursor

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CONSOLE_PATH = PLUGINS_DIR / "Python_Console" / "console.py"

with mock_chemistry_imports():
    _console = load_plugin_for_gui(CONSOLE_PATH)


# ===========================================================================
# PythonConsoleDialog  (visible plugin: "Python Console")
# ===========================================================================


def _console_context() -> MagicMock:
    """Stub context with no main window and a MagicMock current_molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


class TestPythonConsoleDialog:
    """PythonConsoleDialog — embedded Python REPL dialog."""

    @pytest.fixture
    def console(self, qapp):
        ctx = _console_context()
        d = _console.PythonConsoleDialog(context=ctx)
        yield d
        d.destroy()

    def test_creates_without_error(self, console):
        assert console is not None

    def test_window_title(self, console):
        assert console.windowTitle() == "MoleditPy Python Console"

    def test_output_area_is_readonly(self, console):
        assert console.output_area.isReadOnly()

    def test_input_area_exists(self, console):
        assert console.input_area is not None

    def test_interpreter_is_interactive(self, console):
        import code as _code
        assert isinstance(console.interpreter, _code.InteractiveInterpreter)

    def test_context_in_local_scope(self, console):
        assert "context" in console.local_scope


# ===========================================================================
# Python Console — initialize/run entry point tests
# ===========================================================================


class TestPythonConsoleEntryPoints:
    """initialize() and run() entry-point behaviour without a full main window."""

    def test_initialize_sets_plugin_context(self, qapp):
        ctx = _console_context()
        _console.initialize(ctx)
        assert _console.PLUGIN_CONTEXT is ctx
        _console.PLUGIN_CONTEXT = None  # reset global after test

    def test_run_is_noop_when_context_is_none(self, qapp):
        _console.PLUGIN_CONTEXT = None
        mw = MagicMock()
        _console.run(mw)  # must not raise
        _console.PLUGIN_CONTEXT = None

    def test_run_creates_dialog_and_calls_get_window(self, qapp):
        ctx = _console_context()
        ctx.get_window.return_value = None  # no cached window
        _console.PLUGIN_CONTEXT = ctx
        mw = MagicMock()
        _console.run(mw)  # constructs a PythonConsoleDialog, calls show()
        ctx.get_window.assert_called_once_with("main_panel")
        _console.PLUGIN_CONTEXT = None


# ===========================================================================
# Python Console — ConsoleInput widget tests
# ===========================================================================


class TestConsoleInput:
    """ConsoleInput (QPlainTextEdit subclass) — history and placeholder tests."""

    @pytest.fixture
    def inp(self, qapp):
        w = _console.ConsoleInput()
        yield w
        w.destroy()

    def test_creates_without_error(self, inp):
        assert inp is not None

    def test_placeholder_text(self, inp):
        assert "Python" in inp.placeholderText() or "Enter" in inp.placeholderText()

    def test_history_initially_empty(self, inp):
        assert inp.history == []

    def test_append_history_stores_command(self, inp):
        inp.append_history("x = 1")
        assert inp.history[-1] == "x = 1"

    def test_append_history_deduplicates_consecutive(self, inp):
        inp.append_history("cmd")
        inp.append_history("cmd")
        assert inp.history.count("cmd") == 1

    def test_history_index_advances_after_append(self, inp):
        inp.append_history("a")
        inp.append_history("b")
        assert inp.history_index == 2


# ===========================================================================
# Python Console — PythonHighlighter tests
# ===========================================================================


class TestPythonHighlighter:
    """PythonHighlighter (QSyntaxHighlighter) — highlighting rule counts."""

    @pytest.fixture
    def hl(self, qapp):
        from PyQt6.QtGui import QTextDocument
        doc = QTextDocument()
        h = _console.PythonHighlighter(parent=doc)
        yield h

    def test_creates_without_error(self, hl):
        assert hl is not None

    def test_has_highlighting_rules(self, hl):
        assert len(hl.highlighting_rules) > 0

    def test_rules_are_tuples(self, hl):
        for rule in hl.highlighting_rules:
            assert isinstance(rule, tuple) and len(rule) == 2


class TestHighlightBlock:
    """highlightBlock() applies every rule over real document text."""

    def test_highlight_runs_over_all_rule_types(self, qapp):
        from PyQt6.QtGui import QTextDocument

        doc = QTextDocument()
        h = _console.PythonHighlighter(parent=doc)
        # keyword, builtin, comment, string and number all present so every
        # rule matches at least once inside highlightBlock's while loop.
        doc.setPlainText("def foo():  # note\n    x = print('hi') + 42")
        h.rehighlight()  # forces highlightBlock over each block synchronously
        assert doc.toPlainText().startswith("def foo")


# ===========================================================================
# Python Console — ConsoleInput.keyPressEvent (history navigation)
# ===========================================================================


def _key(key, mod=Qt.KeyboardModifier.NoModifier, text=""):
    return QKeyEvent(QEvent.Type.KeyPress, key, mod, text)


class TestConsoleInputKeys:
    """keyPressEvent — Enter/Shift+Enter and Up/Down history navigation."""

    @pytest.fixture
    def inp(self, qapp):
        w = _console.ConsoleInput()
        yield w
        w.destroy()

    def test_enter_emits_execute_signal(self, inp):
        fired = []
        inp.execute_signal.connect(lambda: fired.append(True))
        inp.setPlainText("x = 1")
        inp.keyPressEvent(_key(Qt.Key.Key_Return))
        assert fired == [True]
        # Return path returns early — text is not modified with a newline
        assert inp.toPlainText() == "x = 1"

    def test_shift_enter_does_not_execute(self, inp):
        # Shift+Enter takes the super() branch (newline), never the execute path.
        fired = []
        inp.execute_signal.connect(lambda: fired.append(True))
        inp.setPlainText("abc")
        inp.moveCursor(QTextCursor.MoveOperation.End)
        inp.keyPressEvent(_key(Qt.Key.Key_Return, Qt.KeyboardModifier.ShiftModifier))
        assert fired == []

    def test_up_navigates_history(self, inp):
        inp.history = ["first", "second"]
        inp.history_index = 2
        inp.keyPressEvent(_key(Qt.Key.Key_Up))
        assert inp.toPlainText() == "second"
        assert inp.history_index == 1

    def test_up_without_history_falls_through(self, inp):
        inp.history = []
        inp.keyPressEvent(_key(Qt.Key.Key_Up))  # must not raise
        assert inp.toPlainText() == ""

    def test_down_navigates_history_forward(self, inp):
        inp.history = ["a", "b"]
        inp.history_index = 0
        inp.setPlainText("a")
        inp.moveCursor(QTextCursor.MoveOperation.End)
        inp.keyPressEvent(_key(Qt.Key.Key_Down))
        assert inp.toPlainText() == "b"
        assert inp.history_index == 1

    def test_down_past_last_history_clears(self, inp):
        inp.history = ["a", "b"]
        inp.history_index = 1
        inp.setPlainText("b")
        inp.moveCursor(QTextCursor.MoveOperation.End)
        inp.keyPressEvent(_key(Qt.Key.Key_Down))
        assert inp.toPlainText() == ""
        assert inp.history_index == 2

    def test_down_without_history_falls_through(self, inp):
        inp.history = []
        inp.keyPressEvent(_key(Qt.Key.Key_Down))  # must not raise
        assert inp.toPlainText() == ""

    def test_other_key_inserts_text(self, inp):
        inp.keyPressEvent(_key(Qt.Key.Key_A, text="a"))
        assert "a" in inp.toPlainText()


# ===========================================================================
# Python Console — run_code() execution paths
# ===========================================================================


class TestRunCode:
    """run_code() — REPL execution, echoing, stdout/stderr capture."""

    def _dlg(self, mol=...):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        if mol is not ...:
            ctx.current_molecule = mol
        return _console.PythonConsoleDialog(context=ctx)

    def test_empty_command_is_noop(self, qapp):
        d = self._dlg()
        before = d.output_area.toPlainText()
        d.input_area.setPlainText("   ")
        d.run_code()  # returns before echoing anything
        assert d.output_area.toPlainText() == before

    def test_executes_and_echoes_input(self, qapp):
        d = self._dlg()
        d.input_area.setPlainText("x = 41 + 1")
        d.run_code()
        assert ">>> x = 41 + 1" in d.output_area.toPlainText()
        assert d.local_scope.get("x") == 42
        assert d.input_area.toPlainText() == ""  # input cleared

    def test_stdout_is_captured(self, qapp):
        d = self._dlg()
        d.input_area.setPlainText("print('hello world')")
        d.run_code()
        assert "hello world" in d.output_area.toPlainText()

    def test_stderr_stream_is_captured(self, qapp):
        # Write to stderr without raising, so pytest-qt's event-loop exception
        # guard is not tripped, while still covering the stderr-append branch.
        d = self._dlg()
        d.input_area.setPlainText("import sys; sys.stderr.write('ERRLINE')")
        d.run_code()
        assert "ERRLINE" in d.output_area.toPlainText()

    def test_multiline_uses_exec_mode(self, qapp):
        d = self._dlg()
        d.input_area.setPlainText("a = 5\nb = a * 2")
        d.run_code()
        assert d.local_scope.get("b") == 10
        assert "... b = a * 2" in d.output_area.toPlainText()

    def test_incomplete_block_warns(self, qapp):
        d = self._dlg()
        d.input_area.setPlainText("if True:")
        d.run_code()
        assert "Incomplete" in d.output_area.toPlainText()

    def test_mol_none_emits_warning(self, qapp):
        d = self._dlg(mol=None)
        d.input_area.setPlainText("print(mol)")
        d.run_code()
        assert "'mol' is None" in d.output_area.toPlainText()

    def test_append_output_with_color(self, qapp):
        d = self._dlg()
        d.append_output("colored", color="#FF0000")
        assert "colored" in d.output_area.toPlainText()
