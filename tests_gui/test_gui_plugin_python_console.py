"""
Headless GUI tests for the Python Console plugin.

Covers: PythonConsoleDialog, initialize()/run() entry points, ConsoleInput,
PythonHighlighter.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

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
