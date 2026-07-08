"""
Tests for the Python Console plugin: ConsoleInput.append_history, GUIHelp.__repr__.

Qt-inheriting classes become MagicMock when loaded under mock_optional_imports
(their metaclass is the mocked PyQt6 type). append_history is extracted from
the plugin source via AST and compiled as a standalone function so we can
test pure logic that does not touch self's Qt attributes.
"""

from __future__ import annotations

from pathlib import Path

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
