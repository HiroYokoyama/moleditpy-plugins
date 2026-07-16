"""
GUI tests for the Chat with Molecule Neo (Gemini) plugin.

Covers: ChatMoleculeWindow.

Chemistry/AI libs are mocked; real PyQt6 is used. Gemini is safe with
context.get_main_window()=None.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

GEMINI_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo" / "chat_with_molecule_neo.py"

with mock_chemistry_imports():
    _gemini = load_plugin_for_gui(GEMINI_PATH)


def _no_window_context() -> MagicMock:
    """Context whose get_main_window() returns None — safe for ChatGPT/Gemini."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# Gemini variant
# ===========================================================================


class TestChatMoleculeWindowGemini:
    """ChatMoleculeWindow — Chat with Molecule Neo (Gemini)."""

    @pytest.fixture
    def win(self, qapp, monkeypatch):
        # The 100ms singleShot in __init__ would fire initialize_session on
        # a dead widget during a later test's event loop (and it crashes on
        # mocked markdown anyway) — neutralize it before construction.
        monkeypatch.setattr(
            _gemini.ChatMoleculeWindow, "initialize_session", lambda self: None
        )
        ctx = _no_window_context()
        w = _gemini.ChatMoleculeWindow(context=ctx)
        yield w
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "Chat with Molecule Neo (Gemini)"

    def test_api_key_placeholder_mentions_gemini(self, win):
        assert "Gemini" in win.txt_api_key.placeholderText()

    def test_model_combo_has_item(self, win):
        assert win.combo_model.count() > 0

    def test_chat_display_is_readonly(self, win):
        assert win.chat_display.isReadOnly()

    def test_input_initially_disabled(self, win):
        assert not win.txt_input.isEnabled()

    def test_send_button_initially_disabled(self, win):
        assert not win.btn_send.isEnabled()

    def test_context_label_initial_text(self, win):
        assert "No molecule" in win.lbl_context.text()

    def test_chat_history_log_initially_empty(self, win):
        assert win.chat_history_log == []

    def test_worker_initially_none(self, win):
        assert win.worker is None

    def test_chat_session_initially_none(self, win):
        assert win.chat_session is None
