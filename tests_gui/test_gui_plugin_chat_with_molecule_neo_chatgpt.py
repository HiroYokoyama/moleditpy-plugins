"""
GUI tests for the Chat with Molecule Neo (ChatGPT) plugin.

Covers: ChatMoleculeWindow.

Chemistry/AI libs are mocked; real PyQt6 is used. ChatGPT is safe with
context.get_main_window()=None.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CHATGPT_PATH = (
    PLUGINS_DIR / "Chat_with_Molecule_Neo_ChatGPT" / "chat_with_molecule_neo_chatGPT.py"
)

with mock_chemistry_imports():
    _chatgpt = load_plugin_for_gui(CHATGPT_PATH)


def _no_window_context() -> MagicMock:
    """Context whose get_main_window() returns None — safe for ChatGPT/Gemini."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# ChatGPT variant
# ===========================================================================


class TestChatMoleculeWindowChatGPT:
    """ChatMoleculeWindow — Chat with Molecule Neo (ChatGPT)."""

    @pytest.fixture
    def win(self, qapp, monkeypatch):
        # The 100ms singleShot in __init__ would fire initialize_session on
        # a dead widget during a later test's event loop (and it crashes on
        # mocked markdown anyway) — neutralize it before construction.
        monkeypatch.setattr(
            _chatgpt.ChatMoleculeWindow, "initialize_session", lambda self: None
        )
        ctx = _no_window_context()
        w = _chatgpt.ChatMoleculeWindow(context=ctx)
        yield w
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "Chat with Molecule Neo (ChatGPT)"

    def test_api_key_field_is_password(self, win):
        from PyQt6.QtWidgets import QLineEdit
        assert win.txt_api_key.echoMode() == QLineEdit.EchoMode.Password

    def test_api_key_placeholder(self, win):
        assert "OpenAI" in win.txt_api_key.placeholderText()

    def test_model_combo_default(self, win):
        assert "gpt" in win.combo_model.currentText().lower()

    def test_chat_display_is_readonly(self, win):
        assert win.chat_display.isReadOnly()

    def test_input_initially_disabled(self, win):
        assert not win.txt_input.isEnabled()

    def test_send_button_initially_disabled(self, win):
        assert not win.btn_send.isEnabled()

    def test_context_label_initial_text(self, win):
        assert "No molecule" in win.lbl_context.text()

    def test_loading_bar_initially_hidden(self, win):
        assert win.loading_bar.isHidden()

    def test_thinking_label_initially_hidden(self, win):
        assert win.lbl_thinking.isHidden()

    def test_chat_history_log_initially_empty(self, win):
        assert win.chat_history_log == []

    def test_worker_initially_none(self, win):
        assert win.worker is None
