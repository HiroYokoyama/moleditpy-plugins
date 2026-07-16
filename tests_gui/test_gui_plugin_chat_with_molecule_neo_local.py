"""
GUI tests for the Chat with Molecule Neo (Local) plugin.

Covers: ChatMoleculeWindow.

Chemistry/AI libs are mocked; real PyQt6 is used.

Construction notes:
- Local: accesses self.main_window.init_manager etc. in __init__, so needs a
  real QWidget (not None, not MagicMock) as the main window parent.  We use a
  bare QDialog() with manually attached MagicMock attributes.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

LOCAL_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo_Local" / "chat_with_molecule_neo_local.py"

with mock_chemistry_imports():
    _local = load_plugin_for_gui(LOCAL_PATH)

# markdown is mocked → markdown.markdown() returns MagicMock → re.sub crashes.
# Disable markdown rendering for all tests in this file.
_local.HAS_MARKDOWN = False


def _local_context(qapp):
    """
    Context for the Local variant.

    The Local plugin assigns self.main_window.init_manager etc. in __init__,
    so main_window must be a real QWidget (Qt rejects MagicMock as parent and
    None would crash on attribute access).  We attach MagicMock manager attrs
    to a bare QDialog instance.
    """
    from PyQt6.QtWidgets import QDialog

    fake_mw = QDialog()
    fake_mw.init_manager = MagicMock()
    fake_mw.state_manager = MagicMock()
    fake_mw.edit_actions_manager = MagicMock()
    fake_mw.view_3d_manager = MagicMock()
    ctx = MagicMock()
    ctx.get_main_window.return_value = fake_mw
    ctx.current_molecule = None
    return ctx, fake_mw


# ===========================================================================
# Local variant
# ===========================================================================


class TestChatMoleculeWindowLocal:
    """ChatMoleculeWindow — Chat with Molecule Neo (Local).

    Requires a real QWidget as main_window because __init__ assigns
    self.main_window.init_manager etc. before the try/except can catch it.
    """

    @pytest.fixture
    def win(self, qapp):
        ctx, fake_mw = _local_context(qapp)
        w = _local.ChatMoleculeWindow(context=ctx)
        yield w
        w.destroy()
        fake_mw.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "Chat with Molecule Neo (Local)"

    def test_api_base_url_placeholder(self, win):
        assert "localhost" in win.txt_api_base.placeholderText()

    def test_chat_display_is_readonly(self, win):
        assert win.chat_display.isReadOnly()

    def test_input_initially_disabled(self, win):
        assert not win.txt_input.isEnabled()

    def test_send_button_initially_disabled(self, win):
        assert not win.btn_send.isEnabled()

    def test_client_initially_none(self, win):
        assert win.client is None

    def test_worker_initially_none(self, win):
        assert win.worker is None

    def test_chat_history_log_initially_empty(self, win):
        assert win.chat_history_log == []

    def test_context_label_initial_text(self, win):
        assert "No molecule" in win.lbl_context.text()
