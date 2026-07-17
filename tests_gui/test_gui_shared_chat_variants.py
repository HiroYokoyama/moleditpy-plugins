"""
Shared GUI tests for the three Chat with Molecule Neo variants (Gemini,
ChatGPT, Local). Their window logic is copy-shared, so these behaviors are
parametrized across all three real modules in one file rather than
duplicated per variant.

Covers real-widget behavior the mocked suite doesn't reach: append_message
bookkeeping, the save_settings method (UI -> settings -> persist -> re-init),
and the export_history empty-guard / write paths.

Real PyQt6 (QT_QPA_PLATFORM=offscreen); chemistry/AI libs are mocked, so
markdown.markdown() returns a MagicMock — render_content is stubbed to a
plain string in each built window.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

_PATHS = {
    "gemini": PLUGINS_DIR / "Chat_with_Molecule_Neo" / "chat_with_molecule_neo.py",
    "chatgpt": PLUGINS_DIR
    / "Chat_with_Molecule_Neo_ChatGPT"
    / "chat_with_molecule_neo_chatGPT.py",
    "local": PLUGINS_DIR
    / "Chat_with_Molecule_Neo_Local"
    / "chat_with_molecule_neo_local.py",
}

with mock_chemistry_imports():
    _MODS = {name: load_plugin_for_gui(path) for name, path in _PATHS.items()}

# The Local variant assigns self.main_window.init_manager etc. in __init__,
# so it needs a real QWidget main window; Gemini/ChatGPT accept None.
_NEEDS_REAL_MW = {"gemini": False, "chatgpt": False, "local": True}


def _build(mod, monkeypatch, needs_real_mw):
    monkeypatch.setattr(
        mod.ChatMoleculeWindow, "initialize_session", lambda self: None
    )
    monkeypatch.setattr(mod, "append_log", lambda *a, **k: None)
    ctx = MagicMock()
    mw = None
    if needs_real_mw:
        from PyQt6.QtWidgets import QDialog

        mw = QDialog()
        mw.init_manager = MagicMock()
        mw.state_manager = MagicMock()
        mw.edit_actions_manager = MagicMock()
        mw.view_3d_manager = MagicMock()
        ctx.get_main_window.return_value = mw
    else:
        ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    w = mod.ChatMoleculeWindow(context=ctx)
    w.render_content = lambda text: text  # markdown lib is mocked
    return w, mw


@pytest.fixture(params=list(_MODS), ids=list(_MODS))
def chatwin(request, qapp, monkeypatch):
    name = request.param
    mod = _MODS[name]
    w, mw = _build(mod, monkeypatch, _NEEDS_REAL_MW[name])
    yield w, mod
    w.destroy()
    if mw is not None:
        mw.destroy()


# ===========================================================================
# append_message
# ===========================================================================


class TestSharedAppendMessage:
    def test_appends_to_history_log(self, chatwin):
        w, _ = chatwin
        w.append_message("You", "hello")
        assert w.chat_history_log[-1] == {"sender": "You", "text": "hello"}

    def test_text_shown_in_display(self, chatwin):
        w, _ = chatwin
        w.append_message("Assistant", "a-unique-reply-token")
        assert "a-unique-reply-token" in w.chat_display.toPlainText()

    def test_multiple_messages_accumulate(self, chatwin):
        w, _ = chatwin
        w.append_message("You", "one")
        w.append_message("Assistant", "two")
        assert [e["text"] for e in w.chat_history_log] == ["one", "two"]

    def test_system_message_recorded(self, chatwin):
        w, _ = chatwin
        w.append_message("System", "note", "green")
        assert w.chat_history_log[-1]["sender"] == "System"


# ===========================================================================
# save_settings (method)
# ===========================================================================


class TestSharedSaveSettings:
    def test_collects_api_key_and_model(self, chatwin, monkeypatch):
        w, mod = chatwin
        monkeypatch.setattr(mod, "save_settings", lambda settings: None)
        w.initialize_session = MagicMock()
        w.txt_api_key.setText("  KEY123  ")
        w.combo_model.clear()
        w.combo_model.addItem("model-x")
        w.save_settings()
        assert w.settings["api_key"] == "KEY123"  # trimmed
        assert w.settings["model"] == "model-x"

    def test_persists_via_module_save(self, chatwin, monkeypatch):
        w, mod = chatwin
        captured = {}
        monkeypatch.setattr(mod, "save_settings", captured.update)
        w.initialize_session = MagicMock()
        w.txt_api_key.setText("ABC")
        w.save_settings()
        assert captured.get("api_key") == "ABC"

    def test_reinitializes_session(self, chatwin, monkeypatch):
        w, mod = chatwin
        monkeypatch.setattr(mod, "save_settings", lambda settings: None)
        w.initialize_session = MagicMock()
        w.save_settings()
        w.initialize_session.assert_called_once()


# ===========================================================================
# export_history
# ===========================================================================


class TestSharedExportHistory:
    def test_empty_history_shows_info_and_writes_nothing(
        self, chatwin, monkeypatch, tmp_path
    ):
        w, mod = chatwin
        info = MagicMock()
        save = MagicMock(return_value=("", ""))
        monkeypatch.setattr(mod.QMessageBox, "information", info)
        monkeypatch.setattr(mod.QFileDialog, "getSaveFileName", save)
        w.chat_history_log = []
        w.export_history()
        info.assert_called_once()
        save.assert_not_called()
        assert list(tmp_path.iterdir()) == []

    def test_json_export_writes_history(self, chatwin, monkeypatch, tmp_path):
        w, mod = chatwin
        target = str(tmp_path / "hist.json")
        monkeypatch.setattr(
            mod.QFileDialog, "getSaveFileName", lambda *a, **k: (target, "")
        )
        monkeypatch.setattr(mod.QMessageBox, "information", MagicMock())
        w.chat_history_log = [{"sender": "You", "text": "hi"}]
        w.export_history()
        assert json.loads(Path(target).read_text(encoding="utf-8")) == [
            {"sender": "You", "text": "hi"}
        ]

    def test_md_export_contains_message_text(self, chatwin, monkeypatch, tmp_path):
        w, mod = chatwin
        target = str(tmp_path / "hist.md")
        monkeypatch.setattr(
            mod.QFileDialog, "getSaveFileName", lambda *a, **k: (target, "")
        )
        monkeypatch.setattr(mod.QMessageBox, "information", MagicMock())
        w.chat_history_log = [{"sender": "You", "text": "unique-md-token"}]
        w.export_history()
        assert "unique-md-token" in Path(target).read_text(encoding="utf-8")
