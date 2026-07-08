"""
Tests for the Chat with Molecule Neo (Gemini) plugin: load_settings /
save_settings / latex_to_html fallback.
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CHAT_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo" / "chat_with_molecule_neo.py"

with mock_optional_imports():
    _chat = load_plugin(CHAT_PATH)


class TestChatWithMoleculeNeoSettings:
    def test_round_trip(self, tmp_path, monkeypatch):
        p = str(tmp_path / "chat.json")
        monkeypatch.setattr(_chat, "SETTINGS_FILE", p)
        monkeypatch.setattr(_chat, "DEMO_MODE", False)
        _chat.save_settings({"api_key": "tok123", "model": "gemini-pro"})
        loaded = _chat.load_settings()
        assert loaded["api_key"] == "tok123"
        assert loaded["model"] == "gemini-pro"

    def test_load_missing_returns_empty(self, tmp_path, monkeypatch):
        monkeypatch.setattr(_chat, "SETTINGS_FILE", str(tmp_path / "missing.json"))
        monkeypatch.setattr(_chat, "DEMO_MODE", False)
        assert _chat.load_settings() == {}

    def test_demo_mode_save_does_not_write(self, tmp_path, monkeypatch):
        p = tmp_path / "chat.json"
        monkeypatch.setattr(_chat, "SETTINGS_FILE", str(p))
        monkeypatch.setattr(_chat, "DEMO_MODE", True)
        _chat.save_settings({"key": "val"})
        assert not p.exists()

    def test_demo_mode_load_returns_empty(self, tmp_path, monkeypatch):
        p = tmp_path / "chat.json"
        p.write_text('{"key": "val"}', encoding="utf-8")
        monkeypatch.setattr(_chat, "SETTINGS_FILE", str(p))
        monkeypatch.setattr(_chat, "DEMO_MODE", True)
        assert _chat.load_settings() == {}

    def test_latex_to_html_fallback_returns_italic(self, monkeypatch):
        monkeypatch.setattr(_chat, "HAS_MATPLOTLIB", False)
        latex = r"$E = mc^2$"
        result = _chat.latex_to_html(latex)
        assert result.startswith("<i>")
        assert latex in result
