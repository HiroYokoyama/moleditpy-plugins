"""
Unit tests for Chat with Molecule Neo variants:
  - Chat with Molecule Neo (ChatGPT)
  - Chat with Molecule Neo (Local)

Tests cover pure module-level functions only (no Qt required):
  - load_settings / save_settings round-trips
  - latex_to_html fallback path (HAS_MATPLOTLIB=False)
  - PubChemResolver.resolve_name_to_smiles with empty input (no network)
  - run(main_window) smoke
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CHATGPT_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo_ChatGPT" / "chat_with_molecule_neo_chatGPT.py"
LOCAL_PATH   = PLUGINS_DIR / "Chat_with_Molecule_Neo_Local"   / "chat_with_molecule_neo_local.py"

# ---------------------------------------------------------------------------
# Load both plugins once at collection time (deps mocked during import)
# ---------------------------------------------------------------------------

with mock_optional_imports():
    _chatgpt = load_plugin(CHATGPT_PATH)
    _local   = load_plugin(LOCAL_PATH)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _reset_cache(mod):
    """Clear the per-module LATEX_CACHE between tests."""
    mod.LATEX_CACHE.clear()


# ---------------------------------------------------------------------------
# Chat with Molecule Neo (ChatGPT) — settings
# ---------------------------------------------------------------------------

class TestChatGPTSettings:
    def test_round_trip(self, tmp_path, monkeypatch):
        p = str(tmp_path / "chatgpt.json")
        monkeypatch.setattr(_chatgpt, "SETTINGS_FILE", p)
        monkeypatch.setattr(_chatgpt, "DEMO_MODE", False)
        _chatgpt.save_settings({"api_key": "sk-test", "model": "gpt-4"})
        loaded = _chatgpt.load_settings()
        assert loaded["api_key"] == "sk-test"
        assert loaded["model"] == "gpt-4"

    def test_load_missing_returns_empty(self, tmp_path, monkeypatch):
        monkeypatch.setattr(_chatgpt, "SETTINGS_FILE", str(tmp_path / "missing.json"))
        monkeypatch.setattr(_chatgpt, "DEMO_MODE", False)
        assert _chatgpt.load_settings() == {}

    def test_load_corrupt_returns_empty(self, tmp_path, monkeypatch):
        f = tmp_path / "bad.json"
        f.write_text("not json", encoding="utf-8")
        monkeypatch.setattr(_chatgpt, "SETTINGS_FILE", str(f))
        monkeypatch.setattr(_chatgpt, "DEMO_MODE", False)
        assert _chatgpt.load_settings() == {}

    def test_demo_mode_load_returns_empty(self, tmp_path, monkeypatch):
        f = tmp_path / "settings.json"
        f.write_text('{"api_key": "real"}', encoding="utf-8")
        monkeypatch.setattr(_chatgpt, "SETTINGS_FILE", str(f))
        monkeypatch.setattr(_chatgpt, "DEMO_MODE", True)
        assert _chatgpt.load_settings() == {}

    def test_demo_mode_save_does_not_write(self, tmp_path, monkeypatch):
        p = tmp_path / "settings.json"
        monkeypatch.setattr(_chatgpt, "SETTINGS_FILE", str(p))
        monkeypatch.setattr(_chatgpt, "DEMO_MODE", True)
        _chatgpt.save_settings({"key": "val"})
        assert not p.exists()

    def test_save_written_as_valid_json(self, tmp_path, monkeypatch):
        p = tmp_path / "settings.json"
        monkeypatch.setattr(_chatgpt, "SETTINGS_FILE", str(p))
        monkeypatch.setattr(_chatgpt, "DEMO_MODE", False)
        _chatgpt.save_settings({"model": "gpt-4o", "temperature": 0.7})
        data = json.loads(p.read_text(encoding="utf-8"))
        assert data["model"] == "gpt-4o"
        assert data["temperature"] == pytest.approx(0.7)


# ---------------------------------------------------------------------------
# Chat with Molecule Neo (ChatGPT) — latex_to_html
# ---------------------------------------------------------------------------

class TestChatGPTLatexToHtml:
    def test_fallback_no_matplotlib_returns_italic(self, monkeypatch):
        monkeypatch.setattr(_chatgpt, "HAS_MATPLOTLIB", False)
        _reset_cache(_chatgpt)
        result = _chatgpt.latex_to_html(r"$E = mc^2$")
        assert result.startswith("<i>")
        assert r"$E = mc^2$" in result

    def test_fallback_empty_string_does_not_raise(self, monkeypatch):
        monkeypatch.setattr(_chatgpt, "HAS_MATPLOTLIB", False)
        _reset_cache(_chatgpt)
        result = _chatgpt.latex_to_html("")
        assert isinstance(result, str)

    def test_fallback_plain_text_wrapped_in_italic(self, monkeypatch):
        monkeypatch.setattr(_chatgpt, "HAS_MATPLOTLIB", False)
        _reset_cache(_chatgpt)
        result = _chatgpt.latex_to_html("alpha")
        assert "<i>" in result

    def test_with_matplotlib_mocked_returns_string(self, monkeypatch):
        # HAS_MATPLOTLIB is True (matplotlib mocked at load time); result is
        # either an img tag or the <code> fallback — either way a str.
        _reset_cache(_chatgpt)
        result = _chatgpt.latex_to_html(r"$\alpha$")
        assert isinstance(result, str)

    def test_cache_hit_returns_same_value(self, monkeypatch):
        monkeypatch.setattr(_chatgpt, "HAS_MATPLOTLIB", False)
        _reset_cache(_chatgpt)
        first = _chatgpt.latex_to_html(r"$x$")
        second = _chatgpt.latex_to_html(r"$x$")
        assert first == second


# ---------------------------------------------------------------------------
# Chat with Molecule Neo (ChatGPT) — PubChemResolver
# ---------------------------------------------------------------------------

class TestChatGPTPubChemResolver:
    def test_empty_name_returns_error_immediately(self):
        smiles, err = _chatgpt.PubChemResolver.resolve_name_to_smiles("")
        assert smiles is None
        assert err is not None
        assert "empty" in err.lower() or "name" in err.lower()

    def test_none_name_returns_error_immediately(self):
        smiles, err = _chatgpt.PubChemResolver.resolve_name_to_smiles(None)
        assert smiles is None
        assert err is not None


# ---------------------------------------------------------------------------
# Chat with Molecule Neo (ChatGPT) — run() smoke
# ---------------------------------------------------------------------------

def test_chatgpt_run_does_not_raise():
    with mock_optional_imports():
        mod = load_plugin(CHATGPT_PATH)
        mw = MagicMock()
        mod.run(mw)  # must not raise


# ---------------------------------------------------------------------------
# Chat with Molecule Neo (Local) — settings
# ---------------------------------------------------------------------------

class TestLocalSettings:
    def test_round_trip(self, tmp_path, monkeypatch):
        p = str(tmp_path / "local.json")
        monkeypatch.setattr(_local, "SETTINGS_FILE", p)
        monkeypatch.setattr(_local, "DEMO_MODE", False)
        _local.save_settings({"model": "llama3", "endpoint": "http://localhost:11434"})
        loaded = _local.load_settings()
        assert loaded["model"] == "llama3"
        assert loaded["endpoint"] == "http://localhost:11434"

    def test_load_missing_returns_empty(self, tmp_path, monkeypatch):
        monkeypatch.setattr(_local, "SETTINGS_FILE", str(tmp_path / "missing.json"))
        monkeypatch.setattr(_local, "DEMO_MODE", False)
        assert _local.load_settings() == {}

    def test_load_corrupt_returns_empty(self, tmp_path, monkeypatch):
        f = tmp_path / "bad.json"
        f.write_text("{invalid", encoding="utf-8")
        monkeypatch.setattr(_local, "SETTINGS_FILE", str(f))
        monkeypatch.setattr(_local, "DEMO_MODE", False)
        assert _local.load_settings() == {}

    def test_demo_mode_load_returns_empty(self, tmp_path, monkeypatch):
        f = tmp_path / "settings.json"
        f.write_text('{"endpoint": "real"}', encoding="utf-8")
        monkeypatch.setattr(_local, "SETTINGS_FILE", str(f))
        monkeypatch.setattr(_local, "DEMO_MODE", True)
        assert _local.load_settings() == {}

    def test_demo_mode_save_does_not_write(self, tmp_path, monkeypatch):
        p = tmp_path / "settings.json"
        monkeypatch.setattr(_local, "SETTINGS_FILE", str(p))
        monkeypatch.setattr(_local, "DEMO_MODE", True)
        _local.save_settings({"key": "val"})
        assert not p.exists()

    def test_save_written_as_valid_json(self, tmp_path, monkeypatch):
        p = tmp_path / "settings.json"
        monkeypatch.setattr(_local, "SETTINGS_FILE", str(p))
        monkeypatch.setattr(_local, "DEMO_MODE", False)
        _local.save_settings({"model": "mistral", "temperature": 0.5})
        data = json.loads(p.read_text(encoding="utf-8"))
        assert data["model"] == "mistral"


# ---------------------------------------------------------------------------
# Chat with Molecule Neo (Local) — latex_to_html
# ---------------------------------------------------------------------------

class TestLocalLatexToHtml:
    def test_fallback_no_matplotlib_returns_italic(self, monkeypatch):
        monkeypatch.setattr(_local, "HAS_MATPLOTLIB", False)
        _reset_cache(_local)
        result = _local.latex_to_html(r"$\beta$")
        assert result.startswith("<i>")
        assert r"$\beta$" in result

    def test_fallback_empty_string_does_not_raise(self, monkeypatch):
        monkeypatch.setattr(_local, "HAS_MATPLOTLIB", False)
        _reset_cache(_local)
        result = _local.latex_to_html("")
        assert isinstance(result, str)

    def test_with_matplotlib_mocked_returns_string(self, monkeypatch):
        _reset_cache(_local)
        result = _local.latex_to_html(r"$\gamma$")
        assert isinstance(result, str)

    def test_cache_hit_returns_same_value(self, monkeypatch):
        monkeypatch.setattr(_local, "HAS_MATPLOTLIB", False)
        _reset_cache(_local)
        first = _local.latex_to_html(r"$y$")
        second = _local.latex_to_html(r"$y$")
        assert first == second


# ---------------------------------------------------------------------------
# Chat with Molecule Neo (Local) — PubChemResolver
# ---------------------------------------------------------------------------

class TestLocalPubChemResolver:
    def test_empty_name_returns_error_immediately(self):
        smiles, err = _local.PubChemResolver.resolve_name_to_smiles("")
        assert smiles is None
        assert err is not None

    def test_none_name_returns_error_immediately(self):
        smiles, err = _local.PubChemResolver.resolve_name_to_smiles(None)
        assert smiles is None
        assert err is not None


# ---------------------------------------------------------------------------
# Chat with Molecule Neo (Local) — run() smoke
# ---------------------------------------------------------------------------

def test_local_run_does_not_raise():
    with mock_optional_imports():
        mod = load_plugin(LOCAL_PATH)
        mw = MagicMock()
        mod.run(mw)  # must not raise
