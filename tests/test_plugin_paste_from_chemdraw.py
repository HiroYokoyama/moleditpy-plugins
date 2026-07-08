"""
Tests for the Paste from ChemDraw plugin: reconstruct_from_flat_text
(module-level function).
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CHEMDRAW_PATH = PLUGINS_DIR / "Paste_from_ChemDraw" / "paste_chemdraw.py"

with mock_optional_imports():
    _chemdraw = load_plugin(CHEMDRAW_PATH)


class TestReconstructFromFlatText:
    def test_empty_string_returns_none(self):
        assert _chemdraw.reconstruct_from_flat_text("") is None

    def test_no_v2000_marker_returns_none(self):
        assert _chemdraw.reconstruct_from_flat_text("some ordinary text") is None

    def test_too_few_header_tokens_returns_none(self):
        # V2000 present but fewer than 10 header tokens before it
        result = _chemdraw.reconstruct_from_flat_text("only three tokens V2000 rest")
        assert result is None

    def test_non_printable_chars_dont_raise(self):
        # Control characters should be sanitised; no V2000 → None
        result = _chemdraw.reconstruct_from_flat_text("\x01\x02\x16garbage")
        assert result is None
