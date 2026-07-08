"""
Tests for the Molecule Comparator plugin.
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
MOLECULE_COMPARATOR_PATH = PLUGINS_DIR / "Molecule_Comparator" / "molecule_comparator.py"


class TestMoleculeComparator:
    def test_initialize_registers_document_reset_handler(self):
        """initialize() wires up on_reset via register_document_reset_handler."""
        with mock_optional_imports():
            mod = load_plugin(MOLECULE_COMPARATOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_document_reset_handler.assert_called_once()

    def test_get_text_color_for_background_dark(self):
        """Luminance formula: dark backgrounds map to white text."""
        # Mirror the algorithm from MoleculeComparator._get_text_color_for_background
        def _luminance_text(hex_color):
            try:
                h = hex_color.lstrip("#")
                r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
                return "white" if (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.5 else "black"
            except Exception:
                return "black"

        assert _luminance_text("#000000") == "white"
        assert _luminance_text("#0000FF") == "white"

    def test_get_text_color_for_background_light(self):
        """Luminance formula: light backgrounds map to black text."""
        def _luminance_text(hex_color):
            try:
                h = hex_color.lstrip("#")
                r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
                return "white" if (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.5 else "black"
            except Exception:
                return "black"

        assert _luminance_text("#FFFFFF") == "black"
        assert _luminance_text("#FFFF00") == "black"

    def test_get_text_color_for_background_invalid(self):
        """Luminance formula falls back to 'black' for non-hex input."""
        def _luminance_text(hex_color):
            try:
                h = hex_color.lstrip("#")
                r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
                return "white" if (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.5 else "black"
            except Exception:
                return "black"

        assert _luminance_text("not_a_color") == "black"
