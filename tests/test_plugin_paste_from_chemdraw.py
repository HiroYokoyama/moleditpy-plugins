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




# ---------------------------------------------------------------------------
# reconstruct_from_flat_text MolBlock output
# ---------------------------------------------------------------------------


def _reconstruct(flat_text):
    """Run reconstruct_from_flat_text and capture the MolBlock it builds."""
    with mock_optional_imports():
        mod = load_plugin(CHEMDRAW_PATH)
        captured = {}

        def capture(block):
            captured["block"] = block
            return "MOL_SENTINEL"

        mod.Chem.MolFromMolBlock = capture
        result = mod.reconstruct_from_flat_text(flat_text)
        return result, captured.get("block")


def _flat(n_atoms, n_bonds, body_tokens):
    header = f"molecule ChemDraw 2D  {n_atoms} {n_bonds} 0 0 0 0 0 0 0 0999 V2000"
    return header + " " + " ".join(body_tokens)


class TestChemDrawReconstruct:
    def test_no_v2000_returns_none(self):
        result, block = _reconstruct("just some random text")
        assert result is None
        assert block is None

    def test_short_header_returns_none(self):
        result, block = _reconstruct("1 1 V2000 rest")
        assert result is None
        assert block is None

    def test_non_numeric_counts_return_none(self):
        text = "a b c d e f g h i j k l V2000 stuff"
        result, block = _reconstruct(text)
        assert result is None

    def test_counts_line_format(self):
        tokens = ["0.0", "0.0", "0.0", "C"] + ["0"] * 12
        _, block = _reconstruct(_flat(1, 0, tokens))
        assert "  1  0  0  0  0  0  0  0  0  0999 V2000" in block

    def test_atom_line_formatting(self):
        tokens = ["1.5", "-2.25", "0.0", "N"] + ["0"] * 12
        _, block = _reconstruct(_flat(1, 0, tokens))
        atom_line = block.splitlines()[4]
        assert atom_line.startswith("    1.5000   -2.2500    0.0000 N")

    def test_phantom_f_token_skipped(self):
        tokens = ["F", "1.0", "2.0", "3.0", "O"] + ["0"] * 12
        _, block = _reconstruct(_flat(1, 0, tokens))
        atom_line = block.splitlines()[4]
        assert " O" in atom_line
        assert "1.0000" in atom_line

    def test_bond_line_built(self):
        tokens = (
            ["0.0", "0.0", "0.0", "C"]
            + ["0"] * 12
            + ["1.0", "0.0", "0.0", "C"]
            + ["0"] * 12
            + ["1", "2", "2", "0"]
        )
        _, block = _reconstruct(_flat(2, 1, tokens))
        assert "  1  2  2  0  0  0  0" in block

    def test_m_chg_preserved_before_end(self):
        tokens = (
            ["0.0", "0.0", "0.0", "N"]
            + ["0"] * 12
            + ["M", "CHG", "1", "1", "1", "M", "END"]
        )
        _, block = _reconstruct(_flat(1, 0, tokens))
        lines = block.splitlines()
        chg_idx = next(i for i, l in enumerate(lines) if l.startswith("M  CHG"))
        assert lines[chg_idx + 1] == "M  END"

    def test_block_ends_with_m_end(self):
        tokens = ["0.0", "0.0", "0.0", "C"] + ["0"] * 12
        _, block = _reconstruct(_flat(1, 0, tokens))
        assert block.splitlines()[-1] == "M  END"

    def test_control_characters_sanitized(self):
        tokens = ["0.0", "0.0", "0.0", "C"] + ["0"] * 12
        dirty = "\x01\x16" + _flat(1, 0, tokens) + "\x02"
        result, block = _reconstruct(dirty)
        assert result == "MOL_SENTINEL"

    def test_returns_parser_result(self):
        tokens = ["0.0", "0.0", "0.0", "C"] + ["0"] * 12
        result, _ = _reconstruct(_flat(1, 0, tokens))
        assert result == "MOL_SENTINEL"

    def test_two_atoms_two_lines(self):
        tokens = (
            ["0.0", "0.0", "0.0", "C"]
            + ["0"] * 12
            + ["1.4", "0.0", "0.0", "O"]
            + ["0"] * 12
        )
        _, block = _reconstruct(_flat(2, 0, tokens))
        lines = block.splitlines()
        assert " C" in lines[4]
        assert " O" in lines[5]
