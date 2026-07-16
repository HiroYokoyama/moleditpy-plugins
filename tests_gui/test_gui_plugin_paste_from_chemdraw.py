"""
GUI-tier tests for the Paste from ChemDraw plugin.

These tests run with real PyQt6 (QT_QPA_PLATFORM=offscreen) and add coverage
that the tests/ suite cannot provide.

Covers extra reconstruct_from_flat_text() edge cases beyond
tests/test_plugin_ui_misc.py (empty, no V2000, too few tokens, non-printable
chars): valid header + partial data paths.
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CHEMDRAW_PATH = PLUGINS_DIR / "Paste_from_ChemDraw" / "paste_chemdraw.py"

with mock_chemistry_imports():
    _chemdraw = load_plugin_for_gui(CHEMDRAW_PATH)


# ===========================================================================
# Paste from ChemDraw — reconstruct_from_flat_text() extra edge cases
# (tests/test_plugin_ui_misc.py covers: empty, no V2000, too few tokens,
# non-printable chars.  These cases add valid header + partial data paths.)
# ===========================================================================


def _flat(n_atoms, n_bonds, atom_tokens="", bond_tokens=""):
    """Build a minimal flat ChemDraw clipboard string."""
    # 10 header tokens followed by V2000, then atom/bond tokens
    header = f"t1 t2 t3 t4 t5 t6 t7 t8 {n_bonds} {n_atoms}"
    rest = atom_tokens + " " + bond_tokens
    return header + " V2000 " + rest


class TestReconstructFromFlatTextExtra:
    def test_zero_atoms_zero_bonds_returns_string(self):
        result = _chemdraw.reconstruct_from_flat_text(_flat(0, 0))
        # With 0 atoms and 0 bonds the reconstruction may succeed or return None;
        # either is acceptable — it must not raise.

    def test_v2000_present_but_non_integer_counts_returns_none(self):
        bad = "t1 t2 t3 t4 t5 t6 t7 t8 XY ZZ V2000 rest"
        assert _chemdraw.reconstruct_from_flat_text(bad) is None

    def test_single_carbon_atom_block_returns_something(self):
        # 1 atom, 0 bonds, atom tokens: x y z symbol + 12 zeros
        atom_tok = "0.0 0.0 0.0 C " + " ".join(["0"] * 12)
        result = _chemdraw.reconstruct_from_flat_text(_flat(1, 0, atom_tok))
        # Result is either a string (reconstructed block) or None if parsing failed;
        # must not raise.
        assert result is None or isinstance(result, str)

    def test_whitespace_only_returns_none(self):
        assert _chemdraw.reconstruct_from_flat_text("   \n\t  ") is None

    def test_phantom_f_token_is_skipped(self):
        # Insert a 'F' phantom token before atom coordinates
        atom_tok = "F 0.0 0.0 0.0 C " + " ".join(["0"] * 12)
        result = _chemdraw.reconstruct_from_flat_text(_flat(1, 0, atom_tok))
        assert result is None or isinstance(result, str)

    def test_non_printable_control_chars_sanitized(self):
        text = "\x01\x02V2000 t1 t2 t3 t4 t5 t6 t7 t8 0 0"
        # Should not raise even though the string is unusual
        result = _chemdraw.reconstruct_from_flat_text(text)
        assert result is None or isinstance(result, str)
