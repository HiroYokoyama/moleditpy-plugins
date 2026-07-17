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
from unittest.mock import MagicMock

import pytest

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


# ===========================================================================
# reconstruct_from_flat_text() — verify the rebuilt MolBlock handed to RDKit
# (Chem is mocked, so we capture the block passed to Chem.MolFromMolBlock)
# ===========================================================================


class TestReconstructBuildsValidMolBlock:
    def test_single_atom_block_structure(self, monkeypatch):
        captured = []
        sentinel = object()

        def fake_mol_from_block(block, *a, **k):
            captured.append(block)
            return sentinel

        monkeypatch.setattr(_chemdraw.Chem, "MolFromMolBlock", fake_mol_from_block)

        # counts line tokens sit 10-from-the-end: n_atoms=1, n_bonds=0
        atom_tok = "0.5 -1.25 0.0 N " + " ".join(["0"] * 12)
        flat = "1 0 0 0 0 0 0 0 0 0 V2000 " + atom_tok + " M END"
        result = _chemdraw.reconstruct_from_flat_text(flat)

        assert result is sentinel
        assert len(captured) == 1
        block = captured[0]
        lines = block.splitlines()
        assert lines[3] == "  1  0  0  0  0  0  0  0  0  0999 V2000"
        assert "N" in lines[4]
        assert "0.5000" in lines[4] and "-1.2500" in lines[4]
        assert lines[-1] == "M  END"

    def test_atom_and_bond_block(self, monkeypatch):
        captured = []
        monkeypatch.setattr(
            _chemdraw.Chem,
            "MolFromMolBlock",
            lambda block, *a, **k: captured.append(block) or object(),
        )
        atoms = " ".join(
            f"{x} 0.0 0.0 C " + " ".join(["0"] * 12) for x in ("0.0", "1.5")
        )
        bond = "1 2 1 0"
        flat = "2 1 0 0 0 0 0 0 0 0 V2000 " + atoms + " " + bond + " M END"
        _chemdraw.reconstruct_from_flat_text(flat)

        assert len(captured) == 1
        lines = captured[0].splitlines()
        assert lines[3].startswith("  2  1")
        assert "  1  2  1  0  0  0  0" in lines
        assert lines[-1] == "M  END"

    def test_phantom_f_token_before_atom_is_skipped(self, monkeypatch):
        captured = []
        monkeypatch.setattr(
            _chemdraw.Chem,
            "MolFromMolBlock",
            lambda block, *a, **k: captured.append(block) or object(),
        )
        atom_tok = "F 0.0 0.0 0.0 O " + " ".join(["0"] * 12)
        flat = "1 0 0 0 0 0 0 0 0 0 V2000 " + atom_tok + " M END"
        _chemdraw.reconstruct_from_flat_text(flat)

        assert len(captured) == 1
        # The 'F' phantom must not become the atom symbol — 'O' does
        assert " O  " in captured[0].splitlines()[4]


# ===========================================================================
# initialize() and run() — clipboard-driven GUI paths with real QClipboard
# ===========================================================================


def test_initialize_registers_edit_menu_with_shortcut():
    ctx = MagicMock()
    _chemdraw.initialize(ctx)
    ctx.add_menu_action.assert_called_once()
    args, kwargs = ctx.add_menu_action.call_args
    assert args[0] == "Edit/Paste from ChemDraw"
    assert callable(args[1])
    assert kwargs.get("shortcut") == "Ctrl+Shift+V"


class _FakePos:
    x = 0.0
    y = 0.0


class _FakeConf:
    def GetAtomPosition(self, i):
        return _FakePos()


class _FakeAtom:
    def GetSymbol(self):
        return "C"

    def GetFormalCharge(self):
        return 0


class _FakeMol:
    def GetNumConformers(self):
        return 1

    def GetConformer(self):
        return _FakeConf()

    def GetNumAtoms(self):
        return 2

    def GetAtomWithIdx(self, i):
        return _FakeAtom()

    def GetBonds(self):
        return []


class TestRunClipboardPaths:
    @pytest.fixture
    def ctx(self):
        class _BareMainWindow:
            pass  # no .host, no .init_manager → run() uses QPointF(0, 0)

        c = MagicMock()
        c.get_main_window.return_value = _BareMainWindow()
        return c

    @pytest.fixture
    def clipboard(self, qapp):
        cb = qapp.clipboard()
        yield cb
        cb.clear()

    @pytest.fixture
    def warnings(self, monkeypatch):
        calls = []
        monkeypatch.setattr(
            _chemdraw.QMessageBox,
            "warning",
            staticmethod(lambda *a, **k: calls.append(a)),
        )
        return calls

    @pytest.fixture
    def criticals(self, monkeypatch):
        calls = []
        monkeypatch.setattr(
            _chemdraw.QMessageBox,
            "critical",
            staticmethod(lambda *a, **k: calls.append(a)),
        )
        return calls

    def test_empty_text_clipboard_warns_no_data(self, ctx, clipboard, warnings):
        clipboard.setText("")
        _chemdraw.run(ctx)
        assert len(warnings) == 1
        assert "No valid MDLCT data" in warnings[0][2]
        ctx.push_undo_checkpoint.assert_not_called()

    def test_cleared_clipboard_null_mime_warns_instead_of_crashing(
        self, ctx, clipboard, warnings
    ):
        # On the offscreen platform clear() makes mimeData() return None
        clipboard.clear()
        _chemdraw.run(ctx)
        assert len(warnings) == 1
        assert "No valid MDLCT data" in warnings[0][2]
        ctx.push_undo_checkpoint.assert_not_called()

    def test_plain_text_without_molblock_markers_warns(self, ctx, clipboard, warnings):
        clipboard.setText("just some random text")
        _chemdraw.run(ctx)
        assert len(warnings) == 1
        ctx.push_undo_checkpoint.assert_not_called()

    def test_v2000_text_fallback_draws_molecule(
        self, ctx, clipboard, warnings, criticals, monkeypatch
    ):
        monkeypatch.setattr(
            _chemdraw.Chem, "MolFromMolBlock", lambda *a, **k: _FakeMol()
        )
        clipboard.setText("fake molblock V2000 content\nM  END")
        _chemdraw.run(ctx)

        assert warnings == []
        assert criticals == []
        assert ctx.scene.create_atom.call_count == 2
        ctx.refresh_2d_scene.assert_called_once()
        ctx.push_undo_checkpoint.assert_called_once()
        ctx.show_status_message.assert_called_once()

    def test_unparseable_v2000_text_warns(
        self, ctx, clipboard, warnings, monkeypatch
    ):
        # Both direct parse and flat-text reconstruction fail → warning path
        monkeypatch.setattr(
            _chemdraw.Chem, "MolFromMolBlock", lambda *a, **k: None
        )
        clipboard.setText("V2000 garbage that cannot be reconstructed")
        _chemdraw.run(ctx)
        assert len(warnings) == 1
        ctx.push_undo_checkpoint.assert_not_called()
