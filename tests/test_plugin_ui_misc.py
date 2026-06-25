"""
Unit tests for visible plugins — focused on pure / logic-only functions.

Each test exercises a specific function without launching a GUI.

When a plugin class inherits from a Qt base, the class creation delegates to
the mocked metaclass and becomes a MagicMock rather than a real class.
For those methods we extract the function body from the plugin source via AST
and compile it as a standalone function — no Qt dependency at runtime.

For classes that do NOT inherit from Qt (e.g. GUIHelp) we instantiate them
directly from the loaded module.

Covered plugins:
- Python Console   : ConsoleInput.append_history, GUIHelp.__repr__
- MS Spectrum Neo  : parse_formula_str, to_subscript, to_superscript, get_adduct_delta
- VDW Radii Overlay: load_settings / save_settings
- Settings Saver   : load_library / save_library / get_plugin_data_path
- Chat with Mol Neo: load_settings / save_settings / latex_to_html fallback
- Paste from ChemDraw: reconstruct_from_flat_text  (module-level fn)
"""

from __future__ import annotations

import ast
import json
import os
import textwrap
from pathlib import Path

import pytest

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CONSOLE_PATH = PLUGINS_DIR / "Python_Console" / "console.py"
MS_PATH = PLUGINS_DIR / "MS_Spectrum_Simulation_Neo" / "ms_spectrum_neo.py"
VDW_PATH = PLUGINS_DIR / "VDW_Radii_Overlay" / "vdw_radii_overlay.py"
SAVER_PATH = PLUGINS_DIR / "Settings_Saver" / "settings_saver.py"
CHAT_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo" / "chat_with_molecule_neo.py"
CHEMDRAW_PATH = PLUGINS_DIR / "Paste_from_ChemDraw" / "paste_chemdraw.py"

# ---------------------------------------------------------------------------
# Load each plugin once at collection time; deps mocked only during import
# ---------------------------------------------------------------------------

with mock_optional_imports():
    _console = load_plugin(CONSOLE_PATH)
    _ms = load_plugin(MS_PATH)
    _vdw = load_plugin(VDW_PATH)
    _saver = load_plugin(SAVER_PATH)
    _chat = load_plugin(CHAT_PATH)
    _chemdraw = load_plugin(CHEMDRAW_PATH)


# ---------------------------------------------------------------------------
# AST helper: extract a method from a class and compile as a standalone fn
#
# Qt-inheriting classes become MagicMock when loaded under mock_optional_imports
# (their metaclass is the mocked PyQt6 type).  This helper extracts the raw
# Python source of a method and compiles it without the class context so we
# can test pure logic that does not touch self's Qt attributes.
# ---------------------------------------------------------------------------

def _extract_method(source_path: Path, class_name: str, method_name: str, extra_globals: dict | None = None):
    """Return a callable built from ``class_name.method_name`` in *source_path*."""
    src = source_path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    lines = src.splitlines()

    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == method_name:
                    method_lines = lines[item.lineno - 1 : item.end_lineno]
                    method_src = textwrap.dedent("\n".join(method_lines))
                    ns: dict = dict(extra_globals or {})
                    exec(compile(method_src, f"{class_name}.{method_name}", "exec"), ns)  # noqa: S102
                    return ns[method_name]
    raise ValueError(f"{class_name}.{method_name} not found in {source_path}")


# ---------------------------------------------------------------------------
# Extract pure methods from Qt-inheriting classes
# ---------------------------------------------------------------------------

# ConsoleInput.append_history — only touches self.history / self.history_index
_append_history = _extract_method(CONSOLE_PATH, "ConsoleInput", "append_history")

# MSSpectrumDialog methods — self is unused in all of them
_parse_formula_str = _extract_method(MS_PATH, "MSSpectrumDialog", "parse_formula_str")
_to_subscript = _extract_method(MS_PATH, "MSSpectrumDialog", "to_subscript")
_to_superscript = _extract_method(MS_PATH, "MSSpectrumDialog", "to_superscript")
_get_adduct_delta = _extract_method(MS_PATH, "MSSpectrumDialog", "get_adduct_delta")


# ---------------------------------------------------------------------------
# Python Console — ConsoleInput.append_history
# ---------------------------------------------------------------------------

class _HistoryStub:
    """Minimal stand-in for ConsoleInput — only carries history state."""

    def __init__(self):
        self.history: list = []
        self.history_index: int = 0


class TestAppendHistory:
    def test_adds_first_entry(self):
        stub = _HistoryStub()
        _append_history(stub, "print('hello')")
        assert stub.history == ["print('hello')"]
        assert stub.history_index == 1

    def test_ignores_empty_string(self):
        stub = _HistoryStub()
        _append_history(stub, "")
        assert stub.history == []
        assert stub.history_index == 0

    def test_no_consecutive_duplicate(self):
        stub = _HistoryStub()
        _append_history(stub, "x = 1")
        _append_history(stub, "x = 1")
        assert stub.history == ["x = 1"]
        assert stub.history_index == 1

    def test_non_consecutive_duplicate_is_allowed(self):
        stub = _HistoryStub()
        for cmd in ["x = 1", "y = 2", "x = 1"]:
            _append_history(stub, cmd)
        assert stub.history == ["x = 1", "y = 2", "x = 1"]

    def test_index_always_points_past_end(self):
        stub = _HistoryStub()
        for cmd in ["a", "b", "c"]:
            _append_history(stub, cmd)
        assert stub.history_index == len(stub.history)

    def test_multiple_distinct_entries(self):
        stub = _HistoryStub()
        cmds = ["import math", "math.sqrt(4)", "print(2)"]
        for cmd in cmds:
            _append_history(stub, cmd)
        assert stub.history == cmds
        assert stub.history_index == 3


# ---------------------------------------------------------------------------
# Python Console — GUIHelp (pure Python class, no Qt base)
# ---------------------------------------------------------------------------

class TestGUIHelp:
    def test_repr_mentions_help_object(self):
        h = _console.GUIHelp()
        assert "help(object)" in repr(h)

    def test_repr_warns_about_freezing(self):
        r = repr(_console.GUIHelp())
        lowered = r.lower()
        assert "disabled" in lowered or "prevent" in lowered


# ---------------------------------------------------------------------------
# MS Spectrum Simulation Neo — parse_formula_str
# ---------------------------------------------------------------------------

class TestParseFormulaStr:
    def test_water(self):
        assert _parse_formula_str(None, "H2O") == {"H": 2, "O": 1}

    def test_glucose(self):
        assert _parse_formula_str(None, "C6H12O6") == {"C": 6, "H": 12, "O": 6}

    def test_single_atom_no_count(self):
        assert _parse_formula_str(None, "C") == {"C": 1}

    def test_parentheses_multiplier(self):
        assert _parse_formula_str(None, "(CH2)3") == {"C": 3, "H": 6}

    def test_multichar_element(self):
        assert _parse_formula_str(None, "NaCl") == {"Na": 1, "Cl": 1}

    def test_caffeine(self):
        assert _parse_formula_str(None, "C8H10N4O2") == {"C": 8, "H": 10, "N": 4, "O": 2}

    def test_empty_returns_falsy(self):
        result = _parse_formula_str(None, "")
        assert not result

    def test_invalid_char_returns_none(self):
        assert _parse_formula_str(None, "C6?H6") is None

    def test_charge_symbol_stripped(self):
        # '+' is in the tokenizer — should not produce None
        result = _parse_formula_str(None, "C6H5+")
        assert result is not None
        assert result.get("C") == 6


# ---------------------------------------------------------------------------
# MS Spectrum Simulation Neo — superscript / subscript helpers
# ---------------------------------------------------------------------------

class TestSuperSubscript:
    def test_subscript_digits(self):
        assert _to_subscript(None, "123") == "₁₂₃"

    def test_superscript_positive(self):
        assert _to_superscript(None, "2+") == "²⁺"

    def test_superscript_minus(self):
        assert _to_superscript(None, "1-") == "¹⁻"

    def test_unknown_chars_pass_through(self):
        assert _to_subscript(None, "abc") == "abc"

    def test_empty_string(self):
        assert _to_subscript(None, "") == ""


# ---------------------------------------------------------------------------
# MS Spectrum Simulation Neo — get_adduct_delta
# ---------------------------------------------------------------------------

class TestGetAdductDelta:
    # Positive mode
    def test_positive_m_no_delta(self):
        assert _get_adduct_delta(None, 0, "Positive", 1) == {}

    def test_positive_proton(self):
        assert _get_adduct_delta(None, 1, "Positive", 1) == {"H": 1}

    def test_positive_sodium(self):
        assert _get_adduct_delta(None, 2, "Positive", 1) == {"Na": 1}

    def test_positive_potassium(self):
        assert _get_adduct_delta(None, 3, "Positive", 1) == {"K": 1}

    def test_positive_double_charge_proton(self):
        assert _get_adduct_delta(None, 1, "Positive", 2) == {"H": 2}

    # Negative mode
    def test_negative_m_no_delta(self):
        assert _get_adduct_delta(None, 0, "Negative", 1) == {}

    def test_negative_deprotonation(self):
        assert _get_adduct_delta(None, 1, "Negative", 1) == {"H": -1}

    def test_negative_chloride(self):
        assert _get_adduct_delta(None, 2, "Negative", 1) == {"Cl": 1}

    def test_negative_formate(self):
        assert _get_adduct_delta(None, 3, "Negative", 1) == {"C": 1, "H": 1, "O": 2}


# ---------------------------------------------------------------------------
# VDW Radii Overlay — load_settings / save_settings
# ---------------------------------------------------------------------------

class TestVDWSettings:
    def test_round_trip_occupancy(self, tmp_path, monkeypatch):
        fresh = {"occupancy": 0.3, "resolution": 0.125}
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "vdw.json"))
        monkeypatch.setattr(_vdw, "_vdw_settings", fresh)
        fresh["occupancy"] = 0.75
        _vdw.save_settings()
        fresh["occupancy"] = 0.3  # reset
        _vdw.load_settings()
        assert fresh["occupancy"] == pytest.approx(0.75)

    def test_round_trip_resolution(self, tmp_path, monkeypatch):
        fresh = {"occupancy": 0.3, "resolution": 0.125}
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "vdw.json"))
        monkeypatch.setattr(_vdw, "_vdw_settings", fresh)
        fresh["resolution"] = 0.25
        _vdw.save_settings()
        fresh["resolution"] = 0.125
        _vdw.load_settings()
        assert fresh["resolution"] == pytest.approx(0.25)

    def test_load_missing_file_is_noop(self, tmp_path, monkeypatch):
        fresh = {"occupancy": 0.3, "resolution": 0.125}
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(tmp_path / "nonexistent.json"))
        monkeypatch.setattr(_vdw, "_vdw_settings", fresh)
        _vdw.load_settings()  # must not raise
        assert fresh["occupancy"] == pytest.approx(0.3)

    def test_saved_file_is_valid_json(self, tmp_path, monkeypatch):
        fresh = {"occupancy": 0.5, "resolution": 0.2}
        path = tmp_path / "vdw.json"
        monkeypatch.setattr(_vdw, "SETTINGS_FILE", str(path))
        monkeypatch.setattr(_vdw, "_vdw_settings", fresh)
        _vdw.save_settings()
        data = json.loads(path.read_text())
        assert "occupancy" in data
        assert "resolution" in data


# ---------------------------------------------------------------------------
# Settings Saver — load_library / save_library / get_plugin_data_path
# ---------------------------------------------------------------------------

class TestSettingsSaver:
    def _patch_path(self, tmp_path, monkeypatch) -> str:
        data_path = str(tmp_path / "settings_saver.json")
        monkeypatch.setattr(_saver, "get_plugin_data_path", lambda: data_path)
        return data_path

    def test_round_trip(self, tmp_path, monkeypatch):
        self._patch_path(tmp_path, monkeypatch)
        payload = {
            "preset_A": {"opacity": 0.5},
            "_PLUGIN_CONFIG": {"always_save_to_project": True},
        }
        _saver.save_library(payload)
        loaded = _saver.load_library()
        assert loaded["preset_A"]["opacity"] == pytest.approx(0.5)
        assert loaded["_PLUGIN_CONFIG"]["always_save_to_project"] is True

    def test_load_missing_returns_empty_dict(self, tmp_path, monkeypatch):
        self._patch_path(tmp_path, monkeypatch)
        assert _saver.load_library() == {}

    def test_save_overwrites_existing(self, tmp_path, monkeypatch):
        self._patch_path(tmp_path, monkeypatch)
        _saver.save_library({"old": 1})
        _saver.save_library({"new": 2})
        assert _saver.load_library() == {"new": 2}

    def test_get_plugin_data_path_returns_json(self):
        path = _saver.get_plugin_data_path()
        assert isinstance(path, str)
        assert path.endswith(".json")
        assert os.path.basename(path) == _saver.SETTINGS_FILENAME


# ---------------------------------------------------------------------------
# Chat with Molecule Neo (Gemini) — load_settings / save_settings / latex_to_html
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Paste from ChemDraw — reconstruct_from_flat_text (module-level function)
# ---------------------------------------------------------------------------

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
