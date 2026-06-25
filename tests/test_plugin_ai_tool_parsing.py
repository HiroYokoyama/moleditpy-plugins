"""
Tests for AI tool-call parsing / response-processing logic shared across the
three Chat with Molecule Neo plugins (Gemini, ChatGPT, Local).

No network, no LLM, no Qt required — all heavy deps are mocked via
mock_optional_imports().  Pure-stdlib regex / json sections need no plugin
load at all.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

# ---------------------------------------------------------------------------
# Plugin paths
# ---------------------------------------------------------------------------

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

GEMINI_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo" / "chat_with_molecule_neo.py"
CHATGPT_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo_ChatGPT" / "chat_with_molecule_neo_chatGPT.py"
LOCAL_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo_Local" / "chat_with_molecule_neo_local.py"

# ---------------------------------------------------------------------------
# Load all three plugins once at module scope
# ---------------------------------------------------------------------------

with mock_optional_imports():
    _gemini = load_plugin(GEMINI_PATH)
    _chatgpt = load_plugin(CHATGPT_PATH)
    _local = load_plugin(LOCAL_PATH)

_PLUGINS = [
    ("Gemini", _gemini),
    ("ChatGPT", _chatgpt),
    ("Local", _local),
]

# ---------------------------------------------------------------------------
# Shared regex / logic (no plugin needed)
# ---------------------------------------------------------------------------

TOOL_REGEX = r"```json\s*([\[{].*?[\]}])\s*```"


def _collect_tools(text: str) -> list[dict]:
    """Mirror of the tool-collection logic used in all three chat plugins."""
    json_matches = re.findall(TOOL_REGEX, text, re.DOTALL)
    all_tools: list[dict] = []
    for json_str in json_matches:
        try:
            payload = json.loads(json_str)
            if isinstance(payload, list):
                all_tools.extend(payload)
            elif "tool" in payload and "params" in payload:
                all_tools.append(payload)
        except json.JSONDecodeError:
            pass
    return all_tools


# ---------------------------------------------------------------------------
# 1. Tool-call extraction regex
# ---------------------------------------------------------------------------


class TestToolRegex:
    def _blocks(self, text: str) -> list[str]:
        return re.findall(TOOL_REGEX, text, re.DOTALL)

    def test_single_tool_object_one_match(self):
        text = '```json\n{"tool": "apply_transformation", "params": {"atom_index": 5}}\n```'
        matches = self._blocks(text)
        assert len(matches) == 1
        payload = json.loads(matches[0])
        assert payload["tool"] == "apply_transformation"
        assert payload["params"]["atom_index"] == 5

    def test_array_of_tools_one_match(self):
        text = '```json\n[{"tool": "t1", "params": {}}, {"tool": "t2", "params": {}}]\n```'
        matches = self._blocks(text)
        assert len(matches) == 1
        payload = json.loads(matches[0])
        assert isinstance(payload, list)
        assert len(payload) == 2

    def test_two_separate_blocks_two_matches(self):
        text = (
            '```json\n{"tool": "t1", "params": {}}\n```\n'
            'some text\n'
            '```json\n{"tool": "t2", "params": {}}\n```'
        )
        assert len(self._blocks(text)) == 2

    def test_non_tool_json_found_but_no_tool_key(self):
        text = '```json\n{"key": "value", "other": 42}\n```'
        matches = self._blocks(text)
        assert len(matches) == 1
        payload = json.loads(matches[0])
        assert "tool" not in payload

    def test_invalid_json_block_found_but_unparseable(self):
        # Must start with { and end with } for regex to match, but body is not valid JSON
        text = '```json\n{broken: "json here"}\n```'
        matches = self._blocks(text)
        assert len(matches) == 1
        with pytest.raises(json.JSONDecodeError):
            json.loads(matches[0])

    def test_no_json_blocks_zero_matches(self):
        text = "Here is my analysis. No tools needed. The molecule looks fine."
        assert self._blocks(text) == []

    def test_tool_without_params_found_but_incomplete(self):
        text = '```json\n{"tool": "convert_to_3d"}\n```'
        matches = self._blocks(text)
        assert len(matches) == 1
        payload = json.loads(matches[0])
        assert "tool" in payload
        assert "params" not in payload


# ---------------------------------------------------------------------------
# 2. Tool-call collection logic
# ---------------------------------------------------------------------------


class TestCollectTools:
    def test_single_valid_tool(self):
        text = '```json\n{"tool": "highlight_substructure", "params": {"atom_indices": [1, 2]}}\n```'
        tools = _collect_tools(text)
        assert len(tools) == 1
        assert tools[0]["tool"] == "highlight_substructure"

    def test_array_of_three_tools(self):
        payload = [
            {"tool": "t1", "params": {}},
            {"tool": "t2", "params": {"x": 1}},
            {"tool": "t3", "params": {"y": 2}},
        ]
        text = f"```json\n{json.dumps(payload)}\n```"
        tools = _collect_tools(text)
        assert len(tools) == 3
        assert [t["tool"] for t in tools] == ["t1", "t2", "t3"]

    def test_two_separate_blocks_all_collected(self):
        text = (
            '```json\n{"tool": "apply_transformation", "params": {"reaction_smarts": "[c:1][H]>>[c:1][Cl]", "atom_index": 1}}\n```\n'
            'Converting...\n'
            '```json\n{"tool": "convert_to_3d", "params": {}}\n```'
        )
        tools = _collect_tools(text)
        assert len(tools) == 2
        assert tools[0]["tool"] == "apply_transformation"
        assert tools[1]["tool"] == "convert_to_3d"

    def test_invalid_json_block_ignored_valid_still_collected(self):
        # First block starts/ends with braces (matches regex) but body is not valid JSON
        text = (
            '```json\n{broken: "json here"}\n```\n'
            '```json\n{"tool": "clear_canvas", "params": {}}\n```'
        )
        tools = _collect_tools(text)
        assert len(tools) == 1
        assert tools[0]["tool"] == "clear_canvas"

    def test_non_tool_json_ignored(self):
        text = '```json\n{"result": "success", "value": 42}\n```'
        tools = _collect_tools(text)
        assert tools == []

    def test_tool_without_params_not_collected(self):
        text = '```json\n{"tool": "convert_to_3d"}\n```'
        tools = _collect_tools(text)
        assert tools == []

    def test_empty_text_no_tools(self):
        assert _collect_tools("") == []

    def test_mixed_tool_and_non_tool_blocks(self):
        text = (
            '```json\n{"info": "metadata"}\n```\n'
            '```json\n{"tool": "calculate_descriptors", "params": {"properties": ["MW"]}}\n```'
        )
        tools = _collect_tools(text)
        assert len(tools) == 1
        assert tools[0]["tool"] == "calculate_descriptors"


# ---------------------------------------------------------------------------
# 3. SMILES link pattern → HTML anchor
# ---------------------------------------------------------------------------

SMILES_PATTERN = r"\[([^\]]+)\]\((smiles:[^\)]+)\)"


def _smiles_to_html(text: str) -> str:
    return re.sub(SMILES_PATTERN, r'<a href="\2">\1</a>', text)


class TestSmilesLinkPattern:
    def test_single_smiles_link_converted(self):
        result = _smiles_to_html("[Benzene](smiles:c1ccccc1)")
        assert result == '<a href="smiles:c1ccccc1">Benzene</a>'

    def test_multiple_smiles_links_all_converted(self):
        text = "[Benzene](smiles:c1ccccc1) and [Caffeine](smiles:Cn1cnc2c1c(=O)n(C)c(=O)n2C)"
        result = _smiles_to_html(text)
        assert '<a href="smiles:c1ccccc1">Benzene</a>' in result
        assert "Caffeine" in result
        assert "smiles:Cn1cnc" in result

    def test_https_link_not_converted(self):
        text = "[MoleditPy](https://example.com)"
        result = _smiles_to_html(text)
        assert result == text

    def test_no_links_unchanged(self):
        text = "Plain text with no markdown links at all."
        assert _smiles_to_html(text) == text

    def test_smiles_without_parens_preserved(self):
        # The regex [^\)]+ stops at the first ) so only use paren-free SMILES here.
        # Benzene c1ccccc1 has no parentheses and round-trips cleanly.
        result = _smiles_to_html("[Ethane](smiles:CC)")
        assert 'href="smiles:CC"' in result
        assert ">Ethane</a>" in result

    def test_mixed_smiles_and_regular_links(self):
        text = "[Benzene](smiles:c1ccccc1) see also [docs](https://example.com)"
        result = _smiles_to_html(text)
        assert '<a href="smiles:c1ccccc1">Benzene</a>' in result
        assert "[docs](https://example.com)" in result


# ---------------------------------------------------------------------------
# 4. latex_to_html fallback (matplotlib mocked → HAS_MATPLOTLIB = False)
# ---------------------------------------------------------------------------


class TestLatexToHtmlFallback:
    """
    matplotlib is mocked via mock_optional_imports() but the import succeeds
    (MagicMock, not ImportError), so HAS_MATPLOTLIB ends up True.  We
    monkeypatch it to False for tests that specifically exercise the fallback
    branch, and test the as-loaded behaviour separately.
    """

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_empty_string_returns_string(self, plugin_name, mod):
        result = mod.latex_to_html("")
        assert isinstance(result, str)

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_plain_text_wrapped_in_italic_when_no_matplotlib(self, plugin_name, mod, monkeypatch):
        monkeypatch.setattr(mod, "HAS_MATPLOTLIB", False)
        result = mod.latex_to_html("x")
        assert result == "<i>x</i>"

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_latex_dollar_math_wrapped_when_no_matplotlib(self, plugin_name, mod, monkeypatch):
        monkeypatch.setattr(mod, "HAS_MATPLOTLIB", False)
        result = mod.latex_to_html(r"$\alpha$")
        assert result == r"<i>$\alpha$</i>"

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_has_matplotlib_is_bool(self, plugin_name, mod):
        # matplotlib mock makes the import succeed, so HAS_MATPLOTLIB=True here.
        # The important invariant is that it is a plain bool.
        assert isinstance(mod.HAS_MATPLOTLIB, bool)

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_with_matplotlib_mocked_returns_string(self, plugin_name, mod):
        # HAS_MATPLOTLIB=True path: mocked matplotlib returns MagicMocks;
        # the function should still return a string (img tag or code fallback).
        result = mod.latex_to_html("test")
        assert isinstance(result, str)
        assert len(result) > 0


# ---------------------------------------------------------------------------
# 5. PubChemResolver early-return paths (no network)
# ---------------------------------------------------------------------------


class TestPubChemResolverEarlyReturn:
    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_resolve_inchikey_empty_string(self, plugin_name, mod):
        name, err = mod.PubChemResolver.resolve_inchikey_to_name("")
        assert name is None
        assert "Empty" in err or "empty" in err.lower()

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_resolve_inchikey_none(self, plugin_name, mod):
        name, err = mod.PubChemResolver.resolve_inchikey_to_name(None)
        assert name is None
        assert err  # non-empty error message

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_resolve_name_to_smiles_empty_string(self, plugin_name, mod):
        smiles, err = mod.PubChemResolver.resolve_name_to_smiles("")
        assert smiles is None
        assert "Empty" in err or "empty" in err.lower()

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_resolve_name_to_smiles_none(self, plugin_name, mod):
        smiles, err = mod.PubChemResolver.resolve_name_to_smiles(None)
        assert smiles is None
        assert err


# ---------------------------------------------------------------------------
# 6. Module-level constants
# ---------------------------------------------------------------------------


class TestModuleConstants:
    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_max_history_is_int(self, plugin_name, mod):
        assert isinstance(mod.MAX_HISTORY, int)

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_max_history_is_even(self, plugin_name, mod):
        assert mod.MAX_HISTORY % 2 == 0, "MAX_HISTORY should be even (user+AI pairs)"

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_max_history_at_least_two(self, plugin_name, mod):
        assert mod.MAX_HISTORY >= 2

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_demo_mode_is_false(self, plugin_name, mod):
        assert mod.DEMO_MODE is False

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_system_prompt_non_empty(self, plugin_name, mod):
        assert isinstance(mod.SYSTEM_PROMPT, str)
        assert len(mod.SYSTEM_PROMPT) > 100

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_system_prompt_contains_tool_instructions(self, plugin_name, mod):
        assert "tool" in mod.SYSTEM_PROMPT.lower()

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_generation_config_has_temperature(self, plugin_name, mod):
        if hasattr(mod, "GENERATION_CONFIG"):
            assert "temperature" in mod.GENERATION_CONFIG
            assert isinstance(mod.GENERATION_CONFIG["temperature"], float)

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_plugin_version_present(self, plugin_name, mod):
        assert hasattr(mod, "PLUGIN_VERSION")
        assert mod.PLUGIN_VERSION


# ---------------------------------------------------------------------------
# 7. append_log smoke (all three variants — same pattern)
# ---------------------------------------------------------------------------


class TestAppendLog:
    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_append_log_creates_file(self, plugin_name, mod, tmp_path, monkeypatch):
        log_file = tmp_path / f"{plugin_name}_test.log"
        monkeypatch.setattr(mod, "LOG_FILE", str(log_file))
        mod.append_log("TestSender", "Hello from test")
        assert log_file.exists()

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_append_log_writes_content(self, plugin_name, mod, tmp_path, monkeypatch):
        log_file = tmp_path / f"{plugin_name}_content.log"
        monkeypatch.setattr(mod, "LOG_FILE", str(log_file))
        mod.append_log("Alice", "Test message content")
        content = log_file.read_text(encoding="utf-8")
        assert "Alice" in content
        assert "Test message content" in content

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_append_log_appends(self, plugin_name, mod, tmp_path, monkeypatch):
        log_file = tmp_path / f"{plugin_name}_append.log"
        monkeypatch.setattr(mod, "LOG_FILE", str(log_file))
        mod.append_log("A", "first")
        mod.append_log("B", "second")
        content = log_file.read_text(encoding="utf-8")
        assert "first" in content
        assert "second" in content

    @pytest.mark.parametrize("plugin_name,mod", _PLUGINS)
    def test_append_log_bad_path_does_not_raise(self, plugin_name, mod, monkeypatch):
        monkeypatch.setattr(mod, "LOG_FILE", "/nonexistent_dir/no_such.log")
        mod.append_log("X", "Y")  # must not raise — exceptions are swallowed
