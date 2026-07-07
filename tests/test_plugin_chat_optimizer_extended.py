"""
Extended tests for the Chat with Molecule Neo variants (Gemini / ChatGPT /
Local), the PubChem resolver plugins, and the optimizer plugins
(All-Trans Optimizer, Complex Molecule Untangler, Conformational Search).

Focus areas NOT covered by the existing suites:
- Chat history pruning / usage logging (regression tests for the 2026.07.07
  fixes: Local variant never pruned, ChatGPT variant dropped the system
  prompt when pruning)
- Tool dispatch routing (_dispatch_tool) for all three chat variants
- PubChemResolver mocked-HTTP success / error paths (only the empty-input
  guards were covered before)
- Context-message assembly (_build_context_msg)
- Conformational Search energy-window dedup filter
- UntangleWorker Monte Carlo accept/reject logic
- PubChem Name Resolver run_search() network paths

No network, no LLM, no GUI — heavy deps are mocked via mock_optional_imports();
methods of Qt-derived classes are extracted standalone via AST (the Qt bases
are MagicMock, so the classes themselves are not instantiable).
"""

from __future__ import annotations

import ast
import textwrap
import types as _types
import urllib.error
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

GEMINI_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo" / "chat_with_molecule_neo.py"
CHATGPT_PATH = (
    PLUGINS_DIR / "Chat_with_Molecule_Neo_ChatGPT" / "chat_with_molecule_neo_chatGPT.py"
)
LOCAL_PATH = (
    PLUGINS_DIR / "Chat_with_Molecule_Neo_Local" / "chat_with_molecule_neo_local.py"
)
NAME_RESOLVER_PATH = PLUGINS_DIR / "PubChem_Name_Ressolver" / "pubchem_ressolver.py"
ALLTRANS_PATH = PLUGINS_DIR / "All-Trans_Optimizer" / "all-trans_optimizer.py"
UNTANGLER_PATH = (
    PLUGINS_DIR / "Complex_Molecule_Untangler" / "complex_molecule_untangler.py"
)
CONF_SEARCH_PATH = PLUGINS_DIR / "Conformational_Search" / "conf_search.py"

with mock_optional_imports():
    _gemini = load_plugin(GEMINI_PATH)
    _chatgpt = load_plugin(CHATGPT_PATH)
    _local = load_plugin(LOCAL_PATH)
    _alltrans = load_plugin(ALLTRANS_PATH)

_CHAT_VARIANTS = [
    ("Gemini", _gemini),
    ("ChatGPT", _chatgpt),
    ("Local", _local),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _extract_method_as_fn(
    path: Path,
    class_name: str,
    method_name: str,
    extra_globals: dict | None = None,
):
    """
    Use AST to extract a class method as a standalone callable.

    Needed because Qt base classes are MagicMock instances whose metaclass call
    returns a MagicMock instead of a real type, so the class definition doesn't
    produce a usable type and object.__new__ fails.
    """
    import logging

    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if (
                    isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and item.name == method_name
                ):
                    func_src = ast.get_source_segment(source, item)
                    if func_src:
                        local_ns: dict = {}
                        globs = {"logging": logging, **(extra_globals or {})}
                        exec(textwrap.dedent(func_src), globs, local_ns)
                        return local_ns[method_name]
    return None


class _Stub:
    """Attribute bucket without MagicMock auto-attribute behavior."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


# ---------------------------------------------------------------------------
# 1. History pruning — regression tests for the 2026.07.07 fixes
# ---------------------------------------------------------------------------


def _conv(n):
    """n alternating user/assistant messages."""
    return [
        {"role": "user" if i % 2 == 0 else "assistant", "content": f"m{i}"}
        for i in range(n)
    ]


class TestChatGPTPruneHistory:
    def _prune(self):
        append_log = MagicMock()
        fn = _extract_method_as_fn(
            CHATGPT_PATH,
            "ChatMoleculeWindow",
            "_prune_history",
            extra_globals={"append_log": append_log},
        )
        return fn, append_log

    def test_system_prompt_survives_prune(self):
        # Regression: a plain [-10:] slice dropped the system prompt.
        fn, _ = self._prune()
        history = [{"role": "system", "content": "SYS"}] + _conv(24)
        self_stub = _Stub(chat_history_state=history, append_message=MagicMock())
        fn(self_stub)
        assert self_stub.chat_history_state[0] == {"role": "system", "content": "SYS"}
        assert len(self_stub.chat_history_state) == 11

    def test_keeps_last_ten_conversation_messages(self):
        fn, _ = self._prune()
        history = [{"role": "system", "content": "SYS"}] + _conv(24)
        self_stub = _Stub(chat_history_state=history, append_message=MagicMock())
        fn(self_stub)
        assert self_stub.chat_history_state[1:] == _conv(24)[-10:]

    def test_no_system_message_keeps_last_ten(self):
        fn, _ = self._prune()
        self_stub = _Stub(chat_history_state=_conv(25), append_message=MagicMock())
        fn(self_stub)
        assert self_stub.chat_history_state == _conv(25)[-10:]

    def test_user_notified_and_logged(self):
        fn, append_log = self._prune()
        self_stub = _Stub(chat_history_state=_conv(25), append_message=MagicMock())
        fn(self_stub)
        self_stub.append_message.assert_called_once()
        append_log.assert_called_once_with("System", "History pruned.")

    def test_error_is_swallowed_and_logged(self):
        fn, append_log = self._prune()
        self_stub = _Stub(append_message=MagicMock())  # no chat_history_state
        fn(self_stub)  # must not raise
        assert append_log.call_args[0][0] == "Error"


class TestLocalPruneHistory:
    def _prune(self):
        append_log = MagicMock()
        fn = _extract_method_as_fn(
            LOCAL_PATH,
            "ChatMoleculeWindow",
            "_prune_history",
            extra_globals={"append_log": append_log},
        )
        return fn, append_log

    def test_prunes_chat_history_state(self):
        # Regression: the old code operated on a never-initialized Gemini
        # chat_session via genai.GenerativeModel — history was never pruned.
        fn, _ = self._prune()
        history = [{"role": "system", "content": "SYS"}] + _conv(24)
        self_stub = _Stub(chat_history_state=history, append_message=MagicMock())
        fn(self_stub)
        assert len(self_stub.chat_history_state) == 11
        assert self_stub.chat_history_state[0]["role"] == "system"
        assert self_stub.chat_history_state[1:] == _conv(24)[-10:]

    def test_no_genai_reference_needed(self):
        # The extracted function must run with no genai in its globals.
        fn, append_log = self._prune()
        self_stub = _Stub(chat_history_state=_conv(25), append_message=MagicMock())
        fn(self_stub)
        assert self_stub.chat_history_state == _conv(25)[-10:]
        assert append_log.call_args[0][0] != "Error"


class TestLocalLogUsage:
    def _log_usage(self, max_history=20):
        append_log = MagicMock()
        prune = MagicMock()
        fn = _extract_method_as_fn(
            LOCAL_PATH,
            "ChatMoleculeWindow",
            "log_usage",
            extra_globals={"append_log": append_log, "MAX_HISTORY": max_history},
        )
        return fn, append_log, prune

    def test_prunes_when_history_exceeds_max(self):
        # Regression: old code checked self.chat_session (never set) so the
        # AttributeError was swallowed and pruning never happened.
        fn, _, prune = self._log_usage()
        self_stub = _Stub(chat_history_state=_conv(21), _prune_history=prune)
        fn(self_stub, _Stub())
        prune.assert_called_once()

    def test_no_prune_at_exactly_max(self):
        fn, _, prune = self._log_usage()
        self_stub = _Stub(chat_history_state=_conv(20), _prune_history=prune)
        fn(self_stub, _Stub())
        prune.assert_not_called()

    def test_openai_usage_attribute_logged(self):
        # Regression: old code only looked at Gemini's usage_metadata.
        fn, append_log, prune = self._log_usage()
        self_stub = _Stub(chat_history_state=[], _prune_history=prune)
        fn(self_stub, _Stub(usage="TOKENS=42"))
        append_log.assert_called_once_with("usage", "TOKENS=42")

    def test_gemini_style_usage_metadata_still_logged(self):
        fn, append_log, prune = self._log_usage()
        self_stub = _Stub(chat_history_state=[], _prune_history=prune)
        fn(self_stub, _Stub(usage_metadata="META"))
        append_log.assert_called_once_with("usage", "META")

    def test_no_usage_attribute_logs_nothing(self):
        fn, append_log, prune = self._log_usage()
        self_stub = _Stub(chat_history_state=[], _prune_history=prune)
        fn(self_stub, _Stub())
        append_log.assert_not_called()


class TestGeminiPruneHistory:
    def _prune(self):
        append_log = MagicMock()
        genai = MagicMock()
        types_mock = MagicMock()
        fn = _extract_method_as_fn(
            GEMINI_PATH,
            "ChatMoleculeWindow",
            "_prune_history",
            extra_globals={
                "append_log": append_log,
                "genai": genai,
                "types": types_mock,
                "SYSTEM_PROMPT": "SYS",
                "GENERATION_CONFIG": {
                    "temperature": 0.2,
                    "top_p": 0.9,
                    "top_k": 40,
                    "max_output_tokens": 4096,
                },
            },
        )
        return fn, append_log, genai

    def test_session_recreated_with_last_ten_messages(self):
        fn, _, _ = self._prune()
        client = MagicMock()
        old_session = _Stub(history=list(range(25)))
        self_stub = _Stub(
            chat_session=old_session,
            settings={"model": "gemini-test"},
            client=client,
            append_message=MagicMock(),
        )
        fn(self_stub)
        _, kwargs = client.chats.create.call_args
        assert kwargs["history"] == list(range(25))[-10:]
        assert kwargs["model"] == "gemini-test"
        assert self_stub.chat_session is client.chats.create.return_value

    def test_create_failure_logged_not_raised(self):
        fn, append_log, _ = self._prune()
        client = MagicMock()
        client.chats.create.side_effect = RuntimeError("boom")
        self_stub = _Stub(
            chat_session=_Stub(history=list(range(25))),
            settings={"model": "m"},
            client=client,
            append_message=MagicMock(),
        )
        fn(self_stub)  # must not raise
        assert append_log.call_args[0][0] == "Error"
