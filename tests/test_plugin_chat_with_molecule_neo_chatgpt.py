"""
Plugin-specific tests for Chat with Molecule Neo (ChatGPT) that are NOT
already covered by the shared suites (test_shared_chat_variants.py,
test_shared_ai_tool_parsing.py, test_shared_chat_optimizer.py):

- InitWorker.run(): filters fetched models to those containing "gpt" and
  sorts them (ChatGPT-specific — the Local variant does not filter).
- fetch_models(): API-key-required guard, DEMO_MODE short-circuit, worker
  wiring.
- on_models_fetched(): error dialog vs success list population.
- save_settings(): persists exactly {api_key, model} (no api_base/
  disable_pubchem — those are Local-only keys).
"""

from __future__ import annotations

import ast
import textwrap
from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin, mock_optional_imports

PLUGIN_PATH = (
    Path(__file__).resolve().parents[1]
    / "plugins"
    / "Chat_with_Molecule_Neo_ChatGPT"
    / "chat_with_molecule_neo_chatGPT.py"
)

with mock_optional_imports():
    _mod = load_plugin(PLUGIN_PATH)


def _extract_method_as_fn(class_name, method_name, extra_globals=None):
    source = PLUGIN_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == method_name:
                    func_src = ast.get_source_segment(source, item)
                    ns: dict = {"MagicMock": MagicMock, **(extra_globals or {})}
                    exec(textwrap.dedent(func_src), ns)  # noqa: S102
                    fn = ns[method_name]
                    # If a global with the same name as the method itself was
                    # supplied (e.g. module-level save_settings() vs the
                    # method save_settings()), the `def` statement above
                    # rebound that name in `ns` to the new function, shadowing
                    # the intended global. Restore it now that `fn` already
                    # holds the real reference; `ns` is `fn.__globals__`.
                    ns.update(extra_globals or {})
                    return fn
    raise AssertionError(f"{class_name}.{method_name} not found")


class _Stub:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class _Text:
    """Stand-in for a QLineEdit: .text() returns a fixed string."""

    def __init__(self, s):
        self._s = s

    def text(self):
        return self._s


# ---------------------------------------------------------------------------
# InitWorker.run() — gpt-only filtering, sorted, error handling
# ---------------------------------------------------------------------------


class TestInitWorkerRun:
    def _run(self):
        openai_mock = MagicMock()
        fn = _extract_method_as_fn("InitWorker", "run", extra_globals={"openai": openai_mock})
        return fn, openai_mock

    def test_filters_to_gpt_models_only(self):
        fn, openai_mock = self._run()
        client = openai_mock.OpenAI.return_value
        client.models.list.return_value = [
            MagicMock(id="gpt-4o"),
            MagicMock(id="text-embedding-3-small"),
            MagicMock(id="gpt-3.5-turbo"),
            MagicMock(id="whisper-1"),
        ]
        self_stub = _Stub(api_key="sk-test", finished=MagicMock())
        fn(self_stub)
        models, err = self_stub.finished.emit.call_args[0]
        assert models == ["gpt-3.5-turbo", "gpt-4o"]
        assert err == ""

    def test_passes_api_key_only_no_base_url(self):
        fn, openai_mock = self._run()
        client = openai_mock.OpenAI.return_value
        client.models.list.return_value = []
        self_stub = _Stub(api_key="sk-abc", finished=MagicMock())
        fn(self_stub)
        openai_mock.OpenAI.assert_called_once_with(api_key="sk-abc")

    def test_exception_emits_empty_list_and_message(self):
        fn, openai_mock = self._run()
        openai_mock.OpenAI.side_effect = RuntimeError("network down")
        self_stub = _Stub(api_key="sk-test", finished=MagicMock())
        fn(self_stub)  # must not raise
        models, err = self_stub.finished.emit.call_args[0]
        assert models == []
        assert err == "network down"

    def test_no_gpt_models_returns_empty_list(self):
        fn, openai_mock = self._run()
        client = openai_mock.OpenAI.return_value
        client.models.list.return_value = [MagicMock(id="text-embedding-3-small")]
        self_stub = _Stub(api_key="sk-test", finished=MagicMock())
        fn(self_stub)
        models, err = self_stub.finished.emit.call_args[0]
        assert models == []
        assert err == ""


# ---------------------------------------------------------------------------
# fetch_models()
# ---------------------------------------------------------------------------


class TestFetchModels:
    def _run(self, api_key_text, demo_mode=False):
        qmb = MagicMock()
        init_worker_cls = MagicMock()
        fn = _extract_method_as_fn(
            "ChatMoleculeWindow",
            "fetch_models",
            extra_globals={
                "DEMO_MODE": demo_mode,
                "QMessageBox": qmb,
                "InitWorker": init_worker_cls,
            },
        )
        self_stub = _Stub(
            txt_api_key=_Text(api_key_text),
            combo_model=MagicMock(),
            append_message=MagicMock(),
            on_models_fetched=MagicMock(),
        )
        return fn, self_stub, qmb, init_worker_cls

    def test_empty_api_key_warns_and_does_not_start_worker(self):
        fn, self_stub, qmb, init_worker_cls = self._run("")
        fn(self_stub)
        qmb.warning.assert_called_once()
        init_worker_cls.assert_not_called()

    def test_demo_mode_short_circuits(self):
        fn, self_stub, qmb, init_worker_cls = self._run("", demo_mode=True)
        fn(self_stub)
        self_stub.on_models_fetched.assert_called_once_with(["demo-model"], None)
        init_worker_cls.assert_not_called()

    def test_valid_key_starts_worker(self):
        fn, self_stub, qmb, init_worker_cls = self._run("sk-real")
        worker = init_worker_cls.return_value
        fn(self_stub)
        init_worker_cls.assert_called_once_with("sk-real")
        worker.finished.connect.assert_called_once_with(self_stub.on_models_fetched)
        worker.start.assert_called_once()
        assert self_stub.temp_worker is worker


# ---------------------------------------------------------------------------
# on_models_fetched()
# ---------------------------------------------------------------------------


class TestOnModelsFetched:
    def _run(self):
        qmb = MagicMock()
        fn = _extract_method_as_fn(
            "ChatMoleculeWindow", "on_models_fetched", extra_globals={"QMessageBox": qmb}
        )
        return fn, qmb

    def test_error_shows_critical_dialog(self):
        fn, qmb = self._run()
        self_stub = _Stub(
            combo_model=MagicMock(),
            settings={},
            append_message=MagicMock(),
        )
        fn(self_stub, [], "boom")
        qmb.critical.assert_called_once()
        self_stub.append_message.assert_not_called()

    def test_success_populates_combo_and_sets_current(self):
        fn, qmb = self._run()
        combo = MagicMock()
        self_stub = _Stub(
            combo_model=combo,
            settings={"model": "gpt-4o"},
            append_message=MagicMock(),
        )
        fn(self_stub, ["gpt-3.5-turbo", "gpt-4o"], "")
        combo.clear.assert_called_once()
        combo.addItems.assert_called_once_with(["gpt-3.5-turbo", "gpt-4o"])
        combo.setCurrentText.assert_called_once_with("gpt-4o")
        self_stub.append_message.assert_called_once()

    def test_current_model_not_in_list_leaves_selection_unset(self):
        fn, qmb = self._run()
        combo = MagicMock()
        self_stub = _Stub(
            combo_model=combo,
            settings={"model": "gpt-9000"},
            append_message=MagicMock(),
        )
        fn(self_stub, ["gpt-4o"], "")
        combo.setCurrentText.assert_not_called()


# ---------------------------------------------------------------------------
# save_settings() — persists exactly {api_key, model}
# ---------------------------------------------------------------------------


class TestSaveSettings:
    def test_persists_api_key_and_model_only(self):
        save_settings_mock = MagicMock()
        fn = _extract_method_as_fn(
            "ChatMoleculeWindow", "save_settings", extra_globals={"save_settings": save_settings_mock}
        )
        combo = MagicMock()
        combo.currentText.return_value = "gpt-4o"
        self_stub = _Stub(
            txt_api_key=_Text("sk-new"),
            combo_model=combo,
            settings={},
            initialize_session=MagicMock(),
            append_message=MagicMock(),
        )
        fn(self_stub)
        assert self_stub.settings == {"api_key": "sk-new", "model": "gpt-4o"}
        assert "api_base" not in self_stub.settings
        assert "disable_pubchem" not in self_stub.settings
        save_settings_mock.assert_called_once_with(self_stub.settings)
        self_stub.initialize_session.assert_called_once()
