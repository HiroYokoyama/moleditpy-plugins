"""
Plugin-specific tests for Chat with Molecule Neo (Local) that are NOT
already covered by the shared suites (test_shared_chat_variants.py,
test_shared_ai_tool_parsing.py, test_shared_chat_optimizer.py):

- InitWorker.run(): no model-name filtering (unlike the ChatGPT variant),
  passes both api_key and base_url to the OpenAI-compatible client.
- fetch_models(): falls back to the dummy "lm-studio" API key when the field
  is empty (local servers usually don't require one) instead of blocking
  like the ChatGPT variant does.
- save_settings(): persists {api_key, api_base, model, disable_pubchem} and
  additionally refreshes the privacy warning label and Save-button state —
  none of which exist on the ChatGPT/Gemini variants.
- _get_privacy_details(): Local-only privacy-notice classification based on
  the configured API base URL (loopback / local-network / remote).
- check_settings_changed(): Local-only Save-button enable/label/style logic.
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
    / "Chat_with_Molecule_Neo_Local"
    / "chat_with_molecule_neo_local.py"
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
                    # Restore any global shadowed by the `def` itself (e.g.
                    # a module-level function with the same name as the
                    # extracted method) now that `fn` holds the real ref.
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
# InitWorker.run() — no filtering, base_url wired through
# ---------------------------------------------------------------------------


class TestInitWorkerRun:
    def _run(self):
        openai_mock = MagicMock()
        fn = _extract_method_as_fn("InitWorker", "run", extra_globals={"openai": openai_mock})
        return fn, openai_mock

    def test_all_models_returned_unfiltered(self):
        fn, openai_mock = self._run()
        client = openai_mock.OpenAI.return_value
        client.models.list.return_value = [
            MagicMock(id="llama3"),
            MagicMock(id="text-embedding-ada"),
            MagicMock(id="mistral-7b"),
        ]
        self_stub = _Stub(api_key="lm-studio", base_url="http://localhost:1234/v1", finished=MagicMock())
        fn(self_stub)
        models, err = self_stub.finished.emit.call_args[0]
        assert models == ["llama3", "mistral-7b", "text-embedding-ada"]
        assert err == ""

    def test_passes_api_key_and_base_url(self):
        fn, openai_mock = self._run()
        client = openai_mock.OpenAI.return_value
        client.models.list.return_value = []
        self_stub = _Stub(api_key="k", base_url="http://localhost:11434/v1", finished=MagicMock())
        fn(self_stub)
        openai_mock.OpenAI.assert_called_once_with(api_key="k", base_url="http://localhost:11434/v1")

    def test_exception_emits_empty_list_and_message(self):
        fn, openai_mock = self._run()
        openai_mock.OpenAI.side_effect = RuntimeError("connection refused")
        self_stub = _Stub(api_key="k", base_url="http://x", finished=MagicMock())
        fn(self_stub)  # must not raise
        models, err = self_stub.finished.emit.call_args[0]
        assert models == []
        assert err == "connection refused"


# ---------------------------------------------------------------------------
# fetch_models() — empty key falls back to "lm-studio" dummy
# ---------------------------------------------------------------------------


class TestFetchModels:
    def _run(self, api_key_text, api_base_text="http://localhost:1234/v1", demo_mode=False):
        init_worker_cls = MagicMock()
        fn = _extract_method_as_fn(
            "ChatMoleculeWindow",
            "fetch_models",
            extra_globals={"DEMO_MODE": demo_mode, "InitWorker": init_worker_cls},
        )
        self_stub = _Stub(
            txt_api_key=_Text(api_key_text),
            txt_api_base=_Text(api_base_text),
            combo_model=MagicMock(),
            append_message=MagicMock(),
            on_models_fetched=MagicMock(),
        )
        return fn, self_stub, init_worker_cls

    def test_empty_api_key_uses_lm_studio_dummy(self):
        fn, self_stub, init_worker_cls = self._run("")
        fn(self_stub)
        init_worker_cls.assert_called_once_with("lm-studio", "http://localhost:1234/v1")

    def test_non_empty_api_key_used_as_is(self):
        fn, self_stub, init_worker_cls = self._run("real-key")
        fn(self_stub)
        init_worker_cls.assert_called_once_with("real-key", "http://localhost:1234/v1")

    def test_demo_mode_short_circuits(self):
        fn, self_stub, init_worker_cls = self._run("", demo_mode=True)
        fn(self_stub)
        self_stub.on_models_fetched.assert_called_once_with(["demo-model"], None)
        init_worker_cls.assert_not_called()

    def test_worker_started_and_stored(self):
        fn, self_stub, init_worker_cls = self._run("real-key")
        worker = init_worker_cls.return_value
        fn(self_stub)
        worker.finished.connect.assert_called_once_with(self_stub.on_models_fetched)
        worker.start.assert_called_once()
        assert self_stub.temp_worker is worker


# ---------------------------------------------------------------------------
# save_settings() — persists api_base + disable_pubchem, updates warning UI
# ---------------------------------------------------------------------------


class TestSaveSettings:
    def _run(self, checked=True):
        save_settings_mock = MagicMock()
        fn = _extract_method_as_fn(
            "ChatMoleculeWindow", "save_settings", extra_globals={"save_settings": save_settings_mock}
        )
        combo = MagicMock()
        combo.currentText.return_value = "llama3"
        chk = MagicMock()
        chk.isChecked.return_value = checked
        self_stub = _Stub(
            txt_api_key=_Text("sk-new"),
            txt_api_base=_Text("http://localhost:11434/v1"),
            combo_model=combo,
            chk_enable_pubchem=chk,
            settings={},
            initialize_session=MagicMock(),
            append_message=MagicMock(),
            update_warning_label=MagicMock(),
            _get_privacy_details=MagicMock(return_value=("note", "gray")),
            check_settings_changed=MagicMock(),
        )
        return fn, self_stub, save_settings_mock

    def test_persists_full_settings_dict(self):
        fn, self_stub, save_settings_mock = self._run(checked=True)
        fn(self_stub)
        assert self_stub.settings == {
            "api_key": "sk-new",
            "api_base": "http://localhost:11434/v1",
            "model": "llama3",
            "disable_pubchem": False,  # checked (enabled) -> not disabled
        }
        save_settings_mock.assert_called_once_with(self_stub.settings)

    def test_unchecked_pubchem_box_disables_it(self):
        fn, self_stub, _ = self._run(checked=False)
        fn(self_stub)
        assert self_stub.settings["disable_pubchem"] is True

    def test_refreshes_warning_ui_and_button_state(self):
        fn, self_stub, _ = self._run()
        fn(self_stub)
        self_stub.update_warning_label.assert_called_once()
        self_stub.check_settings_changed.assert_called_once()
        self_stub.initialize_session.assert_called_once()


# ---------------------------------------------------------------------------
# _get_privacy_details() — pure classification logic (Local-only feature)
# ---------------------------------------------------------------------------


class TestGetPrivacyDetails:
    def _run(self, url, disable_pubchem):
        fn = _extract_method_as_fn("ChatMoleculeWindow", "_get_privacy_details")
        self_stub = _Stub(txt_api_base=_Text(url), settings={"disable_pubchem": disable_pubchem})
        return fn(self_stub)

    def test_loopback_pubchem_disabled_is_green(self):
        text, color = self._run("http://localhost:1234/v1", True)
        assert color == "green"
        assert "NOT sent" in text

    def test_loopback_pubchem_enabled_is_orange(self):
        text, color = self._run("http://127.0.0.1:1234/v1", False)
        assert color == "orange"
        assert "sent to PubChem" in text

    def test_local_network_pubchem_disabled_is_blue(self):
        text, color = self._run("http://192.168.1.50:1234/v1", True)
        assert color == "blue"

    def test_local_network_pubchem_enabled_is_orange(self):
        text, color = self._run("http://192.168.1.50:1234/v1", False)
        assert color == "orange"

    def test_empty_url_pubchem_disabled_is_gray(self):
        text, color = self._run("", True)
        assert color == "gray"
        assert "NOT sent" in text

    def test_remote_url_is_red(self):
        text, color = self._run("https://api.example.com/v1", True)
        assert color == "red"
        assert "external server" in text


# ---------------------------------------------------------------------------
# check_settings_changed() — Save button enable/label/style logic
# ---------------------------------------------------------------------------


class TestCheckSettingsChanged:
    def _run(self, current, saved):
        fn = _extract_method_as_fn("ChatMoleculeWindow", "check_settings_changed")
        combo = MagicMock()
        combo.currentText.return_value = current["model"]
        chk = MagicMock()
        chk.isChecked.return_value = current["pubchem_checked"]
        btn = MagicMock()
        self_stub = _Stub(
            txt_api_base=_Text(current["api_base"]),
            txt_api_key=_Text(current["api_key"]),
            combo_model=combo,
            chk_enable_pubchem=chk,
            settings=saved,
            btn_save=btn,
        )
        fn(self_stub)
        return btn

    def test_matching_settings_disables_save_button(self):
        saved = {"api_base": "b", "api_key": "k", "model": "m", "disable_pubchem": False}
        current = {"api_base": "b", "api_key": "k", "model": "m", "pubchem_checked": True}
        btn = self._run(current, saved)
        btn.setEnabled.assert_called_once_with(False)
        btn.setText.assert_called_once_with("Save & Reload")

    def test_changed_api_key_enables_save_button(self):
        saved = {"api_base": "b", "api_key": "old", "model": "m", "disable_pubchem": False}
        current = {"api_base": "b", "api_key": "new", "model": "m", "pubchem_checked": True}
        btn = self._run(current, saved)
        btn.setEnabled.assert_called_once_with(True)
        btn.setText.assert_called_once_with("Save & Reload *")

    def test_changed_pubchem_checkbox_enables_save_button(self):
        saved = {"api_base": "b", "api_key": "k", "model": "m", "disable_pubchem": False}
        current = {"api_base": "b", "api_key": "k", "model": "m", "pubchem_checked": False}
        btn = self._run(current, saved)
        btn.setEnabled.assert_called_once_with(True)
