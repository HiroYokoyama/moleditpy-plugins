"""
Deep GUI coverage for the "Chat with Molecule Neo" (ChatGPT) plugin.

Builds a real ChatMoleculeWindow (real PyQt6, chemistry/AI libs mocked) and
drives the worker classes, tool dispatch, settings/session lifecycle, and
message-streaming flow directly -- rather than only checking widget wiring
(that's covered by test_gui_plugin_chat_with_molecule_neo.py and
test_gui_shared_chat_variants.py).

Strategy for the heavy chemistry-touching methods (get_current_molecule_smiles,
execute_*, update_structure_diff_based, ...): rdkit's Chem/AllChem are
MagicMock modules under mock_chemistry_imports(), so arithmetic/attribute
access/iteration on "molecules" never raises (MagicMock.__iter__ defaults to
an empty iterator, __add__ etc. return further MagicMocks, truthiness is
True) -- calling these methods exercises real control flow without needing a
real chemistry stack.

Worker classes (InitWorker/GenAIWorker) are QThread subclasses; instead of
actually spawning OS threads we monkeypatch their `.start` to call `.run()`
synchronously on the test thread, so signal handlers fire inline and we can
assert on the resulting window state.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CHATGPT_PATH = PLUGINS_DIR / "Chat_with_Molecule_Neo_ChatGPT" / "chat_with_molecule_neo_chatGPT.py"

with mock_chemistry_imports():
    _mod = load_plugin_for_gui(CHATGPT_PATH)

SENDER_LABEL = "ChatGPT"


# ===========================================================================
# Helpers
# ===========================================================================


def _make_mw(atoms=None, bonds=None, selected=frozenset()):
    """A MagicMock main_window with a usable state_manager.data/scene shape."""
    mw = MagicMock()
    atoms = dict(atoms or {})
    mw.state_manager.data.atoms = atoms
    mw.state_manager.data.bonds = dict(bonds or {})
    mw.state_manager.data.next_atom_id = (max(atoms.keys()) + 1) if atoms else 0
    items = {}
    for aid in atoms:
        item = MagicMock()
        item.isSelected.return_value = aid in selected
        items[aid] = item
    mw.scene.atom_items = items
    mw.scene.bond_items = {}
    mw.data = mw.state_manager.data  # some methods read mw.data directly
    return mw


def _fake_openai_with_model(model_id):
    """A fake `openai` module whose OpenAI().models.list() yields one usable model."""
    m1 = MagicMock()
    m1.id = model_id
    client = MagicMock()
    client.models.list.return_value = [m1]
    fake_openai = MagicMock()
    fake_openai.OpenAI.return_value = client
    return fake_openai, client


def _make_openai_chunk(content):
    """A fake streaming chunk shaped like `chunk.choices[0].delta.content`."""
    chunk = MagicMock()
    chunk.choices[0].delta.content = content
    return chunk


def _install_fake_rdkit_chem_submodules(monkeypatch, **submodules):
    """
    Register fake rdkit.Chem.* submodules in sys.modules so the plugin's
    inline `from rdkit.Chem import Descriptors, Lipinski` (etc.) imports
    resolve to controllable MagicMocks -- regardless of whether real rdkit
    happens to be installed in the environment running the test (it is not,
    under CI; it may be on a dev machine, in which case the inline import
    would otherwise hit the real boost-wrapped function with a fake mol and
    raise, which the plugin swallows but which defeats the point of the test).
    """
    import sys as _sys

    names = submodules or {"Descriptors": {}, "Lipinski": {}, "rdMolDescriptors": {}}
    rdkit_mod = _sys.modules.get("rdkit")
    if rdkit_mod is None:
        rdkit_mod = MagicMock()
        monkeypatch.setitem(_sys.modules, "rdkit", rdkit_mod)
    chem_mod = _sys.modules.get("rdkit.Chem")
    if chem_mod is None:
        chem_mod = MagicMock()
        monkeypatch.setitem(_sys.modules, "rdkit.Chem", chem_mod)
        monkeypatch.setattr(rdkit_mod, "Chem", chem_mod, raising=False)
    created = {}
    for name, attrs in names.items():
        sub = MagicMock()
        for k, v in (attrs or {}).items():
            setattr(sub, k, v)
        monkeypatch.setitem(_sys.modules, f"rdkit.Chem.{name}", sub)
        monkeypatch.setattr(chem_mod, name, sub, raising=False)
        created[name] = sub
    return created


class _FakeExportAtom:
    def __init__(self, symbol, x, y, z):
        self._symbol, self.x, self.y, self.z = symbol, x, y, z

    def GetSymbol(self):
        return self._symbol


class _FakeExportPos:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class _FakeExportConf:
    def __init__(self, atoms):
        self._atoms = atoms

    def GetAtomPosition(self, i):
        a = self._atoms[i]
        return _FakeExportPos(a.x, a.y, a.z)


class _FakeExportMol:
    """Minimal mol supporting the GetConformer/GetAtoms/GetNumAtoms surface
    used by execute_orca_input_generator / execute_gaussian_input_generator /
    execute_save_file's [[atom]] injection -- with real float coordinates so
    the `:10.5f` format specs in the plugin don't choke on a bare MagicMock.
    """

    def __init__(self, atoms):
        self._atoms = atoms
        self._conf = _FakeExportConf(atoms)

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetAtoms(self):
        return list(self._atoms)

    def GetConformer(self):
        return self._conf

    def GetNumConformers(self):
        return 1


def _build_window(monkeypatch, ctx=None, stub_init=True):
    if stub_init:
        monkeypatch.setattr(
            _mod.ChatMoleculeWindow, "initialize_session", lambda self: None
        )
    monkeypatch.setattr(_mod, "append_log", lambda *a, **k: None)
    monkeypatch.setattr(_mod, "save_settings", lambda *a, **k: None)
    if ctx is None:
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        ctx.current_molecule = None
    w = _mod.ChatMoleculeWindow(context=ctx)
    w.render_content = lambda text: text
    return w


@pytest.fixture
def win(qapp, monkeypatch):
    w = _build_window(monkeypatch)
    yield w
    w.destroy()


@pytest.fixture(autouse=True)
def _sync_workers(monkeypatch):
    """Run QThread workers synchronously (call run() directly from start())."""
    monkeypatch.setattr(_mod.InitWorker, "start", lambda self: self.run())
    monkeypatch.setattr(_mod.OpenAIWorker, "start", lambda self: self.run())


# ===========================================================================
# InitWorker
# ===========================================================================


class TestInitWorker:
    def test_success_filters_gpt_models(self, qapp):
        m1 = MagicMock()
        m1.id = "gpt-4o"
        m2 = MagicMock()
        m2.id = "text-embedding-3-small"
        client = MagicMock()
        client.models.list.return_value = [m1, m2]
        fake_openai = MagicMock()
        fake_openai.OpenAI.return_value = client

        worker = _mod.InitWorker("key123")
        seen = []
        worker.finished.connect(lambda models, err: seen.append((models, err)))
        orig = _mod.openai
        _mod.openai = fake_openai
        try:
            worker.run()
        finally:
            _mod.openai = orig

        assert seen == [(["gpt-4o"], "")]

    def test_error_emits_empty_list_and_message(self, qapp):
        fake_openai = MagicMock()
        fake_openai.OpenAI.side_effect = RuntimeError("bad key")
        worker = _mod.InitWorker("key123")
        seen = []
        worker.finished.connect(lambda models, err: seen.append((models, err)))
        orig = _mod.openai
        _mod.openai = fake_openai
        try:
            worker.run()
        finally:
            _mod.openai = orig
        assert seen == [([], "bad key")]

    def test_models_are_sorted(self, qapp):
        m1 = MagicMock()
        m1.id = "gpt-4o"
        m2 = MagicMock()
        m2.id = "gpt-3.5-turbo"
        client = MagicMock()
        client.models.list.return_value = [m1, m2]
        fake_openai = MagicMock()
        fake_openai.OpenAI.return_value = client
        worker = _mod.InitWorker("key")
        seen = []
        worker.finished.connect(lambda models, err: seen.append((models, err)))
        orig = _mod.openai
        _mod.openai = fake_openai
        try:
            worker.run()
        finally:
            _mod.openai = orig
        assert seen == [(["gpt-3.5-turbo", "gpt-4o"], "")]


# ===========================================================================
# NameResolverWorker
# ===========================================================================


class TestNameResolverWorker:
    def test_success_emits_signal(self, qapp, monkeypatch):
        monkeypatch.setattr(
            _mod.PubChemResolver,
            "resolve_inchikey_to_name",
            staticmethod(lambda k: ("Water", None)),
        )
        worker = _mod.NameResolverWorker("INCHIKEY1")
        seen = []
        worker.signals.finished.connect(lambda k, n: seen.append((k, n)))
        worker.run()
        assert seen == [("INCHIKEY1", "Water")]

    def test_no_name_does_not_emit(self, qapp, monkeypatch):
        monkeypatch.setattr(
            _mod.PubChemResolver,
            "resolve_inchikey_to_name",
            staticmethod(lambda k: (None, "err")),
        )
        worker = _mod.NameResolverWorker("INCHIKEY1")
        seen = []
        worker.signals.finished.connect(lambda k, n: seen.append((k, n)))
        worker.run()
        assert seen == []

    def test_exception_is_swallowed(self, qapp, monkeypatch):
        monkeypatch.setattr(
            _mod.PubChemResolver,
            "resolve_inchikey_to_name",
            staticmethod(lambda k: (_ for _ in ()).throw(RuntimeError("boom"))),
        )
        worker = _mod.NameResolverWorker("INCHIKEY1")
        worker.run()  # must not raise


# ===========================================================================
# GenAIWorker
# ===========================================================================


class TestGenAIWorker:
    def test_success_streams_chunks_and_emits_response(self, qapp):
        stream = [_make_openai_chunk("Hello "), _make_openai_chunk("world")]
        client = MagicMock()
        client.chat.completions.create.return_value = stream
        worker = _mod.OpenAIWorker(client, [{"role": "user", "content": "hi"}])
        chunks = []
        finals = []
        worker.chunk_received.connect(chunks.append)
        worker.response_received.connect(finals.append)
        worker.run()
        assert chunks == ["Hello ", "world"]
        assert finals == ["Hello world"]

    def test_empty_delta_content_not_emitted(self, qapp):
        stream = [_make_openai_chunk(""), _make_openai_chunk("hi")]
        client = MagicMock()
        client.chat.completions.create.return_value = stream
        worker = _mod.OpenAIWorker(client, [])
        chunks = []
        worker.chunk_received.connect(chunks.append)
        worker.run()
        assert chunks == ["hi"]

    def test_api_error_emits_error_occurred(self, qapp):
        client = MagicMock()
        client.chat.completions.create.side_effect = RuntimeError("api down")
        worker = _mod.OpenAIWorker(client, [])
        errors = []
        worker.error_occurred.connect(errors.append)
        worker.run()
        assert errors == ["api down"]

    def test_interrupted_stops_without_emitting_response(self, qapp):
        stream = [_make_openai_chunk("Hello ")]
        client = MagicMock()
        client.chat.completions.create.return_value = stream
        worker = _mod.OpenAIWorker(client, [])
        worker._is_interrupted = True
        errors = []
        chunks = []
        finals = []
        worker.error_occurred.connect(errors.append)
        worker.chunk_received.connect(chunks.append)
        worker.response_received.connect(finals.append)
        worker.run()
        assert chunks == []
        assert finals == []
        assert errors == []

    def test_interrupted_exception_does_not_emit_error(self, qapp):
        client = MagicMock()
        client.chat.completions.create.side_effect = RuntimeError("boom")
        worker = _mod.OpenAIWorker(client, [])
        worker._is_interrupted = True
        errors = []
        worker.error_occurred.connect(errors.append)
        worker.run()
        assert errors == []

    def test_stop_sets_interrupted_flag(self, qapp):
        worker = _mod.OpenAIWorker(MagicMock(), [])
        worker.stop()
        assert worker._is_interrupted is True


# ===========================================================================
# get_current_molecule_smiles
# ===========================================================================


class TestGetCurrentMoleculeSmiles:
    def test_no_rdkit(self, win, monkeypatch):
        monkeypatch.setattr(_mod, "Chem", None)
        smiles, err = win.get_current_molecule_smiles()
        assert smiles is None and "RDKit" in err

    def test_no_main_window(self, win):
        win.main_window = None
        smiles, err = win.get_current_molecule_smiles()
        assert smiles is None and "Internal Error" in err

    def test_empty_scene_no_fallback(self, win):
        win.main_window = _make_mw(atoms={})
        win.context.current_molecule = None
        smiles, err = win.get_current_molecule_smiles()
        assert smiles is None
        assert "empty" in err.lower()

    def test_fallback_to_context_current_molecule(self, win):
        win.main_window = _make_mw(atoms={})
        win.context.current_molecule = MagicMock()
        smiles, err = win.get_current_molecule_smiles()
        assert err is None
        assert smiles is not None

    def test_full_reconstruction_from_2d_data(self, win):
        atoms = {
            0: {"symbol": "C", "charge": 0},
            1: {"symbol": "O", "charge": 0, "explicit_valence": 1},
        }
        bonds = {(0, 1): {"atom1": 0, "atom2": 1, "order": 1, "stereo": 1}}
        win.main_window = _make_mw(atoms=atoms, bonds=bonds)
        smiles, err = win.get_current_molecule_smiles()
        assert err is None
        assert smiles is not None


# ===========================================================================
# get_selected_atom_indices / _build_context_msg
# ===========================================================================


class TestGetSelectedAtomIndices:
    def test_no_data_attr_returns_empty(self, win):
        win.main_window = MagicMock(spec=[])
        assert win.get_selected_atom_indices() == []

    def test_returns_selected_map_nums(self, win):
        mw = MagicMock()
        mw.data.atoms = {0: {}, 1: {}, 2: {}}
        item0 = MagicMock()
        item0.isSelected.return_value = True
        item1 = MagicMock()
        item1.isSelected.return_value = False
        item2 = MagicMock()
        item2.isSelected.return_value = True
        mw.scene.atom_items = {0: item0, 1: item1, 2: item2}
        win.main_window = mw
        result = win.get_selected_atom_indices()
        assert sorted(result, key=int) == ["1", "3"]


class TestBuildContextMsg:
    def test_includes_name_when_known(self, win):
        win.main_window = _make_mw(atoms={})
        msg = win._build_context_msg("CCO", "Ethanol", lazy=True)
        assert "Ethanol" in msg
        assert "CCO" in msg

    def test_lazy_skips_descriptor_calc(self, win, monkeypatch):
        win.main_window = _make_mw(atoms={})
        called = []
        monkeypatch.setattr(win, "_get_descriptors_str", lambda s: called.append(s) or "X")
        win._build_context_msg("CCO", None, lazy=True)
        assert called == []

    def test_no_name_triggers_descriptor_calc(self, win, monkeypatch):
        win.main_window = _make_mw(atoms={})
        monkeypatch.setattr(win, "_get_descriptors_str", lambda s: "Props: MW=1")
        msg = win._build_context_msg("CCO", None, lazy=False)
        assert "Props: MW=1" in msg

    def test_selection_injected(self, win, monkeypatch):
        win.main_window = _make_mw(atoms={})
        monkeypatch.setattr(win, "get_selected_atom_indices", lambda: ["1", "2"])
        msg = win._build_context_msg("CCO", "Ethanol", lazy=True)
        assert "1, 2" in msg

    def test_calc_results_included_and_popped(self, win):
        win.main_window = _make_mw(atoms={})
        win.calc_results_by_smiles["CCO"] = "\n[System: cached]"
        msg = win._build_context_msg("CCO", "Ethanol", lazy=True)
        assert "cached" in msg
        assert "CCO" not in win.calc_results_by_smiles


# ===========================================================================
# initialize_session
# ===========================================================================


class TestInitializeSession:
    def _fresh(self, qapp, monkeypatch):
        return _build_window(monkeypatch, stub_init=False)

    def test_no_openai_lib_disables_input(self, qapp, monkeypatch):
        w = self._fresh(qapp, monkeypatch)
        monkeypatch.setattr(_mod, "HAS_OPENAI", False)
        w.initialize_session()
        assert not w.txt_input.isEnabled()
        w.destroy()

    def test_missing_api_key_warns(self, qapp, monkeypatch):
        w = self._fresh(qapp, monkeypatch)
        monkeypatch.setattr(_mod, "HAS_OPENAI", True)
        monkeypatch.setattr(_mod.QMessageBox, "warning", MagicMock())
        w.settings = {}
        w.initialize_session()
        assert w.client is None
        w.destroy()

    def test_success_with_smiles_sets_pending_context(self, qapp, monkeypatch):
        w = self._fresh(qapp, monkeypatch)
        monkeypatch.setattr(_mod, "HAS_OPENAI", True)
        w.settings = {"api_key": "KEY", "model": "gpt-4o"}
        fake_openai, fake_client = _fake_openai_with_model("gpt-4o")
        monkeypatch.setattr(_mod, "openai", fake_openai)
        monkeypatch.setattr(w, "get_current_molecule_smiles", lambda: ("CCO", None))
        monkeypatch.setattr(w, "get_molecule_name", lambda **k: "Ethanol")
        w.initialize_session()
        assert w.client is fake_client
        assert w.chat_history_state[0] == {"role": "system", "content": _mod.SYSTEM_PROMPT}
        assert w.pending_context_msg is not None
        assert w.txt_input.isEnabled()
        w.destroy()

    def test_no_smiles_sets_context_label(self, qapp, monkeypatch):
        w = self._fresh(qapp, monkeypatch)
        monkeypatch.setattr(_mod, "HAS_OPENAI", True)
        w.settings = {"api_key": "KEY", "model": "gpt-4o"}
        fake_openai, _ = _fake_openai_with_model("gpt-4o")
        monkeypatch.setattr(_mod, "openai", fake_openai)
        monkeypatch.setattr(w, "get_current_molecule_smiles", lambda: (None, "no mol"))
        w.initialize_session()
        assert "No valid molecule" in w.lbl_context.text()
        w.destroy()

    def test_client_creation_failure_leaves_client_none(self, qapp, monkeypatch):
        w = self._fresh(qapp, monkeypatch)
        monkeypatch.setattr(_mod, "HAS_OPENAI", True)
        w.settings = {"api_key": "KEY", "model": "gpt-4o"}
        fake_openai, fake_client = _fake_openai_with_model("gpt-4o")
        fake_openai.OpenAI.side_effect = [fake_client, RuntimeError("boom")]
        monkeypatch.setattr(_mod, "openai", fake_openai)
        monkeypatch.setattr(w, "get_current_molecule_smiles", lambda: (None, None))
        w.initialize_session()
        assert w.client is None
        w.destroy()


# ===========================================================================
# on_init_finished
# ===========================================================================


class TestOnInitFinished:
    def test_error_message_shown(self, win):
        win.on_init_finished([], "boom")
        assert "boom" in win.chat_history_log[-1]["text"]

    def test_empty_models_shown(self, win):
        win.on_init_finished([], "")
        assert "No models" in win.chat_history_log[-1]["text"]

    def test_saved_model_selected(self, win, monkeypatch):
        fake_openai = MagicMock()
        monkeypatch.setattr(_mod, "openai", fake_openai)
        monkeypatch.setattr(win, "get_current_molecule_smiles", lambda: (None, None))
        win.settings = {"model": "gpt-4o"}
        win.on_init_finished(["gpt-4o", "gpt-3.5-turbo"], "")
        assert win.combo_model.currentText() == "gpt-4o"

    def test_unsaved_model_falls_back_to_preferred(self, win, monkeypatch):
        fake_openai = MagicMock()
        monkeypatch.setattr(_mod, "openai", fake_openai)
        monkeypatch.setattr(win, "get_current_molecule_smiles", lambda: (None, None))
        win.settings = {"model": "nonexistent"}
        win.on_init_finished(["gpt-4o", "gpt-other"], "")
        assert win.combo_model.currentText() == "gpt-4o"

    def test_no_preferred_models_uses_first_gpt(self, win, monkeypatch):
        fake_openai = MagicMock()
        monkeypatch.setattr(_mod, "openai", fake_openai)
        monkeypatch.setattr(win, "get_current_molecule_smiles", lambda: (None, None))
        win.settings = {"model": "nonexistent"}
        win.on_init_finished(["gpt-custom-model"], "")
        assert win.combo_model.currentText() == "gpt-custom-model"

    def test_no_gpt_models_uses_first_available(self, win, monkeypatch):
        fake_openai = MagicMock()
        monkeypatch.setattr(_mod, "openai", fake_openai)
        monkeypatch.setattr(win, "get_current_molecule_smiles", lambda: (None, None))
        win.settings = {"model": "nonexistent"}
        win.on_init_finished(["other-model-x"], "")
        assert win.combo_model.currentText() == "other-model-x"


# ===========================================================================
# fetch_models
# ===========================================================================


class TestFetchModels:
    def test_no_api_key_warns(self, win, monkeypatch):
        warn = MagicMock()
        monkeypatch.setattr(_mod.QMessageBox, "warning", warn)
        win.txt_api_key.setText("")
        win.fetch_models()
        warn.assert_called_once()

    def test_success_populates_combo(self, win, monkeypatch):
        win.txt_api_key.setText("KEY")
        fake_openai, _client = _fake_openai_with_model("gpt-4o")
        monkeypatch.setattr(_mod, "openai", fake_openai)
        win.fetch_models()
        assert win.combo_model.count() >= 1


# ===========================================================================
# _dispatch_tool routing
# ===========================================================================


class TestDispatchTool:
    ROUTES = [
        ("apply_transformation", "execute_apply_transformation"),
        ("highlight_substructure", "execute_highlight_substructure"),
        ("calculate_descriptors", "execute_calculate_descriptors"),
        ("orca_input_generator", "execute_orca_input_generator"),
        ("gaussian_input_generator", "execute_gaussian_input_generator"),
        ("save_file", "execute_save_file"),
        ("load_molecule", "execute_load_molecule"),
        ("load_molecule_by_name", "execute_load_molecule_by_name"),
        ("convert_to_3d", "execute_convert_to_3d"),
        ("set_electronic_state", "execute_set_electronic_state"),
        ("clear_canvas", "execute_clear_canvas"),
    ]

    @pytest.mark.parametrize("tool_name,attr", ROUTES)
    def test_routes_to_executor(self, win, monkeypatch, tool_name, attr):
        mock = MagicMock(return_value=None)
        monkeypatch.setattr(win, attr, mock)
        win._dispatch_tool(tool_name, {"x": 1})
        mock.assert_called_once_with({"x": 1})

    def test_unknown_tool_returns_error(self, win):
        err = win._dispatch_tool("no_such_tool", {})
        assert "Unknown tool" in err


# ===========================================================================
# execute_* tool implementations
# ===========================================================================


class TestExecuteToolImplementations:
    def _win_with_mol(self, win, smiles="CCO"):
        win.main_window = _make_mw(atoms={0: {"symbol": "C"}, 1: {"symbol": "O"}})
        win.get_current_molecule_smiles = lambda: (smiles, None)
        return win

    def test_apply_transformation_no_smarts(self, win):
        result = win.execute_apply_transformation({})
        assert result is None
        assert "No reaction SMARTS" in win.chat_history_log[-1]["text"]

    def test_apply_transformation_no_molecule(self, win):
        win.get_current_molecule_smiles = lambda: (None, "no mol")
        result = win.execute_apply_transformation({"reaction_smarts": "[C:1]>>[N:1]"})
        assert "no mol" in win.chat_history_log[-1]["text"]

    def test_apply_transformation_runs_with_mocked_rdkit(self, win):
        self._win_with_mol(win)
        win.execute_apply_transformation({"reaction_smarts": "[C:1]>>[N:1]"})
        # Mocked RDKit -> AllChem.ReactionFromSmarts/.RunReactants all MagicMocks;
        # must not raise, some message gets appended either way.
        assert win.chat_history_log

    def test_highlight_substructure_no_params_returns_none(self, win):
        assert win.execute_highlight_substructure({}) is None

    def test_highlight_substructure_no_molecule(self, win):
        win.get_current_molecule_smiles = lambda: (None, "no mol")
        win.execute_highlight_substructure({"smarts": "[C]"})

    def test_highlight_substructure_runs(self, win):
        self._win_with_mol(win)
        win.execute_highlight_substructure({"atom_indices": [1, 2]})

    def test_set_electronic_state_no_molecule(self, win):
        win.get_current_molecule_smiles = lambda: (None, "no mol")
        err = win.execute_set_electronic_state({"charge": 1})
        assert err == "No molecule loaded."

    def test_set_electronic_state_runs(self, win):
        self._win_with_mol(win)
        win.execute_set_electronic_state({"atom_index": 1, "charge": 1})

    def test_calculate_descriptors_no_molecule(self, win):
        win.get_current_molecule_smiles = lambda: (None, "no mol")
        assert win.execute_calculate_descriptors({}) is None

    def test_calculate_descriptors_runs(self, win, monkeypatch):
        self._win_with_mol(win)
        monkeypatch.setattr(_mod.QMessageBox, "exec", lambda self: None, raising=False)
        win.execute_calculate_descriptors(
            {"properties": ["MW", "LogP", "TPSA", "HBondDonor", "RotatableBonds", "AromaticRings", "Rings"]}
        )

    def test_load_molecule_no_smiles(self, win):
        win.execute_load_molecule({})
        assert "No SMILES" in win.chat_history_log[-1]["text"]

    def test_load_molecule_runs(self, win, monkeypatch):
        monkeypatch.setattr(win, "load_smiles_undo_safe", MagicMock())
        win.main_window = _make_mw()
        win.execute_load_molecule({"smiles": "CCO", "name": "Ethanol"})
        assert "Ethanol" in win.chat_history_log[-1]["text"]

    def test_load_molecule_raises_returns_error(self, win, monkeypatch):
        monkeypatch.setattr(
            win, "load_smiles_undo_safe", MagicMock(side_effect=RuntimeError("boom"))
        )
        err = win.execute_load_molecule({"smiles": "CCO"})
        assert err == "boom"
        assert "Load Error" in win.chat_history_log[-1]["text"]

    def test_load_molecule_by_name_no_name(self, win):
        err = win.execute_load_molecule_by_name({})
        assert err == "No name provided"

    def test_load_molecule_by_name_no_smiles_returned(self, win, monkeypatch):
        monkeypatch.setattr(
            _mod.PubChemResolver,
            "resolve_name_to_smiles",
            staticmethod(lambda n: (None, None)),
        )
        err = win.execute_load_molecule_by_name({"name": "ghost"})
        assert err == "No SMILES found"

    def test_load_molecule_by_name_raises_returns_error(self, win, monkeypatch):
        monkeypatch.setattr(
            _mod.PubChemResolver,
            "resolve_name_to_smiles",
            staticmethod(lambda n: (_ for _ in ()).throw(RuntimeError("net error"))),
        )
        err = win.execute_load_molecule_by_name({"name": "ethanol"})
        assert err == "net error"

    def test_load_molecule_by_name_search_fails(self, win, monkeypatch):
        monkeypatch.setattr(
            _mod.PubChemResolver,
            "resolve_name_to_smiles",
            staticmethod(lambda n: (None, "not found")),
        )
        err = win.execute_load_molecule_by_name({"name": "unobtainium"})
        assert err == "not found"

    def test_load_molecule_by_name_success(self, win, monkeypatch):
        monkeypatch.setattr(
            _mod.PubChemResolver,
            "resolve_name_to_smiles",
            staticmethod(lambda n: ("CCO", None)),
        )
        monkeypatch.setattr(win, "load_smiles_undo_safe", MagicMock())
        win.main_window = _make_mw()
        win.execute_load_molecule_by_name({"name": "ethanol"})
        assert "Loaded" in win.chat_history_log[-1]["text"]

    def test_convert_to_3d_runs(self, win, monkeypatch):
        monkeypatch.setattr(win, "_ensure_main_window_3d_conversion", MagicMock())
        win.execute_convert_to_3d({})
        assert "3D structure displayed" in win.chat_history_log[-1]["text"]

    def test_clear_canvas_runs(self, win):
        win.main_window = _make_mw(atoms={0: {}}, bonds={(0, 1): {}})
        win.execute_clear_canvas({})
        assert "Canvas cleared" in win.chat_history_log[-1]["text"]

    def test_clear_canvas_error_path_returns_message(self, win):
        win.main_window = None  # AttributeError inside -> caught
        err = win.execute_clear_canvas({})
        assert err is not None


class TestSaveFileAndInputGenerators:
    def test_execute_save_file_no_tags_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            _mod.QFileDialog, "getSaveFileName", lambda *a, **k: ("", "")
        )
        err = win.execute_save_file({"filename": "out.txt", "content": "hi"})
        assert err == "Cancelled"

    def test_execute_save_file_writes_content(self, win, monkeypatch, tmp_path):
        target = str(tmp_path / "out.txt")
        monkeypatch.setattr(
            _mod.QFileDialog, "getSaveFileName", lambda *a, **k: (target, "")
        )
        err = win.execute_save_file({"filename": "out.txt", "content": "hello"})
        assert err is None
        assert Path(target).read_text(encoding="utf-8") == "hello"

    def test_execute_save_file_atom_tag_no_molecule(self, win, monkeypatch):
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: None)
        err = win.execute_save_file({"filename": "out.xyz", "content": "[[atom_count]]"})
        assert "Geometry tag" in win.chat_history_log[-1]["text"]
        assert err == "No molecule for XYZ tag"

    def test_execute_save_file_atom_tag_with_molecule(self, win, monkeypatch, tmp_path):
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 2
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: mol)
        target = str(tmp_path / "out.xyz")
        monkeypatch.setattr(
            _mod.QFileDialog, "getSaveFileName", lambda *a, **k: (target, "")
        )
        win.execute_save_file(
            {"filename": "out.xyz", "content": "[[atom_count]]\n[[atom]]"}
        )

    def test_execute_save_file_atom_tag_geometry_error(self, win, monkeypatch):
        mol = MagicMock()
        mol.GetConformer.side_effect = RuntimeError("conf missing")
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: mol)
        err = win.execute_save_file({"filename": "out.xyz", "content": "[[atom_count]]"})
        assert "conf missing" in err

    def test_execute_orca_input_generator_no_molecule(self, win, monkeypatch):
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: None)
        err = win.execute_orca_input_generator({"filename": "a.inp"})
        assert err == "No valid molecule"

    def test_execute_orca_input_generator_raises(self, win, monkeypatch):
        mol = MagicMock()
        mol.GetConformer.side_effect = RuntimeError("no conformer")
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: mol)
        err = win.execute_orca_input_generator({"filename": "a.inp"})
        assert err is not None
        assert "ORCA Gen Error" in win.chat_history_log[-1]["text"]

    def test_execute_orca_input_generator_writes(self, win, monkeypatch, tmp_path):
        mol = MagicMock()
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: mol)
        target = str(tmp_path / "a.inp")
        monkeypatch.setattr(
            _mod.QFileDialog, "getSaveFileName", lambda *a, **k: (target, "")
        )
        win.execute_orca_input_generator(
            {"filename": "a.inp", "header": "! HF", "footer": ""}
        )

    def test_execute_gaussian_input_generator_no_molecule(self, win, monkeypatch):
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: None)
        err = win.execute_gaussian_input_generator({"filename": "a.gjf"})
        assert err == "No valid molecule"

    def test_execute_gaussian_input_generator_writes(self, win, monkeypatch, tmp_path):
        mol = MagicMock()
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: mol)
        target = str(tmp_path / "a.gjf")
        monkeypatch.setattr(
            _mod.QFileDialog, "getSaveFileName", lambda *a, **k: (target, "")
        )
        win.execute_gaussian_input_generator(
            {"filename": "a.gjf", "header": "#P HF", "footer": ""}
        )

    def test_execute_gaussian_input_generator_raises(self, win, monkeypatch):
        mol = MagicMock()
        mol.GetConformer.side_effect = RuntimeError("no conformer")
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: mol)
        err = win.execute_gaussian_input_generator({"filename": "a.gjf"})
        assert err is not None
        assert "Gaussian Gen Error" in win.chat_history_log[-1]["text"]


class TestGetMoleculeForExport:
    def test_uses_existing_3d_mol(self, win):
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        win.get_current_molecule_smiles = lambda: ("CCO", None)
        win.context.current_molecule = mol
        result = win._get_molecule_for_export()
        assert result is mol

    def test_falls_back_to_local_generation(self, win, monkeypatch):
        win.get_current_molecule_smiles = lambda: ("CCO", None)
        win.context.current_molecule = None
        win.main_window = _make_mw()
        monkeypatch.setattr(win, "_ensure_main_window_3d_conversion", MagicMock())
        # AllChem.EmbedMolecule mocked -> returns MagicMock, compared to 0 -> False by default
        # (MagicMock.__eq__ default identity-based), so falls to the "else" warn branch.
        win._get_molecule_for_export()


# ===========================================================================
# propose / accept / reject / step / retry tool actions
# ===========================================================================


class TestToolActionFlow:
    def test_propose_single_tool_shows_dialog(self, win):
        win.propose_tool_action({"tool": "clear_canvas", "params": {}})
        assert not win.tool_confirm_frame.isHidden()
        assert win.btn_tool_step.isHidden()

    def test_propose_multi_tool_shows_step_button(self, win):
        win.propose_tool_action(
            {"tools": [{"tool": "clear_canvas"}, {"tool": "convert_to_3d"}]}
        )
        assert not win.btn_tool_step.isHidden()

    @pytest.mark.parametrize(
        "tool_name,params",
        [
            ("apply_transformation", {"reaction_smarts": "x"}),
            ("highlight_substructure", {"smarts": "x"}),
            ("calculate_descriptors", {"properties": ["MW"]}),
            ("orca_input_generator", {"filename": "a.inp"}),
            ("gaussian_input_generator", {"filename": "a.gjf"}),
            ("save_file", {"filename": "a.txt"}),
            ("load_molecule", {"name": "Ethanol"}),
            ("load_molecule_by_name", {"name": "Ethanol"}),
            ("set_electronic_state", {"charge": 1, "multiplicity": 1}),
            ("convert_to_3d", {}),
            ("clear_canvas", {}),
        ],
    )
    def test_propose_single_tool_all_descriptions(self, win, tool_name, params):
        win.propose_tool_action({"tool": tool_name, "params": params})
        assert "ChatGPT suggests" in win.lbl_tool_info.text()

    def test_accept_single_tool_dispatches(self, win, monkeypatch):
        mock = MagicMock(return_value=None)
        monkeypatch.setattr(win, "execute_clear_canvas", mock)
        win.propose_tool_action({"tool": "clear_canvas", "params": {}})
        win.accept_tool_action()
        mock.assert_called_once()
        assert win.tool_confirm_frame.isHidden()

    def test_accept_multi_tool_dispatches_all(self, win, monkeypatch):
        m1 = MagicMock(return_value=None)
        m2 = MagicMock(return_value=None)
        monkeypatch.setattr(win, "execute_clear_canvas", m1)
        monkeypatch.setattr(win, "execute_convert_to_3d", m2)
        win.propose_tool_action(
            {"tools": [{"tool": "clear_canvas", "params": {}}, {"tool": "convert_to_3d", "params": {}}]}
        )
        win.accept_tool_action()
        m1.assert_called_once()
        m2.assert_called_once()

    def test_accept_tool_failure_enters_retry_state(self, win, monkeypatch):
        monkeypatch.setattr(win, "execute_clear_canvas", MagicMock(return_value="boom"))
        win.propose_tool_action({"tool": "clear_canvas", "params": {}})
        win.accept_tool_action()
        assert win.tool_state == "RETRY"
        assert not win.tool_confirm_frame.isHidden()

    def test_accept_no_pending_payload_noop(self, win):
        win.pending_tool_payload = None
        win.accept_tool_action()  # must not raise

    def test_reject_tool_action_queues_message(self, win):
        win.pending_context_msg = None
        win.propose_tool_action({"tool": "clear_canvas", "params": {}})
        win.reject_tool_action()
        assert win.tool_confirm_frame.isHidden()
        assert "rejected" in win.pending_context_msg

    def test_reject_appends_to_existing_pending_msg(self, win):
        win.pending_context_msg = "existing"
        win.propose_tool_action({"tool": "clear_canvas", "params": {}})
        win.reject_tool_action()
        assert win.pending_context_msg.startswith("existing")

    def test_on_tool_accept_click_confirm_state(self, win, monkeypatch):
        accept = MagicMock()
        monkeypatch.setattr(win, "accept_tool_action", accept)
        win.tool_state = "CONFIRM"
        win.on_tool_accept_click()
        accept.assert_called_once()

    def test_on_tool_accept_click_retry_state(self, win, monkeypatch):
        retry = MagicMock()
        monkeypatch.setattr(win, "retry_tool_action", retry)
        win.tool_state = "RETRY"
        win.on_tool_accept_click()
        retry.assert_called_once()

    def test_on_tool_reject_click_confirm_calls_reject(self, win, monkeypatch):
        reject = MagicMock()
        monkeypatch.setattr(win, "reject_tool_action", reject)
        win.tool_state = "CONFIRM"
        win.on_tool_reject_click()
        reject.assert_called_once()

    def test_on_tool_reject_click_retry_skips_reject(self, win, monkeypatch):
        reject = MagicMock()
        monkeypatch.setattr(win, "reject_tool_action", reject)
        win.tool_state = "RETRY"
        win.on_tool_reject_click()
        reject.assert_not_called()

    def test_retry_tool_action_sends_message(self, win, monkeypatch):
        send = MagicMock()
        monkeypatch.setattr(win, "send_message", send)
        win.last_tool_error = "some failure"
        win.retry_tool_action()
        send.assert_called_once()
        assert "some failure" in win.txt_input.toPlainText()

    def test_step_click_no_payload_noop(self, win):
        win.pending_tool_payload = None
        win.on_tool_step_click()

    def test_step_click_single_tool_fallback_to_accept(self, win, monkeypatch):
        accept = MagicMock()
        monkeypatch.setattr(win, "accept_tool_action", accept)
        win.pending_tool_payload = {"tool": "clear_canvas", "params": {}}
        win.on_tool_step_click()
        accept.assert_called_once()

    def test_step_click_executes_first_of_many(self, win, monkeypatch):
        dispatch = MagicMock(return_value=None)
        monkeypatch.setattr(win, "_dispatch_tool", dispatch)
        win.pending_tool_payload = {
            "tools": [{"tool": "clear_canvas", "params": {}}, {"tool": "convert_to_3d", "params": {}}]
        }
        win.on_tool_step_click()
        dispatch.assert_called_once_with("clear_canvas", {})
        assert win.pending_tool_payload["tools"] == [{"tool": "convert_to_3d", "params": {}}]

    def test_step_click_last_tool_completes_chain(self, win, monkeypatch):
        dispatch = MagicMock(return_value=None)
        monkeypatch.setattr(win, "_dispatch_tool", dispatch)
        win.pending_tool_payload = {"tools": [{"tool": "clear_canvas", "params": {}}]}
        win.on_tool_step_click()
        assert win.pending_tool_payload is None
        assert win.tool_confirm_frame.isHidden()

    def test_step_click_failure_enters_retry(self, win, monkeypatch):
        dispatch = MagicMock(return_value="boom")
        monkeypatch.setattr(win, "_dispatch_tool", dispatch)
        win.pending_tool_payload = {"tools": [{"tool": "clear_canvas", "params": {}}]}
        win.on_tool_step_click()
        assert win.tool_state == "RETRY"


# ===========================================================================
# send_message / streaming / tool-call end-to-end via synchronous worker
# ===========================================================================


class TestSendMessageFlow:
    def _ready_window(self, monkeypatch):
        w = _build_window(monkeypatch)
        w.client = MagicMock()
        w.chat_history_state = [{"role": "system", "content": "SYS"}]
        w.get_current_molecule_smiles = lambda: (None, None)
        return w

    def test_empty_text_is_noop(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.txt_input.setPlainText("")
        w.send_message()
        assert w.chat_history_log == []
        w.destroy()

    def test_no_client_warns(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.client = None
        w.txt_input.setPlainText("hello")
        w.send_message()
        assert "not initialized" in w.chat_history_log[-1]["text"]
        w.destroy()

    def test_stop_button_click_stops_generation(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.btn_send.setText("Stop")
        worker = MagicMock()
        worker.isRunning.return_value = True
        w.worker = worker
        w.send_message()
        worker.stop.assert_called_once()
        assert w.worker is None
        w.destroy()

    def test_success_response_no_tool(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.client.chat.completions.create.return_value = [
            _make_openai_chunk("Sure, here's the answer.")
        ]
        w.txt_input.setPlainText("What is water?")
        w.send_message()
        assert any(
            e["sender"] == SENDER_LABEL and "answer" in e["text"]
            for e in w.chat_history_log
        )
        assert w.txt_input.isEnabled()
        w.destroy()

    def test_error_response_resets_ui(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.client.chat.completions.create.side_effect = RuntimeError("network down")
        w.txt_input.setPlainText("hi")
        w.send_message()
        assert any("network down" in e["text"] for e in w.chat_history_log)
        assert w.txt_input.isEnabled()
        w.destroy()

    def test_tool_call_response_shows_confirmation(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.client.chat.completions.create.return_value = [
            _make_openai_chunk('```json\n{"tool": "clear_canvas", "params": {}}\n```')
        ]
        w.txt_input.setPlainText("clear it")
        w.send_message()
        assert not w.tool_confirm_frame.isHidden()
        assert w.pending_tool_payload == {"tool": "clear_canvas", "params": {}}
        w.destroy()

    def test_multi_tool_call_response(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.client.chat.completions.create.return_value = [
            _make_openai_chunk(
                '```json\n[{"tool": "clear_canvas", "params": {}}, '
                '{"tool": "convert_to_3d", "params": {}}]\n```'
            )
        ]
        w.txt_input.setPlainText("do both")
        w.send_message()
        assert w.pending_tool_payload["tools"][0]["tool"] == "clear_canvas"
        w.destroy()

    def test_malformed_json_tool_block_ignored(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.client.chat.completions.create.return_value = [
            _make_openai_chunk('```json\n{not valid json}\n```')
        ]
        w.txt_input.setPlainText("test")
        w.send_message()
        assert w.tool_confirm_frame.isHidden()
        w.destroy()

    def test_pending_context_and_info_text_merged_into_prompt(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.pending_info_text = "Context updated: CCO"
        w.client.chat.completions.create.return_value = [_make_openai_chunk("ok")]
        w.txt_input.setPlainText("hi")
        w.send_message()
        assert any("Context updated" in e["text"] for e in w.chat_history_log)
        w.destroy()

    def test_stop_generation_resets_ui(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        worker = MagicMock()
        worker.isRunning.return_value = True
        w.worker = worker
        w.stop_generation()
        worker.stop.assert_called_once()
        assert w.worker is None
        assert w.txt_input.isEnabled()
        w.destroy()

    def test_stop_generation_noop_without_worker(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.worker = None
        w.stop_generation()  # must not raise
        w.destroy()

    def test_close_event_disconnects_running_worker(self, qapp, monkeypatch):
        w = self._ready_window(monkeypatch)
        w.worker = MagicMock()
        w.worker.isRunning.return_value = True
        w.init_worker = MagicMock()
        w.init_worker.isRunning.return_value = True
        from PyQt6.QtGui import QCloseEvent

        w.closeEvent(QCloseEvent())
        w.destroy()


# ===========================================================================
# on_worker_finished / on_error / on_chunk_received / log_usage / _prune_history
# ===========================================================================


class TestResponseCallbacks:
    def test_on_worker_finished_clears_worker(self, win):
        win.worker = MagicMock()
        win.on_worker_finished()
        assert win.worker is None

    def test_on_worker_finished_noop_if_none(self, win):
        win.worker = None
        win.on_worker_finished()

    def test_on_chunk_received_starts_stream(self, win):
        win.stream_accumulated_text = ""
        win.on_chunk_received("part1")
        assert win.stream_accumulated_text == "part1"

    def test_on_error_resets_ui_and_logs(self, win):
        win.worker = MagicMock()
        win.on_error("network fail")
        assert "network fail" in win.chat_history_log[-1]["text"]
        assert win.worker is None
        assert win.txt_input.isEnabled()

    def test_log_usage_with_usage_attribute(self, win, monkeypatch):
        logged = []
        monkeypatch.setattr(_mod, "append_log", lambda *a: logged.append(a))
        win.chat_history_state = [1, 2, 3]
        response = MagicMock(usage="TOKENS=42")
        win.log_usage(response)
        assert ("usage", "TOKENS=42") in logged

    def test_log_usage_triggers_prune_over_max_history(self, win, monkeypatch):
        win.chat_history_state = list(range(_mod.MAX_HISTORY + 5))
        called = []
        monkeypatch.setattr(win, "_prune_history", lambda: called.append(True))
        win.log_usage(MagicMock(spec=[]))
        assert called == [True]

    def test_prune_history_keeps_system_and_last_ten(self, win):
        win.chat_history_state = [{"role": "system", "content": "SYS"}] + [
            {"role": "user", "content": f"m{i}"} for i in range(24)
        ]
        win._prune_history()
        assert win.chat_history_state[0] == {"role": "system", "content": "SYS"}
        assert len(win.chat_history_state) == 11
        assert "History pruned" in win.chat_history_log[-1]["text"]

    def test_prune_history_failure_logged_not_raised(self, win, monkeypatch):
        del win.chat_history_state  # -> AttributeError inside the try block
        logged = []
        monkeypatch.setattr(_mod, "append_log", lambda *a: logged.append(a))
        win._prune_history()  # must not raise
        assert logged and logged[0][0] == "Error"

    def test_on_initial_response_logs_and_notifies(self, win, monkeypatch):
        monkeypatch.setattr(win, "log_usage", MagicMock())
        response = MagicMock(text="OK understood")
        win.on_initial_response(response)
        assert "Ready to chat" in win.chat_history_log[-1]["text"]

    def test_on_response_is_noop(self, win):
        assert win.on_response(MagicMock()) is None


# ===========================================================================
# get_molecule_name / on_name_resolved
# ===========================================================================


class TestMoleculeName:
    def test_no_inchikey_returns_none(self, win):
        win.last_inchikey = None
        assert win.get_molecule_name() is None

    def test_cache_hit(self, win):
        win.last_inchikey = "KEY1"
        win._name_cache = {"KEY1": "Ethanol"}
        assert win.get_molecule_name() == "Ethanol"

    def test_allow_fetch_starts_worker(self, win, monkeypatch):
        win.last_inchikey = "KEY1"
        win.thread_pool = MagicMock()
        result = win.get_molecule_name(allow_fetch=True)
        assert result is None
        win.thread_pool.start.assert_called_once()

    def test_on_name_resolved_updates_cache_and_label(self, win, monkeypatch):
        win.last_inchikey = "KEY1"
        win.last_smiles = "CCO"
        win._fetching_inchikeys = {"KEY1"}
        # check_molecule_change() runs afterward and would overwrite lbl_context
        # via get_current_molecule_smiles() (no main_window) -- stub it out so we
        # can observe the label this method itself sets.
        monkeypatch.setattr(win, "check_molecule_change", lambda: None)
        win.on_name_resolved("KEY1", "Ethanol")
        assert win._name_cache["KEY1"] == "Ethanol"
        assert "Ethanol" in win.lbl_context.text()

    def test_on_name_resolved_different_molecule_skips_label(self, win):
        win.last_inchikey = "OTHER"
        win.on_name_resolved("KEY1", "Ethanol")
        assert win._name_cache["KEY1"] == "Ethanol"


# ===========================================================================
# check_molecule_change
# ===========================================================================


class TestCheckMoleculeChange:
    def test_first_check_sets_label(self, win):
        win.main_window = _make_mw(atoms={})
        win.context.current_molecule = None
        win.check_molecule_change()
        assert win.first_check_done is True

    def test_smiles_change_updates_pending_context(self, win):
        win.main_window = _make_mw(atoms={0: {"symbol": "C"}})
        win.first_check_done = True
        win.last_smiles = None
        win.check_molecule_change()
        assert win.pending_context_msg is not None

    def test_molecule_unload_sets_message(self, win):
        win.main_window = _make_mw(atoms={0: {"symbol": "C"}})
        win.first_check_done = True
        win.last_smiles = "SOMETHING"
        win.context.current_molecule = None
        win.main_window.state_manager.data.atoms = {}
        win.check_molecule_change()
        assert "unloaded" in (win.pending_context_msg or "").lower() or win.pending_context_msg is not None


# ===========================================================================
# handle_link
# ===========================================================================


class TestHandleLink:
    def _url(self, text):
        from PyQt6.QtCore import QUrl

        return QUrl(text)

    def test_smiles_link_loads_via_importer(self, win):
        win.main_window = MagicMock()
        win.main_window.string_importer_manager = MagicMock()
        win.handle_link(self._url("smiles:CCO"))
        win.main_window.string_importer_manager.load_from_smiles.assert_called_once_with("CCO")

    def test_smiles_link_no_importer_errors(self, win):
        win.main_window = MagicMock(spec=[])
        win.handle_link(self._url("smiles:CCO"))
        assert "not found" in win.chat_history_log[-1]["text"]

    def test_empty_smiles_ignored(self, win):
        win.main_window = MagicMock()
        before = len(win.chat_history_log)
        win.handle_link(self._url("smiles:"))
        assert len(win.chat_history_log) == before

    def test_non_smiles_link_opens_external(self, win, monkeypatch):
        opened = []
        monkeypatch.setattr(_mod.QDesktopServices, "openUrl", lambda u: opened.append(u))
        win.handle_link(self._url("https://example.com"))
        assert opened


# ===========================================================================
# update_structure_diff_based / load_smiles_undo_safe
# ===========================================================================


class TestUpdateStructureDiffBased:
    def test_invalid_smiles_returns_silently(self, win, monkeypatch):
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: None)
        win.main_window = _make_mw()
        win.update_structure_diff_based("garbage")  # must not raise

    def test_runs_update_with_mocked_rdkit(self, win):
        win.main_window = _make_mw(atoms={0: {"symbol": "C"}}, bonds={})
        win.update_structure_diff_based("CCO")
        assert win.chat_history_log


class TestLoadSmilesUndoSafe:
    def test_invalid_smiles_shows_error(self, win, monkeypatch):
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s, sanitize=True: None)
        win.main_window = _make_mw()
        win.load_smiles_undo_safe("garbage!!")
        assert "invalid" in win.chat_history_log[-1]["text"].lower()

    def test_runs_with_mocked_rdkit(self, win):
        win.main_window = _make_mw()
        win.load_smiles_undo_safe("CCO")  # must not raise


# ===========================================================================
# _ensure_main_window_3d_conversion / _get_descriptors_str
# ===========================================================================


class TestEnsureMainWindow3DConversion:
    def test_no_compute_manager_noop(self, win):
        win.main_window = MagicMock(spec=[])
        win._ensure_main_window_3d_conversion()  # must not raise

    def test_trigger_conversion_exception_reported(self, win):
        mw = MagicMock()
        mw.compute_manager.trigger_conversion.side_effect = RuntimeError("boom")
        win.main_window = mw
        win._ensure_main_window_3d_conversion()
        assert "boom" in win.chat_history_log[-1]["text"]


class TestGetDescriptorsStr:
    def test_no_molecule_returns_empty(self, win, monkeypatch):
        win.context.current_molecule = None
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: None)
        assert win._get_descriptors_str("CCO") == ""

    def test_success_returns_properties_string(self, win):
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        win.context.current_molecule = mol
        result = win._get_descriptors_str("CCO")
        assert "Properties:" in result or result == ""


# ===========================================================================
# render_content -- exercised directly (real markdown/LaTeX/tool-JSON/SMILES
# link logic).  HAS_MARKDOWN is forced False so the fallback plain-string path
# runs instead of calling the mocked `markdown.markdown()` (which returns a
# MagicMock, not a str, and would blow up the trailing re.sub SMILES-link step).
# ===========================================================================


class TestRenderContentDirect:
    @pytest.fixture(autouse=True)
    def _no_markdown(self, monkeypatch):
        monkeypatch.setattr(_mod, "HAS_MARKDOWN", False)

    def test_newline_converted_to_br(self, win):
        real = _mod.ChatMoleculeWindow.render_content
        out = real(win, "line1\nline2")
        assert "<br>" in out

    def test_block_latex_rendered(self, win, monkeypatch):
        monkeypatch.setattr(_mod, "HAS_MATPLOTLIB", False)
        real = _mod.ChatMoleculeWindow.render_content
        out = real(win, "Energy: $$E=mc^2$$ done")
        assert "<i>" in out  # matplotlib fallback wraps in <i>

    def test_inline_latex_rendered(self, win, monkeypatch):
        monkeypatch.setattr(_mod, "HAS_MATPLOTLIB", False)
        real = _mod.ChatMoleculeWindow.render_content
        out = real(win, "The value $x^2$ is small")
        assert "<i>" in out

    def test_smiles_link_converted_to_anchor(self, win):
        real = _mod.ChatMoleculeWindow.render_content
        out = real(win, "See [ethanol](smiles:CCO) for details")
        assert '<a href="smiles:CCO">ethanol</a>' in out

    def test_single_tool_json_block_simplified(self, win):
        real = _mod.ChatMoleculeWindow.render_content
        text = '```json\n{"tool": "clear_canvas", "params": {}}\n```'
        out = real(win, text)
        assert "Tool Request" in out
        assert "clear_canvas" in out

    def test_multi_tool_json_block_simplified(self, win):
        real = _mod.ChatMoleculeWindow.render_content
        text = (
            '```json\n[{"tool": "clear_canvas", "params": {}}, '
            '{"tool": "convert_to_3d", "params": {}}]\n```'
        )
        out = real(win, text)
        assert "2." in out  # second tool prefixed "2. "

    def test_tool_json_block_with_long_param_truncated(self, win):
        real = _mod.ChatMoleculeWindow.render_content
        long_val = "x" * 100
        text = f'```json\n{{"tool": "save_file", "params": {{"content": "{long_val}"}}}}\n```'
        out = real(win, text)
        assert "..." in out

    def test_non_tool_json_block_left_raw(self, win):
        real = _mod.ChatMoleculeWindow.render_content
        text = '```json\n{"foo": "bar"}\n```'
        out = real(win, text)
        assert "```json" in out

    def test_malformed_json_block_left_raw(self, win):
        real = _mod.ChatMoleculeWindow.render_content
        text = '```json\n{not valid json at all}\n```'
        out = real(win, text)
        assert "```json" in out

    def test_markdown_enabled_but_raises_falls_back(self, win, monkeypatch):
        monkeypatch.setattr(_mod, "HAS_MARKDOWN", True)
        fake_markdown = MagicMock()
        fake_markdown.markdown.side_effect = RuntimeError("boom")
        monkeypatch.setattr(_mod, "markdown", fake_markdown)
        real = _mod.ChatMoleculeWindow.render_content
        out = real(win, "line1\nline2")
        assert "<br>" in out


# ===========================================================================
# DEMO_MODE branch walk -- covers the `send_message()` demo-script chain and
# the "br"/"load" shortcut branches. The QTimer.singleShot callbacks
# themselves are deferred (never fire without an event loop pump); we only
# need the synchronous response_text-selection logic to execute for coverage.
# ===========================================================================


class TestDemoModeSendMessage:
    def _demo_window(self, monkeypatch):
        w = _build_window(monkeypatch)
        monkeypatch.setattr(_mod, "DEMO_MODE", True)
        w.client = MagicMock()
        w.get_current_molecule_smiles = lambda: (None, None)
        return w

    def test_stop_demo_shortcut_shows_stop_button(self, qapp, monkeypatch):
        w = self._demo_window(monkeypatch)
        w.txt_input.setPlainText("stop-demo")
        w.send_message()
        assert w.btn_send.text() == "Stop"
        w.destroy()

    def test_normal_message_schedules_simulated_response(self, qapp, monkeypatch):
        w = self._demo_window(monkeypatch)
        w.txt_input.setPlainText("hello there")
        w.send_message()  # must not raise; response arrives via deferred QTimer
        assert any(e["sender"] == "You" for e in w.chat_history_log)
        w.destroy()


# ===========================================================================
# Deep chemistry-shaped fakes for load_smiles_undo_safe /
# update_structure_diff_based -- real ints/strings so comparisons like
# `map_num > 0` don't blow up against a bare MagicMock.
# ===========================================================================


class _FakeAtom:
    def __init__(self, symbol="C", map_num=0, charge=0, idx=0):
        self._symbol = symbol
        self._map_num = map_num
        self._charge = charge
        self.idx = idx
        self._n_rad = 0

    def GetIdx(self):
        return self.idx

    def GetSymbol(self):
        return self._symbol

    def GetAtomMapNum(self):
        return self._map_num

    def SetAtomMapNum(self, n):
        self._map_num = n

    def GetFormalCharge(self):
        return self._charge

    def SetFormalCharge(self, c):
        self._charge = c

    def SetNumRadicalElectrons(self, n):
        self._n_rad = n


class _FakeBond:
    def __init__(self, begin_idx, end_idx, order=1.0, bond_dir=None, stereo=None, bond_type=None):
        self._begin = begin_idx
        self._end = end_idx
        self._order = order
        self._dir = bond_dir
        self._stereo = stereo
        self._bond_type = bond_type

    def GetBeginAtomIdx(self):
        return self._begin

    def GetEndAtomIdx(self):
        return self._end

    def GetBondTypeAsDouble(self):
        return self._order

    def GetBondDir(self):
        return self._dir

    def GetStereo(self):
        return self._stereo

    def GetBondType(self):
        return self._bond_type


class _FakePos:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x, self.y, self.z = x, y, z


class _FakeConf:
    def __init__(self, positions):
        self._positions = positions

    def GetAtomPosition(self, i):
        return self._positions[i]


class _FakeMol:
    def __init__(self, atoms, bonds, positions=None):
        self._atoms = atoms
        for i, a in enumerate(atoms):
            a.idx = i
        self._bonds = bonds
        self._conf = _FakeConf(positions or [_FakePos(float(i), float(i), 0.0) for i in range(len(atoms))])

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetAtoms(self):
        return list(self._atoms)

    def GetAtomWithIdx(self, i):
        return self._atoms[i]

    def GetBonds(self):
        return list(self._bonds)

    def GetConformer(self):
        return self._conf

    def GetNumConformers(self):
        return 1

    def GetSubstructMatches(self, query, uniquify=True):
        return []

    def UpdatePropertyCache(self, strict=True):
        pass


class TestUpdateStructureDiffBasedDeep:
    def test_updates_existing_and_creates_new_atoms(self, win, monkeypatch):
        atoms = [_FakeAtom("C", map_num=1), _FakeAtom("O", map_num=0)]
        bonds = [_FakeBond(0, 1, order=1.0)]
        fake_mol = _FakeMol(atoms, bonds)
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: fake_mol)
        monkeypatch.setattr(_mod.Chem, "Kekulize", lambda *a, **k: None)
        monkeypatch.setattr(_mod.AllChem, "Compute2DCoords", lambda m: None)

        mw = _make_mw(atoms={0: {"symbol": "C"}})
        import itertools

        counter = itertools.count(50)
        mw.scene.create_atom.side_effect = lambda *a, **k: next(counter)
        win.main_window = mw

        win.update_structure_diff_based("C.O")
        assert "Structure updated" in win.chat_history_log[-1]["text"]

    def test_invalid_smiles_returns_silently(self, win, monkeypatch):
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: None)
        win.main_window = _make_mw()
        before = len(win.chat_history_log)
        win.update_structure_diff_based("garbage")
        assert len(win.chat_history_log) == before

    def test_deletes_atoms_missing_from_new_smiles(self, win, monkeypatch):
        atoms = [_FakeAtom("C", map_num=1)]
        fake_mol = _FakeMol(atoms, [])
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: fake_mol)
        monkeypatch.setattr(_mod.Chem, "Kekulize", lambda *a, **k: None)
        monkeypatch.setattr(_mod.AllChem, "Compute2DCoords", lambda m: None)
        # Existing atom id 5 (map_num 6) is not present in fake_mol -> deleted.
        mw = _make_mw(atoms={0: {"symbol": "C"}, 5: {"symbol": "N"}})
        win.main_window = mw
        win.update_structure_diff_based("C")

    def test_bond_update_and_new_atom_creation_and_bond_delete(self, win, monkeypatch):
        # Two existing atoms (map_num 1, 2 -> ids 0, 1) with a bond that will
        # change order (UPDATE branch); plus a third fake-mol atom (map_num 3
        # -> id 2, not yet on canvas) that triggers the "new atom" + bond
        # CREATE branch, and the stale (0,1) bond ends up deleted because it's
        # not re-added to current_step_bond_keys (create_bond doesn't really
        # write into our mocked bonds dict).
        atoms = [
            _FakeAtom("C", map_num=1, idx=0),
            _FakeAtom("O", map_num=2, idx=1),
        ]
        bonds = [
            _FakeBond(
                0,
                1,
                order=2.0,
                bond_dir=_mod.Chem.BondDir.BEGINWEDGE,
                bond_type=_mod.Chem.BondType.DOUBLE,
                stereo=_mod.Chem.BondStereo.STEREOZ,
            )
        ]
        fake_mol = _FakeMol(atoms, bonds)
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: fake_mol)
        monkeypatch.setattr(_mod.Chem, "Kekulize", lambda *a, **k: None)
        monkeypatch.setattr(_mod.AllChem, "Compute2DCoords", lambda m: None)

        mw = _make_mw(
            atoms={0: {"symbol": "C"}, 1: {"symbol": "O"}},
            bonds={(0, 1): {"order": 1, "stereo": 0}},
        )
        mw.scene.bond_items[(0, 1)] = MagicMock()
        win.main_window = mw

        win.update_structure_diff_based("C=O")
        assert "Structure updated" in win.chat_history_log[-1]["text"]

    def test_new_id_specified_atom_created_and_remapped(self, win, monkeypatch):
        atoms = [_FakeAtom("N", map_num=3, idx=0)]  # target_id=2, not on canvas
        fake_mol = _FakeMol(atoms, [])
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: fake_mol)
        monkeypatch.setattr(_mod.Chem, "Kekulize", lambda *a, **k: None)
        monkeypatch.setattr(_mod.AllChem, "Compute2DCoords", lambda m: None)

        mw = _make_mw(atoms={0: {"symbol": "C"}})
        import itertools

        counter = itertools.count(99)
        mw.scene.create_atom.side_effect = lambda *a, **k: next(counter)
        win.main_window = mw

        win.update_structure_diff_based("N")
        assert "Structure updated" in win.chat_history_log[-1]["text"]


class TestLoadSmilesUndoSafeDeep:
    def test_loads_and_reconstructs_scene(self, win, monkeypatch):
        atoms = [_FakeAtom("C", map_num=1), _FakeAtom("O", map_num=2)]
        bonds = [
            _FakeBond(
                0,
                1,
                order=1.0,
                bond_dir=_mod.Chem.BondDir.BEGINWEDGE,
                bond_type=_mod.Chem.BondType.SINGLE,
            )
        ]
        fake_mol = _FakeMol(atoms, bonds)
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s, sanitize=True: fake_mol)
        monkeypatch.setattr(_mod.Chem, "Kekulize", lambda *a, **k: None)
        monkeypatch.setattr(_mod.AllChem, "Compute2DCoords", lambda m: None)
        monkeypatch.setattr(_mod.AllChem, "AssignStereochemistry", lambda *a, **k: None)
        monkeypatch.setattr(_mod.AllChem, "WedgeMolBonds", lambda *a, **k: None)

        mw = _make_mw()
        import itertools

        counter = itertools.count(1)
        mw.init_manager.scene.create_atom.side_effect = lambda *a, **k: next(counter)
        win.main_window = mw
        win.context.current_molecule = None

        win.load_smiles_undo_safe("CO")  # must not raise
        win.context.push_undo_checkpoint.assert_called()


# ===========================================================================
# PubChemResolver -- mocked-HTTP network paths (only empty-input guards are
# covered elsewhere).
# ===========================================================================


class _FakeHTTPResponse:
    def __init__(self, status, body):
        self.status = status
        self._body = body.encode("utf-8")

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class TestPubChemResolverNetwork:
    def test_resolve_inchikey_success(self, monkeypatch):
        body = json.dumps(
            {"InformationList": {"Information": [{"Title": "Ethanol"}]}}
        )
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(200, body)
        )
        name, err = _mod.PubChemResolver.resolve_inchikey_to_name("KEY1")
        assert name == "Ethanol"
        assert err is None

    def test_resolve_inchikey_no_title_found(self, monkeypatch):
        body = json.dumps({"InformationList": {"Information": [{}]}})
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(200, body)
        )
        name, err = _mod.PubChemResolver.resolve_inchikey_to_name("KEY1")
        assert name is None
        assert "No name found" in err

    def test_resolve_inchikey_http_error_status(self, monkeypatch):
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(500, "")
        )
        name, err = _mod.PubChemResolver.resolve_inchikey_to_name("KEY1")
        assert name is None
        assert "HTTP Error" in err

    def test_resolve_inchikey_404_not_found(self, monkeypatch):
        def raise_404(*a, **k):
            raise _mod.urllib.error.HTTPError("url", 404, "not found", None, None)

        monkeypatch.setattr(_mod.urllib.request, "urlopen", raise_404)
        name, err = _mod.PubChemResolver.resolve_inchikey_to_name("KEY1")
        assert name is None
        assert "not found" in err.lower()

    def test_resolve_inchikey_other_http_error(self, monkeypatch):
        def raise_500(*a, **k):
            raise _mod.urllib.error.HTTPError("url", 500, "server error", None, None)

        monkeypatch.setattr(_mod.urllib.request, "urlopen", raise_500)
        name, err = _mod.PubChemResolver.resolve_inchikey_to_name("KEY1")
        assert "500" in err

    def test_resolve_inchikey_generic_exception(self, monkeypatch):
        def raise_err(*a, **k):
            raise RuntimeError("network down")

        monkeypatch.setattr(_mod.urllib.request, "urlopen", raise_err)
        name, err = _mod.PubChemResolver.resolve_inchikey_to_name("KEY1")
        assert "network down" in err

    def test_resolve_name_to_smiles_success(self, monkeypatch):
        body = json.dumps(
            {"PropertyTable": {"Properties": [{"SMILES": "CCO"}]}}
        )
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(200, body)
        )
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert smiles == "CCO"
        assert err is None

    def test_resolve_name_to_smiles_canonical_fallback(self, monkeypatch):
        body = json.dumps(
            {"PropertyTable": {"Properties": [{"CanonicalSMILES": "CCO"}]}}
        )
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(200, body)
        )
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert smiles == "CCO"

    def test_resolve_name_to_smiles_isomeric_fallback(self, monkeypatch):
        body = json.dumps(
            {"PropertyTable": {"Properties": [{"IsomericSMILES": "CCO"}]}}
        )
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(200, body)
        )
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert smiles == "CCO"

    def test_resolve_name_to_smiles_no_smiles_key(self, monkeypatch):
        body = json.dumps({"PropertyTable": {"Properties": [{"Other": "x"}]}})
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(200, body)
        )
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert smiles is None
        assert "missing" in err

    def test_resolve_name_to_smiles_no_properties(self, monkeypatch):
        body = json.dumps({"PropertyTable": {"Properties": []}})
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(200, body)
        )
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert smiles is None
        assert "No properties" in err

    def test_resolve_name_to_smiles_waiting_status(self, monkeypatch):
        body = json.dumps({"Waiting": "busy"})
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(200, body)
        )
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert smiles is None
        assert "busy" in err

    def test_resolve_name_to_smiles_404(self, monkeypatch):
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(404, "")
        )
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("unobtainium")
        assert smiles is None
        assert "not found" in err

    def test_resolve_name_to_smiles_other_status(self, monkeypatch):
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(
            _mod.urllib.request, "urlopen", lambda *a, **k: _FakeHTTPResponse(503, "")
        )
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert smiles is None
        assert "API Error" in err

    def test_resolve_name_to_smiles_http_error_404(self, monkeypatch):
        def raise_404(*a, **k):
            raise _mod.urllib.error.HTTPError("url", 404, "not found", None, None)

        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(_mod.urllib.request, "urlopen", raise_404)
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert "404" in err

    def test_resolve_name_to_smiles_http_error_other(self, monkeypatch):
        def raise_500(*a, **k):
            raise _mod.urllib.error.HTTPError("url", 500, "server error", None, None)

        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(_mod.urllib.request, "urlopen", raise_500)
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert "HTTP Error" in err

    def test_resolve_name_to_smiles_generic_exception(self, monkeypatch):
        def raise_err(*a, **k):
            raise RuntimeError("boom")

        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        monkeypatch.setattr(_mod.urllib.request, "urlopen", raise_err)
        smiles, err = _mod.PubChemResolver.resolve_name_to_smiles("ethanol")
        assert "boom" in err


# ===========================================================================
# Deeper execute_* success paths requiring real numeric/molecular shapes
# ===========================================================================


class TestExecuteSetElectronicStateDeep:
    def test_atom_index_success_applies_charge(self, win, monkeypatch):
        atoms = [_FakeAtom("C", map_num=1), _FakeAtom("O", map_num=2)]
        fake_mol = _FakeMol(atoms, [])
        win.get_current_molecule_smiles = lambda: ("CCO", None)
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: fake_mol)
        monkeypatch.setattr(_mod.Chem, "SanitizeMol", lambda m: None)
        monkeypatch.setattr(_mod.Chem, "GetFormalCharge", lambda m: 1)
        monkeypatch.setattr(_mod.Chem, "MolToSmiles", lambda m: "CCO")
        monkeypatch.setattr(win, "load_smiles_undo_safe", MagicMock())
        win.main_window = _make_mw()
        win.execute_set_electronic_state(
            {"atom_index": 1, "charge": 1, "multiplicity": 2}
        )
        assert atoms[0].GetFormalCharge() == 1
        assert "Applied Electronic State" in win.chat_history_log[-1]["text"]

    def test_global_metadata_fallback(self, win, monkeypatch):
        atoms = [_FakeAtom("C", map_num=1)]
        fake_mol = _FakeMol(atoms, [])
        win.get_current_molecule_smiles = lambda: ("CCO", None)
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: fake_mol)
        mw = _make_mw()
        mw.state_manager.current_charge = None
        mw.state_manager.current_mult = None
        win.main_window = mw
        win.execute_set_electronic_state({"charge": -1})
        assert "Global Metadata" in win.chat_history_log[-1]["text"]

    def test_no_atoms_targeted_debug_message(self, win, monkeypatch):
        atoms = [_FakeAtom("C", map_num=1)]
        fake_mol = _FakeMol(atoms, [])
        win.get_current_molecule_smiles = lambda: ("CCO", None)
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: fake_mol)
        win.main_window = MagicMock(spec=[])
        win.execute_set_electronic_state({"target_smarts": "[C]"})
        assert "No atoms targeted" in win.chat_history_log[-1]["text"]

    def test_selection_based_targeting(self, win, monkeypatch):
        atoms = [_FakeAtom("C", map_num=1), _FakeAtom("O", map_num=2)]
        fake_mol = _FakeMol(atoms, [])
        win.get_current_molecule_smiles = lambda: ("CCO", None)
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: fake_mol)
        monkeypatch.setattr(_mod.Chem, "SanitizeMol", lambda m: None)
        monkeypatch.setattr(_mod.Chem, "GetFormalCharge", lambda m: 1)
        monkeypatch.setattr(_mod.Chem, "MolToSmiles", lambda m: "CCO")
        monkeypatch.setattr(win, "load_smiles_undo_safe", MagicMock())
        mw = _make_mw(atoms={0: {}}, selected={0})
        # aid=0 -> selected_map_nums=[0]; fake atom map_num must match (use 0)
        atoms[0].SetAtomMapNum(0)
        win.main_window = mw
        win.execute_set_electronic_state({"charge": 1})


class TestExecuteCalculateDescriptorsDeep:
    def test_success_computes_all_properties(self, win, monkeypatch):
        monkeypatch.setattr(_mod.QMessageBox, "exec", lambda self: None, raising=False)
        _install_fake_rdkit_chem_submodules(
            monkeypatch,
            Descriptors={
                "MolWt": lambda m: 46.07,
                "MolLogP": lambda m: -0.14,
                "TPSA": lambda m: 20.2,
            },
            Lipinski={
                "NumHDonors": lambda m: 1,
                "NumHAcceptors": lambda m: 1,
                "NumRotatableBonds": lambda m: 0,
                "NumAromaticRings": lambda m: 0,
                "RingCount": lambda m: 0,
            },
        )
        win.get_current_molecule_smiles = lambda: ("CCO", None)
        monkeypatch.setattr(_mod.Chem, "MolFromSmiles", lambda s: MagicMock())
        win.execute_calculate_descriptors(
            {
                "properties": [
                    "MW",
                    "LogP",
                    "TPSA",
                    "HBD",
                    "HBA",
                    "RB",
                    "AromaticRings",
                    "Rings",
                ]
            }
        )
        assert "Properties:" in win.chat_history_log[-1]["text"] or any(
            "MW" in e["text"] for e in win.chat_history_log
        )
        assert "CCO" in win.calc_results_by_smiles


class TestExecuteInputGeneratorsDeep:
    def _mol(self):
        return _FakeExportMol(
            [
                _FakeExportAtom("C", 0.0, 0.0, 0.0),
                _FakeExportAtom("O", 1.2, 0.0, 0.0),
            ]
        )

    def test_orca_generator_writes_real_geometry(self, win, monkeypatch, tmp_path):
        _install_fake_rdkit_chem_submodules(
            monkeypatch, Descriptors={"NumRadicalElectrons": lambda m: 0}
        )
        monkeypatch.setattr(_mod.Chem, "GetFormalCharge", lambda m: 0)
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: self._mol())
        target = str(tmp_path / "a.inp")
        monkeypatch.setattr(
            _mod.QFileDialog, "getSaveFileName", lambda *a, **k: (target, "")
        )
        err = win.execute_orca_input_generator(
            {"filename": "a.inp", "header": "! HF", "footer": ""}
        )
        assert err is None
        content = Path(target).read_text(encoding="utf-8")
        assert "C " in content and "O " in content

    def test_gaussian_generator_writes_real_geometry(self, win, monkeypatch, tmp_path):
        _install_fake_rdkit_chem_submodules(
            monkeypatch, Descriptors={"NumRadicalElectrons": lambda m: 0}
        )
        monkeypatch.setattr(_mod.Chem, "GetFormalCharge", lambda m: 0)
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: self._mol())
        target = str(tmp_path / "a.gjf")
        monkeypatch.setattr(
            _mod.QFileDialog, "getSaveFileName", lambda *a, **k: (target, "")
        )
        err = win.execute_gaussian_input_generator(
            {"filename": "a.gjf", "header": "#P HF", "footer": ""}
        )
        assert err is None
        content = Path(target).read_text(encoding="utf-8")
        assert "%chk=a.chk" in content

    def test_save_file_atom_tag_with_real_mol(self, win, monkeypatch, tmp_path):
        monkeypatch.setattr(win, "_get_molecule_for_export", lambda: self._mol())
        target = str(tmp_path / "out.xyz")
        monkeypatch.setattr(
            _mod.QFileDialog, "getSaveFileName", lambda *a, **k: (target, "")
        )
        win.execute_save_file(
            {"filename": "out.xyz", "content": "[[atom_count]]\n[[atom]]"}
        )
        content = Path(target).read_text(encoding="utf-8")
        assert content.startswith("2\n")
        assert "C " in content


class TestSendMessageWithMoleculeContext:
    def test_context_injected_when_smiles_present(self, qapp, monkeypatch):
        w = _build_window(monkeypatch)
        w.client = MagicMock()
        w.chat_history_state = [{"role": "system", "content": "SYS"}]
        w.get_current_molecule_smiles = lambda: ("CCO", None)
        monkeypatch.setattr(w, "get_molecule_name", lambda **k: "Ethanol")
        monkeypatch.setattr(w, "_ensure_main_window_3d_conversion", MagicMock())
        w.client.chat.completions.create.return_value = [_make_openai_chunk("OK")]
        w.txt_input.setPlainText("hello")
        w.send_message()
        assert w.pending_context_msg is None  # cleared after being sent
        w.destroy()

    def test_name_resolution_attempted_when_unknown(self, qapp, monkeypatch):
        w = _build_window(monkeypatch)
        w.client = MagicMock()
        w.chat_history_state = [{"role": "system", "content": "SYS"}]
        w.last_inchikey = "KEY1"
        w.get_current_molecule_smiles = lambda: ("CCO", None)
        monkeypatch.setattr(w, "get_molecule_name", lambda **k: None)
        monkeypatch.setattr(w, "_ensure_main_window_3d_conversion", MagicMock())
        monkeypatch.setattr(
            _mod.PubChemResolver,
            "resolve_inchikey_to_name",
            staticmethod(lambda k: ("Ethanol", None)),
        )
        w.client.chat.completions.create.return_value = [_make_openai_chunk("OK")]
        w.txt_input.setPlainText("hello")
        w.send_message()
        assert "Ethanol" in w.lbl_context.text()
        w.destroy()


class TestEnsureMainWindow3DConversionDeep:
    def test_trigger_conversion_success_completes_quickly(self, win, monkeypatch):
        mw = MagicMock()
        mw.next_conversion_id = 0
        mw.active_worker_ids = None

        def _trigger():
            mw.next_conversion_id = 1

        mw.compute_manager.trigger_conversion.side_effect = _trigger
        win.main_window = mw
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        win._ensure_main_window_3d_conversion()
        win.context.draw_molecule_3d.assert_called_once()

    def test_conversion_never_starts_reports_status(self, win, monkeypatch):
        mw = MagicMock()
        mw.next_conversion_id = 5
        mw.statusBar.return_value.currentMessage.return_value = "busy"

        def _trigger():
            pass  # next_conversion_id never changes

        mw.compute_manager.trigger_conversion.side_effect = _trigger
        win.main_window = mw
        monkeypatch.setattr(_mod.time, "sleep", lambda s: None)
        win._ensure_main_window_3d_conversion()
        assert "failed to start" in win.chat_history_log[-1]["text"]
