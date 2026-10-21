"""Exercise real native threads, rather than calling worker.run() on the GUI."""

import gc
import threading
import time
import weakref
from pathlib import Path

import pytest
from PySide6.QtCore import QCoreApplication, QEvent, QThreadPool, Qt

from conftest import load_plugin_for_gui, mock_chemistry_imports


ROOT = Path(__file__).resolve().parents[1] / "plugins"
PATHS = [
    ROOT / "Chat_with_Molecule_Neo" / "chat_with_molecule_neo.py",
    ROOT / "Chat_with_Molecule_Neo_ChatGPT" / "chat_with_molecule_neo_chatGPT.py",
    ROOT / "Chat_with_Molecule_Neo_Local" / "chat_with_molecule_neo_local.py",
]


@pytest.fixture(params=PATHS, ids=lambda path: path.parent.name)
def chat_module(request, qapp):
    with mock_chemistry_imports():
        return load_plugin_for_gui(request.param)


def pump_until(qapp, predicate):
    deadline = time.monotonic() + 5
    while not predicate() and time.monotonic() < deadline:
        qapp.processEvents()
        time.sleep(0.005)
    assert predicate()


def test_model_result_does_not_release_running_native_thread(chat_module, qapp, monkeypatch):
    mod = chat_module
    release = threading.Event()
    results = []

    def run(worker):
        worker.completed.emit(["synthetic-model"], "")
        release.wait(5)

    monkeypatch.setattr(mod.InitWorker, "run", run)
    args = ("synthetic-key", "http://localhost") if "Local" in mod.PLUGIN_NAME else ("synthetic-key",)
    worker = mod.InitWorker(*args)
    worker.completed.connect(lambda models, error: results.append((models, error)))
    ref = weakref.ref(worker)
    worker.start()
    del worker
    gc.collect()
    try:
        pump_until(qapp, lambda: bool(results))
        assert results == [(["synthetic-model"], "")]
        assert ref() is not None
        assert ref().isRunning()
    finally:
        release.set()
        if ref() is not None:
            assert ref().wait(5000)
    pump_until(qapp, lambda: not mod._active_threads)
    gc.collect()
    assert ref() is None


def test_resolver_signal_object_is_destroyed_on_gui_thread(chat_module, qapp, monkeypatch):
    mod = chat_module
    monkeypatch.setattr(mod, "PubChemResolver", None)
    destroyed_on = []
    gui_thread = threading.get_ident()
    pool = QThreadPool()
    worker = mod.NameResolverWorker("synthetic-inchikey")
    worker.signals.destroyed.connect(
        lambda: destroyed_on.append(threading.get_ident()),
        Qt.ConnectionType.DirectConnection,
    )
    pool.start(worker)
    del worker
    assert pool.waitForDone(5000)
    gc.collect()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    pump_until(qapp, lambda: bool(destroyed_on))
    assert destroyed_on == [gui_thread]


def test_old_stream_completion_does_not_clear_replacement(chat_module, qapp, monkeypatch):
    from PySide6.QtWidgets import QDialog

    window = chat_module.ChatMoleculeWindow.__new__(chat_module.ChatMoleculeWindow)
    QDialog.__init__(window)
    previous, replacement = object(), object()
    window.worker = replacement
    monkeypatch.setattr(window, "sender", lambda: previous)
    window.on_worker_finished()
    assert window.worker is replacement
    monkeypatch.setattr(window, "sender", lambda: replacement)
    window.on_worker_finished()
    assert window.worker is None
