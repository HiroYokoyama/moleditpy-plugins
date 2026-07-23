"""
Headless GUI tests for the Gaussian FCHK Loader plugin.

Covers: FCHKLoaderDialog.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

FCHK_PATH = PLUGINS_DIR / "Gaussian_FCHK_Loader" / "gaussian_fchk_loader.py"

with mock_chemistry_imports():
    _fchk = load_plugin_for_gui(FCHK_PATH)


# ===========================================================================
# FCHKLoaderDialog  (Gaussian FCHK Loader)
# ===========================================================================


def _fchk_context() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = MagicMock()
    return ctx


class TestFCHKLoaderDialog:
    """FCHKLoaderDialog with a fake fchk path and MagicMock context."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _fchk_context()
        d = _fchk.FCHKLoaderDialog(
            parent=None, context=ctx, fchk_path="/tmp/test_molecule.fchk"
        )
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Select FCHK Analysis Mode"

    def test_freq_button_exists(self, dlg):
        assert "Frequency" in dlg.btn_freq.text()

    def test_mo_button_exists(self, dlg):
        assert "MO" in dlg.btn_mo.text() or "Orbital" in dlg.btn_mo.text()

    def test_freq_button_disabled_when_plugin_absent(self, dlg):
        # In test environment the freq analyzer plugin is not found via
        # find_file_recursive, so the button should be disabled.
        if dlg.freq_analyzer_path is None:
            assert not dlg.btn_freq.isEnabled()

    def test_fchk_path_stored(self, dlg):
        assert dlg.fchk_path == "/tmp/test_molecule.fchk"

    def test_basename_shown_in_label(self, dlg):
        # The label text is set during init_ui from os.path.basename(fchk_path)
        # Verify we can find text containing the filename in the dialog's children
        from PyQt6.QtWidgets import QLabel
        labels = dlg.findChildren(QLabel)
        label_texts = " ".join(lbl.text() for lbl in labels)
        assert "test_molecule.fchk" in label_texts


# ===========================================================================
# find_freq_analyzer / find_mo_analyzer branch coverage — isolated file copies
# ===========================================================================


class TestFindAnalyzerBranches:
    def test_no_analyzers_disables_buttons_and_sets_tooltips(
        self, qapp, tmp_path, monkeypatch
    ):
        # Fake __file__ under an isolated dir (not named "plugins" at any
        # level) -> else branch of find_freq_analyzer, nothing found.
        fake_file = tmp_path / "lonely" / "gaussian_fchk_loader.py"
        monkeypatch.setattr(_fchk, "__file__", str(fake_file))
        ctx = _fchk_context()
        dlg = _fchk.FCHKLoaderDialog(parent=None, context=ctx, fchk_path="/tmp/x.fchk")
        try:
            assert dlg.freq_analyzer_path is None
            assert dlg.mo_analyzer_pkg is None
            assert not dlg.btn_freq.isEnabled()
            assert not dlg.btn_mo.isEnabled()
            assert dlg.btn_freq.toolTip() == "Plugin not found"
            assert dlg.btn_mo.toolTip() == "Plugin not found"
        finally:
            dlg.destroy()

    def test_root_named_plugins_elif_branch(self, qapp, tmp_path, monkeypatch):
        # Fake __file__ living directly inside a dir literally named
        # "plugins" -> elif branch of find_freq_analyzer.
        fake_file = tmp_path / "plugins" / "gaussian_fchk_loader.py"
        monkeypatch.setattr(_fchk, "__file__", str(fake_file))
        ctx = _fchk_context()
        dlg = _fchk.FCHKLoaderDialog(parent=None, context=ctx, fchk_path="/tmp/x.fchk")
        try:
            assert dlg.freq_analyzer_path is None
        finally:
            dlg.destroy()


# ===========================================================================
# launch_freq / launch_mo / close_existing_analyzers — real Qt execution
# ===========================================================================


def _make_real_mw(qapp):
    from PyQt6.QtWidgets import QMainWindow

    return QMainWindow()


class TestLaunchFreqRealQt:
    def test_missing_path_is_noop(self, qapp):
        ctx = _fchk_context()
        ctx.get_main_window.return_value = _make_real_mw(qapp)
        dlg = _fchk.FCHKLoaderDialog(parent=None, context=ctx, fchk_path="/tmp/x.fchk")
        dlg.freq_analyzer_path = None
        dlg.launch_freq()
        assert not hasattr(dlg, "analyzer")
        dlg.destroy()

    def test_success_creates_dock_and_loads_file(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QWidget

        class _FakeAnalyzer(QWidget):
            def __init__(self, context, dock_widget=None):
                super().__init__()
                self.context = context
                self.dock = dock_widget
                self.loaded_path = None

            def load_file(self, path):
                self.loaded_path = path

        fake_mod = SimpleNamespace(GaussianFCHKFreqAnalyzer=_FakeAnalyzer)
        mw = _make_real_mw(qapp)
        ctx = _fchk_context()
        ctx.get_main_window.return_value = mw
        dlg = _fchk.FCHKLoaderDialog(parent=None, context=ctx, fchk_path="/tmp/mol.fchk")
        dlg.freq_analyzer_path = "/fake/gaussian_fchk_freq_analyzer.py"
        monkeypatch.setattr(_fchk, "load_module_from_path", lambda name, path: fake_mod)
        dlg.launch_freq()
        assert dlg.analyzer.loaded_path == "/tmp/mol.fchk"
        assert dlg.result() == 1  # QDialog.Accepted
        dlg.destroy()

    def test_exception_shows_error_and_reenables_buttons(self, qapp, monkeypatch):
        mw = _make_real_mw(qapp)
        ctx = _fchk_context()
        ctx.get_main_window.return_value = mw
        dlg = _fchk.FCHKLoaderDialog(parent=None, context=ctx, fchk_path="/tmp/mol.fchk")
        dlg.freq_analyzer_path = "/fake/gaussian_fchk_freq_analyzer.py"

        def _boom(name, path):
            raise RuntimeError("boom")

        monkeypatch.setattr(_fchk, "load_module_from_path", _boom)
        critical_calls = []
        monkeypatch.setattr(
            _fchk.QMessageBox,
            "critical",
            staticmethod(lambda *a, **k: critical_calls.append(a)),
        )
        dlg.launch_freq()
        assert len(critical_calls) == 1
        assert dlg.btn_freq.isEnabled()
        assert dlg.btn_mo.isEnabled()
        dlg.destroy()


class TestLaunchMoRealQt:
    def test_missing_pkg_is_noop(self, qapp):
        ctx = _fchk_context()
        ctx.get_main_window.return_value = _make_real_mw(qapp)
        dlg = _fchk.FCHKLoaderDialog(parent=None, context=ctx, fchk_path="/tmp/x.fchk")
        dlg.mo_analyzer_pkg = None
        dlg.launch_mo()
        assert not hasattr(dlg, "mo_widget")
        dlg.destroy()

    def test_success_creates_dock(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QWidget

        class _FakeOrbitalWidget(QWidget):
            def __init__(self, mw, context, fchk_path):
                super().__init__()
                self.mw = mw
                self.context = context
                self.fchk_path = fchk_path

        fake_mod = SimpleNamespace(OrbitalWidget=_FakeOrbitalWidget)

        class _FakeImportlib:
            def import_module(self, name):
                return fake_mod

            def reload(self, mod):
                return mod

        mw = _make_real_mw(qapp)
        ctx = _fchk_context()
        ctx.get_main_window.return_value = mw
        dlg = _fchk.FCHKLoaderDialog(parent=None, context=ctx, fchk_path="/tmp/mol.fchk")
        dlg.mo_analyzer_pkg = "/fake/gaussian_fchk_mo_analyzer"
        monkeypatch.setattr(_fchk, "importlib", _FakeImportlib())
        dlg.launch_mo()
        assert dlg.mo_widget.fchk_path == "/tmp/mol.fchk"
        assert dlg.result() == 1  # QDialog.Accepted
        dlg.destroy()

    def test_exception_shows_error_and_reenables_buttons(self, qapp, monkeypatch):
        class _FakeImportlib:
            def import_module(self, name):
                raise RuntimeError("boom")

        mw = _make_real_mw(qapp)
        ctx = _fchk_context()
        ctx.get_main_window.return_value = mw
        dlg = _fchk.FCHKLoaderDialog(parent=None, context=ctx, fchk_path="/tmp/mol.fchk")
        dlg.mo_analyzer_pkg = "/fake/gaussian_fchk_mo_analyzer"
        monkeypatch.setattr(_fchk, "importlib", _FakeImportlib())
        critical_calls = []
        monkeypatch.setattr(
            _fchk.QMessageBox,
            "critical",
            staticmethod(lambda *a, **k: critical_calls.append(a)),
        )
        dlg.launch_mo()
        assert len(critical_calls) == 1
        assert dlg.btn_freq.isEnabled()
        assert dlg.btn_mo.isEnabled()
        dlg.destroy()


class TestCloseExistingAnalyzersRealQt:
    def test_closes_matching_docks_only(self, qapp):
        from PyQt6.QtWidgets import QDockWidget

        mw = _make_real_mw(qapp)
        ctx = _fchk_context()
        ctx.get_main_window.return_value = mw
        dlg = _fchk.FCHKLoaderDialog(parent=None, context=ctx, fchk_path="/tmp/x.fchk")

        freq_dock = QDockWidget("Gaussian Freq Analyzer", mw)
        mo_dock = QDockWidget("Gaussian MO Analyzer", mw)
        other_dock = QDockWidget("Something Else", mw)
        mw.addDockWidget(_fchk.Qt.DockWidgetArea.RightDockWidgetArea, freq_dock)
        mw.addDockWidget(_fchk.Qt.DockWidgetArea.RightDockWidgetArea, mo_dock)
        mw.addDockWidget(_fchk.Qt.DockWidgetArea.RightDockWidgetArea, other_dock)

        dlg.close_existing_analyzers()

        assert other_dock.windowTitle() == "Something Else"
        dlg.destroy()


# ===========================================================================
# run(mw) legacy entry point — real Qt run_dialog execution
# ===========================================================================


class _FakeFCHKDialogGui:
    can_freq = True
    can_mo = False

    def __init__(self, mw, context, path):
        self.mw, self.context, self.path = mw, context, path
        self.btn_freq = SimpleNamespace(isEnabled=lambda: _FakeFCHKDialogGui.can_freq)
        self.btn_mo = SimpleNamespace(isEnabled=lambda: _FakeFCHKDialogGui.can_mo)
        self.launch_freq_calls = 0
        self.launch_mo_calls = 0
        self.exec_calls = 0

    def launch_freq(self):
        self.launch_freq_calls += 1

    def launch_mo(self):
        self.launch_mo_calls += 1

    def exec(self):
        self.exec_calls += 1


class TestRunFunctionRealQt:
    def _run_with(self, qapp, monkeypatch, tmp_path, can_freq, can_mo):
        from PyQt6.QtWidgets import QFileDialog, QMessageBox

        fchk_file = tmp_path / "mol.fchk"
        fchk_file.write_text("")

        captured = {}

        class _CapDlg(_FakeFCHKDialogGui):
            def __init__(self, mw, context, path):
                super().__init__(mw, context, path)
                captured["dlg"] = self

        _FakeFCHKDialogGui.can_freq = can_freq
        _FakeFCHKDialogGui.can_mo = can_mo
        monkeypatch.setattr(_fchk, "FCHKLoaderDialog", _CapDlg)
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: (str(fchk_file), "")),
        )
        monkeypatch.setattr(
            _fchk.QTimer, "singleShot", staticmethod(lambda ms, fn: fn())
        )
        warn_calls = []
        monkeypatch.setattr(
            QMessageBox, "warning", staticmethod(lambda *a, **k: warn_calls.append(a))
        )
        ctx = _fchk_context()
        _fchk.PLUGIN_CONTEXT = ctx
        mw = _make_real_mw(qapp)
        _fchk.run(mw)
        return captured["dlg"], warn_calls

    def test_both_available_shows_dialog(self, qapp, monkeypatch, tmp_path):
        dlg, warn_calls = self._run_with(qapp, monkeypatch, tmp_path, True, True)
        assert dlg.exec_calls == 1
        assert dlg.launch_freq_calls == 0
        assert dlg.launch_mo_calls == 0

    def test_only_freq_available(self, qapp, monkeypatch, tmp_path):
        dlg, warn_calls = self._run_with(qapp, monkeypatch, tmp_path, True, False)
        assert dlg.launch_freq_calls == 1
        assert dlg.launch_mo_calls == 0
        assert dlg.exec_calls == 0

    def test_only_mo_available(self, qapp, monkeypatch, tmp_path):
        dlg, warn_calls = self._run_with(qapp, monkeypatch, tmp_path, False, True)
        assert dlg.launch_mo_calls == 1
        assert dlg.launch_freq_calls == 0

    def test_neither_available_warns(self, qapp, monkeypatch, tmp_path):
        dlg, warn_calls = self._run_with(qapp, monkeypatch, tmp_path, False, False)
        assert len(warn_calls) == 1
        assert dlg.launch_freq_calls == 0
        assert dlg.launch_mo_calls == 0

    def test_empty_path_returns_without_dialog(self, qapp, monkeypatch, tmp_path):
        from PyQt6.QtWidgets import QFileDialog

        monkeypatch.setattr(
            QFileDialog, "getOpenFileName", staticmethod(lambda *a, **k: ("", ""))
        )
        called = []
        monkeypatch.setattr(
            _fchk.QTimer,
            "singleShot",
            staticmethod(lambda ms, fn: called.append(fn)),
        )
        _fchk.PLUGIN_CONTEXT = _fchk_context()
        mw = _make_real_mw(qapp)
        _fchk.run(mw)
        assert called == []

    def test_no_plugin_context_returns_without_dialog(self, qapp, monkeypatch, tmp_path):
        from PyQt6.QtWidgets import QFileDialog

        fchk_file = tmp_path / "mol2.fchk"
        fchk_file.write_text("")
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: (str(fchk_file), "")),
        )
        called = []
        monkeypatch.setattr(
            _fchk.QTimer,
            "singleShot",
            staticmethod(lambda ms, fn: called.append(fn)),
        )
        _fchk.PLUGIN_CONTEXT = None
        mw = _make_real_mw(qapp)
        _fchk.run(mw)
        assert called == []
