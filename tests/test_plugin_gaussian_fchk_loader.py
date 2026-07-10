"""
Tests for the Gaussian FCHK Loader plugin (plugins/Gaussian_FCHK_Loader/gaussian_fchk_loader.py).

All heavy deps (PyQt6, rdkit, numpy) are stubbed via mock_optional_imports().
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from conftest import extract_function, load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
FCHK_LOADER_PATH = PLUGINS_DIR / "Gaussian_FCHK_Loader" / "gaussian_fchk_loader.py"

with mock_optional_imports():
    _fchk_loader = load_plugin(FCHK_LOADER_PATH)


class TestFindFileRecursive:
    def test_finds_exact_filename(self, tmp_path):
        target = tmp_path / "sub" / "deep"
        target.mkdir(parents=True)
        (target / "result.fchk").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "result.fchk")
        assert found is not None
        assert found.endswith("result.fchk")

    def test_fnmatch_wildcard_pattern(self, tmp_path):
        (tmp_path / "run01.fchk").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "*.fchk")
        assert found is not None
        assert found.endswith(".fchk")

    def test_returns_none_when_not_found(self, tmp_path):
        (tmp_path / "data.txt").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "*.fchk")
        assert found is None

    def test_returns_none_in_empty_directory(self, tmp_path):
        found = _fchk_loader.find_file_recursive(str(tmp_path), "anything.fchk")
        assert found is None

    def test_finds_first_match_in_nested_dirs(self, tmp_path):
        (tmp_path / "a").mkdir()
        (tmp_path / "a" / "mol.fchk").write_text("")
        (tmp_path / "b").mkdir()
        (tmp_path / "b" / "mol.fchk").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "mol.fchk")
        assert found is not None
        assert found.endswith("mol.fchk")

    def test_pattern_case_sensitive(self, tmp_path):
        (tmp_path / "Result.FCHK").write_text("")
        found = _fchk_loader.find_file_recursive(str(tmp_path), "result.fchk")
        # fnmatch on Windows is case-insensitive; on Linux it's case-sensitive.
        # Just verify the function doesn't raise.
        assert isinstance(found, (str, type(None)))


class TestFindMoAnalyzerModule:
    def test_finds_package_with_init(self, tmp_path):
        pkg = tmp_path / "plugins" / "gaussian_fchk_mo_analyzer"
        pkg.mkdir(parents=True)
        (pkg / "__init__.py").write_text("")
        result = _fchk_loader.find_mo_analyzer_module(str(tmp_path))
        assert result is not None
        assert result.endswith("gaussian_fchk_mo_analyzer")

    def test_no_init_py_returns_none(self, tmp_path):
        # Directory exists but has no __init__.py
        pkg = tmp_path / "gaussian_fchk_mo_analyzer"
        pkg.mkdir()
        result = _fchk_loader.find_mo_analyzer_module(str(tmp_path))
        assert result is None

    def test_no_package_returns_none(self, tmp_path):
        result = _fchk_loader.find_mo_analyzer_module(str(tmp_path))
        assert result is None

    def test_nested_package_found(self, tmp_path):
        pkg = tmp_path / "level1" / "level2" / "gaussian_fchk_mo_analyzer"
        pkg.mkdir(parents=True)
        (pkg / "__init__.py").write_text("")
        result = _fchk_loader.find_mo_analyzer_module(str(tmp_path))
        assert result is not None
        assert "gaussian_fchk_mo_analyzer" in result




# ---------------------------------------------------------------------------
# load_module_from_path + initialize contract
# ---------------------------------------------------------------------------

import sys

from conftest import make_context


class TestLoadModuleFromPath:
    def test_loads_valid_module(self, tmp_path):
        mod_file = tmp_path / "mini_plugin_mod.py"
        mod_file.write_text("ANSWER = 42\n")
        mod = _fchk_loader.load_module_from_path("mini_plugin_mod_t1", str(mod_file))
        assert mod is not None
        assert mod.ANSWER == 42

    def test_module_registered_in_sys_modules(self, tmp_path):
        mod_file = tmp_path / "mini_plugin_mod2.py"
        mod_file.write_text("X = 1\n")
        _fchk_loader.load_module_from_path("mini_plugin_mod_t2", str(mod_file))
        assert "mini_plugin_mod_t2" in sys.modules
        del sys.modules["mini_plugin_mod_t2"]

    def test_syntax_error_returns_none(self, tmp_path):
        mod_file = tmp_path / "broken_mod.py"
        mod_file.write_text("def broken(:\n")
        assert _fchk_loader.load_module_from_path("broken_mod_t", str(mod_file)) is None

    def test_missing_file_returns_none(self, tmp_path):
        missing = tmp_path / "ghost.py"
        assert _fchk_loader.load_module_from_path("ghost_mod_t", str(missing)) is None


class TestFCHKLoaderInitialize:
    def _init(self):
        with mock_optional_imports():
            mod = load_plugin(FCHK_LOADER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        return ctx

    def test_registers_three_extensions_priority_100(self):
        ctx = self._init()
        calls = ctx.register_file_opener.call_args_list
        exts = {c[0][0] for c in calls}
        assert exts == {".fchk", ".fck", ".fch"}
        assert all(c.kwargs.get("priority") == 100 for c in calls)

    def test_drop_handler_accepts_fchk_case_insensitive(self):
        ctx = self._init()
        handler = ctx.register_drop_handler.call_args[0][0]
        assert handler("C:/data/RESULT.FCHK") is True

    def test_drop_handler_rejects_other_files(self):
        ctx = self._init()
        handler = ctx.register_drop_handler.call_args[0][0]
        assert handler("C:/data/result.xyz") is False

    def test_drop_handler_priority_100(self):
        ctx = self._init()
        assert ctx.register_drop_handler.call_args.kwargs.get("priority") == 100


# ---------------------------------------------------------------------------
# open_fchk()'s deferred run_dialog() auto-launch decision logic
# ---------------------------------------------------------------------------


class _FakeFCHKDialog:
    """Stand-in for FCHKLoaderDialog controlling can_freq/can_mo."""

    can_freq = True
    can_mo = False

    def __init__(self, mw, context, path):
        self.mw, self.context, self.path = mw, context, path
        self.btn_freq = SimpleNamespace(isEnabled=lambda: _FakeFCHKDialog.can_freq)
        self.btn_mo = SimpleNamespace(isEnabled=lambda: _FakeFCHKDialog.can_mo)
        self.launch_freq = MagicMock()
        self.launch_mo = MagicMock()
        self.exec = MagicMock()


class TestOpenFchkAutoLaunch:
    def test_both_available_shows_dialog(self):
        captured = {}

        class _CapDlg(_FakeFCHKDialog):
            def __init__(self, mw, context, path):
                super().__init__(mw, context, path)
                captured["dlg"] = self

        with mock_optional_imports():
            mod = load_plugin(FCHK_LOADER_PATH)
            mod.FCHKLoaderDialog = _CapDlg
            _FakeFCHKDialog.can_freq = True
            _FakeFCHKDialog.can_mo = True
            ctx = make_context()
            mod.initialize(ctx)
            open_fchk = ctx.register_file_opener.call_args_list[0][0][1]
            open_fchk("C:/data/out.fchk")
            run_dialog = mod.QTimer.singleShot.call_args[0][1]
            run_dialog()

        captured["dlg"].exec.assert_called_once()
        captured["dlg"].launch_freq.assert_not_called()
        captured["dlg"].launch_mo.assert_not_called()

    def test_only_freq_available_calls_launch_freq(self):
        captured = {}

        class _CapDlg(_FakeFCHKDialog):
            def __init__(self, mw, context, path):
                super().__init__(mw, context, path)
                captured["dlg"] = self

        with mock_optional_imports():
            mod = load_plugin(FCHK_LOADER_PATH)
            mod.FCHKLoaderDialog = _CapDlg
            _FakeFCHKDialog.can_freq = True
            _FakeFCHKDialog.can_mo = False
            ctx = make_context()
            mod.initialize(ctx)
            open_fchk = ctx.register_file_opener.call_args_list[0][0][1]
            open_fchk("C:/data/out.fchk")
            run_dialog = mod.QTimer.singleShot.call_args[0][1]
            run_dialog()

        captured["dlg"].launch_freq.assert_called_once()
        captured["dlg"].launch_mo.assert_not_called()
        captured["dlg"].exec.assert_not_called()

    def test_only_mo_available_calls_launch_mo(self):
        captured = {}

        class _CapDlg(_FakeFCHKDialog):
            def __init__(self, mw, context, path):
                super().__init__(mw, context, path)
                captured["dlg"] = self

        with mock_optional_imports():
            mod = load_plugin(FCHK_LOADER_PATH)
            mod.FCHKLoaderDialog = _CapDlg
            _FakeFCHKDialog.can_freq = False
            _FakeFCHKDialog.can_mo = True
            ctx = make_context()
            mod.initialize(ctx)
            open_fchk = ctx.register_file_opener.call_args_list[0][0][1]
            open_fchk("C:/data/out.fchk")
            run_dialog = mod.QTimer.singleShot.call_args[0][1]
            run_dialog()

        captured["dlg"].launch_mo.assert_called_once()
        captured["dlg"].launch_freq.assert_not_called()

    def test_neither_available_shows_warning(self):
        with mock_optional_imports():
            mod = load_plugin(FCHK_LOADER_PATH)
            mod.FCHKLoaderDialog = _FakeFCHKDialog
            _FakeFCHKDialog.can_freq = False
            _FakeFCHKDialog.can_mo = False
            mod.QMessageBox.warning = MagicMock()
            ctx = make_context()
            mod.initialize(ctx)
            open_fchk = ctx.register_file_opener.call_args_list[0][0][1]
            open_fchk("C:/data/out.fchk")
            run_dialog = mod.QTimer.singleShot.call_args[0][1]
            run_dialog()
            mod.QMessageBox.warning.assert_called_once()


# ---------------------------------------------------------------------------
# run(mw) — legacy entry point guards
# ---------------------------------------------------------------------------


class TestRunEntryPoint:
    def test_empty_path_returns_without_queuing(self):
        with mock_optional_imports():
            mod = load_plugin(FCHK_LOADER_PATH)
            mod.PLUGIN_CONTEXT = make_context()
            from PyQt6.QtWidgets import QFileDialog

            QFileDialog.getOpenFileName = MagicMock(return_value=("", ""))
            mod.run(MagicMock())
            mod.QTimer.singleShot.assert_not_called()

    def test_no_plugin_context_returns_without_queuing(self):
        with mock_optional_imports():
            mod = load_plugin(FCHK_LOADER_PATH)
            assert mod.PLUGIN_CONTEXT is None
            from PyQt6.QtWidgets import QFileDialog

            QFileDialog.getOpenFileName = MagicMock(return_value=("C:/x.fchk", ""))
            mod.run(MagicMock())
            mod.QTimer.singleShot.assert_not_called()

    def test_valid_path_and_context_queues_run_dialog(self):
        with mock_optional_imports():
            mod = load_plugin(FCHK_LOADER_PATH)
            mod.PLUGIN_CONTEXT = make_context()
            mod.FCHKLoaderDialog = _FakeFCHKDialog
            _FakeFCHKDialog.can_freq = False
            _FakeFCHKDialog.can_mo = False
            mod.QMessageBox.warning = MagicMock()
            from PyQt6.QtWidgets import QFileDialog

            QFileDialog.getOpenFileName = MagicMock(return_value=("C:/x.fchk", ""))
            mod.run(MagicMock())
            mod.QTimer.singleShot.assert_called_once()


# ---------------------------------------------------------------------------
# FCHKLoaderDialog.launch_freq / launch_mo / close_existing_analyzers
# ---------------------------------------------------------------------------


def _launch_freq_fn():
    globs = {
        "load_module_from_path": lambda name, path: SimpleNamespace(
            GaussianFCHKFreqAnalyzer=lambda context, dock_widget: SimpleNamespace(
                load_file=MagicMock()
            )
        ),
        "QDockWidget": lambda title, mw: SimpleNamespace(
            setAllowedAreas=MagicMock(),
            setAttribute=MagicMock(),
            setWidget=MagicMock(),
            show=MagicMock(),
        ),
        "Qt": SimpleNamespace(
            DockWidgetArea=SimpleNamespace(
                AllDockWidgetAreas=1, RightDockWidgetArea=2
            ),
            WidgetAttribute=SimpleNamespace(WA_DeleteOnClose=3),
        ),
        "QMessageBox": MagicMock(),
    }
    return extract_function(
        FCHK_LOADER_PATH, "FCHKLoaderDialog", "launch_freq", globs
    ), globs


def _make_dialog_self(freq_path="C:/freq.py", mo_pkg=None):
    mw = MagicMock()
    return SimpleNamespace(
        freq_analyzer_path=freq_path,
        mo_analyzer_pkg=mo_pkg,
        btn_freq=SimpleNamespace(setEnabled=MagicMock()),
        btn_mo=SimpleNamespace(setEnabled=MagicMock()),
        close_existing_analyzers=MagicMock(),
        context=MagicMock(),
        mw=mw,
        fchk_path="C:/data/mol.fchk",
        accept=MagicMock(),
    )


class TestLaunchFreq:
    def test_missing_path_is_noop(self):
        fn, globs = _launch_freq_fn()
        self_ = _make_dialog_self(freq_path=None)
        fn(self_)
        self_.close_existing_analyzers.assert_not_called()

    def test_success_calls_load_file_and_accepts(self):
        fn, globs = _launch_freq_fn()
        self_ = _make_dialog_self()
        fn(self_)
        self_.close_existing_analyzers.assert_called_once()
        self_.accept.assert_called_once()

    def test_exception_reenables_buttons_and_shows_error(self):
        globs = {
            "load_module_from_path": lambda name, path: (_ for _ in ()).throw(
                RuntimeError("boom")
            ),
            "QDockWidget": MagicMock(),
            "Qt": SimpleNamespace(
                DockWidgetArea=SimpleNamespace(
                    AllDockWidgetAreas=1, RightDockWidgetArea=2
                ),
                WidgetAttribute=SimpleNamespace(WA_DeleteOnClose=3),
            ),
            "QMessageBox": MagicMock(),
        }
        fn = extract_function(FCHK_LOADER_PATH, "FCHKLoaderDialog", "launch_freq", globs)
        self_ = _make_dialog_self()
        fn(self_)
        self_.btn_freq.setEnabled.assert_called_with(True)
        self_.btn_mo.setEnabled.assert_called_with(True)
        globs["QMessageBox"].critical.assert_called_once()
        self_.accept.assert_not_called()


def _close_existing_analyzers_fn():
    globs = {"QDockWidget": object}
    return extract_function(
        FCHK_LOADER_PATH, "FCHKLoaderDialog", "close_existing_analyzers", globs
    )


class TestCloseExistingAnalyzers:
    def test_closes_matching_docks_only(self):
        fn = _close_existing_analyzers_fn()
        freq_dock = SimpleNamespace(
            windowTitle=lambda: "Gaussian Freq Analyzer",
            close=MagicMock(),
            deleteLater=MagicMock(),
        )
        mo_dock = SimpleNamespace(
            windowTitle=lambda: "Gaussian MO Analyzer",
            close=MagicMock(),
            deleteLater=MagicMock(),
        )
        other_dock = SimpleNamespace(
            windowTitle=lambda: "Something Else",
            close=MagicMock(),
            deleteLater=MagicMock(),
        )
        docks = [freq_dock, mo_dock, other_dock]
        mw = SimpleNamespace(findChildren=lambda cls: docks)
        self_ = SimpleNamespace(mw=mw)
        fn(self_)
        freq_dock.close.assert_called_once()
        mo_dock.close.assert_called_once()
        other_dock.close.assert_not_called()
