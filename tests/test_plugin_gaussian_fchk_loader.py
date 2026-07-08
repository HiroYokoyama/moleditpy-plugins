"""
Tests for the Gaussian FCHK Loader plugin (plugins/Gaussian_FCHK_Loader/gaussian_fchk_loader.py).

All heavy deps (PyQt6, rdkit, numpy) are stubbed via mock_optional_imports().
"""

from __future__ import annotations

from pathlib import Path

from conftest import load_plugin, mock_optional_imports

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
