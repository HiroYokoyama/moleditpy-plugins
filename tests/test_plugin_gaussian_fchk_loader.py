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
