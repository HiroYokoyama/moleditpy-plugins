"""
Headless GUI tests for the Gaussian MO Analyzer plugin (package plugin).

Covers: initialize(), OrbitalWidget.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MO_PKG_DIR = PLUGINS_DIR / "Gaussian_MO_Analyzer"

with mock_chemistry_imports():
    sys.path.insert(0, str(MO_PKG_DIR))
    try:
        _mo = importlib.import_module("gaussian_fchk_mo_analyzer")
    finally:
        sys.path.remove(str(MO_PKG_DIR))


# ===========================================================================
# Gaussian MO Analyzer  (package plugin: OrbitalWidget + initialize)
# ===========================================================================


class TestMOAnalyzerInitialize:
    def test_registers_fchk_openers_and_drop_handler(self, qapp):
        ctx = MagicMock()
        _mo.initialize(ctx)
        exts = [c.args[0] for c in ctx.register_file_opener.call_args_list]
        assert exts == [".fchk", ".fck", ".fch"]
        ctx.register_drop_handler.assert_called_once()

    def test_drop_handler_accepts_only_fchk(self, qapp):
        ctx = MagicMock()
        _mo.initialize(ctx)
        handler = ctx.register_drop_handler.call_args.args[0]
        with patch.object(_mo, "OrbitalWidget") as ow:
            assert handler("mol.fchk") is True
            assert handler("mol.xyz") is False
            ow.assert_called_once()


class TestOrbitalWidget:
    @pytest.fixture
    def win(self, qapp, tmp_path):
        # Empty fchk → reader.data empty; basis prep patched out so the
        # widget builds with n_basis=0 (the "Error reading MOs" path).
        fchk = tmp_path / "water.fchk"
        fchk.write_text("dummy fchk\n", encoding="utf-8")
        ctx = MagicMock()
        with patch.object(
            _mo.analyzer.BasisSetEngine, "_prepare_basis_set", lambda self: None
        ):
            w = _mo.OrbitalWidget(None, ctx, str(fchk))
        yield w
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title_includes_filename(self, win):
        assert "Gaussian MO Analyzer" in win.windowTitle()
        assert "water.fchk" in win.windowTitle()

    def test_default_grid_points(self, win):
        assert win.spin_points.value() == 40

    def test_default_margin(self, win):
        assert win.spin_gap.value() == 4.0

    def test_default_isovalue(self, win):
        assert win.spin_iso.value() == 0.02

    def test_default_opacity(self, win):
        assert win.spin_opacity.value() == 0.4

    def test_smooth_shading_checked_by_default(self, win):
        assert win.check_smooth.isChecked()

    def test_progress_bar_hidden_initially(self, win):
        assert win.pbar.isHidden()

    def test_list_shows_error_without_mo_data(self, win):
        assert win.list_widget.count() == 1
        assert win.list_widget.item(0).text() == "Error reading MOs"

    def test_generate_with_no_selection_keeps_button_enabled(self, win):
        win.list_widget.clearSelection()
        win.generate_selected()
        assert win.btn_gen.isEnabled()

    def test_iso_change_without_cube_is_noop(self, win):
        win.last_cube_path = None
        win.on_iso_changed(0.05)  # must not raise

    def test_cube_path_layout(self, win):
        path = win.get_cube_path(7)
        assert path.endswith("water_MO7.cube")
        assert "water_cubes" in path
