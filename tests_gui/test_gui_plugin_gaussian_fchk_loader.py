"""
Headless GUI tests for the Gaussian FCHK Loader plugin.

Covers: FCHKLoaderDialog.
"""

from __future__ import annotations

from pathlib import Path
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
