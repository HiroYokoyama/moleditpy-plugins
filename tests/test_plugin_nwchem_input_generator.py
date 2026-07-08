"""
Tests for the NWChem Input Generator plugin: scf-block spin-state lines.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
NWCHEM_PATH = PLUGINS_DIR / "NWChem_Input_Generator" / "nwchem_input_generator.py"


class TestNwchemScfSpinLines:
    def _fn(self):
        with mock_optional_imports():
            return load_plugin(NWCHEM_PATH)._scf_spin_lines

    def test_singlet_uses_rhf(self):
        assert self._fn()(1) == ["  rhf", "  singlet"]

    @pytest.mark.parametrize(
        "mult,keyword",
        [
            (2, "doublet"),
            (3, "triplet"),
            (4, "quartet"),
            (5, "quintet"),
            (6, "sextet"),
            (7, "septet"),
            (8, "octet"),
        ],
    )
    def test_open_shell_uses_uhf_keyword(self, mult, keyword):
        assert self._fn()(mult) == ["  uhf", f"  {keyword}"]

    def test_unnamed_high_multiplicity_falls_back_to_comment(self):
        lines = self._fn()(9)
        assert lines[0] == "  uhf"
        assert "9" in lines[1] and lines[1].lstrip().startswith("#")
