"""
Shared regression test: all five legacy QM input generator plugins
(MOPAC, GAMESS, PySCF, Psi4, NWChem) must warn and return early from
run(mw) when no molecule is loaded, without constructing their setup
dialog.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MOPAC_PATH = PLUGINS_DIR / "Mopac_Input_Generator" / "mopac_input_generator.py"
GAMESS_PATH = PLUGINS_DIR / "Gamess_Input_Generator" / "gamess_input_generator.py"
PYSCF_PATH = PLUGINS_DIR / "PySCF_Input_Generator" / "pyscf_input_generator.py"
PSI4_PATH = PLUGINS_DIR / "Psi4_Input_Generator" / "psi4_input_generator.py"
NWCHEM_PATH = PLUGINS_DIR / "NWChem_Input_Generator" / "nwchem_input_generator.py"


class TestRunNoMoleculeGuard:
    @pytest.mark.parametrize(
        "path",
        [MOPAC_PATH, GAMESS_PATH, PYSCF_PATH, PSI4_PATH, NWCHEM_PATH],
        ids=lambda p: p.stem,
    )
    def test_run_warns_and_returns_without_molecule(self, path):
        with mock_optional_imports():
            mod = load_plugin(path)
            mw = MagicMock()
            mw.current_mol = None
            mod.QMessageBox.warning.reset_mock()
            mod.run(mw)
            mod.QMessageBox.warning.assert_called_once()
