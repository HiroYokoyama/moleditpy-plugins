"""Use the real engine while forbidding its import in the GUI process."""

import os
import subprocess
import sys
from pathlib import Path


def test_real_openbabel_import_and_export_are_isolated(tmp_path):
    plugin = (
        Path(__file__).resolve().parents[1]
        / "plugins/OpenBabel_Conversion_Tool/openbabel_conversion_tool.py"
    )
    script = r'''
import importlib.abc
import importlib.util
import sys
from pathlib import Path
from PySide6.QtWidgets import QApplication
from rdkit import Chem

class ForbidEngineImport(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "openbabel" or fullname.startswith("openbabel.") or fullname == "pybel":
            raise AssertionError("The GUI process must not import native Open Babel")

sys.meta_path.insert(0, ForbidEngineImport())
app = QApplication([])
spec = importlib.util.spec_from_file_location("conversion_tool", sys.argv[1])
plugin = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plugin)
assert plugin.OBABEL_AVAILABLE is None
assert plugin._ensure_formats(), "Real engine format discovery failed"
assert "smi" in plugin.pybel.informats
directory = Path(sys.argv[2])
source = directory / "benzene.smi"
source.write_text("c1ccccc1 benzene\n", encoding="utf-8")
records = plugin._read_molecules(str(source), "smi", None)
assert len(records) == 1 and records[0].formula == "C6H6"
block = records[0].write("mol")
assert Chem.MolFromMolBlock(block).GetNumAtoms() == 6
molfile = directory / "benzene.mol"
molfile.write_text(block, encoding="utf-8")
output = directory / "export.smi"
plugin._run_conversion(plugin._CONVERSION_SCRIPT, ["export", "smi", str(molfile), str(output)], None)
assert Chem.MolFromSmiles(output.read_text().split()[0]).GetNumAtoms() == 6
assert not any(name == "openbabel" or name.startswith("openbabel.") for name in sys.modules)
'''
    environment = dict(os.environ, QT_QPA_PLATFORM="offscreen", QT_API="pyside6")
    result = subprocess.run(
        [sys.executable, "-c", script, str(plugin), str(tmp_path)],
        env=environment,
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
