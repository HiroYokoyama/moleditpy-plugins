"""
GUI tests for the OpenBabel Conversion Tool plugin.

Covers: MoleculeSelectionDialog — openbabel/pybel are mocked by
mock_chemistry_imports(), so OBABEL_AVAILABLE=True and pybel is a MagicMock.
The dialog only needs a list of mock "molecule" objects with .title / .formula.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports, patch_mock_attr

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

OPENBABEL_PATH = PLUGINS_DIR / "OpenBabel_Conversion_Tool" / "openbabel_conversion_tool.py"

with mock_chemistry_imports():
    _openbabel = load_plugin_for_gui(OPENBABEL_PATH)

# Native format discovery is exercised separately; dialog tests use the plain
# metadata the discovery worker returns rather than importing native Open Babel.
_openbabel.OBABEL_AVAILABLE = True


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_mock_mols(n: int = 3) -> list:
    mols = []
    for i in range(n):
        m = MagicMock()
        m.title = f"Molecule {i + 1}"
        m.formula = f"C{i + 1}H{2 * (i + 1)}"
        mols.append(m)
    return mols


# ===========================================================================
# MoleculeSelectionDialog  (OpenBabel Conversion Tool)
# ===========================================================================


class TestMoleculeSelectionDialog:
    """MoleculeSelectionDialog with mock pybel molecule objects."""

    @pytest.fixture
    def dlg(self, qapp):
        mols = _make_mock_mols(3)
        d = _openbabel.MoleculeSelectionDialog(mols, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Select Molecule"

    def test_list_widget_has_correct_item_count(self, dlg):
        assert dlg.list_widget.count() == 3

    def test_first_row_selected_by_default(self, dlg):
        assert dlg.list_widget.currentRow() == 0

    def test_get_selected_index_returns_zero_initially(self, dlg):
        assert dlg.get_selected_index() == 0

    def test_get_selected_index_after_selection_change(self, dlg):
        dlg.list_widget.setCurrentRow(2)
        assert dlg.get_selected_index() == 2

    def test_single_selection_mode(self, dlg):
        from PySide6.QtWidgets import QAbstractItemView
        assert dlg.list_widget.selectionMode() == QAbstractItemView.SelectionMode.SingleSelection

    def test_empty_mol_list_creates_empty_widget(self, qapp):
        d = _openbabel.MoleculeSelectionDialog([], parent=None)
        assert d.list_widget.count() == 0
        d.destroy()

    def test_get_selected_index_returns_none_when_nothing_selected(self, qapp):
        d = _openbabel.MoleculeSelectionDialog([], parent=None)
        d.list_widget.clearSelection()
        assert d.get_selected_index() is None
        d.destroy()

    def test_single_mol_list(self, qapp):
        d = _openbabel.MoleculeSelectionDialog(_make_mock_mols(1), parent=None)
        assert d.list_widget.count() == 1
        d.destroy()

    def test_obabel_available_flag(self):
        # Under mock_chemistry_imports, openbabel is mocked → OBABEL_AVAILABLE=True
        assert _openbabel.OBABEL_AVAILABLE is True


class TestIsolatedConversion:
    def test_format_discovery_worker_returns_plain_metadata(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_openbabel, 'OBABEL_AVAILABLE', None)
        # Exercise the worker whether or not openbabel is installed here, and
        # never read or write a format cache beside the plugin source.
        monkeypatch.setattr(_openbabel, '_openbabel_fingerprint', lambda: {'v': 1})
        monkeypatch.setattr(_openbabel, '_formats_cache_path', lambda: str(tmp_path / 'formats.json'))
        monkeypatch.setattr(_openbabel, '_FORMAT_SCRIPT', """
import json, sys
with open(sys.argv[1], 'w', encoding='utf-8') as handle:
    json.dump({'informats': {'pdb': 'PDB'}, 'outformats': {'xyz': 'XYZ'}}, handle)
"""
        )
        assert _openbabel._ensure_formats()
        assert _openbabel.pybel.informats == {'pdb': 'PDB'}
        assert _openbabel.pybel.outformats == {'xyz': 'XYZ'}

    def test_format_discovery_native_crash_does_not_crash_gui(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_openbabel, 'OBABEL_AVAILABLE', None)
        # Exercise the worker whether or not openbabel is installed here, and
        # never read or write a format cache beside the plugin source.
        monkeypatch.setattr(_openbabel, '_openbabel_fingerprint', lambda: {'v': 1})
        monkeypatch.setattr(_openbabel, '_formats_cache_path', lambda: str(tmp_path / 'formats.json'))
        monkeypatch.setattr(_openbabel, '_FORMAT_SCRIPT', 'import os; os._exit(139)')
        assert _openbabel._ensure_formats() is False

    def test_process_success_and_failure(self, qapp):
        _openbabel._run_conversion('print("done")', [], None)
        with pytest.raises(ValueError, match='failed'):
            _openbabel._run_conversion('raise ValueError("unsupported data")', [], None)

    def test_timeout_keeps_gui_responsive(self, qapp, monkeypatch):
        from PySide6.QtCore import QTimer
        monkeypatch.setattr(_openbabel, 'CONVERSION_TIMEOUT_MS', 250)
        ticks = []
        timer = QTimer()
        timer.timeout.connect(lambda: ticks.append(1))
        timer.start(10)
        try:
            with pytest.raises(ValueError, match='timed out'):
                _openbabel._run_conversion('import time; time.sleep(60)', [], None)
        finally:
            timer.stop()
        assert len(ticks) >= 2

    def test_cancel_terminates_native_process(self, qapp):
        from PySide6.QtCore import QTimer
        from PySide6.QtWidgets import QApplication, QProgressDialog
        def cancel():
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, QProgressDialog):
                    widget.canceled.emit()
        QTimer.singleShot(100, cancel)
        with pytest.raises(ValueError, match='cancelled'):
            _openbabel._run_conversion('import time; time.sleep(60)', [], None)

    def test_unsupported_or_large_file_never_launches_parser(self, qapp, tmp_path, monkeypatch):
        monkeypatch.setattr(_openbabel.pybel, 'informats', {'xyz': 'XYZ'})
        runner = MagicMock()
        monkeypatch.setattr(_openbabel, '_run_conversion', runner)
        source = tmp_path / 'large.xyz'
        source.write_bytes(b'123456')
        monkeypatch.setattr(_openbabel, 'MAX_FILE_BYTES', 5)
        with pytest.raises(ValueError, match='not supported'):
            _openbabel._read_molecules(str(source), 'bad', None)
        with pytest.raises(ValueError, match='exceeds'):
            _openbabel._read_molecules(str(source), 'xyz', None)
        runner.assert_not_called()

    def test_unicode_source_copy_and_result_cleanup(self, qapp, tmp_path, monkeypatch):
        import json
        source = tmp_path / '水.xyz'
        source.write_text('1\nwater\nO 0 0 0\n', encoding='utf-8')
        paths = []
        def convert(script, args, parent):
            paths.extend(args[2:])
            assert Path(args[2]).read_bytes() == source.read_bytes()
            Path(args[3]).write_text(json.dumps([dict(title='water', formula='O', block='MOL')]))
        monkeypatch.setattr(_openbabel.pybel, 'informats', {'xyz': 'XYZ'})
        monkeypatch.setattr(_openbabel, '_run_conversion', convert)
        mols = _openbabel._read_molecules(str(source), 'xyz', None)
        assert mols[0].write('mol') == 'MOL'
        assert mols[0].title == 'water'
        assert all(not Path(path).exists() for path in paths)

    def test_child_atom_and_record_limits(self, qapp, tmp_path, monkeypatch):
        import sys
        from types import SimpleNamespace
        mol = MagicMock()
        mol.OBMol.NumAtoms.return_value = 10001
        monkeypatch.setitem(sys.modules, 'openbabel', SimpleNamespace(pybel=SimpleNamespace(readfile=lambda *a: iter([mol]))))
        monkeypatch.setattr(sys, 'argv', ['worker', 'import', 'xyz', 'input', str(tmp_path / 'out')])
        with pytest.raises(ValueError, match='10000 atoms'):
            exec(_openbabel._CONVERSION_SCRIPT, {})
        mol.OBMol.NumAtoms.return_value = 1
        mol.write.return_value = 'MOL'
        mol.title = 'sample'
        mol.formula = 'H'
        monkeypatch.setitem(sys.modules, 'openbabel', SimpleNamespace(pybel=SimpleNamespace(readfile=lambda *a: iter([mol] * 501))))
        with pytest.raises(ValueError, match='500 molecules'):
            exec(_openbabel._CONVERSION_SCRIPT, {})

    def test_failed_export_preserves_existing_file(self, qapp, tmp_path, monkeypatch):
        target = tmp_path / 'existing.xyz'
        target.write_bytes(b'original')
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        ctx.current_molecule.GetNumAtoms.return_value = 2
        monkeypatch.setattr(_openbabel.pybel, 'outformats', {'xyz': 'XYZ'})
        patch_mock_attr(monkeypatch, _openbabel, "Chem", 'MolToMolBlock', lambda *a: 'MOL')
        monkeypatch.setattr(_openbabel.QFileDialog, 'getSaveFileName', lambda *a: (str(target), ''))
        monkeypatch.setattr(_openbabel, '_run_conversion', MagicMock(side_effect=ValueError('timed out')))
        critical = MagicMock()
        monkeypatch.setattr(_openbabel.QMessageBox, 'critical', critical)
        _openbabel.export_with_openbabel(ctx)
        assert target.read_bytes() == b'original'
        critical.assert_called_once()
