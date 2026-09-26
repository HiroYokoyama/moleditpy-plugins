"""
Headless GUI tests for the Elemental Analysis plugin.

Covers: ElementalAnalysisDialog, CopyableTable.

Real PyQt6 (QT_QPA_PLATFORM=offscreen); rdkit is replaced via
mock_chemistry_imports() and a small fake periodic table is patched in so
the real recalc / check_update bodies run to completion.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QKeyEvent, QKeySequence
from PyQt6.QtWidgets import QApplication, QTableWidgetSelectionRange, QWidget

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGIN_PATH = (
    Path(__file__).resolve().parents[1]
    / "plugins"
    / "Elemental_Analysis"
    / "elemental_analysis.py"
)

with mock_chemistry_imports():
    ea = load_plugin_for_gui(PLUGIN_PATH)


_WEIGHTS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.067}
_ISOTOPES = {1: 1.00782503207, 6: 12.0, 7: 14.0030740048, 8: 15.99491461956, 16: 31.97207100}
_ANUM = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16}


class _FakePT:
    def GetAtomicWeight(self, sym):
        if sym not in _WEIGHTS:
            raise RuntimeError(f"Element '{sym}' not found")
        return _WEIGHTS[sym]

    def GetAtomicNumber(self, sym):
        return _ANUM[sym]

    def GetMostCommonIsotope(self, anum):
        return anum

    def GetMassForIsotope(self, anum, _mass_number):
        return _ISOTOPES[anum]


class _FakeChem:
    @staticmethod
    def GetPeriodicTable():
        return _FakePT()

    @staticmethod
    def AddHs(mol):
        return mol


class _FakeRDMolDescriptors:
    @staticmethod
    def CalcMolFormula(mol):
        return mol.formula


@pytest.fixture
def fake_chem(monkeypatch):
    monkeypatch.setattr(ea, "Chem", _FakeChem)
    monkeypatch.setattr(ea, "rdMolDescriptors", _FakeRDMolDescriptors)
    monkeypatch.setattr(ea, "Draw", None)


def _context(mol=None, main_window=None):
    ctx = MagicMock()
    ctx.get_main_window.return_value = main_window
    ctx.current_molecule = mol
    ctx.scene = None
    return ctx


@pytest.fixture
def dlg(qapp, fake_chem):
    d = ea.ElementalAnalysisDialog(_context())
    yield d
    d.done(0)
    d.deleteLater()


class TestDialogBasics:
    def test_title_and_widgets(self, dlg):
        assert dlg.windowTitle() == "Elemental Analysis"
        assert dlg.btn_export_csv.text() == "Export CSV"
        assert dlg.btn_report.text() == "PDF Report"
        assert dlg.table.columnCount() == 3
        assert dlg.solvent_combo.count() == len(ea.SOLVATES)

    def test_sync_disabled_without_main_window(self, dlg):
        assert not dlg.sync_check.isEnabled()
        assert not dlg.use_2d_check.isChecked()

    def test_registers_window(self, qapp, fake_chem):
        ctx = _context()
        d = ea.ElementalAnalysisDialog(ctx)
        ctx.register_window.assert_called_once_with("main_panel", d)
        d.done(0)

    def test_empty_formula_shows_dashes(self, dlg):
        assert dlg.table.rowCount() == 0
        assert dlg.lbl_mw.text().endswith("-")
        assert dlg.lbl_report.text() == ""


class TestRecalc:
    def test_caffeine(self, dlg):
        dlg.formula_input.setText("C8H10N4O2")
        assert dlg.table.rowCount() == 4
        assert [dlg.table.item(r, 0).text() for r in range(4)] == ["C", "H", "N", "O"]
        assert dlg.table.item(0, 2).text() == "49.48"
        assert "194.1940" in dlg.lbl_mw.text()
        assert "194.0804" in dlg.lbl_em.text()
        assert dlg.lbl_report.text() == "Anal. Calcd for C8H10N4O2: C, 49.48; H, 5.19; N, 28.85."

    def test_solvate_changes_composition_but_not_exact_mass(self, dlg):
        dlg.formula_input.setText("C8H10N4O2")
        dlg.solvent_combo.setCurrentIndex(0)  # H2O
        dlg.solvate_spin.setValue(0.5)
        assert dlg.total_formula == "C8H10N4O2·0.5H2O"
        assert "203.2015" in dlg.lbl_mw.text()
        assert "194.0804" in dlg.lbl_em.text()
        assert dlg.table.item(1, 1).text() == "11"
        assert dlg.lbl_report.text().startswith("Anal. Calcd for C8H10N4O2·0.5H2O: C, 47.29")

    def test_one_equivalent_has_no_coefficient(self, dlg):
        dlg.solvate_spin.setValue(1.0)
        assert dlg.solvate_suffix() == "·H2O"
        dlg.solvate_spin.setValue(0.0)
        assert dlg.solvate_suffix() == ""

    def test_invalid_formula_clears(self, dlg):
        dlg.formula_input.setText("C6H6")
        dlg.formula_input.setText("C6?H6")
        assert dlg.table.rowCount() == 0
        assert dlg.lbl_report.text() == ""

    def test_unknown_element_clears(self, dlg):
        dlg.formula_input.setText("Xx2")
        assert dlg.rows == []
        assert dlg.lbl_mw.text().endswith("-")

    def test_formula_input_is_normalised_to_hill(self, dlg):
        dlg.formula_input.setText("OHCH3")
        assert dlg.total_formula == "CH4O"


class TestCopy:
    def test_copy_report_line(self, dlg):
        dlg.formula_input.setText("C8H10N4O2")
        dlg.copy_report_line()
        assert QApplication.clipboard().text() == dlg.lbl_report.text()

    def test_copy_table(self, dlg):
        dlg.formula_input.setText("C8H10N4O2")
        dlg.copy_table()
        assert QApplication.clipboard().text() == ea.rows_to_tsv(dlg.rows)

    def test_copy_nothing_when_empty(self, dlg):
        QApplication.clipboard().setText("keep")
        dlg.copy_report_line()
        dlg.copy_table()
        dlg.table.copy_selection()
        assert QApplication.clipboard().text() == "keep"

    def test_ctrl_c_copies_selected_cells(self, dlg):
        dlg.formula_input.setText("C8H10N4O2")
        dlg.table.setRangeSelected(QTableWidgetSelectionRange(0, 0, 1, 2), True)
        event = QKeyEvent(
            QKeyEvent.Type.KeyPress, Qt.Key.Key_C, Qt.KeyboardModifier.ControlModifier
        )
        assert event.matches(QKeySequence.StandardKey.Copy)
        dlg.table.keyPressEvent(event)
        assert QApplication.clipboard().text() == "C\t8\t49.48\nH\t10\t5.19"

    def test_other_keys_pass_through(self, dlg):
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Down, Qt.KeyboardModifier.NoModifier)
        dlg.table.keyPressEvent(event)  # must not raise


class TestSync:
    def test_check_update_from_2d_state(self, qapp, fake_chem):
        mol = SimpleNamespace(formula="C2H6O")
        mw = QWidget()
        mw.state_manager = SimpleNamespace(
            data=SimpleNamespace(to_rdkit_mol=lambda: mol)
        )
        d = ea.ElementalAnalysisDialog(_context(main_window=mw))
        try:
            assert d.sync_check.isChecked()
            assert d.timer.isActive()
            assert d.formula_input.text() == "C2H6O"
            assert d.total_formula == "C2H6O"
            d.sync_check.setChecked(False)
            assert not d.timer.isActive()
        finally:
            d.done(0)

    def test_manual_edit_turns_sync_off(self, qapp, fake_chem):
        mol = SimpleNamespace(formula="CH4")
        mw = QWidget()
        d = ea.ElementalAnalysisDialog(_context(mol=mol, main_window=mw))
        try:
            assert d.sync_check.isChecked()
            d.formula_input.textEdited.emit("C2H6")
            assert not d.sync_check.isChecked()
        finally:
            d.done(0)

    def test_done_stops_timer(self, qapp, fake_chem):
        mw = QWidget()
        d = ea.ElementalAnalysisDialog(_context(main_window=mw))
        assert d.timer.isActive()
        d.done(0)
        assert not d.timer.isActive()

    def test_check_update_without_rdkit_is_noop(self, dlg, monkeypatch):
        monkeypatch.setattr(ea, "Chem", None)
        dlg.check_update()


class TestExports:
    def test_export_csv(self, dlg, tmp_path, monkeypatch):
        dlg.formula_input.setText("C8H10N4O2")
        out = tmp_path / "ea.csv"
        monkeypatch.setattr(ea.QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), ""))
        monkeypatch.setattr(ea.QMessageBox, "information", MagicMock())
        dlg.export_csv()
        lines = out.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "# Formula,C8H10N4O2"
        assert lines[2].startswith("C,8,96.08800,49.4")

    def test_export_csv_cancelled(self, dlg, monkeypatch):
        monkeypatch.setattr(ea.QFileDialog, "getSaveFileName", lambda *a, **k: ("", ""))
        dlg.export_csv()

    def test_export_csv_error(self, dlg, tmp_path, monkeypatch):
        crit = MagicMock()
        monkeypatch.setattr(
            ea.QFileDialog, "getSaveFileName", lambda *a, **k: (str(tmp_path), "")
        )
        monkeypatch.setattr(ea.QMessageBox, "critical", crit)
        dlg.export_csv()  # a directory cannot be opened for writing
        crit.assert_called_once()

    def test_pdf_report(self, dlg, tmp_path, monkeypatch):
        dlg.formula_input.setText("C8H10N4O2")
        out = tmp_path / "report.pdf"
        info, crit = MagicMock(), MagicMock()
        monkeypatch.setattr(ea.QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), ""))
        monkeypatch.setattr(ea.QMessageBox, "information", info)
        monkeypatch.setattr(ea.QMessageBox, "critical", crit)
        dlg.create_report()
        crit.assert_not_called()
        info.assert_called_once()
        assert out.stat().st_size > 1000

    def test_report_html_contents(self, dlg):
        dlg.formula_input.setText("C8H10N4O2")
        html = dlg.build_report_html("", 0, 0)
        assert "Elemental Analysis Report" in html
        assert "49.48" in html
        assert "No Structure" in html
        assert "Anal. Calcd for C8H10N4O2" in html

    def test_molecule_image_without_mol(self, dlg):
        dlg.mol = None
        assert dlg._molecule_image_b64() == ("", 0, 0)


class TestInitialize:
    def test_toggle_creates_then_reuses(self, qapp, fake_chem, monkeypatch):
        ctx = _context()
        windows = {}
        ctx.register_window.side_effect = windows.__setitem__
        ctx.get_window.side_effect = windows.get
        ea.initialize(ctx)
        toggle = ctx.add_analysis_tool.call_args[0][1]
        toggle()
        first = windows["main_panel"]
        toggle()
        assert windows["main_panel"] is first
        first.done(0)

    def test_toggle_without_rdkit(self, monkeypatch):
        monkeypatch.setattr(ea, "Chem", None)
        ctx = _context()
        ea.initialize(ctx)
        ctx.add_analysis_tool.call_args[0][1]()
        ctx.show_status_message.assert_called_once()
