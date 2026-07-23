"""
Headless GUI tests for the MOPAC Input Generator plugin.

Covers: MopacSetupDialog.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

MOPAC_PATH = PLUGINS_DIR / "Mopac_Input_Generator" / "mopac_input_generator.py"

with mock_chemistry_imports():
    _mopac = load_plugin_for_gui(MOPAC_PATH)


class _FakePos:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class _FakeConf:
    def __init__(self, coords):
        self._coords = coords

    def GetAtomPosition(self, i):
        return _FakePos(*self._coords[i])


class _FakeAtom:
    def __init__(self, symbol, num=6):
        self._symbol = symbol
        self._num = num

    def GetSymbol(self):
        return self._symbol

    def GetAtomicNum(self):
        return self._num


class _FakeMol:
    """Plain duck-typed mol (no real rdkit needed) with a conformer."""

    def __init__(self, atoms, coords):
        self._atoms = atoms
        self._coords = coords

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetAtomWithIdx(self, i):
        return self._atoms[i]

    def GetAtoms(self):
        return list(self._atoms)

    def GetConformer(self):
        return _FakeConf(self._coords)


def _water():
    return _FakeMol(
        [_FakeAtom("O", 8), _FakeAtom("H", 1), _FakeAtom("H", 1)],
        [(0.0, 0.0, 0.117), (0.0, 0.757, -0.469), (0.0, -0.757, -0.469)],
    )


# ===========================================================================
# MopacSetupDialog  (visible plugin: "MOPAC Input Generator")
# ===========================================================================


class TestMopacSetupDialog:
    """MopacSetupDialog with mol=None — exercises the no-molecule UI path."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _mopac.MopacSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "MOPAC Input Generator"

    def test_default_keywords_field(self, dlg):
        assert dlg.keywords_edit.text() == "PM7 PRECISE"

    def test_default_title_field(self, dlg):
        assert dlg.title_edit.text() == "MOPAC Calculation"

    def test_charge_spinbox_range(self, dlg):
        assert dlg.charge_spin.minimum() == -10
        assert dlg.charge_spin.maximum() == 10

    def test_mult_spinbox_minimum_is_one(self, dlg):
        assert dlg.mult_spin.minimum() == 1

    def test_template_combo_populated(self, dlg):
        assert dlg.template_combo.count() > 0

    def test_nosym_checkbox_initially_unchecked(self, dlg):
        assert not dlg.chk_nosym.isChecked()

    def test_preview_area_is_editable(self, dlg):
        assert not dlg.preview_text.isReadOnly()

    def test_preset_combo_initially_empty(self, dlg):
        assert dlg.preset_combo.count() == 0


class TestMopacApplyTemplateReal:
    @pytest.fixture
    def dlg(self, qapp):
        d = _mopac.MopacSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_template_selection_updates_keywords(self, dlg):
        dlg.template_combo.setCurrentText("Single Point (PM7)")
        assert dlg.keywords_edit.text() == "PM7 1SCF"

    def test_placeholder_selection_leaves_keywords(self, dlg):
        dlg.keywords_edit.setText("CUSTOM")
        dlg.template_combo.setCurrentIndex(0)
        assert dlg.keywords_edit.text() == "CUSTOM"

    def test_thermodynamics_template(self, dlg):
        dlg.template_combo.setCurrentText("Thermodynamics (PM7)")
        assert dlg.keywords_edit.text() == "PM7 THERMO"

    def test_transition_state_template(self, dlg):
        dlg.template_combo.setCurrentText("Transition State (PM7 TS)")
        assert dlg.keywords_edit.text() == "PM7 TS PRECISE"


class TestMopacGenerateContentReal:
    @pytest.fixture
    def dlg(self, qapp):
        d = _mopac.MopacSetupDialog(parent=None, mol=_water())
        yield d
        d.destroy()

    def test_geometry_lines_present(self, dlg):
        content = dlg.preview_text.toPlainText()
        lines = content.split("\n")
        assert len(lines) == 6
        assert lines[3].startswith("O")

    def test_charge_and_mult_change_updates_preview(self, dlg):
        dlg.charge_spin.setValue(-1)
        dlg.mult_spin.setValue(2)
        content = dlg.preview_text.toPlainText()
        assert "CHARGE=-1" in content.split("\n")[0]
        assert "DOUBLET" in content.split("\n")[0]

    def test_nosym_checkbox_updates_preview(self, dlg):
        dlg.chk_nosym.setChecked(True)
        assert "NOSYM" in dlg.preview_text.toPlainText().split("\n")[0]

    def test_1scf_keyword_gives_fixed_geometry(self, dlg):
        dlg.keywords_edit.setText("PM7 1SCF")
        atom_line = dlg.preview_text.toPlainText().split("\n")[3]
        assert atom_line.split()[2] == "0"

    def test_calc_initial_charge_mult_sets_spinboxes(self, dlg, monkeypatch):
        monkeypatch.setattr(_mopac.Chem, "GetFormalCharge", lambda m: 0)
        dlg.calc_initial_charge_mult()
        assert dlg.charge_spin.value() == 0
        assert dlg.mult_spin.value() == 1

    def test_calc_initial_charge_mult_odd_electrons(self, dlg, monkeypatch):
        mol = _FakeMol([_FakeAtom("O", 8), _FakeAtom("H", 1)], [(0, 0, 0), (1, 0, 0)])
        dlg.mol = mol
        monkeypatch.setattr(_mopac.Chem, "GetFormalCharge", lambda m: 0)
        dlg.calc_initial_charge_mult()
        assert dlg.mult_spin.value() == 2

    def test_calc_initial_charge_mult_exception_is_silenced(self, dlg, monkeypatch):
        def _raise(_m):
            raise ValueError("boom")

        monkeypatch.setattr(_mopac.Chem, "GetFormalCharge", _raise)
        dlg.calc_initial_charge_mult()  # should not raise


class TestMopacSaveFile:
    @pytest.fixture
    def dlg(self, qapp):
        d = _mopac.MopacSetupDialog(parent=None, mol=_water(), filename="input.xyz")
        yield d
        d.destroy()

    def test_save_writes_file(self, dlg, monkeypatch, tmp_path):
        out = tmp_path / "job.mop"
        monkeypatch.setattr(
            _mopac.QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), "")
        )
        monkeypatch.setattr(_mopac.QMessageBox, "information", MagicMock())
        dlg.save_file()
        assert out.exists()
        assert out.read_text(encoding="utf-8") == dlg.preview_text.toPlainText()

    def test_save_cancelled_writes_nothing(self, dlg, monkeypatch, tmp_path):
        monkeypatch.setattr(
            _mopac.QFileDialog, "getSaveFileName", lambda *a, **k: ("", "")
        )
        dlg.save_file()  # no exception, no file created

    def test_save_write_error_shows_critical(self, dlg, monkeypatch, tmp_path):
        # Missing parent dir fails to open on every OS (a bare "Z:\..." string is
        # a *valid* filename on Linux, so the write would succeed there).
        bad_path = str(tmp_path / "missing_dir" / "job.mop")
        monkeypatch.setattr(
            _mopac.QFileDialog, "getSaveFileName", lambda *a, **k: (bad_path, "")
        )
        crit = MagicMock()
        monkeypatch.setattr(_mopac.QMessageBox, "critical", crit)
        monkeypatch.setattr(_mopac.QMessageBox, "information", MagicMock())
        dlg.save_file()
        crit.assert_called_once()


class TestMopacPresetsReal:
    @pytest.fixture
    def dlg(self, qapp, tmp_path, monkeypatch):
        monkeypatch.setattr(_mopac, "SETTINGS_FILE", str(tmp_path / "mopac.json"))
        d = _mopac.MopacSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_save_and_reload_preset(self, dlg, monkeypatch):
        monkeypatch.setattr(
            _mopac.QInputDialog, "getText", lambda *a, **k: ("MyPreset", True)
        )
        dlg.keywords_edit.setText("PM6 PRECISE")
        dlg.chk_nosym.setChecked(True)
        dlg.save_preset_dialog()
        assert "MyPreset" in dlg.presets_data
        assert dlg.preset_combo.currentText() == "MyPreset"

        dlg.load_presets_from_file()
        assert "MyPreset" in dlg.presets_data
        assert dlg.keywords_edit.text() == "PM6 PRECISE"

    def test_save_preset_cancelled_does_nothing(self, dlg, monkeypatch):
        monkeypatch.setattr(
            _mopac.QInputDialog, "getText", lambda *a, **k: ("", False)
        )
        dlg.save_preset_dialog()
        assert dlg.presets_data == {}

    def test_delete_preset_confirmed(self, dlg, monkeypatch):
        monkeypatch.setattr(
            _mopac.QInputDialog, "getText", lambda *a, **k: ("ToDelete", True)
        )
        dlg.save_preset_dialog()
        assert "ToDelete" in dlg.presets_data

        monkeypatch.setattr(
            _mopac.QMessageBox,
            "question",
            lambda *a, **k: _mopac.QMessageBox.StandardButton.Yes,
        )
        dlg.delete_preset()
        assert "ToDelete" not in dlg.presets_data

    def test_delete_preset_no_selection_is_noop(self, dlg):
        dlg.preset_combo.clear()
        dlg.delete_preset()  # no exception

    def test_load_corrupt_presets_file_gives_empty(self, dlg, monkeypatch, tmp_path):
        bad = tmp_path / "bad.json"
        bad.write_text("{not json", encoding="utf-8")
        monkeypatch.setattr(_mopac, "SETTINGS_FILE", str(bad))
        dlg.load_presets_from_file()
        assert dlg.presets_data == {}
