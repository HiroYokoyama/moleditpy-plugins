"""
Headless GUI tests for the NWChem Input Generator plugin.

Covers: NwchemSetupDialog.

Follows the same pattern as MOPAC/GAMESS/Psi4:
- parent=None, mol=None is safe (calc_initial_charge_mult guards with `if not self.mol`)
- setup_ui() is pure Qt widget construction
- load_presets_from_file() reads a settings file only if it exists

Run locally:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/test_gui_plugin_nwchem_input_generator.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
NWCHEM_PATH = PLUGINS_DIR / "NWChem_Input_Generator" / "nwchem_input_generator.py"

with mock_chemistry_imports():
    _nwchem = load_plugin_for_gui(NWCHEM_PATH)


# ===========================================================================
# NwchemSetupDialog  (visible plugin: "NWChem Input Generator")
# ===========================================================================


class TestNwchemSetupDialog:
    """NwchemSetupDialog with mol=None — exercises the no-molecule UI path."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _nwchem.NwchemSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "NWChem Input Generator"

    def test_default_title_field(self, dlg):
        assert dlg.title_edit.text() == "NWChem Job"

    def test_module_combo_default(self, dlg):
        assert dlg.module_combo.currentText() == "dft"

    def test_functional_combo_default(self, dlg):
        assert dlg.functional_combo.currentText() == "b3lyp"

    def test_task_combo_default(self, dlg):
        assert dlg.task_combo.currentText() == "optimize"

    def test_basis_combo_default(self, dlg):
        assert dlg.basis_combo.currentText() == "6-31G*"

    def test_charge_spinbox_range(self, dlg):
        assert dlg.charge_spin.minimum() == -10
        assert dlg.charge_spin.maximum() == 10

    def test_mult_spinbox_minimum_is_one(self, dlg):
        assert dlg.mult_spin.minimum() == 1

    def test_preset_combo_initially_empty(self, dlg):
        assert dlg.preset_combo.count() == 0

    def test_preview_area_is_editable(self, dlg):
        assert not dlg.preview_text.isReadOnly()

    def test_module_combo_items_include_scf(self, dlg):
        items = [dlg.module_combo.itemText(i) for i in range(dlg.module_combo.count())]
        assert "scf" in items

    def test_basis_combo_includes_cc_pvdz(self, dlg):
        items = [dlg.basis_combo.itemText(i) for i in range(dlg.basis_combo.count())]
        assert "cc-pvdz" in items

    def test_task_combo_includes_freq(self, dlg):
        items = [dlg.task_combo.itemText(i) for i in range(dlg.task_combo.count())]
        assert "freq" in items


# ===========================================================================
# Fake RDKit-like mol objects (real Python, not MagicMock, so real widget
# calls like setValue()/int math work against them)
# ===========================================================================


class _FakePos:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class _FakeAtom:
    def __init__(self, symbol, num):
        self._symbol = symbol
        self._num = num

    def GetSymbol(self):
        return self._symbol

    def GetAtomicNum(self):
        return self._num


class _FakeConformer:
    def __init__(self, coords):
        self._coords = coords

    def GetAtomPosition(self, i):
        return _FakePos(*self._coords[i])


class _FakeMol:
    def __init__(self, atoms, coords):
        self._atoms = atoms
        self._coords = coords

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetAtomWithIdx(self, i):
        return self._atoms[i]

    def GetConformer(self):
        return _FakeConformer(self._coords)

    def GetAtoms(self):
        return self._atoms


def _water():
    return _FakeMol(
        [_FakeAtom("O", 8), _FakeAtom("H", 1), _FakeAtom("H", 1)],
        [(0.0, 0.0, 0.117), (0.0, 0.757, -0.469), (0.0, -0.757, -0.469)],
    )


def _iodine_pair():
    return _FakeMol(
        [_FakeAtom("I", 53), _FakeAtom("H", 1)],
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
    )


def _no_block_msgbox(monkeypatch):
    """QMessageBox.* pop up modal dialogs that block exec() under offscreen
    Qt with no user interaction; replace with recording stubs."""
    calls = {"warning": [], "information": [], "critical": [], "question": []}
    for kind in ("warning", "information", "critical"):
        monkeypatch.setattr(
            _nwchem.QMessageBox, kind,
            staticmethod(lambda *a, _k=kind, **kw: calls[_k].append(a)),
        )
    return calls


# ===========================================================================
# calc_initial_charge_mult (real mol, real GetFormalCharge)
# ===========================================================================


class TestNwchemChargeMultReal:
    def test_neutral_singlet_water(self, qapp, monkeypatch):
        monkeypatch.setattr(_nwchem.Chem, "GetFormalCharge", lambda mol: 0)
        d = _nwchem.NwchemSetupDialog(parent=None, mol=_water())
        try:
            assert d.charge_spin.value() == 0
            assert d.mult_spin.value() == 1
        finally:
            d.destroy()

    def test_charged_radical_doublet(self, qapp, monkeypatch):
        # CH3 radical: 6 + 1 + 1 + 1 = 9 electrons, odd -> doublet
        mol = _FakeMol(
            [_FakeAtom("C", 6), _FakeAtom("H", 1), _FakeAtom("H", 1), _FakeAtom("H", 1)],
            [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0)],
        )
        monkeypatch.setattr(_nwchem.Chem, "GetFormalCharge", lambda mol: 0)
        d = _nwchem.NwchemSetupDialog(parent=None, mol=mol)
        try:
            assert d.mult_spin.value() == 2
        finally:
            d.destroy()

    def test_exception_is_silenced(self, qapp, monkeypatch):
        def boom(mol):
            raise RuntimeError("bad mol")

        monkeypatch.setattr(_nwchem.Chem, "GetFormalCharge", boom)
        d = _nwchem.NwchemSetupDialog(parent=None, mol=_water())
        try:
            assert d is not None
        finally:
            d.destroy()


# ===========================================================================
# check_heavy_atoms (real mol)
# ===========================================================================


class TestNwchemHeavyAtomsReal:
    def test_heavy_atom_sets_warning(self, qapp, monkeypatch):
        monkeypatch.setattr(_nwchem.Chem, "GetFormalCharge", lambda mol: 0)
        d = _nwchem.NwchemSetupDialog(parent=None, mol=_iodine_pair())
        try:
            assert "ECP" in d.ecp_warning.text()
        finally:
            d.destroy()

    def test_light_atoms_no_warning(self, qapp, monkeypatch):
        monkeypatch.setattr(_nwchem.Chem, "GetFormalCharge", lambda mol: 0)
        d = _nwchem.NwchemSetupDialog(parent=None, mol=_water())
        try:
            assert d.ecp_warning.text() == ""
        finally:
            d.destroy()


# ===========================================================================
# generate_content via real widgets — atom lines, module/charge combinations
# ===========================================================================


class TestNwchemPreviewReal:
    @pytest.fixture
    def dlg(self, qapp, monkeypatch):
        monkeypatch.setattr(_nwchem.Chem, "GetFormalCharge", lambda mol: 0)
        d = _nwchem.NwchemSetupDialog(parent=None, mol=_water())
        yield d
        d.destroy()

    def test_atom_coordinates_in_preview(self, dlg):
        text = dlg.preview_text.toPlainText()
        assert "O" in text and "0.117000" in text

    def test_nonzero_charge_shown_in_preview(self, dlg):
        dlg.charge_spin.setValue(-1)
        text = dlg.preview_text.toPlainText()
        assert "charge -1" in text.splitlines()

    def test_scf_module_shows_spin_keyword(self, dlg):
        dlg.module_combo.setCurrentText("scf")
        dlg.mult_spin.setValue(2)
        text = dlg.preview_text.toPlainText()
        assert "  uhf" in text.splitlines()
        assert "  doublet" in text.splitlines()

    def test_scf_module_disables_functional_combo(self, dlg):
        dlg.module_combo.setCurrentText("scf")
        assert not dlg.functional_combo.isEnabled()

    def test_dft_module_reenables_functional_combo(self, dlg):
        dlg.module_combo.setCurrentText("scf")
        dlg.module_combo.setCurrentText("dft")
        assert dlg.functional_combo.isEnabled()


# ===========================================================================
# save_file — real widget, mocked QFileDialog/QMessageBox
# ===========================================================================


class TestNwchemSaveFileReal:
    @pytest.fixture
    def dlg(self, qapp, monkeypatch):
        monkeypatch.setattr(_nwchem.Chem, "GetFormalCharge", lambda mol: 0)
        d = _nwchem.NwchemSetupDialog(parent=None, mol=_water())
        yield d
        d.destroy()

    def test_save_replaces_start_line(self, dlg, monkeypatch, tmp_path):
        target = tmp_path / "water_job.nw"
        calls = _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(
            _nwchem.QFileDialog, "getSaveFileName",
            staticmethod(lambda *a, **k: (str(target), "")),
        )
        dlg.save_file()
        content = target.read_text(encoding="utf-8")
        assert "start water_job" in content
        assert len(calls["information"]) == 1

    def test_save_uses_filename_stem_as_default(self, qapp, monkeypatch):
        monkeypatch.setattr(_nwchem.Chem, "GetFormalCharge", lambda mol: 0)
        d = _nwchem.NwchemSetupDialog(parent=None, mol=_water(), filename="/tmp/foo/bar.xyz")
        try:
            captured = {}

            def fake_get_save(*a, **k):
                captured["default"] = a[2]
                return ("", "")

            monkeypatch.setattr(
                _nwchem.QFileDialog, "getSaveFileName", staticmethod(fake_get_save),
            )
            d.save_file()
            assert captured["default"] == "bar.nw"
        finally:
            d.destroy()

    def test_save_cancelled_writes_nothing(self, dlg, monkeypatch, tmp_path):
        monkeypatch.setattr(
            _nwchem.QFileDialog, "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        dlg.save_file()
        assert list(tmp_path.glob("*.nw")) == []

    def test_save_write_failure_shows_critical(self, dlg, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(
            _nwchem.QFileDialog, "getSaveFileName",
            staticmethod(lambda *a, **k: ("/nonexistent_dir_xyz/job.nw", "")),
        )
        dlg.save_file()
        assert len(calls["critical"]) == 1


# ===========================================================================
# Preset management — real widgets, real settings file on disk
# ===========================================================================


class TestNwchemPresetsReal:
    @pytest.fixture
    def settings_file(self, monkeypatch, tmp_path):
        path = tmp_path / "presets.json"
        monkeypatch.setattr(_nwchem, "SETTINGS_FILE", str(path))
        return path

    @pytest.fixture
    def dlg(self, qapp, settings_file):
        d = _nwchem.NwchemSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_save_preset_persists_and_selects(self, dlg, settings_file, monkeypatch):
        import json

        monkeypatch.setattr(
            _nwchem.QInputDialog, "getText",
            staticmethod(lambda *a, **k: ("My Preset", True)),
        )
        dlg.module_combo.setCurrentText("mp2")
        dlg.title_edit.setText("Custom Title")
        dlg.save_preset_dialog()

        assert dlg.preset_combo.currentText() == "My Preset"
        on_disk = json.loads(settings_file.read_text(encoding="utf-8"))
        assert on_disk["My Preset"]["module"] == "mp2"
        assert on_disk["My Preset"]["title"] == "Custom Title"

    def test_save_preset_cancelled_is_noop(self, dlg, settings_file, monkeypatch):
        monkeypatch.setattr(
            _nwchem.QInputDialog, "getText",
            staticmethod(lambda *a, **k: ("", False)),
        )
        dlg.save_preset_dialog()
        assert dlg.presets_data == {}
        assert not settings_file.exists()

    def test_new_dialog_loads_saved_preset_and_applies(self, qapp, settings_file):
        import json

        settings_file.write_text(
            json.dumps(
                {
                    "Fast": {
                        "module": "scf",
                        "functional": "pbe0",
                        "task": "energy",
                        "basis": "cc-pvtz",
                        "title": "Fast Job",
                        "start": "fast_start",
                    }
                }
            ),
            encoding="utf-8",
        )
        d = _nwchem.NwchemSetupDialog(parent=None, mol=None)
        try:
            assert d.preset_combo.currentText() == "Fast"
            assert d.module_combo.currentText() == "scf"
            assert d.functional_combo.currentText() == "pbe0"
            assert d.task_combo.currentText() == "energy"
            assert d.basis_combo.currentText() == "cc-pvtz"
            assert d.title_edit.text() == "Fast Job"
            assert d.start_edit.text() == "fast_start"
            assert not d.functional_combo.isEnabled()
        finally:
            d.destroy()

    def test_delete_preset_confirmed_removes_entry(self, dlg, settings_file, monkeypatch):
        import json

        dlg.presets_data = {"Doomed": {"module": "dft"}}
        dlg.update_preset_combo()
        assert dlg.preset_combo.currentText() == "Doomed"
        monkeypatch.setattr(
            _nwchem.QMessageBox, "question",
            staticmethod(lambda *a, **k: _nwchem.QMessageBox.StandardButton.Yes),
        )
        dlg.delete_preset()
        assert dlg.presets_data == {}
        assert dlg.preset_combo.count() == 0
        assert json.loads(settings_file.read_text(encoding="utf-8")) == {}

    def test_delete_preset_declined_keeps_entry(self, dlg, monkeypatch):
        dlg.presets_data = {"Kept": {"module": "dft"}}
        dlg.update_preset_combo()
        monkeypatch.setattr(
            _nwchem.QMessageBox, "question",
            staticmethod(lambda *a, **k: _nwchem.QMessageBox.StandardButton.No),
        )
        dlg.delete_preset()
        assert "Kept" in dlg.presets_data

    def test_delete_preset_no_selection_is_noop(self, dlg):
        dlg.delete_preset()  # empty combo -> name is "" -> early return, no raise

    def test_corrupt_settings_file_logs_and_starts_empty(self, qapp, monkeypatch, tmp_path):
        path = tmp_path / "bad.json"
        path.write_text("{not valid json", encoding="utf-8")
        monkeypatch.setattr(_nwchem, "SETTINGS_FILE", str(path))
        d = _nwchem.NwchemSetupDialog(parent=None, mol=None)
        try:
            assert d.presets_data == {}
        finally:
            d.destroy()


# ===========================================================================
# run(mw) — module entry point
# ===========================================================================


def _fake_main_window(current_mol=None, current_file_path=None):
    from PyQt6.QtWidgets import QWidget

    w = QWidget()
    w.current_mol = current_mol
    w.current_file_path = current_file_path
    return w


class TestNwchemRun:
    def test_no_molecule_shows_warning(self, qapp, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        mw = _fake_main_window(current_mol=None)
        _nwchem.run(mw)
        assert len(calls["warning"]) == 1

    def test_with_molecule_opens_dialog(self, qapp, monkeypatch):
        monkeypatch.setattr(_nwchem.Chem, "GetFormalCharge", lambda mol: 0)
        exec_calls = []
        monkeypatch.setattr(
            _nwchem.NwchemSetupDialog, "exec",
            lambda self: exec_calls.append(self) or 0,
        )
        mw = _fake_main_window(current_mol=_water(), current_file_path="water.xyz")
        _nwchem.run(mw)
        assert len(exec_calls) == 1
