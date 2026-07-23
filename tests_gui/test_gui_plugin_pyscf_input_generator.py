"""
Headless GUI tests for the PySCF Input Generator plugin.

Covers: PyscfSetupDialog.

Follows the same pattern as MOPAC/GAMESS/Psi4:
- parent=None, mol=None is safe (calc_initial_charge_mult guards with `if not self.mol`)
- setup_ui() is pure Qt widget construction
- load_presets_from_file() reads a settings file only if it exists

Run locally:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/test_gui_plugin_pyscf_input_generator.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
PYSCF_PATH = PLUGINS_DIR / "PySCF_Input_Generator" / "pyscf_input_generator.py"

with mock_chemistry_imports():
    _pyscf = load_plugin_for_gui(PYSCF_PATH)


# ===========================================================================
# PyscfSetupDialog  (visible plugin: "PySCF Input Generator")
# ===========================================================================


class TestPyscfSetupDialog:
    """PyscfSetupDialog with mol=None — exercises the no-molecule UI path."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _pyscf.PyscfSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "PySCF Input Generator"

    def test_category_combo_default(self, dlg):
        assert dlg.category_combo.currentText() == "Hartree-Fock"

    def test_basis_combo_default(self, dlg):
        assert dlg.basis_combo.currentText() == "def2-svp"

    def test_charge_spinbox_range(self, dlg):
        assert dlg.charge_spin.minimum() == -10
        assert dlg.charge_spin.maximum() == 10

    def test_mult_spinbox_minimum_is_one(self, dlg):
        assert dlg.mult_spin.minimum() == 1

    def test_symmetry_checked_by_default(self, dlg):
        assert dlg.symmetry_check.isChecked()

    def test_functional_combo_initially_disabled(self, dlg):
        assert not dlg.functional_combo.isEnabled()

    def test_post_hf_combo_initially_disabled(self, dlg):
        assert not dlg.post_hf_combo.isEnabled()

    def test_post_hf_ref_combo_initially_disabled(self, dlg):
        assert not dlg.post_hf_ref_combo.isEnabled()

    def test_preset_combo_initially_empty(self, dlg):
        assert dlg.preset_combo.count() == 0

    def test_preview_area_is_editable(self, dlg):
        assert not dlg.preview_text.isReadOnly()

    def test_basis_combo_includes_sto_3g(self, dlg):
        items = [dlg.basis_combo.itemText(i) for i in range(dlg.basis_combo.count())]
        assert "sto-3g" in items

    def test_category_combo_includes_dft(self, dlg):
        items = [dlg.category_combo.itemText(i) for i in range(dlg.category_combo.count())]
        assert "DFT" in items

    def test_save_button_label(self, dlg):
        assert dlg.btn_save.text() == "Save Python Script..."


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


def _no_block_msgbox(monkeypatch):
    """QMessageBox.* pop up modal dialogs that block exec() under offscreen
    Qt with no user interaction; replace with recording stubs."""
    calls = {"warning": [], "information": [], "critical": [], "question": []}
    for kind in ("warning", "information", "critical"):
        monkeypatch.setattr(
            _pyscf.QMessageBox, kind,
            staticmethod(lambda *a, _k=kind, **kw: calls[_k].append(a)),
        )
    return calls


# ===========================================================================
# calc_initial_charge_mult (real mol, real GetFormalCharge)
# ===========================================================================


class TestPyscfChargeMultReal:
    def test_neutral_singlet_water(self, qapp, monkeypatch):
        monkeypatch.setattr(_pyscf.Chem, "GetFormalCharge", lambda mol: 0)
        d = _pyscf.PyscfSetupDialog(parent=None, mol=_water())
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
        monkeypatch.setattr(_pyscf.Chem, "GetFormalCharge", lambda mol: 0)
        d = _pyscf.PyscfSetupDialog(parent=None, mol=mol)
        try:
            assert d.mult_spin.value() == 2
        finally:
            d.destroy()

    def test_exception_is_silenced(self, qapp, monkeypatch):
        def boom(mol):
            raise RuntimeError("bad mol")

        monkeypatch.setattr(_pyscf.Chem, "GetFormalCharge", boom)
        # Must not raise despite calc_initial_charge_mult's internal failure.
        d = _pyscf.PyscfSetupDialog(parent=None, mol=_water())
        try:
            assert d is not None
        finally:
            d.destroy()


# ===========================================================================
# generate_content via real widgets (update_preview) — atom lines, method
# selection across all categories/spins, post-HF variants
# ===========================================================================


class TestPyscfPreviewReal:
    @pytest.fixture
    def dlg(self, qapp, monkeypatch):
        monkeypatch.setattr(_pyscf.Chem, "GetFormalCharge", lambda mol: 0)
        d = _pyscf.PyscfSetupDialog(parent=None, mol=_water())
        yield d
        d.destroy()

    def test_atom_coordinates_in_preview(self, dlg):
        text = dlg.preview_text.toPlainText()
        assert "O 0.00000000 0.00000000 0.11700000" in text

    def test_category_dft_enables_functional_and_updates_preview(self, dlg):
        dlg.category_combo.setCurrentText("DFT")
        assert dlg.functional_combo.isEnabled()
        text = dlg.preview_text.toPlainText()
        assert "mf = dft.RKS(mol)" in text

    def test_category_post_hf_enables_combos(self, dlg):
        dlg.category_combo.setCurrentText("Post-HF (MP2, CCSD)")
        assert dlg.post_hf_combo.isEnabled()
        assert dlg.post_hf_ref_combo.isEnabled()

    def test_open_shell_hf_uses_rohf(self, dlg):
        dlg.mult_spin.setValue(2)
        text = dlg.preview_text.toPlainText()
        assert "mf = scf.ROHF(mol)" in text

    def test_open_shell_dft_uses_roks(self, dlg):
        dlg.category_combo.setCurrentText("DFT")
        dlg.mult_spin.setValue(2)
        text = dlg.preview_text.toPlainText()
        assert "mf = dft.ROKS(mol)" in text

    def test_functional_change_reflected(self, dlg):
        dlg.category_combo.setCurrentText("DFT")
        dlg.functional_combo.setCurrentText("wb97x-d")
        text = dlg.preview_text.toPlainText()
        assert "mf.xc = 'wb97x-d'" in text

    @pytest.mark.parametrize("ref", ["RHF", "UHF", "ROHF"])
    def test_post_hf_reference_reflected(self, dlg, ref):
        dlg.category_combo.setCurrentText("Post-HF (MP2, CCSD)")
        dlg.post_hf_ref_combo.setCurrentText(ref)
        text = dlg.preview_text.toPlainText()
        assert f"mf = scf.{ref}(mol)" in text

    @pytest.mark.parametrize(
        "method,expected",
        [
            ("MP2", "pt = mp.MP2(mf)"),
            ("CCSD", "mycc = cc.CCSD(mf)"),
            ("CCSD(T)", "e_t = mycc.ccsd_t()"),
        ],
    )
    def test_post_hf_method_reflected(self, dlg, method, expected):
        dlg.category_combo.setCurrentText("Post-HF (MP2, CCSD)")
        dlg.post_hf_combo.setCurrentText(method)
        text = dlg.preview_text.toPlainText()
        assert expected in text

    def test_symmetry_toggle_reflected(self, dlg):
        dlg.symmetry_check.setChecked(False)
        assert "symmetry=False" in dlg.preview_text.toPlainText()

    def test_basis_change_reflected(self, dlg):
        dlg.basis_combo.setCurrentText("cc-pvtz")
        assert "basis='cc-pvtz'," in dlg.preview_text.toPlainText()


# ===========================================================================
# save_file — real widget, mocked QFileDialog/QMessageBox
# ===========================================================================


class TestPyscfSaveFileReal:
    @pytest.fixture
    def dlg(self, qapp, monkeypatch):
        monkeypatch.setattr(_pyscf.Chem, "GetFormalCharge", lambda mol: 0)
        d = _pyscf.PyscfSetupDialog(parent=None, mol=_water())
        yield d
        d.destroy()

    def test_save_writes_chk_substituted_content(self, dlg, monkeypatch, tmp_path):
        target = tmp_path / "water_opt.py"
        calls = _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(
            _pyscf.QFileDialog, "getSaveFileName",
            staticmethod(lambda *a, **k: (str(target), "")),
        )
        dlg.save_file()
        content = target.read_text(encoding="utf-8")
        assert "mf.chkfile = 'water_opt.chk'" in content
        assert len(calls["information"]) == 1

    def test_save_cancelled_writes_nothing(self, dlg, monkeypatch, tmp_path):
        monkeypatch.setattr(
            _pyscf.QFileDialog, "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        dlg.save_file()
        assert list(tmp_path.glob("*.py")) == []

    def test_save_error_shows_critical(self, dlg, monkeypatch, tmp_path):
        # A directory path can't be opened for writing -> IOError -> critical box.
        calls = _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(
            _pyscf.QFileDialog, "getSaveFileName",
            staticmethod(lambda *a, **k: (str(tmp_path), "")),
        )
        dlg.save_file()
        assert len(calls["critical"]) == 1

    def test_default_save_name_derived_from_source_filename(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_pyscf.Chem, "GetFormalCharge", lambda mol: 0)
        d = _pyscf.PyscfSetupDialog(parent=None, mol=_water(), filename="mystructure.xyz")
        try:
            captured = {}

            def _fake_get_save_file_name(*a, **k):
                captured["hint"] = a[2]
                return (str(tmp_path / "out.py"), "")

            monkeypatch.setattr(
                _pyscf.QFileDialog, "getSaveFileName",
                staticmethod(_fake_get_save_file_name),
            )
            _no_block_msgbox(monkeypatch)
            d.save_file()
            assert captured["hint"] == "mystructure_pyscf.py"
        finally:
            d.destroy()

    def test_manually_edited_content_missing_chkfile_inserted_before_kernel(self, dlg, monkeypatch, tmp_path):
        target = tmp_path / "abc.py"
        _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(
            _pyscf.QFileDialog, "getSaveFileName",
            staticmethod(lambda *a, **k: (str(target), "")),
        )
        dlg.preview_text.setPlainText("mf = scf.RHF(mol)\nmf.kernel()\n")
        dlg.save_file()
        content = target.read_text(encoding="utf-8")
        assert "mf.chkfile = 'abc.chk'\nmf.kernel()" in content

    def test_manually_edited_content_without_kernel_line_appends_chk(self, dlg, monkeypatch, tmp_path):
        target = tmp_path / "manual.py"
        _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(
            _pyscf.QFileDialog, "getSaveFileName",
            staticmethod(lambda *a, **k: (str(target), "")),
        )
        dlg.preview_text.setPlainText("print('no kernel or chkfile here')")
        dlg.save_file()
        content = target.read_text(encoding="utf-8")
        assert "mf.chkfile = 'manual.chk'" in content


# ===========================================================================
# Presets: load/save/apply/delete via the real widget
# ===========================================================================


class TestPyscfPresetsReal:
    @pytest.fixture
    def dlg(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_pyscf, "SETTINGS_FILE", str(tmp_path / "presets.json"))
        d = _pyscf.PyscfSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_save_preset_persists_and_selects(self, dlg, monkeypatch, tmp_path):
        import json

        monkeypatch.setattr(
            _pyscf.QInputDialog, "getText",
            staticmethod(lambda *a, **k: ("Fast HF", True)),
        )
        dlg.category_combo.setCurrentText("DFT")
        dlg.functional_combo.setCurrentText("pbe")
        dlg.save_preset_dialog()

        assert dlg.preset_combo.currentText() == "Fast HF"
        on_disk = json.loads((tmp_path / "presets.json").read_text(encoding="utf-8"))
        assert on_disk["Fast HF"]["category"] == "DFT"
        assert on_disk["Fast HF"]["functional"] == "pbe"

    def test_save_preset_cancelled_is_noop(self, dlg, monkeypatch):
        monkeypatch.setattr(
            _pyscf.QInputDialog, "getText",
            staticmethod(lambda *a, **k: ("", False)),
        )
        dlg.save_preset_dialog()
        assert dlg.presets_data == {}

    def test_new_dialog_loads_saved_preset_and_applies(self, qapp, monkeypatch, tmp_path):
        import json

        path = tmp_path / "presets.json"
        path.write_text(
            json.dumps(
                {
                    "Saved": {
                        "category": "Post-HF (MP2, CCSD)",
                        "functional": "b3lyp",
                        "post_hf": "CCSD",
                        "post_hf_ref": "UHF",
                        "basis": "cc-pvdz",
                        "symmetry": False,
                    }
                }
            ),
            encoding="utf-8",
        )
        monkeypatch.setattr(_pyscf, "SETTINGS_FILE", str(path))
        d = _pyscf.PyscfSetupDialog(parent=None, mol=None)
        try:
            assert d.preset_combo.currentText() == "Saved"
            assert d.category_combo.currentText() == "Post-HF (MP2, CCSD)"
            assert d.post_hf_combo.currentText() == "CCSD"
            assert d.post_hf_ref_combo.currentText() == "UHF"
            assert d.basis_combo.currentText() == "cc-pvdz"
            assert not d.symmetry_check.isChecked()
        finally:
            d.destroy()

    def test_corrupt_presets_file_is_silenced(self, qapp, monkeypatch, tmp_path):
        path = tmp_path / "presets.json"
        path.write_text("{not valid json", encoding="utf-8")
        monkeypatch.setattr(_pyscf, "SETTINGS_FILE", str(path))
        d = _pyscf.PyscfSetupDialog(parent=None, mol=None)
        try:
            assert d.presets_data == {}
        finally:
            d.destroy()

    def test_update_preset_combo_reselects_current(self, dlg):
        # Second call with an already-current, still-valid preset name exercises
        # the "if current in presets_data: setCurrentText(current)" branch.
        dlg.presets_data = {"Foo": {"category": "Hartree-Fock"}}
        dlg.update_preset_combo()
        assert dlg.preset_combo.currentText() == "Foo"
        dlg.update_preset_combo()
        assert dlg.preset_combo.currentText() == "Foo"

    def test_delete_preset_confirmed_removes_entry(self, dlg, monkeypatch, tmp_path):
        import json

        dlg.presets_data = {"Doomed": {"category": "Hartree-Fock"}}
        dlg.update_preset_combo()
        monkeypatch.setattr(
            _pyscf.QMessageBox, "question",
            staticmethod(lambda *a, **k: _pyscf.QMessageBox.StandardButton.Yes),
        )
        dlg.delete_preset()
        assert dlg.presets_data == {}
        assert json.loads((tmp_path / "presets.json").read_text(encoding="utf-8")) == {}

    def test_delete_preset_declined_keeps_entry(self, dlg, monkeypatch):
        dlg.presets_data = {"Kept": {"category": "Hartree-Fock"}}
        dlg.update_preset_combo()
        monkeypatch.setattr(
            _pyscf.QMessageBox, "question",
            staticmethod(lambda *a, **k: _pyscf.QMessageBox.StandardButton.No),
        )
        dlg.delete_preset()
        assert "Kept" in dlg.presets_data

    def test_delete_preset_no_selection_is_noop(self, dlg):
        assert dlg.preset_combo.currentText() == ""
        dlg.delete_preset()  # must not raise

    def test_save_presets_error_shows_warning(self, dlg, monkeypatch, tmp_path):
        # Parent directory doesn't exist -> open() raises -> warning box.
        bad_path = tmp_path / "missing_dir" / "presets.json"
        monkeypatch.setattr(_pyscf, "SETTINGS_FILE", str(bad_path))
        calls = _no_block_msgbox(monkeypatch)
        dlg.save_presets_to_file()
        assert len(calls["warning"]) == 1


# ===========================================================================
# run(mw) — module entry point
# ===========================================================================


def _fake_main_window(qapp, current_mol=None, current_file_path=None):
    # run(mw) passes mw straight through as the QDialog parent, so it must be
    # a real QWidget (not an arbitrary Python object) under real Qt.
    from PyQt6.QtWidgets import QWidget

    w = QWidget()
    w.current_mol = current_mol
    w.current_file_path = current_file_path
    return w


class TestPyscfRun:
    def test_no_molecule_shows_warning(self, qapp, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        mw = _fake_main_window(qapp, current_mol=None)
        _pyscf.run(mw)
        assert len(calls["warning"]) == 1

    def test_with_molecule_opens_dialog(self, qapp, monkeypatch):
        monkeypatch.setattr(_pyscf.Chem, "GetFormalCharge", lambda mol: 0)
        exec_calls = []
        monkeypatch.setattr(
            _pyscf.PyscfSetupDialog, "exec",
            lambda self: exec_calls.append(self) or 0,
        )
        mw = _fake_main_window(qapp, current_mol=_water(), current_file_path="water.xyz")
        _pyscf.run(mw)
        assert len(exec_calls) == 1
        assert exec_calls[0].filename == "water.xyz"
