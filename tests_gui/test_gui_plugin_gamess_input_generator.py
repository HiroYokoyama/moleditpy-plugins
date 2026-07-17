"""
Headless GUI tests for the GAMESS Input Generator plugin.

Covers: GamessSetupDialog.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

GAMESS_PATH = PLUGINS_DIR / "Gamess_Input_Generator" / "gamess_input_generator.py"

with mock_chemistry_imports():
    _gamess = load_plugin_for_gui(GAMESS_PATH)


# ===========================================================================
# GamessSetupDialog  (visible plugin: "GAMESS Input Generator")
# ===========================================================================


class TestGamessSetupDialog:
    """GamessSetupDialog with mol=None."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _gamess.GamessSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "GAMESS Input Generator"

    def test_run_type_default(self, dlg):
        assert dlg.run_type.currentText() == "OPTIMIZE"

    def test_scf_type_default(self, dlg):
        assert dlg.scf_type.currentText() == "RHF"

    def test_nosym_checked_by_default(self, dlg):
        assert dlg.chk_nosym.isChecked()

    def test_basis_ngauss_default(self, dlg):
        assert dlg.basis_ngauss.value() == 6

    def test_mem_spin_default(self, dlg):
        assert dlg.mem_spin.value() == 100


class TestGamessGenerateContent:
    """generate_content() / update_preview() wiring with mol=None."""

    @pytest.fixture
    def dlg(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_gamess, "SETTINGS_FILE", str(tmp_path / "presets.json"))
        d = _gamess.GamessSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_default_contrl_line(self, dlg):
        content = dlg.generate_content()
        assert "SCFTYP=RHF" in content
        assert "RUNTYP=OPTIMIZE" in content
        assert "NOSYM=1" in content
        assert "DFTTYP" not in content  # None (Hartree-Fock)

    def test_nosym_unchecked_removes_keyword(self, dlg):
        dlg.chk_nosym.setChecked(False)
        assert "NOSYM" not in dlg.generate_content()

    def test_dft_functional_appears_in_contrl(self, dlg):
        dlg.dft_type.setCurrentText("B3LYP")
        assert "DFTTYP=B3LYP" in dlg.generate_content()

    def test_charge_and_mult_propagate(self, dlg):
        dlg.charge_spin.setValue(-1)
        dlg.mult_spin.setValue(2)
        content = dlg.generate_content()
        assert "ICHARG=-1" in content
        assert "MULT=2" in content

    def test_npfunc_zero_omitted_nonzero_included(self, dlg):
        assert "NPFUNC" not in dlg.generate_content()
        dlg.basis_npfunc.setValue(2)
        assert "NPFUNC=2" in dlg.generate_content()

    def test_widget_change_refreshes_preview(self, dlg):
        dlg.run_type.setCurrentText("ENERGY")
        assert "RUNTYP=ENERGY" in dlg.preview_text.toPlainText()

    def test_mol_none_data_block_has_no_atoms(self, dlg):
        content = dlg.generate_content()
        data_block = content.split(" $DATA")[1]
        assert "C1" in data_block
        assert len(data_block.strip().splitlines()) == 3  # title, C1, $END


class TestGamessBasisOptions:
    @pytest.fixture
    def dlg(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_gamess, "SETTINGS_FILE", str(tmp_path / "presets.json"))
        d = _gamess.GamessSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_sto_forces_ngauss_3(self, dlg):
        dlg.basis_gbasis.setCurrentText("STO")
        assert dlg.basis_ngauss.value() == 3
        assert dlg.basis_ngauss.isEnabled()

    def test_back_to_n31_restores_ngauss_6(self, dlg):
        dlg.basis_gbasis.setCurrentText("STO")
        dlg.basis_gbasis.setCurrentText("N31")
        assert dlg.basis_ngauss.value() == 6

    def test_tzv_sets_ngauss_1(self, dlg):
        dlg.basis_gbasis.setCurrentText("TZV")
        assert dlg.basis_ngauss.value() == 1
        assert "GBASIS=TZV NGAUSS=1" in dlg.generate_content()


class TestGamessSaveFile:
    @pytest.fixture
    def dlg(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_gamess, "SETTINGS_FILE", str(tmp_path / "presets.json"))
        d = _gamess.GamessSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_save_writes_preview_content(self, dlg, monkeypatch, tmp_path):
        target = tmp_path / "job.inp"
        infos = []
        monkeypatch.setattr(
            _gamess.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(target), "GAMESS Input (*.inp *.src)")),
        )
        monkeypatch.setattr(
            _gamess.QMessageBox,
            "information",
            staticmethod(lambda *a, **k: infos.append(a)),
        )
        dlg.preview_text.setPlainText("MANUALLY EDITED CONTENT")
        dlg.save_file()
        assert target.read_text(encoding="utf-8") == "MANUALLY EDITED CONTENT"
        assert len(infos) == 1

    def test_save_cancelled_writes_nothing(self, dlg, monkeypatch, tmp_path):
        monkeypatch.setattr(
            _gamess.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        dlg.save_file()
        assert list(tmp_path.glob("*.inp")) == []


class TestGamessPresets:
    @pytest.fixture
    def settings_file(self, monkeypatch, tmp_path):
        path = tmp_path / "presets.json"
        monkeypatch.setattr(_gamess, "SETTINGS_FILE", str(path))
        return path

    @pytest.fixture
    def dlg(self, qapp, settings_file):
        d = _gamess.GamessSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_save_preset_persists_and_selects(self, dlg, settings_file, monkeypatch):
        import json

        monkeypatch.setattr(
            _gamess.QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("My Preset", True)),
        )
        dlg.run_type.setCurrentText("HESSIAN")
        dlg.mem_spin.setValue(250)
        dlg.save_preset_dialog()

        assert dlg.preset_combo.currentText() == "My Preset"
        on_disk = json.loads(settings_file.read_text(encoding="utf-8"))
        assert on_disk["My Preset"]["run_type"] == "HESSIAN"
        assert on_disk["My Preset"]["memory"] == 250

    def test_save_preset_cancelled_is_noop(self, dlg, settings_file, monkeypatch):
        monkeypatch.setattr(
            _gamess.QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("", False)),
        )
        dlg.save_preset_dialog()
        assert dlg.presets_data == {}
        assert not settings_file.exists()

    def test_new_dialog_loads_saved_preset_and_applies(
        self, qapp, settings_file, monkeypatch
    ):
        import json

        settings_file.write_text(
            json.dumps(
                {
                    "Fast": {
                        "run_type": "ENERGY",
                        "scf_type": "UHF",
                        "dft_type": "PBE0",
                        "nosym": False,
                        "gbasis": "N21",
                        "ngauss": 3,
                        "ndfunc": 0,
                        "npfunc": 1,
                        "memory": 42,
                    }
                }
            ),
            encoding="utf-8",
        )
        d = _gamess.GamessSetupDialog(parent=None, mol=None)
        try:
            assert d.preset_combo.currentText() == "Fast"
            assert d.run_type.currentText() == "ENERGY"
            assert d.scf_type.currentText() == "UHF"
            assert d.dft_type.currentText() == "PBE0"
            assert not d.chk_nosym.isChecked()
            assert d.mem_spin.value() == 42
            content = d.generate_content()
            assert "DFTTYP=PBE0" in content
            assert "NOSYM" not in content
        finally:
            d.destroy()

    def test_delete_preset_confirmed_removes_entry(
        self, dlg, settings_file, monkeypatch
    ):
        import json

        dlg.presets_data = {"Doomed": {"run_type": "ENERGY"}}
        dlg.update_preset_combo()
        assert dlg.preset_combo.currentText() == "Doomed"
        monkeypatch.setattr(
            _gamess.QMessageBox,
            "question",
            staticmethod(lambda *a, **k: _gamess.QMessageBox.StandardButton.Yes),
        )
        dlg.delete_preset()
        assert dlg.presets_data == {}
        assert dlg.preset_combo.count() == 0
        assert json.loads(settings_file.read_text(encoding="utf-8")) == {}

    def test_delete_preset_declined_keeps_entry(self, dlg, monkeypatch):
        dlg.presets_data = {"Kept": {"run_type": "ENERGY"}}
        dlg.update_preset_combo()
        monkeypatch.setattr(
            _gamess.QMessageBox,
            "question",
            staticmethod(lambda *a, **k: _gamess.QMessageBox.StandardButton.No),
        )
        dlg.delete_preset()
        assert "Kept" in dlg.presets_data
