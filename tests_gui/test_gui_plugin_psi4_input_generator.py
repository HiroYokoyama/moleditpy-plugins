"""
Headless GUI tests for the Psi4 Input Generator plugin.

Covers: Psi4SetupDialog.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PSI4_PATH = PLUGINS_DIR / "Psi4_Input_Generator" / "psi4_input_generator.py"

with mock_chemistry_imports():
    _psi4 = load_plugin_for_gui(PSI4_PATH)


# ===========================================================================
# Psi4SetupDialog  (visible plugin: "Psi4 Input Generator")
# ===========================================================================


class TestPsi4SetupDialog:
    """Psi4SetupDialog with mol=None."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _psi4.Psi4SetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Psi4 Input Generator"

    def test_method_default(self, dlg):
        assert dlg.method_combo.currentText() == "b3lyp"

    def test_basis_default(self, dlg):
        assert dlg.basis_combo.currentText() == "def2-svp"

    def test_reference_default(self, dlg):
        assert dlg.ref_combo.currentText() == "rks"

    def test_mem_spin_default(self, dlg):
        assert dlg.mem_spin.value() == 2

    def test_thread_spin_default(self, dlg):
        assert dlg.thread_spin.value() == 4


class TestPsi4AutoReference:
    """update_auto_reference(): method/multiplicity drive the reference combo."""

    @pytest.fixture
    def dlg(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_psi4, "SETTINGS_FILE", str(tmp_path / "presets.json"))
        d = _psi4.Psi4SetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_wavefunction_singlet_gets_rhf(self, dlg):
        dlg.method_combo.setCurrentText("mp2")
        assert dlg.ref_combo.currentText() == "rhf"

    def test_wavefunction_open_shell_gets_uhf(self, dlg):
        dlg.method_combo.setCurrentText("mp2")
        dlg.mult_spin.setValue(3)
        assert dlg.ref_combo.currentText() == "uhf"

    def test_dft_open_shell_gets_uks(self, dlg):
        dlg.mult_spin.setValue(2)  # method stays b3lyp (DFT)
        assert dlg.ref_combo.currentText() == "uks"

    def test_dft_back_to_singlet_gets_rks(self, dlg):
        dlg.mult_spin.setValue(2)
        dlg.mult_spin.setValue(1)
        assert dlg.ref_combo.currentText() == "rks"

    def test_ccsd_t_singlet_gets_rhf(self, dlg):
        dlg.method_combo.setCurrentText("ccsd(t)")
        assert dlg.ref_combo.currentText() == "rhf"


class TestPsi4GenerateContent:
    @pytest.fixture
    def dlg(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_psi4, "SETTINGS_FILE", str(tmp_path / "presets.json"))
        d = _psi4.Psi4SetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_resources_lines(self, dlg):
        dlg.mem_spin.setValue(8)
        dlg.thread_spin.setValue(16)
        content = dlg.generate_content()
        assert "memory 8 GB" in content
        assert "set_num_threads(16)" in content

    def test_charge_mult_line_inside_molecule_block(self, dlg):
        dlg.charge_spin.setValue(-2)
        dlg.mult_spin.setValue(3)
        content = dlg.generate_content()
        mol_block = content.split("molecule {")[1].split("}")[0]
        assert "-2 3" in mol_block

    def test_basis_and_reference_in_set_block(self, dlg):
        dlg.basis_combo.setCurrentText("cc-pvtz")
        content = dlg.generate_content()
        assert "basis cc-pvtz" in content
        assert "reference rks" in content

    def test_task_line_combines_task_and_method(self, dlg):
        dlg.task_combo.setCurrentText("optimize")
        assert dlg.generate_content().endswith("optimize('b3lyp')")

    def test_filename_hint_used_for_output_file(self, dlg):
        content = dlg.generate_content(filename_hint="myjob")
        assert "psi4.core.set_output_file('myjob.out', False)" in content

    def test_preview_uses_placeholder_hint(self, dlg):
        assert (
            "psi4.core.set_output_file('[filename].out', False)"
            in dlg.preview_text.toPlainText()
        )


class TestPsi4SaveFile:
    @pytest.fixture
    def dlg(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(_psi4, "SETTINGS_FILE", str(tmp_path / "presets.json"))
        d = _psi4.Psi4SetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_save_rewrites_output_file_to_match_basename(
        self, dlg, monkeypatch, tmp_path
    ):
        target = tmp_path / "benzene_opt.dat"
        infos = []
        monkeypatch.setattr(
            _psi4.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(target), "PSI4 Input (*.dat *.in)")),
        )
        monkeypatch.setattr(
            _psi4.QMessageBox,
            "information",
            staticmethod(lambda *a, **k: infos.append(a)),
        )
        dlg.save_file()
        saved = target.read_text(encoding="utf-8")
        assert "psi4.core.set_output_file('benzene_opt.out', False)" in saved
        assert "[filename]" not in saved
        assert len(infos) == 1

    def test_save_prepends_output_line_when_missing(self, dlg, monkeypatch, tmp_path):
        target = tmp_path / "bare.dat"
        monkeypatch.setattr(
            _psi4.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(target), "")),
        )
        monkeypatch.setattr(
            _psi4.QMessageBox, "information", staticmethod(lambda *a, **k: None)
        )
        dlg.preview_text.setPlainText("energy('scf')")
        dlg.save_file()
        saved = target.read_text(encoding="utf-8")
        assert saved.startswith("psi4.core.set_output_file('bare.out', False)")
        assert "energy('scf')" in saved

    def test_save_cancelled_writes_nothing(self, dlg, monkeypatch, tmp_path):
        monkeypatch.setattr(
            _psi4.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        dlg.save_file()
        assert list(tmp_path.glob("*.dat")) == []


class TestPsi4Presets:
    @pytest.fixture
    def settings_file(self, monkeypatch, tmp_path):
        path = tmp_path / "presets.json"
        monkeypatch.setattr(_psi4, "SETTINGS_FILE", str(path))
        return path

    @pytest.fixture
    def dlg(self, qapp, settings_file):
        d = _psi4.Psi4SetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_save_preset_persists_and_selects(self, dlg, settings_file, monkeypatch):
        import json

        monkeypatch.setattr(
            _psi4.QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("HighAcc", True)),
        )
        dlg.method_combo.setCurrentText("ccsd(t)")
        dlg.thread_spin.setValue(8)
        dlg.save_preset_dialog()

        assert dlg.preset_combo.currentText() == "HighAcc"
        on_disk = json.loads(settings_file.read_text(encoding="utf-8"))
        assert on_disk["HighAcc"]["method"] == "ccsd(t)"
        assert on_disk["HighAcc"]["threads"] == 8

    def test_new_dialog_loads_saved_preset_and_applies(self, qapp, settings_file):
        import json

        settings_file.write_text(
            json.dumps(
                {
                    "Quick": {
                        "method": "scf",
                        "reference": "rhf",
                        "basis": "cc-pvdz",
                        "task": "energy",
                        "memory": 1,
                        "threads": 2,
                    }
                }
            ),
            encoding="utf-8",
        )
        d = _psi4.Psi4SetupDialog(parent=None, mol=None)
        try:
            assert d.preset_combo.currentText() == "Quick"
            assert d.method_combo.currentText() == "scf"
            assert d.basis_combo.currentText() == "cc-pvdz"
            assert d.mem_spin.value() == 1
            content = d.generate_content()
            assert "basis cc-pvdz" in content
            assert content.endswith("energy('scf')")
        finally:
            d.destroy()

    def test_delete_preset_confirmed_removes_entry(
        self, dlg, settings_file, monkeypatch
    ):
        dlg.presets_data = {"Doomed": {"method": "scf"}}
        dlg.update_preset_combo()
        monkeypatch.setattr(
            _psi4.QMessageBox,
            "question",
            staticmethod(lambda *a, **k: _psi4.QMessageBox.StandardButton.Yes),
        )
        dlg.delete_preset()
        assert dlg.presets_data == {}
        assert dlg.preset_combo.count() == 0
