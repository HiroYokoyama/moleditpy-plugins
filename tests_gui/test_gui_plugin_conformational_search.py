"""
Headless GUI tests for the Conformational Search plugin.

Covers: ConformerSearchDialog.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_conformational_search.py
"""

from __future__ import annotations

import contextlib
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

CONF_SEARCH_PATH = PLUGINS_DIR / "Conformational_Search" / "conf_search.py"

with mock_chemistry_imports():
    _conf_search = load_plugin_for_gui(CONF_SEARCH_PATH)


# Import the real rdkit up front so it is already cached in sys.modules
# before _mock_chemistry_keep_real_chem() snapshots it below — otherwise the
# chemistry-blocking meta path finder would intercept the first-ever
# `import rdkit` and hand back a MagicMock instead.
# Guarded so CI's bare test-gui job (only pytest+PyQt6 installed) skips this
# real-rdkit module instead of erroring at collection.
_Chem = pytest.importorskip("rdkit.Chem")
_AllChem = pytest.importorskip("rdkit.Chem.AllChem")
_Point3D = pytest.importorskip("rdkit.Geometry").Point3D


@contextlib.contextmanager
def _mock_chemistry_keep_real_chem():
    """Like mock_chemistry_imports(), but numpy/rdkit resolve to the real
    packages (both installed in this environment).

    Needed so the plugin's real ETKDG embedding / MMFF/UFF optimization
    actually runs instead of chasing MagicMock attribute chains.
    """
    keep_prefixes = ("numpy", "rdkit")
    real_mods = {
        k: v for k, v in sys.modules.items() if k.split(".")[0] in keep_prefixes
    }
    with mock_chemistry_imports():
        sys.modules.update(real_mods)
        yield


# Second module instance with real rdkit, used for tests that drive real
# conformer generation/optimization through the plugin's bound methods
# (separate from `_conf_search` above, which keeps everything mocked for
# the plain widget-construction tests).
with _mock_chemistry_keep_real_chem():
    _conf_rn = load_plugin_for_gui(CONF_SEARCH_PATH)


def _ctx_no_mol() -> MagicMock:
    """Context with no main window and no active molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


def _real_mol_with_conformer(smiles="CCCC", seed=42):
    """Real RDKit mol (with Hs) carrying exactly one embedded conformer."""
    mol = _Chem.AddHs(_Chem.MolFromSmiles(smiles))
    _AllChem.EmbedMolecule(mol, randomSeed=seed)
    return mol


def _no_block_msgbox(monkeypatch, mod):
    """QMessageBox.* pop up modal dialogs that block exec() under offscreen
    Qt with no user interaction; replace with recording stubs."""
    calls = {"warning": [], "information": [], "critical": []}
    for kind in calls:
        monkeypatch.setattr(
            mod.QMessageBox,
            kind,
            staticmethod(lambda *a, _k=kind, **kw: calls[_k].append(a)),
        )
    return calls


# ===========================================================================
# ConformerSearchDialog  (visible plugin: "Conformational Search")
# ===========================================================================


class TestConformerSearchDialog:
    """ConformerSearchDialog with no molecule loaded."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        d = _conf_search.ConformerSearchDialog(context=ctx, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Conformational Search & Preview"

    def test_table_has_two_columns(self, dlg):
        assert dlg.table.columnCount() == 2

    def test_table_column_headers(self, dlg):
        assert dlg.table.horizontalHeaderItem(0).text() == "Rank"
        assert dlg.table.horizontalHeaderItem(1).text() == "Energy (kcal/mol)"

    def test_table_initially_empty(self, dlg):
        assert dlg.table.rowCount() == 0

    def test_force_field_combo_has_options(self, dlg):
        texts = [dlg.combo_ff.itemText(i) for i in range(dlg.combo_ff.count())]
        assert "MMFF94" in texts
        assert "UFF" in texts

    def test_show_all_checkbox_unchecked_by_default(self, dlg):
        assert not dlg.cb_show_all.isChecked()

    def test_run_button_exists(self, dlg):
        assert dlg.btn_run.text() == "Run Search"

    def test_close_button_exists(self, dlg):
        assert dlg.btn_close.text() == "Close"

    def test_conformer_data_initially_empty(self, dlg):
        assert dlg.conformer_data == []


class TestConformerSearchDialogDefaultForceField:
    def test_uff_setting_selects_uff(self, qapp):
        ctx = _ctx_no_mol()
        ctx.get_setting.return_value = "UFF_RDKIT"
        d = _conf_search.ConformerSearchDialog(context=ctx, parent=None)
        assert d.combo_ff.currentText() == "UFF"
        d.destroy()

    def test_mmff_setting_selects_mmff94(self, qapp):
        ctx = _ctx_no_mol()
        ctx.get_setting.return_value = "MMFF_RDKIT"
        d = _conf_search.ConformerSearchDialog(context=ctx, parent=None)
        assert d.combo_ff.currentText() == "MMFF94"
        d.destroy()

    def test_no_setting_leaves_first_item_selected(self, qapp):
        ctx = _ctx_no_mol()
        ctx.get_setting.return_value = None
        d = _conf_search.ConformerSearchDialog(context=ctx, parent=None)
        assert d.combo_ff.currentText() == "MMFF94"
        d.destroy()


# ===========================================================================
# accept / reject / closeEvent — real dialog + real rdkit molecule
# ===========================================================================


class TestAcceptRejectCloseReal:
    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        ctx.current_mol = _real_mol_with_conformer()
        d = _conf_rn.ConformerSearchDialog(context=ctx, parent=None)
        yield d
        d.destroy()

    def test_accept_pushes_undo_checkpoint_and_unregisters_window(self, dlg):
        dlg.accept()
        dlg.context.push_undo_checkpoint.assert_called_once()
        dlg.context.register_window.assert_called_with("main_panel", None)

    def test_reject_restores_original_coords_and_refreshes_view(self, dlg):
        conf = dlg.target_mol.GetConformer()
        original = [conf.GetAtomPosition(i) for i in range(dlg.target_mol.GetNumAtoms())]
        # Mutate the live conformer to simulate a preview having changed it.
        conf.SetAtomPosition(0, _Point3D(9.0, 9.0, 9.0))
        dlg.reject()
        restored = conf.GetAtomPosition(0)
        assert restored.x == pytest.approx(original[0].x)
        assert restored.y == pytest.approx(original[0].y)
        assert restored.z == pytest.approx(original[0].z)
        dlg.context.refresh_3d_view.assert_called_once()
        dlg.context.register_window.assert_called_with("main_panel", None)

    def test_close_event_delegates_to_accept(self, dlg, monkeypatch):
        calls = []
        monkeypatch.setattr(dlg, "accept", lambda: calls.append("accept"))
        event = MagicMock()
        dlg.closeEvent(event)
        assert calls == ["accept"]
        event.ignore.assert_called_once()


class TestAcceptRejectNoMolReal:
    def test_accept_no_mol_does_not_push_undo(self, qapp):
        ctx = _ctx_no_mol()
        d = _conf_rn.ConformerSearchDialog(context=ctx, parent=None)
        d.accept()
        ctx.push_undo_checkpoint.assert_not_called()
        ctx.register_window.assert_called_with("main_panel", None)
        d.destroy()

    def test_reject_no_mol_does_not_raise(self, qapp):
        ctx = _ctx_no_mol()
        d = _conf_rn.ConformerSearchDialog(context=ctx, parent=None)
        d.reject()
        ctx.refresh_3d_view.assert_not_called()
        ctx.register_window.assert_called_with("main_panel", None)
        d.destroy()


# ===========================================================================
# run_search — real RDKit ETKDG embedding + MMFF/UFF optimization
# ===========================================================================


class TestRunSearchReal:
    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        ctx.current_mol = _real_mol_with_conformer()
        d = _conf_rn.ConformerSearchDialog(context=ctx, parent=None)
        yield d
        d.destroy()

    def test_mmff94_search_populates_table_and_sorts_ascending(self, dlg):
        dlg.combo_ff.setCurrentText("MMFF94")
        dlg.run_search()
        assert dlg.btn_run.isEnabled()
        assert len(dlg.results_raw) > 0
        energies = [e for e, _ in dlg.results_raw]
        assert energies == sorted(energies)
        assert dlg.table.rowCount() == len(dlg.conformer_data)
        assert dlg.temp_mol is not None
        assert "Showing" in dlg.lbl_info.text()

    def test_uff_search_populates_table(self, dlg):
        dlg.combo_ff.setCurrentText("UFF")
        dlg.run_search()
        assert len(dlg.results_raw) > 0
        assert dlg.table.rowCount() > 0

    def test_no_molecule_warns_and_skips(self, qapp, monkeypatch):
        ctx = _ctx_no_mol()
        d = _conf_rn.ConformerSearchDialog(context=ctx, parent=None)
        calls = _no_block_msgbox(monkeypatch, _conf_rn)
        d.run_search()
        assert len(calls["warning"]) == 1
        assert d.results_raw == []
        d.destroy()

    def test_embed_failure_warns_and_resets_label(self, dlg, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, _conf_rn)
        monkeypatch.setattr(
            _conf_rn.AllChem, "EmbedMultipleConfs", lambda mol, numConfs, params: []
        )
        dlg.run_search()
        assert len(calls["warning"]) == 1
        assert dlg.lbl_info.text() == "Failed."
        assert dlg.btn_run.isEnabled()

    def test_all_optimizations_fail_warns_and_resets(self, dlg, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, _conf_rn)
        monkeypatch.setattr(
            _conf_rn.AllChem, "MMFFOptimizeMolecule", lambda mol, confId=None: -1
        )
        dlg.combo_ff.setCurrentText("MMFF94")
        dlg.run_search()
        assert len(calls["warning"]) == 1
        assert dlg.results_raw == []
        assert dlg.btn_run.isEnabled()

    def test_exception_during_search_shows_critical(self, dlg, monkeypatch):
        calls = _no_block_msgbox(monkeypatch, _conf_rn)

        def _boom(*a, **k):
            raise RuntimeError("boom")

        monkeypatch.setattr(_conf_rn.AllChem, "ETKDGv3", _boom)
        dlg.run_search()
        assert len(calls["critical"]) == 1
        assert dlg.lbl_info.text() == "Error occurred."
        assert dlg.btn_run.isEnabled()

    def test_target_mol_refreshed_when_current_mol_changes(self, dlg):
        new_mol = _real_mol_with_conformer("CCO", seed=7)
        dlg.context.current_mol = new_mol
        dlg.run_search()
        assert dlg.target_mol is new_mol


# ===========================================================================
# apply_filter_and_update / update_table — real dialog after a real search
# ===========================================================================


class TestApplyFilterAndUpdateReal:
    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        ctx.current_mol = _real_mol_with_conformer()
        d = _conf_rn.ConformerSearchDialog(context=ctx, parent=None)
        d.run_search()
        yield d
        d.destroy()

    def test_default_filter_deduplicates_by_energy(self, dlg):
        assert dlg.cb_show_all.isChecked() is False
        assert len(dlg.conformer_data) <= len(dlg.results_raw)
        energies = [round(e, 4) for e, _ in dlg.conformer_data]
        assert len(energies) == len(set(energies))

    def test_show_all_reveals_every_result(self, dlg):
        dlg.cb_show_all.setChecked(True)
        assert dlg.conformer_data == dlg.results_raw
        assert dlg.table.rowCount() == len(dlg.results_raw)

    def test_table_rows_ranked_from_one(self, dlg):
        assert dlg.table.item(0, 0).text() == "1"

    def test_no_results_leaves_conformer_data_untouched(self, qapp):
        ctx = _ctx_no_mol()
        d = _conf_rn.ConformerSearchDialog(context=ctx, parent=None)
        d.apply_filter_and_update()  # results_raw is still empty -> early return
        assert d.conformer_data == []
        assert d.table.rowCount() == 0
        d.destroy()


# ===========================================================================
# preview_conformer — real dialog, real coordinate copy on table selection
# ===========================================================================


class TestPreviewConformerReal:
    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        ctx.current_mol = _real_mol_with_conformer()
        d = _conf_rn.ConformerSearchDialog(context=ctx, parent=None)
        d.run_search()
        yield d
        d.destroy()

    def test_selecting_row_copies_coordinates_and_refreshes_view(self, dlg):
        assert dlg.table.rowCount() > 0
        dlg.table.setCurrentCell(0, 0)
        assert dlg.context.current_mol is dlg.target_mol
        dlg.context.refresh_3d_view.assert_called()

    def test_atom_count_mismatch_shows_error(self, dlg):
        # Simulate the molecule having changed externally (fewer atoms).
        dlg.target_mol = _real_mol_with_conformer("CC", seed=1)
        dlg.table.setCurrentCell(0, 0)
        assert "Restart search" in dlg.lbl_info.text()
