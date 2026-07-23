"""
Headless GUI tests for the Complex Molecule Untangler plugin.

Covers: UntanglerDialog, run_plugin.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_complex_molecule_untangler.py
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

UNTANGLER_PATH = (
    PLUGINS_DIR / "Complex_Molecule_Untangler" / "complex_molecule_untangler.py"
)

with mock_chemistry_imports():
    _untangler = load_plugin_for_gui(UNTANGLER_PATH)


def _ctx(setting="MMFF_RDKIT") -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.get_setting.return_value = setting
    return ctx


@pytest.fixture
def no_msgbox(monkeypatch):
    """Replace the module's QMessageBox so warnings/info don't block offscreen."""
    box = MagicMock()
    monkeypatch.setattr(_untangler, "QMessageBox", box)
    return {_untangler.__name__: box}


# ===========================================================================
# UntanglerDialog  (visible plugin: "Complex Molecule Untangler")
# ===========================================================================


class TestUntanglerDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _untangler.UntanglerDialog(context=_ctx(), parent=None)
        yield d
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Complex Molecule Untangler"

    def test_force_field_options(self, dlg):
        texts = [dlg.combo_ff.itemText(i) for i in range(dlg.combo_ff.count())]
        assert texts == ["MMFF94", "UFF"]

    def test_default_ff_from_uff_setting(self, qapp):
        d = _untangler.UntanglerDialog(context=_ctx(setting="UFF"), parent=None)
        assert d.combo_ff.currentText() == "UFF"
        d.destroy()

    def test_iterations_defaults(self, dlg):
        assert dlg.spin_iter.value() == 500
        assert (dlg.spin_iter.minimum(), dlg.spin_iter.maximum()) == (100, 10000)
        assert dlg.spin_iter.singleStep() == 100

    def test_progress_bar_initial(self, dlg):
        assert dlg.pbar.value() == 0
        assert not dlg.pbar.isTextVisible()

    def test_run_button(self, dlg):
        assert dlg.btn_run.text() == "Untangle Molecule"
        assert dlg.btn_run.isEnabled()

    def test_run_without_molecule_warns(self, dlg, no_msgbox):
        dlg.run_untangle()
        no_msgbox[_untangler.__name__].warning.assert_called_once()
        assert dlg.btn_run.isEnabled()
        assert dlg.worker is None

    def test_on_finished_error_restores_button(self, dlg, no_msgbox):
        dlg.btn_run.setEnabled(False)
        dlg.btn_run.setText("Processing...")
        dlg.on_finished(None, "No rotatable bonds found.")
        assert dlg.btn_run.isEnabled()
        assert dlg.btn_run.text() == "Untangle Molecule"
        no_msgbox[_untangler.__name__].warning.assert_called_once()

    def test_on_finished_success_updates_molecule(self, dlg, no_msgbox):
        new_mol = MagicMock()
        dlg.on_finished(new_mol, "Processed 5 bonds.")
        assert dlg.context.current_mol is new_mol
        dlg.context.refresh_3d_view.assert_called_once()
        dlg.context.push_undo_checkpoint.assert_called_once()
        no_msgbox[_untangler.__name__].information.assert_called_once()

    def test_worker_stores_parameters(self, qapp):
        mol = MagicMock()
        w = _untangler.UntangleWorker(mol, max_iter=250, force_field="UFF")
        assert w.mol is mol
        assert w.max_iter == 250
        assert w.force_field == "UFF"

    def test_run_plugin_registers_new_dialog(self):
        ctx = _ctx()
        ctx.get_window.return_value = None
        _untangler.run_plugin(ctx)
        args = ctx.register_window.call_args[0]
        assert args[0] == "main_panel"
        assert isinstance(args[1], _untangler.UntanglerDialog)
        args[1].destroy()

    def test_run_plugin_raises_existing_window(self):
        ctx = _ctx()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _untangler.run_plugin(ctx)
        existing.show.assert_called_once()
        ctx.register_window.assert_not_called()

    def test_initialize_then_run_launches(self):
        ctx = _ctx()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        _untangler.initialize(ctx)
        _untangler.run(MagicMock())
        existing.show.assert_called_once()

    def test_run_with_no_launch_fn_uses_plugin_context(self, monkeypatch):
        """run(mw) falls back to PLUGIN_CONTEXT when _launch_fn is unset."""
        ctx = _ctx()
        existing = MagicMock()
        ctx.get_window.return_value = existing
        monkeypatch.setattr(_untangler, "_launch_fn", None)
        monkeypatch.setattr(_untangler, "PLUGIN_CONTEXT", ctx)
        host_mw = SimpleNamespace(host=MagicMock())
        _untangler.run(host_mw)
        existing.show.assert_called_once()

    def test_run_untangle_valid_molecule_starts_worker(self, dlg, monkeypatch):
        """Covers the happy path of run_untangle() with a real molecule (lines 252-265)."""
        fake_worker = MagicMock()
        worker_cls = MagicMock(return_value=fake_worker)
        monkeypatch.setattr(_untangler, "UntangleWorker", worker_cls)
        dlg.context.current_mol = "MOL"
        dlg.spin_iter.setValue(250)
        dlg.combo_ff.setCurrentText("UFF")
        dlg.run_untangle()
        worker_cls.assert_called_once_with("MOL", max_iter=250, force_field="UFF")
        assert not dlg.btn_run.isEnabled()
        assert dlg.btn_run.text() == "Processing..."
        assert dlg.pbar.isTextVisible()
        fake_worker.start.assert_called_once()
        assert dlg.worker is fake_worker


# ===========================================================================
# UntangleWorker.run() — real QThread method, executed synchronously.
# Chemistry libs (Chem/AllChem/rdMolTransforms/random) are monkeypatched on
# the module directly so the real source lines execute (not AST-extracted).
# ===========================================================================


class _FakeAtomGui:
    def __init__(self, idx, nbr_idxs):
        self._idx = idx
        self._nbrs = nbr_idxs

    def GetIdx(self):
        return self._idx

    def GetNeighbors(self):
        return [_FakeAtomGui(i, []) for i in self._nbrs]


class _FakeWorkMolGui:
    """Linear chain 0-1-2-3 with one rotatable central bond (1,2)."""

    def __init__(self, matches=((1, 2),)):
        self.matches = matches
        self.conf = MagicMock(name="conf")

    def GetSubstructMatches(self, patt):
        return self.matches

    def GetAtomWithIdx(self, idx):
        neighbors = {1: [0, 2], 2: [1, 3]}
        return _FakeAtomGui(idx, neighbors.get(idx, []))

    def GetConformer(self):
        return self.conf


def _patch_chem(monkeypatch, work_mol, ff=None, uff=None, props_ok=True):
    chem = MagicMock()
    chem.Mol.return_value = work_mol
    chem.MolFromSmarts.return_value = "SMARTS"

    allchem = MagicMock()
    allchem.MMFFGetMoleculeProperties.return_value = MagicMock() if props_ok else None
    allchem.MMFFGetMoleculeForceField.return_value = ff
    allchem.UFFGetMoleculeForceField.return_value = uff

    rdt = MagicMock()
    rdt.GetDihedralDeg.return_value = 10.0

    monkeypatch.setattr(_untangler, "Chem", chem)
    monkeypatch.setattr(_untangler, "AllChem", allchem)
    monkeypatch.setattr(_untangler, "rdMolTransforms", rdt)
    monkeypatch.setattr(
        _untangler,
        "random",
        SimpleNamespace(choice=lambda seq: seq[0], uniform=lambda a, b: 42.0),
    )
    return chem, allchem, rdt


def _collect_finished(worker):
    results = {}
    worker.finished.connect(lambda new_mol, msg: results.update(new_mol=new_mol, msg=msg))
    return results


class TestUntangleWorkerRunReal:
    def test_mmff_setup_failure_reports_hint(self, qapp, monkeypatch):
        work_mol = _FakeWorkMolGui()
        _patch_chem(monkeypatch, work_mol, ff=None, props_ok=False)
        w = _untangler.UntangleWorker(work_mol, max_iter=3, force_field="MMFF94")
        results = _collect_finished(w)
        w.run()
        assert results["new_mol"] is None
        assert "Try using UFF" in results["msg"]

    def test_uff_setup_failure_reports_message(self, qapp, monkeypatch):
        work_mol = _FakeWorkMolGui()
        _patch_chem(monkeypatch, work_mol, uff=None)
        w = _untangler.UntangleWorker(work_mol, max_iter=3, force_field="UFF")
        results = _collect_finished(w)
        w.run()
        assert results["new_mol"] is None
        assert "Try using UFF" not in results["msg"]

    def test_no_rotatable_bonds(self, qapp, monkeypatch):
        work_mol = _FakeWorkMolGui(matches=())
        ff = MagicMock()
        ff.CalcEnergy.side_effect = [100.0]
        _patch_chem(monkeypatch, work_mol, ff=ff)
        w = _untangler.UntangleWorker(work_mol, max_iter=3, force_field="MMFF94")
        results = _collect_finished(w)
        w.run()
        assert results["new_mol"] is None
        assert results["msg"] == "No rotatable bonds found."

    def test_no_dihedrals(self, qapp, monkeypatch):
        class _TwoAtomWorkMolGui(_FakeWorkMolGui):
            def __init__(self):
                super().__init__(matches=((0, 1),))

            def GetAtomWithIdx(self, idx):
                neighbors = {0: [1], 1: [0]}
                return _FakeAtomGui(idx, neighbors[idx])

        work_mol = _TwoAtomWorkMolGui()
        ff = MagicMock()
        ff.CalcEnergy.side_effect = [100.0]
        _patch_chem(monkeypatch, work_mol, ff=ff)
        w = _untangler.UntangleWorker(work_mol, max_iter=3, force_field="MMFF94")
        results = _collect_finished(w)
        w.run()
        assert results["new_mol"] is None
        assert results["msg"] == "Could not define dihedrals."

    def test_mmff_full_run_accept_and_revert(self, qapp, monkeypatch):
        work_mol = _FakeWorkMolGui()
        ff = MagicMock()
        # initial 100 -> iter1 90 (accept) -> iter2 95 (revert) -> iter3 80 (accept)
        ff.CalcEnergy.side_effect = [100.0, 90.0, 95.0, 80.0]
        _patch_chem(monkeypatch, work_mol, ff=ff)
        w = _untangler.UntangleWorker(work_mol, max_iter=3, force_field="MMFF94")
        results = _collect_finished(w)
        progress_calls = []
        w.progress.connect(progress_calls.append)
        w.run()
        assert results["new_mol"] is work_mol
        assert "Processed 1 bonds" in results["msg"]
        assert "80.00" in results["msg"]
        assert progress_calls == [1, 2, 3]

    def test_uff_full_run(self, qapp, monkeypatch):
        work_mol = _FakeWorkMolGui()
        uff = MagicMock()
        uff.CalcEnergy.side_effect = [100.0, 90.0]
        _, allchem, _ = _patch_chem(monkeypatch, work_mol, uff=uff)
        w = _untangler.UntangleWorker(work_mol, max_iter=1, force_field="UFF")
        results = _collect_finished(w)
        w.run()
        assert results["new_mol"] is work_mol
        allchem.UFFOptimizeMolecule.assert_called_once()

    def test_final_optimization_exception_is_silenced(self, qapp, monkeypatch):
        work_mol = _FakeWorkMolGui()
        ff = MagicMock()
        ff.CalcEnergy.side_effect = [100.0, 90.0]
        _, allchem, _ = _patch_chem(monkeypatch, work_mol, ff=ff)
        allchem.MMFFOptimizeMolecule.side_effect = RuntimeError("opt failed")
        w = _untangler.UntangleWorker(work_mol, max_iter=1, force_field="MMFF94")
        results = _collect_finished(w)
        w.run()
        assert results["new_mol"] is work_mol

    def test_exception_reported_via_finished(self, qapp, monkeypatch):
        work_mol = _FakeWorkMolGui()
        chem, _, _ = _patch_chem(monkeypatch, work_mol)
        chem.Mol.side_effect = RuntimeError("boom")
        w = _untangler.UntangleWorker(work_mol, max_iter=3, force_field="MMFF94")
        results = _collect_finished(w)
        w.run()
        assert results["new_mol"] is None
        assert results["msg"] == "boom"
