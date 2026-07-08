"""
Tests for the Complex Molecule Untangler plugin.

Covers:
  1. initialize() must register at least one menu/export/plugin action
  2. No-molecule guard paths (run_plugin with mol=None)
  3. Dialog accept/reject round-trips
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
UNTANGLER_PATH = PLUGINS_DIR / "Complex_Molecule_Untangler" / "complex_molecule_untangler.py"


class TestComplexMoleculeUntangler:
    def test_initialize_stores_context(self):
        """initialize() must store the context in PLUGIN_CONTEXT and set _launch_fn."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx
            assert mod._launch_fn is not None

    def test_run_plugin_opens_or_raises_dialog(self):
        """run_plugin with a valid context should not raise."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            mod.initialize(ctx)
            # run_plugin creates a dialog — should not raise with mocked Qt
            mod.run_plugin(ctx)

    def test_run_plugin_reuses_existing_window(self):
        """If a window is already registered, run_plugin shows it without creating new."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            existing_win = MagicMock()
            ctx.get_window.return_value = existing_win
            mod.run_plugin(ctx)
            existing_win.show.assert_called_once()
            existing_win.raise_.assert_called_once()

    def test_run_plugin_registers_new_dialog_window(self):
        """run_plugin() registers the created dialog as 'main_panel' when no existing window."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            mod.run_plugin(ctx)
            ctx.register_window.assert_called_once()


class _FakeAtom2:
    def __init__(self, idx, neighbor_idxs):
        self._idx = idx
        self._nbrs = neighbor_idxs

    def GetIdx(self):
        return self._idx

    def GetNeighbors(self):
        return [_FakeAtom2(i, []) for i in self._nbrs]


class _FakeWorkMol:
    """Linear chain 0-1-2-3 with one rotatable central bond (1,2)."""

    def __init__(self, matches=((1, 2),)):
        self.matches = matches
        self.conf = MagicMock(name="conf")

    def GetSubstructMatches(self, patt):
        return self.matches

    def GetAtomWithIdx(self, idx):
        neighbors = {1: [0, 2], 2: [1, 3]}
        return _FakeAtom2(idx, neighbors.get(idx, []))

    def GetConformer(self):
        return self.conf


def _untangle_env(work_mol, energies, props_ok=True, ff=None):
    """Build globals for UntangleWorker.run with scripted force-field energies."""
    import logging

    chem = MagicMock()
    chem.Mol.return_value = work_mol
    chem.MolFromSmarts.return_value = "SMARTS"

    if ff is None:
        ff = MagicMock()
        ff.CalcEnergy.side_effect = list(energies)

    allchem = MagicMock()
    allchem.MMFFGetMoleculeProperties.return_value = (
        MagicMock() if props_ok else None
    )
    allchem.MMFFGetMoleculeForceField.return_value = ff

    rdt = MagicMock()
    rdt.GetDihedralDeg.return_value = 10.0

    rnd = SimpleNamespace(
        choice=lambda seq: seq[0],
        uniform=lambda a, b: 42.0,
    )

    globs = {
        "Chem": chem,
        "AllChem": allchem,
        "rdMolTransforms": rdt,
        "random": rnd,
        "logging": logging,
    }
    fn = extract_function(UNTANGLER_PATH, "UntangleWorker", "run", globs)
    return fn, globs, ff, rdt


def _untangle_self(mol="MOL", max_iter=3, force_field="MMFF94"):
    return SimpleNamespace(
        mol=mol,
        max_iter=max_iter,
        force_field=force_field,
        progress=MagicMock(),
        finished=MagicMock(),
    )


class TestUntangleWorkerRun:
    def test_ff_setup_failure_reports_error(self):
        fn, globs, _, _ = _untangle_env(_FakeWorkMol(), [], props_ok=False)
        globs["AllChem"].MMFFGetMoleculeForceField.return_value = None
        s = _untangle_self()
        fn(s)
        new_mol, msg = s.finished.emit.call_args[0]
        assert new_mol is None
        assert "Could not setup Force Field" in msg
        assert "Try using UFF" in msg  # MMFF-specific hint

    def test_no_rotatable_bonds_reports_message(self):
        fn, _, _, _ = _untangle_env(_FakeWorkMol(matches=()), [100.0])
        s = _untangle_self()
        fn(s)
        new_mol, msg = s.finished.emit.call_args[0]
        assert new_mol is None
        assert msg == "No rotatable bonds found."

    def test_improvement_accepted_worsening_reverted(self):
        work_mol = _FakeWorkMol()
        # initial 100 -> iter1 90 (accept) -> iter2 95 (revert) -> iter3 80 (accept)
        fn, globs, ff, rdt = _untangle_env(work_mol, [100.0, 90.0, 95.0, 80.0])
        s = _untangle_self(max_iter=3)
        fn(s)

        set_calls = rdt.SetDihedralDeg.call_args_list
        # 3 rotations + 1 revert (iter2 back to old angle 10.0)
        assert len(set_calls) == 4
        assert set_calls[0].args == (work_mol.conf, 0, 1, 2, 3, 42.0)
        assert set_calls[2].args == (work_mol.conf, 0, 1, 2, 3, 10.0)  # revert

        new_mol, msg = s.finished.emit.call_args[0]
        assert new_mol is work_mol
        assert "Processed 1 bonds" in msg
        assert "80.00" in msg  # best energy reported

    def test_progress_emitted_each_iteration(self):
        fn, _, _, _ = _untangle_env(_FakeWorkMol(), [100.0, 90.0, 80.0, 70.0])
        s = _untangle_self(max_iter=3)
        fn(s)
        assert s.progress.emit.call_args_list == [((1,),), ((2,),), ((3,),)]

    def test_final_optimization_uses_selected_ff(self):
        fn, globs, _, _ = _untangle_env(_FakeWorkMol(), [100.0, 90.0])
        s = _untangle_self(max_iter=1)
        fn(s)
        globs["AllChem"].MMFFOptimizeMolecule.assert_called_once()

    def test_uff_path_builds_uff_force_field(self):
        work_mol = _FakeWorkMol()
        fn, globs, _, _ = _untangle_env(work_mol, [])
        uff = MagicMock()
        uff.CalcEnergy.side_effect = [100.0, 90.0]
        globs["AllChem"].UFFGetMoleculeForceField.return_value = uff
        s = _untangle_self(max_iter=1, force_field="UFF")
        fn(s)
        globs["AllChem"].UFFGetMoleculeForceField.assert_called_once()
        globs["AllChem"].UFFOptimizeMolecule.assert_called_once()
        new_mol, _ = s.finished.emit.call_args[0]
        assert new_mol is work_mol

    def test_exception_reported_via_finished(self):
        fn, globs, _, _ = _untangle_env(_FakeWorkMol(), [])
        globs["Chem"].Mol.side_effect = RuntimeError("boom")
        s = _untangle_self()
        fn(s)
        new_mol, msg = s.finished.emit.call_args[0]
        assert new_mol is None
        assert msg == "boom"
