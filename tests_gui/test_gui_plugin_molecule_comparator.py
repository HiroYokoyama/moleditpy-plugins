"""
Headless GUI tests for the Molecule Comparator plugin.

Covers: MoleculeComparator.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_molecule_comparator.py
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

COMPARATOR_PATH = PLUGINS_DIR / "Molecule_Comparator" / "molecule_comparator.py"

with mock_chemistry_imports():
    _comparator = load_plugin_for_gui(COMPARATOR_PATH)


def _ctx_no_mol() -> MagicMock:
    """Context with no main window and no active molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# MoleculeComparator  (visible plugin: "Molecule Comparator")
# ===========================================================================


class TestMoleculeComparator:
    """MoleculeComparator with no main window."""

    @pytest.fixture
    def comp(self, qapp):
        ctx = _ctx_no_mol()
        w = _comparator.MoleculeComparator(context=ctx)
        yield w
        w.destroy()

    def test_creates_without_error(self, comp):
        assert comp is not None

    def test_window_title(self, comp):
        assert comp.windowTitle() == "Molecule Comparator"

    def test_mol_list_initially_empty(self, comp):
        assert comp.mol_list.count() == 0

    def test_molecules_list_initially_empty(self, comp):
        assert comp.molecules == []

    def test_style_combo_default_sticks(self, comp):
        assert comp.combo_style.currentText() == "Sticks"

    def test_style_combo_has_four_options(self, comp):
        assert comp.combo_style.count() == 4

    def test_align_method_default(self, comp):
        assert comp.combo_align_method.currentText() == "Substructure (MCS)"

    def test_ignore_hydrogens_unchecked_by_default(self, comp):
        assert not comp.check_ignore_hs.isChecked()

    def test_wireframe_lighting_unchecked_by_default(self, comp):
        assert not comp.check_wireframe_lighting.isChecked()

    def test_add_current_button_exists(self, comp):
        assert comp.btn_add_current.text() == "Add Current"

    def test_align_button_exists(self, comp):
        assert comp.btn_align.text() == "Align & Calculate RMSD"


# ===========================================================================
# Helpers: comparator with visualization side effects stubbed out
# ===========================================================================


def _make_comp(qapp_unused=None):
    ctx = _ctx_no_mol()
    w = _comparator.MoleculeComparator(context=ctx)
    w.update_visualization = MagicMock()
    w.update_wireframe_lighting = MagicMock()
    w.reset_view = MagicMock()
    return w, ctx


def _fake_entry_mol(name=None):
    mol = MagicMock()
    mol.HasProp.return_value = name is not None
    mol.GetProp.return_value = name
    return mol


def _add_mols(w, n):
    for i in range(n):
        w.molecules.append(
            {
                "name": f"M{i}",
                "mol": MagicMock(),
                "color": "#ff0000",
                "scope": "Carbon Only",
                "rms": None,
            }
        )
    w.update_list()


# ===========================================================================
# Molecule list management
# ===========================================================================


class TestComparatorListManagement:
    def test_add_current_without_molecule_warns(self, qapp):
        w, ctx = _make_comp()
        with pytest.MonkeyPatch.context() as mp:
            warned = MagicMock()
            mp.setattr(_comparator.QMessageBox, "warning", warned)
            w.add_current_molecule()
            warned.assert_called_once()
        assert w.molecules == []
        w.destroy()

    def test_add_current_uses_mol_name_prop(self, qapp):
        w, ctx = _make_comp()
        ctx.current_molecule = _fake_entry_mol("benzene")
        w.add_current_molecule()
        assert w.molecules[0]["name"] == "benzene"
        w.destroy()

    def test_add_current_falls_back_to_numbered_name(self, qapp):
        w, ctx = _make_comp()
        ctx.current_molecule = _fake_entry_mol(None)
        w.add_current_molecule()
        assert w.molecules[0]["name"] == "Mol 1"
        w.destroy()

    def test_add_current_cycles_default_colors(self, qapp):
        w, ctx = _make_comp()
        ctx.current_molecule = _fake_entry_mol(None)
        w.add_current_molecule()
        w.add_current_molecule()
        colors = _comparator.DEFAULT_COLORS
        assert w.molecules[0]["color"] == colors[0]
        assert w.molecules[1]["color"] == colors[1]
        w.destroy()

    def test_update_list_marks_reference_and_selects_last(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 3)
        assert w.mol_list.item(0).text() == "M0 (Ref)"
        assert w.mol_list.item(1).text() == "M1"
        assert w.mol_list.currentRow() == 2
        w.destroy()

    def test_set_as_reference_moves_selected_to_top(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 3)
        w.mol_list.setCurrentRow(2)
        w.set_as_reference()
        assert [m["name"] for m in w.molecules] == ["M2", "M0", "M1"]
        assert w.mol_list.item(0).text() == "M2 (Ref)"
        w.destroy()

    def test_set_as_reference_noop_for_reference_row(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 2)
        w.mol_list.setCurrentRow(0)
        w.set_as_reference()
        assert [m["name"] for m in w.molecules] == ["M0", "M1"]
        w.destroy()

    def test_remove_molecule_removes_current_row(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 2)
        w.mol_list.setCurrentRow(0)
        w.remove_molecule()
        assert [m["name"] for m in w.molecules] == ["M1"]
        w.destroy()

    def test_remove_molecule_noop_without_selection(self, qapp):
        w, ctx = _make_comp()
        w.remove_molecule()
        assert w.molecules == []
        w.destroy()


# ===========================================================================
# Selection-driven UI state
# ===========================================================================


class TestComparatorSelectionState:
    def test_invalid_row_disables_color_and_scope(self, qapp):
        w, ctx = _make_comp()
        w.on_selection_changed(-1)
        assert not w.btn_color.isEnabled()
        assert not w.combo_scope.isEnabled()
        w.destroy()

    def test_valid_row_enables_and_reflects_scope(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 1)
        w.molecules[0]["scope"] = "All Atoms"
        w.on_selection_changed(0)
        assert w.btn_color.isEnabled()
        assert w.combo_scope.isEnabled()
        assert w.combo_scope.currentText() == "All Atoms"
        w.destroy()

    def test_change_scope_updates_entry(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 1)
        w.mol_list.setCurrentRow(0)
        w.combo_scope.setCurrentText("All Atoms")
        w.change_scope()
        assert w.molecules[0]["scope"] == "All Atoms"
        w.destroy()

    def test_text_color_black_on_light_white_on_dark(self, qapp):
        w, ctx = _make_comp()
        assert w._get_text_color_for_background("#ffffff") == "black"
        assert w._get_text_color_for_background("#000000") == "white"
        assert w._get_text_color_for_background("garbage") == "black"
        w.destroy()


# ===========================================================================
# Results table / clipboard
# ===========================================================================


class TestComparatorResults:
    def test_results_table_formats_rms_states(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 3)
        w.molecules[0]["rms"] = 0.0
        w.molecules[1]["rms"] = -1.0
        w.molecules[2]["rms"] = None
        w.update_results_table()
        assert w.table_results.item(0, 1).text() == "0.0000"
        assert w.table_results.item(1, 1).text() == "N/A"
        assert w.table_results.item(2, 1).text() == "-"
        w.destroy()

    def test_copy_results_noop_when_empty(self, qapp):
        from PyQt6.QtWidgets import QApplication

        w, ctx = _make_comp()
        QApplication.clipboard().setText("sentinel")
        w.copy_results_to_clipboard()
        assert QApplication.clipboard().text() == "sentinel"
        w.destroy()

    def test_copy_results_writes_tsv(self, qapp):
        from PyQt6.QtWidgets import QApplication

        w, ctx = _make_comp()
        _add_mols(w, 2)
        w.molecules[0]["rms"] = 0.0
        w.molecules[1]["rms"] = 1.23456
        w.copy_results_to_clipboard()
        text = QApplication.clipboard().text()
        assert text.splitlines() == [
            "Molecule\tRMSD (Å)",
            "M0\t0.0000",
            "M1\t1.2346",
        ]
        ctx.show_status_message.assert_called_with(
            "Results copied to clipboard", 3000
        )
        w.destroy()

    def test_run_alignment_needs_two_molecules(self, qapp):
        w, ctx = _make_comp()
        _add_mols(w, 1)
        w.run_alignment()
        assert not hasattr(w, "progress_dialog")
        w.destroy()


# ===========================================================================
# Close / cleanup behavior
# ===========================================================================


class TestComparatorCleanup:
    def _closable(self):
        w, ctx = _make_comp()
        w.mw = MagicMock()  # cleanup path touches mw.view_3d_manager
        w.exit_3d_only_mode = MagicMock()
        ctx.draw_molecule_3d.reset_mock()  # construction already draws once
        return w, ctx

    def test_cleanup_redraws_current_molecule(self, qapp):
        w, ctx = self._closable()
        ctx.current_molecule = MagicMock()
        w.cleanup_and_close()
        ctx.draw_molecule_3d.assert_called_once_with(ctx.current_molecule)
        w.exit_3d_only_mode.assert_called_once()
        w.destroy()

    def test_cleanup_clears_plotter_without_molecule(self, qapp):
        w, ctx = self._closable()
        ctx.current_molecule = None
        w.cleanup_and_close()
        ctx.plotter.clear.assert_called_once()
        w.destroy()

    def test_cleanup_resets_color_overrides(self, qapp):
        w, ctx = self._closable()
        w.mw.view_3d_manager._plugin_color_overrides = {"stale": 1}
        w.cleanup_and_close()
        assert w.mw.view_3d_manager._plugin_color_overrides == {}
        w.destroy()

    def test_close_plugin_hides_window(self, qapp):
        w, ctx = self._closable()
        w.show()
        w.close_plugin()
        assert w.isHidden()
        w.destroy()


# ===========================================================================
# initialize() — document reset closes the window
# ===========================================================================


class TestComparatorInitialize:
    def test_registers_document_reset_handler(self, qapp):
        ctx = _ctx_no_mol()
        _comparator.initialize(ctx)
        ctx.register_document_reset_handler.assert_called_once()

    def test_reset_closes_open_comparator_window(self, qapp):
        ctx = _ctx_no_mol()
        handlers = []
        ctx.register_document_reset_handler.side_effect = handlers.append
        _comparator.initialize(ctx)
        mw = MagicMock()
        ctx.get_main_window.return_value = mw
        handlers[0]()
        mw.molecule_comparator_window.close.assert_called_once()


# ===========================================================================
# run(mw) module function
# ===========================================================================


class TestRunFunction:
    def test_run_noop_without_context(self, qapp):
        _comparator.PLUGIN_CONTEXT = None
        mw = MagicMock(spec=[])
        _comparator.run(mw)  # no window attribute created

    def _real_mw_namespace(self):
        # Must be a real QWidget: MoleculeComparator.__init__ passes it as the
        # Qt parent. enter_3d_only_mode() unconditionally touches
        # mw.ui_manager / mw.init_manager (no hasattr guard on those two), so
        # both must exist too.
        from PyQt6.QtWidgets import QWidget

        mw = QWidget()
        mw.ui_manager = SimpleNamespace()
        mw.init_manager = SimpleNamespace()
        mw.view_3d_manager = MagicMock()  # touched by cleanup_and_close/change_style
        return mw

    def test_run_unwraps_host(self, qapp):
        ctx = _ctx_no_mol()
        ctx.get_main_window.return_value = self._real_mw_namespace()
        _comparator.PLUGIN_CONTEXT = ctx
        inner_mw = SimpleNamespace()
        wrapper = SimpleNamespace(host=inner_mw)
        _comparator.run(wrapper)
        assert hasattr(inner_mw, "molecule_comparator_window")
        win = inner_mw.molecule_comparator_window
        assert win.isVisible()
        win.destroy()

    def test_run_toggles_existing_window_closed(self, qapp):
        ctx = _ctx_no_mol()
        ctx.get_main_window.return_value = self._real_mw_namespace()
        _comparator.PLUGIN_CONTEXT = ctx
        mw = SimpleNamespace()
        _comparator.run(mw)  # creates + shows
        win = mw.molecule_comparator_window
        assert win.isVisible()
        _comparator.run(mw)  # second call closes it
        assert win.isHidden()
        win.destroy()


# ===========================================================================
# AlignmentWorker.run() — real QThread instance, run() called synchronously.
# rdkit is fully MagicMocked at the module level here, so Chem/AllChem/rdFMCS
# attributes are monkeypatched per-test with plain callables returning real
# python values so the arithmetic comparisons (`rms < best_rms`) behave.
# ===========================================================================


class _WAtom:
    def __init__(self, idx, z=6):
        self._idx = idx
        self._z = z

    def GetIdx(self):
        return self._idx

    def GetAtomicNum(self):
        return self._z


class _WMol:
    """Minimal fake rdkit Mol for AlignmentWorker.run()."""

    def __init__(self, atoms, matches=()):
        self._atoms = list(atoms)
        self._matches = matches

    def GetAtoms(self):
        return list(self._atoms)

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetSubstructMatches(self, patt, uniquify=False):
        return self._matches


def _run_worker(worker):
    received = {}
    errors = {}
    worker.finished_signal.connect(lambda r: received.setdefault("r", r))
    worker.error_signal.connect(lambda m: errors.setdefault("m", m))
    worker.run()
    return received.get("r"), errors.get("m")


class TestAlignmentWorkerAtomIds:
    def test_success_updates_mol_and_rms(self, qapp, monkeypatch):
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", lambda *a, **k: 0.42)
        ref = _WMol([_WAtom(0), _WAtom(1), _WAtom(2)])
        probe = _WMol([_WAtom(0), _WAtom(1), _WAtom(2)])
        worker = _comparator.AlignmentWorker(ref, [(1, probe)], "Atom IDs", False)
        results, err = _run_worker(worker)
        assert err is None
        assert results[0]["rms"] == 0.42
        assert results[0]["index"] == 1

    def test_ignore_hs_uses_heavy_atom_map(self, qapp, monkeypatch):
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", lambda *a, **k: 0.1)
        ref = _WMol([_WAtom(0, 6), _WAtom(1, 1), _WAtom(2, 8)])
        probe = _WMol([_WAtom(0, 6), _WAtom(1, 1), _WAtom(2, 8)])
        worker = _comparator.AlignmentWorker(ref, [(0, probe)], "Atom IDs", True)
        results, _err = _run_worker(worker)
        assert results[0]["rms"] == 0.1

    def test_mismatched_count_skips_alignment(self, qapp, monkeypatch):
        align = MagicMock()
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", align)
        ref = _WMol([_WAtom(0), _WAtom(1)])
        probe = _WMol([_WAtom(0), _WAtom(1), _WAtom(2)])
        worker = _comparator.AlignmentWorker(ref, [(0, probe)], "Atom IDs", False)
        results, _err = _run_worker(worker)
        align.assert_not_called()
        assert results[0]["rms"] == -1.0

    def test_runtime_error_is_silenced(self, qapp, monkeypatch):
        def _boom(*a, **k):
            raise RuntimeError("no match")

        monkeypatch.setattr(_comparator.AllChem, "AlignMol", _boom)
        ref = _WMol([_WAtom(0), _WAtom(1)])
        probe = _WMol([_WAtom(0), _WAtom(1)])
        worker = _comparator.AlignmentWorker(ref, [(0, probe)], "Atom IDs", False)
        results, err = _run_worker(worker)
        assert err is None
        assert results[0]["rms"] == -1.0


class TestAlignmentWorkerMCS:
    def _patch_fmcs(self, monkeypatch, num_atoms, ref_matches, probe_matches):
        # AlignmentWorker.run() does a *local* `from rdkit.Chem import rdFMCS`
        # inside the MCS branch. Since rdkit is genuinely installed in this
        # environment and mock_chemistry_imports() only shadows it during the
        # plugin's module-load, that statement re-resolves against the real
        # installed `rdkit.Chem` package at call time -- patching the frozen
        # `_comparator.Chem` mock (used for the module-level `Chem` global
        # elsewhere) has no effect on it. Patch the real package instead.
        # importorskip so CI test-gui (no rdkit installed) skips these MCS
        # tests instead of failing; local runs still exercise them.
        _real_rdkit_chem = pytest.importorskip("rdkit.Chem")

        monkeypatch.setattr(
            _real_rdkit_chem,
            "rdFMCS",
            SimpleNamespace(
                FindMCS=lambda mols, **kw: SimpleNamespace(
                    numAtoms=num_atoms, smartsString="C"
                )
            ),
            raising=False,
        )
        monkeypatch.setattr(_comparator.Chem, "MolFromSmarts", lambda s: "PATT")

    def test_no_mcs_match_leaves_rms_sentinel(self, qapp, monkeypatch):
        self._patch_fmcs(monkeypatch, num_atoms=0, ref_matches=(), probe_matches=())
        ref = _WMol([_WAtom(0), _WAtom(1)], matches=())
        probe = _WMol([_WAtom(0), _WAtom(1)], matches=())
        worker = _comparator.AlignmentWorker(
            ref, [(0, probe)], "Substructure (MCS)", False
        )
        results, _err = _run_worker(worker)
        assert results[0]["rms"] == -1.0

    def test_best_match_selected_without_ignore_hs(self, qapp, monkeypatch):
        self._patch_fmcs(monkeypatch, num_atoms=3, ref_matches=None, probe_matches=None)
        ref = _WMol([_WAtom(0), _WAtom(1), _WAtom(2)], matches=((0, 1, 2),))
        probe = _WMol(
            [_WAtom(0), _WAtom(1), _WAtom(2)], matches=((0, 1, 2), (2, 1, 0))
        )
        rms_values = iter([0.9, 0.3])
        monkeypatch.setattr(
            _comparator.AllChem, "AlignMol", lambda *a, **k: next(rms_values)
        )
        worker = _comparator.AlignmentWorker(
            ref, [(0, probe)], "Substructure (MCS)", False
        )
        results, _err = _run_worker(worker)
        assert results[0]["rms"] == pytest.approx(0.3)

    def test_ignore_hs_full_map_transform_applied(self, qapp, monkeypatch):
        self._patch_fmcs(monkeypatch, num_atoms=2, ref_matches=None, probe_matches=None)
        # Match indices are positions within the heavy-atom-only substructure
        # (as RemoveHs would produce), i.e. 0-based positions into p_heavy /
        # r_heavy — not original full-molecule atom indices.
        ref = _WMol(
            [_WAtom(0, 6), _WAtom(1, 1), _WAtom(2, 8)], matches=((0, 1),)
        )
        probe = _WMol(
            [_WAtom(0, 6), _WAtom(1, 1), _WAtom(2, 8)], matches=((0, 1),)
        )
        monkeypatch.setattr(_comparator.Chem, "RemoveHs", lambda m: m)
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", lambda *a, **k: 0.05)
        worker = _comparator.AlignmentWorker(
            ref, [(0, probe)], "Substructure (MCS)", True
        )
        results, _err = _run_worker(worker)
        assert results[0]["rms"] == pytest.approx(0.05)

    def test_ignore_hs_transform_exception_keeps_rms_but_not_mol(
        self, qapp, monkeypatch
    ):
        self._patch_fmcs(monkeypatch, num_atoms=2, ref_matches=None, probe_matches=None)
        ref = _WMol([_WAtom(0, 6), _WAtom(1, 1)], matches=((0,),))
        probe = _WMol([_WAtom(0, 6), _WAtom(1, 1)], matches=((0,),))

        calls = {"n": 0}

        def _align(probe_arg, ref_arg, atomMap, reflect=False):
            calls["n"] += 1
            if calls["n"] == 1:
                return 0.5
            raise ValueError("full map alignment failed")

        monkeypatch.setattr(_comparator.Chem, "RemoveHs", lambda m: m)
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", _align)
        worker = _comparator.AlignmentWorker(
            ref, [(0, probe)], "Substructure (MCS)", True
        )
        results, err = _run_worker(worker)
        assert err is None
        # best_rms was still recorded even though the transform application
        # raised and best_transform stayed None.
        assert results[0]["rms"] == pytest.approx(0.5)

    def test_max_combinations_limit_aborts_search(self, qapp, monkeypatch, capsys):
        self._patch_fmcs(monkeypatch, num_atoms=2, ref_matches=None, probe_matches=None)
        ref = _WMol([_WAtom(0), _WAtom(1)], matches=((0, 1), (1, 0)))
        probe = _WMol(
            [_WAtom(0), _WAtom(1)], matches=((0, 1), (1, 0), (0, 0))
        )
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", lambda *a, **k: 1.0)
        worker = _comparator.AlignmentWorker(
            ref, [(0, probe)], "Substructure (MCS)", False
        )
        worker.MAX_COMBINATIONS = 2
        results, _err = _run_worker(worker)
        assert results[0]["rms"] == pytest.approx(1.0)
        out = capsys.readouterr().out
        assert "Too many MCS matches" in out

    def test_progress_emitted_for_many_combinations(self, qapp, monkeypatch):
        self._patch_fmcs(monkeypatch, num_atoms=2, ref_matches=None, probe_matches=None)
        ref_matches = tuple((0, 1) for _ in range(5))
        probe_matches = tuple((0, 1) for _ in range(12))  # 60 combos > 50
        ref = _WMol([_WAtom(0), _WAtom(1)], matches=ref_matches)
        probe = _WMol([_WAtom(0), _WAtom(1)], matches=probe_matches)
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", lambda *a, **k: 1.0)
        worker = _comparator.AlignmentWorker(
            ref, [(0, probe)], "Substructure (MCS)", False
        )
        progress_msgs = []
        worker.progress.connect(lambda i, t, msg: progress_msgs.append(msg))
        worker.run()
        assert any("%" in m for m in progress_msgs)

    def test_interruption_breaks_inner_and_outer_loop(self, qapp, monkeypatch):
        self._patch_fmcs(monkeypatch, num_atoms=2, ref_matches=None, probe_matches=None)
        ref = _WMol([_WAtom(0), _WAtom(1)], matches=((0, 1),))
        probe = _WMol([_WAtom(0), _WAtom(1)], matches=((0, 1), (1, 0)))
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", lambda *a, **k: 1.0)
        worker = _comparator.AlignmentWorker(
            ref, [(0, probe)], "Substructure (MCS)", False
        )

        state = {"n": 0}

        def _fake_check():
            state["n"] += 1
            return state["n"] > 2

        worker.isInterruptionRequested = _fake_check
        results, _err = _run_worker(worker)
        assert results[0]["rms"] == pytest.approx(1.0)

    def test_outer_interruption_stops_before_targets(self, qapp, monkeypatch):
        ref = _WMol([_WAtom(0)])
        probe = _WMol([_WAtom(0)])
        worker = _comparator.AlignmentWorker(
            ref, [(0, probe)], "Substructure (MCS)", False
        )
        worker.isInterruptionRequested = lambda: True
        results, _err = _run_worker(worker)
        assert results == []

    def test_exception_reports_via_error_signal(self, qapp):
        ref = _WMol([_WAtom(0)])
        worker = _comparator.AlignmentWorker(ref, [(0, None)], "Atom IDs", False)
        results, err = _run_worker(worker)
        assert results is None
        assert err is not None


# ===========================================================================
# run_alignment / update_progress / on_alignment_finished / on_alignment_error
# — full plumbing, worker.start() monkeypatched to run() synchronously so no
# real QThread is spawned (repo convention, see symmetry_analyzer tests).
# ===========================================================================


class TestRunAlignmentIntegration:
    def _comp_with_two_mols(self, monkeypatch):
        monkeypatch.setattr(
            _comparator.AlignmentWorker, "start", _comparator.AlignmentWorker.run
        )
        w, ctx = _make_comp()
        w.molecules = [
            {
                "name": "ref",
                "mol": _WMol([_WAtom(0), _WAtom(1)]),
                "color": "#ff0000",
                "scope": "Carbon Only",
                "rms": None,
            },
            {
                "name": "probe",
                "mol": _WMol([_WAtom(0), _WAtom(1)]),
                "color": "#00ff00",
                "scope": "Carbon Only",
                "rms": None,
            },
        ]
        w.combo_align_method.setCurrentText("Atom IDs")
        return w, ctx

    def test_needs_two_molecules_early_return(self, qapp):
        w, _ctx = _make_comp()
        _add_mols(w, 1)
        w.run_alignment()
        assert not hasattr(w, "progress_dialog")
        w.destroy()

    def test_full_run_updates_results_table(self, qapp, monkeypatch):
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", lambda *a, **k: 0.77)
        w, ctx = self._comp_with_two_mols(monkeypatch)
        w.run_alignment()
        assert w.molecules[0]["rms"] == 0.0
        assert w.molecules[1]["rms"] == 0.77
        assert w.table_results.item(1, 1).text() == "0.7700"
        w.destroy()

    def test_error_signal_shows_critical_dialog(self, qapp, monkeypatch):
        monkeypatch.setattr(_comparator.AllChem, "AlignMol", lambda *a, **k: 0.5)
        w, ctx = self._comp_with_two_mols(monkeypatch)
        w.molecules[1]["mol"] = None  # forces AttributeError -> error_signal
        critical = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "critical", critical)
        w.run_alignment()
        critical.assert_called_once()
        w.destroy()

    def test_update_progress_ignored_when_cancelled(self, qapp):
        w, _ctx = _make_comp()
        _add_mols(w, 2)
        w.progress_dialog = MagicMock()
        w.progress_dialog.wasCanceled.return_value = True
        w.update_progress(1, 2, "msg")
        w.progress_dialog.setValue.assert_not_called()
        w.destroy()

    def test_update_progress_updates_dialog(self, qapp):
        w, _ctx = _make_comp()
        w.progress_dialog = MagicMock()
        w.progress_dialog.wasCanceled.return_value = False
        w.update_progress(1, 5, "halfway")
        w.progress_dialog.setValue.assert_called_once_with(1)
        w.progress_dialog.setLabelText.assert_called_once_with("halfway")
        w.destroy()


# ===========================================================================
# change_color / change_scope early-return / redraw_visualization
# ===========================================================================


class TestChangeColorAndScope:
    def test_change_color_noop_without_selection(self, qapp):
        w, _ctx = _make_comp()
        w.change_color()  # no selection, no molecules -> currentRow() == -1
        w.destroy()

    def test_change_color_updates_entry_on_valid_pick(self, qapp, monkeypatch):
        from PyQt6.QtGui import QColor

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        w.mol_list.setCurrentRow(0)
        monkeypatch.setattr(
            _comparator.QColorDialog, "getColor", lambda *a, **k: QColor("#123456")
        )
        w.change_color()
        assert w.molecules[0]["color"] == "#123456"
        w.update_visualization.assert_called()
        w.destroy()

    def test_change_color_ignores_invalid_pick(self, qapp, monkeypatch):
        from PyQt6.QtGui import QColor

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        w.mol_list.setCurrentRow(0)
        original = w.molecules[0]["color"]
        monkeypatch.setattr(
            _comparator.QColorDialog, "getColor", lambda *a, **k: QColor()
        )
        w.change_color()
        assert w.molecules[0]["color"] == original
        w.destroy()

    def test_change_scope_noop_without_selection(self, qapp):
        w, _ctx = _make_comp()
        w.change_scope()  # no molecules -> currentRow() == -1
        w.destroy()

    def test_redraw_visualization_calls_both_updates(self, qapp):
        w, _ctx = _make_comp()
        w.redraw_visualization()
        w.update_visualization.assert_called_once()
        w.update_wireframe_lighting.assert_called_once()
        w.destroy()


# ===========================================================================
# load_from_file / add_molecule_from_path / drag & drop
# ===========================================================================


class TestLoadAndDrop:
    def test_load_from_file_cancel_is_noop(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog

        w, _ctx = _make_comp()
        monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: ("", ""))
        w.add_molecule_from_path = MagicMock()
        w.load_from_file()
        w.add_molecule_from_path.assert_not_called()
        w.destroy()

    def test_load_from_file_delegates_to_add(self, qapp, monkeypatch, tmp_path):
        from PyQt6.QtWidgets import QFileDialog

        w, _ctx = _make_comp()
        path = str(tmp_path / "mol.mol")
        monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (path, ""))
        w.add_molecule_from_path = MagicMock()
        w.load_from_file()
        w.add_molecule_from_path.assert_called_once_with(path)
        w.destroy()

    def test_add_molecule_from_path_mol_success(self, qapp, monkeypatch, tmp_path):
        w, _ctx = _make_comp()
        f = tmp_path / "a.mol"
        f.write_text("dummy")
        fake_mol = _WMol([_WAtom(0)])
        monkeypatch.setattr(_comparator.Chem, "MolFromMolFile", lambda *a, **k: fake_mol)
        monkeypatch.setattr(_comparator.Chem, "SanitizeMol", lambda m: None)
        w.add_molecule_from_path(str(f))
        assert len(w.molecules) == 1
        assert w.molecules[0]["name"] == "a.mol"
        w.update_visualization.assert_called_once()
        w.destroy()

    def test_add_molecule_from_path_pdb_success(self, qapp, monkeypatch, tmp_path):
        w, _ctx = _make_comp()
        f = tmp_path / "b.pdb"
        f.write_text("dummy")
        fake_mol = _WMol([_WAtom(0)])
        monkeypatch.setattr(_comparator.Chem, "MolFromPDBFile", lambda *a, **k: fake_mol)
        monkeypatch.setattr(_comparator.Chem, "SanitizeMol", lambda m: None)
        w.add_molecule_from_path(str(f))
        assert len(w.molecules) == 1
        w.destroy()

    def test_add_molecule_from_path_xyz_success_with_bond_perception(
        self, qapp, monkeypatch, tmp_path
    ):
        w, _ctx = _make_comp()
        f = tmp_path / "c.xyz"
        f.write_text("dummy")
        fake_mol = MagicMock()
        fake_mol.GetNumBonds.return_value = 0
        monkeypatch.setattr(_comparator.Chem, "MolFromXYZBlock", lambda s: fake_mol)
        monkeypatch.setattr(_comparator.Chem, "SanitizeMol", lambda m: None)
        determine = MagicMock()
        monkeypatch.setattr(_comparator, "rdDetermineBonds", determine)
        w.add_molecule_from_path(str(f))
        assert len(w.molecules) == 1
        determine.DetermineConnectivity.assert_called_once_with(fake_mol)
        determine.DetermineBondOrders.assert_called_once_with(fake_mol)
        w.destroy()

    def test_add_molecule_from_path_xyz_bond_perception_failure_silenced(
        self, qapp, monkeypatch, tmp_path
    ):
        w, _ctx = _make_comp()
        f = tmp_path / "d.xyz"
        f.write_text("dummy")
        fake_mol = MagicMock()
        fake_mol.GetNumBonds.return_value = 0
        monkeypatch.setattr(_comparator.Chem, "MolFromXYZBlock", lambda s: fake_mol)
        monkeypatch.setattr(_comparator.Chem, "SanitizeMol", lambda m: None)
        determine = MagicMock()
        determine.DetermineConnectivity.side_effect = ValueError("boom")
        monkeypatch.setattr(_comparator, "rdDetermineBonds", determine)
        w.add_molecule_from_path(str(f))
        assert len(w.molecules) == 1  # failure is swallowed, mol still added
        w.destroy()

    def test_add_molecule_from_path_load_failure_warns(self, qapp, monkeypatch, tmp_path):
        w, _ctx = _make_comp()
        f = tmp_path / "e.mol"
        f.write_text("dummy")
        monkeypatch.setattr(_comparator.Chem, "MolFromMolFile", lambda *a, **k: None)
        warned = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "warning", warned)
        w.add_molecule_from_path(str(f))
        warned.assert_called_once()
        assert w.molecules == []
        w.destroy()

    def test_add_molecule_from_path_sanitize_failure_silenced(
        self, qapp, monkeypatch, tmp_path
    ):
        w, _ctx = _make_comp()
        f = tmp_path / "f.mol"
        f.write_text("dummy")
        fake_mol = _WMol([_WAtom(0)])
        monkeypatch.setattr(_comparator.Chem, "MolFromMolFile", lambda *a, **k: fake_mol)

        def _sanitize_boom(m):
            raise ValueError("bad valence")

        monkeypatch.setattr(_comparator.Chem, "SanitizeMol", _sanitize_boom)
        w.add_molecule_from_path(str(f))
        assert len(w.molecules) == 1  # sanitize failure is swallowed
        w.destroy()

    def test_add_molecule_from_path_outer_exception_shows_critical(
        self, qapp, monkeypatch, tmp_path
    ):
        w, _ctx = _make_comp()
        f = tmp_path / "g.mol"
        f.write_text("dummy")

        def _boom(*a, **k):
            raise ValueError("catastrophic")

        monkeypatch.setattr(_comparator.Chem, "MolFromMolFile", _boom)
        critical = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "critical", critical)
        w.add_molecule_from_path(str(f))
        critical.assert_called_once()
        assert w.molecules == []
        w.destroy()

    def test_add_molecule_from_path_unsupported_extension_warns(
        self, qapp, monkeypatch, tmp_path
    ):
        w, _ctx = _make_comp()
        f = tmp_path / "h.foobar"
        f.write_text("dummy")
        warned = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "warning", warned)
        w.add_molecule_from_path(str(f))
        warned.assert_called_once()
        w.destroy()

    def test_drag_enter_accepts_urls(self, qapp):
        w, _ctx = _make_comp()
        event = MagicMock()
        event.mimeData.return_value.hasUrls.return_value = True
        w.dragEnterEvent(event)
        event.acceptProposedAction.assert_called_once()
        w.destroy()

    def test_drag_enter_ignores_without_urls(self, qapp):
        w, _ctx = _make_comp()
        event = MagicMock()
        event.mimeData.return_value.hasUrls.return_value = False
        w.dragEnterEvent(event)
        event.acceptProposedAction.assert_not_called()
        w.destroy()

    def test_drop_event_loads_local_files(self, qapp, monkeypatch, tmp_path):
        f = tmp_path / "dropped.mol"
        f.write_text("dummy")
        w, _ctx = _make_comp()
        w.add_molecule_from_path = MagicMock()
        url = MagicMock()
        url.toLocalFile.return_value = str(f)
        event = MagicMock()
        event.mimeData.return_value.urls.return_value = [url]
        w.dropEvent(event)
        w.add_molecule_from_path.assert_called_once_with(str(f))
        w.destroy()

    def test_drop_event_skips_missing_files(self, qapp):
        w, _ctx = _make_comp()
        w.add_molecule_from_path = MagicMock()
        url = MagicMock()
        url.toLocalFile.return_value = "/nonexistent/path/x.mol"
        event = MagicMock()
        event.mimeData.return_value.urls.return_value = [url]
        w.dropEvent(event)
        w.add_molecule_from_path.assert_not_called()
        w.destroy()


# ===========================================================================
# show_results_context_menu / export_molecule_coords
# ===========================================================================


class TestResultsContextMenu:
    def test_no_item_at_pos_is_noop(self, qapp):
        from PyQt6.QtCore import QPoint

        w, _ctx = _make_comp()
        assert w.table_results.itemAt(QPoint(5, 5)) is None
        w.show_results_context_menu(QPoint(5, 5))  # must not raise
        w.destroy()

    def test_menu_shown_for_valid_row(self, qapp, monkeypatch):
        from PyQt6.QtCore import QPoint

        w, _ctx = _make_comp()
        _add_mols(w, 2)
        w.update_results_table()
        item = MagicMock()
        item.row.return_value = 0
        monkeypatch.setattr(w.table_results, "itemAt", lambda pos: item)
        exec_mock = MagicMock()
        monkeypatch.setattr(_comparator.QMenu, "exec", exec_mock)
        w.show_results_context_menu(QPoint(1, 1))
        exec_mock.assert_called_once()
        w.destroy()

    def test_menu_out_of_range_row_is_noop(self, qapp, monkeypatch):
        from PyQt6.QtCore import QPoint

        w, _ctx = _make_comp()
        item = MagicMock()
        item.row.return_value = 5  # no molecules loaded
        monkeypatch.setattr(w.table_results, "itemAt", lambda pos: item)
        exec_mock = MagicMock()
        monkeypatch.setattr(_comparator.QMenu, "exec", exec_mock)
        w.show_results_context_menu(QPoint(1, 1))
        exec_mock.assert_not_called()
        w.destroy()

    def test_export_molecule_coords_out_of_range_is_noop(self, qapp):
        w, _ctx = _make_comp()
        w.export_molecule_coords(3, "mol")  # no molecules
        w.destroy()

    def test_export_molecule_coords_mol_writes_file(self, qapp, monkeypatch, tmp_path):
        from PyQt6.QtWidgets import QFileDialog

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        out = tmp_path / "out.mol"
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), "")
        )
        monkeypatch.setattr(
            _comparator.Chem, "MolToMolBlock", lambda m: "MOLBLOCK-DATA"
        )
        info = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "information", info)
        w.export_molecule_coords(0, "mol")
        assert "MOLBLOCK-DATA" in out.read_text()
        info.assert_called_once()
        w.destroy()

    def test_export_molecule_coords_xyz_writes_file(self, qapp, monkeypatch, tmp_path):
        from PyQt6.QtWidgets import QFileDialog

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        out = tmp_path / "out.xyz"
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), "")
        )
        monkeypatch.setattr(_comparator.Chem, "MolToXYZBlock", lambda m: "XYZ-DATA")
        info = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "information", info)
        w.export_molecule_coords(0, "xyz")
        assert "XYZ-DATA" in out.read_text()
        w.destroy()

    def test_export_molecule_coords_cancel_is_noop(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: ("", ""))
        w.export_molecule_coords(0, "mol")
        w.destroy()

    def test_export_molecule_coords_write_error_shows_critical(
        self, qapp, monkeypatch, tmp_path
    ):
        from PyQt6.QtWidgets import QFileDialog

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        bad_path = str(tmp_path / "nonexistent_dir" / "out.mol")
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: (bad_path, "")
        )
        monkeypatch.setattr(_comparator.Chem, "MolToMolBlock", lambda m: "X")
        critical = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "critical", critical)
        w.export_molecule_coords(0, "mol")
        critical.assert_called_once()
        w.destroy()


# ===========================================================================
# save_results_to_file / export_as_png
# ===========================================================================


class TestSaveResultsAndExportPng:
    def test_save_results_no_molecules_warns(self, qapp, monkeypatch):
        w, _ctx = _make_comp()
        warned = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "warning", warned)
        w.save_results_to_file()
        warned.assert_called_once()
        w.destroy()

    def test_save_results_cancel_is_noop(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QFileDialog

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: ("", ""))
        w.save_results_to_file()
        w.destroy()

    def test_save_results_writes_csv(self, qapp, monkeypatch, tmp_path):
        from PyQt6.QtWidgets import QFileDialog

        w, _ctx = _make_comp()
        _add_mols(w, 2)
        w.molecules[0]["rms"] = 0.5
        w.molecules[1]["name"] = "M1, extra"
        out = tmp_path / "results.csv"
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), "")
        )
        info = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "information", info)
        w.save_results_to_file()
        content = out.read_text(encoding="utf-8")
        assert "0.5000" in content
        assert '"M1, extra"' in content
        info.assert_called_once()
        w.destroy()

    def test_save_results_write_error_shows_critical(self, qapp, monkeypatch, tmp_path):
        from PyQt6.QtWidgets import QFileDialog

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        bad_path = str(tmp_path / "missing_dir" / "out.csv")
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: (bad_path, "")
        )
        critical = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "critical", critical)
        w.save_results_to_file()
        critical.assert_called_once()
        w.destroy()

    def test_export_png_no_molecules_warns(self, qapp, monkeypatch):
        w, _ctx = _make_comp()
        warned = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "warning", warned)
        w.export_as_png()
        warned.assert_called_once()
        w.destroy()

    def test_export_png_dialog_rejected_is_noop(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QDialog

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.DialogCode.Rejected)
        w.export_as_png()
        w.destroy()

    def test_export_png_cancel_after_settings_is_noop(self, qapp, monkeypatch):
        from PyQt6.QtWidgets import QDialog, QFileDialog

        w, _ctx = _make_comp()
        _add_mols(w, 1)
        monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.DialogCode.Accepted)
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: ("", ""))
        w.export_as_png()
        w.destroy()

    def test_export_png_success_appends_extension(self, qapp, monkeypatch, tmp_path):
        from PyQt6.QtWidgets import QDialog, QFileDialog

        w, ctx = _make_comp()
        _add_mols(w, 1)
        out = tmp_path / "shot"  # no .png extension
        monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.DialogCode.Accepted)
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), "")
        )
        info = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "information", info)
        w.export_as_png()
        ctx.plotter.screenshot.assert_called_once()
        args, kwargs = ctx.plotter.screenshot.call_args
        assert args[0].endswith(".png")
        info.assert_called_once()
        w.destroy()

    def test_export_png_screenshot_error_shows_critical(
        self, qapp, monkeypatch, tmp_path
    ):
        from PyQt6.QtWidgets import QDialog, QFileDialog

        w, ctx = _make_comp()
        _add_mols(w, 1)
        out = tmp_path / "shot.png"
        monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.DialogCode.Accepted)
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), "")
        )
        ctx.plotter.screenshot.side_effect = RuntimeError("no plotter")
        critical = MagicMock()
        monkeypatch.setattr(_comparator.QMessageBox, "critical", critical)
        w.export_as_png()
        critical.assert_called_once()
        w.destroy()


# ===========================================================================
# update_visualization — combined-molecule build + color application
# ===========================================================================


class _VAtom:
    def __init__(self, idx, z):
        self._idx = idx
        self._z = z

    def GetIdx(self):
        return self._idx

    def GetAtomicNum(self):
        return self._z


class _VMol:
    def __init__(self, atoms):
        self._atoms = list(atoms)

    def GetAtoms(self):
        return list(self._atoms)

    def GetNumAtoms(self):
        return len(self._atoms)


def _patch_chem_combine(monkeypatch):
    def fake_mol_ctor(*args):
        if not args:
            return _VMol([])
        return _VMol(args[0].GetAtoms())

    def fake_combine(a, b):
        return _VMol(a.GetAtoms() + b.GetAtoms())

    monkeypatch.setattr(_comparator.Chem, "Mol", fake_mol_ctor)
    monkeypatch.setattr(_comparator.Chem, "CombineMols", fake_combine)


class TestUpdateVisualization:
    def test_no_molecules_draws_none(self, qapp):
        ctx = _ctx_no_mol()
        w = _comparator.MoleculeComparator(context=ctx)
        ctx.draw_molecule_3d.reset_mock()
        w.update_visualization()
        ctx.draw_molecule_3d.assert_called_once_with(None)
        w.destroy()

    def test_controller_path_applies_colors_by_scope(self, qapp, monkeypatch):
        _patch_chem_combine(monkeypatch)
        w, ctx = _make_comp()
        w.update_visualization = _comparator.MoleculeComparator.update_visualization.__get__(w)
        w.molecules = [
            {
                "name": "A",
                "mol": _VMol([_VAtom(0, 6), _VAtom(1, 1)]),
                "color": "#ff0000",
                "scope": "Carbon Only",
                "rms": None,
            },
            {
                "name": "B",
                "mol": _VMol([_VAtom(0, 8)]),
                "color": "#00ff00",
                "scope": "All Atoms",
                "rms": None,
            },
        ]
        controller = MagicMock()
        ctx.get_3d_controller.return_value = controller
        w.update_visualization()
        # Carbon-only entry colors only the carbon (global idx 0); All-Atoms
        # entry colors its single oxygen atom (global idx 2).
        calls = [c.args for c in controller.set_atom_color.call_args_list]
        assert (0, "#ff0000") in calls
        assert (2, "#00ff00") in calls
        assert len(calls) == 2
        w.destroy()

    def test_controller_set_atom_color_exception_is_silenced(self, qapp, monkeypatch):
        _patch_chem_combine(monkeypatch)
        w, ctx = _make_comp()
        w.update_visualization = _comparator.MoleculeComparator.update_visualization.__get__(w)
        w.molecules = [
            {
                "name": "A",
                "mol": _VMol([_VAtom(0, 6)]),
                "color": "#ff0000",
                "scope": "All Atoms",
                "rms": None,
            }
        ]
        controller = MagicMock()
        controller.set_atom_color.side_effect = RuntimeError("no atom")
        ctx.get_3d_controller.return_value = controller
        w.update_visualization()  # must not raise
        w.destroy()

    def test_fallback_path_sets_color_overrides(self, qapp, monkeypatch):
        _patch_chem_combine(monkeypatch)
        w, ctx = _make_comp()
        w.update_visualization = _comparator.MoleculeComparator.update_visualization.__get__(w)
        w.mw = MagicMock()
        ctx.get_3d_controller.return_value = None
        w.molecules = [
            {
                "name": "A",
                "mol": _VMol([_VAtom(0, 6)]),
                "color": "#ff0000",
                "scope": "All Atoms",
                "rms": None,
            }
        ]
        w.update_visualization()
        assert w.mw.view_3d_manager._plugin_color_overrides == {0: "#ff0000"}
        w.destroy()


# ===========================================================================
# _find_style_tool_button / _find_style_actions / change_style — real
# QToolButton/QAction children so findChildren() works without mocking.
# ===========================================================================


class TestStyleHelpers:
    def _real_mw_with_actions(self):
        from PyQt6.QtWidgets import QWidget, QToolButton

        mw = QWidget()
        btn = QToolButton(mw)
        btn.setText("3D Style")
        actions = [
            _comparator.QAction("CPK", mw),
            _comparator.QAction("Sticks (Stick)", mw),
            _comparator.QAction("Ball && Stick", mw),
            _comparator.QAction("Wireframe", mw),
        ]
        for a in actions:
            a.setCheckable(True)  # setChecked() is a no-op otherwise
        return mw, btn, actions

    def test_find_style_tool_button_found(self, qapp):
        w, _ctx = _make_comp()
        mw, btn, _actions = self._real_mw_with_actions()
        w.mw = mw
        found = w._find_style_tool_button()
        assert found is btn
        w.destroy()

    def test_find_style_tool_button_missing_returns_none(self, qapp):
        from PyQt6.QtWidgets import QWidget

        w, _ctx = _make_comp()
        w.mw = QWidget()
        assert w._find_style_tool_button() is None
        w.destroy()

    def test_find_style_actions_matches_all_kinds(self, qapp):
        w, _ctx = _make_comp()
        mw, _btn, actions = self._real_mw_with_actions()
        w.mw = mw
        found = w._find_style_actions()
        assert len(found) == len(actions)
        w.destroy()

    def test_change_style_checks_matching_action_and_updates(self, qapp):
        w, _ctx = _make_comp()
        mw, _btn, actions = self._real_mw_with_actions()
        w.mw = mw
        w.mw.view_3d_manager = MagicMock()
        w.change_style("CPK")
        assert actions[0].isChecked()
        w.update_visualization.assert_called_once()
        w.mw.view_3d_manager.set_3d_style.assert_called_once_with("cpk")
        w.destroy()

    def test_change_style_ball_and_stick_matches(self, qapp):
        w, _ctx = _make_comp()
        mw, _btn, actions = self._real_mw_with_actions()
        w.mw = mw
        w.mw.view_3d_manager = MagicMock()
        w.change_style("Ball and Stick")
        assert actions[2].isChecked()
        w.destroy()

    def test_change_style_without_view_3d_manager(self, qapp):
        from PyQt6.QtWidgets import QWidget

        w, _ctx = _make_comp()
        w.mw = QWidget()  # no view_3d_manager attribute
        w.change_style("Wireframe")  # must not raise
        w.destroy()


# ===========================================================================
# update_wireframe_lighting / reset_view
# ===========================================================================


class _FakeProp:
    def __init__(self):
        self.lighting = None

    def SetLighting(self, v):
        self.lighting = v


class _FakeActor:
    def __init__(self):
        self._prop = _FakeProp()

    def GetProperty(self):
        return self._prop


class _FakeActorCollection:
    def __init__(self, actors):
        self._actors = list(actors)
        self._i = 0

    def __bool__(self):
        return bool(self._actors)

    def InitTraversal(self):
        self._i = 0

    def GetNumberOfItems(self):
        return len(self._actors)

    def GetNextItem(self):
        item = self._actors[self._i]
        self._i += 1
        return item


def _raw_comp():
    """MoleculeComparator with NO methods stubbed out (unlike _make_comp())."""
    ctx = _ctx_no_mol()
    w = _comparator.MoleculeComparator(context=ctx)
    return w, ctx


class TestUpdateWireframeLighting:
    def test_wireframe_style_uses_checkbox_state(self, qapp):
        w, _ctx = _raw_comp()
        w.combo_style.setCurrentText("Wireframe")
        w.check_wireframe_lighting.setChecked(True)
        actor = _FakeActor()
        w.context.plotter.renderer.GetActors.return_value = _FakeActorCollection(
            [actor]
        )
        w.context.plotter.render.reset_mock()
        w.update_wireframe_lighting()
        assert actor._prop.lighting is True
        w.context.plotter.render.assert_called_once()
        w.destroy()

    def test_non_wireframe_style_forces_lighting_on(self, qapp):
        w, _ctx = _raw_comp()
        w.combo_style.setCurrentText("Sticks")
        actor = _FakeActor()
        w.context.plotter.renderer.GetActors.return_value = _FakeActorCollection(
            [actor]
        )
        w.update_wireframe_lighting()
        assert actor._prop.lighting is True
        w.destroy()

    def test_no_actors_is_noop(self, qapp):
        w, _ctx = _raw_comp()
        w.context.plotter.renderer.GetActors.return_value = _FakeActorCollection([])
        w.update_wireframe_lighting()  # must not raise
        w.destroy()

    def test_plotter_exception_is_silenced(self, qapp):
        w, _ctx = _raw_comp()
        w.context.plotter.renderer.GetActors.side_effect = RuntimeError("no vtk")
        w.update_wireframe_lighting()  # must not raise
        w.destroy()


class TestResetView:
    def test_reset_view_calls_camera_reset(self, qapp, monkeypatch):
        w, ctx = _raw_comp()
        ctx.reset_3d_camera.reset_mock()
        ctx.plotter.render.reset_mock()
        monkeypatch.setattr(_comparator.QTimer, "singleShot", lambda ms, fn: fn())
        w.reset_view()
        ctx.reset_3d_camera.assert_called_once()
        ctx.plotter.render.assert_called_once()
        w.destroy()

    def test_reset_view_exception_is_silenced(self, qapp, monkeypatch):
        w, ctx = _raw_comp()
        ctx.reset_3d_camera.side_effect = RuntimeError("no camera")
        monkeypatch.setattr(_comparator.QTimer, "singleShot", lambda ms, fn: fn())
        w.reset_view()  # must not raise
        w.destroy()


# ===========================================================================
# enter_3d_only_mode / exit_3d_only_mode
# ===========================================================================


class _FakeSplitter:
    def __init__(self, sizes):
        self._sizes = sizes
        self.set_calls = []

    def sizes(self):
        return self._sizes

    def setSizes(self, s):
        self.set_calls.append(list(s))
        self._sizes = s


class TestEnter3DOnlyMode:
    def _mw_with_splitter(self, initial_sizes):
        mw = SimpleNamespace()
        mw.splitter = True  # only presence checked via hasattr
        splitter = _FakeSplitter(initial_sizes)
        edit_action = MagicMock()
        mw.init_manager = SimpleNamespace(splitter=splitter, edit_3d_action=edit_action)
        mw.ui_manager = SimpleNamespace(toggle_3d_edit_mode=MagicMock())
        return mw, splitter, edit_action

    def test_enter_collapses_visible_left_pane(self, qapp):
        w, _ctx = _make_comp()
        mw, splitter, edit_action = self._mw_with_splitter([200, 400])
        w.mw = mw
        w.enter_3d_only_mode()
        assert splitter.set_calls[-1] == [0, 10000]
        assert w.saved_splitter_sizes == [200, 400]
        mw.ui_manager.toggle_3d_edit_mode.assert_called_once_with(False)
        edit_action.setEnabled.assert_called_once_with(False)
        w.destroy()

    def test_enter_noop_when_already_collapsed(self, qapp):
        w, _ctx = _make_comp()
        mw, splitter, _edit_action = self._mw_with_splitter([0, 600])
        w.mw = mw
        w.enter_3d_only_mode()
        assert splitter.set_calls == []
        w.destroy()

    def test_enter_without_splitter_attribute(self, qapp):
        w, _ctx = _make_comp()
        # No `splitter` attr -> the collapse branch is skipped entirely, but
        # ui_manager/init_manager are touched unconditionally (no hasattr
        # guard on those two) so they must still be present.
        mw = SimpleNamespace(
            ui_manager=SimpleNamespace(), init_manager=SimpleNamespace()
        )
        w.mw = mw
        w.enter_3d_only_mode()  # must not raise
        w.destroy()

    def test_exit_restores_saved_splitter_sizes(self, qapp):
        w, _ctx = _make_comp()
        mw, splitter, edit_action = self._mw_with_splitter([200, 400])
        w.mw = mw
        w.enter_3d_only_mode()
        w.exit_3d_only_mode()
        assert splitter.set_calls[-1] == [200, 400]
        assert w.saved_splitter_sizes is None
        edit_action.setEnabled.assert_called_with(True)
        w.destroy()

    def test_exit_without_saved_sizes_is_noop(self, qapp):
        w, _ctx = _make_comp()
        mw, _splitter, edit_action = self._mw_with_splitter([200, 400])
        w.mw = mw
        w.exit_3d_only_mode()  # no enter() call first -> saved_splitter_sizes unset
        edit_action.setEnabled.assert_called_once_with(True)
        w.destroy()
