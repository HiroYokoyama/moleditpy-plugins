"""
Tests for the Conformational Search plugin.

Covers:
  1. initialize() must register at least one menu/export/plugin action
  2. No-molecule guard paths (run_plugin with mol=None)
  3. Dialog accept/reject round-trips
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
CONF_SEARCH_PATH = PLUGINS_DIR / "Conformational_Search" / "conf_search.py"


def _menu_registered(ctx: MagicMock) -> bool:
    """Return True if initialize() called any recognised registration method."""
    return (
        ctx.add_menu_action.called
        or ctx.add_export_action.called
        or ctx.add_plugin_menu.called
        or ctx.add_analysis_tool.called
        or ctx.add_toolbar_action.called
    )


class TestConformationalSearch:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert _menu_registered(ctx), (
                "Conformational Search initialize() must call add_menu_action()"
            )

    def test_initialize_menu_path_contains_conformational(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            call_args = ctx.add_menu_action.call_args
            assert call_args is not None
            path = call_args[0][0]
            assert "Conformational" in path or "conformational" in path.lower()

    def test_run_plugin_no_mol_warns(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            mock_warn.assert_called_once()

    def test_dialog_accept_does_not_raise(self):
        """ConformerSearchDialog.accept() must not raise when target_mol is None."""
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.current_mol = None
            # QDialog base is mocked — create a real instance normally
            dialog = mod.ConformerSearchDialog(ctx)
            dialog.accept()  # super().accept() calls mocked QDialog.accept → OK

    def test_run_plugin_registers_dialog_window(self):
        """run_plugin() when a mol is present creates and registers the conformer dialog."""
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            ctx.current_mol = MagicMock()
            mod.run_plugin(ctx)
            ctx.register_window.assert_called_once()


def _conf_filter_fn():
    return extract_function(
        CONF_SEARCH_PATH, "ConformerSearchDialog", "apply_filter_and_update", {}
    )


def _conf_self(results_raw, show_all=False):
    s = SimpleNamespace()
    s.results_raw = results_raw
    s.conformer_data = []
    s.cb_show_all = MagicMock()
    s.cb_show_all.isChecked.return_value = show_all
    s.update_table = MagicMock()
    s.lbl_info = MagicMock()
    return s


class TestConfSearchFilter:
    def test_empty_results_returns_early(self):
        fn = _conf_filter_fn()
        s = _conf_self([])
        fn(s)
        s.update_table.assert_not_called()
        assert s.conformer_data == []

    def test_show_all_keeps_everything(self):
        fn = _conf_filter_fn()
        raw = [(1.0, 0), (1.0, 1), (2.0, 2)]
        s = _conf_self(raw, show_all=True)
        fn(s)
        assert s.conformer_data == raw
        s.update_table.assert_called_once()

    def test_duplicate_energies_deduplicated(self):
        fn = _conf_filter_fn()
        # 1.0 and 1.00005 are within the 1e-4 window -> one survivor
        raw = [(1.0, 0), (1.00005, 1), (2.0, 2)]
        s = _conf_self(raw)
        fn(s)
        assert s.conformer_data == [(1.0, 0), (2.0, 2)]

    def test_distinct_energies_all_kept(self):
        fn = _conf_filter_fn()
        raw = [(1.0, 0), (1.5, 1), (2.0, 2)]
        s = _conf_self(raw)
        fn(s)
        assert s.conformer_data == raw

    def test_info_label_shows_counts(self):
        fn = _conf_filter_fn()
        raw = [(1.0, 0), (1.00001, 1), (3.0, 2)]
        s = _conf_self(raw)
        fn(s)
        msg = s.lbl_info.setText.call_args[0][0]
        assert "Showing 2 conformers" in msg
        assert "Total found: 3" in msg


class _FakeTable:
    """Recorder stand-in for the installer's QTableWidget."""

    def __init__(self):
        self.rows = 0
        self.items = {}  # (row, col) -> _Item
        self.cell_widgets = {}  # (row, col) -> widget
        self.hidden = {}  # row -> bool

    def setRowCount(self, n):
        self.rows = n
        if n == 0:
            self.items.clear()
            self.cell_widgets.clear()

    def rowCount(self):
        return self.rows

    def insertRow(self, row):
        self.rows += 1

    def setItem(self, row, col, item):
        self.items[(row, col)] = item

    def item(self, row, col):
        return self.items.get((row, col))

    def setCellWidget(self, row, col, widget):
        self.cell_widgets[(row, col)] = widget

    def setUpdatesEnabled(self, flag):
        pass

    def setRowHidden(self, row, hidden):
        self.hidden[row] = hidden


class TestConfSearchUpdateTable:
    def _run(self, conformer_data):
        items = []

        class _RecItem:
            def __init__(self, text):
                self.text = text
                self.user_data = None
                items.append(self)

            def setData(self, role, value):
                self.user_data = value

        globs = {"QTableWidgetItem": _RecItem, "Qt": MagicMock()}
        fn = extract_function(
            CONF_SEARCH_PATH, "ConformerSearchDialog", "update_table", globs
        )
        s = SimpleNamespace()
        s.conformer_data = conformer_data
        s.table = _FakeTable()
        fn(s)
        return s, items

    def test_rows_ranked_and_energy_formatted(self):
        s, _ = self._run([(1.23456789, 7), (2.5, 3)])
        assert s.table.rows == 2
        assert s.table.items[(0, 0)].text == "1"
        assert s.table.items[(0, 1)].text == "1.2346"  # 4 decimal places
        assert s.table.items[(1, 0)].text == "2"
        assert s.table.items[(1, 1)].text == "2.5000"

    def test_conformer_id_stored_as_user_data(self):
        s, _ = self._run([(1.0, 42)])
        assert s.table.items[(0, 0)].user_data == 42


# ---------------------------------------------------------------------------
# run_search: RDKit embedding orchestration
# ---------------------------------------------------------------------------


class _FakeFF:
    def __init__(self, energy):
        self._energy = energy
        self.initialized = 0

    def Initialize(self):
        self.initialized += 1

    def CalcEnergy(self):
        return self._energy


class _FakeAllChem:
    """Records calls and lets tests control embedding/optimization outcomes."""

    def __init__(self):
        self.embed_calls = []
        self.mmff_optimize_calls = []
        self.uff_optimize_calls = []
        self.embed_return = [0, 1, 2]
        self.mmff_optimize_return = 0  # 0 == success, -1 == failure
        self.uff_optimize_return = 0
        self.mmff_energies = {0: 5.0, 1: 3.0, 2: 8.0}
        self.uff_energies = {0: 1.0, 1: 2.0, 2: 0.5}
        self.etkdg_raises = False

    def ETKDGv3(self):
        if self.etkdg_raises:
            raise RuntimeError("boom")
        return SimpleNamespace(useSmallRingTorsions=None, clearConfs=True)

    def EmbedMultipleConfs(self, mol, numConfs=None, params=None):
        self.embed_calls.append((numConfs, params))
        # Embedding runs in batches; hand out the next slice each call.
        start = sum(n for n, _ in self.embed_calls[:-1])
        return self.embed_return[start : start + (numConfs or 0)]

    def MMFFOptimizeMolecule(self, mol, confId=None):
        self.mmff_optimize_calls.append(confId)
        return self.mmff_optimize_return

    def MMFFGetMoleculeProperties(self, mol):
        return object()

    def MMFFGetMoleculeForceField(self, mol, prop, confId=None):
        return _FakeFF(self.mmff_energies[confId])

    def UFFOptimizeMolecule(self, mol, confId=None):
        self.uff_optimize_calls.append(confId)
        return self.uff_optimize_return

    def UFFGetMoleculeForceField(self, mol, confId=None):
        return _FakeFF(self.uff_energies[confId])


class _FakeConfSearchMol:
    def __init__(self, n=2):
        self.n = n

    def GetNumAtoms(self):
        return self.n

    def GetConformer(self):
        return SimpleNamespace(GetAtomPosition=lambda i: (0.0, 0.0, 0.0))


class _SignalRecorder:
    """Stands in for a pyqtSignal: records every emit()."""

    def __init__(self):
        self.emissions = []

    def emit(self, *args):
        self.emissions.append(args)


def _worker_run_fn(allchem, chem=None):
    return extract_function(
        CONF_SEARCH_PATH,
        "ConformerSearchWorker",
        "run",
        {
            "Chem": chem
            or SimpleNamespace(AssignStereochemistryFrom3D=lambda m: None),
            "AllChem": allchem,
        },
    )


def _optimize_fn(allchem):
    return extract_function(
        CONF_SEARCH_PATH, "ConformerSearchWorker", "_optimize", {"AllChem": allchem}
    )


def _worker_self(allchem, mol=None, ff="MMFF94", num_confs=30, stop_after=None):
    """Fake worker instance.

    ``stop_after``: flip the stop flag once that many conformers have been
    optimized, standing in for the user hitting Stop mid-run.
    """
    optimize = _optimize_fn(allchem)
    self_ = SimpleNamespace(
        mol=mol if mol is not None else _FakeConfSearchMol(),
        num_confs=num_confs,
        force_field=ff,
        _stop=False,
        EMBED_BATCH=5,
        progress=_SignalRecorder(),
        failed=_SignalRecorder(),
        completed=_SignalRecorder(),
    )
    optimized = []

    def _optimize(mol_calc, cid):
        optimized.append(cid)
        if stop_after is not None and len(optimized) >= stop_after:
            self_._stop = True
        return optimize(self_, mol_calc, cid)

    self_._optimize = _optimize
    self_.optimized = optimized
    return self_


class TestConfSearchWorkerRun:
    def test_embed_called_in_batches_with_expected_params(self):
        allchem = _FakeAllChem()
        fn = _worker_run_fn(allchem)
        self_ = _worker_self(allchem)
        fn(self_)
        num_confs, params = allchem.embed_calls[0]
        assert num_confs == 5  # EMBED_BATCH, so Stop is honoured mid-embedding
        assert params.useSmallRingTorsions is True
        assert params.clearConfs is False  # batches accumulate

    def test_mmff_energies_sorted_ascending(self):
        allchem = _FakeAllChem()
        fn = _worker_run_fn(allchem)
        mol = _FakeConfSearchMol()
        self_ = _worker_self(allchem, mol=mol, ff="MMFF94")
        fn(self_)
        assert len(self_.completed.emissions) == 1
        emitted_mol, results, was_stopped = self_.completed.emissions[0]
        assert emitted_mol is mol
        assert [e for e, _ in results] == [3.0, 5.0, 8.0]
        assert was_stopped is False

    def test_uff_branch_uses_uff_functions_only(self):
        allchem = _FakeAllChem()
        fn = _worker_run_fn(allchem)
        self_ = _worker_self(allchem, ff="UFF")
        fn(self_)
        assert allchem.uff_optimize_calls == [0, 1, 2]
        assert allchem.mmff_optimize_calls == []
        _, results, _ = self_.completed.emissions[0]
        assert [e for e, _ in results] == [0.5, 1.0, 2.0]

    def test_embed_failure_reports_failure(self):
        allchem = _FakeAllChem()
        allchem.embed_return = []
        fn = _worker_run_fn(allchem)
        self_ = _worker_self(allchem)
        fn(self_)
        assert self_.completed.emissions == []
        assert self_.failed.emissions == [("Failed to generate conformers.",)]

    def test_all_optimizations_fail_reports_force_field(self):
        allchem = _FakeAllChem()
        allchem.mmff_optimize_return = -1
        fn = _worker_run_fn(allchem)
        self_ = _worker_self(allchem)
        fn(self_)
        assert self_.completed.emissions == []
        assert "MMFF94" in self_.failed.emissions[0][0]

    def test_exception_is_reported_not_raised(self):
        allchem = _FakeAllChem()
        allchem.etkdg_raises = True
        fn = _worker_run_fn(allchem)
        self_ = _worker_self(allchem)
        fn(self_)  # must not raise
        assert self_.failed.emissions[0][0].startswith("Error during search:")

    def test_progress_is_emitted_per_conformer(self):
        allchem = _FakeAllChem()
        fn = _worker_run_fn(allchem)
        self_ = _worker_self(allchem)
        fn(self_)
        optimizing = [e for e in self_.progress.emissions if "Optimizing" in e[2]]
        assert [done for done, _, _ in optimizing] == [1, 2, 3]
        assert all(total == 3 for _, total, _ in optimizing)


class TestConfSearchWorkerStop:
    def test_stop_before_embedding_completes_with_no_results(self):
        allchem = _FakeAllChem()
        fn = _worker_run_fn(allchem)
        self_ = _worker_self(allchem)
        self_._stop = True
        fn(self_)
        assert allchem.embed_calls == []
        assert self_.failed.emissions == []
        assert self_.completed.emissions == [(self_.mol, [], True)]

    def test_stop_mid_optimization_keeps_conformers_found_so_far(self):
        allchem = _FakeAllChem()
        fn = _worker_run_fn(allchem)
        self_ = _worker_self(allchem, stop_after=2)
        fn(self_)
        # Third conformer never optimized; the first two survive.
        assert self_.optimized == [0, 1]
        _, results, was_stopped = self_.completed.emissions[0]
        assert [cid for _, cid in results] == [1, 0]  # sorted by energy 3.0, 5.0
        assert was_stopped is True

    def test_stop_with_zero_optimized_reports_no_failure(self):
        allchem = _FakeAllChem()
        fn = _worker_run_fn(allchem)
        self_ = _worker_self(allchem, stop_after=1)
        allchem.mmff_optimize_return = -1  # the one attempt also fails
        fn(self_)
        assert self_.failed.emissions == []
        _, results, was_stopped = self_.completed.emissions[0]
        assert results == []
        assert was_stopped is True


# ---------------------------------------------------------------------------
# run_search: dialog-side guards and worker hand-off
# ---------------------------------------------------------------------------


def _dialog_run_search_fn(qmessage=None, worker_cls=None):
    return extract_function(
        CONF_SEARCH_PATH,
        "ConformerSearchDialog",
        "run_search",
        {
            "QMessageBox": qmessage if qmessage is not None else MagicMock(),
            "copy": __import__("copy"),
            "PLUGIN_NAME": "Conformational Search",
            "ConformerSearchWorker": worker_cls or MagicMock(),
        },
    )


def _dialog_run_search_self(mol, ff="MMFF94"):
    return SimpleNamespace(
        context=SimpleNamespace(current_mol=mol),
        target_mol=mol,
        original_coords=[],
        worker=None,
        lbl_info=MagicMock(),
        combo_ff=MagicMock(currentText=lambda: ff),
        spin_confs=MagicMock(value=lambda: 30),
        _set_running=MagicMock(),
        on_progress=MagicMock(),
        on_failed=MagicMock(),
        on_completed=MagicMock(),
    )


class TestConfSearchRunSearch:
    def test_worker_started_with_dialog_settings(self):
        worker_cls = MagicMock()
        fn = _dialog_run_search_fn(worker_cls=worker_cls)
        mol = _FakeConfSearchMol()
        self_ = _dialog_run_search_self(mol, ff="UFF")
        fn(self_)
        args, kwargs = worker_cls.call_args
        mol_calc, num_confs, force_field = args
        assert mol_calc is not mol  # deep-copied, not the molecule on screen
        assert num_confs == 30
        assert force_field == "UFF"
        assert kwargs["parent"] is self_
        self_.worker.start.assert_called_once()
        self_._set_running.assert_called_once_with(True)

    def test_second_click_while_running_is_ignored(self):
        worker_cls = MagicMock()
        fn = _dialog_run_search_fn(worker_cls=worker_cls)
        self_ = _dialog_run_search_self(_FakeConfSearchMol())
        self_.worker = MagicMock()
        fn(self_)
        worker_cls.assert_not_called()

    def test_target_mol_refreshed_when_current_mol_changed(self):
        fn = _dialog_run_search_fn()
        old_mol = _FakeConfSearchMol()
        new_mol = _FakeConfSearchMol(n=3)
        self_ = _dialog_run_search_self(old_mol)
        self_.context.current_mol = new_mol
        fn(self_)
        assert self_.target_mol is new_mol


class TestConfSearchNoMoleculeWarning:
    def test_no_molecule_warns(self):
        qmsg = MagicMock()
        worker_cls = MagicMock()
        fn = _dialog_run_search_fn(qmessage=qmsg, worker_cls=worker_cls)
        self_ = _dialog_run_search_self(None)
        fn(self_)
        qmsg.warning.assert_called_once()
        worker_cls.assert_not_called()


# ---------------------------------------------------------------------------
# stop_search / progress / completion handlers
# ---------------------------------------------------------------------------


def _dialog_method(name, extra_globals=None):
    return extract_function(
        CONF_SEARCH_PATH, "ConformerSearchDialog", name, extra_globals or {}
    )


class TestConfSearchStopHandlers:
    def test_stop_search_flags_worker_and_disables_button(self):
        fn = _dialog_method("stop_search")
        self_ = SimpleNamespace(
            worker=MagicMock(), btn_stop=MagicMock(), lbl_info=MagicMock()
        )
        fn(self_)
        self_.worker.stop.assert_called_once()
        self_.btn_stop.setEnabled.assert_called_once_with(False)

    def test_stop_search_without_worker_is_a_no_op(self):
        fn = _dialog_method("stop_search")
        self_ = SimpleNamespace(
            worker=None, btn_stop=MagicMock(), lbl_info=MagicMock()
        )
        fn(self_)  # must not raise
        self_.btn_stop.setEnabled.assert_not_called()

    def test_shutdown_worker_stops_and_waits(self):
        fn = _dialog_method("_shutdown_worker")
        worker = MagicMock()
        self_ = SimpleNamespace(worker=worker)
        fn(self_)
        worker.stop.assert_called_once()
        worker.wait.assert_called_once()
        assert self_.worker is None

    def test_on_failed_clears_worker_and_warns(self):
        qmsg = MagicMock()
        fn = _dialog_method("on_failed", {"QMessageBox": qmsg, "PLUGIN_NAME": "P"})
        self_ = SimpleNamespace(
            worker=MagicMock(), _set_running=MagicMock(), lbl_info=MagicMock()
        )
        fn(self_, "boom")
        assert self_.worker is None
        self_._set_running.assert_called_once_with(False)
        qmsg.warning.assert_called_once()

    def test_on_completed_publishes_results(self):
        fn = _dialog_method("on_completed")
        results = [(1.0, 0)]
        self_ = SimpleNamespace(
            worker=MagicMock(),
            _set_running=MagicMock(),
            lbl_info=MagicMock(),
            temp_mol=None,
            results_raw=[],
            apply_filter_and_update=MagicMock(),
        )
        mol = object()
        fn(self_, mol, results, False)
        assert self_.temp_mol is mol
        assert self_.results_raw == results
        self_.apply_filter_and_update.assert_called_once()
        assert self_.worker is None

    def test_on_completed_with_no_results_keeps_previous_preview(self):
        fn = _dialog_method("on_completed")
        self_ = SimpleNamespace(
            worker=MagicMock(),
            _set_running=MagicMock(),
            lbl_info=MagicMock(),
            temp_mol="previous",
            results_raw=[],
            apply_filter_and_update=MagicMock(),
        )
        fn(self_, object(), [], True)
        assert self_.temp_mol == "previous"
        self_.apply_filter_and_update.assert_not_called()

    def test_on_progress_scales_percentage(self):
        fn = _dialog_method("on_progress")
        self_ = SimpleNamespace(progress=MagicMock(), lbl_info=MagicMock())
        fn(self_, 3, 4, "Optimizing 3/4...")
        self_.progress.setValue.assert_called_once_with(75)
        self_.lbl_info.setText.assert_called_once_with("Optimizing 3/4...")

    def test_on_progress_ignores_zero_total(self):
        fn = _dialog_method("on_progress")
        self_ = SimpleNamespace(progress=MagicMock(), lbl_info=MagicMock())
        fn(self_, 0, 0, "Embedding...")
        self_.progress.setValue.assert_not_called()

    def test_set_running_toggles_controls(self):
        fn = _dialog_method("_set_running")
        self_ = SimpleNamespace(
            btn_run=MagicMock(),
            btn_stop=MagicMock(),
            combo_ff=MagicMock(),
            spin_confs=MagicMock(),
            progress=MagicMock(),
        )
        fn(self_, True)
        self_.btn_run.setEnabled.assert_called_once_with(False)
        self_.btn_stop.setEnabled.assert_called_once_with(True)
        self_.spin_confs.setEnabled.assert_called_once_with(False)
        self_.progress.setVisible.assert_called_once_with(True)


# ---------------------------------------------------------------------------
# accept / reject / closeEvent: real dialog instances
# ---------------------------------------------------------------------------


class _RSMol:
    """Fake rdkit mol usable both as mol and as its own conformer."""

    def __init__(self, n, positions):
        self.n = n
        self.positions = list(positions)
        self.set_calls = []

    def GetNumAtoms(self):
        return self.n

    def GetConformer(self):
        return self

    def GetAtomPosition(self, i):
        return self.positions[i]

    def SetAtomPosition(self, i, pos):
        self.set_calls.append((i, pos))


def _accept_self(target_mol, ctx):
    return SimpleNamespace(
        context=ctx, target_mol=target_mol, _shutdown_worker=lambda: None
    )


# accept()/reject() call super().accept()/reject(); zero-arg super() needs a
# __class__ cell that only exists inside a real class body, so we stub the
# global `super` builtin when extracting these methods standalone.
_FAKE_SUPER_GLOBALS = {
    "super": lambda *a: SimpleNamespace(accept=lambda: None, reject=lambda: None)
}


class TestConfSearchAcceptReject:
    def test_accept_pushes_undo_checkpoint_when_mol_present(self):
        fn = extract_function(
            CONF_SEARCH_PATH, "ConformerSearchDialog", "accept", _FAKE_SUPER_GLOBALS
        )
        ctx = make_context()
        self_ = _accept_self(_RSMol(1, [(0.0, 0.0, 0.0)]), ctx)
        fn(self_)
        ctx.push_undo_checkpoint.assert_called_once()
        ctx.register_window.assert_called_with("main_panel", None)

    def test_accept_no_undo_checkpoint_when_no_mol(self):
        fn = extract_function(
            CONF_SEARCH_PATH, "ConformerSearchDialog", "accept", _FAKE_SUPER_GLOBALS
        )
        ctx = make_context()
        self_ = _accept_self(None, ctx)
        fn(self_)
        ctx.push_undo_checkpoint.assert_not_called()
        ctx.register_window.assert_called_with("main_panel", None)

    def test_reject_restores_original_coordinates(self):
        fn = extract_function(
            CONF_SEARCH_PATH, "ConformerSearchDialog", "reject", _FAKE_SUPER_GLOBALS
        )
        ctx = make_context()
        mol = _RSMol(2, ["p0_new", "p1_new"])
        self_ = SimpleNamespace(
            context=ctx,
            target_mol=mol,
            original_coords=["p0", "p1"],
            _shutdown_worker=lambda: None,
        )
        fn(self_)
        assert mol.set_calls == [(0, "p0"), (1, "p1")]
        ctx.refresh_3d_view.assert_called_once()
        ctx.register_window.assert_called_with("main_panel", None)

    def test_reject_without_mol_does_not_raise(self):
        fn = extract_function(
            CONF_SEARCH_PATH, "ConformerSearchDialog", "reject", _FAKE_SUPER_GLOBALS
        )
        ctx = make_context()
        self_ = SimpleNamespace(
            context=ctx,
            target_mol=None,
            original_coords=[],
            _shutdown_worker=lambda: None,
        )
        fn(self_)  # must not raise
        ctx.register_window.assert_called_with("main_panel", None)

    def test_close_event_delegates_to_accept_and_ignores_event(self):
        fn = extract_function(CONF_SEARCH_PATH, "ConformerSearchDialog", "closeEvent")
        ctx = make_context()
        calls = []
        self_ = SimpleNamespace(accept=lambda: calls.append("accept"))
        event = MagicMock()
        fn(self_, event)
        assert calls == ["accept"]
        event.ignore.assert_called_once()


# ---------------------------------------------------------------------------
# preview_conformer
# ---------------------------------------------------------------------------


def _preview_fn():
    qt = SimpleNamespace(ItemDataRole=SimpleNamespace(UserRole="UR"))
    return extract_function(
        CONF_SEARCH_PATH, "ConformerSearchDialog", "preview_conformer", {"Qt": qt}
    )


class _PreviewItem:
    def __init__(self, cid):
        self._cid = cid

    def data(self, role):
        return self._cid


class _PreviewTable:
    def __init__(self, cid):
        self._item = _PreviewItem(cid)

    def item(self, row, col):
        return self._item


class TestConfSearchPreviewConformer:
    def test_matching_atom_counts_copies_coordinates(self):
        fn = _preview_fn()
        set_calls = []
        temp_mol = SimpleNamespace(
            GetNumAtoms=lambda: 2,
            GetConformer=lambda cid: SimpleNamespace(
                GetAtomPosition=lambda i: f"p{i}"
            ),
        )
        target_mol = SimpleNamespace(
            GetNumAtoms=lambda: 2,
            GetConformer=lambda: SimpleNamespace(
                SetAtomPosition=lambda i, pos: set_calls.append((i, pos))
            ),
        )
        self_ = SimpleNamespace(
            temp_mol=temp_mol,
            target_mol=target_mol,
            table=_PreviewTable(5),
            lbl_info=MagicMock(),
            context=MagicMock(),
        )
        current = SimpleNamespace(row=lambda: 0)
        fn(self_, current, None)
        assert set_calls == [(0, "p0"), (1, "p1")]
        assert self_.context.current_mol is target_mol
        self_.context.refresh_3d_view.assert_called_once()

    def test_atom_count_mismatch_shows_error_and_skips_copy(self):
        fn = _preview_fn()
        temp_mol = SimpleNamespace(
            GetNumAtoms=lambda: 2,
            GetConformer=lambda cid: SimpleNamespace(GetAtomPosition=lambda i: f"p{i}"),
        )
        target_mol = SimpleNamespace(
            GetNumAtoms=lambda: 3,
            GetConformer=lambda: SimpleNamespace(SetAtomPosition=MagicMock()),
        )
        self_ = SimpleNamespace(
            temp_mol=temp_mol,
            target_mol=target_mol,
            table=_PreviewTable(5),
            lbl_info=MagicMock(),
            context=MagicMock(),
        )
        current = SimpleNamespace(row=lambda: 0)
        fn(self_, current, None)
        msg = self_.lbl_info.setText.call_args[0][0]
        assert "Restart search" in msg
        self_.context.refresh_3d_view.assert_not_called()

    def test_no_current_selection_returns_early(self):
        fn = _preview_fn()
        self_ = SimpleNamespace(
            temp_mol=SimpleNamespace(),
            target_mol=SimpleNamespace(),
            table=None,
            lbl_info=MagicMock(),
            context=MagicMock(),
        )
        fn(self_, None, None)  # current=None -> early return, no crash
        self_.context.refresh_3d_view.assert_not_called()

    def test_no_temp_mol_returns_early(self):
        fn = _preview_fn()
        self_ = SimpleNamespace(
            temp_mol=None,
            target_mol=SimpleNamespace(),
            table=None,
            lbl_info=MagicMock(),
            context=MagicMock(),
        )
        current = SimpleNamespace(row=lambda: 0)
        fn(self_, current, None)  # must not raise
        self_.context.refresh_3d_view.assert_not_called()


# ---------------------------------------------------------------------------
# run_plugin: reuse existing window
# ---------------------------------------------------------------------------


class TestConfSearchRunPluginReuse:
    def test_reuses_existing_window_without_creating_new_dialog(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            ctx.current_mol = MagicMock()
            win = MagicMock()
            ctx.get_window.return_value = win
            mod.run_plugin(ctx)
            win.show.assert_called_once()
            win.raise_.assert_called_once()
            win.activateWindow.assert_called_once()
            ctx.register_window.assert_not_called()


# ---------------------------------------------------------------------------
# legacy run(mw) entry point
# ---------------------------------------------------------------------------


class TestConfSearchLegacyRun:
    def test_run_invokes_launch_fn_when_set(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            mod._launch_fn = MagicMock()
            mod.run(SimpleNamespace())
            mod._launch_fn.assert_called_once()

    def test_run_noop_when_launch_fn_none(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            mod._launch_fn = None
            mod.run(SimpleNamespace())  # must not raise

    def test_run_unwraps_host_attribute(self):
        with mock_optional_imports():
            mod = load_plugin(CONF_SEARCH_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            mod._launch_fn = MagicMock()
            mw = SimpleNamespace(host=SimpleNamespace())
            mod.run(mw)
            mod._launch_fn.assert_called_once()


# ---------------------------------------------------------------------------
# dialog default force-field selection from context.get_setting
# ---------------------------------------------------------------------------


def _init_ui_fn():
    return extract_function(
        CONF_SEARCH_PATH,
        "ConformerSearchDialog",
        "init_ui",
        {
            "QVBoxLayout": MagicMock(),
            "QHBoxLayout": MagicMock(),
            "QLabel": MagicMock(),
            "QComboBox": MagicMock(),
            "QCheckBox": MagicMock(),
            "QTableWidget": MagicMock(),
            "QHeaderView": SimpleNamespace(ResizeMode=SimpleNamespace(Stretch=1)),
            "QAbstractItemView": SimpleNamespace(
                SelectionBehavior=SimpleNamespace(SelectRows=1),
                SelectionMode=SimpleNamespace(SingleSelection=1),
                EditTrigger=SimpleNamespace(NoEditTriggers=1),
            ),
            "QPushButton": MagicMock(),
            "QSpinBox": MagicMock(),
            "QProgressBar": MagicMock(),
        },
    )


class TestConfSearchDefaultForceField:
    def test_uff_default_setting_sets_combo_to_uff(self):
        fn = _init_ui_fn()
        self_ = SimpleNamespace(
            context=SimpleNamespace(get_setting=MagicMock(return_value="UFF_RDKIT")),
            combo_ff=MagicMock(),
            apply_filter_and_update=MagicMock(),
            run_search=MagicMock(),
            stop_search=MagicMock(),
            accept=MagicMock(),
            preview_conformer=MagicMock(),
        )
        fn(self_)
        self_.combo_ff.setCurrentText.assert_called_with("UFF")

    def test_mmff_default_setting_sets_combo_to_mmff94(self):
        fn = _init_ui_fn()
        self_ = SimpleNamespace(
            context=SimpleNamespace(get_setting=MagicMock(return_value="MMFF_RDKIT")),
            combo_ff=MagicMock(),
            apply_filter_and_update=MagicMock(),
            run_search=MagicMock(),
            stop_search=MagicMock(),
            accept=MagicMock(),
            preview_conformer=MagicMock(),
        )
        fn(self_)
        self_.combo_ff.setCurrentText.assert_called_with("MMFF94")

    def test_no_default_setting_skips_combo_override(self):
        fn = _init_ui_fn()
        self_ = SimpleNamespace(
            context=SimpleNamespace(get_setting=MagicMock(return_value=None)),
            combo_ff=MagicMock(),
            apply_filter_and_update=MagicMock(),
            run_search=MagicMock(),
            stop_search=MagicMock(),
            accept=MagicMock(),
            preview_conformer=MagicMock(),
        )
        fn(self_)
        self_.combo_ff.setCurrentText.assert_not_called()


class TestStereochemistryPreservation:
    """ETKDG embeds from the graph, not from the displayed coordinates.

    Without re-perceiving stereo from 3D first, a molecule with unset chiral
    tags produces a mix of both enantiomers, and one whose tags went stale
    after a 3D edit produces conformers that are all the mirror image of what
    the user is looking at.
    """

    def _drive(self, assign):
        calls = []
        mol = MagicMock(name="mol")
        chem = SimpleNamespace(AssignStereochemistryFrom3D=lambda m: assign(calls, m))

        def embed(m, numConfs=0, params=None):
            calls.append(("embed", m))
            return []  # no conformers -> the worker stops after one batch

        allchem = SimpleNamespace(
            ETKDGv3=lambda: SimpleNamespace(
                useSmallRingTorsions=False, clearConfs=True
            ),
            EmbedMultipleConfs=embed,
        )
        fn = _worker_run_fn(allchem, chem=chem)
        self_ = _worker_self(MagicMock(), mol=mol)
        fn(self_)
        return calls, mol

    @staticmethod
    def _record(calls, m):
        calls.append(("assign_from_3d", m))

    def test_stereo_is_reperceived_from_3d_before_embedding(self):
        calls, mol = self._drive(self._record)
        names = [name for name, _ in calls]
        assert "assign_from_3d" in names, (
            "the worker must call Chem.AssignStereochemistryFrom3D; without it "
            "ETKDG follows stale/absent chiral tags and can invert the molecule"
        )
        assert names.index("assign_from_3d") < names.index("embed"), (
            "stereo must be re-perceived BEFORE EmbedMultipleConfs, not after"
        )

    def test_stereo_is_reperceived_on_the_calculation_copy(self):
        """It must run on the working copy, never mutate the user's molecule."""
        calls, mol = self._drive(self._record)
        assigned = [m for name, m in calls if name == "assign_from_3d"]
        embedded = [m for name, m in calls if name == "embed"]
        assert assigned == embedded, "the embedded mol must be the one re-perceived"

    def test_stereo_failure_does_not_abort_the_search(self):
        """A molecule with no conformer must not kill the whole search."""

        def boom(calls, m):
            calls.append(("assign_from_3d", m))
            raise ValueError("no conformer")

        calls, mol = self._drive(boom)
        assert [n for n, _ in calls] == ["assign_from_3d", "embed"]
