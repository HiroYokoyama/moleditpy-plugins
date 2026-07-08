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
        return SimpleNamespace(useSmallRingTorsions=None)

    def EmbedMultipleConfs(self, mol, numConfs=None, params=None):
        self.embed_calls.append((numConfs, params))
        return self.embed_return

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


def _run_search_fn(allchem, qmessage=None):
    return extract_function(
        CONF_SEARCH_PATH,
        "ConformerSearchDialog",
        "run_search",
        {
            "AllChem": allchem,
            "QMessageBox": qmessage if qmessage is not None else MagicMock(),
            "QApplication": SimpleNamespace(processEvents=lambda: None),
            "copy": __import__("copy"),
            "PLUGIN_NAME": "Conformational Search",
        },
    )


def _run_search_self(mol, ff="MMFF94"):
    return SimpleNamespace(
        context=SimpleNamespace(current_mol=mol),
        target_mol=mol,
        original_coords=[],
        btn_run=MagicMock(),
        lbl_info=MagicMock(),
        combo_ff=MagicMock(currentText=lambda: ff),
        results_raw=[],
        temp_mol=None,
        apply_filter_and_update=MagicMock(),
    )


class TestConfSearchRunSearch:
    def test_embed_called_with_expected_params(self):
        allchem = _FakeAllChem()
        fn = _run_search_fn(allchem)
        mol = _FakeConfSearchMol()
        self_ = _run_search_self(mol)
        fn(self_)
        assert len(allchem.embed_calls) == 1
        num_confs, params = allchem.embed_calls[0]
        assert num_confs == 30
        assert params.useSmallRingTorsions is True

    def test_mmff_energies_sorted_ascending(self):
        allchem = _FakeAllChem()
        fn = _run_search_fn(allchem)
        mol = _FakeConfSearchMol()
        self_ = _run_search_self(mol, ff="MMFF94")
        fn(self_)
        energies = [e for e, _ in self_.results_raw]
        assert energies == sorted(energies)
        assert energies == [3.0, 5.0, 8.0]
        self_.apply_filter_and_update.assert_called_once()
        assert self_.temp_mol is not None
        assert self_.temp_mol is not mol  # deep-copied, not the original

    def test_uff_branch_uses_uff_functions_only(self):
        allchem = _FakeAllChem()
        fn = _run_search_fn(allchem)
        mol = _FakeConfSearchMol()
        self_ = _run_search_self(mol, ff="UFF")
        fn(self_)
        assert allchem.uff_optimize_calls == [0, 1, 2]
        assert allchem.mmff_optimize_calls == []
        energies = [e for e, _ in self_.results_raw]
        assert energies == [0.5, 1.0, 2.0]

    def test_embed_failure_warns_and_stops(self):
        allchem = _FakeAllChem()
        allchem.embed_return = []
        fn = _run_search_fn(allchem)
        mol = _FakeConfSearchMol()
        self_ = _run_search_self(mol)
        fn(self_)
        self_.apply_filter_and_update.assert_not_called()
        assert self_.lbl_info.setText.call_args_list[-1][0][0] == "Failed."
        self_.btn_run.setEnabled.assert_any_call(True)

    def test_all_optimizations_fail_warns_and_stops(self):
        allchem = _FakeAllChem()
        allchem.mmff_optimize_return = -1
        fn = _run_search_fn(allchem)
        mol = _FakeConfSearchMol()
        self_ = _run_search_self(mol)
        fn(self_)
        self_.apply_filter_and_update.assert_not_called()
        assert self_.results_raw == []

    def test_exception_during_search_shows_critical(self):
        allchem = _FakeAllChem()
        allchem.etkdg_raises = True
        fn = _run_search_fn(allchem)
        mol = _FakeConfSearchMol()
        self_ = _run_search_self(mol)
        fn(self_)  # must not raise
        assert self_.lbl_info.setText.call_args_list[-1][0][0] == "Error occurred."
        self_.btn_run.setEnabled.assert_any_call(True)

    def test_target_mol_refreshed_when_current_mol_changed(self):
        allchem = _FakeAllChem()
        fn = _run_search_fn(allchem)
        old_mol = _FakeConfSearchMol()
        new_mol = _FakeConfSearchMol(n=3)
        self_ = _run_search_self(old_mol)
        self_.context.current_mol = new_mol
        fn(self_)
        assert self_.target_mol is new_mol


class TestConfSearchNoMoleculeWarning:
    def test_no_molecule_warns(self):
        qmsg = MagicMock()
        fn = _run_search_fn(_FakeAllChem(), qmessage=qmsg)
        self_ = _run_search_self(None)
        fn(self_)
        qmsg.warning.assert_called_once()


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
    return SimpleNamespace(context=ctx, target_mol=target_mol)


class TestConfSearchAcceptReject:
    def test_accept_pushes_undo_checkpoint_when_mol_present(self):
        fn = extract_function(CONF_SEARCH_PATH, "ConformerSearchDialog", "accept")
        ctx = make_context()
        self_ = _accept_self(_RSMol(1, [(0.0, 0.0, 0.0)]), ctx)
        fn(self_)
        ctx.push_undo_checkpoint.assert_called_once()
        ctx.register_window.assert_called_with("main_panel", None)

    def test_accept_no_undo_checkpoint_when_no_mol(self):
        fn = extract_function(CONF_SEARCH_PATH, "ConformerSearchDialog", "accept")
        ctx = make_context()
        self_ = _accept_self(None, ctx)
        fn(self_)
        ctx.push_undo_checkpoint.assert_not_called()
        ctx.register_window.assert_called_with("main_panel", None)

    def test_reject_restores_original_coordinates(self):
        fn = extract_function(CONF_SEARCH_PATH, "ConformerSearchDialog", "reject")
        ctx = make_context()
        mol = _RSMol(2, ["p0_new", "p1_new"])
        self_ = SimpleNamespace(
            context=ctx, target_mol=mol, original_coords=["p0", "p1"]
        )
        fn(self_)
        assert mol.set_calls == [(0, "p0"), (1, "p1")]
        ctx.refresh_3d_view.assert_called_once()
        ctx.register_window.assert_called_with("main_panel", None)

    def test_reject_without_mol_does_not_raise(self):
        fn = extract_function(CONF_SEARCH_PATH, "ConformerSearchDialog", "reject")
        ctx = make_context()
        self_ = SimpleNamespace(context=ctx, target_mol=None, original_coords=[])
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
            accept=MagicMock(),
            preview_conformer=MagicMock(),
        )
        fn(self_)
        self_.combo_ff.setCurrentText.assert_not_called()
