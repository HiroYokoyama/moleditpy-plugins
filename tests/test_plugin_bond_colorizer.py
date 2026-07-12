"""
Tests for the Bond Colorizer plugin: save/load/reset handler payloads,
apply_color logic (bond-index input AND atom-pair input).

BondColorizerWindow inherits a Qt base, so apply_color/get_selection_from_viewer
are extracted from the plugin source via AST and invoked with a fake self object.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
COLORIZER_PATH = PLUGINS_DIR / "Bond_Colorizer" / "bond_colorizer.py"


class TestBondColorizerHandlers:
    def _handlers(self, ctx):
        with mock_optional_imports():
            mod = load_plugin(COLORIZER_PATH)
            mod.initialize(ctx)
        save = ctx.register_save_handler.call_args[0][0]
        load = ctx.register_load_handler.call_args[0][0]
        reset = ctx.register_document_reset_handler.call_args[0][0]
        return save, load, reset

    def test_initialize_registers_all_three_handlers(self):
        ctx = make_context()
        self._handlers(ctx)
        ctx.register_save_handler.assert_called_once()
        ctx.register_load_handler.assert_called_once()
        ctx.register_document_reset_handler.assert_called_once()

    def test_save_returns_empty_dict_without_manager(self):
        ctx = make_context()
        save, _, _ = self._handlers(ctx)
        ctx.get_main_window.return_value.view_3d_manager = None
        assert save() == {}

    def test_save_stringifies_override_keys(self):
        ctx = make_context()
        save, _, _ = self._handlers(ctx)
        mw = ctx.get_main_window.return_value
        mw.view_3d_manager._plugin_bond_color_overrides = {0: "#ff0000", 2: "#00ff00"}
        assert save() == {"bond_colors": {"0": "#ff0000", "2": "#00ff00"}}

    def test_load_with_empty_data_is_noop(self):
        ctx = make_context()
        _, load, _ = self._handlers(ctx)
        load({})
        ctx.get_3d_controller.assert_not_called()
        ctx.refresh_3d_view.assert_not_called()

    def test_load_applies_each_color(self):
        ctx = make_context()
        _, load, _ = self._handlers(ctx)
        load({"bond_colors": {"1": "#00ff00", "3": "#0000ff"}})
        controller = ctx.get_3d_controller.return_value
        controller.set_bond_color.assert_any_call(1, "#00ff00")
        controller.set_bond_color.assert_any_call(3, "#0000ff")
        ctx.refresh_3d_view.assert_called_once()

    def test_load_with_bad_index_does_not_raise(self):
        ctx = make_context()
        _, load, _ = self._handlers(ctx)
        load({"bond_colors": {"abc": "#00ff00"}})
        ctx.refresh_3d_view.assert_called_once()

    def test_reset_clears_overrides_and_both_fields(self):
        ctx = make_context()
        _, _, reset = self._handlers(ctx)
        mw = ctx.get_main_window.return_value
        overrides = {0: "#ff0000"}
        mw.view_3d_manager._plugin_bond_color_overrides = overrides
        reset()
        assert overrides == {}
        ctx.get_window.return_value.le_bond_ids.clear.assert_called_once()
        ctx.get_window.return_value.le_atom_pairs.clear.assert_called_once()

    def test_run_without_context_is_noop(self):
        with mock_optional_imports():
            mod = load_plugin(COLORIZER_PATH)
        assert mod.PLUGIN_CONTEXT is None
        mod.run(MagicMock())


class TestBondColorizerApplyColor:
    def _apply(self, qmb=None):
        return extract_function(
            COLORIZER_PATH, "BondColorizerWindow", "apply_color",
            extra_globals={"QMessageBox": qmb or MagicMock()},
        )

    def _self(self, bond_ids_text="", atom_pairs_text="", num_bonds=5, bond_between=None):
        s = MagicMock()
        s.le_bond_ids.text.return_value = bond_ids_text
        s.le_atom_pairs.text.return_value = atom_pairs_text
        s.current_color.name.return_value = "#112233"
        s.context.current_molecule.GetNumBonds.return_value = num_bonds

        def get_bond_between(a1, a2):
            if bond_between is None:
                return None
            return bond_between.get((a1, a2))

        s.context.current_molecule.GetBondBetweenAtoms.side_effect = get_bond_between
        return s

    def test_both_fields_empty_warns(self):
        qmb = MagicMock()
        fn = self._apply(qmb)
        s = self._self("", "")
        fn(s)
        qmb.warning.assert_called_once()
        s.context.get_3d_controller.assert_not_called()

    def test_bond_ids_applied(self):
        fn = self._apply()
        s = self._self(bond_ids_text="0,2")
        fn(s)
        controller = s.context.get_3d_controller.return_value
        controller.set_bond_color.assert_any_call(0, "#112233")
        controller.set_bond_color.assert_any_call(2, "#112233")
        assert controller.set_bond_color.call_count == 2
        s.context.refresh_3d_view.assert_called_once()

    def test_atom_pair_resolved_to_bond_idx(self):
        fake_bond = MagicMock()
        fake_bond.GetIdx.return_value = 4
        fn = self._apply()
        s = self._self(atom_pairs_text="0-1", bond_between={(0, 1): fake_bond})
        fn(s)
        controller = s.context.get_3d_controller.return_value
        controller.set_bond_color.assert_called_once_with(4, "#112233")

    def test_non_bonded_pair_skipped_and_warned_valid_still_applied(self):
        qmb = MagicMock()
        fake_bond = MagicMock()
        fake_bond.GetIdx.return_value = 4
        fn = self._apply(qmb)
        s = self._self(
            atom_pairs_text="0-1,2-3", bond_between={(0, 1): fake_bond}
        )
        fn(s)
        controller = s.context.get_3d_controller.return_value
        controller.set_bond_color.assert_called_once_with(4, "#112233")
        qmb.warning.assert_called_once()

    def test_malformed_bond_id_warns_and_returns(self):
        qmb = MagicMock()
        fn = self._apply(qmb)
        s = self._self(bond_ids_text="1, x, 2")
        fn(s)
        qmb.warning.assert_called_once()
        s.context.get_3d_controller.assert_not_called()

    def test_malformed_atom_pair_warns_and_returns(self):
        qmb = MagicMock()
        fn = self._apply(qmb)
        s = self._self(atom_pairs_text="0-1-2")
        fn(s)
        qmb.warning.assert_called_once()
        s.context.get_3d_controller.assert_not_called()

    def test_malformed_atom_pair_non_int_warns_and_returns(self):
        qmb = MagicMock()
        fn = self._apply(qmb)
        s = self._self(atom_pairs_text="a-b")
        fn(s)
        qmb.warning.assert_called_once()
        s.context.get_3d_controller.assert_not_called()

    def test_no_molecule_warns(self):
        qmb = MagicMock()
        fn = self._apply(qmb)
        s = self._self(bond_ids_text="0")
        s.context.current_molecule = None
        fn(s)
        qmb.warning.assert_called_once()

    def test_out_of_range_bond_id_skipped_and_warned(self):
        qmb = MagicMock()
        fn = self._apply(qmb)
        s = self._self(bond_ids_text="1, 50", num_bonds=3)
        fn(s)
        controller = s.context.get_3d_controller.return_value
        controller.set_bond_color.assert_called_once_with(1, "#112233")
        qmb.warning.assert_called_once()


class TestBondColorizerResetColors:
    def _fn(self):
        return extract_function(
            COLORIZER_PATH, "BondColorizerWindow", "reset_colors",
            extra_globals={"QMessageBox": MagicMock()},
        )

    def test_resets_each_bond_and_clears_fields(self):
        fn = self._fn()
        s = MagicMock()
        s.context.current_molecule.GetNumBonds.return_value = 3
        fn(s)
        controller = s.context.get_3d_controller.return_value
        controller.set_bond_color.assert_any_call(0, None)
        controller.set_bond_color.assert_any_call(1, None)
        controller.set_bond_color.assert_any_call(2, None)
        assert controller.set_bond_color.call_count == 3
        s.context.refresh_3d_view.assert_called_once()
        s.le_bond_ids.clear.assert_called_once()
        s.le_atom_pairs.clear.assert_called_once()

    def test_no_molecule_is_noop(self):
        fn = self._fn()
        s = MagicMock()
        s.context.current_molecule = None
        fn(s)
        s.context.get_3d_controller.assert_not_called()


class _FakeLineEdit:
    def __init__(self, text="", focused=False):
        self._text = text
        self._focused = focused
        self.set_calls = []

    def text(self):
        return self._text

    def setText(self, value):
        self.set_calls.append(value)
        self._text = value

    def hasFocus(self):
        return self._focused


class _FakeBond:
    def __init__(self, a1, a2, idx):
        self._a1 = a1
        self._a2 = a2
        self._idx = idx

    def GetBeginAtomIdx(self):
        return self._a1

    def GetEndAtomIdx(self):
        return self._a2

    def GetIdx(self):
        return self._idx


class TestGetSelectionFromViewer:
    """get_selection_from_viewer merges 2D + 3D + measurement atom picks,
    then derives the bonds fully inside that selection into le_atom_pairs."""

    def _fn(self):
        return extract_function(
            COLORIZER_PATH, "BondColorizerWindow", "get_selection_from_viewer"
        )

    def _self(self, sel_2d, sel_3d=None, picked=None, field_text="", bonds=None, has_mol=True):
        fake = MagicMock()
        fake.context.get_selected_atom_indices.return_value = sel_2d
        mw = MagicMock()
        mw.edit_3d_manager.selected_atoms_3d = sel_3d if sel_3d is not None else set()
        mw.edit_3d_manager.selected_atoms_for_measurement = picked if picked is not None else []
        fake.context.get_main_window.return_value = mw
        fake.le_atom_pairs = _FakeLineEdit(field_text)
        if has_mol:
            fake.context.current_molecule.GetBonds.return_value = bonds or []
        else:
            fake.context.current_molecule = None
        return fake

    def test_derives_bonds_fully_inside_selection(self):
        bonds = [
            _FakeBond(0, 1, 10),
            _FakeBond(1, 2, 11),
            _FakeBond(2, 5, 12),  # 5 not selected
        ]
        fake = self._self([0, 1, 2], bonds=bonds)
        self._fn()(fake)
        assert fake.le_atom_pairs.text() == "0-1,1-2"

    def test_merges_2d_3d_and_measurement_sources(self):
        bonds = [_FakeBond(2, 5, 1), _FakeBond(5, 7, 2)]
        fake = self._self([0, 2], sel_3d={5}, picked=(7,), bonds=bonds)
        self._fn()(fake)
        assert fake.le_atom_pairs.text() == "2-5,5-7"

    def test_no_settext_when_selection_unchanged(self):
        bonds = [_FakeBond(1, 3, 0)]
        fake = self._self([1, 3], field_text="1-3", bonds=bonds)
        self._fn()(fake)
        assert fake.le_atom_pairs.set_calls == []

    def test_no_molecule_writes_nothing(self):
        fake = self._self([0, 1], has_mol=False)
        self._fn()(fake)
        assert fake.le_atom_pairs.set_calls == []

    def test_no_main_window_uses_2d_only(self):
        bonds = [_FakeBond(1, 9, 0)]
        fake = self._self([9, 1], bonds=bonds)
        fake.context.get_main_window.return_value = None
        self._fn()(fake)
        assert fake.le_atom_pairs.text() == "1-9"


class TestAutoUpdateSelection:
    def _fn(self):
        return extract_function(
            COLORIZER_PATH, "BondColorizerWindow", "_auto_update_selection"
        )

    def test_skips_update_while_bond_ids_focused(self):
        fake = MagicMock()
        fake.le_bond_ids.hasFocus.return_value = True
        fake.le_atom_pairs.hasFocus.return_value = False
        self._fn()(fake)
        fake.get_selection_from_viewer.assert_not_called()

    def test_skips_update_while_atom_pairs_focused(self):
        fake = MagicMock()
        fake.le_bond_ids.hasFocus.return_value = False
        fake.le_atom_pairs.hasFocus.return_value = True
        self._fn()(fake)
        fake.get_selection_from_viewer.assert_not_called()

    def test_updates_when_neither_focused(self):
        fake = MagicMock()
        fake.le_bond_ids.hasFocus.return_value = False
        fake.le_atom_pairs.hasFocus.return_value = False
        self._fn()(fake)
        fake.get_selection_from_viewer.assert_called_once()


class TestBondColorizerLoadRefresh:
    def _handlers(self, ctx):
        with mock_optional_imports():
            mod = load_plugin(COLORIZER_PATH)
            mod.initialize(ctx)
        return (
            ctx.register_save_handler.call_args[0][0],
            ctx.register_load_handler.call_args[0][0],
            ctx.register_document_reset_handler.call_args[0][0],
        )

    def test_load_refreshes_view_even_with_empty_color_map(self):
        ctx = make_context()
        _, load, _ = self._handlers(ctx)
        load({"bond_colors": {}})
        ctx.refresh_3d_view.assert_called_once()

    def test_load_none_data_does_not_refresh(self):
        ctx = make_context()
        _, load, _ = self._handlers(ctx)
        load(None)
        ctx.refresh_3d_view.assert_not_called()

    def test_reset_without_window_does_not_raise(self):
        ctx = make_context()
        _, _, reset = self._handlers(ctx)
        ctx.get_window.return_value = None
        mw = ctx.get_main_window.return_value
        mw.view_3d_manager._plugin_bond_color_overrides = {0: "#ff0000"}
        reset()
        assert mw.view_3d_manager._plugin_bond_color_overrides == {}
