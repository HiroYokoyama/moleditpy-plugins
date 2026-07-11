"""
Tests for the Atom Colorizer plugin: save/load/reset handler payloads,
apply_color logic.

AtomColorizerWindow inherits a Qt base, so apply_color is extracted from the
plugin source via AST and invoked with a fake self object.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
COLORIZER_PATH = PLUGINS_DIR / "Atom_Colorizer" / "atom_colorizer.py"


class TestAtomColorizerHandlers:
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
        mw.view_3d_manager._plugin_color_overrides = {0: "#ff0000", 2: "#00ff00"}
        assert save() == {"atom_colors": {"0": "#ff0000", "2": "#00ff00"}}

    def test_load_with_empty_data_is_noop(self):
        ctx = make_context()
        _, load, _ = self._handlers(ctx)
        load({})
        ctx.get_3d_controller.assert_not_called()
        ctx.refresh_3d_view.assert_not_called()

    def test_load_applies_each_color(self):
        ctx = make_context()
        _, load, _ = self._handlers(ctx)
        load({"atom_colors": {"1": "#00ff00", "3": "#0000ff"}})
        controller = ctx.get_3d_controller.return_value
        controller.set_atom_color.assert_any_call(1, "#00ff00")
        controller.set_atom_color.assert_any_call(3, "#0000ff")
        ctx.refresh_3d_view.assert_called_once()

    def test_load_with_bad_index_does_not_raise(self):
        ctx = make_context()
        _, load, _ = self._handlers(ctx)
        load({"atom_colors": {"abc": "#00ff00"}})
        ctx.refresh_3d_view.assert_called_once()

    def test_reset_clears_overrides_and_indices_field(self):
        ctx = make_context()
        _, _, reset = self._handlers(ctx)
        mw = ctx.get_main_window.return_value
        overrides = {0: "#ff0000"}
        mw.view_3d_manager._plugin_color_overrides = overrides
        reset()
        assert overrides == {}
        ctx.get_window.return_value.le_indices.clear.assert_called_once()

    def test_run_without_context_is_noop(self):
        with mock_optional_imports():
            mod = load_plugin(COLORIZER_PATH)
        assert mod.PLUGIN_CONTEXT is None
        mod.run(MagicMock())


class TestAtomColorizerApplyColor:
    def _apply(self, qmb=None):
        return extract_function(
            COLORIZER_PATH, "AtomColorizerWindow", "apply_color",
            extra_globals={"QMessageBox": qmb or MagicMock()},
        )

    def _self(self, text, num_atoms=3):
        s = MagicMock()
        s.le_indices.text.return_value = text
        s.current_color.name.return_value = "#112233"
        s.context.current_molecule.GetNumAtoms.return_value = num_atoms
        return s

    def test_empty_selection_warns(self):
        qmb = MagicMock()
        fn = self._apply(qmb)
        s = self._self("")
        fn(s)
        qmb.warning.assert_called_once()
        s.context.get_3d_controller.assert_not_called()

    def test_invalid_indices_warn(self):
        qmb = MagicMock()
        fn = self._apply(qmb)
        s = self._self("1, x, 2")
        fn(s)
        qmb.warning.assert_called_once()
        s.context.get_3d_controller.assert_not_called()

    def test_no_molecule_warns(self):
        qmb = MagicMock()
        fn = self._apply(qmb)
        s = self._self("0")
        s.context.current_molecule = None
        fn(s)
        qmb.warning.assert_called_once()

    def test_applies_color_to_valid_indices(self):
        fn = self._apply()
        s = self._self("0, 2", num_atoms=3)
        fn(s)
        controller = s.context.get_3d_controller.return_value
        controller.set_atom_color.assert_any_call(0, "#112233")
        controller.set_atom_color.assert_any_call(2, "#112233")
        assert controller.set_atom_color.call_count == 2
        s.context.refresh_3d_view.assert_called_once()

    def test_out_of_range_indices_skipped(self):
        fn = self._apply()
        s = self._self("1, 5, -1", num_atoms=3)
        fn(s)
        controller = s.context.get_3d_controller.return_value
        controller.set_atom_color.assert_called_once_with(1, "#112233")


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


class TestGetSelectionFromViewer:
    """get_selection_from_viewer merges 2D + 3D + measurement picks,
    deduplicated and sorted, into the indices field."""

    def _fn(self):
        return extract_function(
            COLORIZER_PATH, "AtomColorizerWindow", "get_selection_from_viewer"
        )

    def _self(self, sel_2d, sel_3d=None, picked=None, field_text=""):
        fake = MagicMock()
        fake.context.get_selected_atom_indices.return_value = sel_2d
        mw = MagicMock()
        mw.edit_3d_manager.selected_atoms_3d = sel_3d if sel_3d is not None else set()
        mw.edit_3d_manager.selected_atoms_for_measurement = picked if picked is not None else []
        fake.context.get_main_window.return_value = mw
        fake.le_indices = _FakeLineEdit(field_text)
        return fake

    def test_merges_sorts_and_dedups_all_sources(self):
        fake = self._self([0, 2], sel_3d={2, 5}, picked=(7, 5))
        self._fn()(fake)
        assert fake.le_indices.text() == "0,2,5,7"

    def test_no_settext_when_selection_unchanged(self):
        fake = self._self([1, 3], field_text="1,3")
        self._fn()(fake)
        assert fake.le_indices.set_calls == []

    def test_non_collection_3d_attribute_ignored(self):
        fake = self._self([4], sel_3d="garbage", picked=None)
        self._fn()(fake)
        assert fake.le_indices.text() == "4"

    def test_no_main_window_uses_2d_only(self):
        fake = self._self([9, 1])
        fake.context.get_main_window.return_value = None
        self._fn()(fake)
        assert fake.le_indices.text() == "1,9"


class TestAutoUpdateSelection:
    def _fn(self):
        return extract_function(
            COLORIZER_PATH, "AtomColorizerWindow", "_auto_update_selection"
        )

    def test_skips_update_while_field_focused(self):
        fake = MagicMock()
        fake.le_indices.hasFocus.return_value = True
        self._fn()(fake)
        fake.get_selection_from_viewer.assert_not_called()

    def test_updates_when_not_focused(self):
        fake = MagicMock()
        fake.le_indices.hasFocus.return_value = False
        self._fn()(fake)
        fake.get_selection_from_viewer.assert_called_once()


class TestAtomColorizerLoadRefresh:
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
        load({"atom_colors": {}})
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
        mw.view_3d_manager._plugin_color_overrides = {0: "#ff0000"}
        reset()
        assert mw.view_3d_manager._plugin_color_overrides == {}
