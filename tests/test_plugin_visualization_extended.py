"""
Extended tests for visualization / rendering plugins:
  - Dark Mode Theme        (stylesheet content; autorun settings mutation)
  - Atom Colorizer         (save/load/reset handler payloads; apply_color logic)
  - Vector Viewer          (get_com; update_visualization parse/draw logic)
  - High Resolution Imager (run() guards; rejected dialog short-circuit)
  - VDW Radii Overlay      (partial/typed settings; initialize style registration)
  - Cube File Viewer (+Advanced, Mapped): MO-header cubes, Bohr->Angstrom
    conversion, drop handlers, stale-structure detach regression,
    FlexibleDoubleSpinBox formatting
  - Advanced Rendering     (gather/apply settings dict; sync_style_ui; re-init cleanup)

Qt / PyVista / RDKit stay mocked; real numpy is injected for parser math.
Methods on Qt-derived classes are extracted via AST (see
test_plugin_advanced_and_misc.py for the rationale).
"""

from __future__ import annotations

import ast
import json
import logging
import textwrap
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

DARK_PATH = PLUGINS_DIR / "Dark_Mode_Theme" / "dark_mode_plugin.py"
COLORIZER_PATH = PLUGINS_DIR / "Atom_Colorizer" / "atom_colorizer.py"
VECTOR_PATH = PLUGINS_DIR / "Vector_Viewer" / "vector_viewer.py"
HIRES_PATH = PLUGINS_DIR / "High_Resolution_Imager" / "high_res_imager.py"
VDW_PATH = PLUGINS_DIR / "VDW_Radii_Overlay" / "vdw_radii_overlay.py"
CUBE_PATH = PLUGINS_DIR / "Cube_File_Viewer" / "cube_viewer.py"
CUBE_ADV_PATH = PLUGINS_DIR / "Cube_File_Viewer_Advanced" / "cube_viewer_advanced.py"
MAPPED_PATH = PLUGINS_DIR / "Mapped_Cube_Viewer" / "mapped_cube_viewer.py"
ADVREND_PATH = PLUGINS_DIR / "Advanced_Rendering" / "advanced_rendering.py"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _extract_method_as_fn(path: Path, class_name: str, method_name: str,
                          extra_globals: dict | None = None):
    """
    Extract a class method as a standalone callable via AST.

    Qt base classes are MagicMock instances, so plugin classes deriving from
    them don't produce real types; pure methods must be exec'd in isolation.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if (isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
                        and item.name == method_name):
                    func_src = ast.get_source_segment(source, item)
                    if func_src:
                        local_ns: dict = {}
                        globs = {"logging": logging, **(extra_globals or {})}
                        exec(textwrap.dedent(func_src), globs, local_ns)
                        return local_ns[method_name]
    raise AssertionError(f"{class_name}.{method_name} not found in {path.name}")


def _load_real_np(plugin_path: Path):
    """Load plugin with heavy deps mocked, then inject real numpy."""
    with mock_optional_imports():
        mod = load_plugin(plugin_path)
    mod.np = np
    return mod


# ---------------------------------------------------------------------------
# Dark Mode Theme
# ---------------------------------------------------------------------------


class TestDarkModeStylesheet:
    @pytest.fixture(scope="class")
    def mod(self):
        with mock_optional_imports():
            return load_plugin(DARK_PATH)

    @pytest.mark.parametrize("selector", [
        "QWidget", "QMainWindow", "QMenuBar", "QMenu", "QTabBar::tab",
        "QToolBar", "QPushButton", "QLineEdit", "QCheckBox", "QSlider",
        "QScrollBar:vertical", "QHeaderView::section", "QGroupBox",
        "QStatusBar", "QSplitter::handle", "QTableView",
    ])
    def test_stylesheet_covers_widget(self, mod, selector):
        assert selector in mod.DARK_STYLESHEET

    def test_stylesheet_uses_dark_background(self, mod):
        assert "#2b2b2b" in mod.DARK_STYLESHEET

    def test_stylesheet_uses_accent_color(self, mod):
        assert "#3a6ea5" in mod.DARK_STYLESHEET

    def test_stylesheet_has_balanced_braces(self, mod):
        qss = mod.DARK_STYLESHEET
        assert qss.count("{") == qss.count("}")
        assert qss.count("{") > 20  # sanity: rich stylesheet


class TestDarkModeAutorun:
    def _fresh(self):
        with mock_optional_imports():
            return load_plugin(DARK_PATH)

    def test_applies_stylesheet_to_main_window(self):
        mod = self._fresh()
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        mw.setStyleSheet.assert_called_once_with(mod.DARK_STYLESHEET)

    def test_sets_dark_3d_background_in_settings(self):
        mod = self._fresh()
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        assert mw.init_manager.settings["background_color"] == "#2b2b2b"

    def test_sets_white_icon_foreground(self):
        mod = self._fresh()
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        assert mw.init_manager.settings["icon_foreground"] == "#FFFFFF"

    def test_triggers_apply_3d_settings(self):
        mod = self._fresh()
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        mw.view_3d_manager.apply_3d_settings.assert_called_once_with()

    def test_settings_error_does_not_raise(self):
        mod = self._fresh()
        mw = MagicMock()
        broken = MagicMock()
        broken.__setitem__.side_effect = RuntimeError("boom")
        mw.init_manager.settings = broken
        mod.autorun(mw)  # must be swallowed
        mw.setStyleSheet.assert_called_once()

    def test_status_message_via_context_after_initialize(self):
        mod = self._fresh()
        ctx = make_context()
        mod.initialize(ctx)
        ctx.get_main_window.return_value.setStyleSheet.assert_called_once()
        ctx.show_status_message.assert_called_once()
        msg = ctx.show_status_message.call_args[0][0]
        assert "Dark Mode" in msg

    def test_status_message_falls_back_to_statusbar(self):
        mod = self._fresh()
        assert mod._CONTEXT is None
        mw = MagicMock()
        mw.init_manager.settings = {}
        mod.autorun(mw)
        mw.statusBar.return_value.showMessage.assert_called_once()


# ---------------------------------------------------------------------------
# Atom Colorizer — persistence handlers
# ---------------------------------------------------------------------------


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
        mod.run(MagicMock())  # must not raise


class TestAtomColorizerApplyColor:
    def _apply(self, qmb=None):
        return _extract_method_as_fn(
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


# ---------------------------------------------------------------------------
# Vector Viewer
# ---------------------------------------------------------------------------


class TestVectorViewerGetCom:
    def _get_com(self):
        return _extract_method_as_fn(
            VECTOR_PATH, "VectorViewerPlugin", "get_com",
            extra_globals={"np": np},
        )

    def _mol(self, coords):
        mol = MagicMock()
        mol.GetNumAtoms.return_value = len(coords)
        conf = mol.GetConformer.return_value
        conf.GetAtomPosition.side_effect = [
            SimpleNamespace(x=c[0], y=c[1], z=c[2]) for c in coords
        ]
        return mol

    def test_mean_of_positions(self):
        fn = self._get_com()
        s = MagicMock()
        s.context.current_molecule = self._mol([(0, 0, 0), (2, 4, 6)])
        com = fn(s)
        assert np.allclose(com, [1.0, 2.0, 3.0])

    def test_no_molecule_returns_origin(self):
        fn = self._get_com()
        s = MagicMock()
        s.context.current_molecule = None
        assert np.allclose(fn(s), [0.0, 0.0, 0.0])

    def test_conformer_error_returns_origin(self):
        fn = self._get_com()
        s = MagicMock()
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 1
        mol.GetConformer.side_effect = RuntimeError("no conformer")
        s.context.current_molecule = mol
        assert np.allclose(fn(s), [0.0, 0.0, 0.0])


class TestVectorViewerUpdateVisualization:
    def _fn(self, pv_mock):
        return _extract_method_as_fn(
            VECTOR_PATH, "VectorViewerPlugin", "update_visualization",
            extra_globals={"np": np, "pv": pv_mock},
        )

    def _self(self, text, scale=1.0, reverse=False, com=(0.0, 0.0, 0.0)):
        s = MagicMock()
        s.vec_input.text.return_value = text
        s.reverse_chk.isChecked.return_value = reverse
        s.scale_spin.value.return_value = scale
        s.res_spin.value.return_value = 20
        s.opacity_spin.value.return_value = 0.5
        s.arrow_color.name.return_value = "#00ff00"
        s.get_com.return_value = np.array(com)
        s.vis_actor = None
        return s

    def test_no_plotter_returns_early(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1 0 0")
        s.context.plotter = None
        fn(s)
        pv_mock.Arrow.assert_not_called()

    def test_empty_text_returns_early(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("   ")
        fn(s)
        pv_mock.Arrow.assert_not_called()

    def test_two_components_rejected(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        fn(self._self("1.0 2.0"))
        pv_mock.Arrow.assert_not_called()

    def test_non_numeric_rejected(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        fn(self._self("a b c"))
        pv_mock.Arrow.assert_not_called()

    def test_zero_vector_rejected(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        fn(self._self("0 0 0"))
        pv_mock.Arrow.assert_not_called()

    def test_comma_separated_vector_drawn(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1, 0, 0", scale=2.0)
        fn(s)
        kwargs = pv_mock.Arrow.call_args.kwargs
        assert np.allclose(kwargs["direction"], [2.0, 0.0, 0.0])
        # Arrow centred on COM: start = com - scaled/2
        assert np.allclose(kwargs["start"], [-1.0, 0.0, 0.0])
        s.context.plotter.add_mesh.assert_called_once()
        assert s.context.plotter.add_mesh.call_args.kwargs["name"] == \
            "vector_viewer_arrow"
        s.context.plotter.render.assert_called_once()

    def test_reverse_flips_direction(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1 2 3", reverse=True)
        fn(s)
        kwargs = pv_mock.Arrow.call_args.kwargs
        assert np.allclose(kwargs["direction"], [-1.0, -2.0, -3.0])

    def test_old_actor_removed_before_redraw(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1 0 0")
        old_actor = MagicMock()
        s.vis_actor = old_actor
        fn(s)
        s.context.plotter.remove_actor.assert_called_once_with(old_actor)

    def test_actor_stored_after_draw(self):
        pv_mock = MagicMock()
        fn = self._fn(pv_mock)
        s = self._self("1 0 0")
        fn(s)
        assert s.vis_actor is s.context.plotter.add_mesh.return_value


class TestVectorViewerEntryPoints:
    def test_run_before_initialize_is_noop(self):
        with mock_optional_imports():
            mod = load_plugin(VECTOR_PATH)
        assert mod._launch_fn is None
        mod.run(MagicMock())  # must not raise

    def test_initialize_sets_launch_and_status(self):
        with mock_optional_imports():
            mod = load_plugin(VECTOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        assert callable(mod._launch_fn)
        ctx.show_status_message.assert_called_once()


# ---------------------------------------------------------------------------
# High Resolution Imager
# ---------------------------------------------------------------------------


class TestHighResImagerRun:
    def test_run_without_plugin_manager_returns(self):
        with mock_optional_imports():
            mod = load_plugin(HIRES_PATH)
            mod.PLUGIN_CONTEXT = MagicMock()
            mw = MagicMock(spec=[])  # no plugin_manager attribute
            mod.run(mw)  # must not raise, must not open dialog
        mod.PLUGIN_CONTEXT.get_main_window.assert_not_called()

    def test_run_without_context_returns(self):
        with mock_optional_imports():
            mod = load_plugin(HIRES_PATH)
            assert mod.PLUGIN_CONTEXT is None
            mod.run(MagicMock())  # must not raise

    def test_rejected_dialog_skips_file_dialog(self):
        with mock_optional_imports():
            mod = load_plugin(HIRES_PATH)
            ctx = make_context()
            # QDialog is a MagicMock: exec() returns a MagicMock which never
            # equals DialogCode.Accepted, i.e. the "user cancelled" path.
            mod.take_screenshot(ctx)
            assert mod.QFileDialog.getSaveFileName.call_count == 0
            ctx.plotter.screenshot.assert_not_called()


# ---------------------------------------------------------------------------
# VDW Radii Overlay — settings edge cases + initialize
# ---------------------------------------------------------------------------


class TestVDWSettingsEdgeCases:
    def _mod(self, tmp_path):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        mod.SETTINGS_FILE = str(tmp_path / "vdw.json")
        return mod

    def test_partial_file_keeps_other_default(self, tmp_path):
        mod = self._mod(tmp_path)
        Path(mod.SETTINGS_FILE).write_text(json.dumps({"occupancy": 0.7}))
        mod.load_settings()
        assert mod._vdw_settings["occupancy"] == pytest.approx(0.7)
        assert mod._vdw_settings["resolution"] == pytest.approx(0.125)

    def test_string_values_coerced_to_float(self, tmp_path):
        mod = self._mod(tmp_path)
        Path(mod.SETTINGS_FILE).write_text(
            json.dumps({"occupancy": "0.55", "resolution": "0.25"})
        )
        mod.load_settings()
        assert mod._vdw_settings["occupancy"] == pytest.approx(0.55)
        assert mod._vdw_settings["resolution"] == pytest.approx(0.25)

    def test_unknown_keys_ignored(self, tmp_path):
        mod = self._mod(tmp_path)
        Path(mod.SETTINGS_FILE).write_text(
            json.dumps({"bogus": 1, "occupancy": 0.4})
        )
        mod.load_settings()
        assert "bogus" not in mod._vdw_settings
        assert mod._vdw_settings["occupancy"] == pytest.approx(0.4)

    def test_save_writes_valid_json_with_both_keys(self, tmp_path):
        mod = self._mod(tmp_path)
        mod._vdw_settings["occupancy"] = 0.9
        mod._vdw_settings["resolution"] = 0.2
        mod.save_settings()
        on_disk = json.loads(Path(mod.SETTINGS_FILE).read_text())
        assert on_disk == {"occupancy": 0.9, "resolution": 0.2}

    def test_defaults(self, tmp_path):
        mod = self._mod(tmp_path)
        assert mod._vdw_settings["occupancy"] == pytest.approx(0.3)
        assert mod._vdw_settings["resolution"] == pytest.approx(0.125)


class TestVDWInitialize:
    def test_registers_vdw_overlay_style(self):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        ctx.register_3d_style.assert_called_once()
        name, drawer = ctx.register_3d_style.call_args[0]
        assert name == "vdw_overlay"
        assert drawer is mod.draw_vdw_overlay

    def test_run_without_context_is_noop(self):
        with mock_optional_imports():
            mod = load_plugin(VDW_PATH)
        assert mod.PLUGIN_CONTEXT is None
        mod.run(MagicMock())  # must not raise
