"""
Headless GUI tests for:
  - Encrypted Project      → PmeencPlugin (patch/unpatch with real QActions)
  - Metadata Saver         → MetadataSaverDialog + save/load handlers
  - PubChem Name Resolver  → MoleculeResolverDialog
  - Structural Updater     → SettingsDialog + StructuralUpdaterPlugin
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

ENC_PATH = PLUGINS_DIR / "Encrypted_Project" / "encrypted_project.py"
META_PATH = PLUGINS_DIR / "Metadata_Saver" / "metadata_saver.py"
PUBCHEM_PATH = PLUGINS_DIR / "PubChem_Name_Ressolver" / "pubchem_ressolver.py"
UPDATER_PATH = PLUGINS_DIR / "Structural_Updater" / "structural_updater.py"

with mock_chemistry_imports():
    _enc = load_plugin_for_gui(ENC_PATH)
    _meta = load_plugin_for_gui(META_PATH)
    _pubchem = load_plugin_for_gui(PUBCHEM_PATH)
    _updater = load_plugin_for_gui(UPDATER_PATH)


# ===========================================================================
# PmeencPlugin  (Encrypted Project)
# ===========================================================================


def _enc_plugin(mw=None):
    ctx = MagicMock()
    ctx.get_main_window.return_value = mw if mw is not None else MagicMock()
    plugin = _enc.PmeencPlugin(ctx)
    return plugin, ctx


class TestPmeencPlugin:
    def test_cryptography_flag_true_in_test_env(self):
        assert _enc.CRYPTOGRAPHY_AVAILABLE

    def test_initialize_registers_hooks(self, qapp):
        plugin, ctx = _enc_plugin()
        plugin.initialize()
        ctx.register_file_opener.assert_called_once()
        assert ctx.register_file_opener.call_args.args[0] == ".pmeenc"
        ctx.add_export_action.assert_called_once()
        ctx.register_document_reset_handler.assert_called_once()
        ctx.register_drop_handler.assert_called_once()

    def test_drop_accepts_pmeenc(self, qapp):
        plugin, _ = _enc_plugin()
        plugin.on_import = MagicMock()  # keep the deferred QTimer call harmless
        assert plugin.on_drop("project.pmeenc") is True

    def test_drop_rejects_other_extensions(self, qapp):
        plugin, _ = _enc_plugin()
        assert plugin.on_drop("project.pmeprj") is False

    def test_document_reset_clears_keys(self, qapp):
        plugin, ctx = _enc_plugin()
        plugin.current_key = b"key"
        plugin.current_salt = b"salt"
        ctx.current_molecule = None
        plugin.on_document_reset()
        assert plugin.current_key is None
        assert plugin.current_salt is None

    def test_document_reset_clears_colorizer_overrides(self, qapp):
        plugin, ctx = _enc_plugin()
        plugin.mw._plugin_color_overrides = {0: "#ff0000"}
        plugin.on_document_reset()
        assert plugin.mw._plugin_color_overrides == {}
        ctx.draw_molecule_3d.assert_called_once()

    @pytest.fixture
    def patched_widget(self, qapp):
        from PyQt6.QtGui import QAction
        from PyQt6.QtWidgets import QWidget

        w = QWidget()
        action = QAction("Save Project", w)
        w.addAction(action)
        calls = []
        w.save_project = lambda: calls.append("orig")
        plugin, _ = _enc_plugin(mw=w)
        yield plugin, w, action, calls
        w.destroy()

    def test_patch_method_replaces_and_reconnects(self, patched_widget):
        plugin, w, action, calls = patched_widget
        plugin._patch_method(
            "save_project",
            lambda mw, *a, **k: calls.append("patched"),
            ["Save Project"],
        )
        assert getattr(w.save_project, "_is_pmeenc_patch", False)
        action.trigger()
        assert calls == ["patched"]

    def test_patch_method_is_idempotent(self, patched_widget):
        plugin, w, action, calls = patched_widget
        patch_fn = lambda mw, *a, **k: calls.append("patched")
        plugin._patch_method("save_project", patch_fn, ["Save Project"])
        first = w.save_project
        plugin._patch_method("save_project", patch_fn, ["Save Project"])
        assert w.save_project is first

    def test_unpatch_method_restores_original(self, patched_widget):
        plugin, w, action, calls = patched_widget
        original = w.save_project
        plugin._patch_method(
            "save_project",
            lambda mw, *a, **k: calls.append("patched"),
            ["Save Project"],
        )
        plugin._unpatch_method("save_project", ["Save Project"])
        assert w.save_project is original
        action.trigger()
        assert calls == ["orig"]

    def test_patched_clear_all_forwards_to_original(self, qapp):
        plugin, _ = _enc_plugin()
        called = []
        plugin._original_clear_all = lambda *a, **k: called.append(a)
        plugin._patched_clear_all(plugin.mw, True)
        assert called == [(True,)]


# ===========================================================================
# Metadata Saver — handlers and MetadataSaverDialog
# ===========================================================================


class TestMetadataHandlers:
    def test_default_config_shape(self):
        cfg = _meta._default_config()
        assert cfg["enabled"] is True
        assert cfg["silence_notice"] is False
        assert set(cfg["enabled_fields"]) == {f["key"] for f in _meta.ALL_FIELDS}

    def test_collect_metadata_skips_disabled_fields(self):
        cfg = _meta._default_config()
        cfg["enabled_fields"] = {k: False for k in cfg["enabled_fields"]}
        cfg["enabled_fields"]["os_name"] = True
        meta = _meta.collect_metadata(cfg, None)
        assert list(meta) == ["os_name"]

    def test_collect_metadata_skips_empty_note(self):
        cfg = _meta._default_config()
        cfg["enabled_fields"]["note"] = True
        cfg["custom_note"] = "   "
        assert "note" not in _meta.collect_metadata(cfg, None)

    def test_on_save_project_returns_versioned_payload(self):
        old_ctx = _meta._CONTEXT
        _meta._CONTEXT = MagicMock()
        try:
            data = _meta.on_save_project()
        finally:
            _meta._CONTEXT = old_ctx
        assert data["plugin_version"] == _meta.PLUGIN_VERSION
        assert "saved_at" in data["metadata"]

    def test_on_load_project_roundtrip_and_reset(self):
        _meta.on_load_project({"metadata": {"username": "alice"}})
        assert _meta._LOADED_METADATA == {"username": "alice"}
        _meta.on_document_reset()
        assert _meta._LOADED_METADATA == {}

    def test_on_load_project_ignores_malformed_data(self):
        _meta.on_load_project({"metadata": {"username": "alice"}})
        _meta.on_load_project("not a dict")
        assert _meta._LOADED_METADATA == {}

    def test_initialize_registers_menu_and_handlers(self, qapp):
        old_ctx = _meta._CONTEXT
        ctx = MagicMock()
        try:
            _meta.initialize(ctx)
        finally:
            _meta._CONTEXT = old_ctx
        assert (
            ctx.add_menu_action.call_args.args[0] == "Settings/Metadata Saver..."
        )
        ctx.register_save_handler.assert_called_once()
        ctx.register_load_handler.assert_called_once()
        ctx.register_document_reset_handler.assert_called_once()


class TestMetadataSaverDialog:
    @pytest.fixture
    def dlg(self, qapp):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        d = _meta.MetadataSaverDialog(ctx, parent=None)
        yield d
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Metadata Saver — Settings"

    def test_master_toggle_checked_by_default(self, dlg):
        assert dlg._chk_master.isChecked()

    def test_one_checkbox_per_field(self, dlg):
        assert set(dlg._checkboxes) == {f["key"] for f in _meta.ALL_FIELDS}

    def test_checkbox_defaults_match_registry(self, dlg):
        for field in _meta.ALL_FIELDS:
            assert dlg._checkboxes[field["key"]].isChecked() == field["default"]

    def test_master_toggle_disables_fields(self, dlg):
        dlg._chk_master.setChecked(False)
        assert not dlg._fields_container.isEnabled()
        dlg._chk_master.setChecked(True)
        assert dlg._fields_container.isEnabled()

    def test_build_cfg_reflects_ui_state(self, dlg):
        dlg._checkboxes["username"].setChecked(False)
        dlg._note_edit.setText("hello")
        cfg = dlg._build_cfg_from_ui()
        assert cfg["enabled_fields"]["username"] is False
        assert cfg["custom_note"] == "hello"

    def test_preview_fills_viewer(self, dlg):
        dlg._on_preview()
        text = dlg._metadata_view.toPlainText()
        assert "os_name:" in text

    def test_preview_with_nothing_selected(self, dlg):
        for chk in dlg._checkboxes.values():
            chk.setChecked(False)
        dlg._on_preview()
        assert dlg._metadata_view.toPlainText() == "(no fields selected)"


# ===========================================================================
# MoleculeResolverDialog  (PubChem Name Resolver)
# ===========================================================================


class TestMoleculeResolverDialog:
    @pytest.fixture
    def dlg(self, qapp):
        ctx = MagicMock()
        ctx.get_main_window.return_value = None
        d = _pubchem.MoleculeResolverDialog(ctx)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "PubChem Name Resolver"

    def test_search_type_options(self, dlg):
        items = [dlg.combo_type.itemText(i) for i in range(dlg.combo_type.count())]
        assert items == ["Auto (Name/CAS)", "SMILES"]

    def test_table_columns(self, dlg):
        assert dlg.table.columnCount() == 3
        headers = [
            dlg.table.horizontalHeaderItem(i).text() for i in range(3)
        ]
        assert headers == ["Name/Synonym", "Formula", "SMILES"]

    def test_registered_as_main_panel(self, dlg):
        dlg.context.register_window.assert_called_with("main_panel", dlg)

    def test_empty_search_is_noop(self, dlg):
        dlg.line_input.setText("   ")
        dlg.run_search()
        assert dlg.btn_search.isEnabled()
        assert dlg.lbl_info.text() == "Enter a chemical identifier and click Search."

    def test_smiles_search_populates_table(self, dlg):
        dlg.combo_type.setCurrentText("SMILES")
        dlg.line_input.setText("CCO")
        dlg.run_search()
        assert dlg.table.rowCount() == 1
        assert dlg.table.item(0, 0).text() == "User Input SMILES"
        assert dlg.table.item(0, 2).text() == "CCO"
        assert "Found 1 candidates" in dlg.lbl_info.text()
        assert dlg.btn_search.isEnabled()

    def test_invalid_smiles_shows_error(self, dlg):
        dlg.combo_type.setCurrentText("SMILES")
        dlg.line_input.setText("not-a-smiles")

        # The mocked requests module's RequestException is not a real exception
        # class; the except clause needs one to let the ValueError pass through.
        class _FakeRequestException(Exception):
            pass

        with patch.object(_pubchem.Chem, "MolFromSmiles", return_value=None), \
                patch.object(
                    _pubchem.requests.exceptions,
                    "RequestException",
                    _FakeRequestException,
                ), \
                patch.object(_pubchem, "QMessageBox") as mb:
            dlg.run_search()
        mb.critical.assert_called_once()
        assert "Invalid SMILES" in mb.critical.call_args.args[2]
        assert dlg.lbl_info.text() == "Error occurred."

    def test_load_without_selection_warns(self, dlg):
        with patch.object(_pubchem, "QMessageBox") as mb:
            dlg.load_molecule()
        mb.warning.assert_called_once()

    def test_load_selected_row_uses_string_importer(self, dlg):
        dlg.candidates_data = [
            {"name": "Ethanol", "smiles": "CCO", "formula": "C2H6O"}
        ]
        dlg.update_table()
        dlg.table.selectRow(0)
        mw = MagicMock()
        dlg.context.get_main_window.return_value = mw
        with patch.object(_pubchem, "QMessageBox"):
            dlg.load_molecule()
        mw.string_importer_manager.load_from_smiles.assert_called_once_with("CCO")


# ===========================================================================
# Structural Updater — SettingsDialog + StructuralUpdaterPlugin
# ===========================================================================


class TestUpdaterSettingsDialog:
    @pytest.fixture
    def dlg(self, qapp):
        d = _updater.SettingsDialog(parent=None, current_enabled=True)
        yield d
        d.destroy()

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Structural Updater Settings"

    def test_checkbox_reflects_current_enabled(self, dlg):
        assert dlg.chk_enable.isChecked()

    def test_disabled_state_propagates(self, qapp):
        d = _updater.SettingsDialog(parent=None, current_enabled=False)
        assert not d.chk_enable.isChecked()
        d.destroy()

    def test_accept_picks_up_checkbox_state(self, dlg):
        dlg.chk_enable.setChecked(False)
        dlg.accept()
        assert dlg.enabled is False

    def test_reject_keeps_previous_value(self, dlg):
        dlg.chk_enable.setChecked(False)
        dlg.reject()
        assert dlg.enabled is True


class TestStructuralUpdaterPlugin:
    @pytest.fixture
    def plugin(self, qapp):
        from PyQt6.QtWidgets import QWidget

        # QTimer needs a real QObject parent; managers stay MagicMock.
        w = QWidget()
        w.init_manager = MagicMock()
        w.init_manager.convert_button.text.return_value = "Convert 2D to 3D"
        w.compute_manager = MagicMock()
        ctx = MagicMock()
        ctx.get_main_window.return_value = w
        ctx.current_molecule = None
        saved_originals = dict(_updater._ORIGINAL_METHODS)
        _updater._ORIGINAL_METHODS.clear()
        p = _updater.StructuralUpdaterPlugin(ctx)
        yield p
        p.timer.stop()
        _updater._ORIGINAL_METHODS.clear()
        _updater._ORIGINAL_METHODS.update(saved_originals)
        w.destroy()

    def test_enabled_by_default(self, plugin):
        assert plugin.enabled is True

    def test_menu_action_registered(self, plugin):
        args = plugin.context.add_menu_action.call_args.args
        assert args[0] == "Settings/Structural Updater..."

    def test_timer_running_at_one_second(self, plugin):
        assert plugin.timer.isActive()
        assert plugin.timer.interval() == 1000

    def test_patch_stores_originals(self, plugin):
        assert "trigger_conversion" in _updater._ORIGINAL_METHODS
        assert "on_calculation_finished" in _updater._ORIGINAL_METHODS

    def test_check_state_enters_apply_mode(self, plugin):
        atom = MagicMock()
        atom.HasProp.return_value = True
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 3
        mol.GetAtoms.return_value = [atom]
        plugin.context.current_molecule = mol
        plugin.check_state()
        assert plugin.apply_mode_active is True
        plugin.mw.init_manager.convert_button.setText.assert_called_with(
            "Apply 2D Changes to 3D"
        )

    def test_check_state_exits_apply_mode_without_molecule(self, plugin):
        plugin.apply_mode_active = True
        plugin.context.current_molecule = None
        plugin.check_state()
        assert plugin.apply_mode_active is False
        plugin.mw.init_manager.convert_button.setText.assert_called_with(
            "Convert 2D to 3D"
        )

    def test_check_state_skips_during_running_calculation(self, plugin):
        plugin.mw.init_manager.convert_button.text.return_value = "Halt conversion"
        plugin.apply_mode_active = False
        atom = MagicMock()
        atom.HasProp.return_value = True
        mol = MagicMock()
        mol.GetNumAtoms.return_value = 3
        mol.GetAtoms.return_value = [atom]
        plugin.context.current_molecule = mol
        plugin.check_state()
        assert plugin.apply_mode_active is False

    def test_trigger_falls_back_to_original_when_not_apply_mode(self, plugin):
        plugin.apply_mode_active = False
        plugin.new_trigger_conversion()
        _updater._ORIGINAL_METHODS["trigger_conversion"].assert_called_once()

    def test_trigger_respects_temp_mode_override(self, plugin):
        plugin.apply_mode_active = True
        plugin.mw._temp_conv_mode = "force_full"
        plugin.new_trigger_conversion()
        _updater._ORIGINAL_METHODS["trigger_conversion"].assert_called_once()
