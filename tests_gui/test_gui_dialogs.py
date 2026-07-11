"""
Headless GUI tests for plugin dialog widgets.

These tests use real PyQt6 (QT_QPA_PLATFORM=offscreen) rather than mocking
all Qt classes.  Chemistry/scientific libraries (rdkit, numpy, …) are still
replaced with MagicMock via mock_chemistry_imports() so no installed chemistry
stack is required.

Only visible plugins (registry visible=true) are tested here.

Run locally:
    QT_QPA_PLATFORM=offscreen pytest tests_gui/ -v
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PASTE_XYZ_PATH = PLUGINS_DIR / "Paste_XYZ" / "paste_xyz.py"

MS_NEO_PATH = PLUGINS_DIR / "MS_Spectrum_Simulation_Neo" / "ms_spectrum_neo.py"
MOPAC_PATH = PLUGINS_DIR / "Mopac_Input_Generator" / "mopac_input_generator.py"
GAMESS_PATH = PLUGINS_DIR / "Gamess_Input_Generator" / "gamess_input_generator.py"
PSI4_PATH = PLUGINS_DIR / "Psi4_Input_Generator" / "psi4_input_generator.py"
ATOM_COLORIZER_PATH = PLUGINS_DIR / "Atom_Colorizer" / "atom_colorizer.py"
PLUGIN_INSTALLER_PATH = PLUGINS_DIR / "Plugin_Installer" / "plugin_installer.py"
CONSOLE_PATH = PLUGINS_DIR / "Python_Console" / "console.py"


# ---------------------------------------------------------------------------
# Load each plugin once at collection time.
# Real Qt classes are used; chemistry deps are MagicMock.
# ---------------------------------------------------------------------------

with mock_chemistry_imports():
    _paste_xyz = load_plugin_for_gui(PASTE_XYZ_PATH)
    _ms_neo = load_plugin_for_gui(MS_NEO_PATH)
    _mopac = load_plugin_for_gui(MOPAC_PATH)
    _gamess = load_plugin_for_gui(GAMESS_PATH)
    _psi4 = load_plugin_for_gui(PSI4_PATH)
    _atom_colorizer = load_plugin_for_gui(ATOM_COLORIZER_PATH)
    _plugin_installer = load_plugin_for_gui(PLUGIN_INSTALLER_PATH)
    _console = load_plugin_for_gui(CONSOLE_PATH)


def _ms_context() -> MagicMock:
    """Minimal stub context for MSSpectrumDialog (Neo).

    get_main_window() returns None → dialog has no parent, sync disabled.
    current_molecule = None  → formula_input stays blank, no Chem calls.
    """
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# PasteXYZDialog  (visible plugin: "Paste XYZ")
# ===========================================================================


class TestPasteXYZDialog:
    """PasteXYZDialog — pure-Qt dialog; rdkit guarded by try/except."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _paste_xyz.PasteXYZDialog(parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Paste XYZ"

    def test_default_size_wide_enough(self, dlg):
        assert dlg.width() >= 400

    def test_get_data_initially_empty(self, dlg):
        assert dlg.get_data() == ""

    def test_get_data_returns_typed_text(self, dlg):
        dlg.text_edit.setPlainText("C 0.0 0.0 0.0\nH 1.0 0.0 0.0")
        result = dlg.get_data()
        assert "C 0.0 0.0 0.0" in result
        assert "H 1.0 0.0 0.0" in result

    def test_multiline_xyz_round_trip(self, dlg):
        xyz = "C 0.0 0.0 0.0\nH 1.089 0.0 0.0\nH -0.363 1.027 0.0"
        dlg.text_edit.setPlainText(xyz)
        assert dlg.get_data().strip() == xyz

    def test_load_button_label(self, dlg):
        assert dlg.load_btn.text() == "Load"

    def test_cancel_button_label(self, dlg):
        assert dlg.cancel_btn.text() == "Cancel"

    def test_text_edit_placeholder_mentions_xyz(self, dlg):
        ph = dlg.text_edit.placeholderText().lower()
        assert "xyz" in ph

    def test_clear_resets_get_data(self, dlg):
        dlg.text_edit.setPlainText("something")
        dlg.text_edit.clear()
        assert dlg.get_data() == ""


# ===========================================================================
# RouteBuilderDialog  (visible plugin: "Gaussian Input Generator Neo")
# ===========================================================================


@pytest.mark.skip(reason="Gaussian Input Generator Neo retired; replaced by Gaussian Input Generator Pro")
class TestRouteBuilderDialog:
    """RouteBuilderDialog — tabbed QDialog with combo-driven route preview."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _gaussian.RouteBuilderDialog(parent=None, current_route="")
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Route Builder"

    def test_has_four_tabs(self, dlg):
        assert dlg.tabs.count() == 4

    def test_tab_method_basis_present(self, dlg):
        names = [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())]
        assert any("Method" in n or "Basis" in n for n in names)

    def test_tab_job_type_present(self, dlg):
        names = [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())]
        assert any("Job" in n for n in names)

    def test_tab_properties_present(self, dlg):
        names = [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())]
        assert any("Propert" in n for n in names)

    def test_preview_label_exists(self, dlg):
        assert hasattr(dlg, "preview_label")

    def test_preview_non_empty_after_init(self, dlg):
        assert dlg.preview_label.text() != ""

    def test_preview_contains_default_basis_set(self, dlg):
        assert "6-31G" in dlg.preview_label.text()

    def test_ok_button_label(self, dlg):
        assert dlg.btn_ok.text() == "Apply to Job"

    def test_cancel_button_label(self, dlg):
        assert dlg.btn_cancel.text() == "Cancel"

    def test_method_type_combo_has_dft(self, dlg):
        items = [dlg.method_type.itemText(i) for i in range(dlg.method_type.count())]
        assert any("DFT" in it for it in items)

    def test_dft_method_list_contains_b3lyp(self, dlg):
        dlg.method_type.setCurrentText("DFT")
        items = [dlg.method_name.itemText(i) for i in range(dlg.method_name.count())]
        assert "B3LYP" in items

    def test_mp2_method_list_contains_ccsd(self, dlg):
        dlg.method_type.setCurrentText("MP2")
        items = [dlg.method_name.itemText(i) for i in range(dlg.method_name.count())]
        assert "CCSD" in items

    def test_hf_method_list_contains_hf(self, dlg):
        dlg.method_type.setCurrentText("Hartree-Fock")
        items = [dlg.method_name.itemText(i) for i in range(dlg.method_name.count())]
        assert "HF" in items

    def test_switching_method_updates_preview(self, dlg):
        dlg.method_type.setCurrentText("Hartree-Fock")
        dlg.method_name.setCurrentText("HF")
        assert "HF" in dlg.preview_label.text()

    def test_basis_set_change_updates_preview(self, dlg):
        dlg.method_type.setCurrentText("DFT")
        dlg.basis_set.setCurrentText("def2TZVP")
        assert "def2TZVP" in dlg.preview_label.text()

    def test_job_opt_only_shows_opt_group_hides_freq(self, dlg):
        # "Optimization Only (Opt)" is index 1
        dlg.job_type.setCurrentIndex(1)
        # Use isHidden() — isVisible() requires the window to be .show()n first
        assert not dlg.opt_group.isHidden()
        assert dlg.freq_group.isHidden()

    def test_job_freq_only_shows_freq_group_hides_opt(self, dlg):
        # "Frequency Only (Freq)" is index 2
        dlg.job_type.setCurrentIndex(2)
        assert dlg.opt_group.isHidden()
        assert not dlg.freq_group.isHidden()

    def test_job_sp_hides_both_groups(self, dlg):
        # "Single Point Energy (SP)" is index 3
        dlg.job_type.setCurrentIndex(3)
        assert dlg.opt_group.isHidden()
        assert dlg.freq_group.isHidden()

    def test_opt_freq_job_shows_both_groups(self, dlg):
        # "Optimization + Freq (Opt Freq)" is index 0
        dlg.job_type.setCurrentIndex(0)
        assert not dlg.opt_group.isHidden()
        assert not dlg.freq_group.isHidden()


# ===========================================================================
# MSSpectrumDialog  (visible plugin: "MS Spectrum Simulation Neo")
# ===========================================================================


class TestMSSpectrumDialogNeo:
    """MSSpectrumDialog (Neo) with mol=None context — exercises no-molecule UI."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ms_context()
        d = _ms_neo.MSSpectrumDialog(context=ctx)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "MS Spectrum Simulation Neo"

    def test_default_size_tall_enough(self, dlg):
        assert dlg.height() >= 400

    def test_formula_input_initially_empty(self, dlg):
        assert dlg.formula_input.text() == ""

    def test_export_button_exists(self, dlg):
        assert dlg.export_btn.text() == "Export to Image"

    def test_csv_button_exists(self, dlg):
        assert dlg.btn_export_csv.text() == "Export CSV"

    def test_sync_check_disabled_when_no_main_window(self, dlg):
        # With get_main_window()=None, sync_check must be disabled
        assert not dlg.sync_check.isEnabled()

    def test_charge_spinbox_default_positive(self, dlg):
        assert dlg.charge_spin.value() == 1

    def test_adduct_combo_populated(self, dlg):
        assert dlg.adduct_combo.count() > 0

    def test_formula_input_accepts_text(self, dlg):
        dlg.formula_input.setText("C6H12O6")
        assert dlg.formula_input.text() == "C6H12O6"


# ===========================================================================
# MopacSetupDialog  (visible plugin: "MOPAC Input Generator")
# ===========================================================================


class TestMopacSetupDialog:
    """MopacSetupDialog with mol=None — exercises the no-molecule UI path."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _mopac.MopacSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "MOPAC Input Generator"

    def test_default_keywords_field(self, dlg):
        assert dlg.keywords_edit.text() == "PM7 PRECISE"

    def test_default_title_field(self, dlg):
        assert dlg.title_edit.text() == "MOPAC Calculation"

    def test_charge_spinbox_range(self, dlg):
        assert dlg.charge_spin.minimum() == -10
        assert dlg.charge_spin.maximum() == 10

    def test_mult_spinbox_minimum_is_one(self, dlg):
        assert dlg.mult_spin.minimum() == 1

    def test_template_combo_populated(self, dlg):
        assert dlg.template_combo.count() > 0

    def test_nosym_checkbox_initially_unchecked(self, dlg):
        assert not dlg.chk_nosym.isChecked()

    def test_preview_area_is_editable(self, dlg):
        assert not dlg.preview_text.isReadOnly()

    def test_preset_combo_initially_empty(self, dlg):
        assert dlg.preset_combo.count() == 0


# ===========================================================================
# GamessSetupDialog  (visible plugin: "GAMESS Input Generator")
# ===========================================================================


class TestGamessSetupDialog:
    """GamessSetupDialog with mol=None."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _gamess.GamessSetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "GAMESS Input Generator"

    def test_run_type_default(self, dlg):
        assert dlg.run_type.currentText() == "OPTIMIZE"

    def test_scf_type_default(self, dlg):
        assert dlg.scf_type.currentText() == "RHF"

    def test_nosym_checked_by_default(self, dlg):
        assert dlg.chk_nosym.isChecked()

    def test_basis_ngauss_default(self, dlg):
        assert dlg.basis_ngauss.value() == 6

    def test_mem_spin_default(self, dlg):
        assert dlg.mem_spin.value() == 100


# ===========================================================================
# Psi4SetupDialog  (visible plugin: "Psi4 Input Generator")
# ===========================================================================


class TestPsi4SetupDialog:
    """Psi4SetupDialog with mol=None."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _psi4.Psi4SetupDialog(parent=None, mol=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Psi4 Input Generator"

    def test_method_default(self, dlg):
        assert dlg.method_combo.currentText() == "b3lyp"

    def test_basis_default(self, dlg):
        assert dlg.basis_combo.currentText() == "def2-svp"

    def test_reference_default(self, dlg):
        assert dlg.ref_combo.currentText() == "rks"

    def test_mem_spin_default(self, dlg):
        assert dlg.mem_spin.value() == 2

    def test_thread_spin_default(self, dlg):
        assert dlg.thread_spin.value() == 4


# ===========================================================================
# AtomColorizerWindow  (visible plugin: "Atom Colorizer")
# ===========================================================================


def _colorizer_context() -> MagicMock:
    """Minimal stub: get_main_window() returns None so no QObject parent issues."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


class TestAtomColorizerWindow:
    """AtomColorizerWindow with a no-main-window context."""

    @pytest.fixture
    def win(self, qapp):
        ctx = _colorizer_context()
        w = _atom_colorizer.AtomColorizerWindow(context=ctx)
        yield w
        w.sel_timer.stop()
        w.destroy()

    def test_creates_without_error(self, win):
        assert win is not None

    def test_window_title(self, win):
        assert win.windowTitle() == "Atom Colorizer"

    def test_indices_placeholder(self, win):
        assert win.le_indices.placeholderText() == "e.g. 0, 1, 5"

    def test_choose_color_button_exists(self, win):
        assert win.btn_color.text() == "Choose Color"

    def test_timer_is_active(self, win):
        assert win.sel_timer.isActive()

    def test_is_non_modal(self, win):
        assert not win.isModal()


# ===========================================================================
# PluginDetailsDialog  (inside Plugin Installer)
# ===========================================================================


def _make_details_dialog(*, dependencies=None, local_info=None) -> object:
    """Create a PluginDetailsDialog with safe defaults."""
    return _plugin_installer.PluginDetailsDialog(
        None,
        "Test Plugin",
        "Test Author",
        "2026.01.01",
        "A headless test plugin.",
        dependencies if dependencies is not None else [],
        local_info,
        None,
        ">=4.0.0",
    )


class TestPluginDetailsDialog:
    """PluginDetailsDialog — lightweight details popup, no network calls."""

    @pytest.fixture
    def dlg(self, qapp):
        d = _make_details_dialog()
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title_includes_plugin_name(self, dlg):
        assert "Test Plugin" in dlg.windowTitle()

    def test_no_error_with_empty_dependencies(self, qapp):
        d = _make_details_dialog(dependencies=[])
        assert d is not None
        d.destroy()

    def test_no_error_with_missing_dependency(self, qapp):
        d = _make_details_dialog(dependencies=["fake-nonexistent-pkg-xyz"])
        assert d is not None
        d.destroy()

    def test_no_error_with_installed_dependency(self, qapp):
        d = _make_details_dialog(dependencies=["pytest"])
        assert d is not None
        d.destroy()

    def test_target_file_stored(self, qapp):
        d = _plugin_installer.PluginDetailsDialog(
            None, "X", "A", "1.0", "Desc", [], None, "/some/path.py"
        )
        assert d.target_file == "/some/path.py"
        d.destroy()


# ===========================================================================
# PluginInstallerWindow  (visible plugin: "Plugin Installer")
# ===========================================================================


class TestPluginInstallerWindow:
    """PluginInstallerWindow — main window of the Plugin Installer.

    check_updates() is patched to prevent any network calls during __init__.
    """

    @pytest.fixture
    def installer(self, qapp):
        from unittest.mock import patch

        with patch.object(
            _plugin_installer.PluginInstallerWindow,
            "check_updates",
            lambda self: None,
        ):
            d = _plugin_installer.PluginInstallerWindow(None, auto_check=False)
        yield d
        d.destroy()

    def test_creates_without_error(self, installer):
        assert installer is not None

    def test_window_title(self, installer):
        assert installer.windowTitle() == "Plugin Installer"

    def test_table_has_six_columns(self, installer):
        assert installer.table.columnCount() == 6

    def test_search_input_exists(self, installer):
        assert installer.search_input.placeholderText() != ""

    def test_startup_checkbox_exists(self, installer):
        assert installer.chk_startup is not None


# ===========================================================================
# PythonConsoleDialog  (visible plugin: "Python Console")
# ===========================================================================


def _console_context() -> MagicMock:
    """Stub context with no main window and a MagicMock current_molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    return ctx


class TestPythonConsoleDialog:
    """PythonConsoleDialog — embedded Python REPL dialog."""

    @pytest.fixture
    def console(self, qapp):
        ctx = _console_context()
        d = _console.PythonConsoleDialog(context=ctx)
        yield d
        d.destroy()

    def test_creates_without_error(self, console):
        assert console is not None

    def test_window_title(self, console):
        assert console.windowTitle() == "MoleditPy Python Console"

    def test_output_area_is_readonly(self, console):
        assert console.output_area.isReadOnly()

    def test_input_area_exists(self, console):
        assert console.input_area is not None

    def test_interpreter_is_interactive(self, console):
        import code as _code
        assert isinstance(console.interpreter, _code.InteractiveInterpreter)

    def test_context_in_local_scope(self, console):
        assert "context" in console.local_scope


# ===========================================================================
# Plugin Installer — pure-function unit tests (no Qt needed)
# ===========================================================================


class TestPluginInstallerVersionCompat:
    """is_app_version_compatible(): pure Python, no Qt required."""

    def test_within_range(self):
        assert _plugin_installer.is_app_version_compatible("4.2.0", ">=4.0.0,<5.0.0")

    def test_equal_lower_bound(self):
        assert _plugin_installer.is_app_version_compatible("4.0.0", ">=4.0.0")

    def test_below_lower_bound(self):
        assert not _plugin_installer.is_app_version_compatible("3.9.9", ">=4.0.0")

    def test_equal_upper_bound_exclusive(self):
        assert not _plugin_installer.is_app_version_compatible("5.0.0", "<5.0.0")

    def test_wildcard_always_passes(self):
        assert _plugin_installer.is_app_version_compatible("1.0.0", "*")

    def test_empty_specifier_always_passes(self):
        assert _plugin_installer.is_app_version_compatible("1.0.0", "")

    def test_not_equal(self):
        assert not _plugin_installer.is_app_version_compatible("4.0.0", "!=4.0.0")

    def test_compatible_release_same_major(self):
        assert _plugin_installer.is_app_version_compatible("4.2.0", "~=4.0")

    def test_compatible_release_different_major(self):
        assert not _plugin_installer.is_app_version_compatible("5.0.0", "~=4.0")


class TestPluginInstallerParseDependency:
    """parse_dependency(): pure Python, no Qt required."""

    def test_name_only(self):
        name, spec = _plugin_installer.parse_dependency("numpy")
        assert name == "numpy"
        assert spec == ""

    def test_name_with_specifier(self):
        name, spec = _plugin_installer.parse_dependency("rdkit>=2022.03.5")
        assert name == "rdkit"
        assert "2022" in spec

    def test_empty_string_returns_empty_pair(self):
        name, spec = _plugin_installer.parse_dependency("")
        assert name == ""
        assert spec == ""

    def test_url_with_colon_is_rejected(self):
        name, spec = _plugin_installer.parse_dependency("http://example.com/package")
        assert name == ""

    def test_dash_in_name(self):
        name, spec = _plugin_installer.parse_dependency("my-package>=1.0")
        assert name == "my-package"


class TestPluginInstallerCheckDependency:
    """check_dependency_satisfied(): queries importlib.metadata, pure Python."""

    def test_installed_package_returns_true(self):
        assert _plugin_installer.check_dependency_satisfied("pytest")

    def test_nonexistent_package_returns_false(self):
        assert not _plugin_installer.check_dependency_satisfied(
            "fake-nonexistent-dep-xyz123"
        )

    def test_empty_string_returns_false(self):
        assert not _plugin_installer.check_dependency_satisfied("")


# ===========================================================================
# Plugin Installer — PluginInstallerWindow more widget tests
# ===========================================================================


class TestPluginInstallerWindowWidgets:
    """Additional widget-state tests for PluginInstallerWindow."""

    @pytest.fixture
    def installer(self, qapp):
        from unittest.mock import patch

        with patch.object(
            _plugin_installer.PluginInstallerWindow,
            "check_updates",
            lambda self: None,
        ):
            d = _plugin_installer.PluginInstallerWindow(None, auto_check=False)
        yield d
        d.destroy()

    def test_table_header_first_column(self, installer):
        assert installer.table.horizontalHeaderItem(0).text() == "Plugin Name"

    def test_table_header_action_column(self, installer):
        assert installer.table.horizontalHeaderItem(5).text() == "Action"

    def test_progress_bar_initially_hidden(self, installer):
        assert installer._progress_bar.isHidden()

    def test_upgrade_button_initially_hidden(self, installer):
        assert installer.btn_upgrade_app.isHidden()

    def test_update_all_initially_disabled(self, installer):
        assert not installer.btn_update_all.isEnabled()

    def test_table_empty_before_fetch(self, installer):
        assert installer.table.rowCount() == 0

    def test_remote_data_initially_empty(self, installer):
        assert installer.remote_data == []


# ===========================================================================
# Python Console — initialize/run entry point tests
# ===========================================================================


class TestPythonConsoleEntryPoints:
    """initialize() and run() entry-point behaviour without a full main window."""

    def test_initialize_sets_plugin_context(self, qapp):
        ctx = _console_context()
        _console.initialize(ctx)
        assert _console.PLUGIN_CONTEXT is ctx
        _console.PLUGIN_CONTEXT = None  # reset global after test

    def test_run_is_noop_when_context_is_none(self, qapp):
        _console.PLUGIN_CONTEXT = None
        mw = MagicMock()
        _console.run(mw)  # must not raise
        _console.PLUGIN_CONTEXT = None

    def test_run_creates_dialog_and_calls_get_window(self, qapp):
        ctx = _console_context()
        ctx.get_window.return_value = None  # no cached window
        _console.PLUGIN_CONTEXT = ctx
        mw = MagicMock()
        _console.run(mw)  # constructs a PythonConsoleDialog, calls show()
        ctx.get_window.assert_called_once_with("main_panel")
        _console.PLUGIN_CONTEXT = None


# ===========================================================================
# Python Console — ConsoleInput widget tests
# ===========================================================================


class TestConsoleInput:
    """ConsoleInput (QPlainTextEdit subclass) — history and placeholder tests."""

    @pytest.fixture
    def inp(self, qapp):
        w = _console.ConsoleInput()
        yield w
        w.destroy()

    def test_creates_without_error(self, inp):
        assert inp is not None

    def test_placeholder_text(self, inp):
        assert "Python" in inp.placeholderText() or "Enter" in inp.placeholderText()

    def test_history_initially_empty(self, inp):
        assert inp.history == []

    def test_append_history_stores_command(self, inp):
        inp.append_history("x = 1")
        assert inp.history[-1] == "x = 1"

    def test_append_history_deduplicates_consecutive(self, inp):
        inp.append_history("cmd")
        inp.append_history("cmd")
        assert inp.history.count("cmd") == 1

    def test_history_index_advances_after_append(self, inp):
        inp.append_history("a")
        inp.append_history("b")
        assert inp.history_index == 2


# ===========================================================================
# Python Console — PythonHighlighter tests
# ===========================================================================


class TestPythonHighlighter:
    """PythonHighlighter (QSyntaxHighlighter) — highlighting rule counts."""

    @pytest.fixture
    def hl(self, qapp):
        from PyQt6.QtGui import QTextDocument
        doc = QTextDocument()
        h = _console.PythonHighlighter(parent=doc)
        yield h

    def test_creates_without_error(self, hl):
        assert hl is not None

    def test_has_highlighting_rules(self, hl):
        assert len(hl.highlighting_rules) > 0

    def test_rules_are_tuples(self, hl):
        for rule in hl.highlighting_rules:
            assert isinstance(rule, tuple) and len(rule) == 2
