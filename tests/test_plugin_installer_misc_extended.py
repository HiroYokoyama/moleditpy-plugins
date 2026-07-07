"""
Second-wave extended tests: Plugin Installer table/version logic plus
PubChem Name Resolver, Python Console, Conformational Search,
All-Trans Optimizer, Complex Molecule Untangler and Advanced Rendering.

All tests run headlessly — PyQt6/rdkit/network fully mocked.
"""

import ast
import sys
import textwrap
import types
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

# Reuse the Plugin Installer module loaded by the existing suite so both files
# share one PyQt6 stub (QDialog/QThread are real classes there, which lets
# object.__new__ build bare window instances).
import test_plugin_installer as _tpi

PI = _tpi.PI

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

PUBCHEM_PATH = PLUGINS_DIR / "PubChem_Name_Ressolver" / "pubchem_ressolver.py"
CONSOLE_PATH = PLUGINS_DIR / "Python_Console" / "console.py"
CONF_SEARCH_PATH = PLUGINS_DIR / "Conformational_Search" / "conf_search.py"
ALL_TRANS_PATH = PLUGINS_DIR / "All-Trans_Optimizer" / "all-trans_optimizer.py"
UNTANGLER_PATH = (
    PLUGINS_DIR / "Complex_Molecule_Untangler" / "complex_molecule_untangler.py"
)
ADV_RENDER_PATH = PLUGINS_DIR / "Advanced_Rendering" / "advanced_rendering.py"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _extract_fn(path: Path, class_name: str | None, fn_name: str, extra_globals=None):
    """
    Extract a function or method from source via AST and exec it standalone.

    Needed for methods of Qt-derived classes: the mocked Qt bases make the
    class object a MagicMock, so methods can't be reached via the class.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    scope = tree
    if class_name is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == class_name:
                scope = node
                break
        else:
            raise AssertionError(f"class {class_name} not found in {path}")
    for node in scope.body:
        if isinstance(node, ast.FunctionDef) and node.name == fn_name:
            func_src = ast.get_source_segment(source, node)
            assert func_src is not None
            local_ns: dict = {}
            globs: dict = {}
            globs.update(extra_globals or {})
            exec(textwrap.dedent(func_src), globs, local_ns)
            return local_ns[fn_name]
    raise AssertionError(f"{fn_name} not found in {path}")


class _Item:
    """Recorder stand-in for QTableWidgetItem."""

    def __init__(self, text=""):
        self._text = str(text)
        self.tooltip = None
        self.background = None

    def setBackground(self, color):
        self.background = color

    def setToolTip(self, tip):
        self.tooltip = tip

    def text(self):
        return self._text


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


# ---------------------------------------------------------------------------
# Plugin Installer — populate_table statuses and compatibility filtering
# ---------------------------------------------------------------------------


class TestPopulateTable:
    def _make_installer(self, remote_data, installed, app_version="4.2.0"):
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.table = _FakeTable()
        inst.latest_app_version = "Unknown"  # skips the app-upgrade block
        inst.remote_data = remote_data
        inst.main_window = MagicMock()
        inst.main_window.plugin_manager.plugins = installed
        inst._pending_installs = {}
        inst.btn_upgrade_app = MagicMock()
        inst.btn_update_all = None
        inst.get_app_version = lambda: app_version
        inst._update_status_label = MagicMock()
        return inst

    def _statuses(self, inst):
        """Return {plugin name: status text} from the fake table."""
        out = {}
        for (row, col), item in inst.table.items.items():
            if col == 4:
                name_item = inst.table.items.get((row, 0))
                out[name_item.text()] = item.text()
        return out

    def _run(self, inst, monkeypatch):
        monkeypatch.setattr(PI, "QTableWidgetItem", _Item)
        inst.populate_table(silent=True)

    def test_invisible_uninstalled_plugin_skipped(self, monkeypatch):
        remote = [
            {"name": "Hidden", "version": "1.0.0", "visible": False,
             "downloadUrl": "https://x/h.py"},
        ]
        inst = self._make_installer(remote, installed=[])
        self._run(inst, monkeypatch)
        assert inst.table.rows == 0
        inst._update_status_label.assert_called_once_with(0)

    def test_incompatible_uninstalled_shows_incompatible_status(self, monkeypatch):
        remote = [
            {"name": "Future", "version": "9.0.0", "visible": True,
             "supported_moleditpy_version": ">=5.0.0",
             "downloadUrl": "https://x/f.py"},
        ]
        inst = self._make_installer(remote, installed=[], app_version="4.2.0")
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Future": "Incompatible"}
        status_item = inst.table.items[(0, 4)]
        assert ">=5.0.0" in status_item.tooltip
        assert "4.2.0" in status_item.tooltip

    def test_compatible_uninstalled_gets_install_button(self, monkeypatch):
        remote = [
            {"name": "Fresh", "version": "1.0.0", "visible": True,
             "supported_moleditpy_version": ">=4.0.0, <5.0.0",
             "downloadUrl": "https://x/fresh.py"},
        ]
        inst = self._make_installer(remote, installed=[])
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Fresh": "Not Installed"}
        assert (0, 5) in inst.table.cell_widgets  # Install button placed

    def test_update_available_counted(self, monkeypatch):
        remote = [
            {"name": "Old", "version": "2.0.0", "visible": True,
             "downloadUrl": "https://x/old.py"},
        ]
        installed = [{"name": "Old", "version": "1.0.0", "filepath": None}]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Old": "Update Available"}
        assert inst.updates_found is True
        inst._update_status_label.assert_called_once_with(1)

    def test_up_to_date_not_counted(self, monkeypatch):
        remote = [
            {"name": "Same", "version": "1.0.0", "visible": True,
             "downloadUrl": "https://x/same.py"},
        ]
        installed = [{"name": "Same", "version": "1.0.0", "filepath": None}]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Same": "Up to date"}
        inst._update_status_label.assert_called_once_with(0)

    def test_local_newer_than_registry(self, monkeypatch):
        remote = [
            {"name": "Dev", "version": "1.0.0", "visible": True,
             "downloadUrl": "https://x/dev.py"},
        ]
        installed = [{"name": "Dev", "version": "2.0.0", "filepath": None}]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Dev": "Newer"}

    def test_installed_but_not_in_registry(self, monkeypatch):
        installed = [{"name": "Private", "version": "0.1.0", "filepath": None}]
        inst = self._make_installer([], installed)
        self._run(inst, monkeypatch)
        assert self._statuses(inst) == {"Private": "Not in Registry"}
        assert (0, 5) not in inst.table.cell_widgets  # no action button

    def test_local_version_read_from_file(self, monkeypatch, tmp_path):
        f = tmp_path / "p.py"
        f.write_text('PLUGIN_VERSION = "3.3.3"\n', encoding="utf-8")
        remote = [
            {"name": "OnDisk", "version": "3.3.3", "visible": True,
             "downloadUrl": "https://x/p.py"},
        ]
        installed = [{"name": "OnDisk", "version": "0.0.1", "filepath": str(f)}]
        inst = self._make_installer(remote, installed)
        self._run(inst, monkeypatch)
        # AST-read version (3.3.3) wins over the stale metadata (0.0.1)
        assert self._statuses(inst) == {"OnDisk": "Up to date"}
        assert inst.table.items[(0, 2)].text() == "3.3.3"


# ---------------------------------------------------------------------------
# Plugin Installer — filter_plugins
# ---------------------------------------------------------------------------


class TestFilterPlugins:
    def _make_installer(self, rows):
        """rows: list of (name, author)."""
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.table = _FakeTable()
        inst.table.rows = len(rows)
        for r, (name, author) in enumerate(rows):
            inst.table.items[(r, 0)] = _Item(name)
            inst.table.items[(r, 1)] = _Item(author)
        inst.search_input = MagicMock()
        return inst

    def test_empty_filter_shows_all(self):
        inst = self._make_installer([("Alpha", "Ann"), ("Beta", "Bob")])
        inst.search_input.text.return_value = ""
        inst.filter_plugins()
        assert inst.table.hidden == {0: False, 1: False}

    def test_matches_name_case_insensitive(self):
        inst = self._make_installer([("Alpha", "Ann"), ("Beta", "Bob")])
        inst.search_input.text.return_value = "ALPHA"
        inst.filter_plugins()
        assert inst.table.hidden == {0: False, 1: True}

    def test_matches_author(self):
        inst = self._make_installer([("Alpha", "Ann"), ("Beta", "Bob")])
        inst.search_input.text.return_value = "bob"
        inst.filter_plugins()
        assert inst.table.hidden == {0: True, 1: False}

    def test_no_match_hides_all(self):
        inst = self._make_installer([("Alpha", "Ann"), ("Beta", "Bob")])
        inst.search_input.text.return_value = "zzz"
        inst.filter_plugins()
        assert inst.table.hidden == {0: True, 1: True}

    def test_missing_items_treated_as_empty(self):
        inst = object.__new__(PI.PluginInstallerWindow)
        inst.table = _FakeTable()
        inst.table.rows = 1  # row exists but has no items
        inst.search_input = MagicMock()
        inst.search_input.text.return_value = "x"
        inst.filter_plugins()  # must not raise
        assert inst.table.hidden == {0: True}


# ---------------------------------------------------------------------------
# Plugin Installer — app version / package name detection
# ---------------------------------------------------------------------------


def _fake_pkg(monkeypatch, top, version):
    """Install a fake <top>.utils.constants module tree with VERSION."""
    pkg = types.ModuleType(top)
    pkg.__path__ = []
    utils = types.ModuleType(f"{top}.utils")
    utils.__path__ = []
    constants = types.ModuleType(f"{top}.utils.constants")
    constants.VERSION = version
    pkg.utils = utils
    utils.constants = constants
    monkeypatch.setitem(sys.modules, top, pkg)
    monkeypatch.setitem(sys.modules, f"{top}.utils", utils)
    monkeypatch.setitem(sys.modules, f"{top}.utils.constants", constants)


def _block_pkg(monkeypatch, top):
    """Make ``import <top>`` raise ImportError even if really installed."""
    monkeypatch.setitem(sys.modules, top, None)


class TestAppVersionDetection:
    def _inst(self):
        return object.__new__(PI.PluginInstallerWindow)

    def test_get_app_version_from_moleditpy(self, monkeypatch):
        _fake_pkg(monkeypatch, "moleditpy", "4.9.9")
        assert self._inst().get_app_version() == "4.9.9"

    def test_get_app_version_from_linux_variant(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _fake_pkg(monkeypatch, "moleditpy_linux", "4.8.0")
        assert self._inst().get_app_version() == "4.8.0"

    def test_get_app_version_fallback_zero(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        assert self._inst().get_app_version() == "0.0.0"

    def test_package_name_moleditpy(self, monkeypatch):
        _fake_pkg(monkeypatch, "moleditpy", "4.9.9")
        assert self._inst()._get_package_name() == "moleditpy"

    def test_package_name_linux(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _fake_pkg(monkeypatch, "moleditpy_linux", "4.8.0")
        assert self._inst()._get_package_name() == "moleditpy-linux"

    def test_package_name_default_when_neither(self, monkeypatch):
        _block_pkg(monkeypatch, "moleditpy")
        _block_pkg(monkeypatch, "moleditpy_linux")
        assert self._inst()._get_package_name() == "moleditpy"


# ---------------------------------------------------------------------------
# PubChem Name Resolver — run_search / load_molecule
# ---------------------------------------------------------------------------


def _make_requests_stub():
    """Stub requests module with a real RequestException class."""
    m = types.ModuleType("requests")

    class RequestException(Exception):
        pass

    m.exceptions = SimpleNamespace(RequestException=RequestException)
    m.get = MagicMock()
    return m


def _pubchem_search_env(requests_stub, chem=None):
    globs = {
        "requests": requests_stub,
        "Chem": chem or MagicMock(),
        "QMessageBox": MagicMock(),
        "QApplication": MagicMock(),
        "Qt": MagicMock(),
        "PLUGIN_NAME": "PubChem Name Resolver",
    }
    fn = _extract_fn(PUBCHEM_PATH, "MoleculeResolverDialog", "run_search", globs)
    return fn, globs


def _pubchem_self(query, search_type="Auto (Name/CAS)"):
    s = SimpleNamespace()
    s.line_input = MagicMock()
    s.line_input.text.return_value = query
    s.combo_type = MagicMock()
    s.combo_type.currentText.return_value = search_type
    s.lbl_info = MagicMock()
    s.btn_search = MagicMock()
    s.table = MagicMock()
    s.context = MagicMock()
    s.candidates_data = []
    s.update_table = MagicMock()
    return s


def _http_response(status_code, payload=None):
    return SimpleNamespace(status_code=status_code, json=lambda: payload or {})


class TestPubChemRunSearch:
    def test_empty_query_returns_immediately(self):
        req = _make_requests_stub()
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("   ")
        fn(s)
        s.btn_search.setEnabled.assert_not_called()
        req.get.assert_not_called()

    def test_http_200_populates_candidates(self):
        req = _make_requests_stub()
        req.get.return_value = _http_response(200, {
            "PropertyTable": {"Properties": [
                {"IsomericSMILES": "CCO", "Title": "Ethanol",
                 "MolecularFormula": "C2H6O"},
            ]}
        })
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("ethanol")
        fn(s)
        assert s.candidates_data == [
            {"name": "Ethanol", "smiles": "CCO", "formula": "C2H6O"}
        ]
        s.update_table.assert_called_once()
        assert "Found 1 candidates" in s.lbl_info.setText.call_args[0][0]

    def test_canonical_smiles_fallback(self):
        req = _make_requests_stub()
        req.get.return_value = _http_response(200, {
            "PropertyTable": {"Properties": [
                {"CanonicalSMILES": "CC", "Title": "Ethane",
                 "MolecularFormula": "C2H6"},
            ]}
        })
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("ethane")
        fn(s)
        assert s.candidates_data[0]["smiles"] == "CC"

    def test_any_smiles_key_fallback(self):
        req = _make_requests_stub()
        req.get.return_value = _http_response(200, {
            "PropertyTable": {"Properties": [
                {"ConnectivitySMILES": "CN", "Title": "X",
                 "MolecularFormula": "CH5N"},
            ]}
        })
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("methylamine")
        fn(s)
        assert s.candidates_data[0]["smiles"] == "CN"

    def test_http_404_reports_not_found(self):
        req = _make_requests_stub()
        req.get.return_value = _http_response(404)
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("nosuchcompound")
        fn(s)
        assert s.lbl_info.setText.call_args[0][0] == "Not found in PubChem search."
        s.update_table.assert_not_called()

    def test_network_error_shows_dialog(self):
        req = _make_requests_stub()
        req.get.side_effect = req.exceptions.RequestException("offline")
        fn, globs = _pubchem_search_env(req)
        s = _pubchem_self("water")
        fn(s)
        globs["QMessageBox"].critical.assert_called_once()
        assert s.lbl_info.setText.call_args[0][0] == "Network error."

    def test_search_button_reenabled_after_error(self):
        req = _make_requests_stub()
        req.get.side_effect = req.exceptions.RequestException("offline")
        fn, _ = _pubchem_search_env(req)
        s = _pubchem_self("water")
        fn(s)
        s.btn_search.setEnabled.assert_called_with(True)

    def test_smiles_mode_valid_input(self):
        req = _make_requests_stub()
        chem = MagicMock()
        chem.rdMolDescriptors.CalcMolFormula.return_value = "C6H6"
        fn, _ = _pubchem_search_env(req, chem=chem)
        s = _pubchem_self("c1ccccc1", search_type="SMILES")
        fn(s)
        assert s.candidates_data == [
            {"name": "User Input SMILES", "smiles": "c1ccccc1", "formula": "C6H6"}
        ]
        req.get.assert_not_called()  # no network in SMILES mode

    def test_smiles_mode_invalid_input_reports_error(self):
        req = _make_requests_stub()
        chem = MagicMock()
        chem.MolFromSmiles.return_value = None
        fn, globs = _pubchem_search_env(req, chem=chem)
        s = _pubchem_self("not_smiles", search_type="SMILES")
        fn(s)
        globs["QMessageBox"].critical.assert_called_once()
        assert "Invalid SMILES" in globs["QMessageBox"].critical.call_args[0][2]
        assert s.lbl_info.setText.call_args[0][0] == "Error occurred."


class TestPubChemLoadMolecule:
    def _env(self):
        globs = {
            "QMessageBox": MagicMock(),
            "QApplication": MagicMock(),
            "Qt": MagicMock(),
            "PLUGIN_NAME": "PubChem Name Resolver",
        }
        fn = _extract_fn(
            PUBCHEM_PATH, "MoleculeResolverDialog", "load_molecule", globs
        )
        return fn, globs

    def _self(self, selected_row=0, smiles="CCO"):
        s = SimpleNamespace()
        item = MagicMock()
        item.row.return_value = selected_row
        s.table = MagicMock()
        s.table.selectedItems.return_value = [item] if selected_row is not None else []
        s.candidates_data = [{"name": "Ethanol", "smiles": smiles}]
        s.lbl_info = MagicMock()
        s.context = MagicMock()
        s.accept = MagicMock()
        return s

    def test_no_selection_warns(self):
        fn, globs = self._env()
        s = self._self(selected_row=None)
        fn(s)
        globs["QMessageBox"].warning.assert_called_once()
        s.accept.assert_not_called()

    def test_empty_smiles_warns(self):
        fn, globs = self._env()
        s = self._self(smiles="")
        fn(s)
        globs["QMessageBox"].warning.assert_called_once()
        s.accept.assert_not_called()

    def test_success_loads_via_string_importer(self):
        fn, globs = self._env()
        s = self._self()
        mw = MagicMock()
        s.context.get_main_window.return_value = mw
        fn(s)
        mw.string_importer_manager.load_from_smiles.assert_called_once_with("CCO")
        s.accept.assert_called_once()
        globs["QMessageBox"].information.assert_called_once()

    def test_missing_importer_reports_error(self):
        fn, globs = self._env()
        s = self._self()
        mw = SimpleNamespace()  # no string_importer_manager attribute
        s.context.get_main_window.return_value = mw
        fn(s)
        globs["QMessageBox"].critical.assert_called_once()
        s.accept.assert_not_called()


# ---------------------------------------------------------------------------
# Python Console — run_code execution and output capture
# ---------------------------------------------------------------------------


def _console_self(command, mol="MOL"):
    import code as code_mod

    s = SimpleNamespace()
    s.input_area = MagicMock()
    s.input_area.toPlainText.return_value = command
    s.outputs = []  # (text, color)
    s.append_output = lambda text, color=None: s.outputs.append((text, color))
    s.output_area = MagicMock()
    s.context = MagicMock()
    s._get_best_mol = lambda: mol
    s.local_scope = {}
    s.interpreter = code_mod.InteractiveInterpreter(s.local_scope)
    return s


def _run_code_fn():
    import io
    import traceback
    from contextlib import redirect_stderr, redirect_stdout

    globs = {
        "io": io,
        "traceback": traceback,
        "redirect_stdout": redirect_stdout,
        "redirect_stderr": redirect_stderr,
    }
    return _extract_fn(CONSOLE_PATH, "PythonConsoleDialog", "run_code", globs)


class TestConsoleRunCode:
    def test_empty_command_is_ignored(self):
        fn = _run_code_fn()
        s = _console_self("   \n  ")
        fn(s)
        s.input_area.append_history.assert_not_called()
        assert s.outputs == []

    def test_print_output_captured(self):
        fn = _run_code_fn()
        s = _console_self("print(21 * 2)")
        fn(s)
        s.output_area.append.assert_any_call("42")

    def test_expression_result_echoed_in_single_mode(self):
        fn = _run_code_fn()
        s = _console_self("1 + 1")
        fn(s)
        s.output_area.append.assert_any_call("2")

    # code.InteractiveInterpreter routes the traceback through sys.excepthook,
    # which pytest-qt's exception capture would report as a Qt-loop error.
    @pytest.mark.qt_no_exception_capture
    def test_exception_written_in_error_color(self):
        fn = _run_code_fn()
        s = _console_self("1/0")
        fn(s)
        err = [t for t, c in s.outputs if c == "#FF5252"]
        assert err and "ZeroDivisionError" in err[0]

    def test_incomplete_block_warns(self):
        fn = _run_code_fn()
        s = _console_self("def f():")
        fn(s)
        assert any("Incomplete" in t for t, _ in s.outputs)

    def test_multiline_runs_in_exec_mode(self):
        fn = _run_code_fn()
        s = _console_self("a = 6\nprint(a * 7)")
        fn(s)
        s.output_area.append.assert_any_call("42")

    def test_command_stored_in_history(self):
        fn = _run_code_fn()
        s = _console_self("x = 1")
        fn(s)
        s.input_area.append_history.assert_called_once_with("x = 1")
        s.input_area.clear.assert_called_once()

    def test_mol_and_mw_synced_into_scope(self):
        fn = _run_code_fn()
        s = _console_self("x = 1", mol="THEMOL")
        mw = object()
        s.context.get_main_window.return_value = mw
        fn(s)
        assert s.local_scope["mol"] == "THEMOL"
        assert s.local_scope["mw"] is mw

    def test_none_mol_warning_when_command_uses_mol(self):
        fn = _run_code_fn()
        s = _console_self("mol", mol=None)
        fn(s)
        assert any("'mol' is None" in t for t, _ in s.outputs)

    def test_no_mol_warning_for_unrelated_command(self):
        fn = _run_code_fn()
        s = _console_self("x = 5", mol=None)
        fn(s)
        assert not any("'mol' is None" in t for t, _ in s.outputs)

    def test_input_echoed_with_prompt_markers(self):
        fn = _run_code_fn()
        s = _console_self("a = 1\nb = 2")
        fn(s)
        texts = [t for t, _ in s.outputs]
        assert ">>> a = 1" in texts
        assert "... b = 2" in texts


# ---------------------------------------------------------------------------
# Conformational Search — energy-window dedup and table content
# ---------------------------------------------------------------------------


def _conf_filter_fn():
    return _extract_fn(
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
        fn = _extract_fn(
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
# All-Trans Optimizer — torsion application paths
# ---------------------------------------------------------------------------


class TestAllTransRunPlugin:
    def _load(self):
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
        return mod

    def _context(self, mol):
        ctx = MagicMock()
        ctx.get_main_window.return_value = MagicMock()
        ctx.current_mol = mol
        return ctx

    def test_no_conformer_warns_and_stops(self):
        mod = self._load()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 0
        ctx = self._context(mol)
        mod.QMessageBox.reset_mock()
        mod.rdMolTransforms.reset_mock()
        mod.run_plugin(ctx)
        mod.QMessageBox.warning.assert_called_once()
        mod.rdMolTransforms.SetDihedralDeg.assert_not_called()

    def test_matches_set_to_180_and_view_refreshed(self):
        mod = self._load()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        conf = MagicMock()
        mol.GetConformer.return_value = conf
        mol.GetSubstructMatches.return_value = [(0, 1, 2, 3), (4, 5, 6, 7)]
        ctx = self._context(mol)
        mod.rdMolTransforms.reset_mock()
        mod.QMessageBox.reset_mock()

        mod.run_plugin(ctx)

        calls = mod.rdMolTransforms.SetDihedralDeg.call_args_list
        assert len(calls) == 2
        assert calls[0].args == (conf, 0, 1, 2, 3, 180.0)
        assert calls[1].args == (conf, 4, 5, 6, 7, 180.0)
        ctx.refresh_3d_view.assert_called_once()
        ctx.push_undo_checkpoint.assert_called_once()
        assert "2" in ctx.show_status_message.call_args[0][0]

    def test_no_matches_informs_without_refresh(self):
        mod = self._load()
        mol = MagicMock()
        mol.GetNumConformers.return_value = 1
        mol.GetSubstructMatches.return_value = []
        ctx = self._context(mol)
        mod.rdMolTransforms.reset_mock()
        mod.QMessageBox.reset_mock()

        mod.run_plugin(ctx)

        mod.QMessageBox.information.assert_called_once()
        mod.rdMolTransforms.SetDihedralDeg.assert_not_called()
        ctx.refresh_3d_view.assert_not_called()
        ctx.push_undo_checkpoint.assert_not_called()


# ---------------------------------------------------------------------------
# Complex Molecule Untangler — UntangleWorker.run Monte Carlo loop
# ---------------------------------------------------------------------------


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
    fn = _extract_fn(UNTANGLER_PATH, "UntangleWorker", "run", globs)
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


# ---------------------------------------------------------------------------
# Advanced Rendering — registered style drawers
# ---------------------------------------------------------------------------


class TestAdvancedRenderingStyleDrawers:
    def _drawers(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        return {
            call.args[0]: call.args[1]
            for call in ctx.register_3d_style.call_args_list
        }

    def test_four_styles_registered(self):
        drawers = self._drawers()
        assert set(drawers) == {
            "Ball & Stick (Advanced Rendering)",
            "CPK (Advanced Rendering)",
            "Wireframe (Advanced Rendering)",
            "Stick (Advanced Rendering)",
        }

    def test_drawer_uses_normalized_style_key(self):
        drawers = self._drawers()
        viewer = MagicMock()
        mw_obj = SimpleNamespace(
            view_3d_manager=MagicMock(), _adv_rendering_viewer=viewer
        )
        mw_obj.view_3d_manager.current_3d_style = "Ball & Stick (Advanced Rendering)"

        drawers["Ball & Stick (Advanced Rendering)"](mw_obj, "MOL")

        mw_obj.view_3d_manager.draw_standard_3d_style.assert_called_once_with(
            "MOL", style_override="ball_and_stick"
        )
        viewer.apply_pbr_forced.assert_called_once()
        viewer.update_lights.assert_called_once()
        viewer.sync_style_ui.assert_called_once_with(
            "Ball & Stick (Advanced Rendering)"
        )

    def test_cpk_drawer_key(self):
        drawers = self._drawers()
        mw_obj = SimpleNamespace(
            view_3d_manager=MagicMock(), _adv_rendering_viewer=MagicMock()
        )
        drawers["CPK (Advanced Rendering)"](mw_obj, "MOL")
        mw_obj.view_3d_manager.draw_standard_3d_style.assert_called_once_with(
            "MOL", style_override="cpk"
        )

    def test_drawer_without_viewer_does_not_crash(self):
        drawers = self._drawers()
        mw_obj = SimpleNamespace(view_3d_manager=MagicMock())  # no viewer attr
        drawers["Stick (Advanced Rendering)"](mw_obj, "MOL")
        mw_obj.view_3d_manager.draw_standard_3d_style.assert_called_once()

    def test_drawer_without_view_manager_does_not_crash(self):
        # Regression (2026.07.08): sync_style_ui read mw_obj.view_3d_manager
        # unguarded even though the draw call above is hasattr-guarded for it.
        drawers = self._drawers()
        viewer = MagicMock()
        mw_obj = SimpleNamespace(_adv_rendering_viewer=viewer)
        drawers["Wireframe (Advanced Rendering)"](mw_obj, "MOL")
        viewer.sync_style_ui.assert_called_once_with("")
