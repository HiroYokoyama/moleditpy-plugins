"""
Tests for the xTB Optimizer plugin.

All heavy dependencies (PyQt6, rdkit, tblite, ase, numpy, …) are mocked out
via conftest.mock_optional_imports(), so the suite runs headlessly with no
chemistry libraries required.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pytest

from conftest import (
    FakeAtom,
    FakeBond,
    FakeConf,
    FakeMol,
    extract_function,
    load_plugin,
    make_context,
    mock_optional_imports,
)

PLUGIN_PATH = (
    Path(__file__).resolve().parents[1] / "plugins" / "XTB_Optimizer" / "xtb_optimizer.py"
)


# ---------------------------------------------------------------------------
# Module-level fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def plugin_mod():
    """Load the plugin once for the whole module (deps mocked)."""
    with mock_optional_imports():
        return load_plugin(PLUGIN_PATH)


# ---------------------------------------------------------------------------
# Metadata / constants
# ---------------------------------------------------------------------------


def test_plugin_name(plugin_mod):
    assert plugin_mod.PLUGIN_NAME == "xTB Optimizer"


def test_plugin_version_format(plugin_mod):
    """Version must be YYYY.MM.DD."""
    v = plugin_mod.PLUGIN_VERSION
    parts = v.split(".")
    assert len(parts) == 3, f"Expected YYYY.MM.DD, got {v!r}"
    year, month, day = parts
    assert len(year) == 4 and year.isdigit()
    assert 1 <= int(month) <= 12
    assert 1 <= int(day) <= 31


def test_plugin_dependencies(plugin_mod):
    deps = plugin_mod.PLUGIN_DEPENDENCIES
    assert "tblite" in deps
    assert "ase" in deps


def test_required_constants(plugin_mod):
    for attr in (
        "PLUGIN_NAME",
        "PLUGIN_VERSION",
        "PLUGIN_SUPPORTED_MOLEDITPY_VERSION",
        "PLUGIN_AUTHOR",
        "PLUGIN_DESCRIPTION",
    ):
        assert hasattr(plugin_mod, attr), f"Missing constant: {attr}"


# ---------------------------------------------------------------------------
# initialize() registration
# ---------------------------------------------------------------------------


def test_initialize_registers_menu_action():
    with mock_optional_imports():
        mod = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        mod.initialize(ctx)
        ctx.add_menu_action.assert_called_once()
        path_arg = ctx.add_menu_action.call_args[0][0]
        assert "xTB" in path_arg, (
            f"Expected 'xTB' in menu path, got {path_arg!r}"
        )


def test_initialize_sets_launch_fn():
    """initialize() must store a callable in the module-level _launch_fn."""
    with mock_optional_imports():
        mod = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        mod.initialize(ctx)
        assert callable(mod._launch_fn)


def test_run_after_initialize_calls_launch_fn():
    """run(mw) must delegate to _launch_fn."""
    with mock_optional_imports():
        mod = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        mod.initialize(ctx)

        called = []
        original_fn = mod._launch_fn
        mod._launch_fn = lambda: called.append(True)

        mw = MagicMock()
        mod.run(mw)

        assert called, "run(mw) did not call _launch_fn"
        mod._launch_fn = original_fn  # restore


# ---------------------------------------------------------------------------
# run_plugin() guard — no molecule
# ---------------------------------------------------------------------------


def test_run_plugin_warns_when_no_molecule():
    """run_plugin must show a warning and return early if no molecule loaded."""
    with mock_optional_imports():
        mod = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        ctx.current_mol = None

        # Patch QMessageBox so we don't need a real Qt
        import sys
        qmb = sys.modules.get("PyQt6.QtWidgets")
        if qmb is None:
            pytest.skip("PyQt6 not available (expected when mocked)")

        # With mocked imports, QMessageBox.warning is a MagicMock
        mod.run_plugin(ctx)
        # Should NOT have opened a window (get_window not called with intent to show)
        ctx.get_window.assert_not_called()


# ---------------------------------------------------------------------------
# XtbWorker class shape
# ---------------------------------------------------------------------------


def _source_has_class_method(class_name: str, method_name: str) -> bool:
    """Return True if class_name defines method_name in the plugin source (AST)."""
    import ast
    source = PLUGIN_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in ast.walk(node):
                if isinstance(item, ast.FunctionDef) and item.name == method_name:
                    return True
    return False


def _source_has_class_attr(class_name: str, attr_name: str) -> bool:
    """Return True if class_name assigns attr_name at class body level (AST)."""
    import ast
    source = PLUGIN_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, ast.Assign):
                    for t in item.targets:
                        if isinstance(t, ast.Name) and t.id == attr_name:
                            return True
    return False


def test_worker_has_required_signals():
    """XtbWorker class must define the three worker signals at class level."""
    for signal in ("log_message", "step_update", "finished"):
        assert _source_has_class_attr("XtbWorker", signal), (
            f"XtbWorker is missing signal assignment: {signal}"
        )


def test_worker_has_cancel_method():
    assert _source_has_class_method("XtbWorker", "cancel")


def test_worker_has_run_method():
    assert _source_has_class_method("XtbWorker", "run")


# ---------------------------------------------------------------------------
# XtbOptimizerDialog class shape (AST-based — base class is MagicMock)
# ---------------------------------------------------------------------------


def test_dialog_class_exists():
    import ast
    source = PLUGIN_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    class_names = {
        n.name for n in ast.walk(tree) if isinstance(n, ast.ClassDef)
    }
    assert "XtbOptimizerDialog" in class_names


def test_dialog_has_on_run_method():
    assert _source_has_class_method("XtbOptimizerDialog", "_on_run")


def test_dialog_has_on_cancel_method():
    assert _source_has_class_method("XtbOptimizerDialog", "_on_cancel")


def test_dialog_has_on_finished_method():
    assert _source_has_class_method("XtbOptimizerDialog", "_on_finished")


# ---------------------------------------------------------------------------
# _on_finished: coordinate application logic (extracted via AST)
# ---------------------------------------------------------------------------


def test_on_finished_failure_does_not_update_mol():
    """
    When success=False, _on_finished must NOT set context.current_mol or
    call push_undo_checkpoint.
    """
    ctx = MagicMock()
    ctx.current_mol = MagicMock()

    fn = extract_function(PLUGIN_PATH, "XtbOptimizerDialog", "_on_finished")

    self_mock = MagicMock()
    self_mock.context = ctx
    self_mock._set_running = MagicMock()
    self_mock.lbl_status = MagicMock()
    self_mock._worker = None
    self_mock._step_count = 0

    fn(self_mock, False, "Cancelled by user.")

    ctx.push_undo_checkpoint.assert_not_called()


# ---------------------------------------------------------------------------
# Singleton dialog pattern
# ---------------------------------------------------------------------------


def test_run_plugin_reuses_existing_window():
    """run_plugin should raise an existing window instead of creating a new one."""
    with mock_optional_imports():
        mod = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        mol_mock = MagicMock()
        ctx.current_mol = mol_mock

        existing_win = MagicMock()
        ctx.get_window.return_value = existing_win

        mod.run_plugin(ctx)

        existing_win.show.assert_called_once()
        existing_win.raise_.assert_called_once()
        # No new dialog should be registered
        ctx.register_window.assert_not_called()


# ---------------------------------------------------------------------------
# Menu path
# ---------------------------------------------------------------------------


def test_menu_path_is_under_3d_edit():
    """The menu action must be registered under the '3D Edit/' prefix."""
    with mock_optional_imports():
        mod = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        mod.initialize(ctx)
        path_arg = ctx.add_menu_action.call_args[0][0]
        assert path_arg.startswith("3D Edit/"), (
            f"Expected path to start with '3D Edit/', got {path_arg!r}"
        )


# ---------------------------------------------------------------------------
# run(mw) host unwrapping
# ---------------------------------------------------------------------------


def test_run_unwraps_host_attribute():
    """run(mw) must use mw.host when the attribute exists."""
    with mock_optional_imports():
        mod = load_plugin(PLUGIN_PATH)
        ctx = make_context()
        mod.initialize(ctx)

        calls = []
        mod._launch_fn = lambda: calls.append(True)

        host = MagicMock()
        host.host = MagicMock()   # mw.host exists
        mod.run(host)

        assert calls, "run(mw) did not call _launch_fn when mw.host present"


def test_run_without_launch_fn_does_not_raise():
    """run(mw) with _launch_fn=None must be a silent no-op."""
    with mock_optional_imports():
        mod = load_plugin(PLUGIN_PATH)
        mod._launch_fn = None
        mw = MagicMock()
        mod.run(mw)   # must not raise


# ---------------------------------------------------------------------------
# PLUGIN_SUPPORTED_MOLEDITPY_VERSION format
# ---------------------------------------------------------------------------


def test_supported_version_specifier(plugin_mod):
    """Version specifier must include '>=4' to target the V4 API."""
    spec = plugin_mod.PLUGIN_SUPPORTED_MOLEDITPY_VERSION
    assert ">=4" in spec, (
        f"Expected '>=4' in PLUGIN_SUPPORTED_MOLEDITPY_VERSION, got {spec!r}"
    )


# ---------------------------------------------------------------------------
# XtbWorker — cancel sets internal flag (AST + extract_function)
# ---------------------------------------------------------------------------


def test_worker_cancel_sets_flag():
    """cancel() must set self._cancelled = True."""
    fn = extract_function(PLUGIN_PATH, "XtbWorker", "cancel")
    self_mock = MagicMock()
    self_mock._cancelled = False
    fn(self_mock)
    assert self_mock._cancelled is True


# ---------------------------------------------------------------------------
# _on_finished — success path calls push_undo_checkpoint
# ---------------------------------------------------------------------------


def test_on_finished_success_calls_undo_checkpoint():
    """
    When success=True, _on_finished must apply coordinates and call
    push_undo_checkpoint exactly once.
    """
    fn = extract_function(
        PLUGIN_PATH, "XtbOptimizerDialog", "_on_finished",
        extra_globals={"PLUGIN_NAME": "xTB Optimizer", "logger": MagicMock()},
    )
    fake_conf = MagicMock()
    fake_mol = MagicMock()
    fake_mol.GetConformer.return_value = fake_conf

    ctx = MagicMock()
    ctx.current_mol = fake_mol

    self_mock = MagicMock()
    self_mock.context = ctx
    self_mock._set_running = MagicMock()
    self_mock.lbl_status = MagicMock()
    self_mock._worker = None
    self_mock._step_count = 3
    self_mock.combo_method = MagicMock()
    self_mock.combo_method.currentText.return_value = "GFN2-xTB"

    # Two atoms, two positions
    new_positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]
    fn(self_mock, True, new_positions)

    ctx.push_undo_checkpoint.assert_called_once()
    ctx.refresh_3d_view.assert_called_once()


def test_on_finished_success_sets_each_atom_position():
    """All positions from the worker payload must be applied to the conformer."""
    fn = extract_function(
        PLUGIN_PATH, "XtbOptimizerDialog", "_on_finished",
        extra_globals={"PLUGIN_NAME": "xTB Optimizer", "logger": MagicMock()},
    )

    fake_conf = MagicMock()
    fake_mol = MagicMock()
    fake_mol.GetConformer.return_value = fake_conf

    ctx = MagicMock()
    ctx.current_mol = fake_mol

    self_mock = MagicMock()
    self_mock.context = ctx
    self_mock._set_running = MagicMock()
    self_mock.lbl_status = MagicMock()
    self_mock._worker = None
    self_mock._step_count = 1
    self_mock.combo_method = MagicMock()
    self_mock.combo_method.currentText.return_value = "GFN2-xTB"

    positions = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    fn(self_mock, True, positions)

    assert fake_conf.SetAtomPosition.call_count == len(positions)


# ---------------------------------------------------------------------------
# XtbOptimizerDialog — additional lifecycle methods (AST)
# ---------------------------------------------------------------------------


def test_dialog_has_close_event():
    assert _source_has_class_method("XtbOptimizerDialog", "closeEvent")


def test_dialog_has_build_ui():
    assert _source_has_class_method("XtbOptimizerDialog", "_build_ui")


def test_dialog_has_set_running():
    assert _source_has_class_method("XtbOptimizerDialog", "_set_running")


def test_dialog_has_append_log():
    assert _source_has_class_method("XtbOptimizerDialog", "_append_log")


def test_dialog_has_on_step_update():
    assert _source_has_class_method("XtbOptimizerDialog", "_on_step_update")


# ---------------------------------------------------------------------------
# PLUGIN_TAGS
# ---------------------------------------------------------------------------


def test_plugin_tags_contains_optimization(plugin_mod):
    assert "Optimization" in plugin_mod.PLUGIN_TAGS


def test_plugin_tags_is_minimal(plugin_mod):
    """Tags should be minimal — only ['Optimization']."""
    assert plugin_mod.PLUGIN_TAGS == ["Optimization"]

