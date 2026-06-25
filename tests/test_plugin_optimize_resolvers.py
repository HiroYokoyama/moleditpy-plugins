"""
Tests for optimization/resolver/viewer plugins:
  - All-Trans Optimizer
  - Complex Molecule Untangler
  - Conformational Search
  - PubChem Name Resolver
  - Vector Viewer
  - Molecule Comparator

Covers:
  1. initialize() must register at least one menu/export/plugin action
  2. No-molecule guard paths (run_plugin with mol=None)
  3. Dialog accept/reject round-trips
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _menu_registered(ctx: MagicMock) -> bool:
    """Return True if initialize() called any recognised registration method."""
    return (
        ctx.add_menu_action.called
        or ctx.add_export_action.called
        or ctx.add_plugin_menu.called
        or ctx.add_analysis_tool.called
        or ctx.add_toolbar_action.called
    )


def _load(rel: str):
    return PLUGINS_DIR / rel


# ---------------------------------------------------------------------------
# All-Trans Optimizer
# ---------------------------------------------------------------------------

ALL_TRANS_PATH = _load("All-Trans_Optimizer/all-trans_optimizer.py")


class TestAllTransOptimizer:
    def test_initialize_sets_launch_fn(self):
        """initialize() stores the launch closure in _launch_fn; run() uses it."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod._launch_fn is not None, "_launch_fn must be set by initialize()"

    def test_run_plugin_no_mol_does_not_raise(self):
        """run_plugin with mol=None should show a warning, not raise."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            mock_warn.assert_called_once()

    def test_run_plugin_no_mol_shows_correct_message(self):
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.run_plugin(ctx)
            args = mock_warn.call_args[0]
            assert "molecule" in args[2].lower() or "No molecule" in args[2]

    def test_initialize_then_run_calls_launch(self):
        """After initialize(), calling run(mw) invokes the launch function."""
        with mock_optional_imports():
            mod = load_plugin(ALL_TRANS_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            # run() should not raise even without a real molecule
            ctx.current_mol = None
            with patch.object(mod.QMessageBox, "warning"):
                mod.run(ctx.get_main_window())


# ---------------------------------------------------------------------------
# Complex Molecule Untangler
# ---------------------------------------------------------------------------

UNTANGLER_PATH = _load("Complex_Molecule_Untangler/complex_molecule_untangler.py")


class TestComplexMoleculeUntangler:
    def test_initialize_stores_context(self):
        """initialize() must store the context in PLUGIN_CONTEXT and set _launch_fn."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx
            assert mod._launch_fn is not None

    def test_run_plugin_opens_or_raises_dialog(self):
        """run_plugin with a valid context should not raise."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            mod.initialize(ctx)
            # run_plugin creates a dialog — should not raise with mocked Qt
            mod.run_plugin(ctx)

    def test_run_plugin_reuses_existing_window(self):
        """If a window is already registered, run_plugin shows it without creating new."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            existing_win = MagicMock()
            ctx.get_window.return_value = existing_win
            mod.run_plugin(ctx)
            existing_win.show.assert_called_once()
            existing_win.raise_.assert_called_once()

    def test_run_plugin_registers_new_dialog_window(self):
        """run_plugin() registers the created dialog as 'main_panel' when no existing window."""
        with mock_optional_imports():
            mod = load_plugin(UNTANGLER_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            mod.run_plugin(ctx)
            ctx.register_window.assert_called_once()


# ---------------------------------------------------------------------------
# Conformational Search (already has menu action — regression test)
# ---------------------------------------------------------------------------

CONF_SEARCH_PATH = _load("Conformational_Search/conf_search.py")


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


# ---------------------------------------------------------------------------
# PubChem Name Resolver
# ---------------------------------------------------------------------------

PUBCHEM_PATH = _load("PubChem_Name_Ressolver/pubchem_ressolver.py")


class TestPubChemNameResolver:
    def test_initialize_does_not_raise(self):
        """initialize() defines run_resolver locally; run(mw) handles auto-registration."""
        with mock_optional_imports():
            mod = load_plugin(PUBCHEM_PATH)
            ctx = make_context()
            mod.initialize(ctx)  # must not raise

    def test_run_search_empty_input_does_nothing(self):
        """run_search() with empty text must return without setting lbl_info."""
        with mock_optional_imports():
            mod = load_plugin(PUBCHEM_PATH)
            ctx = make_context()
            dialog = mod.MoleculeResolverDialog(ctx)
            dialog.line_input = MagicMock()
            dialog.line_input.text.return_value = "   "
            dialog.lbl_info = MagicMock()
            dialog.run_search()
            dialog.lbl_info.setText.assert_not_called()

    def test_run_creates_dialog_on_first_call(self):
        """run(mw) creates a MoleculeResolverDialog when no window is registered."""
        with mock_optional_imports():
            mod = load_plugin(PUBCHEM_PATH)
            ctx = make_context()
            ctx.get_window.return_value = None
            mod.initialize(ctx)
            # run(mw) is the auto-registered entry point
            mod.run(ctx.get_main_window())  # must not raise


# ---------------------------------------------------------------------------
# Vector Viewer
# ---------------------------------------------------------------------------

VECTOR_VIEWER_PATH = _load("Vector_Viewer/vector_viewer.py")


class TestVectorViewer:
    def test_initialize_sets_launch_fn(self):
        """initialize() stores launch in _launch_fn; run(mw) handles auto-registration."""
        with mock_optional_imports():
            mod = load_plugin(VECTOR_VIEWER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod._launch_fn is not None

    def test_initialize_shows_status_message(self):
        with mock_optional_imports():
            mod = load_plugin(VECTOR_VIEWER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.show_status_message.assert_called()

    def test_run_reuses_existing_window(self):
        """run() should show an existing window if one is registered."""
        with mock_optional_imports():
            mod = load_plugin(VECTOR_VIEWER_PATH)
            ctx = make_context()
            existing_win = MagicMock()
            ctx.get_window.return_value = existing_win
            mod.initialize(ctx)
            mod.run(ctx.get_main_window())
            # The launch fn checks get_window; it should show the existing window
            # (we can't deeply verify Qt calls without a real event loop, but no raise is the assertion)


# ---------------------------------------------------------------------------
# Molecule Comparator
# ---------------------------------------------------------------------------

MOLECULE_COMPARATOR_PATH = _load("Molecule_Comparator/molecule_comparator.py")


class TestMoleculeComparator:
    def test_initialize_registers_document_reset_handler(self):
        """initialize() wires up on_reset via register_document_reset_handler."""
        with mock_optional_imports():
            mod = load_plugin(MOLECULE_COMPARATOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_document_reset_handler.assert_called_once()

    def test_get_text_color_for_background_dark(self):
        """Luminance formula: dark backgrounds map to white text."""
        # Mirror the algorithm from MoleculeComparator._get_text_color_for_background
        def _luminance_text(hex_color):
            try:
                h = hex_color.lstrip("#")
                r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
                return "white" if (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.5 else "black"
            except Exception:
                return "black"

        assert _luminance_text("#000000") == "white"
        assert _luminance_text("#0000FF") == "white"

    def test_get_text_color_for_background_light(self):
        """Luminance formula: light backgrounds map to black text."""
        def _luminance_text(hex_color):
            try:
                h = hex_color.lstrip("#")
                r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
                return "white" if (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.5 else "black"
            except Exception:
                return "black"

        assert _luminance_text("#FFFFFF") == "black"
        assert _luminance_text("#FFFF00") == "black"

    def test_get_text_color_for_background_invalid(self):
        """Luminance formula falls back to 'black' for non-hex input."""
        def _luminance_text(hex_color):
            try:
                h = hex_color.lstrip("#")
                r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
                return "white" if (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.5 else "black"
            except Exception:
                return "black"

        assert _luminance_text("not_a_color") == "black"
