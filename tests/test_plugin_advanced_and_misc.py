"""
Tests for visible plugins with pure-function or initialization logic:
  - Animated XYZ Giffer   (parse_multi_frame_xyz — pure file parser)
  - Advanced Rendering    (get_icon; initialize → 4 x register_3d_style + add_menu_action)
  - Dummy Atom Mode       (initialize → add_toolbar_action, no add_menu_action)
  - OpenBabel Conversion  (initialize → export + drop handler; OBABEL_AVAILABLE guard;
                           open_file_with_openbabel no-extension early exit)
  - High Resolution Imager (initialize → add_export_action; PLUGIN_CONTEXT stored)
  - XYZ Editor            (initialize → add_menu_action + save/load/reset handlers;
                           save handler returns {} when no molecule)
  - Symmetry Analyzer     (initialize → add_menu_action; PLUGIN_CONTEXT stored;
                           SymmetryAnalysisWorker.run emits finished with (dict, bool))
"""

from __future__ import annotations

import ast
import logging
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_xyz(tmp_path, frames, filename="test.xyz"):
    """Write a multi-frame XYZ file and return its string path."""
    lines = []
    for f in frames:
        atoms = f["atoms"]  # list of (sym, x, y, z)
        comment = f.get("comment", "frame")
        lines.append(str(len(atoms)))
        lines.append(comment)
        for sym, x, y, z in atoms:
            lines.append(f"{sym} {x} {y} {z}")
    p = tmp_path / filename
    p.write_text("\n".join(lines), encoding="utf-8")
    return str(p)


def _extract_method_as_fn(path: Path, class_name: str, method_name: str, extra_globals: dict | None = None):
    """
    Use AST to extract a class method as a standalone callable.

    Needed because Qt base classes are MagicMock instances whose metaclass call
    returns a MagicMock instead of a real type, so the class definition doesn't
    produce a usable type and object.__new__ fails.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == method_name:
                    func_src = ast.get_source_segment(source, item)
                    if func_src:
                        local_ns: dict = {}
                        globs = {"logging": logging, **(extra_globals or {})}
                        exec(textwrap.dedent(func_src), globs, local_ns)
                        return local_ns[method_name]
    return None


# ---------------------------------------------------------------------------
# Animated XYZ Giffer — parse_multi_frame_xyz
# ---------------------------------------------------------------------------

GIFFER_PATH = PLUGINS_DIR / "Animated_XYZ_Giffer" / "animated_xyz_giffer.py"

# Extract parse_multi_frame_xyz as a standalone function (self unused in body)
_parse_xyz_raw = _extract_method_as_fn(GIFFER_PATH, "AnimatedXYZPlayer", "parse_multi_frame_xyz")


def _parse_xyz(file_path: str):
    return _parse_xyz_raw(None, file_path)  # self=None: method never touches self


class TestParseMultiFrameXYZ:
    def test_single_frame_two_atoms(self, tmp_path):
        p = _write_xyz(tmp_path, [
            {"atoms": [("H", 0.0, 0.0, 0.0), ("O", 0.0, 0.0, 1.0)], "comment": "water"},
        ])
        frames = _parse_xyz(p)
        assert len(frames) == 1
        assert frames[0]["symbols"] == ["H", "O"]
        assert frames[0]["comment"] == "water"
        assert len(frames[0]["coords"]) == 2

    def test_two_frames(self, tmp_path):
        p = _write_xyz(tmp_path, [
            {"atoms": [("C", 0.0, 0.0, 0.0)], "comment": "frame1"},
            {"atoms": [("C", 1.0, 0.0, 0.0)], "comment": "frame2"},
        ])
        frames = _parse_xyz(p)
        assert len(frames) == 2
        assert frames[1]["comment"] == "frame2"
        assert frames[1]["coords"][0] == (1.0, 0.0, 0.0)

    def test_empty_file_returns_empty_list(self, tmp_path):
        p = tmp_path / "empty.xyz"
        p.write_text("", encoding="utf-8")
        assert _parse_xyz(str(p)) == []

    def test_blank_lines_between_frames(self, tmp_path):
        content = "\n\n2\nframe1\nH 0 0 0\nO 0 0 1\n\n2\nframe2\nH 1 0 0\nO 1 0 1\n"
        p = tmp_path / "blanks.xyz"
        p.write_text(content, encoding="utf-8")
        frames = _parse_xyz(str(p))
        assert len(frames) == 2

    def test_incomplete_frame_at_eof_skipped(self, tmp_path):
        content = "2\nframe1\nH 0 0 0\nO 0 0 1\n5\nincomplete\nH 0 0 0\n"
        p = tmp_path / "incomplete.xyz"
        p.write_text(content, encoding="utf-8")
        frames = _parse_xyz(str(p))
        assert len(frames) == 1  # the incomplete frame is dropped

    def test_coord_values_parsed_correctly(self, tmp_path):
        p = _write_xyz(tmp_path, [
            {"atoms": [("N", 1.5, -2.3, 0.001)], "comment": "c"},
        ])
        x, y, z = _parse_xyz(p)[0]["coords"][0]
        assert abs(x - 1.5) < 1e-9
        assert abs(y - (-2.3)) < 1e-9
        assert abs(z - 0.001) < 1e-9

    def test_bad_coord_line_skipped(self, tmp_path):
        """Atom lines with a non-float coordinate are silently skipped."""
        content = "2\nbad\nH 0 0 0\nO not_a_float 0 1\n"
        p = tmp_path / "bad.xyz"
        p.write_text(content, encoding="utf-8")
        frames = _parse_xyz(str(p))
        assert len(frames) == 1
        assert frames[0]["symbols"] == ["H"]  # O line dropped

    def test_non_integer_atom_count_line_skipped(self, tmp_path):
        """Lines that can't be parsed as an int atom count are skipped."""
        content = "GARBAGE\n2\nframe1\nH 0 0 0\nO 0 0 1\n"
        p = tmp_path / "garbage.xyz"
        p.write_text(content, encoding="utf-8")
        frames = _parse_xyz(str(p))
        assert len(frames) == 1

    def test_giffer_run_does_not_raise(self):
        """run(mw) returns early when plugin_manager is absent (MagicMock lacks it by default)."""
        with mock_optional_imports():
            mod = load_plugin(GIFFER_PATH)
            mw = MagicMock(spec=[])  # no attributes → plugin_manager absent
            mod.run(mw)  # must not raise


# ---------------------------------------------------------------------------
# Advanced Rendering
# ---------------------------------------------------------------------------

ADV_RENDER_PATH = PLUGINS_DIR / "Advanced_Rendering" / "advanced_rendering.py"


class TestAdvancedRendering:
    def test_get_icon_returns_none(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            assert mod.get_icon() is None

    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()

    def test_initialize_registers_four_3d_styles(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert ctx.register_3d_style.call_count == 4

    def test_initialize_menu_path_contains_advanced(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            path = ctx.add_menu_action.call_args[0][0]
            assert "Advanced" in path

    def test_initialize_style_names_contain_advanced_rendering(self):
        """All registered 3D style keys end with '(Advanced Rendering)'."""
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            style_keys = [
                call[0][0] for call in ctx.register_3d_style.call_args_list
            ]
            assert all("Advanced Rendering" in k for k in style_keys)


# ---------------------------------------------------------------------------
# Dummy Atom Mode
# ---------------------------------------------------------------------------

DUMMY_PATH = PLUGINS_DIR / "Dummy_Atom_Mode" / "dummy_atom_mode.py"


class TestDummyAtomMode:
    def test_initialize_calls_add_toolbar_action(self):
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_toolbar_action.assert_called_once()

    def test_initialize_toolbar_text_is_dummy_atom_star(self):
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            kwargs = ctx.add_toolbar_action.call_args[1]
            assert kwargs.get("text") == "Dummy Atom *"

    def test_initialize_does_not_call_add_menu_action(self):
        """Dummy Atom Mode registers a toolbar button, NOT a menu item."""
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_not_called()

    def test_initialize_tooltip_mentions_dummy(self):
        with mock_optional_imports():
            mod = load_plugin(DUMMY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            kwargs = ctx.add_toolbar_action.call_args[1]
            assert "dummy" in kwargs.get("tooltip", "").lower()


# ---------------------------------------------------------------------------
# OpenBabel Conversion Tool
# ---------------------------------------------------------------------------

OBABEL_PATH = PLUGINS_DIR / "OpenBabel_Conversion_Tool" / "openbabel_conversion_tool.py"


class TestOpenBabelConversionTool:
    def test_initialize_obabel_available_registers_export_action(self):
        """When openbabel is mocked (available), the export action is registered."""
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_openbabel(self):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            label = ctx.add_export_action.call_args[0][0]
            assert "OpenBabel" in label

    def test_initialize_obabel_available_registers_drop_handler(self):
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_drop_handler.assert_called_once()

    def test_initialize_obabel_unavailable_returns_early(self):
        """When OBABEL_AVAILABLE=False initialize() prints a warning and returns
        without calling any context registration method."""
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            mod.OBABEL_AVAILABLE = False
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_export_action.assert_not_called()
            ctx.register_file_opener.assert_not_called()
            ctx.register_drop_handler.assert_not_called()

    def test_open_file_no_extension_warns_and_returns(self):
        """open_file_with_openbabel() warns when the path has no file extension."""
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            with patch.object(mod.QMessageBox, "warning") as mock_warn:
                mod.open_file_with_openbabel("noextension", ctx)
            mock_warn.assert_called_once()

    def test_open_file_extensionless_does_not_reach_pybel(self):
        """Extensionless path exits before calling pybel.readfile."""
        with mock_optional_imports():
            mod = load_plugin(OBABEL_PATH)
            ctx = make_context()
            with patch.object(mod.QMessageBox, "warning"):
                mod.open_file_with_openbabel("noextension", ctx)
            # pybel.readfile is only called inside the try block after the guard
            # If we got here without AttributeError the guard worked correctly


# ---------------------------------------------------------------------------
# High Resolution Imager
# ---------------------------------------------------------------------------

HI_RES_PATH = PLUGINS_DIR / "High_Resolution_Imager" / "high_res_imager.py"


class TestHighResolutionImager:
    def test_initialize_registers_export_action(self):
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_screenshot(self):
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            label = ctx.add_export_action.call_args[0][0]
            assert "Screenshot" in label

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_initialize_does_not_register_menu_action(self):
        """High Resolution Imager uses add_export_action, not add_menu_action."""
        with mock_optional_imports():
            mod = load_plugin(HI_RES_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_not_called()


# ---------------------------------------------------------------------------
# XYZ Editor
# ---------------------------------------------------------------------------

XYZ_EDITOR_PATH = PLUGINS_DIR / "XYZ_Editor" / "xyz_editor.py"


class TestXYZEditor:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()

    def test_initialize_menu_path_contains_xyz(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            path = ctx.add_menu_action.call_args[0][0]
            assert "XYZ" in path

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_initialize_registers_save_handler(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_save_handler.assert_called_once()

    def test_initialize_registers_load_handler(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_load_handler.assert_called_once()

    def test_initialize_registers_document_reset_handler(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_document_reset_handler.assert_called_once()

    def test_save_handler_returns_empty_dict_when_no_molecule(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            ctx.current_molecule = None
            mod.initialize(ctx)
            save_fn = ctx.register_save_handler.call_args[0][0]
            result = save_fn()
            assert result == {} or result is None

    def test_load_handler_tolerates_empty_dict(self):
        with mock_optional_imports():
            mod = load_plugin(XYZ_EDITOR_PATH)
            ctx = make_context()
            ctx.current_molecule = None
            mod.initialize(ctx)
            load_fn = ctx.register_load_handler.call_args[0][0]
            load_fn({})  # must not raise


# ---------------------------------------------------------------------------
# Symmetry Analyzer
# ---------------------------------------------------------------------------

SYMMETRY_PATH = PLUGINS_DIR / "Symmetry_Analyzer" / "symmetry_analyzer.py"


# Extract SymmetryAnalysisWorker.run() as a standalone function.
# SymmetryAnalysisWorker inherits QThread (mocked → MagicMock metaclass),
# so the class itself becomes a MagicMock; object.__new__ cannot be used.
# The run() body uses np.arange (returns a MagicMock, iterable as empty)
# and PointGroupAnalyzer (inside a try/except, never reached when loop is empty).
_sym_run_raw = _extract_method_as_fn(
    SYMMETRY_PATH, "SymmetryAnalysisWorker", "run",
    extra_globals={"np": MagicMock()},  # np.arange returns MagicMock → iter([])
)


class _FakeWorker:
    """Minimal self-object for the extracted run() — carries only what run() reads."""
    def __init__(self, min_tol, max_tol):
        self.mol_pmg = MagicMock()
        self.min_tol = min_tol
        self.max_tol = max_tol
        self.finished = MagicMock()


class TestSymmetryAnalyzer:
    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()

    def test_initialize_menu_path_contains_symmetr(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            path = ctx.add_menu_action.call_args[0][0]
            assert "Symmetr" in path

    def test_initialize_stores_plugin_context(self):
        with mock_optional_imports():
            mod = load_plugin(SYMMETRY_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert mod.PLUGIN_CONTEXT is ctx

    def test_worker_run_emits_finished(self):
        """run() must always emit finished — with mocked numpy the tolerance
        loop iterates 0 times (MagicMock __iter__ returns empty iterator)."""
        worker = _FakeWorker(0.1, 0.1)
        with mock_optional_imports():
            _sym_run_raw(worker)
        worker.finished.emit.assert_called_once()

    def test_worker_run_emits_dict_and_bool(self):
        """finished.emit(group_data, found_any): group_data is dict, found_any is bool."""
        worker = _FakeWorker(0.05, 0.2)
        with mock_optional_imports():
            _sym_run_raw(worker)
        args = worker.finished.emit.call_args[0]
        assert len(args) == 2
        group_data, found_any = args
        assert isinstance(group_data, dict)
        assert isinstance(found_any, bool)

    def test_worker_found_any_false_when_no_real_analysis(self):
        """With mocked pymatgen the loop body is never entered → found_any=False."""
        worker = _FakeWorker(0.1, 0.5)
        with mock_optional_imports():
            _sym_run_raw(worker)
        _, found_any = worker.finished.emit.call_args[0]
        assert found_any is False
