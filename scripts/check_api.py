#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_api.py -- moleditpy-plugins repo-specific API checker.
Version: 2026.04.24

Pre-configured for the DEV_MAIN workspace layout:

  DEV_MAIN/
    python_molecular_editor/   <-- main app (default --app)
    moleditpy-plugins/
      REGISTRY/plugins.json
      plugins/                 <-- all plugin files (default --plugin)
      scripts/
        check_api.py           <-- this script

Usage (run from repo root or anywhere):
  python scripts/check_api.py                                    # all plugins, strict
  python scripts/check_api.py --registry                         # visible plugins only
  python scripts/check_api.py --registry --default-allowlist     # suppress manager false positives
  python scripts/check_api.py --registry --default-allowlist --mw-allowlist  # also suppress mw.X compat bridges
  python scripts/check_api.py --plugin plugins/Atom_Colorizer/atom_colorizer.py
  python scripts/check_api.py --check-context --show-api

Allowlist flags:
  --default-allowlist   Suppress manager attrs set via self.host.manager.X = ... (AST-invisible).
                        Safe to use routinely; these are confirmed false positives.
  --mw-allowlist        Also suppress direct mw.X legacy compat bridge attrs (mw.host, mw.view3d …).
                        Off by default — enabling it hides real V3 migration bugs.

All flags from the generic plugin_api_checker.py are supported.
Override defaults with --app / --plugin as needed.
"""

import ast
import argparse
import io
import json
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# Ensure Unicode output works on Windows terminals with narrow code pages.
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf-16"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ---------------------------------------------------------------------------
# Repo-relative defaults (computed from this script's location)
# ---------------------------------------------------------------------------
_SCRIPT_DIR  = Path(__file__).resolve().parent          # .../moleditpy-plugins/scripts/
_REPO_ROOT   = _SCRIPT_DIR.parent                       # .../moleditpy-plugins/
_DEFAULT_PLUGINS = (_REPO_ROOT / "plugins").resolve()
_DEFAULT_REGISTRY = (_REPO_ROOT / "REGISTRY" / "plugins.json").resolve()


# ---------------------------------------------------------------------------
# Registry loader
# ---------------------------------------------------------------------------

def load_visible_plugin_files(registry_path: Path, plugins_root: Path) -> list[Path]:
    """
    Parse REGISTRY/plugins.json, keep entries where visible==true,
    resolve relative downloadUrls to absolute paths, skip external URLs.
    Returns a sorted list of existing .py files.
    """
    try:
        entries = json.loads(registry_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        _die(f"Cannot read registry {registry_path}: {exc}")

    files: list[Path] = []
    skipped_external = 0
    skipped_missing = 0

    for entry in entries:
        if not entry.get("visible", False):
            continue
        url: str = entry.get("downloadUrl", "")
        if url.startswith("http://") or url.startswith("https://"):
            skipped_external += 1
            continue
        # Resolve relative path: downloadUrl is relative to REGISTRY/
        resolved = (registry_path.parent / url).resolve()
        if not resolved.exists():
            # Try resolving relative to plugins_root parent as fallback
            resolved2 = (plugins_root.parent / url.lstrip("./")).resolve()
            if resolved2.exists():
                resolved = resolved2
            else:
                skipped_missing += 1
                continue
        # Only scan .py files (zips / other formats are skipped)
        if resolved.suffix == ".py":
            files.append(resolved)

    if skipped_external:
        print(f"  Registry: skipped {skipped_external} external URL(s) (not local)")
    if skipped_missing:
        print(f"  Registry: skipped {skipped_missing} entry/entries (file not found)")

    return sorted(set(files))


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Issue:
    file: str
    line: int
    code: str
    message: str
    in_try: bool = False   # True if the access is inside a try: body

    def __str__(self) -> str:
        tag = "[try]" if self.in_try else "     "
        return f"  {tag}[{self.code}] line {self.line}: {self.message}"

    def key(self) -> tuple:
        return (self.file, self.line, self.code, self.message)


@dataclass
class APIInfo:
    """Public API surface detected from the main app."""
    mw_members: dict[str, str] = field(default_factory=dict)
    manager_members: dict[str, set[str]] = field(default_factory=dict)
    manager_class_names: dict[str, str] = field(default_factory=dict)
    context_members: set[str] = field(default_factory=set)


# ---------------------------------------------------------------------------
# Qt inherited methods -- always valid on QMainWindow, suppress false positives.
# ---------------------------------------------------------------------------
_QT_INHERITED = frozenset({
    "statusBar", "menuBar", "toolBar", "centralWidget", "setCentralWidget",
    "addToolBar", "addDockWidget", "removeDockWidget",
    "show", "hide", "close", "resize", "move", "setWindowTitle", "windowTitle", "setWindowIcon",
    "setMinimumSize", "setMaximumSize", "setFixedSize", "setSizePolicy",
    "update", "repaint", "raise_", "lower", "activateWindow",
    "setEnabled", "setDisabled", "isEnabled", "isVisible", "setVisible",
    "width", "height", "size", "pos", "geometry", "setGeometry",
    "parentWidget", "parent", "children", "findChild", "findChildren",
    "setAttribute", "testAttribute", "setStyleSheet",
    "setToolTip", "setStatusTip", "setWhatsThis",
    "installEventFilter", "removeEventFilter",
    "grabKeyboard", "releaseKeyboard",
    "connect", "disconnect", "blockSignals", "signalsBlocked",
    "deleteLater", "objectName", "setObjectName",
    "property", "setProperty", "metaObject",
    "thread", "moveToThread",
    "accept", "reject", "done",
})


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Allowlists
# ---------------------------------------------------------------------------
# _MANAGER_ALLOWLIST  -- safe, confirmed false positives: attrs assigned via
#   self.host.manager.X = ... patterns that the AST scanner cannot see.
#   Activated by --default-allowlist.
#
# _MW_ALLOWLIST       -- opt-in suppression for direct mw.X accesses that are
#   known legacy compat bridges (mw.host, mw.view3d, …).  NOT included in
#   --default-allowlist because hiding mw.X accesses masks real V3 migration
#   bugs.  Activated by --mw-allowlist.
# ---------------------------------------------------------------------------

_MANAGER_ALLOWLIST: dict[str, dict | set] = {
    "manager": {
        "state_manager": {"data"},
        "init_manager": {
            "scene",
            "view_2d",
            "current_file_path",
            "measurement_action",
            "analysis_action",
            "splitter",
            "edit_3d_action",
            "convert_button",
            "optimize_3d_button",
            "mode_actions",
        },
        "view_3d_manager": {
            "plotter",   # CustomQtInteractor -- set via self.host.view_3d_manager.plotter = ...
        },
    },
}

_MW_ALLOWLIST: dict[str, dict | set] = {
    "mw": {
        "host", "view3d", "string_importers", "apply_3d_settings",
        "main_window_ui_manager", "main_window_string_importers",
        "toggle_atom_info_display", "halt_all_calculations", "close_all_3d_edit_dialogs",
        "create_json_data", "load_from_json_data", "get_current_state", "set_state_from_data",
        "clear_2d_editor",
    },
}

def _load_site_allowlist(plugin_path: Path) -> dict:
    """
    Walk up from plugin_path looking for .moleditpy-api-allowlist.

    Each value can be a list (attrs only) or an object (attr -> reason string).
    The reason strings are ignored by the checker — they serve as inline comments.

      List form (no comments):
        { "mw": ["attr1", "attr2"] }

      Object form (with per-attr reasons):
        { "mw": { "attr1": "why it is safe to skip", "attr2": "..." } }

      Manager attrs follow the same pattern per manager:
        { "manager": { "state_manager": { "data": "injected by loader" } } }

    Returns a merged allowlist dict, or {} if no file is found.
    """
    search = plugin_path if plugin_path.is_dir() else plugin_path.parent
    for directory in [search, *search.parents]:
        candidate = directory / ".moleditpy-api-allowlist"
        if candidate.exists():
            try:
                data = json.loads(candidate.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError) as exc:
                print(f"Warning: could not parse {candidate}: {exc}", file=sys.stderr)
                return {}
            result: dict = {}
            if "mw" in data:
                mw_val = data["mw"]
                result["mw"] = set(mw_val.keys() if isinstance(mw_val, dict) else mw_val)
            if "manager" in data:
                result["manager"] = {
                    k: set(v.keys() if isinstance(v, dict) else v)
                    for k, v in data["manager"].items()
                }
            return result
        # Stop at repo root or filesystem root — don't leak into parent repos.
        if (directory / ".git").exists() or directory == directory.parent:
            break
    return {}


# Keep a combined alias for callers that want both (e.g. --default-allowlist + --mw-allowlist)
def _merge_allowlists(*lists) -> dict:
    merged: dict = {}
    for al in lists:
        for key, val in al.items():
            if key == "manager":
                mgr = merged.setdefault("manager", {})
                for mgr_name, members in val.items():
                    mgr.setdefault(mgr_name, set()).update(members)
            elif key == "mw":
                merged.setdefault("mw", set()).update(val)
    return merged


# ---------------------------------------------------------------------------
# Phase 1 -- Main App API extraction
# ---------------------------------------------------------------------------

class AppAPIExtractor:
    MANAGER_CLASS_HINTS = frozenset({
        "StateManager", "IOManager", "View3DManager", "Edit3DManager",
        "EditActionsManager", "ComputeManager", "DialogManager", "UIManager",
        "ExportManager", "MainInitManager", "StringImporterManager",
    })

    def __init__(self, app_root: Path, verbose: bool = False):
        self.app_root = app_root
        self.verbose = verbose
        self.api = APIInfo()

    def extract(self) -> APIInfo:
        mw_file = self._find_file_containing("class MainWindow")
        if not mw_file:
            _die(f"Could not find 'class MainWindow' under {self.app_root}")
        self._parse_main_window(mw_file)

        # Also scan all app .py files for `self.host.X = ...` patterns.
        # Managers set attributes on MainWindow this way (e.g. self.host.plugin_manager = ...).
        self._scan_host_assignments()

        ctx_file = self._find_file_containing("class PluginContext")
        if ctx_file:
            self._parse_plugin_context(ctx_file)

        if self.verbose:
            self._print_api_summary()

        return self.api

    def _scan_host_assignments(self):
        """
        Walk every .py file under app_root, find `self.host.X = ...` assignments,
        and add X to mw_members as a dynamically-set attribute.
        These are valid at runtime but invisible to MainWindow's own class body.
        """
        for f in sorted(self.app_root.rglob("*.py")):
            if "__pycache__" in str(f):
                continue
            try:
                tree = self._parse(f)
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                if not isinstance(node, ast.Assign):
                    continue
                for target in node.targets:
                    # Match self.host.X = ...
                    if (
                        isinstance(target, ast.Attribute)
                        and isinstance(target.value, ast.Attribute)
                        and isinstance(target.value.value, ast.Name)
                        and target.value.value.id == "self"
                        and target.value.attr == "host"
                        and not target.attr.startswith("_")
                    ):
                        self.api.mw_members.setdefault(target.attr, "dynamic")

    def _parse_main_window(self, path: Path):
        tree = self._parse(path)
        mw_class = self._find_class(tree, "MainWindow")
        if not mw_class:
            return
        for node in mw_class.body:
            if isinstance(node, ast.Assign):
                for t in node.targets:
                    if isinstance(t, ast.Name):
                        self.api.mw_members.setdefault(t.id, "signal/classattr")
            elif isinstance(node, ast.FunctionDef):
                self.api.mw_members[node.name] = _method_kind(node)
                if node.name == "__init__":
                    self._mine_mw_init(node, path.parent)

    def _mine_mw_init(self, init_fn: ast.FunctionDef, search_dir: Path):
        for stmt in ast.walk(init_fn):
            if not isinstance(stmt, ast.Assign):
                continue
            for target in stmt.targets:
                if not _is_self_attr(target):
                    continue
                attr_name = target.attr  # type: ignore[attr-defined]
                self.api.mw_members.setdefault(attr_name, "attr")
                if isinstance(stmt.value, ast.Call):
                    class_name = _call_class_name(stmt.value)
                    if class_name and class_name in self.MANAGER_CLASS_HINTS:
                        self.api.mw_members[attr_name] = f"manager:{class_name}"
                        self.api.manager_class_names[attr_name] = class_name
                        self._load_manager_class(attr_name, class_name, search_dir)

    def _load_manager_class(self, attr_name: str, class_name: str, search_dir: Path):
        py_file = self._find_file_containing(f"class {class_name}", under=search_dir)
        if not py_file:
            py_file = self._find_file_containing(f"class {class_name}")
        if not py_file:
            if self.verbose:
                print(f"  [warn] Could not find source for {class_name}")
            return
        tree = self._parse(py_file)
        cls = self._find_class(tree, class_name)
        if not cls:
            return
        members: set[str] = set()
        for node in cls.body:
            if isinstance(node, ast.FunctionDef):
                members.add(node.name)
                self._collect_self_attrs(node, members)
            elif isinstance(node, ast.Assign):
                for t in node.targets:
                    if isinstance(t, ast.Name):
                        members.add(t.id)
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                members.add(node.target.id)
        self.api.manager_members[attr_name] = members
        if self.verbose:
            print(f"  Loaded {class_name} ({attr_name}): {len(members)} members")

    @staticmethod
    def _collect_self_attrs(func: ast.FunctionDef, out: set[str]):
        for stmt in ast.walk(func):
            if isinstance(stmt, ast.Assign):
                for t in stmt.targets:
                    if _is_self_attr(t):
                        out.add(t.attr)  # type: ignore[attr-defined]
            elif isinstance(stmt, ast.AnnAssign):
                if _is_self_attr(stmt.target):
                    out.add(stmt.target.attr)  # type: ignore[attr-defined]

    def _parse_plugin_context(self, path: Path):
        tree = self._parse(path)
        cls = self._find_class(tree, "PluginContext")
        if not cls:
            return
        for node in cls.body:
            if isinstance(node, ast.FunctionDef):
                self.api.context_members.add(node.name)
            elif isinstance(node, ast.Assign):
                for t in node.targets:
                    if isinstance(t, ast.Name):
                        self.api.context_members.add(t.id)

    def _find_file_containing(self, needle: str, under: Optional[Path] = None) -> Optional[Path]:
        root = under or self.app_root
        for f in sorted(root.rglob("*.py")):
            try:
                if needle in f.read_text(encoding="utf-8", errors="ignore"):
                    return f
            except OSError:
                pass
        return None

    @staticmethod
    def _parse(path: Path) -> ast.Module:
        return ast.parse(path.read_text(encoding="utf-8", errors="ignore"), filename=str(path))

    @staticmethod
    def _find_class(tree: ast.Module, name: str) -> Optional[ast.ClassDef]:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == name:
                return node
        return None

    def _print_api_summary(self):
        print("\n=== Detected MainWindow API surface ===")
        for name, kind in sorted(self.api.mw_members.items()):
            print(f"  mw.{name}  [{kind}]")
        print("\n  Managers with parsed members:")
        for mgr, members in self.api.manager_members.items():
            cls = self.api.manager_class_names.get(mgr, "?")
            print(f"    mw.{mgr}  ({cls}, {len(members)} members)")
        if self.api.context_members:
            print(f"\n  PluginContext members: {len(self.api.context_members)}")
        print()


# ---------------------------------------------------------------------------
# Phase 2 -- Plugin scanning
# ---------------------------------------------------------------------------

_MW_VAR_NAMES = frozenset({
    "mw", "main_window", "parent_window", "app_window", "mainwindow",
})
_MW_GETTER_ATTRS = frozenset({"get_main_window"})
_SAFE_CALL_FUNCS = frozenset({"hasattr", "getattr", "setattr", "isinstance", "type"})
_CONTEXT_PARAM_NAMES = frozenset({"context", "ctx"})


class PluginFileChecker:
    def __init__(
        self,
        filepath: Path,
        api: APIInfo,
        check_context: bool = False,
        allowlist: Optional[dict] = None,
    ):
        self.filepath = filepath
        self.api = api
        self.check_context = check_context
        self._allowlist = allowlist or {}
        self.issues: list[Issue] = []
        self._mw_refs: set[str] = set()
        self._ctx_refs: set[str] = set()
        self._seen: set[tuple] = set()
        self._dynamic_mw_attrs: set[str] = set()
        self._dynamic_manager_attrs: dict[str, set[str]] = {}

    def check(self) -> list[Issue]:
        try:
            source = self.filepath.read_text(encoding="utf-8", errors="ignore")
            tree = ast.parse(source, filename=str(self.filepath))
        except SyntaxError as exc:
            self.issues.append(Issue(str(self.filepath), 0, "SYNTAX_ERROR", str(exc)))
            return self.issues
        self._try_body_lines = _collect_try_body_lines(tree)
        self._pass1_collect_aliases(tree)
        self._pass1_5_collect_dynamic_attrs(tree)
        self._pass2_check_accesses(tree)
        return self.issues

    def _pass1_5_collect_dynamic_attrs(self, tree: ast.Module):
        """
        AST Pass 1.5: Collect attributes dynamically assigned by the plugin to MainWindow 
        or its managers to suppress them as false positives (e.g. self.mw.my_custom_var = 1).
        """
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Attribute):
                        if self._is_mw_ref(target.value):
                            self._dynamic_mw_attrs.add(target.attr)
                        else:
                            mgr = self._is_mw_manager_ref(target.value)
                            if mgr:
                                self._dynamic_manager_attrs.setdefault(mgr, set()).add(target.attr)
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Attribute):
                if self._is_mw_ref(node.target.value):
                    self._dynamic_mw_attrs.add(node.target.attr)
                else:
                    mgr = self._is_mw_manager_ref(node.target.value)
                    if mgr:
                        self._dynamic_manager_attrs.setdefault(mgr, set()).add(node.target.attr)
            elif isinstance(node, ast.Call):
                if isinstance(node.func, ast.Name) and node.func.id == "setattr":
                    if len(node.args) >= 2:
                        arg0 = node.args[0]
                        arg1 = node.args[1]
                        attr_name = None
                        if isinstance(arg1, ast.Constant) and isinstance(arg1.value, str):
                            attr_name = arg1.value
                            
                        if attr_name:
                            if self._is_mw_ref(arg0):
                                self._dynamic_mw_attrs.add(attr_name)
                            else:
                                mgr = self._is_mw_manager_ref(arg0)
                                if mgr:
                                    self._dynamic_manager_attrs.setdefault(mgr, set()).add(attr_name)

    def _pass1_collect_aliases(self, tree: ast.Module):
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                is_mw = self._rhs_is_mw(node.value)
                is_ctx = self._rhs_is_ctx(node.value)
                for target in node.targets:
                    if is_mw:
                        self._record_mw(target)
                    if is_ctx:
                        self._record_ctx(target)
            elif isinstance(node, ast.AnnAssign) and node.value:
                if self._rhs_is_mw(node.value):
                    self._record_mw(node.target)
                if self._rhs_is_ctx(node.value):
                    self._record_ctx(node.target)
            elif isinstance(node, ast.FunctionDef):
                for arg in node.args.args:
                    if arg.arg in _MW_VAR_NAMES:
                        self._mw_refs.add(arg.arg)
                    if arg.arg in _CONTEXT_PARAM_NAMES:
                        self._ctx_refs.add(arg.arg)

    def _rhs_is_mw(self, node: ast.expr) -> bool:
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr in _MW_GETTER_ATTRS:
                return True
        if isinstance(node, ast.Name) and node.id in _MW_VAR_NAMES:
            return True
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
            if node.value.id == "self" and node.attr in _MW_VAR_NAMES:
                return True
        return False

    def _rhs_is_ctx(self, node: ast.expr) -> bool:
        return isinstance(node, ast.Name) and node.id in _CONTEXT_PARAM_NAMES

    def _record_mw(self, target: ast.expr):
        if isinstance(target, ast.Name):
            self._mw_refs.add(target.id)
        elif _is_self_attr(target):
            self._mw_refs.add(f"self.{target.attr}")  # type: ignore[attr-defined]

    def _record_ctx(self, target: ast.expr):
        if isinstance(target, ast.Name):
            self._ctx_refs.add(target.id)
        elif _is_self_attr(target):
            self._ctx_refs.add(f"self.{target.attr}")  # type: ignore[attr-defined]

    def _pass2_check_accesses(self, tree: ast.Module):
        safe = _collect_safe_positions(tree)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Attribute):
                continue
            if id(node) in safe:
                continue
            obj, attr = node.value, node.attr
            if attr.startswith("_") or attr in _QT_INHERITED:
                continue

            in_try = node.lineno in self._try_body_lines

            if self._is_mw_ref(obj):
                mw_allow = self._allowlist.get("mw", set())
                if attr not in self.api.mw_members and attr not in mw_allow and attr not in self._dynamic_mw_attrs:
                    self._add_issue(node.lineno, "UNKNOWN_MW_ATTR",
                        f"`{_repr(obj)}.{attr}` -- '{attr}' not found on MainWindow",
                        in_try)

            elif isinstance(obj, ast.Attribute) and self._is_mw_ref(obj.value):
                mgr_attr = obj.attr
                if mgr_attr in self.api.manager_members:
                    mgr_allow = self._allowlist.get("manager", {}).get(mgr_attr, set())
                    dyn_mgr_allow = self._dynamic_manager_attrs.get(mgr_attr, set())
                    if attr not in self.api.manager_members[mgr_attr] and attr not in mgr_allow and attr not in dyn_mgr_allow:
                        cls = self.api.manager_class_names.get(mgr_attr, mgr_attr)
                        self._add_issue(node.lineno, "UNKNOWN_MANAGER_ATTR",
                            f"`{_repr(obj.value)}.{mgr_attr}.{attr}` -- "
                            f"'{attr}' not found in {cls}",
                            in_try)

            elif self.check_context and self._is_ctx_ref(obj):
                if self.api.context_members and attr not in self.api.context_members:
                    self._add_issue(node.lineno, "UNKNOWN_CONTEXT_ATTR",
                        f"`{_repr(obj)}.{attr}` -- '{attr}' not found on PluginContext",
                        in_try)

    def _is_mw_ref(self, node: ast.expr) -> bool:
        if isinstance(node, ast.Name):
            return node.id in self._mw_refs or node.id in _MW_VAR_NAMES
        if _is_self_attr(node):
            return f"self.{node.attr}" in self._mw_refs  # type: ignore[attr-defined]
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr in _MW_GETTER_ATTRS:
                return True
        return False

    def _is_mw_manager_ref(self, node: ast.expr) -> Optional[str]:
        if isinstance(node, ast.Attribute) and self._is_mw_ref(node.value):
            return node.attr
        return None

    def _is_ctx_ref(self, node: ast.expr) -> bool:
        if isinstance(node, ast.Name):
            return node.id in self._ctx_refs or node.id in _CONTEXT_PARAM_NAMES
        if _is_self_attr(node):
            return f"self.{node.attr}" in self._ctx_refs  # type: ignore[attr-defined]
        return False

    def _add_issue(self, line: int, code: str, message: str, in_try: bool = False):
        issue = Issue(str(self.filepath), line, code, message, in_try)
        key = issue.key()
        if key not in self._seen:
            self._seen.add(key)
            self.issues.append(issue)


# ---------------------------------------------------------------------------
# AST helpers
# ---------------------------------------------------------------------------

def _is_self_attr(node: ast.expr) -> bool:
    return (
        isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == "self"
    )


def _method_kind(func: ast.FunctionDef) -> str:
    for dec in func.decorator_list:
        if isinstance(dec, ast.Name) and dec.id == "property":
            return "property"
        if isinstance(dec, ast.Attribute) and dec.attr in ("setter", "deleter"):
            return "property"
        if isinstance(dec, ast.Name) and dec.id in ("staticmethod", "classmethod"):
            return dec.id
    return "method"


def _call_class_name(call: ast.Call) -> Optional[str]:
    func = call.func
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return None


def _collect_try_body_lines(tree: ast.Module) -> set[int]:
    """
    Return the set of line numbers that fall inside a try: body (not except/finally).
    Issues on these lines may be caught at runtime -- severity = WARNING.
    Issues outside try bodies will definitely crash -- severity = CRITICAL.
    """
    lines: set[int] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Try):
            for stmt in node.body:
                for child in ast.walk(stmt):
                    if hasattr(child, "lineno"):
                        lines.add(child.lineno)
    return lines


def _collect_safe_positions(tree: ast.Module) -> set[int]:
    safe: set[int] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        func_name = (
            func.id if isinstance(func, ast.Name) else
            func.attr if isinstance(func, ast.Attribute) else None
        )
        if func_name in _SAFE_CALL_FUNCS and node.args:
            first = node.args[0]
            safe.add(id(first))
            if isinstance(first, ast.Attribute):
                safe.add(id(first))
    return safe


def _repr(node: ast.expr) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return f"{_repr(node.value)}.{node.attr}"
    if isinstance(node, ast.Call):
        func = node.func
        if isinstance(func, ast.Attribute):
            return f"{_repr(func.value)}.{func.attr}()"
        if isinstance(func, ast.Name):
            return f"{func.id}()"
    return "?"


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def collect_plugin_files(plugin_path: Path) -> list[Path]:
    if plugin_path.is_file():
        return [plugin_path]
    return [
        f for f in sorted(plugin_path.rglob("*.py"))
        if "__pycache__" not in str(f)
    ]


def run(args) -> int:
    app_root    = Path(args.app).expanduser().resolve()
    plugin_path = Path(args.plugin).expanduser().resolve()

    if not app_root.exists():
        _die(f"App path does not exist: {app_root}")
    if not plugin_path.exists():
        _die(f"Plugin path does not exist: {plugin_path}")

    print(f"App root  : {app_root}")
    print(f"Plugin(s) : {plugin_path}")

    # Phase 1 -- API extraction
    print("\nExtracting MainWindow API surface...")
    extractor = AppAPIExtractor(app_root, verbose=args.show_api)
    api = extractor.extract()
    print(
        f"  {len(api.mw_members)} MainWindow members, "
        f"{len(api.manager_members)} managers parsed"
        + (f", {len(api.context_members)} PluginContext members" if api.context_members else "")
    )

    # Phase 2 -- Collect files
    print()
    if args.registry:
        registry_path = Path(args.registry).expanduser().resolve()
        if not registry_path.exists():
            _die(f"Registry not found: {registry_path}")
        print(f"Registry  : {registry_path}")
        plugin_files = load_visible_plugin_files(registry_path, plugin_path)
        print(f"Scanning {len(plugin_files)} visible registered plugin file(s)...\n")
    else:
        plugin_files = collect_plugin_files(plugin_path)
        print(f"Scanning {len(plugin_files)} plugin file(s)...\n")

    parts = []
    if args.default_allowlist:
        parts.append(_MANAGER_ALLOWLIST)
    if args.mw_allowlist:
        parts.append(_MW_ALLOWLIST)

    site_allowlist = _load_site_allowlist(plugin_path)
    if site_allowlist:
        n_s_mw  = len(site_allowlist.get("mw", set()))
        n_s_mgr = sum(len(v) for v in site_allowlist.get("manager", {}).values())
        site_desc = []
        if n_s_mw:  site_desc.append(f"{n_s_mw} mw attr(s)")
        if n_s_mgr: site_desc.append(f"{n_s_mgr} manager attr(s)")
        print(f"Site allowlist (.moleditpy-api-allowlist): {', '.join(site_desc)}\n")
        parts.append(site_allowlist)

    allowlist = _merge_allowlists(*parts) if parts else {}
    if allowlist:
        n_mgr = sum(len(v) for v in allowlist.get("manager", {}).values())
        n_mw  = len(allowlist.get("mw", set()))
        suppressed = []
        if n_mgr: suppressed.append(f"{n_mgr} manager attr(s)")
        if n_mw:  suppressed.append(f"{n_mw} mw attr(s)")
        print(f"Allowlist active: {', '.join(suppressed)} suppressed\n")

    all_issues: list[tuple[Path, list[Issue]]] = []
    for pf in plugin_files:
        checker = PluginFileChecker(pf, api, check_context=args.check_context, allowlist=allowlist)
        issues = checker.check()
        if args.skip_try:
            issues = [i for i in issues if not i.in_try]
        if issues:
            all_issues.append((pf, issues))

    # Report
    if not all_issues:
        print("No issues found.")
        return 0

    total = sum(len(iss) for _, iss in all_issues)
    n_try = sum(1 for _, iss in all_issues for i in iss if i.in_try)

    suffix = f"  ({n_try} inside try blocks, shown with [try])" if n_try and not args.skip_try else ""
    print(f"Found {total} issue(s) in {len(all_issues)} file(s){suffix}:")
    print("=" * 72)

    display_root = plugin_path.parent if plugin_path.is_file() else plugin_path

    for filepath, issues in all_issues:
        try:
            display = filepath.relative_to(display_root)
        except ValueError:
            display = filepath
        print(f"\n{display}  ({len(issues)} issue(s))")
        print("-" * 60)
        for issue in sorted(issues, key=lambda i: i.line):
            print(issue)

    print("\n" + "=" * 72)
    print(f"Total: {total} issue(s)\n")

    code_counts = Counter(i.code for _, iss in all_issues for i in iss)
    for code, count in code_counts.most_common():
        label = {
            "UNKNOWN_MW_ATTR":      "Unknown MainWindow attribute",
            "UNKNOWN_MANAGER_ATTR": "Unknown manager attribute",
            "UNKNOWN_CONTEXT_ATTR": "Unknown PluginContext attribute",
            "SYNTAX_ERROR":         "Syntax error in plugin file",
        }.get(code, code)
        print(f"  {label:<40} {count}")

    return 1


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(2)


def main():
    parser = argparse.ArgumentParser(
        description="Check moleditpy-plugins for API disconnections against the main app.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--app",
        required=True,
        help="Path to the main app root (e.g. ../python_molecular_editor/)",
    )
    parser.add_argument(
        "--plugin",
        default=str(_DEFAULT_PLUGINS),
        help=f"Plugin file or directory (default: {_DEFAULT_PLUGINS})",
    )
    parser.add_argument(
        "--registry",
        nargs="?",
        const=str(_DEFAULT_REGISTRY),
        default=None,
        metavar="PATH",
        help=(
            "Scan only visible plugins from REGISTRY/plugins.json. "
            f"Defaults to {_DEFAULT_REGISTRY} when flag is given without a path."
        ),
    )
    parser.add_argument(
        "--default-allowlist", action="store_true",
        help=(
            "Suppress known manager false positives: runtime attrs assigned via "
            "self.host.manager.X patterns invisible to AST "
            "(init_manager.scene, state_manager.data, etc.)"
        ),
    )
    parser.add_argument(
        "--mw-allowlist", action="store_true",
        help=(
            "Also suppress known direct mw.X legacy compat attrs "
            "(mw.host, mw.view3d, etc.).  Off by default so mw.X issues remain visible."
        ),
    )
    parser.add_argument(
        "--skip-try", action="store_true",
        help=(
            "Hide issues whose access is inside a try: block. "
            "By default all issues are shown; those inside try blocks are tagged [try]."
        ),
    )
    parser.add_argument(
        "--check-context", action="store_true",
        help="Also check context.xxx accesses against the PluginContext API",
    )
    parser.add_argument(
        "--show-api", action="store_true",
        help="Print the detected MainWindow API surface before scanning",
    )
    args = parser.parse_args()
    sys.exit(run(args))


if __name__ == "__main__":
    main()
