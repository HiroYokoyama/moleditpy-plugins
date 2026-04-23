#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plugin_api_checker.py -- AST-based MoleditPy plugin/main-app API disconnection finder.

Detects AttributeError-prone accesses like `mw.nonexistent_method()` or
`mw.view_3d_manager.bad_method()` in plugin code -- without running anything.

Usage
-----
  # Check all plugins in a directory against the main app:
  python scripts/plugin_api_checker.py --app ../python_molecular_editor/ \\
                                        --plugin plugins/

  # Check a single plugin file:
  python scripts/plugin_api_checker.py --app ../python_molecular_editor/ \\
                                        --plugin plugins/Atom_Colorizer/atom_colorizer.py

  # Also check context.xxx accesses against the PluginContext API:
  python scripts/plugin_api_checker.py --app ../python_molecular_editor/ \\
                                        --plugin plugins/ --check-context

  # Show the full API surface that was detected:
  python scripts/plugin_api_checker.py --app ../python_molecular_editor/ \\
                                        --plugin plugins/ --show-api
"""

import ast
import argparse
import io
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# Ensure Unicode output works on Windows terminals with narrow code pages.
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf-16"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Issue:
    file: str
    line: int
    code: str
    message: str

    def __str__(self) -> str:
        return f"  [{self.code}] line {self.line}: {self.message}"

    def key(self) -> tuple:
        return (self.file, self.line, self.code, self.message)


@dataclass
class APIInfo:
    """Public API surface detected from the main app."""
    # MainWindow members: name -> kind ("method", "property", "attr", "signal/classattr")
    mw_members: dict[str, str] = field(default_factory=dict)
    # Per-manager attr name -> set of its public member names
    manager_members: dict[str, set[str]] = field(default_factory=dict)
    # Manager attr name -> source class name (for readable error messages)
    manager_class_names: dict[str, str] = field(default_factory=dict)
    # PluginContext public members (used when --check-context is set)
    context_members: set[str] = field(default_factory=set)


# ---------------------------------------------------------------------------
# Qt inherited methods -- these are always valid on a QMainWindow/QWidget.
# Accesses to these on an MW variable are suppressed to avoid false positives.
# ---------------------------------------------------------------------------
_QT_INHERITED = frozenset({
    # QMainWindow
    "statusBar", "menuBar", "toolBar", "centralWidget", "setCentralWidget",
    "addToolBar", "addDockWidget", "removeDockWidget",
    # QWidget
    "show", "hide", "close", "resize", "move", "setWindowTitle", "setWindowIcon",
    "setMinimumSize", "setMaximumSize", "setFixedSize", "setSizePolicy",
    "update", "repaint", "raise_", "lower", "activateWindow",
    "setEnabled", "setDisabled", "isEnabled", "isVisible",
    "width", "height", "size", "pos", "geometry", "setGeometry",
    "parentWidget", "parent", "children", "findChild", "findChildren",
    "setAttribute", "testAttribute", "setStyleSheet",
    "setToolTip", "setStatusTip", "setWhatsThis",
    "installEventFilter", "removeEventFilter",
    "grabKeyboard", "releaseKeyboard",
    # QObject
    "connect", "disconnect", "blockSignals", "signalsBlocked",
    "deleteLater", "objectName", "setObjectName",
    "property", "setProperty", "metaObject",
    "thread", "moveToThread",
    # Common slots/signals that may be called directly
    "accept", "reject", "done",
})


# ---------------------------------------------------------------------------
# Phase 1 -- Main App API extraction
# ---------------------------------------------------------------------------

class AppAPIExtractor:
    """
    Walks the main app source tree and builds APIInfo by parsing:
      - main_window.py       -> MainWindow class members + manager attr names
      - Each manager module  -> all public methods + self.xxx assignments (all methods)
      - plugin_interface.py  -> PluginContext public members
    """

    # Class names that, when found in MainWindow.__init__ as `self.X = SomeClass(...)`,
    # indicate X is a manager whose internals plugins might call.
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

        ctx_file = self._find_file_containing("class PluginContext")
        if ctx_file:
            self._parse_plugin_context(ctx_file)

        if self.verbose:
            self._print_api_summary()

        return self.api

    # ------------------------------------------------------------------ #
    # MainWindow parsing
    # ------------------------------------------------------------------ #

    def _parse_main_window(self, path: Path):
        tree = self._parse(path)
        mw_class = self._find_class(tree, "MainWindow")
        if not mw_class:
            return

        for node in mw_class.body:
            # Class-level signals / plain assignments
            if isinstance(node, ast.Assign):
                for t in node.targets:
                    if isinstance(t, ast.Name):
                        self.api.mw_members.setdefault(t.id, "signal/classattr")

            # Methods (including __init__) and properties
            elif isinstance(node, ast.FunctionDef):
                kind = _method_kind(node)
                self.api.mw_members[node.name] = kind
                if node.name == "__init__":
                    self._mine_mw_init(node, path.parent)

    def _mine_mw_init(self, init_fn: ast.FunctionDef, search_dir: Path):
        """Extract self.X = SomeManager(...) assignments from MainWindow.__init__."""
        for stmt in ast.walk(init_fn):
            if not isinstance(stmt, ast.Assign):
                continue
            for target in stmt.targets:
                if not _is_self_attr(target):
                    continue
                attr_name = target.attr  # type: ignore[attr-defined]
                self.api.mw_members.setdefault(attr_name, "attr")

                # If RHS is a constructor call for a known manager class, load it.
                if isinstance(stmt.value, ast.Call):
                    class_name = _call_class_name(stmt.value)
                    if class_name and class_name in self.MANAGER_CLASS_HINTS:
                        self.api.mw_members[attr_name] = f"manager:{class_name}"
                        self.api.manager_class_names[attr_name] = class_name
                        self._load_manager_class(attr_name, class_name, search_dir)

    def _load_manager_class(self, attr_name: str, class_name: str, search_dir: Path):
        """Find the file that defines class_name and extract all its public members."""
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
                # Collect self.xxx from ALL methods (not just __init__),
                # because managers often set attributes in setup/init helpers.
                self._collect_self_attrs(node, members)

            elif isinstance(node, ast.Assign):
                # Class-level attributes / signals
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
        """Add every `self.xxx = ...` assignment found anywhere in func to out."""
        for stmt in ast.walk(func):
            if isinstance(stmt, ast.Assign):
                for t in stmt.targets:
                    if _is_self_attr(t):
                        out.add(t.attr)  # type: ignore[attr-defined]
            elif isinstance(stmt, ast.AnnAssign):
                if _is_self_attr(stmt.target):
                    out.add(stmt.target.attr)  # type: ignore[attr-defined]

    # ------------------------------------------------------------------ #
    # PluginContext parsing
    # ------------------------------------------------------------------ #

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

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #

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
        print(f"\n  Managers with parsed members:")
        for mgr, members in self.api.manager_members.items():
            cls = self.api.manager_class_names.get(mgr, "?")
            print(f"    mw.{mgr}  ({cls}, {len(members)} members)")
        if self.api.context_members:
            print(f"\n  PluginContext members: {len(self.api.context_members)}")
        print()


# ---------------------------------------------------------------------------
# Phase 2 -- Plugin scanning
# ---------------------------------------------------------------------------

# Variable names that strongly suggest the variable holds a MainWindow instance.
# Deliberately conservative: "window" is excluded because it's used for any QWidget.
_MW_VAR_NAMES = frozenset({
    "mw", "main_window", "parent_window", "app_window", "mainwindow",
})

# Attribute names on a call target that return a MainWindow.
_MW_GETTER_ATTRS = frozenset({"get_main_window"})

# Builtins that take an object as first arg but DON'T actually access its attributes
# -- suppress false positives for defensive checks like hasattr(mw, "x").
_SAFE_CALL_FUNCS = frozenset({"hasattr", "getattr", "setattr", "isinstance", "type"})

# Variable/parameter names that strongly suggest a PluginContext.
_CONTEXT_PARAM_NAMES = frozenset({"context", "ctx"})


class PluginFileChecker:
    """
    Two-pass AST visitor for a single plugin file.

    Pass 1 -- Alias collection
        Records which variable names / self.X attributes hold MainWindow or
        PluginContext references, based on assignments and function parameters.

    Pass 2 -- Access checking
        For every `X.attr` where X is an MW reference, checks `attr` against
        the known MainWindow API.  For chained `X.manager.attr`, checks the
        manager's member set.  Accesses inside hasattr/getattr are skipped.
    """

    def __init__(self, filepath: Path, api: APIInfo, check_context: bool = False):
        self.filepath = filepath
        self.api = api
        self.check_context = check_context
        self.issues: list[Issue] = []

        # Tracks names/self-attrs known to hold a MainWindow ref.
        # Simple names stored as-is; self.X attrs stored as "self.X".
        self._mw_refs: set[str] = set()

        # Tracks names/self-attrs known to hold a PluginContext ref.
        self._ctx_refs: set[str] = set()

        # De-dupe: (file, line, code, msg)
        self._seen: set[tuple] = set()

    # ------------------------------------------------------------------ #
    # Public entry point
    # ------------------------------------------------------------------ #

    def check(self) -> list[Issue]:
        try:
            source = self.filepath.read_text(encoding="utf-8", errors="ignore")
            tree = ast.parse(source, filename=str(self.filepath))
        except SyntaxError as exc:
            self.issues.append(Issue(str(self.filepath), 0, "SYNTAX_ERROR", str(exc)))
            return self.issues

        self._pass1_collect_aliases(tree)
        self._pass2_check_accesses(tree)
        return self.issues

    # ------------------------------------------------------------------ #
    # Pass 1 -- alias collection
    # ------------------------------------------------------------------ #

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
        # context.get_main_window() or self.context.get_main_window()
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr in _MW_GETTER_ATTRS:
                return True
        # Variable whose name strongly implies it's an MW
        if isinstance(node, ast.Name) and node.id in _MW_VAR_NAMES:
            return True
        # self.mw / self.main_window assigned from another mw-ref
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
            if node.value.id == "self" and node.attr in _MW_VAR_NAMES:
                return True
        return False

    def _rhs_is_ctx(self, node: ast.expr) -> bool:
        if isinstance(node, ast.Name) and node.id in _CONTEXT_PARAM_NAMES:
            return True
        return False

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

    # ------------------------------------------------------------------ #
    # Pass 2 -- access checking
    # ------------------------------------------------------------------ #

    def _pass2_check_accesses(self, tree: ast.Module):
        safe = _collect_safe_positions(tree)

        for node in ast.walk(tree):
            if not isinstance(node, ast.Attribute):
                continue
            if id(node) in safe:
                continue

            obj = node.value
            attr = node.attr

            # Skip private / dunder names
            if attr.startswith("_"):
                continue

            # Skip known Qt-inherited methods (valid on any QMainWindow)
            if attr in _QT_INHERITED:
                continue

            # ---- mw.something ------------------------------------------
            if self._is_mw_ref(obj):
                if attr not in self.api.mw_members:
                    self._add_issue(
                        node.lineno, "UNKNOWN_MW_ATTR",
                        f"`{_repr(obj)}.{attr}` -- '{attr}' not found on MainWindow"
                    )

            # ---- mw.manager_attr.something -----------------------------
            elif (
                isinstance(obj, ast.Attribute)
                and self._is_mw_ref(obj.value)
            ):
                mgr_attr = obj.attr
                if mgr_attr in self.api.manager_members:
                    if attr not in self.api.manager_members[mgr_attr]:
                        cls = self.api.manager_class_names.get(mgr_attr, mgr_attr)
                        self._add_issue(
                            node.lineno, "UNKNOWN_MANAGER_ATTR",
                            f"`{_repr(obj.value)}.{mgr_attr}.{attr}` -- "
                            f"'{attr}' not found in {cls}"
                        )

            # ---- context.something ------------------------------------
            elif self.check_context and self._is_ctx_ref(obj):
                if self.api.context_members and attr not in self.api.context_members:
                    self._add_issue(
                        node.lineno, "UNKNOWN_CONTEXT_ATTR",
                        f"`{_repr(obj)}.{attr}` -- '{attr}' not found on PluginContext"
                    )

    # ------------------------------------------------------------------ #
    # Reference detection
    # ------------------------------------------------------------------ #

    def _is_mw_ref(self, node: ast.expr) -> bool:
        if isinstance(node, ast.Name):
            return node.id in self._mw_refs or node.id in _MW_VAR_NAMES
        if _is_self_attr(node):
            return f"self.{node.attr}" in self._mw_refs  # type: ignore[attr-defined]
        # Inline: context.get_main_window()
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr in _MW_GETTER_ATTRS:
                return True
        return False

    def _is_ctx_ref(self, node: ast.expr) -> bool:
        if isinstance(node, ast.Name):
            return node.id in self._ctx_refs or node.id in _CONTEXT_PARAM_NAMES
        if _is_self_attr(node):
            return f"self.{node.attr}" in self._ctx_refs  # type: ignore[attr-defined]
        return False

    # ------------------------------------------------------------------ #
    # Issue recording (de-duplicated)
    # ------------------------------------------------------------------ #

    def _add_issue(self, line: int, code: str, message: str):
        issue = Issue(str(self.filepath), line, code, message)
        key = issue.key()
        if key not in self._seen:
            self._seen.add(key)
            self.issues.append(issue)


# ---------------------------------------------------------------------------
# Shared AST helpers
# ---------------------------------------------------------------------------

def _is_self_attr(node: ast.expr) -> bool:
    """True if node is `self.something`."""
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
    """Extract class name from SomeClass(...) or pkg.SomeClass(...)."""
    func = call.func
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return None


def _collect_safe_positions(tree: ast.Module) -> set[int]:
    """
    Return id()s of AST nodes that are the FIRST argument of hasattr/getattr/etc.
    These are defensive checks -- we skip reporting issues for them.
    """
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
            # Also suppress the inner attribute of e.g. hasattr(mw.view_3d_manager, "x")
            if isinstance(first, ast.Attribute):
                safe.add(id(first))
    return safe


def _repr(node: ast.expr) -> str:
    """Best-effort textual representation of an AST expression."""
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
    app_root = Path(args.app).expanduser().resolve()
    plugin_path = Path(args.plugin).expanduser().resolve()

    if not app_root.exists():
        _die(f"App path does not exist: {app_root}")
    if not plugin_path.exists():
        _die(f"Plugin path does not exist: {plugin_path}")

    print(f"App root  : {app_root}")
    print(f"Plugin(s) : {plugin_path}")
    print()

    # Phase 1 -- API extraction
    print("Extracting MainWindow API surface...")
    extractor = AppAPIExtractor(app_root, verbose=args.show_api)
    api = extractor.extract()
    n_mgrs = len(api.manager_members)
    n_mw = len(api.mw_members)
    print(f"  {n_mw} MainWindow members, {n_mgrs} managers parsed", end="")
    if api.context_members:
        print(f", {len(api.context_members)} PluginContext members", end="")
    print("\n")

    # Phase 2 -- Plugin scanning
    plugin_files = collect_plugin_files(plugin_path)
    print(f"Scanning {len(plugin_files)} plugin file(s)...\n")

    all_issues: list[tuple[Path, list[Issue]]] = []
    for pf in plugin_files:
        checker = PluginFileChecker(pf, api, check_context=args.check_context)
        issues = checker.check()
        if issues:
            all_issues.append((pf, issues))

    # Report
    if not all_issues:
        print("No issues found.")
        return 0

    total = sum(len(iss) for _, iss in all_issues)
    print(f"Found {total} issue(s) in {len(all_issues)} file(s):")
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
        description="Find API disconnections between MoleditPy plugins and the main app.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--app", required=True,
        help="Path to the main app root (e.g. ../python_molecular_editor/)",
    )
    parser.add_argument(
        "--plugin", required=True,
        help="Path to a plugin .py file or a directory of plugins (searched recursively)",
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
