"""
Scan plugins/ for silent error patterns.

Checks:
  1. Bare except: clauses (catches everything incl. KeyboardInterrupt/SystemExit)
  2. hasattr(subject, ...) guards that may silently no-op
  3. Bare hasattr(...) expressions whose result is discarded
  4. .get(key) without a default value (returns None silently)

Usage:
  python scripts/scan_code_issues.py
      Run with defaults (hasattr subjects: self, mw, main_window)

  python scripts/scan_code_issues.py --hasattr-subjects=self
      Only scan hasattr(self, ...) patterns

  python scripts/scan_code_issues.py --hasattr-subjects=self,mw,main_window,context
      Add extra subjects to the hasattr scan
"""

import ast
import sys
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGINS_DIR = ROOT / "plugins"


# ---------------------------------------------------------------------------
# 1. Bare except: (text-based)
# ---------------------------------------------------------------------------

def scan_bare_except(path: Path) -> list[tuple[int, str]]:
    hits = []
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return hits
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if "script_lines.append" in line and "except" in line:
            continue
        if stripped.startswith("except:") and "except Exception" not in stripped:
            hits.append((i, line.rstrip()))
    return hits


# ---------------------------------------------------------------------------
# 2. hasattr(self/mw/main_window, ...) guards (AST-based)
# ---------------------------------------------------------------------------

# Default subjects — override via --hasattr-subjects=self,mw,main_window
DEFAULT_HASATTR_SUBJECTS = {"self", "mw", "main_window"}

def scan_self_hasattr(path: Path, subjects: set[str] | None = None) -> list[tuple[int, str, str, str]]:
    if subjects is None:
        subjects = DEFAULT_HASATTR_SUBJECTS
    hits = []
    try:
        source = path.read_text(encoding="utf-8", errors="ignore")
        tree = ast.parse(source, filename=str(path))
    except (OSError, SyntaxError):
        return hits
    lines = source.splitlines()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Name) and func.id == "hasattr"):
            continue
        if len(node.args) < 2:
            continue
        first, second = node.args[0], node.args[1]
        if not (isinstance(first, ast.Name) and first.id in subjects):
            continue
        if not isinstance(second, ast.Constant) or not isinstance(second.value, str):
            continue
        attr = second.value
        if attr.startswith("_"):
            continue
        subject = first.id
        hits.append((node.lineno, subject, attr, lines[node.lineno - 1].rstrip()))
    return hits


# ---------------------------------------------------------------------------
# 3. Bare hasattr(...) expressions — result discarded (AST-based)
# ---------------------------------------------------------------------------

def scan_bare_hasattr_expr(path: Path) -> list[tuple[int, str]]:
    """hasattr() called as a standalone statement — return value thrown away."""
    hits = []
    try:
        source = path.read_text(encoding="utf-8", errors="ignore")
        tree = ast.parse(source, filename=str(path))
    except (OSError, SyntaxError):
        return hits
    lines = source.splitlines()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Expr):
            continue
        val = node.value
        if not isinstance(val, ast.Call):
            continue
        func = val.func
        if isinstance(func, ast.Name) and func.id == "hasattr":
            hits.append((node.lineno, lines[node.lineno - 1].rstrip()))
    return hits


# ---------------------------------------------------------------------------
# 4. .get(key) without default — returns None silently (AST-based)
# ---------------------------------------------------------------------------

def scan_get_no_default(path: Path) -> list[tuple[int, str]]:
    """
    Find dict.get(key) calls with only one argument (no default).
    These return None silently when the key is absent.
    Only flags cases where the result is actually used (not bare expressions).
    """
    hits = []
    try:
        source = path.read_text(encoding="utf-8", errors="ignore")
        tree = ast.parse(source, filename=str(path))
    except (OSError, SyntaxError):
        return hits
    lines = source.splitlines()

    # Collect all .get() calls that are bare expressions (result discarded)
    bare_get_lines = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            call = node.value
            if (isinstance(call.func, ast.Attribute) and
                    call.func.attr == "get" and
                    len(call.args) == 1 and not call.keywords):
                bare_get_lines.add(node.lineno)

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == "get"):
            continue
        # Only 1 positional arg, no keywords => no default provided
        if len(node.args) != 1 or node.keywords:
            continue
        # Skip bare expressions (result unused entirely — different issue)
        if node.lineno in bare_get_lines:
            continue
            
        # Apply the same intentional-skip rules used by auto_fix_issues:
        line_text = lines[node.lineno - 1]
        
        # 1. Skip if inside an f-string
        if re.search(r'f["\'].*\.get\(', line_text):
            continue
            
        # We need the prefix to check for context
        # In AST we only have column offsets, but since we deduplicate by line,
        # we can just use simple regex checks on the line string itself.
        # Check if any .get is part of a dict/list key subscript:
        left_bracket_count = line_text.split('.get(')[0].count('[')
        right_bracket_count = line_text.split('.get(')[0].count(']')
        if left_bracket_count > right_bracket_count:
            continue
        
        # Check if wrapped in int()/float()/str()/bool()
        if re.search(r'\b(?:int|float|str|bool)\s*\([^)]*\.get\(', line_text):
            continue
            
        # Check if chained (e.g. .get().something)
        if '.get(' in line_text and ').' in line_text.split('.get(')[-1]:
            # Basic approximation. AST is cleaner:
            parent = getattr(node, "parent", None) # AST doesn't have parent out of box
            pass # We'll just rely on the regex if we really needed to, but let's do a strict AST check below:
        
        # Simpler AST check for chaining: if this node is inside an Attribute access
        is_chained = False
        for pnode in ast.walk(tree):
            if isinstance(pnode, ast.Attribute) and pnode.value is node:
                is_chained = True
        if is_chained:
            continue

        hits.append((node.lineno, line_text.rstrip()))

    # Deduplicate by line (multiple .get() on same line)
    seen = set()
    deduped = []
    for lineno, code in hits:
        if lineno not in seen:
            seen.add(lineno)
            deduped.append((lineno, code))
    return sorted(deduped)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def report_section(title: str, results: dict[str, list], fmt):
    total = sum(len(v) for v in results.values())
    print("=" * 70)
    print(f"{title} -- {total} occurrences in {len(results)} files")
    print("=" * 70)
    for rel, hits in results.items():
        print(f"\n  {rel}")
        for hit in hits:
            print(f"    L{hit[0]:>5}: {fmt(hit)}")
    print()


def main():
    # Parse --hasattr-subjects=a,b,c  (optional override)
    subjects = DEFAULT_HASATTR_SUBJECTS
    for arg in sys.argv[1:]:
        if arg.startswith("--hasattr-subjects="):
            subjects = set(arg.split("=", 1)[1].split(","))

    registry_path = ROOT / "REGISTRY" / "plugins.json"
    visible_paths = set()
    if registry_path.exists():
        try:
            reg_data = json.loads(registry_path.read_text(encoding="utf-8-sig"))
            for p in reg_data:
                if p.get("visible", False):
                    # the url is like "../plugins/Name/plugin.py" or ".../plugin.zip"
                    url = p.get("downloadUrl", "")
                    if "/plugins/" in url:
                        folder = url.split("/plugins/")[-1].split("/")[0].replace(".zip", "")
                        visible_paths.add(folder.lower())
        except Exception:
            pass

    py_files = sorted(PLUGINS_DIR.rglob("*.py"))

    bare_except   = {}
    self_hasattr  = {}
    bare_hasattr  = {}
    get_no_default = {}

    for path in py_files:
        if "__pycache__" in path.parts:
            continue
            
        # Only scan if part of a visible plugin (or root / core files)
        rel_parts = path.relative_to(PLUGINS_DIR).parts
        if rel_parts and visible_paths:
            plugin_folder = rel_parts[0].lower()
            if plugin_folder not in visible_paths:
                continue

        rel = str(path.relative_to(ROOT))

        r = scan_bare_except(path)
        if r: bare_except[rel] = r

        r = scan_self_hasattr(path, subjects)
        if r: self_hasattr[rel] = r

        r = scan_bare_hasattr_expr(path)
        if r: bare_hasattr[rel] = r

        r = scan_get_no_default(path)
        if r: get_no_default[rel] = r

    report_section("1. BARE except:",       bare_except,    lambda h: h[1])
    subj_label = "/".join(sorted(subjects))
    report_section(f"2. hasattr({subj_label}, ...)",  self_hasattr,   lambda h: f"[{h[1]}.{h[2]}]  {h[3]}")
    report_section("3. BARE hasattr() expr (result discarded)", bare_hasattr, lambda h: h[1])
    report_section("4. .get(key) no default (silent None)", get_no_default, lambda h: h[1])

    total = sum(sum(len(v) for v in d.values()) for d in [bare_except, self_hasattr, bare_hasattr, get_no_default])
    print(f"Total issues found: {total}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
