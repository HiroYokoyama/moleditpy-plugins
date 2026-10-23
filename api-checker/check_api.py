"""Registry-aware wrapper around the shared API checker engine."""

import argparse
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent
_DEFAULT_PLUGINS = (_REPO_ROOT / "plugins").resolve()
_DEFAULT_REGISTRY = (_REPO_ROOT / "REGISTRY/plugins.json").resolve()
_spec = importlib.util.spec_from_file_location(
    "_common_api_engine", _SCRIPT_DIR / "plugin_api_checker.py"
)
_engine = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_engine)
Issue = _engine.Issue
AppAPIExtractor = _engine.AppAPIExtractor
PluginFileChecker = _engine.PluginFileChecker
collect_plugin_files = _engine.collect_plugin_files
_merge_allowlists = _engine._merge_allowlists
_load_site_allowlist = _engine._load_site_allowlist
_MANAGER_ALLOWLIST = _engine._MANAGER_ALLOWLIST
_MW_ALLOWLIST = _engine._MW_ALLOWLIST


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


def run(args) -> int:
    app_root = Path(args.app).expanduser().resolve()
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
        + (
            f", {len(api.context_members)} PluginContext members"
            if api.context_members
            else ""
        )
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
        n_s_mw = len(site_allowlist.get("mw", set()))
        n_s_mgr = sum(len(v) for v in site_allowlist.get("manager", {}).values())
        site_desc = []
        if n_s_mw:
            site_desc.append(f"{n_s_mw} mw attr(s)")
        if n_s_mgr:
            site_desc.append(f"{n_s_mgr} manager attr(s)")
        print(f"Site allowlist (.moleditpy-api-allowlist): {', '.join(site_desc)}\n")
        parts.append(site_allowlist)

    allowlist = _merge_allowlists(*parts) if parts else {}
    if allowlist:
        n_mgr = sum(len(v) for v in allowlist.get("manager", {}).values())
        n_mw = len(allowlist.get("mw", set()))
        suppressed = []
        if n_mgr:
            suppressed.append(f"{n_mgr} manager attr(s)")
        if n_mw:
            suppressed.append(f"{n_mw} mw attr(s)")
        print(f"Allowlist active: {', '.join(suppressed)} suppressed\n")

    all_issues: list[tuple[Path, list[Issue]]] = []
    for pf in plugin_files:
        checker = PluginFileChecker(
            pf, api, check_context=args.check_context, allowlist=allowlist
        )
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

    suffix = (
        f"  ({n_try} inside try blocks, shown with [try])"
        if n_try and not args.skip_try
        else ""
    )
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
            "UNKNOWN_MW_ATTR": "Unknown MainWindow attribute",
            "UNKNOWN_MANAGER_ATTR": "Unknown manager attribute",
            "UNKNOWN_CONTEXT_ATTR": "Unknown PluginContext attribute",
            "SYNTAX_ERROR": "Syntax error in plugin file",
        }.get(code, code)
        print(f"  {label:<40} {count}")

    return 1


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
        "--default-allowlist",
        action="store_true",
        help=(
            "Suppress known manager false positives: runtime attrs assigned via "
            "self.host.manager.X patterns invisible to AST "
            "(init_manager.scene, state_manager.data, etc.)"
        ),
    )
    parser.add_argument(
        "--mw-allowlist",
        action="store_true",
        help=(
            "Also suppress known direct mw.X legacy compat attrs "
            "(mw.host, mw.view3d, etc.).  Off by default so mw.X issues remain visible."
        ),
    )
    parser.add_argument(
        "--skip-try",
        action="store_true",
        help=(
            "Hide issues whose access is inside a try: block. "
            "By default all issues are shown; those inside try blocks are tagged [try]."
        ),
    )
    parser.add_argument(
        "--check-context",
        action="store_true",
        help="Also check context.xxx accesses against the PluginContext API",
    )
    parser.add_argument(
        "--show-api",
        action="store_true",
        help="Print the detected MainWindow API surface before scanning",
    )
    args = parser.parse_args()
    sys.exit(run(args))


if __name__ == "__main__":
    main()
