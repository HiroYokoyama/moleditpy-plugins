"""
Static API-compatibility check for all visible plugins.

Uses ``scripts/check_api.py`` to verify that every visible plugin only accesses
MainWindow / manager attributes that actually exist in the main app.

Both allowlists are active:

- ``--default-allowlist``: suppresses manager attrs set via
  ``self.host.manager.X = ...`` patterns that are AST-invisible.
- ``--mw-allowlist``: suppresses V2 legacy bridge attrs (``mw.host``,
  ``mw.view3d``) that are always guarded by ``hasattr`` and are intentional.

With both allowlists the expected issue count is zero.  Any new plugin that
calls a genuinely non-existent attribute fails the test.

The test is automatically skipped if the main app is not present at the
expected sibling path (``../python_molecular_editor``).
"""

from __future__ import annotations

import importlib.util
import types
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
_REGISTRY_PATH = Path(__file__).resolve().parents[1] / "REGISTRY" / "plugins.json"
_DEFAULT_APP = Path(__file__).resolve().parents[2] / "python_molecular_editor"

# ---------------------------------------------------------------------------
# Load check_api as a module (it lives in scripts/, not a package)
# ---------------------------------------------------------------------------


def _load_check_api():
    script_path = _SCRIPTS_DIR / "check_api.py"
    spec = importlib.util.spec_from_file_location("check_api", script_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    mod.__file__ = str(script_path)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


# ---------------------------------------------------------------------------
# Skip condition
# ---------------------------------------------------------------------------

_APP_MISSING = not (_DEFAULT_APP / "moleditpy").exists()


@pytest.mark.skipif(
    _APP_MISSING,
    reason=f"Main app not found at {_DEFAULT_APP} — skipping API check",
)
def test_no_unknown_api_accesses() -> None:
    """
    All visible plugins must only access MainWindow / manager attributes that
    exist in the main app.

    Both ``_MANAGER_ALLOWLIST`` and ``_MW_ALLOWLIST`` are active, so known
    false positives (AST-invisible manager attrs, V2 compat bridge attrs) are
    suppressed.  Any genuinely unknown attribute causes the test to fail.
    """
    check_api = _load_check_api()

    # Build the APIInfo from the real main app
    extractor = check_api.AppAPIExtractor(_DEFAULT_APP, verbose=False)
    api = extractor.extract()

    # Load visible plugin files from the registry
    plugin_files = check_api.load_visible_plugin_files(
        _REGISTRY_PATH,
        _DEFAULT_APP.parent / "moleditpy-plugins" / "plugins",
    )

    # Both allowlists active — same as running with --default-allowlist --mw-allowlist
    allowlist = check_api._merge_allowlists(  # noqa: SLF001
        check_api._MANAGER_ALLOWLIST,          # noqa: SLF001
        check_api._MW_ALLOWLIST,               # noqa: SLF001
    )

    all_issues = []
    for pf in plugin_files:
        checker = check_api.PluginFileChecker(
            pf, api, check_context=False, allowlist=allowlist
        )
        all_issues.extend(checker.check())

    if all_issues:
        lines = [
            f"  [{i.code}] {Path(i.file).name} line {i.line}: {i.message}"
            for i in all_issues
        ]
        pytest.fail(
            f"{len(all_issues)} unknown API access(es) found in visible plugins:\n"
            + "\n".join(lines)
        )
