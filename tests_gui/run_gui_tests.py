"""
Runner script for headless GUI tests.

Sets QT_QPA_PLATFORM=offscreen automatically so you don't need to export it
manually.  On Windows it also pins the PyQt6 DLL directory before launching
pytest so the correct Qt runtime is loaded.

Usage (from any directory):
    python tests_gui/run_gui_tests.py                              # all GUI tests
    python tests_gui/run_gui_tests.py test_gui_plugin_xyz_editor.py  # one file (fast)
    python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_xyz_editor.py::TestX
    python tests_gui/run_gui_tests.py -k xyz                       # filter by keyword
    python tests_gui/run_gui_tests.py -v -s                        # verbose + no capture

Passing a test file / node-id runs ONLY that target (bare filenames are
resolved against tests_gui/). With no target the whole suite runs. Output
is quiet by default — add ``-v`` for the per-test listing.
"""

from __future__ import annotations

import os
import sys

# ---------------------------------------------------------------------------
# 1. Set the Qt platform plugin before any Qt code is imported.
# ---------------------------------------------------------------------------
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# ---------------------------------------------------------------------------
# 2. On Windows, add PyQt6's own Qt6/bin to the DLL search path so the
#    correct Qt runtime is found when both PyQt6 and PySide6 are installed.
# ---------------------------------------------------------------------------
if sys.platform == "win32":
    try:
        import PyQt6 as _pyqt6_pkg  # noqa: PLC0415

        _qt6_bin = os.path.join(os.path.dirname(_pyqt6_pkg.__file__), "Qt6", "bin")
        if os.path.isdir(_qt6_bin):
            os.add_dll_directory(_qt6_bin)
    except Exception:
        pass

# ---------------------------------------------------------------------------
# 3. Run pytest programmatically, forwarding any extra CLI arguments.
# ---------------------------------------------------------------------------
import pytest  # noqa: E402

TESTS_GUI_DIR = os.path.dirname(os.path.abspath(__file__))


def _is_target(arg: str) -> bool:
    """True if *arg* names a test file / node-id (not an option)."""
    if arg.startswith("-"):
        return False
    base = arg.split("::", 1)[0]
    return base.endswith(".py") or "::" in arg or os.path.exists(base)


def _resolve_target(arg: str) -> str:
    """Let bare filenames/node-ids resolve against tests_gui/ from any cwd."""
    base = arg.split("::", 1)[0]
    if os.path.exists(base):
        return arg
    candidate = os.path.join(TESTS_GUI_DIR, base)
    if os.path.exists(candidate):
        return candidate + arg[len(base):]  # keep any ::node-id suffix
    return arg


if __name__ == "__main__":
    extra_args = sys.argv[1:]
    # Only collect the whole directory when the caller named no explicit
    # target — otherwise pytest would union the dir with the target and run
    # everything (the historical "why is one file so slow" bug).
    if any(_is_target(a) for a in extra_args):
        args = [_resolve_target(a) if _is_target(a) else a for a in extra_args]
    else:
        args = [TESTS_GUI_DIR, *extra_args]
    raise SystemExit(pytest.main([*args, "--tb=short"]))
