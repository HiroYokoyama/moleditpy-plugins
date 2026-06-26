"""
Runner script for headless GUI tests.

Sets QT_QPA_PLATFORM=offscreen automatically so you don't need to export it
manually.  On Windows it also pins the PyQt6 DLL directory before launching
pytest so the correct Qt runtime is loaded.

Usage (from any directory):
    python tests_gui/run_gui_tests.py           # all GUI tests
    python tests_gui/run_gui_tests.py -k xyz    # filter by keyword
    python tests_gui/run_gui_tests.py -v -s     # verbose + no capture
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

if __name__ == "__main__":
    extra_args = sys.argv[1:]
    raise SystemExit(pytest.main([TESTS_GUI_DIR, "-v", "--tb=short", *extra_args]))
