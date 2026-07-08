"""
Tests for the Vector Viewer plugin.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
VECTOR_VIEWER_PATH = PLUGINS_DIR / "Vector_Viewer" / "vector_viewer.py"


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
