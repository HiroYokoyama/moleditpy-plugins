"""Esc must not skip a dialog's cleanup.

A QDialog leaves through ``reject() -> done()`` and raises no QCloseEvent, so
anything kept in ``closeEvent`` alone is skipped entirely when the window is
dismissed with Esc -- or by a Close button wired to ``reject``, which is the
usual spelling.

That went wrong in eight plugins at once, in every flavour the mistake has:

* the colorizers force the 3D viewer into select mode and restore it on close,
  so Esc left the viewer stuck and unusable, with the picking timer still
  polling and the window still registered -- which meant reopening the plugin
  handed back the same dead dialog;
* the xTB optimizer asks before abandoning a running optimization, so Esc
  closed the window with the worker thread still running;
* the chat windows disconnect a streaming reply, so Esc left the worker
  emitting chunks into a window that had gone;
* the plugin installer stops its fetch thread, the animated player stops two
  timers and hands the displayed frame back, the MS spectrum window stops a
  live-sync timer, and the advanced-rendering window saves its settings.

The rule is therefore structural, and cheap to keep: a QDialog that defines
``closeEvent`` must also define ``done`` (or ``reject``), so its cleanup sits
somewhere every close route passes through. Overriding ``done`` is the better
of the two -- ``accept()`` and ``reject()`` both call it.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"


def _dialog_classes():
    """(plugin, file, class, method names) for every QDialog subclass."""
    found = []
    for path in sorted(PLUGINS_DIR.rglob("*.py")):
        if "_old" in path.parts:
            # Retired plugins are frozen; never reported on.
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            if not any("QDialog" in ast.unparse(base) for base in node.bases):
                continue
            methods = {n.name for n in node.body if isinstance(n, ast.FunctionDef)}
            found.append((path.relative_to(PLUGINS_DIR).parts[0], path.name, node.name, methods))
    return found


DIALOGS = _dialog_classes()


def test_the_scan_finds_dialogs_at_all():
    # A rename of the base class, or a broken parse, would otherwise turn every
    # test below into a vacuous pass.
    assert len(DIALOGS) >= 10


@pytest.mark.parametrize(
    "plugin,filename,cls,methods",
    DIALOGS,
    ids=[f"{p}:{c}" for p, _f, c, _m in DIALOGS],
)
def test_cleanup_is_reachable_from_every_close_route(plugin, filename, cls, methods):
    if "closeEvent" not in methods:
        pytest.skip("no close-time cleanup to reach")
    assert methods & {"done", "reject"}, (
        f"{plugin}/{filename}:{cls} cleans up in closeEvent but overrides neither "
        "done() nor reject(). Esc leaves a QDialog through reject() -> done() "
        "and raises no QCloseEvent, so that cleanup would be skipped entirely."
    )
