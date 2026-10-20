"""Exercise the actual Explorer renderer, including untrusted license strings."""

from pathlib import Path
import shutil
import subprocess

import pytest


def test_explorer_license_rendering():
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is needed to run the Explorer JavaScript renderer")
    root = Path(__file__).resolve().parents[1]
    code = r"""
const assert = require('node:assert/strict');
const {renderLicenseMetadata: render} = require('./explorer/license.js');
assert.match(render('LGPL-3.0-or-later'), /LGPL-3.0-or-later/);
for (const absent of [undefined, null, '', 'NOASSERTION', 42]) {
    assert.match(render(absent), /Not specified/);
}
const html = render('<img src=x onerror="alert(1)">');
assert.ok(!html.includes('<img'));
assert.match(html, /&lt;img/);
assert.match(html, /&quot;/);
assert.match(html, /aria-hidden="true"/);
assert.match(html, />\s*License\s*</);
"""
    subprocess.run([node, "-e", code], cwd=root, capture_output=True, text=True, check=True)
