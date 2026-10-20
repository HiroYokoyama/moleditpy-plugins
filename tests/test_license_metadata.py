"""License precedence, published asset parsing and registry-only maintenance."""

import io
import json
from pathlib import Path
import sys
import zipfile

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import license_metadata as licenses
import register_remote_plugin as registration
import update_intra_repo_metadata as intra
import update_plugin_metadata as metadata
from update_registry_licenses import update_licenses


@pytest.mark.parametrize("code,manual,stored,expected", [
    ("MIT", "BSD-3-Clause", "GPL-3.0-only", "MIT"),
    (None, "BSD-3-Clause", "GPL-3.0-only", "BSD-3-Clause"),
    (None, None, "GPL-3.0-only", "GPL-3.0-only"),
    (" ", "", "LGPL-3.0-or-later", "LGPL-3.0-or-later"),
    (None, None, None, "NOASSERTION"),
    ("(MIT OR Apache-2.0) AND BSD-3-Clause", None, None,
     "(MIT OR Apache-2.0) AND BSD-3-Clause"),
])
def test_precedence(code, manual, stored, expected):
    assert licenses.resolve_license(code, manual, stored) == expected


@pytest.mark.parametrize("value", [[], 42, "MIT\nOR GPL", "<script>", "x" * 257])
def test_rejects_invalid_license_values(value):
    with pytest.raises(ValueError):
        licenses.normalize_license(value)


def test_reads_license_from_python_and_zip_without_importing(tmp_path):
    body = 'PLUGIN_NAME = "Demo"\nPLUGIN_VERSION = "1.1.0"\nPLUGIN_LICENSE = "MIT"\nraise RuntimeError("must not run")\n'
    source = tmp_path / "demo.py"
    source.write_text(body, encoding="utf-8")
    assert registration.extract_metadata_from_file(source)["license"] == "MIT"
    archive = tmp_path / "demo.zip"
    with zipfile.ZipFile(archive, "w") as z:
        z.writestr("demo/__init__.py", body)
    assert registration.extract_metadata_from_file(archive)["license"] == "MIT"
    assert registration.extract_metadata_from_code_regex(body)["license"] == "MIT"


@pytest.fixture
def registry(tmp_path, monkeypatch):
    folder = tmp_path / "REGISTRY"
    folder.mkdir()
    path = folder / "plugins.json"
    entry = {"id": "demo", "visible": True, "name": "Demo", "version": "1.0.0",
             "supported_moleditpy_version": ">=5.0.0, <6.0.0", "license": "GPL-3.0-only",
             "sha256": "original-hash", "downloadUrl": "https://github.com/HiroYokoyama/demo/releases/download/v1.0.0/demo.py",
             "projectUrl": "https://github.com/HiroYokoyama/demo"}
    path.write_text(json.dumps([entry]), encoding="utf-8")
    for module in (registration, metadata):
        monkeypatch.setattr(module, "__file__", str(tmp_path / "scripts" / "script.py"))
    monkeypatch.delenv("GITHUB_OUTPUT", raising=False)
    return path


@pytest.mark.parametrize("declared,manual,expected", [
    ('PLUGIN_LICENSE = "MIT"', "BSD-3-Clause", "MIT"),
    ("", "BSD-3-Clause", "BSD-3-Clause"),
    ("", None, "GPL-3.0-only"),
])
@pytest.mark.parametrize("mode", ["release", "metadata"])
def test_remote_workflows_preserve_or_update_license(registry, monkeypatch, declared, manual, expected, mode):
    body = ('PLUGIN_NAME = "Demo"\nPLUGIN_VERSION = "1.1.0"\n'
            'PLUGIN_AUTHOR = "HiroYokoyama"\nPLUGIN_DESCRIPTION = "Demo plugin"\n' + declared)
    monkeypatch.setattr(registration.urllib.request, "urlopen", lambda *a, **kw: io.BytesIO(body.encode()))
    url = "https://github.com/HiroYokoyama/demo/releases/download/v1.1.0/demo.py"
    if mode == "release":
        args = ["register_remote_plugin.py", url]
        if manual:
            args += ["--license", manual]
        monkeypatch.setattr(sys, "argv", args)
        registration.main()
    else:
        metadata.update_metadata_only(url, plugin_license=manual)
    entry = json.loads(registry.read_text(encoding="utf-8"))[0]
    assert entry["license"] == expected
    if mode == "metadata":
        assert entry["version"] == "1.0.0"
        assert entry["sha256"] == "original-hash"
        assert entry["downloadUrl"].endswith("v1.0.0/demo.py")


def test_new_registration_has_explicit_unknown_license(registry, monkeypatch):
    registry.write_text("[]", encoding="utf-8")
    body = (b'PLUGIN_NAME = "Demo"\nPLUGIN_VERSION = "1.0.0"\nPLUGIN_AUTHOR = "HiroYokoyama"\n'
            b'PLUGIN_DESCRIPTION = "Demo"\nPLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=5.0.0"\n')
    monkeypatch.setattr(registration.urllib.request, "urlopen", lambda *a, **kw: io.BytesIO(body))
    monkeypatch.setattr(sys, "argv", ["register_remote_plugin.py", "https://github.com/HiroYokoyama/demo/releases/download/v1.0.0/demo.py"])
    registration.main()
    assert json.loads(registry.read_text())[0]["license"] == "NOASSERTION"


@pytest.mark.parametrize("declaration,expected", [("", "GPL-3.0-only"), ('PLUGIN_LICENSE = "MIT"', "MIT")])
def test_intra_repo_sync_keeps_registry_only_license(registry, declaration, expected):
    data = json.loads(registry.read_text())
    data[0]["downloadUrl"] = "../plugins/Demo/demo.py"
    registry.write_text(json.dumps(data), encoding="utf-8")
    source = registry.parents[1] / "plugins/Demo/demo.py"
    source.parent.mkdir(parents=True)
    source.write_text('PLUGIN_VERSION = "1.0.0"\n' + declaration, encoding="utf-8")
    intra.update_single_json(registry)
    assert json.loads(registry.read_text())[0]["license"] == expected
    assert source.read_text() == 'PLUGIN_VERSION = "1.0.0"\n' + declaration


def test_registry_only_update_preserves_release_identity(registry):
    before = json.loads(registry.read_text())[0]
    update_licenses(registry, ["demo=MIT"])
    after = json.loads(registry.read_text())[0]
    assert after == dict(before, license="MIT")


def test_registry_only_update_is_atomic_on_invalid_id(registry):
    original = registry.read_bytes()
    with pytest.raises(ValueError):
        update_licenses(registry, ["demo=MIT", "missing=MIT"])
    assert registry.read_bytes() == original
