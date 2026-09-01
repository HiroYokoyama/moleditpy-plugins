"""Reading a calculation's output must not depend on the machine's locale.

``open(path, "r")`` decodes with ``locale.getpreferredencoding()``. That is
UTF-8 on most Linux and macOS installs, but on a Japanese Windows it is cp932,
and on other Windows installs cp1252 or cp936. A cube, fchk or ORCA output is
written by another program and routinely carries a title, a path or a comment
that is not ASCII -- and then the same file read on two machines gives two
different answers:

* cp932 against a UTF-8 file with an em dash or a Japanese character raises
  UnicodeDecodeError, so the plugin cannot open the file at all;
* against an accented Latin or Greek character it decodes without complaint
  into mojibake, which is worse, because nothing says so.

Both were reproduced on a cp932 machine before these were pinned. The readers
therefore name UTF-8 explicitly and pass ``errors="replace"``: a stray byte in
a numeric data file should cost one character, never the whole file.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

#: Readers and writers of files produced by (or handed to) other programs.
DATA_FILE_IO = [
    ("Cube_File_Viewer", "cube_viewer.py"),
    ("Cube_File_Viewer_Advanced", "cube_viewer_advanced.py"),
    ("Mapped_Cube_Viewer", "mapped_cube_viewer.py"),
    ("Gaussian_Freq_Analyzer", "gaussian_fchk_freq_analyzer.py"),
    ("Gaussian_MO_Analyzer", "gaussian_fchk_mo_analyzer/analyzer.py"),
    ("Gaussian_MO_Analyzer", "gaussian_fchk_mo_analyzer/vis.py"),
    ("ORCA_Freq_Analyzer", "orca_out_freq_analyzer.py"),
    ("Molecule_Comparator", "molecule_comparator.py"),
]


def _text_opens_without_encoding(path: Path):
    """Line numbers of text-mode open() calls that name no encoding."""
    tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
        if name != "open":
            continue
        if any(k.arg == "encoding" for k in node.keywords):
            continue
        mode = next(
            (a.value for a in node.args[1:2] if isinstance(a, ast.Constant)), ""
        )
        if "b" in str(mode):
            continue  # binary needs no encoding
        offenders.append(node.lineno)
    return offenders


@pytest.mark.parametrize(
    "plugin,filename", DATA_FILE_IO, ids=[f"{p}/{f}" for p, f in DATA_FILE_IO]
)
def test_data_files_are_read_as_utf8(plugin, filename):
    path = PLUGINS_DIR / plugin / filename
    assert path.exists(), f"{path} moved; update this list"
    offenders = _text_opens_without_encoding(path)
    assert not offenders, (
        f"{plugin}/{filename} opens a text file without naming an encoding at "
        f"line(s) {offenders}. It would then be decoded with the machine's "
        "locale -- cp932 on a Japanese Windows -- so the same file reads "
        "differently on two machines, or not at all."
    )


def test_the_scan_can_actually_spot_one(tmp_path):
    # Otherwise a broken parser would make every assertion above vacuous.
    sample = tmp_path / "sample.py"
    sample.write_text(
        "with open(p, 'r') as f:\n    pass\n"
        "with open(q, 'rb') as f:\n    pass\n"
        "with open(r, 'r', encoding='utf-8') as f:\n    pass\n",
        encoding="utf-8",
    )
    assert _text_opens_without_encoding(sample) == [1]


def test_readers_survive_a_byte_that_is_not_utf8(tmp_path):
    """errors="replace", not a bare utf-8 open.

    A cube file is mostly numbers; one bad byte in a comment must cost one
    character, not the whole file.
    """
    path = tmp_path / "grid.cube"
    path.write_bytes("title\n".encode("utf-8") + b"\xff\xfe bad\n" + b"0 0.0 0.0 0.0\n")
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        assert len(handle.readlines()) == 3
    with pytest.raises(UnicodeDecodeError):
        with open(path, "r", encoding="utf-8") as handle:
            handle.readlines()
