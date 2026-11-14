"""Submit to Cluster in the single-file input generators.

The button hands the saved input to the Job Manager plugin's public
``submit_file`` entry point, the same handoff the Pro generators use.
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import extract_function, load_plugin, mock_optional_imports

PLUGINS = Path(__file__).resolve().parent.parent / "plugins"

GENERATORS = [
    (PLUGINS / "Mopac_Input_Generator" / "mopac_input_generator.py", "MopacSetupDialog", ".mop"),
    (PLUGINS / "NWChem_Input_Generator" / "nwchem_input_generator.py", "NwchemSetupDialog", ".nw"),
    (PLUGINS / "Psi4_Input_Generator" / "psi4_input_generator.py", "Psi4SetupDialog", ".dat"),
    (PLUGINS / "PySCF_Input_Generator" / "pyscf_input_generator.py", "PyscfSetupDialog", ".py"),
]
IDS = [p.stem for p, _, _ in GENERATORS]


@pytest.fixture(params=GENERATORS, ids=IDS)
def generator(request):
    path, cls, ext = request.param
    return SimpleNamespace(path=path, cls=cls, ext=ext)


def _main_window(*records):
    return SimpleNamespace(plugin_manager=SimpleNamespace(plugins=list(records)))


class TestFindJobManager:
    def _find(self, generator):
        with mock_optional_imports():
            mod = load_plugin(generator.path)
        return mod.find_job_manager

    def test_found_by_name_and_submit_file(self, generator):
        job_manager = SimpleNamespace(submit_file=lambda *a, **k: True)
        mw = _main_window(
            {"name": "Other", "module": SimpleNamespace(submit_file=lambda *a: True)},
            {"name": "Job Manager", "module": job_manager},
        )
        assert self._find(generator)(mw) is job_manager

    def test_an_old_job_manager_without_submit_file_is_absent(self, generator):
        mw = _main_window({"name": "Job Manager", "module": SimpleNamespace()})
        assert self._find(generator)(mw) is None

    @pytest.mark.parametrize(
        "mw",
        [None, SimpleNamespace(), _main_window(), _main_window("not a record")],
    )
    def test_no_plugin_list_means_absent(self, generator, mw):
        assert self._find(generator)(mw) is None


class _Preview:
    def __init__(self, text):
        self.text = text

    def toPlainText(self):
        return self.text


def _dialog(generator, tmp_path, job_manager, save_to="job", cancel=False):
    """A fake dialog whose save_file writes the preview like the real one."""
    funcs = {
        name: extract_function(
            generator.path,
            generator.cls,
            name,
            {
                "os": os,
                "logging": logging,
                "PLUGIN_NAME": "plugin",
                "QMessageBox": dlg_box,
                "find_job_manager": lambda mw: job_manager,
            },
        )
        for name in ("submit_to_cluster", "_saved_is_current")
    }
    dlg = SimpleNamespace(
        saved_path=None,
        _saved_preview=None,
        _quiet_save=False,
        preview_text=_Preview("input text"),
        parent=lambda: "main window",
        accept=MagicMock(),
        quiet_seen=[],
        saves=0,
    )

    def save_file():
        dlg.saves += 1
        dlg.quiet_seen.append(dlg._quiet_save)
        if cancel:
            return
        path = str(tmp_path / (save_to + generator.ext))
        Path(path).write_text(dlg.preview_text.toPlainText(), encoding="utf-8")
        dlg.saved_path = path
        dlg._saved_preview = dlg.preview_text.toPlainText()

    dlg.save_file = save_file
    for name, fn in funcs.items():
        setattr(dlg, name, fn.__get__(dlg))
    return dlg


dlg_box = MagicMock()


@pytest.fixture(autouse=True)
def _reset_box():
    dlg_box.reset_mock()


class TestSubmitToCluster:
    def test_first_submit_saves_quietly_then_hands_over(self, generator, tmp_path):
        jm = MagicMock()
        jm.submit_file.return_value = True
        dlg = _dialog(generator, tmp_path, jm, save_to="benzene")
        dlg.submit_to_cluster()
        assert dlg.saves == 1 and dlg.quiet_seen == [True] and dlg._quiet_save is False
        jm.submit_file.assert_called_once_with(dlg.saved_path, name="benzene")
        dlg.accept.assert_called_once()

    def test_an_unchanged_saved_input_is_not_saved_again(self, generator, tmp_path):
        jm = MagicMock()
        jm.submit_file.return_value = True
        dlg = _dialog(generator, tmp_path, jm)
        dlg.save_file()
        dlg.submit_to_cluster()
        assert dlg.saves == 1
        jm.submit_file.assert_called_once()

    def test_an_edit_after_saving_saves_again(self, generator, tmp_path):
        jm = MagicMock()
        jm.submit_file.return_value = True
        dlg = _dialog(generator, tmp_path, jm)
        dlg.save_file()
        dlg.preview_text.text = "edited"
        dlg.submit_to_cluster()
        assert dlg.saves == 2
        assert Path(dlg.saved_path).read_text(encoding="utf-8") == "edited"

    def test_a_deleted_file_is_saved_again(self, generator, tmp_path):
        jm = MagicMock()
        jm.submit_file.return_value = True
        dlg = _dialog(generator, tmp_path, jm)
        dlg.save_file()
        os.remove(dlg.saved_path)
        dlg.submit_to_cluster()
        assert dlg.saves == 2 and os.path.isfile(dlg.saved_path)

    def test_cancelling_the_save_submits_nothing(self, generator, tmp_path):
        jm = MagicMock()
        dlg = _dialog(generator, tmp_path, jm, cancel=True)
        dlg.submit_to_cluster()
        jm.submit_file.assert_not_called()
        dlg.accept.assert_not_called()
        dlg_box.warning.assert_not_called()
        assert dlg._quiet_save is False

    def test_cancelling_a_resave_does_not_submit_the_old_file(self, generator, tmp_path):
        jm = MagicMock()
        dlg = _dialog(generator, tmp_path, jm)
        dlg.save_file()
        dlg.preview_text.text = "edited"
        dlg.save_file = lambda: None  # the user cancels the file dialog
        dlg.submit_to_cluster()
        jm.submit_file.assert_not_called()

    @pytest.mark.parametrize("outcome", [False, RuntimeError("boom")])
    def test_a_refusal_warns_and_keeps_the_dialog(self, generator, tmp_path, outcome):
        jm = MagicMock()
        if isinstance(outcome, Exception):
            jm.submit_file.side_effect = outcome
        else:
            jm.submit_file.return_value = outcome
        dlg = _dialog(generator, tmp_path, jm)
        dlg.submit_to_cluster()
        dlg_box.warning.assert_called_once()
        dlg.accept.assert_not_called()

    def test_job_manager_gone_since_the_dialog_opened(self, generator, tmp_path):
        dlg = _dialog(generator, tmp_path, None)
        dlg.submit_to_cluster()
        dlg_box.warning.assert_called_once()
        dlg.accept.assert_not_called()


class TestSaveFileRecordsTheFile:
    """The real save_file remembers what it wrote and stays quiet when asked."""

    def _save(self, generator, tmp_path, quiet):
        target = str(tmp_path / ("job" + generator.ext))
        box = MagicMock()
        dialog = MagicMock()
        dialog.getSaveFileName.return_value = (target, "")
        save = extract_function(
            generator.path,
            generator.cls,
            "save_file",
            {"os": os, "re": re, "QMessageBox": box, "QFileDialog": dialog},
        )
        text = {
            ".nw": "start old\ngeometry\nend\n",
            ".dat": "psi4.core.set_output_file('old.out', False)\n",
            ".py": "mf.chkfile = 'old.chk'\n",
        }.get(generator.ext, "PM7\n\n\n")
        dlg = SimpleNamespace(
            _generation_error=None,
            filename=None,
            preview_text=MagicMock(),
            saved_path=None,
            _saved_preview=None,
            _quiet_save=quiet,
        )
        dlg.preview_text.toPlainText.return_value = text
        save(dlg)
        return dlg, box, target

    def test_records_the_path_and_preview(self, generator, tmp_path):
        dlg, box, target = self._save(generator, tmp_path, quiet=False)
        assert dlg.saved_path == target and os.path.isfile(target)
        assert dlg._saved_preview == dlg.preview_text.toPlainText.return_value
        box.information.assert_called_once()

    def test_quiet_save_shows_no_box(self, generator, tmp_path):
        dlg, box, _ = self._save(generator, tmp_path, quiet=True)
        assert dlg.saved_path
        box.information.assert_not_called()


def test_the_button_is_wired_and_hidden_without_job_manager(generator):
    source = generator.path.read_text(encoding="utf-8")
    assert 'QPushButton("Submit to Cluster...")' in source
    assert "self.btn_submit.clicked.connect(self.submit_to_cluster)" in source
    assert re.search(
        r"self\.btn_submit\.setVisible\(find_job_manager\(self\.parent\(\)\) is not None\)",
        source,
    )
