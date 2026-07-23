"""
Headless GUI tests for the Compound Info Report plugin.

Covers: ReportDialog.

generate_initial_report() calls build_html() which calls re.sub() on a
MagicMock rdkit formula — that raises TypeError.  Patch it to a no-op so
setup_ui() runs (creating all widgets) without triggering rdkit calls.
"""

from __future__ import annotations

import contextlib
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest
from PyQt6.QtCore import QRectF
from PyQt6.QtWidgets import QDialog, QWidget

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

COMPOUND_PATH = PLUGINS_DIR / "Compound_Info_Report" / "compound_info_report.py"

with mock_chemistry_imports():
    _compound = load_plugin_for_gui(COMPOUND_PATH)


# The `test-gui` CI job installs only pytest + PyQt6 -- rdkit/PIL are NOT
# present there (unlike this local dev environment, where they are real
# installed packages). The real-chemistry test classes below are gated on
# HAS_REAL_CHEM so they run (and contribute real coverage) locally, but skip
# cleanly in CI instead of raising ModuleNotFoundError at import time.
try:
    import PIL  # noqa: F401
    import PIL.Image  # noqa: F401
    import rdkit  # noqa: F401
    import rdkit.Chem  # noqa: F401
    from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors  # noqa: F401

    HAS_REAL_CHEM = True
except ImportError:
    HAS_REAL_CHEM = False

requires_real_chem = pytest.mark.skipif(
    not HAS_REAL_CHEM, reason="requires real rdkit + PIL (not installed in CI test-gui job)"
)


@contextlib.contextmanager
def _mock_chemistry_keep_real_rdkit():
    """Like mock_chemistry_imports(), but rdkit resolves to the real package.

    Needed so build_html()/capture_scene_image()/fetch_pubchem_data() exercise
    real RDKit descriptor math instead of chasing MagicMock attribute chains.
    """
    real_mods = {
        k: v
        for k, v in sys.modules.items()
        if k == "rdkit"
        or k.startswith("rdkit.")
        or k == "PIL"
        or k.startswith("PIL.")
    }
    with mock_chemistry_imports():
        sys.modules.update(real_mods)
        yield


# Second module instance with real rdkit, used for tests that drive real
# chemistry through the plugin's bound methods (build_html, capture image,
# fetch_pubchem_data, generate_report, print/save).
if HAS_REAL_CHEM:
    with _mock_chemistry_keep_real_rdkit():
        _compound_real = load_plugin_for_gui(COMPOUND_PATH)
else:
    _compound_real = None


def _ethanol_mol():
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    return mol


# ===========================================================================
# ReportDialog  (Compound Info Report)
# ===========================================================================


class TestReportDialog:
    """ReportDialog with mol=None and generate_initial_report patched."""

    @pytest.fixture
    def dlg(self, qapp):
        with patch.object(
            _compound.ReportDialog, "generate_initial_report", lambda self: None
        ):
            d = _compound.ReportDialog(mol=None, parent=None)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Compound Info Report"

    def test_default_size_tall_enough(self, dlg):
        assert dlg.height() >= 400

    def test_pubchem_checkbox_initially_unchecked(self, dlg):
        assert not dlg.chk_pubchem.isChecked()

    def test_preview_is_text_browser(self, dlg):
        from PyQt6.QtWidgets import QTextBrowser
        assert isinstance(dlg.preview, QTextBrowser)

    def test_preview_has_external_links_enabled(self, dlg):
        assert dlg.preview.openExternalLinks()

    def test_print_button_exists(self, dlg):
        assert dlg.btn_print.text() == "Print..."

    def test_pdf_button_exists(self, dlg):
        assert dlg.btn_pdf.text() == "Save PDF..."

    def test_close_button_exists(self, dlg):
        assert dlg.btn_close.text() == "Close"

    def test_calculate_adducts_returns_list(self, dlg):
        adducts = dlg.calculate_adducts(180.0634)
        assert isinstance(adducts, list)
        assert len(adducts) > 0

    def test_adduct_mh_plus_correct(self, dlg):
        adducts = dlg.calculate_adducts(100.0)
        names = [a[0] for a in adducts]
        assert any("M+H" in n for n in names)

    def test_urllib_disabled_disables_checkbox(self, qapp, monkeypatch):
        monkeypatch.setattr(_compound, "URLLIB_AVAILABLE", False)
        with patch.object(
            _compound.ReportDialog, "generate_initial_report", lambda self: None
        ):
            d = _compound.ReportDialog(mol=None, parent=None)
        try:
            assert not d.chk_pubchem.isEnabled()
            assert not d.chk_pubchem.isChecked()
            assert "Disabled" in d.chk_pubchem.text()
        finally:
            d.destroy()


# ===========================================================================
# ReportDialog with real RDKit chemistry (module `_compound_real`)
# ===========================================================================


class _AtomItem:
    """Stand-in for the main app's QGraphicsItem-derived AtomItem."""

    def __init__(self, rect, visible=True):
        self._rect = rect
        self._visible = visible

    def isVisible(self):
        return self._visible

    def sceneBoundingRect(self):
        return self._rect


class _FakeScene:
    def __init__(self, items=None, bounding_rect=None, raise_on_items=False):
        self._items = items or []
        self._bounding_rect = bounding_rect or QRectF()
        self._raise_on_items = raise_on_items

    def items(self):
        if self._raise_on_items:
            raise RuntimeError("scene boom")
        return self._items

    def itemsBoundingRect(self):
        return self._bounding_rect

    def render(self, painter, target, source):
        pass


def _real_dlg(mol=None, parent=None):
    with patch.object(
        _compound_real.ReportDialog, "generate_initial_report", lambda self: None
    ):
        return _compound_real.ReportDialog(mol=mol, parent=parent)


@requires_real_chem
class TestCaptureSceneImage:
    def test_no_parent_scene_falls_back_to_rdkit(self, qapp):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol, parent=None)
        try:
            img_b64, w, h = d.capture_scene_image()
            assert img_b64
            assert (w, h) == (400, 300)
        finally:
            d.destroy()

    def test_empty_scene_falls_back_to_rdkit(self, qapp):
        mol = _ethanol_mol()
        parent = QWidget()
        parent.scene = _FakeScene(items=[])
        d = _real_dlg(mol=mol, parent=parent)
        try:
            img_b64, w, h = d.capture_scene_image()
            assert img_b64
            assert (w, h) == (400, 300)
        finally:
            d.destroy()
            parent.deleteLater()

    def test_scene_with_visible_atom_items_captured(self, qapp):
        mol = _ethanol_mol()
        parent = QWidget()
        rect = QRectF(0, 0, 50, 40)
        item = _AtomItem(rect, visible=True)
        # type(item).__name__ must equal "AtomItem" for the scene branch
        item.__class__.__name__ = "AtomItem"
        parent.scene = _FakeScene(items=[item])
        d = _real_dlg(mol=mol, parent=parent)
        try:
            img_b64, w, h = d.capture_scene_image()
            assert img_b64
            assert w > 0 and h > 0
        finally:
            d.destroy()
            parent.deleteLater()

    def test_invisible_items_ignored_uses_items_bounding_rect_fallback(self, qapp):
        mol = _ethanol_mol()
        parent = QWidget()
        item = _AtomItem(QRectF(0, 0, 50, 40), visible=False)
        item.__class__.__name__ = "AtomItem"
        # itemsBoundingRect also empty -> fallback to RDKit
        parent.scene = _FakeScene(items=[item], bounding_rect=QRectF())
        d = _real_dlg(mol=mol, parent=parent)
        try:
            img_b64, w, h = d.capture_scene_image()
            assert img_b64
            assert (w, h) == (400, 300)
        finally:
            d.destroy()
            parent.deleteLater()

    def test_scene_exception_falls_back_to_rdkit(self, qapp):
        mol = _ethanol_mol()
        parent = QWidget()
        parent.scene = _FakeScene(raise_on_items=True)
        d = _real_dlg(mol=mol, parent=parent)
        try:
            img_b64, w, h = d.capture_scene_image()
            assert img_b64
            assert (w, h) == (400, 300)
        finally:
            d.destroy()
            parent.deleteLater()


@requires_real_chem
class TestBuildHtml:
    def _empty_pubchem(self):
        return {
            "common_name": "",
            "density": None,
            "phys_desc": None,
            "cas_numbers": [],
        }

    def test_formula_subscripted_and_no_pubchem_section(self, qapp):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            html = d.build_html(self._empty_pubchem())
            assert "C<sub>2</sub>H<sub>6</sub>O" in html or "<sub>" in html
            assert "Online Data (PubChem)" not in html
            assert "Structure" not in html
        finally:
            d.destroy()

    def test_pubchem_data_renders_all_rows(self, qapp):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            pd = {
                "common_name": "Ethanol",
                "density": "0.789 g/cm3",
                "phys_desc": "Colorless liquid",
                "cas_numbers": ["64-17-5"],
            }
            html = d.build_html(pd)
            assert "Ethanol" in html
            assert "0.789 g/cm3" in html
            assert "Colorless liquid" in html
            assert "64-17-5" in html
            assert "Online Data (PubChem)" in html
        finally:
            d.destroy()

    def test_image_scaled_down_when_oversized(self, qapp):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            html = d.build_html(
                self._empty_pubchem(), img_b64="ZmFrZQ==", img_w=800, img_h=600
            )
            assert "Structure" in html
            # scale = min(400/800, 250/600, 0.5) = 0.41666... -> disp_w = 333
            assert 'width="333"' in html
            assert 'height="250"' in html
        finally:
            d.destroy()


@requires_real_chem
class TestFetchPubchemData:
    def test_unchecked_returns_default(self, qapp):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            data = d.fetch_pubchem_data()
            assert data == {
                "common_name": "",
                "density": None,
                "phys_desc": None,
                "cas_numbers": [],
            }
        finally:
            d.destroy()

    def test_checked_full_pipeline_with_progress_callback(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            d.chk_pubchem.setChecked(True)
            monkeypatch.setattr(
                _compound_real.PubChemFetcher,
                "get_synonyms",
                staticmethod(lambda key: ["Ethanol", "64-17-5"]),
            )
            monkeypatch.setattr(
                _compound_real.PubChemFetcher, "get_cid", staticmethod(lambda key: 702)
            )
            monkeypatch.setattr(
                _compound_real.PubChemFetcher,
                "fetch_experimental_properties",
                staticmethod(lambda cid: ("0.789 g/cm3", "Colorless liquid")),
            )
            calls = []
            data = d.fetch_pubchem_data(lambda msg, val: calls.append((msg, val)))
            assert data["common_name"] == "Ethanol"
            assert "64-17-5" in data["cas_numbers"]
            assert data["density"] == "0.789 g/cm3"
            assert data["phys_desc"] == "Colorless liquid"
            assert len(calls) == 3
        finally:
            d.destroy()

    def test_checked_exception_returns_default(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            d.chk_pubchem.setChecked(True)

            def _boom(*args, **kwargs):
                raise RuntimeError("inchikey boom")

            monkeypatch.setattr(_compound_real.Chem, "MolToInchiKey", _boom)
            data = d.fetch_pubchem_data()
            assert data["common_name"] == ""
            assert data["cas_numbers"] == []
        finally:
            d.destroy()

    def test_checked_no_cid_skips_experimental_properties(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            d.chk_pubchem.setChecked(True)
            monkeypatch.setattr(
                _compound_real.PubChemFetcher, "get_synonyms", staticmethod(lambda k: [])
            )
            monkeypatch.setattr(
                _compound_real.PubChemFetcher, "get_cid", staticmethod(lambda k: None)
            )
            called = []
            monkeypatch.setattr(
                _compound_real.PubChemFetcher,
                "fetch_experimental_properties",
                staticmethod(lambda cid: called.append(cid) or (None, None)),
            )
            data = d.fetch_pubchem_data()
            assert not called
            assert data["density"] is None
        finally:
            d.destroy()


@requires_real_chem
class TestGenerateReport:
    def test_unchecked_generates_preview_and_caches(self, qapp):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            d.generate_report()
            assert d.preview.toHtml() != ""
            assert hasattr(d, "last_pubchem_data")
            assert hasattr(d, "last_img")
        finally:
            d.destroy()

    def test_checked_runs_progress_dialog_pipeline(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            d.chk_pubchem.setChecked(True)
            monkeypatch.setattr(
                _compound_real.PubChemFetcher,
                "get_synonyms",
                staticmethod(lambda k: ["Ethanol"]),
            )
            monkeypatch.setattr(
                _compound_real.PubChemFetcher, "get_cid", staticmethod(lambda k: None)
            )
            d.generate_report()
            assert d.last_pubchem_data["common_name"] == "Ethanol"
        finally:
            d.destroy()


@requires_real_chem
class TestPrintReport:
    def test_accepted_calls_preview_print(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            monkeypatch.setattr(
                _compound_real.QPrintDialog,
                "exec",
                lambda self: QDialog.DialogCode.Accepted,
            )
            printed = []
            monkeypatch.setattr(d.preview, "print", lambda printer: printed.append(printer))
            d.print_report()
            assert len(printed) == 1
        finally:
            d.destroy()

    def test_rejected_does_not_print(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            monkeypatch.setattr(
                _compound_real.QPrintDialog,
                "exec",
                lambda self: QDialog.DialogCode.Rejected,
            )
            printed = []
            monkeypatch.setattr(d.preview, "print", lambda printer: printed.append(printer))
            d.print_report()
            assert printed == []
        finally:
            d.destroy()


@requires_real_chem
class TestSavePdf:
    def test_cancel_dialog_does_nothing(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            monkeypatch.setattr(
                _compound_real.QFileDialog, "getSaveFileName", lambda *a, **k: ("", "")
            )
            d.save_pdf()  # should return silently
        finally:
            d.destroy()

    def test_success_writes_pdf_using_cached_image(self, qapp, monkeypatch, tmp_path):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            d.generate_report()  # populate last_img / last_pubchem_data
            out_file = tmp_path / "report.pdf"
            monkeypatch.setattr(
                _compound_real.QFileDialog,
                "getSaveFileName",
                lambda *a, **k: (str(out_file), "PDF Files (*.pdf)"),
            )
            infos = []
            monkeypatch.setattr(
                _compound_real.QMessageBox,
                "information",
                staticmethod(lambda *a, **k: infos.append(a)),
            )
            d.save_pdf()
            assert out_file.exists()
            assert infos
        finally:
            d.destroy()

    def test_fetches_pubchem_data_when_not_cached(self, qapp, monkeypatch, tmp_path):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            d.last_pubchem_data = None
            out_file = tmp_path / "report2.pdf"
            monkeypatch.setattr(
                _compound_real.QFileDialog,
                "getSaveFileName",
                lambda *a, **k: (str(out_file), "PDF Files (*.pdf)"),
            )
            monkeypatch.setattr(
                _compound_real.QMessageBox, "information", staticmethod(lambda *a, **k: None)
            )
            d.save_pdf()
            assert d.last_pubchem_data is not None
            assert out_file.exists()
        finally:
            d.destroy()

    def test_progress_canceled_returns_before_building(self, qapp, monkeypatch, tmp_path):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            d.generate_report()
            out_file = tmp_path / "report3.pdf"
            monkeypatch.setattr(
                _compound_real.QFileDialog,
                "getSaveFileName",
                lambda *a, **k: (str(out_file), "PDF Files (*.pdf)"),
            )
            monkeypatch.setattr(
                _compound_real.QProgressDialog, "wasCanceled", lambda self: True
            )
            d.save_pdf()
            assert not out_file.exists()
        finally:
            d.destroy()

    def test_exception_shows_critical_message(self, qapp, monkeypatch, tmp_path):
        mol = _ethanol_mol()
        d = _real_dlg(mol=mol)
        try:
            d.generate_report()
            out_file = tmp_path / "report4.pdf"
            monkeypatch.setattr(
                _compound_real.QFileDialog,
                "getSaveFileName",
                lambda *a, **k: (str(out_file), "PDF Files (*.pdf)"),
            )

            def _boom(*a, **k):
                raise RuntimeError("build boom")

            monkeypatch.setattr(d, "build_html", _boom)
            errors = []
            monkeypatch.setattr(
                _compound_real.QMessageBox,
                "critical",
                staticmethod(lambda *a, **k: errors.append(a)),
            )
            d.save_pdf()
            assert errors
        finally:
            d.destroy()


# ===========================================================================
# show_report(context)  (module-level function)
# ===========================================================================


class _DummyDialog:
    last_instance = None

    def __init__(self, mol, parent=None):
        self.mol = mol
        self.parent_ = parent
        self.exec_called = False
        _DummyDialog.last_instance = self

    def exec(self):
        self.exec_called = True


@requires_real_chem
class TestShowReport:
    def test_no_rdkit_shows_warning(self, qapp, monkeypatch):
        monkeypatch.setattr(_compound_real, "Chem", None)
        mw = SimpleNamespace()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = None
        warned = []
        monkeypatch.setattr(
            _compound_real.QMessageBox,
            "warning",
            staticmethod(lambda *a, **k: warned.append(a)),
        )
        _compound_real.show_report(ctx)
        assert warned

    def test_no_molecule_no_fallback_shows_warning(self, qapp, monkeypatch):
        mw = SimpleNamespace()  # no state_manager attr
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = None
        warned = []
        monkeypatch.setattr(
            _compound_real.QMessageBox,
            "warning",
            staticmethod(lambda *a, **k: warned.append(a)),
        )
        monkeypatch.setattr(_compound_real, "ReportDialog", _DummyDialog)
        _compound_real.show_report(ctx)
        assert warned
        assert _DummyDialog.last_instance is None or not _DummyDialog.last_instance.exec_called

    def test_fallback_state_manager_used(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        mw = SimpleNamespace(
            state_manager=SimpleNamespace(data=SimpleNamespace(to_rdkit_mol=lambda: mol))
        )
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = None
        monkeypatch.setattr(_compound_real, "ReportDialog", _DummyDialog)
        _DummyDialog.last_instance = None
        _compound_real.show_report(ctx)
        assert _DummyDialog.last_instance is not None
        assert _DummyDialog.last_instance.mol is mol
        assert _DummyDialog.last_instance.exec_called

    def test_molecule_without_conformer_gets_2d_coords(self, qapp, monkeypatch):
        from rdkit import Chem

        mol = Chem.MolFromSmiles("CCO")  # no conformer yet
        mw = SimpleNamespace()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = mol
        monkeypatch.setattr(_compound_real, "ReportDialog", _DummyDialog)
        _DummyDialog.last_instance = None
        _compound_real.show_report(ctx)
        assert _DummyDialog.last_instance is not None
        assert _DummyDialog.last_instance.mol.GetNumConformers() >= 1

    def test_molecule_with_conformer_skips_2d_coords(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        n_conf_before = mol.GetNumConformers()
        assert n_conf_before >= 1
        mw = SimpleNamespace()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = mol
        monkeypatch.setattr(_compound_real, "ReportDialog", _DummyDialog)
        _DummyDialog.last_instance = None
        _compound_real.show_report(ctx)
        assert _DummyDialog.last_instance.mol.GetNumConformers() == n_conf_before

    def test_compute2dcoords_exception_is_silenced(self, qapp, monkeypatch):
        mol = _ethanol_mol()
        mw = SimpleNamespace()
        ctx = MagicMock()
        ctx.get_main_window.return_value = mw
        ctx.current_molecule = mol

        def _boom(m):
            raise RuntimeError("2d boom")

        # Force GetNumConformers() == 0 path but AllChem.Compute2DCoords raises.
        monkeypatch.setattr(mol, "GetNumConformers", lambda: 0)
        monkeypatch.setattr(_compound_real.AllChem, "Compute2DCoords", _boom)
        monkeypatch.setattr(_compound_real, "ReportDialog", _DummyDialog)
        _DummyDialog.last_instance = None
        _compound_real.show_report(ctx)  # must not raise
        assert _DummyDialog.last_instance is not None
