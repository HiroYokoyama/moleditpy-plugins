"""
Tests for the Molecule Comparator plugin.
"""

from __future__ import annotations

import math
import numpy as _ensure_real_numpy_imported  # noqa: F401  (see mocks_with_real_numpy note below)
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from conftest import (
    extract_function,
    load_plugin,
    make_context,
    mock_optional_imports,
    mocks_with_real_numpy,
    FakeConf,
)

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
MOLECULE_COMPARATOR_PATH = PLUGINS_DIR / "Molecule_Comparator" / "molecule_comparator.py"


class TestMoleculeComparator:
    def test_initialize_registers_document_reset_handler(self):
        """initialize() wires up on_reset via register_document_reset_handler."""
        with mock_optional_imports():
            mod = load_plugin(MOLECULE_COMPARATOR_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.register_document_reset_handler.assert_called_once()

    def test_get_text_color_for_background_dark(self):
        """Luminance formula: dark backgrounds map to white text."""
        # Mirror the algorithm from MoleculeComparator._get_text_color_for_background
        def _luminance_text(hex_color):
            try:
                h = hex_color.lstrip("#")
                r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
                return "white" if (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.5 else "black"
            except Exception:
                return "black"

        assert _luminance_text("#000000") == "white"
        assert _luminance_text("#0000FF") == "white"

    def test_get_text_color_for_background_light(self):
        """Luminance formula: light backgrounds map to black text."""
        def _luminance_text(hex_color):
            try:
                h = hex_color.lstrip("#")
                r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
                return "white" if (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.5 else "black"
            except Exception:
                return "black"

        assert _luminance_text("#FFFFFF") == "black"
        assert _luminance_text("#FFFF00") == "black"

    def test_get_text_color_for_background_invalid(self):
        """Luminance formula falls back to 'black' for non-hex input."""
        def _luminance_text(hex_color):
            try:
                h = hex_color.lstrip("#")
                r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
                return "white" if (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.5 else "black"
            except Exception:
                return "black"

        assert _luminance_text("not_a_color") == "black"




# ---------------------------------------------------------------------------
# document-reset handler behaviour
# ---------------------------------------------------------------------------


class TestComparatorResetHandler:
    def _handler(self, mw):
        with mock_optional_imports():
            mod = load_plugin(MOLECULE_COMPARATOR_PATH)
            ctx = make_context()
            ctx.get_main_window.return_value = mw
            mod.initialize(ctx)
        return ctx.register_document_reset_handler.call_args[0][0]

    def test_reset_closes_open_window(self):
        mw = SimpleNamespace(molecule_comparator_window=MagicMock())
        handler = self._handler(mw)
        handler()
        mw.molecule_comparator_window.close.assert_called_once()

    def test_reset_noop_without_window(self):
        mw = SimpleNamespace()  # no molecule_comparator_window attribute
        handler = self._handler(mw)
        handler()


# ---------------------------------------------------------------------------
# _get_text_color_for_background — real extracted method (pure string/math,
# no Qt/self attributes touched other than the color_hex argument)
# ---------------------------------------------------------------------------


def _text_color_fn():
    return extract_function(
        MOLECULE_COMPARATOR_PATH, "MoleculeComparator", "_get_text_color_for_background"
    )


class TestGetTextColorForBackgroundReal:
    def test_dark_background_gives_white_text(self):
        fn = _text_color_fn()
        assert fn(None, "#000000") == "white"
        assert fn(None, "#0000FF") == "white"

    def test_light_background_gives_black_text(self):
        fn = _text_color_fn()
        assert fn(None, "#FFFFFF") == "black"
        assert fn(None, "#FFFF00") == "black"

    def test_invalid_hex_falls_back_to_black(self):
        fn = _text_color_fn()
        assert fn(None, "not_a_color") == "black"

    def test_boundary_luminance_value(self):
        # Pure grey ~0.5 luminance is right at the threshold; exercise both sides.
        fn = _text_color_fn()
        # Luminance(#7F7F7F) < 0.5 -> white text
        assert fn(None, "#7F7F7F") == "white"
        # Luminance(#FFFFFF) = 1.0 -> black text (already covered above)


# ---------------------------------------------------------------------------
# AlignmentWorker.run — real Kabsch-alignment RMSD computed against
# hand-verifiable geometries, via a monkeypatched AllChem.AlignMol that does
# a genuine best-fit rigid alignment using real numpy (no rdkit involved).
# ---------------------------------------------------------------------------


class _RMSDAtom:
    def __init__(self, idx, atomic_num=6):
        self._idx = idx
        self._atomic_num = atomic_num

    def GetIdx(self):
        return self._idx

    def GetAtomicNum(self):
        return self._atomic_num


class _RMSDMol:
    """Minimal fake rdkit Mol exposing exactly what AlignmentWorker.run needs."""

    def __init__(self, coords, atomic_nums=None):
        self.conf = FakeConf(coords)
        self.atoms = [
            _RMSDAtom(i, (atomic_nums[i] if atomic_nums else 6))
            for i in range(len(coords))
        ]

    def GetAtoms(self):
        return list(self.atoms)

    def GetNumAtoms(self):
        return len(self.atoms)


def _kabsch_align_mol(probe, ref, atomMap, reflect=False):
    """Real best-fit (Kabsch) alignment RMSD, used as a stand-in for
    rdkit.Chem.AllChem.AlignMol under mocks_with_real_numpy()."""
    import numpy as np

    probe_pts = np.array(
        [
            (
                probe.conf.coords[p].x,
                probe.conf.coords[p].y,
                probe.conf.coords[p].z,
            )
            for p, _r in atomMap
        ]
    )
    ref_pts = np.array(
        [
            (ref.conf.coords[r].x, ref.conf.coords[r].y, ref.conf.coords[r].z)
            for _p, r in atomMap
        ]
    )

    probe_centroid = probe_pts.mean(axis=0)
    ref_centroid = ref_pts.mean(axis=0)
    pc = probe_pts - probe_centroid
    rc = ref_pts - ref_centroid

    H = pc.T @ rc
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T

    aligned = (R @ pc.T).T + ref_centroid
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - ref_pts) ** 2, axis=1))))
    return rmsd


def _run_fn():
    import copy as _copy
    import logging as _logging

    globs = {
        "Chem": SimpleNamespace(RemoveHs=lambda m: m),
        "AllChem": SimpleNamespace(AlignMol=_kabsch_align_mol),
        "copy": _copy,
        "logging": _logging,
    }
    return extract_function(MOLECULE_COMPARATOR_PATH, "AlignmentWorker", "run", globs)


class _WorkerStub:
    def __init__(self, ref_work, targets, method="Atom IDs", ignore_hs=False):
        self.ref_work = ref_work
        self.targets = targets
        self.method = method
        self.ignore_hs = ignore_hs
        self.MAX_COMBINATIONS = 100000
        self.progress = MagicMock()
        self.finished_signal = MagicMock()
        self.error_signal = MagicMock()
        self._interrupted = False

    def isInterruptionRequested(self):
        return self._interrupted


class TestAlignmentWorkerRun:
    def test_identical_structures_rmsd_zero(self):
        with mocks_with_real_numpy():
            fn = _run_fn()
            coords = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
            ref = _RMSDMol(coords)
            probe = _RMSDMol(coords)
            stub = _WorkerStub(ref, [(1, probe)])
            fn(stub)

        stub.finished_signal.emit.assert_called_once()
        results = stub.finished_signal.emit.call_args[0][0]
        assert len(results) == 1
        assert math.isclose(results[0]["rms"], 0.0, abs_tol=1e-6)

    def test_translated_copy_rmsd_zero_after_alignment(self):
        with mocks_with_real_numpy():
            fn = _run_fn()
            coords = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
            shift = (5.0, -3.0, 2.0)
            translated = [(x + shift[0], y + shift[1], z + shift[2]) for x, y, z in coords]
            ref = _RMSDMol(coords)
            probe = _RMSDMol(translated)
            stub = _WorkerStub(ref, [(0, probe)])
            fn(stub)

        results = stub.finished_signal.emit.call_args[0][0]
        assert math.isclose(results[0]["rms"], 0.0, abs_tol=1e-6)

    def test_mismatched_atom_count_skips_alignment(self):
        with mocks_with_real_numpy():
            fn = _run_fn()
            ref = _RMSDMol([(0, 0, 0), (1, 0, 0)])
            probe = _RMSDMol([(0, 0, 0), (1, 0, 0), (2, 0, 0)])
            stub = _WorkerStub(ref, [(0, probe)])
            fn(stub)

        results = stub.finished_signal.emit.call_args[0][0]
        # best_rms never set -> "rms" left at its initial -1.0 sentinel
        assert results[0]["rms"] == -1.0

    def test_ignore_hs_uses_only_heavy_atoms(self):
        with mocks_with_real_numpy():
            fn = _run_fn()
            # atom 1 is Hydrogen (atomic num 1) in both; only heavy atoms (0, 2)
            # need to line up for the RMSD to come out exactly zero.
            coords = [(0, 0, 0), (9, 9, 9), (1, 1, 1)]
            ref = _RMSDMol(coords, atomic_nums=[6, 1, 8])
            probe = _RMSDMol(coords, atomic_nums=[6, 1, 8])
            stub = _WorkerStub(ref, [(0, probe)], ignore_hs=True)
            fn(stub)

        results = stub.finished_signal.emit.call_args[0][0]
        assert math.isclose(results[0]["rms"], 0.0, abs_tol=1e-6)

    def test_interruption_stops_before_any_target_processed(self):
        with mocks_with_real_numpy():
            fn = _run_fn()
            coords = [(0, 0, 0), (1, 0, 0)]
            ref = _RMSDMol(coords)
            probe = _RMSDMol(coords)
            stub = _WorkerStub(ref, [(0, probe)])
            stub._interrupted = True
            fn(stub)

        results = stub.finished_signal.emit.call_args[0][0]
        assert results == []

    def test_progress_emitted_with_target_index(self):
        with mocks_with_real_numpy():
            fn = _run_fn()
            coords = [(0, 0, 0), (1, 0, 0)]
            ref = _RMSDMol(coords)
            probe = _RMSDMol(coords)
            stub = _WorkerStub(ref, [(2, probe)])
            fn(stub)

        stub.progress.emit.assert_called_once_with(0, 1, "Aligning Mol 3...")

    def test_exception_reports_via_error_signal(self):
        with mocks_with_real_numpy():
            fn = _run_fn()
            ref = _RMSDMol([(0, 0, 0)])
            # A None target: copy.deepcopy(None) -> None, then None.GetNumAtoms()
            # raises AttributeError inside the loop -> caught by run()'s outer
            # try/except and reported via error_signal instead of crashing.
            stub = _WorkerStub(ref, [(0, None)])
            fn(stub)

        stub.error_signal.emit.assert_called_once()
        assert "NoneType" in stub.error_signal.emit.call_args[0][0]
        stub.finished_signal.emit.assert_not_called()


# ---------------------------------------------------------------------------
# update_results_table / copy_results_to_clipboard / save_results_to_file —
# real extracted methods exercising the RMSD string-formatting logic
# ---------------------------------------------------------------------------


class _ResultItem:
    def __init__(self, text):
        self._text = text

    def text(self):
        return self._text


class _ResultsTableStub:
    def __init__(self):
        self.rows = 0
        self.items = {}

    def setRowCount(self, n):
        self.rows = n

    def setItem(self, row, col, item):
        self.items[(row, col)] = item


def _update_results_table_fn():
    globs = {"QTableWidgetItem": _ResultItem}
    return extract_function(
        MOLECULE_COMPARATOR_PATH, "MoleculeComparator", "update_results_table", globs
    )


class TestUpdateResultsTable:
    def test_formats_rms_values(self):
        fn = _update_results_table_fn()
        table = _ResultsTableStub()
        stub = SimpleNamespace(
            table_results=table,
            molecules=[
                {"name": "A", "rms": 0.12345},
                {"name": "B", "rms": -1.0},
                {"name": "C", "rms": None},
            ],
        )
        fn(stub)
        assert table.items[(0, 1)].text() == "0.1235"
        assert table.items[(1, 1)].text() == "N/A"
        assert table.items[(2, 1)].text() == "-"
        assert table.items[(0, 0)].text() == "A"


def _copy_results_fn():
    globs = {"QApplication": MagicMock()}
    return extract_function(
        MOLECULE_COMPARATOR_PATH, "MoleculeComparator", "copy_results_to_clipboard", globs
    ), globs["QApplication"]


class TestCopyResultsToClipboard:
    def test_empty_molecules_is_noop(self):
        fn, qapp = _copy_results_fn()
        stub = SimpleNamespace(molecules=[], context=MagicMock())
        fn(stub)
        qapp.clipboard.assert_not_called()

    def test_builds_tab_separated_text(self):
        # copy_results_to_clipboard does a *local* `from PyQt6.QtWidgets import
        # QApplication` — extra_globals can't intercept that, so run inside
        # mock_optional_imports() and configure the mocked QApplication that
        # import resolves to.
        with mock_optional_imports():
            import PyQt6.QtWidgets as _qtw

            fn, _unused = _copy_results_fn()
            qapp = _qtw.QApplication
            clipboard = MagicMock()
            qapp.clipboard.return_value = clipboard
            stub = SimpleNamespace(
                molecules=[
                    {"name": "A", "rms": 1.5},
                    {"name": "B", "rms": -1.0},
                ],
                context=MagicMock(),
            )
            fn(stub)

        text = clipboard.setText.call_args[0][0]
        assert "A\t1.5000" in text
        assert "B\tN/A" in text
        stub.context.show_status_message.assert_called_once()


class TestSaveResultsToFile:
    def test_no_molecules_warns_and_returns(self):
        fake_qmb = MagicMock()
        fn = extract_function(
            MOLECULE_COMPARATOR_PATH,
            "MoleculeComparator",
            "save_results_to_file",
            {"QMessageBox": fake_qmb, "QFileDialog": MagicMock()},
        )
        stub = SimpleNamespace(molecules=[], mw=MagicMock())
        fn(stub)
        fake_qmb.warning.assert_called_once()

    def test_writes_csv_with_escaped_commas(self, tmp_path):
        out_file = tmp_path / "results.csv"
        fake_qfd = SimpleNamespace(getSaveFileName=lambda *a, **k: (str(out_file), ""))
        fake_qmb = MagicMock()
        fn = extract_function(
            MOLECULE_COMPARATOR_PATH,
            "MoleculeComparator",
            "save_results_to_file",
            {"QFileDialog": fake_qfd, "QMessageBox": fake_qmb},
        )
        stub = SimpleNamespace(
            molecules=[
                {"name": "Mol, A", "rms": 2.0},
                {"name": "Mol B", "rms": None},
            ],
            mw=MagicMock(),
        )
        fn(stub)

        content = out_file.read_text(encoding="utf-8")
        assert '"Mol, A",2.0000' in content
        assert "Mol B,-" in content
        fake_qmb.information.assert_called_once()


# ---------------------------------------------------------------------------
# add_molecule_from_path — real rdkit parsing (.mol/.pdb/.xyz + error paths)
#
# Regression (2026.07.10): the file-open dialog advertises "*.xyz" in its
# filter, but add_molecule_from_path had no ".xyz" branch at all, so every
# XYZ file failed to load with "Failed to load molecule from ...". Fixed by
# adding an .xyz branch that parses via Chem.MolFromXYZBlock + rdDetermineBonds.
# ---------------------------------------------------------------------------

from rdkit import Chem as _real_Chem
from rdkit.Chem import rdDetermineBonds as _real_rdDetermineBonds

DEFAULT_COLORS = None
with mock_optional_imports():
    _mc_mod_for_colors = load_plugin(MOLECULE_COMPARATOR_PATH)
DEFAULT_COLORS = _mc_mod_for_colors.DEFAULT_COLORS

_WATER_XYZ = """3
water
O 0.000000 0.000000 0.000000
H 0.000000 0.000000 0.960000
H 0.930000 0.000000 -0.240000
"""


def _add_molecule_from_path_fn():
    globs = {
        "os": __import__("os"),
        "Chem": _real_Chem,
        "rdDetermineBonds": _real_rdDetermineBonds,
        "QMessageBox": MagicMock(),
        "logging": __import__("logging"),
        "DEFAULT_COLORS": DEFAULT_COLORS,
    }
    fn = extract_function(
        MOLECULE_COMPARATOR_PATH, "MoleculeComparator", "add_molecule_from_path", globs
    )
    return fn, globs


def _comparator_stub():
    return SimpleNamespace(
        molecules=[],
        mw=MagicMock(),
        update_list=MagicMock(),
        update_visualization=MagicMock(),
        update_wireframe_lighting=MagicMock(),
        reset_view=MagicMock(),
    )


class TestAddMoleculeFromPath:
    def test_loads_xyz_file_and_perceives_bonds(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text(_WATER_XYZ)
        stub = _comparator_stub()

        fn(stub, str(xyz_file))

        globs["QMessageBox"].warning.assert_not_called()
        globs["QMessageBox"].critical.assert_not_called()
        assert len(stub.molecules) == 1
        entry = stub.molecules[0]
        assert entry["name"] == "water.xyz"
        assert entry["mol"].GetNumAtoms() == 3
        assert entry["mol"].GetNumBonds() == 2  # O-H, O-H perceived
        stub.update_list.assert_called_once()
        stub.update_visualization.assert_called_once()

    def test_missing_file_shows_critical_not_warning(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        stub = _comparator_stub()
        fn(stub, str(tmp_path / "does_not_exist.xyz"))
        globs["QMessageBox"].critical.assert_called_once()
        globs["QMessageBox"].warning.assert_not_called()
        assert stub.molecules == []

    def test_unparseable_mol_file_shows_warning(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        bad_mol = tmp_path / "bad.mol"
        bad_mol.write_text("not a real mol file\ngarbage\n")
        stub = _comparator_stub()
        fn(stub, str(bad_mol))
        globs["QMessageBox"].warning.assert_called_once()
        assert stub.molecules == []

    def test_unsupported_extension_shows_warning(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        odd_file = tmp_path / "data.foobar"
        odd_file.write_text("irrelevant")
        stub = _comparator_stub()
        fn(stub, str(odd_file))
        globs["QMessageBox"].warning.assert_called_once()
        assert stub.molecules == []

    def test_color_cycles_with_existing_molecules(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text(_WATER_XYZ)
        stub = _comparator_stub()
        stub.molecules = [{"name": "existing", "mol": None, "color": "x", "rms": None}]
        fn(stub, str(xyz_file))
        assert stub.molecules[1]["color"] == DEFAULT_COLORS[1 % len(DEFAULT_COLORS)]
