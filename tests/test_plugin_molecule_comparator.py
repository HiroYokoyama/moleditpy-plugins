"""
Tests for the Molecule Comparator plugin.
"""

from __future__ import annotations

import math
import types
import numpy as _ensure_real_numpy_imported  # noqa: F401  (see mocks_with_real_numpy note below)
import pytest
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
        "Chem": SimpleNamespace(RemoveAllHs=lambda m: m),
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


def _mcs_text_fn():
    return extract_function(
        MOLECULE_COMPARATOR_PATH,
        None,
        "mcs_coverage_text",
        {"MIN_MCS_COVERAGE": 0.5},
    )


def _update_results_table_fn():
    globs = {"QTableWidgetItem": _ResultItem, "mcs_coverage_text": _mcs_text_fn()}
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
        assert table.items[(0, 2)].text() == "-"

    def test_matched_atoms_column(self):
        fn = _update_results_table_fn()
        table = _ResultsTableStub()
        stub = SimpleNamespace(
            table_results=table,
            molecules=[
                {"name": "A", "rms": 0.1, "matched": 72, "of": 74, "loose": True},
            ],
        )
        fn(stub)
        assert table.items[(0, 2)].text() == "72/74 (loose)"


def _copy_results_fn():
    globs = {"QApplication": MagicMock(), "mcs_coverage_text": _mcs_text_fn()}
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
        # copy_results_to_clipboard does a *local* `from PySide6.QtWidgets import
        # QApplication` — extra_globals can't intercept that, so run inside
        # mock_optional_imports() and configure the mocked QApplication that
        # import resolves to.
        with mock_optional_imports():
            import PySide6.QtWidgets as _qtw

            fn, qapp = _copy_results_fn()
            clipboard = MagicMock()
            qapp.clipboard.return_value = clipboard
            stub = SimpleNamespace(
                molecules=[
                    {"name": "A", "rms": 1.5, "matched": 4, "of": 74},
                    {"name": "B", "rms": -1.0},
                ],
                context=MagicMock(),
            )
            fn(stub)

        text = clipboard.setText.call_args[0][0]
        assert text.splitlines()[0] == "Molecule\tRMSD (\u00c5)\tMatched atoms"
        assert "A\t1.5000\t4/74 \u26a0" in text
        assert "B\tN/A\t-" in text
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
            {
                "QFileDialog": fake_qfd,
                "QMessageBox": fake_qmb,
                "csv": __import__("csv"),
                "mcs_coverage_text": _mcs_text_fn(),
            },
        )
        stub = SimpleNamespace(
            molecules=[
                {"name": "Mol, A", "rms": 2.0, "matched": 72, "of": 74},
                {"name": "Mol B", "rms": None},
            ],
            mw=MagicMock(),
        )
        fn(stub)

        content = out_file.read_text(encoding="utf-8")
        assert content.splitlines()[0] == "Molecule,RMSD (\u00c5),Matched atoms"
        assert '"Mol, A",2.0000,72/74' in content
        assert "Mol B,-,-" in content
        fake_qmb.information.assert_called_once()


# ---------------------------------------------------------------------------
# add_molecule_from_path — real rdkit parsing (.mol/.pdb/.xyz + error paths)
#
# Regression (2026.07.10): the file-open dialog advertises "*.xyz" in its
# filter, but add_molecule_from_path had no ".xyz" branch at all, so every
# XYZ file failed to load with "Failed to load molecule from ...". Fixed by
# adding an .xyz branch that parses via Chem.MolFromXYZBlock + rdDetermineBonds.
# ---------------------------------------------------------------------------

try:
    from rdkit import Chem as _real_Chem
    from rdkit.Chem import rdDetermineBonds as _real_rdDetermineBonds
except ImportError:  # CI runs without rdkit installed
    _real_Chem = None
    _real_rdDetermineBonds = None

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
    """Rebuild the real loading path (helpers included) against real rdkit."""
    try:
        from rdkit.Geometry import Point3D as _real_Point3D
    except ImportError:
        _real_Point3D = None

    globs = {
        "os": __import__("os"),
        "math": math,
        "Chem": _real_Chem,
        "rdDetermineBonds": _real_rdDetermineBonds,
        "Point3D": _real_Point3D,
        "QMessageBox": MagicMock(),
        "logging": __import__("logging"),
        "DEFAULT_COLORS": DEFAULT_COLORS,
        "_TEXT_ENCODINGS": _mc_mod_for_colors._TEXT_ENCODINGS,
    }
    # Extracted in dependency order: extract_function snapshots the globals it
    # is given, so each helper has to be in place before the next one is built.
    for helper in (
        "read_text_flexible",
        "_xyz_atom",
        "parse_xyz_rows",
        "mol_from_xyz_lines",
        "perceive_bonds",
        "mol_from_xyz_text",
        "has_3d_coords",
    ):
        globs[helper] = extract_function(
            MOLECULE_COMPARATOR_PATH, None, helper, globs
        )

    globs["_read_molecules"] = extract_function(
        MOLECULE_COMPARATOR_PATH, "MoleculeComparator", "_read_molecules", globs
    )
    fn = extract_function(
        MOLECULE_COMPARATOR_PATH,
        "MoleculeComparator",
        "add_molecules_from_paths",
        globs,
    )

    def add_one(stub, file_path):
        """Call the real worker the way add_molecule_from_path does."""
        return fn(stub, [file_path])

    return add_one, globs


def _comparator_stub(globs=None):
    stub = SimpleNamespace(
        molecules=[],
        # A plain namespace, not a MagicMock: the loader asks the host for
        # io_manager.load_xyz_block, and a MagicMock would answer with a mock
        # molecule instead of letting the rdkit path run.
        mw=SimpleNamespace(io_manager=None),
        update_list=MagicMock(),
        update_visualization=MagicMock(),
        update_wireframe_lighting=MagicMock(),
        reset_view=MagicMock(),
    )
    if globs is not None:
        stub._read_molecules = types.MethodType(globs["_read_molecules"], stub)
    return stub


@pytest.mark.skipif(_real_Chem is None, reason="requires real rdkit")
class TestAddMoleculeFromPath:
    def test_loads_xyz_file_and_perceives_bonds(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text(_WATER_XYZ)
        stub = _comparator_stub(globs)

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

    def test_missing_file_is_reported_with_its_name(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        stub = _comparator_stub(globs)
        fn(stub, str(tmp_path / "does_not_exist.xyz"))
        warned = globs["QMessageBox"].warning
        warned.assert_called_once()
        assert "does_not_exist.xyz" in warned.call_args[0][2]
        assert stub.molecules == []

    def test_unknown_extension_holding_xyz_still_loads(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        odd = tmp_path / "water.txt"
        odd.write_text(_WATER_XYZ)
        stub = _comparator_stub(globs)
        fn(stub, str(odd))
        assert len(stub.molecules) == 1
        assert stub.molecules[0]["mol"].GetNumAtoms() == 3

    def test_headerless_xyz_loads(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        f = tmp_path / "bare.xyz"
        f.write_text(
            "O 0.0 0.0 0.0\nH 0.0 0.0 0.96\nH 0.93 0.0 -0.24\n"
        )
        stub = _comparator_stub(globs)
        fn(stub, str(f))
        assert len(stub.molecules) == 1
        assert stub.molecules[0]["mol"].GetNumAtoms() == 3

    def test_trajectory_xyz_loads_only_its_first_frame(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        f = tmp_path / "traj.xyz"
        f.write_text(_WATER_XYZ + _WATER_XYZ)
        stub = _comparator_stub(globs)
        fn(stub, str(f))
        assert len(stub.molecules) == 1
        assert stub.molecules[0]["mol"].GetNumAtoms() == 3

    def test_two_files_land_as_two_entries(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        first = tmp_path / "a.xyz"
        second = tmp_path / "b.xyz"
        first.write_text(_WATER_XYZ)
        second.write_text(_WATER_XYZ)
        stub = _comparator_stub(globs)
        add_all = extract_function(
            MOLECULE_COMPARATOR_PATH,
            "MoleculeComparator",
            "add_molecules_from_paths",
            globs,
        )
        add_all(stub, [str(first), str(second)])
        assert [entry["name"] for entry in stub.molecules] == ["a.xyz", "b.xyz"]
        stub.update_list.assert_called_once()

    def test_unparseable_mol_file_shows_warning(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        bad_mol = tmp_path / "bad.mol"
        bad_mol.write_text("not a real mol file\ngarbage\n")
        stub = _comparator_stub(globs)
        fn(stub, str(bad_mol))
        globs["QMessageBox"].warning.assert_called_once()
        assert stub.molecules == []

    def test_file_that_holds_no_structure_shows_warning(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        odd_file = tmp_path / "data.foobar"
        odd_file.write_text("irrelevant")
        stub = _comparator_stub(globs)
        fn(stub, str(odd_file))
        globs["QMessageBox"].warning.assert_called_once()
        assert stub.molecules == []

    def test_color_cycles_with_existing_molecules(self, tmp_path):
        fn, globs = _add_molecule_from_path_fn()
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text(_WATER_XYZ)
        stub = _comparator_stub(globs)
        stub.molecules = [{"name": "existing", "mol": None, "color": "x", "rms": None}]
        fn(stub, str(xyz_file))
        assert stub.molecules[1]["color"] == DEFAULT_COLORS[1 % len(DEFAULT_COLORS)]


# ---------------------------------------------------------------------------
# MCS selection and coverage text (2026.11.13)
#
# Regression: structures without hydrogens (X-ray models before H placement,
# XYZ imports) come out of bond perception with radicals and different bond
# orders, so the strict MCS collapsed to ~4 atoms and the plugin reported a
# ~0 A RMSD for a badly displaced overlay.
# ---------------------------------------------------------------------------


class _FakeMol:
    def __init__(self, n):
        self._n = n

    def GetNumAtoms(self):
        return self._n


def _fake_fmcs(strict_atoms, loose_atoms):
    """rdFMCS stand-in answering the strict search first, then the loose one."""
    calls = []

    def find_mcs(mols, **kwargs):
        calls.append(kwargs)
        loose = "atomCompare" in kwargs
        return SimpleNamespace(numAtoms=loose_atoms if loose else strict_atoms)

    fake = SimpleNamespace(
        FindMCS=find_mcs,
        AtomCompare=SimpleNamespace(CompareElements="elements"),
        BondCompare=SimpleNamespace(CompareAny="any"),
    )
    return fake, calls


def _find_mcs_fn(fmcs):
    return extract_function(
        MOLECULE_COMPARATOR_PATH,
        None,
        "find_mcs",
        {"rdFMCS": fmcs, "MIN_MCS_COVERAGE": 0.5},
    )


class TestFindMcs:
    def test_adequate_strict_match_is_kept_without_a_second_search(self):
        fmcs, calls = _fake_fmcs(strict_atoms=72, loose_atoms=80)
        res, loose = _find_mcs_fn(fmcs)(_FakeMol(82), _FakeMol(74))
        assert (res.numAtoms, loose) == (72, False)
        assert len(calls) == 1
        assert calls[0]["matchValences"] is True

    def test_exactly_half_coverage_is_adequate(self):
        fmcs, calls = _fake_fmcs(strict_atoms=3, loose_atoms=6)
        res, loose = _find_mcs_fn(fmcs)(_FakeMol(6), _FakeMol(6))
        assert (res.numAtoms, loose) == (3, False)
        assert len(calls) == 1

    def test_small_strict_match_falls_back_to_loose(self):
        fmcs, calls = _fake_fmcs(strict_atoms=4, loose_atoms=75)
        res, loose = _find_mcs_fn(fmcs)(_FakeMol(82), _FakeMol(75))
        assert (res.numAtoms, loose) == (75, True)
        assert len(calls) == 2
        second = calls[1]
        assert second["matchValences"] is False
        assert second["atomCompare"] == "elements"
        assert second["bondCompare"] == "any"

    def test_loose_search_that_is_no_larger_keeps_strict(self):
        fmcs, _calls = _fake_fmcs(strict_atoms=4, loose_atoms=4)
        res, loose = _find_mcs_fn(fmcs)(_FakeMol(82), _FakeMol(75))
        assert (res.numAtoms, loose) == (4, False)

    def test_empty_molecule_still_returns_a_result(self):
        fmcs, _calls = _fake_fmcs(strict_atoms=0, loose_atoms=0)
        res, loose = _find_mcs_fn(fmcs)(_FakeMol(0), _FakeMol(5))
        assert (res.numAtoms, loose) == (0, False)


class TestMcsCoverageText:
    def test_unknown_when_nothing_matched(self):
        fn = _mcs_text_fn()
        assert fn({}) == "-"
        assert fn({"matched": 0, "of": 74}) == "-"
        assert fn({"matched": 4, "of": 0}) == "-"

    def test_plain_ratio(self):
        assert _mcs_text_fn()({"matched": 72, "of": 74}) == "72/74"

    def test_loose_is_labelled(self):
        text = _mcs_text_fn()({"matched": 72, "of": 74, "loose": True})
        assert text == "72/74 (loose)"

    def test_low_coverage_is_flagged(self):
        assert _mcs_text_fn()({"matched": 4, "of": 74}).endswith("\u26a0")

    def test_exactly_half_is_not_flagged(self):
        assert _mcs_text_fn()({"matched": 3, "of": 6}) == "3/6"


# Real RDKit: the case from the bug report. A probe whose bonds are all single
# and whose atoms carry no hydrogens (what bond perception makes of a
# hydrogen-free model) shares no strict MCS with the proper structure, yet is
# the same graph.
try:
    from rdkit import Chem as _rk_Chem
    from rdkit.Chem import rdFMCS as _rk_rdFMCS
except ImportError:  # CI runs without rdkit installed
    _rk_Chem = _rk_rdFMCS = None


@pytest.mark.skipif(_rk_rdFMCS is None, reason="rdkit not installed")
class TestFindMcsRealRdkit:
    SMILES = "c1ccc(cc1)C(=O)Nc1ccccc1C(=O)O"

    def _pair(self):
        ref = _rk_Chem.RemoveHs(_rk_Chem.AddHs(_rk_Chem.MolFromSmiles(self.SMILES)))
        rw = _rk_Chem.RWMol(ref)
        for bond in rw.GetBonds():
            bond.SetBondType(_rk_Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
        for atom in rw.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
        probe = rw.GetMol()
        _rk_Chem.SanitizeMol(
            probe, _rk_Chem.SANITIZE_ALL ^ _rk_Chem.SANITIZE_PROPERTIES
        )
        return ref, probe

    def test_hydrogen_free_probe_matches_in_full_only_when_loose(self):
        ref, probe = self._pair()
        strict = _rk_rdFMCS.FindMCS(
            [ref, probe],
            matchValences=True,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            timeout=5,
        )
        assert strict.numAtoms < 0.5 * probe.GetNumAtoms()
        res, loose = _find_mcs_fn(_rk_rdFMCS)(ref, probe)
        assert loose is True
        assert res.numAtoms == probe.GetNumAtoms()

    def test_identical_structures_stay_strict(self):
        ref, _probe = self._pair()
        res, loose = _find_mcs_fn(_rk_rdFMCS)(ref, _rk_Chem.Mol(ref))
        assert loose is False
        assert res.numAtoms == ref.GetNumAtoms()
