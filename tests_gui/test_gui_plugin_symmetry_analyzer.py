"""
Headless GUI tests for the Symmetry Analyzer plugin.

Covers: SymmetryAnalysisPlugin.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_symmetry_analyzer.py
"""

from __future__ import annotations

import contextlib
import logging
import math
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtTest import QTest

# Guarded so CI's bare test-gui job (only pytest+PyQt6 installed) skips this
# real-numpy module instead of erroring at collection.
np = pytest.importorskip("numpy")

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

SYMMETRY_PATH = PLUGINS_DIR / "Symmetry_Analyzer" / "symmetry_analyzer.py"

with mock_chemistry_imports():
    _symmetry = load_plugin_for_gui(SYMMETRY_PATH)


@contextlib.contextmanager
def _mock_chemistry_keep_numpy():
    """Like mock_chemistry_imports(), but numpy resolves to the real package.

    Needed so the plugin's real math (np.linalg.eig, np.arange, ...) actually
    runs instead of chasing MagicMock attribute chains.
    """
    real_np_mods = {
        k: v for k, v in sys.modules.items() if k == "numpy" or k.startswith("numpy.")
    }
    with mock_chemistry_imports():
        sys.modules.update(real_np_mods)
        yield


# Second module instance with real numpy, used for tests that drive real
# vector math through the plugin's bound methods (separate from `_symmetry`
# above, which keeps numpy mocked for the plain widget-construction tests).
with _mock_chemistry_keep_numpy():
    _symnp = load_plugin_for_gui(SYMMETRY_PATH)


def _ctx_no_mol() -> MagicMock:
    """Context with no main window and no active molecule."""
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_mol = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# SymmetryAnalysisPlugin  (visible plugin: "Symmetry Analyzer")
# ===========================================================================


class TestSymmetryAnalysisPlugin:
    """SymmetryAnalysisPlugin — pymatgen is mocked so HAS_PYMATGEN=True."""

    @pytest.fixture
    def dlg(self, qapp):
        ctx = _ctx_no_mol()
        d = _symmetry.SymmetryAnalysisPlugin(context=ctx)
        yield d
        d.destroy()

    def test_creates_without_error(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert "Symmetry" in dlg.windowTitle() or dlg.windowTitle() == ""

    def test_groups_list_is_empty_initially(self, dlg):
        assert dlg.groups_list.count() == 0

    def test_ops_list_is_empty_initially(self, dlg):
        assert dlg.ops_list.count() == 0

    def test_point_group_label_default(self, dlg):
        assert dlg.selected_group_label.text() == "Point Group: -"

    def test_op_details_is_readonly(self, dlg):
        assert dlg.op_details.isReadOnly()

    def test_max_tol_spin_default(self, dlg):
        # Beyond ~0.7 A distinct atoms merge and every reported group is an
        # artefact, so the default scan stops well short of that.
        assert dlg.max_tol_spin.value() == pytest.approx(0.5)

    def test_min_tol_spin_default(self, dlg):
        """Issue #9: the floor used to be hard-coded at 0.05 A, which reported
        D2 for a C2 molecule whose true C2 window was 0.02-0.04 A."""
        assert dlg.min_tol_spin.value() == pytest.approx(0.01)
        assert dlg.min_tol_spin.minimum() > 0.0

    def test_symmetrize_button_initially_disabled(self, dlg):
        assert not dlg.sym_btn.isEnabled()

    def test_analyze_button_exists(self, dlg):
        assert dlg.calc_btn is not None


# ===========================================================================
# Shared fakes for the real-numpy tests below (module `_symnp`)
# ===========================================================================


class _P3:
    """Minimal RDKit-Point3D-like: .x/.y/.z + sequence unpacking."""

    def __init__(self, x, y, z):
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __iter__(self):
        return iter((self.x, self.y, self.z))


class _NPAtom:
    #: Average atomic weights by Z, for the mass-weighted centre of mass.
    _MASSES = {1: 1.008, 6: 12.011, 7: 14.007, 8: 15.999, 9: 18.998}

    def __init__(self, z):
        self._z = z

    def GetAtomicNum(self):
        return self._z

    def GetNumRadicalElectrons(self):
        # Real RDKit atoms always expose this; without it the
        # derivation raised into the dialog's silent except.
        return getattr(self, '_radicals', 0)

    def GetMass(self):
        return self._MASSES.get(self._z, 12.011)


class _NPConf:
    def __init__(self, coords):
        self._pts = [_P3(*c) for c in coords]

    def GetAtomPosition(self, i):
        return self._pts[i]

    def SetAtomPosition(self, i, p):
        self._pts[i] = p


class _NPMol:
    """rdkit-Mol-like stand-in built from atomic numbers + 3D coords."""

    def __init__(self, atomic_nums, coords):
        self._atoms = [_NPAtom(z) for z in atomic_nums]
        self._conf = _NPConf(coords)

    def GetAtoms(self):
        return list(self._atoms)

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetConformer(self):
        return self._conf


class _NPMolNoConformer(_NPMol):
    def GetConformer(self):
        raise ValueError("no conformer")


class _Species:
    _MASSES = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998}

    def __init__(self, symbol):
        self.symbol = symbol
        self.mass = self._MASSES.get(symbol, 12.011)


class _FakePlotter:
    def __init__(self):
        self.added = []
        self.removed = []

    def add_mesh(self, mesh, **kwargs):
        actor = SimpleNamespace(mesh=mesh, kwargs=kwargs)
        self.added.append(actor)
        return actor

    def remove_actor(self, actor):
        self.removed.append(actor)

    def render(self):
        pass


def _op(matrix, translation=None):
    tv = translation if translation is not None else np.zeros(3)
    return SimpleNamespace(rotation_matrix=matrix, translation_vector=tv, as_xyz_str=lambda: "x, y, z")


_S3 = math.sqrt(3.0) / 2.0
_MAT_E = np.eye(3)
_MAT_C2 = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
_MAT_C3 = np.array([[-0.5, -_S3, 0.0], [_S3, -0.5, 0.0], [0.0, 0.0, 1.0]])
_MAT_SIGMA = np.diag([1.0, 1.0, -1.0])
_MAT_INV = -np.eye(3)
_MAT_S4 = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])


class _SymOp:
    """pymatgen SymmOp stand-in: real rotation about the origin."""

    def __init__(self, matrix):
        self.matrix = np.array(matrix, dtype=float)

    def operate(self, p):
        return self.matrix @ np.asarray(p, dtype=float)


def _rot_z(deg):
    th = np.radians(deg)
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _refl(axis):
    d = [1.0, 1.0, 1.0]
    d[axis] = -1.0
    return np.diag(d)


_WATER_COORDS = [(0.0, 0.0, 0.3), (0.8, 0.0, -0.3), (-0.8, 0.0, -0.3)]
_WATER_SPECIES = ["O", "H", "H"]
_WATER_OPS = [
    _SymOp(np.eye(3)),
    _SymOp(_rot_z(180)),
    _SymOp(_refl(1)),
    _SymOp(_refl(0)),
]


@pytest.fixture
def dlgnp(qapp):
    """SymmetryAnalysisPlugin built from the real-numpy module `_symnp`."""
    ctx = _ctx_no_mol()
    d = _symnp.SymmetryAnalysisPlugin(context=ctx)
    yield d
    d.destroy()


def _no_block_msgbox(monkeypatch):
    """QMessageBox.* pop up modal dialogs that block exec() under offscreen
    Qt with no user interaction; replace with recording stubs."""
    calls = {"warning": [], "information": [], "critical": []}
    for kind in calls:
        monkeypatch.setattr(
            _symnp.QMessageBox, kind,
            staticmethod(lambda *a, _k=kind, **kw: calls[_k].append(a)),
        )

    # question() blocks the same way, and unlike the others its RETURN VALUE
    # drives control flow, so it must answer rather than just record. Default
    # to Yes so the code under test proceeds and stays observable.
    calls["question"] = []

    def _question(*a, **kw):
        calls["question"].append(a)
        return _symnp.QMessageBox.StandardButton.Yes

    monkeypatch.setattr(_symnp.QMessageBox, "question", staticmethod(_question))
    return calls


# ===========================================================================
# get_pymatgen_molecule
# ===========================================================================


class TestGetPymatgenMoleculeReal:
    def test_no_current_molecule_returns_none(self, dlgnp):
        dlgnp.context.current_molecule = None
        assert dlgnp.get_pymatgen_molecule() is None

    def test_no_conformer_returns_none(self, dlgnp):
        dlgnp.context.current_molecule = _NPMolNoConformer([8, 1, 1], _WATER_COORDS)
        assert dlgnp.get_pymatgen_molecule() is None

    def test_species_and_coords_passed_to_molecule(self, dlgnp, monkeypatch):
        captured = {}

        def fake_molecule(species, coords):
            captured["species"] = species
            captured["coords"] = coords
            return "MOL"

        monkeypatch.setattr(_symnp, "Molecule", fake_molecule)
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        result = dlgnp.get_pymatgen_molecule()
        assert result == "MOL"
        assert captured["species"] == [8, 1, 1]
        assert [list(c) for c in captured["coords"]] == [list(c) for c in _WATER_COORDS]


# ===========================================================================
# SymmetryAnalysisWorker.run() — real numpy tolerance scan
# ===========================================================================


class _FakeAnalyzer:
    def __init__(self, tol):
        self.tol = tol
        self.sch_symbol = "C2v" if tol < 0.15 else "Cs"

    def get_symmetry_operations(self):
        return []


class TestWorkerRunReal:
    def test_run_emits_real_group_data(self, monkeypatch):
        monkeypatch.setattr(_symnp, "PointGroupAnalyzer", lambda mol, tolerance: _FakeAnalyzer(tolerance))
        worker = _symnp.SymmetryAnalysisWorker(mol_pmg=MagicMock(), min_tol=0.0, max_tol=0.2)
        received = {}
        worker.analysis_finished.connect(lambda gd, fa: received.update(group_data=gd, found_any=fa))
        worker.run()
        assert received["found_any"] is True
        assert set(received["group_data"].keys()) == {"C2v", "Cs"}
        # The refinement pass bisects the C2v/Cs transition, so the boundary is
        # located to ~1e-4 A rather than to the nearest grid step.
        highest_c2v = max(received["group_data"]["C2v"]["tols"])
        assert 0.1499 <= highest_c2v < 0.15
        assert min(received["group_data"]["Cs"]["tols"]) == pytest.approx(0.15)
        assert received["group_data"]["C2v"]["bands"] == [(0.0, highest_c2v)]
        assert len(received["group_data"]["Cs"]["bands"]) == 1

    def test_run_tolerates_analyzer_exceptions(self, monkeypatch):
        def _boom(mol, tolerance):
            raise ValueError("bad tolerance")

        monkeypatch.setattr(_symnp, "PointGroupAnalyzer", _boom)
        worker = _symnp.SymmetryAnalysisWorker(mol_pmg=MagicMock(), min_tol=0.0, max_tol=0.1)
        received = {}
        worker.analysis_finished.connect(lambda gd, fa: received.update(group_data=gd, found_any=fa))
        worker.run()
        assert received["found_any"] is False
        assert received["group_data"] == {}


# ===========================================================================
# analyze_symmetry / on_analysis_finished — real UI wiring, synchronous worker
# ===========================================================================


class TestAnalyzeSymmetryReal:
    def test_no_molecule_warns(self, dlgnp, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        dlgnp.context.current_molecule = None
        dlgnp.analyze_symmetry()
        assert len(calls["warning"]) == 1

    def test_full_scan_populates_groups_and_selects_row(self, dlgnp, monkeypatch):
        _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(_symnp, "PointGroupAnalyzer", lambda mol, tolerance: _FakeAnalyzer(tolerance))
        # Run the worker synchronously (in-thread) instead of spawning a real
        # QThread, so the test is deterministic (see repo test-writing gotchas).
        monkeypatch.setattr(_symnp.SymmetryAnalysisWorker, "start", _symnp.SymmetryAnalysisWorker.run)

        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        dlgnp.max_tol_spin.setValue(0.2)
        dlgnp.analyze_symmetry()

        assert dlgnp.calc_btn.isEnabled()
        assert dlgnp.calc_btn.text() == "Analyze (Scan)"
        assert dlgnp.groups_list.count() == 2  # C2v, Cs
        # on_group_selected fires via the real itemSelectionChanged signal.
        assert dlgnp.sym_btn.isEnabled()


# ===========================================================================
# on_group_selected / update_ops_list / _get_op_sort_key / _get_op_label /
# _format_symmetry_symbol — direct calls on the real bound methods
# ===========================================================================


class TestOnGroupSelectedReal:
    def test_selects_group_sorts_ops_and_labels_them(self, dlgnp):
        analyzer = MagicMock()
        analyzer.get_symmetry_operations.return_value = [
            _op(_MAT_S4), _op(_MAT_E), _op(_MAT_C2),
        ]
        dlgnp.group_data = {"Td": {"analyzer": analyzer}}
        dlgnp.groups_list.addItem("Td  (Tol: 0.10 - 1.00 Å)")
        dlgnp.groups_list.setCurrentRow(0)

        assert dlgnp.analyzer is analyzer
        mats = [op.rotation_matrix for op in dlgnp.symmetry_ops]
        assert np.allclose(mats[0], _MAT_E)
        assert np.allclose(mats[1], _MAT_C2)
        assert np.allclose(mats[2], _MAT_S4)
        assert dlgnp.sym_btn.isEnabled()
        assert dlgnp.ops_list.count() == 3
        assert dlgnp.ops_list.item(0).text() == "#1: E (Identity)"
        assert dlgnp.ops_list.item(1).text() == "#2: C₂ (Rotation)"
        assert dlgnp.ops_list.item(2).text() == "#3: S₄ (Improper Rotation)"
        assert "Point Group:" in dlgnp.selected_group_label.text()
        assert "<i>T</i>" in dlgnp.selected_group_label.text()

    def test_no_current_item_is_noop(self, dlgnp):
        dlgnp.on_group_selected()  # nothing selected -> should not raise

    def test_symbol_not_in_group_data(self, dlgnp):
        dlgnp.group_data = {}
        dlgnp.groups_list.addItem("No point groups found.")
        dlgnp.groups_list.setCurrentRow(0)
        assert not dlgnp.sym_btn.isEnabled()


class TestFormatSymmetrySymbolReal:
    def test_c2v(self, dlgnp):
        assert dlgnp._format_symmetry_symbol("C2v") == "<html><i>C</i><sub>2v</sub></html>"

    def test_infinity_axis(self, dlgnp):
        assert dlgnp._format_symmetry_symbol("C*v") == "<html><i>C</i><sub>∞v</sub></html>"

    def test_non_matching_symbol(self, dlgnp):
        assert dlgnp._format_symmetry_symbol("1?") == "1?"


# ===========================================================================
# on_op_selection_changed / _display_single_op_details / visualize_ops /
# _add_op_visualization — real numpy vector math
# ===========================================================================


class TestOnOpSelectionChangedReal:
    def _setup(self, dlgnp, n_ops):
        dlgnp.context.plotter = None  # skip real visualize_ops draw path here
        dlgnp.symmetry_ops = [_op(_MAT_E) for _ in range(n_ops)]
        for i in range(n_ops):
            dlgnp.ops_list.addItem(f"#{i + 1}: E (Identity)")

    def test_zero_selected_clears_details(self, dlgnp):
        self._setup(dlgnp, 3)
        dlgnp.on_op_selection_changed()
        assert dlgnp.op_details.toPlainText() == ""

    def test_one_selected_shows_matrix_details(self, dlgnp):
        self._setup(dlgnp, 3)
        dlgnp.ops_list.item(1).setSelected(True)
        dlgnp.on_op_selection_changed()
        assert "Rotation Matrix" in dlgnp.op_details.toPlainText()

    def test_multiple_selected_shows_count(self, dlgnp):
        self._setup(dlgnp, 3)
        dlgnp.ops_list.item(0).setSelected(True)
        dlgnp.ops_list.item(2).setSelected(True)
        dlgnp.on_op_selection_changed()
        assert "2 operations selected" in dlgnp.op_details.toPlainText()


class TestVisualizeOpsReal:
    # visualize_ops builds real pyvista meshes (Sphere/Line/Disc); with pyvista
    # mocked those calls return MagicMocks and the drawn-actor assertions no
    # longer hold. These tests are only meaningful with a real pyvista/vtk build
    # installed (CI test-gui deliberately omits them), so skip otherwise.
    @pytest.fixture(autouse=True)
    def _need_real_pyvista(self):
        pytest.importorskip("pyvista")
        pytest.importorskip("vtk")

    def test_no_plotter_returns_immediately(self, dlgnp):
        dlgnp.context.plotter = None
        dlgnp.visualize_ops([_op(_MAT_E)])  # should not raise

    def test_empty_ops_clears_and_renders(self, dlgnp):
        plotter = _FakePlotter()
        dlgnp.context.plotter = plotter
        dlgnp.vis_actors = ["stub"]
        dlgnp.visualize_ops([])
        assert dlgnp.vis_actors == []
        assert plotter.removed == ["stub"]

    def test_identity_inversion_rotation_and_mirror_all_draw(self, dlgnp):
        # visualize_ops does `import pyvista as pv` lazily at call time, so
        # whatever real pyvista build is installed is used for mesh creation;
        # only self.context.plotter (our fake) is under test-control here.
        plotter = _FakePlotter()
        dlgnp.context.plotter = plotter
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        ops = [_op(_MAT_E), _op(_MAT_INV), _op(_rot_z(180)), _op(np.diag([1.0, 1.0, -1.0]))]
        dlgnp.visualize_ops(ops)
        # Identity contributes nothing; inversion=sphere(orange), rotation=line(cyan), mirror=disc(magenta).
        colors = [a.kwargs.get("color") for a in plotter.added]
        assert "orange" in colors
        assert "cyan" in colors
        assert "magenta" in colors
        assert len(plotter.added) == 3

    def test_improper_rotation_draws_its_axis(self, dlgnp):
        """S_n is neither a plane nor an inversion; it used to fall through
        both branches undrawn and methane lost all 6 of its S4 operations."""
        plotter = _FakePlotter()
        dlgnp.context.plotter = plotter
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        s4 = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])
        dlgnp.visualize_ops([_op(s4)])
        assert [a.kwargs.get("color") for a in plotter.added] == ["yellow"]
        assert dlgnp.vis_actors == plotter.added

    def test_a_second_selection_replaces_the_previous_actors(self, dlgnp):
        plotter = _FakePlotter()
        dlgnp.context.plotter = plotter
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        dlgnp.visualize_ops([_op(_rot_z(180))])
        first = list(dlgnp.vis_actors)
        dlgnp.visualize_ops([_op(_MAT_INV)])
        assert plotter.removed == first
        assert dlgnp.vis_actors != first

    def test_visualization_without_a_molecule_still_draws(self, dlgnp):
        """No molecule means no centre of mass; the elements fall back to the
        origin rather than raising."""
        plotter = _FakePlotter()
        dlgnp.context.plotter = plotter
        dlgnp.context.current_molecule = None
        dlgnp.visualize_ops([_op(_rot_z(180))])
        assert len(plotter.added) == 1


# ===========================================================================
# symmetrize_structure / update_rdkit_coords
# ===========================================================================


def _mol_pmg_stub(coords, species):
    """Bypasses get_pymatgen_molecule (pymatgen itself stays MagicMock even
    in `_symnp`) with a stand-in exposing real numpy cart_coords/species."""
    arr = np.array(coords, dtype=float)
    masses = np.array([_Species(s).mass for s in species], dtype=float)
    return SimpleNamespace(
        cart_coords=arr,
        species=[_Species(s) for s in species],
        # pymatgen centres on the centre of mass, which is the frame its
        # symmetry operations are expressed in.
        center_of_mass=(arr * masses[:, None]).sum(axis=0) / masses.sum(),
    )


class TestSymmetrizeStructureReal:
    def test_no_analyzer_is_noop(self, dlgnp):
        dlgnp.analyzer = None
        dlgnp.symmetrize_structure()  # should not raise / not touch anything

    def test_no_symmetry_ops_is_noop(self, dlgnp, monkeypatch):
        dlgnp.analyzer = SimpleNamespace(get_symmetry_operations=lambda: [], sch_symbol="C1")
        recorded = []
        monkeypatch.setattr(dlgnp, "update_rdkit_coords", lambda c: recorded.append(c))
        monkeypatch.setattr(dlgnp, "get_pymatgen_molecule", lambda: _mol_pmg_stub(_WATER_COORDS, _WATER_SPECIES))
        dlgnp.symmetrize_structure()
        assert recorded == []

    def test_exact_symmetric_geometry_reproduced(self, dlgnp, monkeypatch):
        _no_block_msgbox(monkeypatch)
        dlgnp.analyzer = SimpleNamespace(
            get_symmetry_operations=lambda: _WATER_OPS, sch_symbol="C2v",
        )
        monkeypatch.setattr(dlgnp, "get_pymatgen_molecule", lambda: _mol_pmg_stub(_WATER_COORDS, _WATER_SPECIES))
        recorded = []
        monkeypatch.setattr(dlgnp, "update_rdkit_coords", lambda c: recorded.append(c))
        dlgnp.symmetrize_structure()
        assert len(recorded) == 1
        assert np.allclose(recorded[0], np.array(_WATER_COORDS), atol=1e-8)

    def test_symmetrizes_without_scipy(self, dlgnp, monkeypatch):
        _no_block_msgbox(monkeypatch)
        monkeypatch.setitem(sys.modules, "scipy", None)
        monkeypatch.setitem(sys.modules, "scipy.optimize", None)
        dlgnp.analyzer = SimpleNamespace(
            get_symmetry_operations=lambda: _WATER_OPS, sch_symbol="C2v",
        )
        monkeypatch.setattr(dlgnp, "get_pymatgen_molecule", lambda: _mol_pmg_stub(_WATER_COORDS, _WATER_SPECIES))
        recorded = []
        monkeypatch.setattr(dlgnp, "update_rdkit_coords", lambda c: recorded.append(c))
        dlgnp.symmetrize_structure()
        assert len(recorded) == 1
        assert np.allclose(recorded[0], np.array(_WATER_COORDS), atol=1e-8)

    def test_missing_scipy_is_logged_not_printed(self, dlgnp, monkeypatch, caplog):
        """A GUI plugin printing to stdout tells the user nothing."""
        _no_block_msgbox(monkeypatch)
        monkeypatch.setitem(sys.modules, "scipy", None)
        monkeypatch.setitem(sys.modules, "scipy.optimize", None)
        dlgnp.analyzer = SimpleNamespace(
            get_symmetry_operations=lambda: _WATER_OPS, sch_symbol="C2v",
        )
        monkeypatch.setattr(dlgnp, "get_pymatgen_molecule", lambda: _mol_pmg_stub(_WATER_COORDS, _WATER_SPECIES))
        monkeypatch.setattr(dlgnp, "update_rdkit_coords", MagicMock())
        with caplog.at_level(logging.WARNING):
            dlgnp.symmetrize_structure()
        assert any("scipy" in r.message for r in caplog.records)

    def test_greedy_fallback_flags_an_atom_it_cannot_map(self, dlgnp, monkeypatch):
        """Without scipy an atom with no partner within 1 A is left unmapped,
        which averages in the wrong site -- it must reach the confirmation."""
        calls = _no_block_msgbox(monkeypatch)
        monkeypatch.setitem(sys.modules, "scipy", None)
        monkeypatch.setitem(sys.modules, "scipy.optimize", None)
        spread = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (0.0, 3.0, 0.0)]
        dlgnp.analyzer = SimpleNamespace(
            get_symmetry_operations=lambda: [_SymOp(np.eye(3)), _SymOp(_rot_z(180))],
            sch_symbol="C2",
        )
        monkeypatch.setattr(dlgnp, "get_pymatgen_molecule", lambda: _mol_pmg_stub(spread, ["H", "H", "H"]))
        monkeypatch.setattr(dlgnp, "update_rdkit_coords", MagicMock())
        dlgnp.symmetrize_structure()
        assert len(calls["question"]) == 1

    def test_mismatched_op_on_chiral_geometry_warns(self, dlgnp, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        chiral_coords = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.5), (0.0, 2.0, 1.5), (-2.0, -0.8, 2.0)]
        chiral_species = ["C", "H", "H", "H"]
        bad_ops = [_SymOp(np.eye(3)), _SymOp(_rot_z(180))]
        dlgnp.analyzer = SimpleNamespace(get_symmetry_operations=lambda: bad_ops, sch_symbol="C1")
        monkeypatch.setattr(dlgnp, "get_pymatgen_molecule", lambda: _mol_pmg_stub(chiral_coords, chiral_species))
        applied = MagicMock()
        monkeypatch.setattr(dlgnp, "update_rdkit_coords", applied)
        dlgnp.symmetrize_structure()
        # Unmapped atoms now prompt before overwriting the geometry, rather
        # than warning and applying a possibly scrambled structure anyway.
        assert len(calls["question"]) == 1
        applied.assert_called_once()  # the stub answers Yes

    def test_declining_the_prompt_leaves_the_geometry_alone(self, dlgnp, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(
            _symnp.QMessageBox, "question",
            staticmethod(
                lambda *a, **kw: (
                    calls["question"].append(a),
                    _symnp.QMessageBox.StandardButton.No,
                )[1]
            ),
        )
        chiral_coords = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.5), (0.0, 2.0, 1.5), (-2.0, -0.8, 2.0)]
        dlgnp.analyzer = SimpleNamespace(
            get_symmetry_operations=lambda: [_SymOp(np.eye(3)), _SymOp(_rot_z(180))],
            sch_symbol="C1",
        )
        monkeypatch.setattr(
            dlgnp, "get_pymatgen_molecule",
            lambda: _mol_pmg_stub(chiral_coords, ["C", "H", "H", "H"]),
        )
        applied = MagicMock()
        monkeypatch.setattr(dlgnp, "update_rdkit_coords", applied)
        dlgnp.symmetrize_structure()
        assert len(calls["question"]) == 1
        applied.assert_not_called()


class TestUpdateRdkitCoordsReal:
    def test_no_molecule_is_noop(self, dlgnp):
        dlgnp.context.current_molecule = None
        dlgnp.update_rdkit_coords(np.zeros((0, 3)))  # should not raise

    def test_updates_conformer_and_refreshes_view(self, dlgnp, monkeypatch):
        # Point3D comes from mocked rdkit even in `_symnp`; swap in a real
        # x/y/z-capturing stand-in so the written-back positions are assertable.
        monkeypatch.setattr(_symnp, "Point3D", _P3)
        mol = _NPMol([8, 1], _WATER_COORDS[:2])
        dlgnp.context.current_molecule = mol
        refreshed = []
        dlgnp.context.refresh_3d_view = lambda: refreshed.append(True)
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        dlgnp.update_rdkit_coords(coords)
        assert list(mol._conf.GetAtomPosition(0)) == [1.0, 2.0, 3.0]
        assert list(mol._conf.GetAtomPosition(1)) == [4.0, 5.0, 6.0]
        assert dlgnp.context.current_mol is mol
        assert refreshed == [True]


# ===========================================================================
# closeEvent
# ===========================================================================


class TestCloseEventReal:
    def test_close_resets_ui_state(self, dlgnp):
        dlgnp.group_data = {"x": 1}
        dlgnp.symmetry_ops = [1, 2]
        dlgnp.analyzer = object()
        dlgnp.close()
        assert dlgnp.group_data == {}
        assert dlgnp.symmetry_ops == []
        assert dlgnp.analyzer is None
        assert not dlgnp.sym_btn.isEnabled()

    def test_close_stops_running_worker(self, dlgnp):
        worker = MagicMock()
        worker.isRunning.side_effect = [True, False]
        dlgnp.worker = worker
        dlgnp.close()
        worker.abort.assert_called_once()
        worker.wait.assert_called_once_with(5000)
        worker.terminate.assert_not_called()

    def test_escape_runs_the_same_cleanup_as_close(self, dlgnp):
        """Issue #10: Esc rejects a QDialog without delivering a close event,
        so the 3D overlay and the scan thread used to survive it."""
        plotter = _FakePlotter()
        dlgnp.context.plotter = plotter
        dlgnp.vis_actors = ["axis", "plane"]
        dlgnp.group_data = {"D2": {}}
        dlgnp.analyzer = object()
        worker = MagicMock()
        worker.isRunning.side_effect = [True, False]
        dlgnp.worker = worker

        dlgnp.show()
        QTest.keyClick(dlgnp, Qt.Key.Key_Escape)

        assert not dlgnp.isVisible()
        assert plotter.removed == ["axis", "plane"]
        assert dlgnp.group_data == {}
        assert dlgnp.analyzer is None
        worker.abort.assert_called_once()

    def test_close_after_escape_does_not_recurse(self, dlgnp):
        """QDialog.closeEvent routes back through reject(); the two paths must
        not call each other (ORCA Result Analyzer v3.10.2 recursed here)."""
        dlgnp.show()
        QTest.keyClick(dlgnp, Qt.Key.Key_Escape)
        dlgnp.close()  # must simply return
        assert not dlgnp.isVisible()

    def test_close_terminates_worker_still_running_after_wait(self, dlgnp):
        worker = MagicMock()
        worker.isRunning.return_value = True
        dlgnp.worker = worker
        dlgnp.close()
        worker.terminate.assert_called_once()

    def test_close_removes_visualization_actors(self, dlgnp):
        plotter = _FakePlotter()
        dlgnp.context.plotter = plotter
        dlgnp.vis_actors = ["a1", "a2"]
        dlgnp.close()
        assert plotter.removed == ["a1", "a2"]


# ===========================================================================
# Analysis freshness: fingerprint guard + document reset
# ===========================================================================


class TestAnalysisFreshness:
    def test_document_reset_handler_is_registered(self, dlgnp):
        """File->New must not leave the previous molecule's groups on screen."""
        handlers = [
            c.args[0]
            for c in dlgnp.context.register_document_reset_handler.call_args_list
        ]
        assert dlgnp._invalidate_analysis in handlers

    def test_invalidate_clears_results_and_overlay(self, dlgnp):
        plotter = _FakePlotter()
        dlgnp.context.plotter = plotter
        dlgnp.vis_actors = ["axis"]
        dlgnp.group_data = {"C2": {}}
        dlgnp.analyzer = object()
        dlgnp.symmetry_ops = [1]
        dlgnp._scanned_fingerprint = "something"
        dlgnp.sym_btn.setEnabled(True)

        dlgnp._invalidate_analysis()

        assert plotter.removed == ["axis"]
        assert dlgnp.group_data == {}
        assert dlgnp.symmetry_ops == []
        assert dlgnp.analyzer is None
        assert dlgnp._scanned_fingerprint is None
        assert not dlgnp.sym_btn.isEnabled()
        assert dlgnp.selected_group_label.text() == "Point Group: -"

    def test_fingerprint_tracks_coordinates(self, dlgnp):
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        first = dlgnp._molecule_fingerprint()
        moved = [(0.0, 0.0, 0.5)] + list(_WATER_COORDS[1:])
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], moved)
        assert dlgnp._molecule_fingerprint() != first

    def test_fingerprint_tracks_elements(self, dlgnp):
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        first = dlgnp._molecule_fingerprint()
        dlgnp.context.current_molecule = _NPMol([16, 1, 1], _WATER_COORDS)
        assert dlgnp._molecule_fingerprint() != first

    def test_fingerprint_is_none_without_a_molecule(self, dlgnp):
        dlgnp.context.current_molecule = None
        assert dlgnp._molecule_fingerprint() is None

    def test_scan_records_the_fingerprint_it_ran_against(self, dlgnp, monkeypatch):
        _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(
            _symnp, "PointGroupAnalyzer", lambda mol, tolerance: _FakeAnalyzer(tolerance)
        )
        monkeypatch.setattr(
            _symnp.SymmetryAnalysisWorker, "start", _symnp.SymmetryAnalysisWorker.run
        )
        mol = _NPMol([8, 1, 1], _WATER_COORDS)
        dlgnp.context.current_molecule = mol
        dlgnp.analyze_symmetry()
        assert dlgnp._scanned_fingerprint == dlgnp._molecule_fingerprint()

    def test_symmetrize_refuses_a_molecule_the_scan_never_saw(
        self, dlgnp, monkeypatch
    ):
        """Issue #11: water's C2v used to be applied to whatever was loaded."""
        calls = _no_block_msgbox(monkeypatch)
        dlgnp.analyzer = SimpleNamespace(
            sch_symbol="C2v", get_symmetry_operations=lambda: _WATER_OPS
        )
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        dlgnp._scanned_fingerprint = "a molecule that is no longer loaded"
        dlgnp.context.push_undo_checkpoint.reset_mock()

        dlgnp.symmetrize_structure()

        assert len(calls["warning"]) == 1
        assert not calls["information"]
        dlgnp.context.push_undo_checkpoint.assert_not_called()
        assert dlgnp.analyzer is None
        assert not dlgnp.sym_btn.isEnabled()


# ===========================================================================
# Atomless molecule / placeholder row
# ===========================================================================


class TestDegenerateInputs:
    def test_atomless_molecule_is_not_handed_to_pymatgen(self, dlgnp):
        """pymatgen's Molecule([], []) raises IndexError."""
        dlgnp.context.current_molecule = _NPMol([], [])
        assert dlgnp.get_pymatgen_molecule() is None

    def test_atomless_molecule_warns_instead_of_scanning(self, dlgnp, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        dlgnp.context.current_molecule = _NPMol([], [])
        dlgnp.analyze_symmetry()
        assert len(calls["warning"]) == 1
        assert dlgnp.worker is None

    def test_selecting_the_placeholder_row_is_ignored(self, dlgnp, monkeypatch):
        _no_block_msgbox(monkeypatch)
        dlgnp.on_analysis_finished({}, False)
        dlgnp.groups_list.setCurrentRow(0)
        assert dlgnp.groups_list.item(0).text() == "No point groups found."
        assert dlgnp.selected_group_label.text() == "Point Group: -"
        assert not dlgnp.sym_btn.isEnabled()


# ===========================================================================
# Stopping a running scan
# ===========================================================================


class TestStopScan:
    def test_button_starts_a_scan_when_idle(self, dlgnp, monkeypatch):
        started = []
        monkeypatch.setattr(dlgnp, "analyze_symmetry", lambda: started.append(1))
        dlgnp.worker = None
        dlgnp.on_analyze_clicked()
        assert started == [1]

    def test_button_stops_the_scan_while_it_runs(self, dlgnp, monkeypatch):
        monkeypatch.setattr(
            dlgnp, "analyze_symmetry", lambda: pytest.fail("started a second scan")
        )
        worker = MagicMock()
        worker.isRunning.return_value = True
        dlgnp.worker = worker

        dlgnp.on_analyze_clicked()

        worker.abort.assert_called_once()
        assert dlgnp.calc_btn.text() == "Stopping..."
        assert not dlgnp.calc_btn.isEnabled()

    def test_stop_is_a_no_op_without_a_running_scan(self, dlgnp):
        worker = MagicMock()
        worker.isRunning.return_value = False
        dlgnp.worker = worker
        dlgnp.stop_analysis()
        worker.abort.assert_not_called()
        assert dlgnp.calc_btn.text() == "Analyze (Scan)"

    def test_button_offers_stop_while_scanning(self, dlgnp, monkeypatch):
        _no_block_msgbox(monkeypatch)
        monkeypatch.setattr(_symnp.SymmetryAnalysisWorker, "start", lambda self: None)
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        dlgnp.analyze_symmetry()
        assert dlgnp.calc_btn.text() == "Stop"
        assert dlgnp.calc_btn.isEnabled()

    def test_a_stopped_scan_still_reports_what_it_found(self, dlgnp, monkeypatch):
        """Aborting is not cancelling: the partial bands are still useful."""
        _no_block_msgbox(monkeypatch)

        def pga(mol, tolerance):
            worker_ref["w"]._abort = True
            return _FakeAnalyzer(tolerance)

        worker_ref = {}
        real_start = _symnp.SymmetryAnalysisWorker.run

        def start(self):
            worker_ref["w"] = self
            real_start(self)

        monkeypatch.setattr(_symnp, "PointGroupAnalyzer", pga)
        monkeypatch.setattr(_symnp.SymmetryAnalysisWorker, "start", start)
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        dlgnp.analyze_symmetry()

        assert dlgnp.groups_list.count() >= 1
        assert dlgnp.calc_btn.text() == "Analyze (Scan)"
        assert dlgnp.calc_btn.isEnabled()

    def test_tolerance_range_reaches_the_worker(self, dlgnp, monkeypatch):
        _no_block_msgbox(monkeypatch)
        captured = {}
        monkeypatch.setattr(
            _symnp.SymmetryAnalysisWorker,
            "start",
            lambda self: captured.update(lo=self.min_tol, hi=self.max_tol),
        )
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        dlgnp.min_tol_spin.setValue(0.005)
        dlgnp.max_tol_spin.setValue(0.3)
        dlgnp.analyze_symmetry()
        assert captured == pytest.approx({"lo": 0.005, "hi": 0.3})

    def test_no_step_control_is_offered(self, dlgnp):
        """The refinement pass makes the result step-independent, so a Step
        setting only looked like it did nothing."""
        assert not hasattr(dlgnp, "step_spin")


# ===========================================================================
# pymatgen names groups whose operations it cannot produce
# ===========================================================================


class TestSymbolReconciliation:
    def _analyzer(self, matrices):
        return SimpleNamespace(
            get_symmetry_operations=lambda: [
                SimpleNamespace(rotation_matrix=np.array(m, dtype=float))
                for m in matrices
            ]
        )

    def test_a_scan_reports_c2_when_only_two_operations_exist(
        self, dlgnp, monkeypatch
    ):
        """pymatgen 2025.10.7 names tetramethylhydrazine D2 while its own
        get_symmetry_operations() returns just E and one C2."""
        _no_block_msgbox(monkeypatch)
        eye = np.eye(3)
        c2z = np.diag([-1.0, -1.0, 1.0])

        class _Claims:
            sch_symbol = "D2"

            def get_symmetry_operations(self):
                return [
                    SimpleNamespace(rotation_matrix=eye),
                    SimpleNamespace(rotation_matrix=c2z),
                ]

        monkeypatch.setattr(_symnp, "PointGroupAnalyzer", lambda mol, tolerance: _Claims())
        monkeypatch.setattr(
            _symnp.SymmetryAnalysisWorker, "start", _symnp.SymmetryAnalysisWorker.run
        )
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        dlgnp.analyze_symmetry()

        rows = [dlgnp.groups_list.item(i).text() for i in range(dlgnp.groups_list.count())]
        assert rows and all(r.startswith("C2 ") for r in rows), rows
        assert "D2" not in " ".join(rows)

    def test_a_real_d2_survives_the_scan(self, dlgnp, monkeypatch):
        _no_block_msgbox(monkeypatch)
        mats = [
            np.eye(3),
            np.diag([-1.0, -1.0, 1.0]),
            np.diag([1.0, -1.0, -1.0]),
            np.diag([-1.0, 1.0, -1.0]),
        ]

        class _Claims:
            sch_symbol = "D2"

            def get_symmetry_operations(self):
                return [SimpleNamespace(rotation_matrix=m) for m in mats]

        monkeypatch.setattr(_symnp, "PointGroupAnalyzer", lambda mol, tolerance: _Claims())
        monkeypatch.setattr(
            _symnp.SymmetryAnalysisWorker, "start", _symnp.SymmetryAnalysisWorker.run
        )
        dlgnp.context.current_molecule = _NPMol([8, 1, 1], _WATER_COORDS)
        dlgnp.analyze_symmetry()

        rows = [dlgnp.groups_list.item(i).text() for i in range(dlgnp.groups_list.count())]
        assert rows and all(r.startswith("D2 ") for r in rows), rows


# ===========================================================================
# _group_from_operations on the real worker (executes the module's own lines)
# ===========================================================================


def _ops(*matrices):
    return [
        SimpleNamespace(rotation_matrix=np.array(m, dtype=float)) for m in matrices
    ]


_EYE = np.eye(3)
_C2Z = np.diag([-1.0, -1.0, 1.0])
_C2X = np.diag([1.0, -1.0, -1.0])
_C2Y = np.diag([-1.0, 1.0, -1.0])
_SGZ = np.diag([1.0, 1.0, -1.0])
_SGX = np.diag([-1.0, 1.0, 1.0])
_SGY = np.diag([1.0, -1.0, 1.0])
_INVERSION = -np.eye(3)
_S4Z = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])
_C3Z = np.array(
    [[-0.5, -(3**0.5) / 2, 0.0], [(3**0.5) / 2, -0.5, 0.0], [0.0, 0.0, 1.0]]
)
_C4Z = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])


class TestGroupFromOperationsReal:
    @pytest.fixture
    def worker(self):
        return _symnp.SymmetryAnalysisWorker(MagicMock(), 0.01, 0.5)

    @pytest.mark.parametrize(
        "matrices,expected",
        [
            ((_EYE,), "C1"),
            ((_EYE, _SGZ), "Cs"),
            ((_EYE, _INVERSION), "Ci"),
            ((_EYE, _C2Z), "C2"),
            ((_EYE, _C2X, _C2Y, _C2Z), "D2"),
            ((_EYE, _C2Z, _SGX, _SGY), "C2v"),
            ((_EYE, _C2Z, _SGZ, _INVERSION), "C2h"),
            ((_EYE, _C2Z, _S4Z, _S4Z.T), "S4"),
            ((_EYE, _C3Z, _C3Z.T), "C3"),
            ((_EYE, _C4Z, _C4Z.T, _C2Z), "C4"),
            ((_EYE, _C2X, _C2Y, _C2Z, _SGZ, _SGX, _SGY, _INVERSION), "D2h"),
        ],
    )
    def test_symbol_derived_from_operations(self, worker, matrices, expected):
        assert worker._group_from_operations(_ops(*matrices)) == expected

    def test_cubic_is_left_to_pymatgen(self, worker):
        body = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        other = np.array([[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        assert worker._group_from_operations(_ops(_EYE, body, body.T, other, other.T)) is None

    def test_a_repeated_axis_counts_once(self, worker):
        """A rotation and its inverse share an axis."""
        assert worker._group_from_operations(_ops(_EYE, _C2Z, _C2Z)) == "C2"

    def test_reconcile_downgrades_an_unsupported_name(self, worker):
        analyzer = SimpleNamespace(get_symmetry_operations=lambda: _ops(_EYE, _C2Z))
        assert worker._reconcile_symbol("D2", analyzer) == "C2"

    def test_reconcile_keeps_a_supported_name(self, worker):
        analyzer = SimpleNamespace(
            get_symmetry_operations=lambda: _ops(_EYE, _C2X, _C2Y, _C2Z)
        )
        assert worker._reconcile_symbol("D2", analyzer) == "D2"

    def test_reconcile_ignores_groups_outside_the_table(self, worker):
        analyzer = SimpleNamespace(get_symmetry_operations=lambda: _ops(_EYE))
        assert worker._reconcile_symbol("Td", analyzer) == "Td"

    def test_reconcile_survives_an_analyzer_without_operations(self, worker):
        assert worker._reconcile_symbol("D2", SimpleNamespace()) == "D2"
