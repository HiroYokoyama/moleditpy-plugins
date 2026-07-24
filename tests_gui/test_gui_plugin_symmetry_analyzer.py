"""
Headless GUI tests for the Symmetry Analyzer plugin.

Covers: SymmetryAnalysisPlugin.

Real PyQt6 is used (QT_QPA_PLATFORM=offscreen); chemistry libs are MagicMocked.
Run via: python tests_gui/run_gui_tests.py tests_gui/test_gui_plugin_symmetry_analyzer.py
"""

from __future__ import annotations

import contextlib
import math
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

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
        assert dlg.max_tol_spin.value() == pytest.approx(1.0)

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
    def __init__(self, z):
        self._z = z

    def GetAtomicNum(self):
        return self._z


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
    def __init__(self, symbol):
        self.symbol = symbol


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
        worker.finished.connect(lambda gd, fa: received.update(group_data=gd, found_any=fa))
        worker.run()
        assert received["found_any"] is True
        assert set(received["group_data"].keys()) == {"C2v", "Cs"}
        # tol 0.0, 0.05, 0.1 -> C2v; 0.15, 0.2 -> Cs (arange step 0.05)
        assert received["group_data"]["C2v"]["tols"] == pytest.approx([0.0, 0.05, 0.1])
        assert received["group_data"]["Cs"]["tols"] == pytest.approx([0.15, 0.2])

    def test_run_tolerates_analyzer_exceptions(self, monkeypatch):
        def _boom(mol, tolerance):
            raise ValueError("bad tolerance")

        monkeypatch.setattr(_symnp, "PointGroupAnalyzer", _boom)
        worker = _symnp.SymmetryAnalysisWorker(mol_pmg=MagicMock(), min_tol=0.0, max_tol=0.1)
        received = {}
        worker.finished.connect(lambda gd, fa: received.update(group_data=gd, found_any=fa))
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


# ===========================================================================
# symmetrize_structure / update_rdkit_coords
# ===========================================================================


def _mol_pmg_stub(coords, species):
    """Bypasses get_pymatgen_molecule (pymatgen itself stays MagicMock even
    in `_symnp`) with a stand-in exposing real numpy cart_coords/species."""
    return SimpleNamespace(
        cart_coords=np.array(coords, dtype=float),
        species=[_Species(s) for s in species],
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

    def test_mismatched_op_on_chiral_geometry_warns(self, dlgnp, monkeypatch):
        calls = _no_block_msgbox(monkeypatch)
        chiral_coords = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.5), (0.0, 2.0, 1.5), (-2.0, -0.8, 2.0)]
        chiral_species = ["C", "H", "H", "H"]
        bad_ops = [_SymOp(np.eye(3)), _SymOp(_rot_z(180))]
        dlgnp.analyzer = SimpleNamespace(get_symmetry_operations=lambda: bad_ops, sch_symbol="C1")
        monkeypatch.setattr(dlgnp, "get_pymatgen_molecule", lambda: _mol_pmg_stub(chiral_coords, chiral_species))
        monkeypatch.setattr(dlgnp, "update_rdkit_coords", MagicMock())
        dlgnp.symmetrize_structure()
        assert len(calls["warning"]) == 1


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
        worker.quit.assert_called_once()
        worker.wait.assert_called_once_with(1000)
        worker.terminate.assert_not_called()

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
