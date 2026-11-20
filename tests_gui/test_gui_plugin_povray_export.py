"""
Headless GUI tests for the POV-Ray Export plugin.

Covers module constants, initialize(), and run().

All tests use real PySide6 (QT_QPA_PLATFORM=offscreen).
Chemistry libraries are mocked via mock_chemistry_imports().
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from conftest import load_plugin_for_gui, mock_chemistry_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"

POVRAY_PATH = PLUGINS_DIR / "POV-Ray_Export" / "povray_export.py"

with mock_chemistry_imports():
    _povray = load_plugin_for_gui(POVRAY_PATH)


def _make_ctx() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# POV-Ray Export — module-level and initialize()/run() tests
# ===========================================================================


class TestPovRayExportModule:
    """Module constants, initialize(), and run() for POV-Ray Export."""

    def test_plugin_name(self):
        assert _povray.PLUGIN_NAME == "POV-Ray Export"

    def test_plugin_version_is_string(self):
        assert isinstance(_povray.PLUGIN_VERSION, str)

    def test_plugin_context_initially_none(self):
        assert _povray.PLUGIN_CONTEXT is None

    def test_initialize_registers_export_action(self):
        ctx = _make_ctx()
        _povray.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_povray(self):
        ctx = _make_ctx()
        _povray.initialize(ctx)
        label = ctx.add_export_action.call_args.args[0]
        assert "POV" in label

    def test_initialize_callback_is_callable(self):
        ctx = _make_ctx()
        _povray.initialize(ctx)
        cb = ctx.add_export_action.call_args.args[1]
        assert callable(cb)

    def test_run_noop_without_plugin_manager(self):
        mw = MagicMock(spec=[])
        _povray.run(mw)  # must not raise

    def test_run_noop_when_context_is_none(self):
        _povray.PLUGIN_CONTEXT = None
        mw = MagicMock()
        _povray.run(mw)  # PLUGIN_CONTEXT is None → early return

    def test_run_calls_export_to_povray_when_context_set(self, monkeypatch):
        ctx = _make_ctx()
        _povray.PLUGIN_CONTEXT = ctx
        export_mock = MagicMock()
        monkeypatch.setattr(_povray, "export_to_povray", export_mock)
        mw = MagicMock()  # has plugin_manager (auto-attribute)
        _povray.run(mw)
        export_mock.assert_called_once_with(ctx)


# ===========================================================================
# export_to_povray — control flow (guards / cancel / save / error)
# ===========================================================================


def _mol(n_atoms=3, n_confs=1):
    mol = MagicMock()
    mol.GetNumAtoms.return_value = n_atoms
    mol.GetNumConformers.return_value = n_confs
    return mol


class TestPovRayExportFlow:
    @pytest.fixture(autouse=True)
    def _qt(self, qapp, monkeypatch):
        # export_to_povray imports these from PySide6.QtWidgets at call time.
        from PySide6.QtWidgets import QFileDialog, QMessageBox

        self.warning = MagicMock()
        self.info = MagicMock()
        self.critical = MagicMock()
        monkeypatch.setattr(QMessageBox, "warning", self.warning)
        monkeypatch.setattr(QMessageBox, "information", self.info)
        monkeypatch.setattr(QMessageBox, "critical", self.critical)
        self._save_ret = ("", "")
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: self._save_ret
        )
        # Isolate the flow from the heavy geometry generator.
        monkeypatch.setattr(
            _povray, "generate_povray_scene", MagicMock(return_value=("SCENE", 800, 600))
        )

    def _ctx(self, mol):
        ctx = MagicMock()
        ctx.get_main_window.return_value = MagicMock()
        ctx.current_molecule = mol
        _povray.PLUGIN_CONTEXT = ctx
        return ctx

    def test_no_molecule_warns_and_skips_generate(self):
        ctx = self._ctx(None)
        _povray.export_to_povray(ctx)
        self.warning.assert_called_once()
        _povray.generate_povray_scene.assert_not_called()

    def test_empty_molecule_warns(self):
        ctx = self._ctx(_mol(n_atoms=0))
        _povray.export_to_povray(ctx)
        self.warning.assert_called_once()

    def test_no_conformer_warns_about_3d(self):
        ctx = self._ctx(_mol(n_confs=0))
        _povray.export_to_povray(ctx)
        assert "3D" in self.warning.call_args.args[1]

    def test_cancel_writes_nothing(self, tmp_path):
        ctx = self._ctx(_mol())
        self._save_ret = ("", "")
        _povray.export_to_povray(ctx)
        _povray.generate_povray_scene.assert_not_called()
        assert list(tmp_path.iterdir()) == []

    def test_save_appends_pov_extension_and_writes(self, tmp_path):
        ctx = self._ctx(_mol())
        target = tmp_path / "scene"
        self._save_ret = (str(target), "")
        _povray.export_to_povray(ctx)
        written = tmp_path / "scene.pov"
        assert written.exists()
        assert written.read_text(encoding="utf-8") == "SCENE"
        self.info.assert_called_once()
        ctx.show_status_message.assert_called()

    def test_generation_error_shows_critical(self, tmp_path):
        ctx = self._ctx(_mol())
        self._save_ret = (str(tmp_path / "x.pov"), "")
        _povray.generate_povray_scene.side_effect = RuntimeError("boom")
        _povray.export_to_povray(ctx)
        self.critical.assert_called_once()

    def test_bad_current_file_path_falls_back_to_default_name(self, tmp_path):
        # current_file_path is a non-str value -> os.path.dirname/splitext raise
        # inside the try block, exercising the silenced except branch.
        ctx = self._ctx(_mol())
        ctx.get_main_window.return_value.init_manager.current_file_path = 12345
        target = tmp_path / "scene"
        self._save_ret = (str(target), "")
        _povray.export_to_povray(ctx)
        assert (tmp_path / "scene.pov").exists()


# ===========================================================================
# Regression: dative / zero bonds through the real export, with real RDKit
# ===========================================================================

# H3N->BH3 core as MoleditPy reads it from a V3000 MOL file: bond type 9 is a
# coordination (DATIVE) bond.
_DATIVE_MOLBLOCK = """
  MoleditPy

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.0 0.0 0.0 0
M  V30 2 B 1.6 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 9 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
"""


def _real_mol(kind):
    """N-B molecule with one *kind* bond and a 3D conformer, in real RDKit."""
    pytest.importorskip("numpy")
    Chem = pytest.importorskip("rdkit.Chem")
    if kind == "DATIVE":
        mol = Chem.MolFromMolBlock(_DATIVE_MOLBLOCK, removeHs=False)
    else:
        from rdkit.Geometry import Point3D

        rw = Chem.RWMol()
        rw.AddAtom(Chem.Atom("N"))
        rw.AddAtom(Chem.Atom("B"))
        rw.AddBond(0, 1, getattr(Chem.BondType, kind))
        conf = Chem.Conformer(2)
        conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
        conf.SetAtomPosition(1, Point3D(1.6, 0.0, 0.0))
        mol = rw.GetMol()
        mol.AddConformer(conf, assignId=True)
    assert mol.GetBondWithIdx(0).GetBondType() == getattr(Chem.BondType, kind)
    return mol


def _export_mw(color_map=None):
    return SimpleNamespace(
        init_manager=SimpleNamespace(settings={}, current_file_path=None),
        view_3d_manager=SimpleNamespace(
            _3d_color_map=color_map or {},
            current_3d_style="ball_and_stick",
            glyph_source=None,
        ),
        plotter=None,
    )


class TestPovRayExportDativeBonds:
    @pytest.fixture(autouse=True)
    def _qt(self, qapp, monkeypatch):
        from PySide6.QtWidgets import QFileDialog, QMessageBox

        self.critical = MagicMock()
        monkeypatch.setattr(QMessageBox, "information", MagicMock())
        monkeypatch.setattr(QMessageBox, "warning", MagicMock())
        monkeypatch.setattr(QMessageBox, "critical", self.critical)
        self._save_ret = ("", "")
        monkeypatch.setattr(
            QFileDialog, "getSaveFileName", lambda *a, **k: self._save_ret
        )

    def _export(self, tmp_path, mol, color_map=None):
        ctx = MagicMock()
        ctx.get_main_window.return_value = _export_mw(color_map)
        ctx.current_molecule = mol
        _povray.PLUGIN_CONTEXT = ctx
        self._save_ret = (str(tmp_path / "scene.pov"), "")
        _povray.export_to_povray(ctx)
        self.critical.assert_not_called()
        return (tmp_path / "scene.pov").read_text(encoding="utf-8")

    @pytest.mark.parametrize("kind", ["DATIVE", "ZERO"])
    def test_bond_exported_as_single_cylinder(self, tmp_path, kind):
        scene = self._export(tmp_path, _real_mol(kind))
        assert scene.count("cylinder {") == 1
        assert "// Bond 0 (fallback)" in scene

    def test_dative_bond_uses_host_split_colors(self, tmp_path):
        scene = self._export(
            tmp_path,
            _real_mol("DATIVE"),
            color_map={"bond_0_start": (255, 0, 0), "bond_0_end": (0, 0, 255)},
        )
        assert scene.count("cylinder {") == 2
        assert "// Bond 0 start" in scene
        assert "// Bond 0 end" in scene
