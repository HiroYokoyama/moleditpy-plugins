"""
Headless GUI tests for the Blender Export plugin.

Covers module constants, initialize(), and run() (no dialog class; export
logic is callable-only).

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

BLENDER_PATH = PLUGINS_DIR / "Blender_Export" / "blender_export.py"

with mock_chemistry_imports():
    _blender = load_plugin_for_gui(BLENDER_PATH)


def _make_ctx() -> MagicMock:
    ctx = MagicMock()
    ctx.get_main_window.return_value = None
    ctx.current_molecule = None
    return ctx


# ===========================================================================
# Blender Export — module-level and initialize()/run() tests
# (no dialog class; export logic is callable-only)
# ===========================================================================


class TestBlenderExportModule:
    """Module constants, initialize(), and run() for Blender Export."""

    def test_plugin_name(self):
        assert _blender.PLUGIN_NAME == "Blender Export"

    def test_plugin_version_is_string(self):
        assert isinstance(_blender.PLUGIN_VERSION, str)

    def test_plugin_context_initially_none(self):
        assert _blender.PLUGIN_CONTEXT is None

    def test_initialize_registers_export_action(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        ctx.add_export_action.assert_called_once()

    def test_initialize_export_label_contains_blender(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        label = ctx.add_export_action.call_args.args[0]
        assert "Blender" in label

    def test_initialize_callback_is_callable(self):
        ctx = _make_ctx()
        _blender.initialize(ctx)
        cb = ctx.add_export_action.call_args.args[1]
        assert callable(cb)

    def test_run_noop_without_plugin_manager(self):
        mw = MagicMock(spec=[])  # no plugin_manager attribute
        _blender.run(mw)  # must not raise

    def test_run_noop_when_context_is_none(self):
        _blender.PLUGIN_CONTEXT = None
        mw = MagicMock()  # has plugin_manager (MagicMock auto-attribute)
        _blender.run(mw)  # PLUGIN_CONTEXT is None → early return, no error


# ===========================================================================
# export_to_blender — control flow (guards / cancel / save / error)
# ===========================================================================


def _mol(n_atoms=3, n_confs=1):
    mol = MagicMock()
    mol.GetNumAtoms.return_value = n_atoms
    mol.GetNumConformers.return_value = n_confs
    return mol


class TestBlenderExportFlow:
    @pytest.fixture(autouse=True)
    def _qt(self, qapp, monkeypatch):
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
        monkeypatch.setattr(
            _blender,
            "generate_blender_script",
            MagicMock(return_value="# blender script"),
        )

    def _ctx(self, mol):
        ctx = MagicMock()
        ctx.get_main_window.return_value = MagicMock()
        ctx.current_molecule = mol
        _blender.PLUGIN_CONTEXT = ctx
        return ctx

    def test_no_molecule_warns_and_skips_generate(self):
        ctx = self._ctx(None)
        _blender.export_to_blender(ctx)
        self.warning.assert_called_once()
        _blender.generate_blender_script.assert_not_called()

    def test_empty_molecule_warns(self):
        ctx = self._ctx(_mol(n_atoms=0))
        _blender.export_to_blender(ctx)
        self.warning.assert_called_once()

    def test_no_conformer_warns_about_3d(self):
        ctx = self._ctx(_mol(n_confs=0))
        _blender.export_to_blender(ctx)
        assert "3D" in self.warning.call_args.args[1]

    def test_cancel_writes_nothing(self, tmp_path):
        ctx = self._ctx(_mol())
        self._save_ret = ("", "")
        _blender.export_to_blender(ctx)
        _blender.generate_blender_script.assert_not_called()
        assert list(tmp_path.iterdir()) == []

    def test_save_appends_py_extension_and_writes(self, tmp_path):
        ctx = self._ctx(_mol())
        target = tmp_path / "scene"
        self._save_ret = (str(target), "")
        _blender.export_to_blender(ctx)
        written = tmp_path / "scene.py"
        assert written.exists()
        assert written.read_text(encoding="utf-8") == "# blender script"
        self.info.assert_called_once()
        ctx.show_status_message.assert_called()

    def test_generation_error_shows_critical(self, tmp_path):
        ctx = self._ctx(_mol())
        self._save_ret = (str(tmp_path / "x.py"), "")
        _blender.generate_blender_script.side_effect = RuntimeError("boom")
        _blender.export_to_blender(ctx)
        self.critical.assert_called_once()


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


class TestBlenderExportDativeBonds:
    @pytest.fixture(autouse=True)
    def _qt(self, qapp, monkeypatch):
        self._monkeypatch = monkeypatch
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
        self._monkeypatch.setattr(_blender, "PLUGIN_CONTEXT", ctx)
        self._save_ret = (str(tmp_path / "scene.py"), "")
        _blender.export_to_blender(ctx)
        self.critical.assert_not_called()
        return (tmp_path / "scene.py").read_text(encoding="utf-8")

    @pytest.mark.parametrize("kind", ["DATIVE", "ZERO"])
    def test_bond_exported_as_single_cylinder(self, tmp_path, kind):
        script = self._export(tmp_path, _real_mol(kind))
        assert script.count("obj = create_bond(") == 1
        assert "color=(0.6, 0.6, 0.6)," in script

    def test_dative_bond_uses_host_split_colors(self, tmp_path):
        script = self._export(
            tmp_path,
            _real_mol("DATIVE"),
            color_map={"bond_0_start": (255, 0, 0), "bond_0_end": (0, 0, 255)},
        )
        assert script.count("obj = create_bond(") == 2
        assert "color=(1.0, 0.0, 0.0)," in script
        assert "color=(0.0, 0.0, 1.0)," in script
