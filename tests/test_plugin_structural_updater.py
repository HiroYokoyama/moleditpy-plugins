"""
Tests for the Structural Updater plugin.

Visible: true in REGISTRY/plugins.json.
All tests run headlessly — PyQt6, rdkit, etc. are mocked.
"""
from __future__ import annotations

import json
from pathlib import Path

from conftest import load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
STRUCTURAL_PATH = PLUGINS_DIR / "Structural_Updater" / "structural_updater.py"

# Load module once at module scope
with mock_optional_imports():
    SU = load_plugin(STRUCTURAL_PATH)


# ---------------------------------------------------------------------------
# Structural Updater — initialize / finalize
# ---------------------------------------------------------------------------

class TestStructuralUpdaterInitialize:
    def test_initialize_does_not_raise(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)

    def test_initialize_sets_plugin_instance(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)
        assert SU._PLUGIN_INSTANCE is not None

    def test_initialize_adds_menu_action(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)
        ctx.add_menu_action.assert_called()
        path = ctx.add_menu_action.call_args[0][0]
        assert "Structural Updater" in path

    def test_finalize_does_not_raise(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)
            SU.finalize()


# ---------------------------------------------------------------------------
# Structural Updater — settings load / save
# ---------------------------------------------------------------------------

class TestStructuralUpdaterSettings:
    def test_default_enabled_is_true(self):
        with mock_optional_imports():
            ctx = make_context()
            SU.initialize(ctx)
        assert SU._PLUGIN_INSTANCE.enabled is True

    def test_save_creates_file(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        plugin.save_settings()
        assert (tmp_path / "su.json").exists()

    def test_save_writes_enabled_true(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        plugin.enabled = True
        plugin.save_settings()
        data = json.loads((tmp_path / "su.json").read_text())
        assert data == {"enabled": True}

    def test_save_writes_enabled_false(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        plugin.enabled = False
        plugin.save_settings()
        data = json.loads((tmp_path / "su.json").read_text())
        assert data == {"enabled": False}

    def test_load_reads_enabled_false(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        (tmp_path / "su.json").write_text(json.dumps({"enabled": False}))
        plugin.load_settings()
        assert plugin.enabled is False

    def test_load_reads_enabled_true(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        (tmp_path / "su.json").write_text(json.dumps({"enabled": True}))
        plugin.load_settings()
        assert plugin.enabled is True

    def test_load_missing_file_creates_default(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        settings_file = tmp_path / "su.json"
        plugin.settings_file = str(settings_file)
        plugin.load_settings()
        assert settings_file.exists()  # save_settings() called automatically

    def test_round_trip(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.settings_file = str(tmp_path / "su.json")
        plugin.enabled = False
        plugin.save_settings()
        plugin.enabled = True  # change in memory
        plugin.load_settings()  # reload from file
        assert plugin.enabled is False


# ---------------------------------------------------------------------------
# Structural Updater — check_state
# ---------------------------------------------------------------------------

class TestStructuralUpdaterCheckState:
    def test_check_state_disabled_exits_early(self):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.enabled = False
        plugin.check_state()  # must not raise

    def test_check_state_enabled_does_not_raise(self):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        plugin.enabled = True
        # MagicMock() > 0 raises TypeError — configure GetNumAtoms to return int
        plugin.context.current_molecule.GetNumAtoms.return_value = 0
        plugin.check_state()  # must not raise




# ---------------------------------------------------------------------------
# trigger dispatch, apply-mode state, finalize
# ---------------------------------------------------------------------------

from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as real_numpy


class FakePos:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class FakeAtom:
    def __init__(self, symbol, num=6, props=None, radicals=0, idx=0):
        self._symbol = symbol
        self._num = num
        self._props = props or {}
        self._radicals = radicals
        self._idx = idx

    def GetSymbol(self):
        return self._symbol

    def GetAtomicNum(self):
        return self._num

    def GetNumRadicalElectrons(self):
        return self._radicals

    def GetIdx(self):
        return self._idx

    def HasProp(self, key):
        return key in self._props

    def GetProp(self, key):
        return self._props[key]

    def SetProp(self, key, val):
        self._props[key] = val


class FakeConf:
    def __init__(self, coords):
        self._coords = coords

    def GetAtomPosition(self, i):
        return FakePos(*self._coords[i])

    def GetPositions(self):
        return real_numpy.array(self._coords, dtype=float)


class FakeMol:
    def __init__(self, atoms, coords, bonds=0):
        self._atoms = atoms
        self._coords = coords
        self._bonds = bonds
        for i, a in enumerate(self._atoms):
            a._idx = i

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetNumBonds(self):
        return self._bonds

    def GetAtoms(self):
        return list(self._atoms)

    def GetAtomWithIdx(self, i):
        return self._atoms[i]

    def GetConformer(self):
        return FakeConf(self._coords)


def _updater_instance(mod, enabled=True, apply_mode=False, temp_mode=None):
    inst = object.__new__(mod.StructuralUpdaterPlugin)
    inst.enabled = enabled
    inst.apply_mode_active = apply_mode
    inst.context = MagicMock()
    inst.mw = MagicMock()
    if temp_mode is None:
        # getattr(self.mw, "_temp_conv_mode", None) must yield None
        inst.mw._temp_conv_mode = None
    else:
        inst.mw._temp_conv_mode = temp_mode
    inst.apply_changes_to_3d = MagicMock()
    return inst


class TestUpdaterTriggerDispatch:
    def _mod(self):
        with mock_optional_imports():
            return load_plugin(STRUCTURAL_PATH)

    def test_apply_mode_uses_plugin_logic(self):
        mod = self._mod()
        orig = MagicMock()
        mod._ORIGINAL_METHODS["trigger_conversion"] = orig
        inst = _updater_instance(mod, enabled=True, apply_mode=True)
        inst.new_trigger_conversion()
        inst.apply_changes_to_3d.assert_called_once()
        orig.assert_not_called()

    def test_disabled_falls_through_to_original(self):
        mod = self._mod()
        orig = MagicMock()
        mod._ORIGINAL_METHODS["trigger_conversion"] = orig
        inst = _updater_instance(mod, enabled=False, apply_mode=True)
        inst.new_trigger_conversion()
        orig.assert_called_once()
        inst.apply_changes_to_3d.assert_not_called()

    def test_no_apply_mode_falls_through(self):
        mod = self._mod()
        orig = MagicMock()
        mod._ORIGINAL_METHODS["trigger_conversion"] = orig
        inst = _updater_instance(mod, enabled=True, apply_mode=False)
        inst.new_trigger_conversion()
        orig.assert_called_once()

    def test_temp_mode_override_bypasses_plugin(self):
        mod = self._mod()
        orig = MagicMock()
        mod._ORIGINAL_METHODS["trigger_conversion"] = orig
        inst = _updater_instance(
            mod, enabled=True, apply_mode=True, temp_mode="fast"
        )
        inst.new_trigger_conversion()
        orig.assert_called_once()
        inst.apply_changes_to_3d.assert_not_called()

    def test_force_full_conversion_resets_apply_mode(self):
        mod = self._mod()
        orig = MagicMock()
        mod._ORIGINAL_METHODS["trigger_conversion"] = orig
        inst = _updater_instance(mod, apply_mode=True)
        inst.force_full_conversion()
        assert inst.apply_mode_active is False
        orig.assert_called_once()


class TestUpdaterCalculationFinished:
    def _run(self, mol, enabled=True):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
            orig = MagicMock()
            mod._ORIGINAL_METHODS["on_calculation_finished"] = orig
            inst = _updater_instance(mod, enabled=enabled)
            inst.new_trigger_conversion = MagicMock()
            inst.context.current_molecule = mol
            inst.new_on_calculation_finished("RESULT")
            return inst, orig

    def test_original_called_with_result(self):
        _, orig = self._run(mol=None)
        orig.assert_called_once_with("RESULT")

    def test_molecule_present_enables_apply_mode(self):
        mol = FakeMol([FakeAtom("C")], [(0, 0, 0)])
        inst, _ = self._run(mol=mol)
        assert inst.apply_mode_active is True
        inst.mw.init_manager.convert_button.setText.assert_called_with(
            "Apply 2D Changes to 3D"
        )

    def test_no_molecule_disables_apply_mode(self):
        inst, _ = self._run(mol=None)
        assert inst.apply_mode_active is False
        inst.mw.init_manager.convert_button.setText.assert_called_with(
            "Convert 2D to 3D"
        )

    def test_button_reconnected_to_wrapper(self):
        inst, _ = self._run(mol=None)
        inst.mw.init_manager.convert_button.clicked.connect.assert_called_with(
            inst.new_trigger_conversion
        )


class TestUpdaterCheckState:
    def _inst(self, mod, mol, apply_mode, halt=False):
        inst = _updater_instance(mod, enabled=True, apply_mode=apply_mode)
        inst.context.current_molecule = mol
        inst.mw.init_manager.convert_button.text.return_value = (
            "Halt conversion" if halt else "Convert 2D to 3D"
        )
        return inst

    def test_3d_props_switch_to_apply_mode(self):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
            mol = FakeMol(
                [FakeAtom("C", props={"_original_atom_id": 1})], [(0, 0, 0)]
            )
            inst = self._inst(mod, mol, apply_mode=False)
            inst.check_state()
            assert inst.apply_mode_active is True

    def test_no_props_switch_back(self):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
            mol = FakeMol([FakeAtom("C")], [(0, 0, 0)])
            inst = self._inst(mod, mol, apply_mode=True)
            inst.check_state()
            assert inst.apply_mode_active is False

    def test_halt_button_freezes_state(self):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
            mol = FakeMol(
                [FakeAtom("C", props={"_original_atom_id": 1})], [(0, 0, 0)]
            )
            inst = self._inst(mod, mol, apply_mode=False, halt=True)
            inst.check_state()
            assert inst.apply_mode_active is False

    def test_disabled_is_noop(self):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
            inst = _updater_instance(mod, enabled=False, apply_mode=False)
            inst.check_state()
            assert inst.apply_mode_active is False


class TestUpdaterFinalize:
    """Regression: finalize() restored originals onto mw, but the patches
    live on mw.compute_manager — so the patch was never actually undone."""

    def test_finalize_restores_compute_manager_methods(self):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
            orig_trigger = MagicMock(name="orig_trigger")
            orig_finished = MagicMock(name="orig_finished")
            mw = MagicMock()
            mod._ORIGINAL_METHODS.clear()
            mod._ORIGINAL_METHODS["trigger_conversion"] = orig_trigger
            mod._ORIGINAL_METHODS["on_calculation_finished"] = orig_finished
            mod._PLUGIN_INSTANCE = SimpleNamespace(mw=mw)

            mod.finalize()

            assert mw.compute_manager.trigger_conversion is orig_trigger
            assert mw.compute_manager.on_calculation_finished is orig_finished

    def test_finalize_without_instance_is_noop(self):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
            mod._PLUGIN_INSTANCE = None
            mod.finalize()  # must not raise

    def test_version_bumped(self):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
            assert mod.PLUGIN_VERSION > "2026.06.27"


# ---------------------------------------------------------------------------
# load_settings / save_settings — exception (warning-log) branches
# ---------------------------------------------------------------------------


class TestUpdaterSettingsExceptionPaths:
    def test_load_settings_bad_json_is_silenced(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        settings_file = tmp_path / "bad.json"
        settings_file.write_text("{not valid json")
        plugin.settings_file = str(settings_file)
        plugin.enabled = True
        plugin.load_settings()  # must not raise; exception is caught & logged
        assert plugin.enabled is True  # unchanged since load failed

    def test_save_settings_bad_path_is_silenced(self, tmp_path):
        with mock_optional_imports():
            ctx = make_context()
            plugin = SU.StructuralUpdaterPlugin(ctx)
        # A directory path cannot be opened for writing -> IsADirectoryError
        plugin.settings_file = str(tmp_path)
        plugin.save_settings()  # must not raise


# ---------------------------------------------------------------------------
# patch_mainwindow / new_on_calculation_finished — disconnect() exception
# ---------------------------------------------------------------------------


class TestUpdaterDisconnectExceptionPaths:
    def test_patch_mainwindow_disconnect_exception_is_silenced(self):
        with mock_optional_imports():
            ctx = make_context()
            mw = ctx.get_main_window()
            mw.init_manager.convert_button.clicked.disconnect.side_effect = TypeError(
                "no connections"
            )
            plugin = SU.StructuralUpdaterPlugin(ctx)  # must not raise
        # Connect still happens after the silenced exception
        mw.init_manager.convert_button.clicked.connect.assert_called_with(
            plugin.new_trigger_conversion
        )

    def test_calculation_finished_disconnect_exception_is_silenced(self):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
            orig = MagicMock()
            mod._ORIGINAL_METHODS["on_calculation_finished"] = orig
            inst = _updater_instance(mod, enabled=False)
            inst.mw.init_manager.convert_button.clicked.disconnect.side_effect = Exception(
                "boom"
            )
            inst.new_on_calculation_finished("RESULT")  # must not raise
            inst.mw.init_manager.convert_button.clicked.connect.assert_called_with(
                inst.new_trigger_conversion
            )


# ---------------------------------------------------------------------------
# apply_changes_to_3d — full core-logic coverage
# ---------------------------------------------------------------------------


class RAtom:
    """Rich fake RDKit atom exposing what apply_changes_to_3d needs."""

    def __init__(self, num, oid=None, idx=0):
        self._num = num
        self._oid = oid
        self._idx = idx

    def GetIdx(self):
        return self._idx

    def GetAtomicNum(self):
        return self._num

    def HasProp(self, key):
        return key == "_original_atom_id" and self._oid is not None

    def GetIntProp(self, key):
        return self._oid

    def GetBonds(self):
        return []

    def GetNeighbors(self):
        return []


class RConf:
    def __init__(self, coords):
        self._coords = list(coords)
        self.set_calls = []

    def GetAtomPosition(self, i):
        return FakePos(*self._coords[i])

    def SetAtomPosition(self, idx, pos):
        self.set_calls.append((idx, pos))


class RMol:
    def __init__(self, atoms, coords, conf_raises=False, num_conformers=1):
        self._atoms = atoms
        self._coords = coords
        self._conf_raises = conf_raises
        self._num_conformers = num_conformers
        for i, a in enumerate(atoms):
            a._idx = i

    def GetAtoms(self):
        return list(self._atoms)

    def GetAtomWithIdx(self, i):
        return self._atoms[i]

    def GetNumAtoms(self):
        return len(self._atoms)

    def GetConformer(self):
        if self._conf_raises:
            raise RuntimeError("no conformer")
        return RConf(self._coords)

    def GetNumConformers(self):
        return self._num_conformers


class TestApplyChangesToThreeD:
    def _mod(self):
        with mock_optional_imports():
            mod = load_plugin(STRUCTURAL_PATH)
        mod.Chem.AddHs.side_effect = lambda mol: mol
        return mod

    def _make_mols(self, n=4, matched=True):
        new_atoms = [RAtom(6, oid=(i if matched else 100 + i)) for i in range(n)]
        old_atoms = [RAtom(6, oid=i) for i in range(n)]
        coords = [(float(i), 0.0, 0.0) for i in range(n)]
        return RMol(new_atoms, coords), RMol(old_atoms, coords)

    def _inst(self, mod, new_mol, old_mol, mw=None):
        inst = object.__new__(mod.StructuralUpdaterPlugin)
        inst.context = MagicMock()
        inst.context.current_molecule = old_mol
        if mw is None:
            mw = MagicMock()
            mw.state_manager.data.to_rdkit_mol.return_value = new_mol
        inst.mw = mw
        inst.force_full_conversion = MagicMock()
        return inst

    def test_invalid_2d_structure(self):
        mod = self._mod()
        inst = self._inst(mod, None, MagicMock())
        inst.apply_changes_to_3d()
        inst.context.show_status_message.assert_any_call(
            "Error: Invalid 2D structure."
        )

    def test_no_old_mol_forces_full_conversion(self):
        mod = self._mod()
        new_mol, _ = self._make_mols()
        inst = self._inst(mod, new_mol, None)
        inst.apply_changes_to_3d()
        inst.force_full_conversion.assert_called_once()

    def test_conformer_exception_forces_full_conversion(self):
        mod = self._mod()
        new_mol, old_mol = self._make_mols()
        old_mol._conf_raises = True
        inst = self._inst(mod, new_mol, old_mol)
        inst.apply_changes_to_3d()
        inst.force_full_conversion.assert_called_once()

    def test_too_few_matches_forces_full_conversion(self):
        mod = self._mod()
        new_mol, old_mol = self._make_mols(matched=False)
        inst = self._inst(mod, new_mol, old_mol)
        inst.apply_changes_to_3d()
        inst.force_full_conversion.assert_called_once()
        inst.context.show_status_message.assert_any_call(
            "Notice: Too few matching atoms (0). Performing full conversion."
        )

    def test_embed_success_first_attempt_mmff(self):
        mod = self._mod()
        new_mol, old_mol = self._make_mols()
        mod.AllChem.EmbedMolecule.return_value = 0
        mod.AllChem.MMFFHasAllMoleculeParams.return_value = True
        ff = MagicMock()
        mod.AllChem.MMFFGetMoleculeForceField.return_value = ff
        inst = self._inst(mod, new_mol, old_mol)
        inst.apply_changes_to_3d()
        inst.mw.compute_manager.on_calculation_finished.assert_called_once_with(
            new_mol
        )
        ff.AddFixedPoint.assert_called()
        ff.Minimize.assert_called_once_with(maxIts=200)
        inst.context.show_status_message.assert_any_call(
            "Applied 2D changes to 3D structure."
        )

    def test_embed_success_uff_path(self):
        mod = self._mod()
        new_mol, old_mol = self._make_mols()
        mod.AllChem.EmbedMolecule.return_value = 0
        mod.AllChem.MMFFHasAllMoleculeParams.return_value = False
        ff = MagicMock()
        mod.AllChem.UFFGetMoleculeForceField.return_value = ff
        inst = self._inst(mod, new_mol, old_mol)
        inst.apply_changes_to_3d()
        mod.AllChem.UFFGetMoleculeForceField.assert_called_once()
        inst.mw.compute_manager.on_calculation_finished.assert_called_once()

    def test_ff_none_skips_minimize(self):
        mod = self._mod()
        new_mol, old_mol = self._make_mols()
        mod.AllChem.EmbedMolecule.return_value = 0
        mod.AllChem.MMFFHasAllMoleculeParams.return_value = True
        mod.AllChem.MMFFGetMoleculeForceField.return_value = None
        inst = self._inst(mod, new_mol, old_mol)
        inst.apply_changes_to_3d()
        inst.mw.compute_manager.on_calculation_finished.assert_called_once()

    def test_minimize_exception_forces_full_conversion(self):
        mod = self._mod()
        new_mol, old_mol = self._make_mols()
        mod.AllChem.EmbedMolecule.return_value = 0
        mod.AllChem.MMFFHasAllMoleculeParams.return_value = True
        ff = MagicMock()
        ff.Minimize.side_effect = RuntimeError("explode")
        mod.AllChem.MMFFGetMoleculeForceField.return_value = ff
        inst = self._inst(mod, new_mol, old_mol)
        inst.apply_changes_to_3d()
        inst.force_full_conversion.assert_called_once()
        inst.context.show_status_message.assert_any_call(
            "Notice: Structure optimization unstable. Re-generating full 3D structure...",
            5000,
        )

    def test_fallback1_success(self):
        mod = self._mod()
        new_mol, old_mol = self._make_mols()
        mod.AllChem.EmbedMolecule.side_effect = [1, 0]
        mod.AllChem.MMFFHasAllMoleculeParams.return_value = True
        mod.AllChem.MMFFGetMoleculeForceField.return_value = MagicMock()
        inst = self._inst(mod, new_mol, old_mol)
        inst.apply_changes_to_3d()
        assert mod.AllChem.EmbedMolecule.call_count == 2
        inst.mw.compute_manager.on_calculation_finished.assert_called_once()

    def test_fallback2_success(self):
        mod = self._mod()
        new_mol, old_mol = self._make_mols()
        mod.AllChem.EmbedMolecule.side_effect = [1, 1, 0]
        mod.AllChem.MMFFHasAllMoleculeParams.return_value = True
        mod.AllChem.MMFFGetMoleculeForceField.return_value = MagicMock()
        inst = self._inst(mod, new_mol, old_mol)
        inst.apply_changes_to_3d()
        assert mod.AllChem.EmbedMolecule.call_count == 3
        inst.mw.compute_manager.on_calculation_finished.assert_called_once()

    def test_all_embeds_fail_forces_full_conversion(self):
        mod = self._mod()
        new_mol, old_mol = self._make_mols()
        mod.AllChem.EmbedMolecule.side_effect = [1, 1, 1]
        inst = self._inst(mod, new_mol, old_mol)
        inst.apply_changes_to_3d()
        inst.force_full_conversion.assert_called_once()

    def test_legacy_data_fallback_uses_mw_data(self):
        """Regression: elif branch used mw.state_manager.data instead of
        mw.data, so it crashed for hosts without a state_manager attr."""
        mod = self._mod()
        new_mol, old_mol = self._make_mols()
        mod.AllChem.EmbedMolecule.return_value = 0
        mod.AllChem.MMFFHasAllMoleculeParams.return_value = True
        mod.AllChem.MMFFGetMoleculeForceField.return_value = MagicMock()

        class LegacyMw:
            pass

        legacy_mw = LegacyMw()
        legacy_mw.data = MagicMock()
        legacy_mw.data.to_rdkit_mol.return_value = new_mol
        legacy_mw.compute_manager = MagicMock()

        inst = self._inst(mod, None, old_mol, mw=legacy_mw)
        inst.apply_changes_to_3d()
        legacy_mw.data.to_rdkit_mol.assert_called_once()
        legacy_mw.compute_manager.on_calculation_finished.assert_called_once_with(
            new_mol
        )
