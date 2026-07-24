"""
Tests for the Advanced Rendering plugin (get_icon; initialize -> 4 x register_3d_style + add_menu_action).
"""

from __future__ import annotations

import json
import logging
import math
import os

import numpy as np
import pytest

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from conftest import extract_function, load_plugin, make_context, mock_optional_imports

PLUGINS_DIR = Path(__file__).resolve().parents[1] / "plugins"
ADV_RENDER_PATH = PLUGINS_DIR / "Advanced_Rendering" / "advanced_rendering.py"


class TestAdvancedRendering:
    def test_get_icon_returns_none(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            assert mod.get_icon() is None

    def test_initialize_registers_menu_action(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            ctx.add_menu_action.assert_called_once()

    def test_initialize_registers_four_3d_styles(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            assert ctx.register_3d_style.call_count == 4

    def test_initialize_menu_path_contains_advanced(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            path = ctx.add_menu_action.call_args[0][0]
            assert "Advanced" in path

    def test_initialize_style_names_contain_advanced_rendering(self):
        """All registered 3D style keys end with '(Advanced Rendering)'."""
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
            style_keys = [
                call[0][0] for call in ctx.register_3d_style.call_args_list
            ]
            assert all("Advanced Rendering" in k for k in style_keys)


class TestAdvancedRenderingStyleDrawers:
    def _drawers(self):
        with mock_optional_imports():
            mod = load_plugin(ADV_RENDER_PATH)
            ctx = make_context()
            mod.initialize(ctx)
        return {
            call.args[0]: call.args[1]
            for call in ctx.register_3d_style.call_args_list
        }

    def test_four_styles_registered(self):
        drawers = self._drawers()
        assert set(drawers) == {
            "Ball & Stick (Advanced Rendering)",
            "CPK (Advanced Rendering)",
            "Wireframe (Advanced Rendering)",
            "Stick (Advanced Rendering)",
        }

    def test_drawer_uses_normalized_style_key(self):
        drawers = self._drawers()
        viewer = MagicMock()
        mw_obj = SimpleNamespace(
            view_3d_manager=MagicMock(), _adv_rendering_viewer=viewer
        )
        mw_obj.view_3d_manager.current_3d_style = "Ball & Stick (Advanced Rendering)"

        drawers["Ball & Stick (Advanced Rendering)"](mw_obj, "MOL")

        mw_obj.view_3d_manager.draw_standard_3d_style.assert_called_once_with(
            "MOL", style_override="ball_and_stick"
        )
        # Drawer must NOT force PBR — it honours the user's saved atom-PBR
        # preference through update_atoms_pbr(); forcing it rendered the model
        # flat gray (PBR with no environment texture).
        viewer.apply_pbr_forced.assert_not_called()
        viewer.update_atoms_pbr.assert_called_once()
        viewer.update_lights.assert_called_once()
        viewer.sync_style_ui.assert_called_once_with(
            "Ball & Stick (Advanced Rendering)"
        )

    def test_cpk_drawer_key(self):
        drawers = self._drawers()
        mw_obj = SimpleNamespace(
            view_3d_manager=MagicMock(), _adv_rendering_viewer=MagicMock()
        )
        drawers["CPK (Advanced Rendering)"](mw_obj, "MOL")
        mw_obj.view_3d_manager.draw_standard_3d_style.assert_called_once_with(
            "MOL", style_override="cpk"
        )

    def test_drawer_without_viewer_does_not_crash(self):
        drawers = self._drawers()
        mw_obj = SimpleNamespace(view_3d_manager=MagicMock())  # no viewer attr
        drawers["Stick (Advanced Rendering)"](mw_obj, "MOL")
        mw_obj.view_3d_manager.draw_standard_3d_style.assert_called_once()

    def test_drawer_without_view_manager_does_not_crash(self):
        # Regression (2026.07.08): sync_style_ui read mw_obj.view_3d_manager
        # unguarded even though the draw call above is hasattr-guarded for it.
        drawers = self._drawers()
        viewer = MagicMock()
        mw_obj = SimpleNamespace(_adv_rendering_viewer=viewer)
        drawers["Wireframe (Advanced Rendering)"](mw_obj, "MOL")
        viewer.sync_style_ui.assert_called_once_with("")

    def test_drawer_tears_down_effects_before_redraw(self):
        # Regression (2026.07.25): plotter.clear() during the redraw left the
        # EDL pass attached, tripping "vtkEDLShading: ColorTexture should have
        # been deleted in ReleaseGraphicsResources()". The drawer must release
        # the pass-based effects BEFORE the host clears/redraws the plotter.
        drawers = self._drawers()
        order = []
        viewer = MagicMock()
        viewer.teardown_scene_effects.side_effect = lambda: order.append("teardown")
        vm = MagicMock()
        vm.current_3d_style = "CPK (Advanced Rendering)"
        vm.draw_standard_3d_style.side_effect = lambda *a, **k: order.append("draw")
        mw_obj = SimpleNamespace(view_3d_manager=vm, _adv_rendering_viewer=viewer)

        drawers["CPK (Advanced Rendering)"](mw_obj, "MOL")

        viewer.teardown_scene_effects.assert_called_once()
        assert order == ["teardown", "draw"]


def _teardown_scene_effects_fn():
    return extract_function(
        ADV_RENDER_PATH,
        "AdvancedGraphicsWidget",
        "teardown_scene_effects",
        {"logging": logging},
    )


class TestTeardownSceneEffects:
    def test_disables_all_passes_and_renders(self):
        fn = _teardown_scene_effects_fn()
        plotter = MagicMock()
        fn(SimpleNamespace(plotter=plotter))
        plotter.disable_eye_dome_lighting.assert_called_once()
        plotter.disable_shadows.assert_called_once()
        plotter.disable_ssao.assert_called_once()
        plotter.disable_depth_peeling.assert_called_once()
        # A render is required so VTK actually runs ReleaseGraphicsResources.
        plotter.render.assert_called_once()

    def test_noop_without_plotter(self):
        fn = _teardown_scene_effects_fn()
        fn(SimpleNamespace(plotter=None))  # must not raise

    def test_render_failure_is_logged_not_raised(self):
        fn = _teardown_scene_effects_fn()
        plotter = MagicMock()
        plotter.render.side_effect = RuntimeError("gl context gone")
        fn(SimpleNamespace(plotter=plotter))  # must not raise


def _sync_style_ui_fn():
    return extract_function(
        ADV_RENDER_PATH, "AdvancedGraphicsWidget", "sync_style_ui", {"logging": logging}
    )


def _sync_stub(use_atom_pbr):
    return SimpleNamespace(
        check_atom_pbr=MagicMock(),
        slider_atom_metallic=MagicMock(),
        slider_atom_roughness=MagicMock(),
        use_atom_pbr=use_atom_pbr,
        update_atoms_pbr=MagicMock(),
        _enforce_scene_state=MagicMock(),
    )


class TestSyncStyleUiPbrPreference:
    def test_advanced_style_does_not_force_pbr_on(self):
        # Regression (2026.07.25): switching to an Advanced style used to force
        # use_atom_pbr=True, rendering the model flat gray. It must respect the
        # saved preference (here: OFF).
        fn = _sync_style_ui_fn()
        stub = _sync_stub(use_atom_pbr=False)
        fn(stub, "Ball & Stick (Advanced Rendering)")
        assert stub.use_atom_pbr is False
        stub.check_atom_pbr.setChecked.assert_called_with(False)
        stub._enforce_scene_state.assert_called_once_with(enable=True)

    def test_advanced_style_reflects_saved_pbr_on(self):
        fn = _sync_style_ui_fn()
        stub = _sync_stub(use_atom_pbr=True)
        fn(stub, "CPK (Advanced Rendering)")
        assert stub.use_atom_pbr is True
        stub.check_atom_pbr.setChecked.assert_called_with(True)
        stub.slider_atom_metallic.setEnabled.assert_called_with(True)
        stub.slider_atom_roughness.setEnabled.assert_called_with(True)

    def test_standard_style_disables_pbr_but_preserves_pref(self):
        fn = _sync_style_ui_fn()
        stub = _sync_stub(use_atom_pbr=True)
        fn(stub, "CPK")  # not an Advanced Rendering style
        # Saved preference must survive a trip through a standard style.
        assert stub.use_atom_pbr is True
        stub.check_atom_pbr.setEnabled.assert_called_with(False)
        stub.check_atom_pbr.setChecked.assert_called_with(False)
        stub._enforce_scene_state.assert_called_once_with(enable=False)


def _enforce_scene_state_fn():
    return extract_function(
        ADV_RENDER_PATH,
        "AdvancedGraphicsWidget",
        "_enforce_scene_state",
        {"logging": logging},
    )


def _enforce_stub(**overrides):
    stub = SimpleNamespace(
        plotter=MagicMock(),
        use_shadows=False,
        use_ssao=False,
        use_edl=False,
        use_aa=False,
        use_depth_peeling=False,
        _clean_render_pipeline=MagicMock(),
        update_lights=MagicMock(),
    )
    for k, v in overrides.items():
        setattr(stub, k, v)
    return stub


class TestEnforceSceneState:
    def test_cleans_pipeline_before_re_enabling_edl(self):
        # Regression (2026.07.25): the redraw re-apply path enabled passes
        # without resetting the pipeline first (unlike the manual toggle
        # handlers), so effects silently failed to re-attach.
        fn = _enforce_scene_state_fn()
        stub = _enforce_stub(use_edl=True)
        fn(stub, enable=True)
        stub._clean_render_pipeline.assert_called_once()
        stub.plotter.enable_eye_dome_lighting.assert_called_once()

    def test_enabling_shadows_reapplies_lights(self):
        fn = _enforce_scene_state_fn()
        stub = _enforce_stub(use_shadows=True)
        fn(stub, enable=True)
        stub.plotter.enable_shadows.assert_called_once()
        stub.update_lights.assert_called_once()

    def test_disable_path_does_not_clean_pipeline(self):
        fn = _enforce_scene_state_fn()
        stub = _enforce_stub(use_edl=True)
        fn(stub, enable=False)
        stub._clean_render_pipeline.assert_not_called()
        stub.plotter.disable_eye_dome_lighting.assert_called_once()

    def test_no_effects_enabled_does_not_clean_pipeline(self):
        fn = _enforce_scene_state_fn()
        stub = _enforce_stub()  # all effects off
        fn(stub, enable=True)
        stub._clean_render_pipeline.assert_not_called()


# ---------------------------------------------------------------------------
# Pure-logic method extraction: settings dict round-trip, lighting math,
# effect-exclusivity logic, and settings file persistence.
# ---------------------------------------------------------------------------


def _gather_settings_dict_fn():
    return extract_function(
        ADV_RENDER_PATH, "AdvancedGraphicsWidget", "gather_settings_dict"
    )


def _apply_settings_dict_fn():
    globs = {"os": os, "logging": logging}
    return extract_function(
        ADV_RENDER_PATH, "AdvancedGraphicsWidget", "apply_settings_dict", globs
    )


def _make_widget_stub(**overrides):
    stub = SimpleNamespace(
        use_shadows=False,
        light_intensity=1.0,
        light_azimuth=0,
        light_elevation=45,
        use_ssao=False,
        use_depth_peeling=False,
        use_aa=False,
        use_edl=False,
        env_texture_path="",
        use_atom_pbr=False,
        atom_metallic=0.0,
        atom_roughness=0.5,
        check_shadows=MagicMock(),
        slider_light=MagicMock(),
        slider_azi=MagicMock(),
        slider_ele=MagicMock(),
        check_ssao=MagicMock(),
        check_depth=MagicMock(),
        check_aa=MagicMock(),
        check_edl=MagicMock(),
        check_atom_pbr=MagicMock(),
        slider_atom_metallic=MagicMock(),
        slider_atom_roughness=MagicMock(),
        apply_texture=MagicMock(),
        update_atoms_pbr=MagicMock(),
        update_lights=MagicMock(),
    )
    for k, v in overrides.items():
        setattr(stub, k, v)
    return stub


class TestGatherSettingsDict:
    def test_gathers_all_expected_keys(self):
        fn = _gather_settings_dict_fn()
        stub = _make_widget_stub(atom_metallic=0.25, light_azimuth=-45)
        data = fn(stub)
        assert data["atom_metallic"] == 0.25
        assert data["light_azimuth"] == -45
        assert set(data.keys()) == {
            "use_shadows",
            "light_intensity",
            "light_azimuth",
            "light_elevation",
            "use_ssao",
            "use_depth_peeling",
            "use_aa",
            "use_edl",
            "env_texture_path",
            "use_atom_pbr",
            "atom_metallic",
            "atom_roughness",
        }


class TestApplySettingsDict:
    def test_round_trips_gather_output(self):
        gather = _gather_settings_dict_fn()
        apply_ = _apply_settings_dict_fn()
        stub = _make_widget_stub(
            use_shadows=True,
            light_intensity=2.5,
            light_azimuth=30,
            light_elevation=-10,
            use_atom_pbr=True,
            atom_metallic=0.7,
            atom_roughness=0.3,
        )
        data = gather(stub)

        fresh = _make_widget_stub()
        apply_(fresh, data)

        assert fresh.light_intensity == pytest.approx(2.5)
        assert fresh.light_azimuth == 30
        assert fresh.light_elevation == -10
        assert fresh.atom_metallic == pytest.approx(0.7)
        assert fresh.atom_roughness == pytest.approx(0.3)
        fresh.check_shadows.setChecked.assert_called_with(True)
        fresh.check_atom_pbr.setChecked.assert_called_with(True)
        fresh.update_atoms_pbr.assert_called_once()
        fresh.update_lights.assert_called_once()

    def test_missing_keys_are_left_untouched(self):
        apply_ = _apply_settings_dict_fn()
        stub = _make_widget_stub(light_intensity=9.9)
        apply_(stub, {})
        # Nothing in `data`, so light_intensity should be unchanged.
        assert stub.light_intensity == 9.9
        stub.check_shadows.setChecked.assert_not_called()

    def test_env_texture_path_triggers_apply_texture(self):
        apply_ = _apply_settings_dict_fn()
        stub = _make_widget_stub()
        apply_(stub, {"env_texture_path": "C:/tex.png"})
        assert stub.env_texture_path == "C:/tex.png"
        stub.apply_texture.assert_called_once()


class _FakeLight:
    def __init__(self):
        self.position = None
        self.focal_point = None
        self.intensity = None
        self.switch = None
        self.positional = None
        self._type = "CameraLight"

    def GetLightType(self):
        return self._type

    def SetLightTypeToCameraLight(self):
        self._type = "CameraLight"

    def SetPosition(self, x, y, z):
        self.position = (x, y, z)

    def SetFocalPoint(self, x, y, z):
        self.focal_point = (x, y, z)

    def SetIntensity(self, v):
        self.intensity = v

    def SetSwitch(self, v):
        self.switch = v

    def SetPositional(self, v):
        self.positional = v


class _FakeLights:
    def __init__(self, lights):
        self._lights = lights

    def GetNumberOfItems(self):
        return len(self._lights)

    def __iter__(self):
        return iter(self._lights)


class _FakeRenderer:
    def __init__(self, lights):
        self._lights_collection = _FakeLights(lights)
        self.added = []

    def GetLights(self):
        return self._lights_collection

    def AddLight(self, light):
        self.added.append(light)
        self._lights_collection = _FakeLights(self._lights_collection._lights + [light])


class _FakeVtkModule:
    VTK_LIGHT_TYPE_CAMERA_LIGHT = "CameraLight"

    def vtkLight(self):
        return _FakeLight()


def _update_lights_fn():
    globs = {"math": math, "logging": logging, "vtk": _FakeVtkModule()}
    return extract_function(
        ADV_RENDER_PATH, "AdvancedGraphicsWidget", "update_lights", globs
    )


class TestUpdateLights:
    def test_noop_when_not_advanced_style(self):
        fn = _update_lights_fn()
        plotter = MagicMock()
        stub = SimpleNamespace(
            mw=SimpleNamespace(view_3d_manager=SimpleNamespace(current_3d_style="CPK")),
            safe_plotter=plotter,
            light_azimuth=0,
            light_elevation=45,
            light_intensity=1.0,
        )
        fn(stub)
        plotter.render.assert_not_called()

    def test_computes_spherical_to_cartesian_position(self):
        fn = _update_lights_fn()
        light = _FakeLight()
        renderer = _FakeRenderer([light])
        plotter = SimpleNamespace(renderer=renderer, render=MagicMock())
        stub = SimpleNamespace(
            mw=SimpleNamespace(
                view_3d_manager=SimpleNamespace(
                    current_3d_style="Ball & Stick (Advanced Rendering)"
                )
            ),
            safe_plotter=plotter,
            light_azimuth=0,
            light_elevation=45,
            light_intensity=1.0,
        )
        fn(stub)

        r = 100.0
        expected_y = r * math.sin(math.radians(45))
        expected_h = r * math.cos(math.radians(45))
        expected_x = expected_h * math.sin(0)
        expected_z = expected_h * math.cos(0)
        assert light.position == pytest.approx((expected_x, expected_y, expected_z))
        assert light.focal_point == (0, 0, 0)
        assert light.intensity == 1.0
        assert light.switch is True
        plotter.render.assert_called_once()

    def test_creates_light_when_none_exist(self):
        fn = _update_lights_fn()
        renderer = _FakeRenderer([])
        plotter = SimpleNamespace(renderer=renderer, render=MagicMock())
        stub = SimpleNamespace(
            mw=SimpleNamespace(
                view_3d_manager=SimpleNamespace(
                    current_3d_style="CPK (Advanced Rendering)"
                )
            ),
            safe_plotter=plotter,
            light_azimuth=90,
            light_elevation=0,
            light_intensity=2.0,
        )
        fn(stub)
        assert len(renderer.added) == 1
        light = renderer.added[0]
        r = 100.0
        assert light.position == pytest.approx((r, 0.0, 0.0), abs=1e-6)

    def test_second_light_is_switched_off(self):
        fn = _update_lights_fn()
        main = _FakeLight()
        extra = _FakeLight()
        extra.switch = True
        renderer = _FakeRenderer([main, extra])
        plotter = SimpleNamespace(renderer=renderer, render=MagicMock())
        stub = SimpleNamespace(
            mw=SimpleNamespace(
                view_3d_manager=SimpleNamespace(
                    current_3d_style="Stick (Advanced Rendering)"
                )
            ),
            safe_plotter=plotter,
            light_azimuth=0,
            light_elevation=45,
            light_intensity=1.0,
        )
        fn(stub)
        assert extra.switch is False


def _disable_conflicting_effects_fn():
    globs = {"logging": logging}
    return extract_function(
        ADV_RENDER_PATH,
        "AdvancedGraphicsWidget",
        "_disable_conflicting_effects",
        globs,
    )


class TestDisableConflictingEffects:
    def test_exclude_depth_disables_edl_shadows_ssao(self):
        fn = _disable_conflicting_effects_fn()
        plotter = MagicMock()
        stub = SimpleNamespace(
            blockSignals=MagicMock(),
            use_edl=True,
            check_edl=MagicMock(),
            use_shadows=True,
            check_shadows=MagicMock(),
            use_ssao=True,
            check_ssao=MagicMock(),
            use_depth_peeling=True,
            check_depth=MagicMock(),
            plotter=plotter,
        )
        fn(stub, exclude="depth")
        assert stub.use_edl is False
        assert stub.use_shadows is False
        assert stub.use_ssao is False
        # depth itself is not touched by its own exclusion
        assert stub.use_depth_peeling is True
        plotter.disable_eye_dome_lighting.assert_called_once()
        plotter.disable_shadows.assert_called_once()
        plotter.disable_ssao.assert_called_once()

    def test_exclude_shadows_disables_depth_peeling_only(self):
        fn = _disable_conflicting_effects_fn()
        plotter = MagicMock()
        stub = SimpleNamespace(
            blockSignals=MagicMock(),
            use_edl=True,
            check_edl=MagicMock(),
            use_shadows=True,
            check_shadows=MagicMock(),
            use_ssao=True,
            check_ssao=MagicMock(),
            use_depth_peeling=True,
            check_depth=MagicMock(),
            plotter=plotter,
        )
        fn(stub, exclude="shadows")
        # exclude="shadows" only handles the depth-peeling branch
        assert stub.use_depth_peeling is False
        assert stub.use_edl is True
        assert stub.use_shadows is True
        plotter.disable_depth_peeling.assert_called_once()


class _DummyQColor:
    """Stand-in for PyQt6.QtGui.QColor — never instantiated by these tests,
    only used as an isinstance() target inside save_settings' sanitize()."""


def _save_settings_fn():
    globs = {
        "os": os,
        "json": json,
        "logging": logging,
        "np": np,
        "QColor": _DummyQColor,
    }
    return extract_function(
        ADV_RENDER_PATH, "AdvancedGraphicsWidget", "save_settings", globs
    )


def _load_settings_fn():
    globs = {"os": os, "json": json, "logging": logging}
    return extract_function(
        ADV_RENDER_PATH, "AdvancedGraphicsWidget", "load_settings", globs
    )


class TestSaveLoadSettings:
    def test_save_then_load_round_trip(self, tmp_path):
        settings_path = str(tmp_path / "adv_settings.json")
        save_fn = _save_settings_fn()
        gather = _gather_settings_dict_fn()

        save_stub = SimpleNamespace(
            gather_settings_dict=lambda: gather(
                _make_widget_stub(atom_metallic=0.42, light_azimuth=17)
            ),
            presets={"Default": {}, "MyPreset": {"atom_metallic": 0.9}},
            get_settings_path=lambda: settings_path,
        )
        save_fn(save_stub)

        assert os.path.exists(settings_path)
        with open(settings_path, "r", encoding="utf-8") as f:
            on_disk = json.load(f)
        assert on_disk["atom_metallic"] == pytest.approx(0.42)
        assert on_disk["presets"]["MyPreset"]["atom_metallic"] == pytest.approx(0.9)

        load_fn = _load_settings_fn()
        applied = {}
        load_stub = SimpleNamespace(
            get_settings_path=lambda: settings_path,
            presets={"Default": {}},
            default_preset_names={"Default"},
            update_preset_combo=MagicMock(),
            apply_settings_dict=lambda data: applied.update(data),
        )
        load_fn(load_stub)

        assert applied["atom_metallic"] == pytest.approx(0.42)
        assert load_stub.presets["MyPreset"]["atom_metallic"] == pytest.approx(0.9)
        load_stub.update_preset_combo.assert_called_once()

    def test_load_missing_file_does_not_crash(self, tmp_path):
        load_fn = _load_settings_fn()
        stub = SimpleNamespace(
            get_settings_path=lambda: str(tmp_path / "does_not_exist.json"),
            presets={},
            default_preset_names={"Default"},
            update_preset_combo=MagicMock(),
            apply_settings_dict=MagicMock(),
        )
        load_fn(stub)  # must not raise
        stub.apply_settings_dict.assert_not_called()

    def test_load_malformed_json_does_not_crash(self, tmp_path):
        bad_path = tmp_path / "bad.json"
        bad_path.write_text("{not valid json", encoding="utf-8")
        load_fn = _load_settings_fn()
        stub = SimpleNamespace(
            get_settings_path=lambda: str(bad_path),
            presets={},
            default_preset_names={"Default"},
            update_preset_combo=MagicMock(),
            apply_settings_dict=MagicMock(),
        )
        load_fn(stub)  # must not raise
        stub.apply_settings_dict.assert_not_called()

    def test_save_failure_is_logged_not_raised(self, tmp_path):
        save_fn = _save_settings_fn()
        gather = _gather_settings_dict_fn()
        stub = SimpleNamespace(
            gather_settings_dict=lambda: gather(_make_widget_stub()),
            presets={},
            # A directory that doesn't exist -> open() raises, caught internally.
            get_settings_path=lambda: str(tmp_path / "nope" / "settings.json"),
        )
        save_fn(stub)  # must not raise


# ---------------------------------------------------------------------------
# Environment texture: (re)apply for PBR image-based lighting, caching, and
# the UI preview label. (2026.07.25)
# ---------------------------------------------------------------------------


def _apply_texture_fn(read_texture):
    fake_pv = SimpleNamespace(read_texture=read_texture)
    return extract_function(
        ADV_RENDER_PATH,
        "AdvancedGraphicsWidget",
        "apply_texture",
        {"os": os, "logging": logging, "pv": fake_pv},
    )


class TestApplyTexture:
    def test_sets_environment_texture_when_path_exists(self, tmp_path):
        p = tmp_path / "env.png"
        p.write_bytes(b"x")
        fn = _apply_texture_fn(lambda pth: f"TEX:{pth}")
        plotter = MagicMock()
        stub = SimpleNamespace(
            plotter=plotter,
            env_texture_path=str(p),
            _update_texture_label=MagicMock(),
        )
        assert fn(stub) is True
        plotter.set_environment_texture.assert_called_once_with(f"TEX:{p}")
        stub._update_texture_label.assert_called_once()

    def test_caches_texture_across_calls(self, tmp_path):
        p = tmp_path / "env.png"
        p.write_bytes(b"x")
        reads = []
        fn = _apply_texture_fn(lambda pth: reads.append(pth) or "TEX")
        stub = SimpleNamespace(
            plotter=MagicMock(),
            env_texture_path=str(p),
            _update_texture_label=MagicMock(),
        )
        fn(stub)
        fn(stub)
        assert len(reads) == 1  # file read once, then served from cache

    def test_removes_texture_when_path_empty(self):
        fn = _apply_texture_fn(lambda pth: "TEX")
        plotter = MagicMock()
        stub = SimpleNamespace(
            plotter=plotter, env_texture_path="", _update_texture_label=MagicMock()
        )
        assert fn(stub) is True
        plotter.remove_environment_texture.assert_called_once()

    def test_read_failure_returns_false(self, tmp_path):
        p = tmp_path / "bad.png"
        p.write_bytes(b"x")

        def boom(pth):
            raise RuntimeError("cannot decode")

        fn = _apply_texture_fn(boom)
        stub = SimpleNamespace(
            plotter=MagicMock(),
            env_texture_path=str(p),
            _update_texture_label=MagicMock(),
        )
        assert fn(stub) is False
        stub._update_texture_label.assert_called_once()

    def test_noop_without_plotter_still_updates_label(self):
        fn = _apply_texture_fn(lambda pth: "TEX")
        stub = SimpleNamespace(
            plotter=None, env_texture_path="x", _update_texture_label=MagicMock()
        )
        assert fn(stub) is True
        stub._update_texture_label.assert_called_once()


class _FakePixmapValid:
    def __init__(self, *a):
        self._path = a[0] if a else ""

    def isNull(self):
        return not self._path

    def scaledToHeight(self, h, mode):
        return f"scaled:{self._path}:{h}"


class _FakePixmapNull:
    def __init__(self, *a):
        pass

    def isNull(self):
        return True

    def scaledToHeight(self, h, mode):  # pragma: no cover - never reached
        return self


_FAKE_QT = SimpleNamespace(
    TransformationMode=SimpleNamespace(SmoothTransformation=0)
)


def _update_texture_label_fn(qpixmap):
    return extract_function(
        ADV_RENDER_PATH,
        "AdvancedGraphicsWidget",
        "_update_texture_label",
        {"os": os, "QPixmap": qpixmap, "Qt": _FAKE_QT},
    )


class TestUpdateTextureLabel:
    def test_shows_thumbnail_for_decodable_image(self, tmp_path):
        p = tmp_path / "env.png"
        p.write_bytes(b"x")
        fn = _update_texture_label_fn(_FakePixmapValid)
        label = MagicMock()
        fn(SimpleNamespace(env_texture_path=str(p), label_tex=label))
        label.setPixmap.assert_called_once_with(f"scaled:{p}:60")
        label.setToolTip.assert_called_once_with(str(p))

    def test_shows_filename_when_pixmap_cannot_decode(self, tmp_path):
        p = tmp_path / "scene.hdr"
        p.write_bytes(b"x")
        fn = _update_texture_label_fn(_FakePixmapNull)
        label = MagicMock()
        fn(SimpleNamespace(env_texture_path=str(p), label_tex=label))
        label.setText.assert_called_once_with("scene.hdr")

    def test_shows_placeholder_when_no_texture(self):
        fn = _update_texture_label_fn(_FakePixmapValid)
        label = MagicMock()
        fn(SimpleNamespace(env_texture_path="", label_tex=label))
        label.setText.assert_called_once_with("(no texture)")

    def test_noop_without_label(self):
        fn = _update_texture_label_fn(_FakePixmapValid)
        fn(SimpleNamespace(env_texture_path="x", label_tex=None))  # must not raise


# ---------------------------------------------------------------------------
# PBR toggle applies the texture; EDL is enable/disable only (no strength API).
# ---------------------------------------------------------------------------


class TestOnAtomPbrToggled:
    def _fn(self):
        return extract_function(
            ADV_RENDER_PATH,
            "AdvancedGraphicsWidget",
            "on_atom_pbr_toggled",
            {"logging": logging},
        )

    def _stub(self, use):
        return SimpleNamespace(
            use_atom_pbr=use,
            slider_atom_metallic=MagicMock(),
            slider_atom_roughness=MagicMock(),
            apply_texture=MagicMock(),
            update_atoms_pbr=MagicMock(),
            save_settings=MagicMock(),
        )

    def test_enabling_pbr_applies_texture(self):
        stub = self._stub(False)
        self._fn()(stub, True)
        assert stub.use_atom_pbr is True
        stub.apply_texture.assert_called_once()
        stub.update_atoms_pbr.assert_called_once()

    def test_disabling_pbr_does_not_apply_texture(self):
        stub = self._stub(True)
        self._fn()(stub, False)
        stub.apply_texture.assert_not_called()


class TestUpdateAtomsPbrTexture:
    def _fn(self):
        return extract_function(
            ADV_RENDER_PATH,
            "AdvancedGraphicsWidget",
            "update_atoms_pbr",
            {"logging": logging, "vtk": None},
        )

    def _stub(self, use_atom_pbr, env_path):
        renderer = SimpleNamespace(actors={})  # no actors -> skip the prop loop
        plotter = SimpleNamespace(renderer=renderer, render=MagicMock())
        return SimpleNamespace(
            safe_plotter=plotter,
            mw=SimpleNamespace(
                view_3d_manager=SimpleNamespace(
                    current_3d_style="CPK (Advanced Rendering)"
                )
            ),
            use_atom_pbr=use_atom_pbr,
            env_texture_path=env_path,
            check_atom_pbr=MagicMock(),
            atom_metallic=0.0,
            atom_roughness=0.5,
            apply_texture=MagicMock(),
        )

    def test_reapplies_texture_when_pbr_active(self):
        stub = self._stub(use_atom_pbr=True, env_path="/env.png")
        self._fn()(stub)
        stub.apply_texture.assert_called_once()

    def test_no_texture_reapply_when_pbr_off(self):
        stub = self._stub(use_atom_pbr=False, env_path="/env.png")
        self._fn()(stub)
        stub.apply_texture.assert_not_called()

    def test_no_texture_reapply_when_no_path(self):
        stub = self._stub(use_atom_pbr=True, env_path="")
        self._fn()(stub)
        stub.apply_texture.assert_not_called()


class TestEdlControls:
    # The EDL strength slider was removed — vtkEDLShading has no strength API,
    # so EDL is enable/disable only.
    def test_no_edl_strength_slider_or_handler(self):
        src = ADV_RENDER_PATH.read_text(encoding="utf-8")
        assert "slider_edl" not in src
        assert "on_edl_strength_changed" not in src
        assert "edl_strength" not in src

    def _apply_edl_fn(self):
        return extract_function(
            ADV_RENDER_PATH, "AdvancedGraphicsWidget", "apply_edl", {"logging": logging}
        )

    def test_apply_edl_enable(self):
        plotter = MagicMock()
        self._apply_edl_fn()(
            SimpleNamespace(plotter=plotter, use_edl=True, use_aa=False)
        )
        plotter.enable_eye_dome_lighting.assert_called_once()

    def test_apply_edl_disable(self):
        plotter = MagicMock()
        self._apply_edl_fn()(
            SimpleNamespace(plotter=plotter, use_edl=False, use_aa=False)
        )
        plotter.disable_eye_dome_lighting.assert_called_once()
