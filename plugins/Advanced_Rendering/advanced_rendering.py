#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Advanced Rendering Plugin (Fixed & Optimized)
Provides detailed control over Scene and Atomic rendering (PBR, Shadows, etc.).
Fixed issues with window lifecycle, actor validation, and render pipeline conflicts.
"""

import os
import json
import logging
import math
import numpy as np
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
    QCheckBox, QSlider, QComboBox, QPushButton, 
    QGroupBox, QTabWidget, QFileDialog, QDialog,
    QMessageBox
)
from PyQt6.QtCore import Qt, QTimer, QCoreApplication
from PyQt6.QtGui import QColor, QCloseEvent

# PyVista / VTK Import Check
try:
    import pyvista as pv
    import vtk
except ImportError:
    pv = None
    vtk = None

PLUGIN_NAME = "Advanced Rendering"
PLUGIN_VERSION = "2026.02.03-fix"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Fine-grained control over Scene lighting, shadows, and PBR effects. (Stability Fixed)"

# --- HELPER CLASSES ---

class HideOnCloseDialog(QDialog):
    """
    A Dialog that hides instead of closes, keeping the internal widget alive.
    Essential for preserving plugin state when the user closes the window.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
    
    def closeEvent(self, event: QCloseEvent):
        # Do not destroy the widget, just hide it
        event.ignore()
        self.hide()

# --- PLUGIN ENTRY POINTS ---

def get_icon():
    return None

def load_plugin(main_window):
    """
    Shows the plugin dialog. Creates it only if necessary.
    """
    # 1. Retrieve the persistent dialog instance
    dialog = getattr(main_window, '_adv_graphics_dialog', None)
    
    # 2. If it doesn't exist or was accidentally destroyed, recreate it
    if not dialog or not isinstance(dialog, QDialog):
        # Ensure the viewer widget exists
        viewer = getattr(main_window, '_adv_rendering_viewer', None)
        if not viewer:
            # Should have been created in initialize, but fallback just in case
            viewer = AdvancedGraphicsWidget(main_window)
            main_window._adv_rendering_viewer = viewer
        
        # Create the container dialog
        dialog = HideOnCloseDialog(main_window)
        dialog.setWindowTitle("Advanced Graphics Settings")
        dialog.setMinimumWidth(400)
        
        dlg_layout = QVBoxLayout(dialog)
        dlg_layout.setContentsMargins(0, 0, 0, 0)
        
        # Reparent viewer to this dialog
        viewer.setParent(dialog)
        dlg_layout.addWidget(viewer)
        
        main_window._adv_graphics_dialog = dialog

    # 3. Show the dialog
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()

def initialize(context):
    """
    Plugin entry point. Initializes backend widget and hooks into menus.
    """
    mw = context.get_main_window()
    
    # 1. Create Background Widget (State Holder)
    if not hasattr(mw, '_adv_rendering_viewer'):
        viewer = AdvancedGraphicsWidget(mw)
        mw._adv_rendering_viewer = viewer
    else:
        viewer = mw._adv_rendering_viewer

    # 2. Register Menu Action
    context.add_menu_action("Settings/Advanced Graphics Settings", lambda: load_plugin(mw))

    # 3. Register PBR Styles & Draw Overrides
    def make_style_drawer(base_style_key, force_pbr=False):
        def drawer(mw_obj, mol):
            # A. Draw standard version using built-in method
            if hasattr(mw_obj, 'main_window_view_3d'):
                mw_obj.main_window_view_3d.draw_standard_3d_style(mol, style_override=base_style_key)
            
            # B. Apply PBR/Effects
            v = getattr(mw_obj, '_adv_rendering_viewer', None)
            if v:
                if force_pbr:
                    # Temporarily force PBR ON for this specific style render
                    v.apply_pbr_forced()
                else:
                    # Use current checkbox state
                    v.update_atoms_pbr()
                
                # Sync UI to match the new style
                v.sync_style_ui(getattr(mw_obj, 'current_3d_style', ''))
                    
        return drawer

    # Register PBR variants
    context.register_3d_style("Ball & Stick (PBR)", make_style_drawer('ball_and_stick', force_pbr=True))
    context.register_3d_style("CPK (PBR)", make_style_drawer('cpk', force_pbr=True))
    context.register_3d_style("Wireframe (PBR)", make_style_drawer('wireframe', force_pbr=True))
    context.register_3d_style("Stick (PBR)", make_style_drawer('stick', force_pbr=True))


# --- MAIN WIDGET ---

class AdvancedGraphicsWidget(QWidget):
    def __init__(self, parent_window):
        super().__init__(parent_window)
        self.mw = parent_window
        self.plotter = getattr(self.mw, 'plotter', None)
        
        # State Variables
        self.env_texture_path = ""
        self.use_shadows = False
        self.light_intensity = 1.0
        self.light_azimuth = 0
        self.light_elevation = 45
        
        self.use_ssao = False
        self.use_depth_peeling = False
        self.use_aa = False
        self.use_edl = False
        self.edl_strength = 0.2
        
        # Atom PBR
        self.atom_metallic = 0.0
        self.atom_roughness = 0.5
        self.use_atom_pbr = False
        self.use_atom_silhouette = False
        self.atom_actor_sil = None
        
        # Presets
        self.presets = {}
        self.default_preset_names = {"Default"}
        self._init_default_presets()
        
        self.init_ui()
        self.update_preset_combo()
        
        # Delay load to ensure UI is ready
        QTimer.singleShot(100, self.load_settings)
        
        # Save on exit
        QCoreApplication.instance().aboutToQuit.connect(self.save_settings)
        
        # Polling for external style changes
        self._last_polled_style = ""
        self._sync_timer = QTimer(self)
        self._sync_timer.timeout.connect(self._poll_style_change)
        self._sync_timer.start(500)

    def _init_default_presets(self):
        self.presets["Default"] = {
            "use_shadows": False,
            "light_intensity": 1.0,
            "light_azimuth": 0,
            "light_elevation": 45,
            "use_ssao": False,
            "use_depth_peeling": False,
            "use_aa": False,
            "use_edl": False,
            "edl_strength": 0.2,
            "env_texture_path": "",
            "use_atom_pbr": False,
            "atom_metallic": 0.0,
            "atom_roughness": 0.5
        }

    def get_settings_path(self):
        # Prefer user home dir or app config dir, but falling back to plugin dir for portability
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "advanced_graphics_settings.json")

    def init_ui(self):
        layout = QVBoxLayout()
        
        # --- PRESETS ---
        group_preset = QGroupBox("Presets")
        preset_layout = QHBoxLayout()
        preset_layout.addWidget(QLabel("Name:"))
        
        self.combo_presets = QComboBox()
        self.combo_presets.setEditable(False)
        self.combo_presets.setPlaceholderText("Select preset...")
        self.combo_presets.setMinimumWidth(150)
        self.combo_presets.activated.connect(self.on_preset_activated)
        preset_layout.addWidget(self.combo_presets)
        
        btn_save = QPushButton("Save")
        btn_save.clicked.connect(self.save_preset)
        preset_layout.addWidget(btn_save)
        
        btn_del = QPushButton("Del")
        btn_del.clicked.connect(self.delete_preset)
        preset_layout.addWidget(btn_del)
        
        group_preset.setLayout(preset_layout)
        layout.addWidget(group_preset)
        
        # --- TABS ---
        tabs = QTabWidget()
        
        # TAB 1: SCENE
        tab_scene = QWidget()
        scene_layout = QVBoxLayout()
        
        # Texture
        tex_layout = QHBoxLayout()
        tex_layout.addWidget(QLabel("Env Texture:"))
        btn_browse_tex = QPushButton("Browse...")
        btn_browse_tex.clicked.connect(self.on_load_env_texture)
        tex_layout.addWidget(btn_browse_tex)
        btn_clear_tex = QPushButton("X")
        btn_clear_tex.setFixedWidth(30)
        btn_clear_tex.clicked.connect(self.on_remove_env_texture)
        tex_layout.addWidget(btn_clear_tex)
        scene_layout.addLayout(tex_layout)
        
        # Shadows
        self.check_shadows = QCheckBox("Enable Shadows (Requires PBR)")
        self.check_shadows.toggled.connect(self.on_shadows_toggled)
        scene_layout.addWidget(self.check_shadows)
        
        # --- Light Control ---
        light_group = QGroupBox("Lighting")
        light_group_layout = QVBoxLayout()

        # Light Intensity
        light_layout = QHBoxLayout()
        light_layout.addWidget(QLabel("Intensity:"))
        self.slider_light = QSlider(Qt.Orientation.Horizontal)
        self.slider_light.setRange(0, 500) # 0.0 - 5.0
        self.slider_light.setValue(100)
        self.slider_light.valueChanged.connect(self.on_light_changed)
        light_layout.addWidget(self.slider_light)
        light_group_layout.addLayout(light_layout)

        # Light Azimuth
        azi_layout = QHBoxLayout()
        azi_layout.addWidget(QLabel("Azimuth:"))
        self.slider_azi = QSlider(Qt.Orientation.Horizontal)
        self.slider_azi.setRange(-180, 180)
        self.slider_azi.setValue(0)
        self.slider_azi.valueChanged.connect(self.on_light_pos_changed)
        azi_layout.addWidget(self.slider_azi)
        light_group_layout.addLayout(azi_layout)

        # Light Elevation
        ele_layout = QHBoxLayout()
        ele_layout.addWidget(QLabel("Elevation:"))
        self.slider_ele = QSlider(Qt.Orientation.Horizontal)
        self.slider_ele.setRange(-90, 90)
        self.slider_ele.setValue(45)
        self.slider_ele.valueChanged.connect(self.on_light_pos_changed)
        ele_layout.addWidget(self.slider_ele)
        light_group_layout.addLayout(ele_layout)

        scene_layout.addWidget(light_group)
        
        # Effects
        self.check_ssao = QCheckBox("SSAO (Ambient Occlusion)")
        self.check_ssao.toggled.connect(self.on_ssao_toggled)
        scene_layout.addWidget(self.check_ssao)
        
        self.check_depth = QCheckBox("Depth Peeling (Transparency Fix)")
        self.check_depth.toggled.connect(self.on_depth_peeling_toggled)
        scene_layout.addWidget(self.check_depth)
        
        self.check_aa = QCheckBox("Anti-Aliasing (MSAA)")
        self.check_aa.toggled.connect(self.on_aa_toggled)
        scene_layout.addWidget(self.check_aa)
        
        self.check_edl = QCheckBox("Eye Dome Lighting (Depth Shading)")
        self.check_edl.toggled.connect(self.on_edl_toggled)
        scene_layout.addWidget(self.check_edl)
        
        edl_str_layout = QHBoxLayout()
        edl_str_layout.addWidget(QLabel("EDL Strength:"))
        self.slider_edl = QSlider(Qt.Orientation.Horizontal)
        self.slider_edl.setRange(1, 100)
        self.slider_edl.setValue(20)
        self.slider_edl.valueChanged.connect(self.on_edl_strength_changed)
        self.slider_edl.setEnabled(False)
        edl_str_layout.addWidget(self.slider_edl)
        scene_layout.addLayout(edl_str_layout)
        
        scene_layout.addStretch()
        tab_scene.setLayout(scene_layout)
        tabs.addTab(tab_scene, "Scene")
        
        # TAB 2: ATOMS
        tab_atoms = QWidget()
        atoms_layout = QVBoxLayout()
        
        self.check_atom_pbr = QCheckBox("Enable Atom PBR")
        self.check_atom_pbr.toggled.connect(self.on_atom_pbr_toggled)
        atoms_layout.addWidget(self.check_atom_pbr)
        
        # Apply/Clear
        button_layout = QHBoxLayout()
        btn_apply = QPushButton("Force Apply")
        btn_apply.setToolTip("Re-apply settings to current scene actors")
        btn_apply.clicked.connect(self.update_atoms_pbr)
        button_layout.addWidget(btn_apply)
        
        btn_clear = QPushButton("Reset Actors")
        btn_clear.clicked.connect(self.clear_atom_settings)
        button_layout.addWidget(btn_clear)
        atoms_layout.addLayout(button_layout)
        
        # Metallic
        met_layout = QHBoxLayout()
        met_layout.addWidget(QLabel("Metallic:"))
        self.slider_atom_metallic = QSlider(Qt.Orientation.Horizontal)
        self.slider_atom_metallic.setRange(0, 100)
        self.slider_atom_metallic.setValue(0)
        self.slider_atom_metallic.valueChanged.connect(self.on_atom_metallic_changed)
        self.slider_atom_metallic.setEnabled(False)
        met_layout.addWidget(self.slider_atom_metallic)
        atoms_layout.addLayout(met_layout)
        
        # Roughness
        rough_layout = QHBoxLayout()
        rough_layout.addWidget(QLabel("Roughness:"))
        self.slider_atom_roughness = QSlider(Qt.Orientation.Horizontal)
        self.slider_atom_roughness.setRange(0, 100)
        self.slider_atom_roughness.setValue(50)
        self.slider_atom_roughness.valueChanged.connect(self.on_atom_roughness_changed)
        self.slider_atom_roughness.setEnabled(False)
        rough_layout.addWidget(self.slider_atom_roughness)
        atoms_layout.addLayout(rough_layout)
        
        # Silhouette
        self.check_atom_sil = QCheckBox("Silhouette (Outline)")
        self.check_atom_sil.toggled.connect(self.on_atom_silhouette_toggled)
        atoms_layout.addWidget(self.check_atom_sil)
        
        atoms_layout.addStretch()
        tab_atoms.setLayout(atoms_layout)
        tabs.addTab(tab_atoms, "Atoms")
        
        layout.addWidget(tabs)
        self.setLayout(layout)

    # --- CORE LOGIC ---

    def _poll_style_change(self):
        """Monitors current_3d_style and updates UI if changed."""
        if not self.mw: return
        curr = getattr(self.mw, 'current_3d_style', '')
        if curr != self._last_polled_style:
            self._last_polled_style = curr
            self.sync_style_ui(curr)

    def apply_pbr_forced(self):
        """Used by PBR-specific styles to force PBR on without changing UI check state permanently."""
        old_state = self.use_atom_pbr
        self.use_atom_pbr = True
        self.update_atoms_pbr()
        self.use_atom_pbr = old_state # Restore internal state

    def update_atoms_pbr(self):
        """Applies PBR settings to relevant actors."""
        if not self.plotter or not hasattr(self.plotter, 'renderer'): return

        try:
            curr_style = getattr(self.mw, 'current_3d_style', '')
            is_pbr_style = "(PBR)" in curr_style
            
            # If manually enabled via checkbox OR it's a PBR style, we use PBR
            should_use_pbr = self.use_atom_pbr or is_pbr_style
            
            # Sync UI Checkbox if it's a PBR style
            if is_pbr_style and hasattr(self, 'check_atom_pbr'):
                self.check_atom_pbr.blockSignals(True)
                self.check_atom_pbr.setChecked(True)
                self.check_atom_pbr.blockSignals(False)

            # --- FIX: Use 'actors' property (dict) or GetActors() for PyVista compatibility ---
            if hasattr(self.plotter.renderer, 'actors'):
                # PyVista standard: actors is a dict {name: actor}
                actors_list = list(self.plotter.renderer.actors.values())
            elif hasattr(self.plotter.renderer, 'GetActors'):
                # VTK fallback
                actors_list = self.plotter.renderer.GetActors()
            else:
                logging.warning("Advanced Rendering: Could not retrieve actors from renderer.")
                return
            
            for actor in actors_list:
                if not actor: continue
                # Validate it is a VTK Actor (Geometry)
                if vtk and not actor.IsA("vtkActor"):
                    continue

                prop = actor.GetProperty()
                if not prop: continue

                if should_use_pbr:
                    if hasattr(prop, 'SetInterpolationToPBR'):
                        prop.SetInterpolationToPBR()
                        prop.SetMetallic(self.atom_metallic)
                        prop.SetRoughness(self.atom_roughness)
                else:
                    if hasattr(prop, 'SetInterpolationToPhong'):
                        prop.SetInterpolationToPhong()
                        # Phong defaults
                        prop.SetMetallic(0.0)
                        prop.SetRoughness(0.5) # Not strictly used in Phong but safe reset
            
            self.plotter.render()
        except Exception as e:
            logging.error(f"PBR Update Error: {e}")

    # --- EVENT HANDLERS ---

    def on_atom_pbr_toggled(self, checked):
        self.use_atom_pbr = checked
        self.slider_atom_metallic.setEnabled(checked)
        self.slider_atom_roughness.setEnabled(checked)
        self.update_atoms_pbr()
        
    def on_atom_metallic_changed(self, val):
        self.atom_metallic = val / 100.0
        if self.use_atom_pbr: self.update_atoms_pbr()
        
    def on_atom_roughness_changed(self, val):
        self.atom_roughness = val / 100.0
        if self.use_atom_pbr: self.update_atoms_pbr()
        
    def on_atom_silhouette_toggled(self, checked):
        self.use_atom_silhouette = checked
        if not self.plotter: return

        try:
             # Remove old silhouette if exists
             if self.atom_actor_sil:
                 self.plotter.remove_actor(self.atom_actor_sil)
                 self.atom_actor_sil = None
             
             if checked:
                 # Try to find a valid target actor
                 target_actor = None
                 if hasattr(self.mw, 'atom_actor') and self.mw.atom_actor:
                     target_actor = self.mw.atom_actor
                 elif hasattr(self.mw, 'main_window_view_3d'):
                      if hasattr(self.mw.main_window_view_3d, 'atom_actor'):
                          target_actor = self.mw.main_window_view_3d.atom_actor
                 
                 if target_actor and hasattr(self.plotter, 'add_silhouette'):
                      self.atom_actor_sil = self.plotter.add_silhouette(
                          target_actor, color='black', line_width=2.0
                      )
             self.plotter.render()
        except Exception as e:
            logging.error(f"Silhouette Toggle Error: {e}")

    def clear_atom_settings(self):
        """Reset atom settings to defaults without clearing the whole scene."""
        self.check_atom_pbr.blockSignals(True)
        self.check_atom_pbr.setChecked(False)
        self.check_atom_pbr.blockSignals(False)
        self.use_atom_pbr = False
        
        self.slider_atom_metallic.setValue(0)
        self.slider_atom_roughness.setValue(50)
        
        # Force Phong on everything
        try:
            # --- FIX: Use correct actors property ---
            if hasattr(self.plotter.renderer, 'actors'):
                actors_list = list(self.plotter.renderer.actors.values())
            elif hasattr(self.plotter.renderer, 'GetActors'):
                actors_list = self.plotter.renderer.GetActors()
            else:
                return

            for actor in actors_list:
                 if vtk and not actor.IsA("vtkActor"): continue
                 prop = actor.GetProperty()
                 if hasattr(prop, 'SetInterpolationToPhong'):
                     prop.SetInterpolationToPhong()
            self.plotter.render()
        except Exception as e:
            logging.error(f"Clear Settings Error: {e}")

    # SCENE HANDLERS
    def on_load_env_texture(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Open Texture", "", "Images (*.png *.jpg *.jpeg *.hdr)")
        if fname:
            self.env_texture_path = fname
            self.apply_texture()
            
    def on_remove_env_texture(self):
        self.env_texture_path = ""
        self.apply_texture()
        
    def apply_texture(self):
        if not self.plotter: return
        try:
            if self.env_texture_path and os.path.exists(self.env_texture_path):
                # Ensure read_texture is used for robustness
                texture = pv.read_texture(self.env_texture_path)
                self.plotter.set_environment_texture(texture)
            else:
                self.plotter.remove_environment_texture()
        except Exception as e:
            logging.error(f"Failed to load environment texture: {e}")
            QMessageBox.warning(self, "Texture Load Error", f"Could not load texture:\n{e}")
        
    def on_shadows_toggled(self, checked):
        if not self.plotter: return
        
        curr_style = getattr(self.mw, 'current_3d_style', '')
        if checked and "(PBR)" not in curr_style and not self.use_atom_pbr:
             # Warn user or auto-enable PBR? For now just auto-disable shadows
             self.check_shadows.blockSignals(True)
             self.check_shadows.setChecked(False)
             self.check_shadows.blockSignals(False)
             return

        self.use_shadows = checked
        if checked:
            self._disable_conflicting_effects(exclude="shadows")

        try:
            if checked:
                self.plotter.enable_shadows()
            else:
                self.plotter.disable_shadows()
            self.plotter.render()
        except Exception as e:
            logging.warning(f"Shadow error: {e}")

    def _disable_conflicting_effects(self, exclude=""):
        """
        Mutually exclusive effects logic.
        Uses blockSignals to prevent recursion.
        """
        self.blockSignals(True) 
        
        if exclude == "depth":
            # Depth Peeling is incompatible with shadows/EDL/SSAO in some VTK versions
            if self.use_edl:
                self.check_edl.setChecked(False)
                self.use_edl = False
                try: 
                    self.plotter.disable_eye_dome_lighting()
                    self.plotter.render() # Added to clean up EDL resources
                except: pass
            
            if self.use_shadows:
                self.check_shadows.setChecked(False)
                self.use_shadows = False
                try: self.plotter.disable_shadows()
                except: pass

            if self.use_ssao:
                self.check_ssao.setChecked(False)
                self.use_ssao = False
                try: self.plotter.disable_ssao()
                except: pass

        elif exclude in ["edl", "shadows", "ssao"]:
            if self.use_depth_peeling:
                self.check_depth.setChecked(False)
                self.use_depth_peeling = False
                try: self.plotter.disable_depth_peeling()
                except: pass
        
        self.blockSignals(False)
    
    def on_light_pos_changed(self):
        self.light_azimuth = self.slider_azi.value()
        self.light_elevation = self.slider_ele.value()
        self.update_lights()

    def update_lights(self):
        if not self.plotter or not hasattr(self.plotter, 'renderer'): return
        
        try:
            # Calculate new position (Spherical -> Cartesian)
            # Assumption: Center is (0,0,0) or Camera Focal Point
            # Assuming camera focal point is better if we want to orbit the subject
            center = (0,0,0)
            if hasattr(self.plotter, 'camera'):
                center = self.plotter.camera.focal_point

            r = 100.0 # Distance
            
            rad_azi = math.radians(self.light_azimuth)
            rad_ele = math.radians(self.light_elevation)
            
            # Y-up world? VTK default is often Y-up or Z-up depending on view.
            # Standard Math: Z-up
            # x = r * cos(ele) * cos(azi)
            # y = r * cos(ele) * sin(azi)
            # z = r * sin(ele)
            
            # Let's assume generic orbit
            x = center[0] + r * math.cos(rad_ele) * math.sin(rad_azi)
            y = center[1] + r * math.sin(rad_ele)
            z = center[2] + r * math.cos(rad_ele) * math.cos(rad_azi)
            
            lights = self.plotter.renderer.GetLights()
            # Modify the first light (Headlight/Keylight)
            found = False
            for light in lights:
                if light.GetSwitch(): # Only modify active lights
                    light.SetPosition(x, y, z)
                    light.SetFocalPoint(center)
                    # If it was a headlight, it might be locked to camera. 
                    # We might need to ensure it's not.
                    found = True
                    # Only do one for now to avoid chaotic lighting
                    break 
            
            if not found:
                 # Create a light if none?
                 pass
            
            self.plotter.render()
        except Exception as e:
            pass

    def on_light_changed(self, val):
        self.light_intensity = val / 100.0
        if not self.plotter: return
        try:
            if hasattr(self.plotter, 'renderer'):
                lights = self.plotter.renderer.GetLights()
                for light in lights:
                    light.SetIntensity(self.light_intensity)
            self.plotter.render()
        except Exception: pass

    def on_ssao_toggled(self, checked):
        self.use_ssao = checked
        if checked:
            self._disable_conflicting_effects(exclude="ssao")
        try:
            if checked: self.plotter.enable_ssao()
            else: self.plotter.disable_ssao()
            self.plotter.render()
        except Exception: pass

    def on_depth_peeling_toggled(self, checked):
        self.use_depth_peeling = checked
        if checked:
            self._disable_conflicting_effects(exclude="depth")
        try:
            if checked: self.plotter.enable_depth_peeling()
            else: self.plotter.disable_depth_peeling()
            self.plotter.render()
        except Exception: pass

    def on_aa_toggled(self, checked):
        self.use_aa = checked
        try:
            if checked: self.plotter.enable_anti_aliasing()
            else: self.plotter.disable_anti_aliasing()
            self.plotter.render()
        except: pass

    def on_edl_toggled(self, checked):
        self.use_edl = checked
        if checked:
            self._disable_conflicting_effects(exclude="edl")
        
        self.slider_edl.setEnabled(checked)
        self.apply_edl()
        
    def on_edl_strength_changed(self, val):
        self.edl_strength = val / 100.0
        if self.use_edl: self.apply_edl()
        
    def apply_edl(self):
        if not self.plotter: return
        try:
            if self.use_edl:
                self.plotter.enable_eye_dome_lighting()
                if hasattr(self.plotter.renderer, '_edl_pass'):
                    edl_pass = self.plotter.renderer._edl_pass
                    if edl_pass and hasattr(edl_pass, 'SetEDLStrength'):
                        edl_pass.SetEDLStrength(self.edl_strength)
            else:
                self.plotter.disable_eye_dome_lighting()
            self.plotter.render()
        except Exception: pass

    # --- PRESETS & PERSISTENCE ---

    def update_preset_combo(self):
        self.combo_presets.blockSignals(True)
        curr = self.combo_presets.currentText()
        self.combo_presets.clear()
        self.combo_presets.addItems(sorted(self.presets.keys()))
        if curr in self.presets:
            self.combo_presets.setCurrentText(curr)
        self.combo_presets.blockSignals(False)

    def on_preset_activated(self, idx):
        name = self.combo_presets.currentText()
        if name in self.presets:
            self.apply_settings_dict(self.presets[name])
            
    def save_preset(self):
        from PyQt6.QtWidgets import QInputDialog
        name, ok = QInputDialog.getText(self, "Save Preset", "Enter preset name:", text=self.combo_presets.currentText())
        
        if not ok or not name.strip(): return
        name = name.strip()
        
        if name in self.default_preset_names:
            QMessageBox.warning(self, "Locked", f"'{name}' is a built-in preset and cannot be overwritten.")
            return
        
        self.presets[name] = self.gather_settings_dict()
        self.update_preset_combo()
        self.combo_presets.setCurrentText(name)
        self.save_settings()
        
    def delete_preset(self):
        name = self.combo_presets.currentText().strip()
        if name in self.default_preset_names:
            QMessageBox.warning(self, "Locked", f"'{name}' is a built-in preset.")
            return
        if name in self.presets:
            del self.presets[name]
            self.update_preset_combo()
            self.save_settings()

    def gather_settings_dict(self):
        return {
            "use_shadows": self.use_shadows,
            "light_intensity": self.light_intensity,
            "light_azimuth": self.light_azimuth,
            "light_elevation": self.light_elevation,
            "use_ssao": self.use_ssao,
            "use_depth_peeling": self.use_depth_peeling,
            "use_aa": self.use_aa,
            "use_edl": self.use_edl,
            "edl_strength": self.edl_strength,
            "env_texture_path": self.env_texture_path,
            "use_atom_pbr": self.use_atom_pbr,
            "atom_metallic": self.atom_metallic,
            "atom_roughness": self.atom_roughness
        }

    def apply_settings_dict(self, data):
        # Update UI first
        if "use_shadows" in data: self.check_shadows.setChecked(data["use_shadows"])
        
        # Light Settings
        if "light_intensity" in data: 
            self.light_intensity = float(data["light_intensity"])
            self.slider_light.setValue(int(self.light_intensity * 100))
        if "light_azimuth" in data:
            self.light_azimuth = int(data["light_azimuth"])
            self.slider_azi.setValue(self.light_azimuth)
        if "light_elevation" in data:
            self.light_elevation = int(data["light_elevation"])
            self.slider_ele.setValue(self.light_elevation)

        if "use_ssao" in data: self.check_ssao.setChecked(data["use_ssao"])
        if "use_depth_peeling" in data: self.check_depth.setChecked(data["use_depth_peeling"])
        if "use_aa" in data: self.check_aa.setChecked(data["use_aa"])
        if "use_edl" in data: self.check_edl.setChecked(data["use_edl"])
        if "edl_strength" in data: 
            self.edl_strength = float(data["edl_strength"])
            self.slider_edl.setValue(int(self.edl_strength * 100))
        
        if "env_texture_path" in data:
            self.env_texture_path = data["env_texture_path"]
            self.apply_texture()
            
        if "use_atom_pbr" in data: self.check_atom_pbr.setChecked(data["use_atom_pbr"])
        if "atom_metallic" in data:
            self.atom_metallic = float(data["atom_metallic"])
            self.slider_atom_metallic.setValue(int(self.atom_metallic * 100))
        if "atom_roughness" in data:
            self.atom_roughness = data.get("atom_roughness", 0.5)
            self.slider_atom_roughness.setValue(int(self.atom_roughness * 100))
        
        self.update_atoms_pbr() 
        self.update_lights() # Apply new light positions

    def sync_style_ui(self, active_style_name):
        """Syncs widget state with active style (simplified)."""
        is_active_pbr = "(PBR)" in active_style_name
        
        # Sync Widget State
        if hasattr(self, 'check_atom_pbr'):
            self.check_atom_pbr.blockSignals(True)
            self.check_atom_pbr.setChecked(is_active_pbr)
            if not is_active_pbr:
                self.use_atom_pbr = False
            self.check_atom_pbr.blockSignals(False)

    # --- FILE IO ---
    def save_settings(self):
        try:
            data = self.gather_settings_dict()
            data["presets"] = self.presets
            
            # Simple sanitization for JSON
            def sanitize(obj):
                if isinstance(obj, np.generic): return obj.item()
                if isinstance(obj, np.ndarray): return obj.tolist()
                if isinstance(obj, (QColor,)): return [obj.red(), obj.green(), obj.blue()]
                if isinstance(obj, dict): return {k: sanitize(v) for k, v in obj.items()}
                if isinstance(obj, list): return [sanitize(v) for v in obj]
                return obj
            
            with open(self.get_settings_path(), 'w') as f:
                json.dump(sanitize(data), f, indent=4)
        except Exception as e:
            logging.error(f"Save adv gfx error: {e}")
            
    def load_settings(self):
        try:
            path = self.get_settings_path()
            if os.path.exists(path):
                with open(path, 'r') as f:
                    data = json.load(f)
                
                if "presets" in data:
                    user_presets = data["presets"]
                    for name, preset_data in user_presets.items():
                        if name not in self.default_preset_names:
                            self.presets[name] = preset_data
                    self.update_preset_combo()
                    
                self.apply_settings_dict(data)
        except Exception as e:
            logging.error(f"Load adv gfx error: {e}")