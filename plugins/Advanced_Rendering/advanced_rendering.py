#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Advanced Rendering Plugin (Fixed & Optimized)
Provides detailed control over Scene and Atomic rendering (PBR, Shadows, etc.).
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
try:
    from PyQt6 import sip
except ImportError:
    try:
        import sip
    except ImportError:
        sip = None
from PyQt6.QtCore import Qt, QTimer, QCoreApplication
from PyQt6.QtGui import QColor, QCloseEvent, QAction

# PyVista / VTK Import Check
try:
    import pyvista as pv
    import vtk
except ImportError:
    pv = None
    vtk = None

PLUGIN_NAME = "Advanced Rendering"
PLUGIN_VERSION = "2026.02.04"
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
        self.viewer = None 
    
    def set_viewer(self, viewer):
        self.viewer = viewer

    def closeEvent(self, event: QCloseEvent):
        # Save settings before hiding
        if self.viewer and hasattr(self.viewer, 'save_settings'):
            self.viewer.save_settings()
        
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
    
    # 2. If it doesn't exist or was accidentally destroyed from C++, recreate it
    # Check if dialog is a disconnected C++ object using sip (if available) or simple try/except later
    is_deleted = False
    try:
        if dialog and sip.isdeleted(dialog):
            is_deleted = True
    except: pass

    if not dialog or not isinstance(dialog, QDialog) or is_deleted:
        # Ensure the viewer widget exists
        viewer = getattr(main_window, '_adv_rendering_viewer', None)
        try:
            if viewer and sip.isdeleted(viewer):
                viewer = None
        except: pass
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
        viewer.show()  # Make widget visible within the dialog
        
        # Link viewer for save on close
        dialog.set_viewer(viewer)
        
        main_window._adv_graphics_dialog = dialog

    # 3. Show the dialog
    try:
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()
    except RuntimeError:
        # Detected deleted C++ object, force recreate
        main_window._adv_graphics_dialog = None
        load_plugin(main_window)
    except Exception as e:
        logging.error(f"Error showing plugin dialog: {e}")

def initialize(context):
    """
    Plugin entry point. Initializes backend widget and hooks into menus.
    """
    mw = context.get_main_window()
    
    # 1. Cleanup Old Instances (Fix Reloading)
    if hasattr(mw, '_adv_rendering_viewer') and mw._adv_rendering_viewer:
        try:
            old_viewer = mw._adv_rendering_viewer
            if hasattr(old_viewer, '_sync_timer'):
                old_viewer._sync_timer.stop()
            old_viewer.deleteLater()
        except: pass
    
    if hasattr(mw, '_adv_graphics_dialog') and mw._adv_graphics_dialog:
        try:
            mw._adv_graphics_dialog.close()
            mw._adv_graphics_dialog.deleteLater()
        except: pass
        mw._adv_graphics_dialog = None  # Ensure ID is cleared
    
    # 2. Create Background Widget (State Holder)
    viewer = AdvancedGraphicsWidget(mw)
    mw._adv_rendering_viewer = viewer

    # Ensure settings are loaded to check toolbar preference
    viewer.load_settings()

    # 2. Register Menu Action - DISABLED (causes menu artifacts)
    context.add_menu_action("Settings/Advanced Graphics Settings", lambda: load_plugin(mw))


    # 4. Register PBR Styles & Draw Overrides
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
                
                # Apply Custom Lighting Override
                v.update_lights()
                
                # Sync UI to match the new style
                v.sync_style_ui(getattr(mw_obj, 'current_3d_style', ''))
                    
        return drawer

    # Register PBR variants
    styles = [
        ("Ball & Stick (Advanced Rendering)", "Ball & Stick (Advanced Rendering)"),
        ("CPK (Advanced Rendering)", "CPK (Advanced Rendering)"),
        ("Wireframe (Advanced Rendering)", "Wireframe (Advanced Rendering)"),
        ("Stick (Advanced Rendering)", "Stick (Advanced Rendering)")
    ]

    for display_name, style_key in styles:
        context.register_3d_style(style_key, make_style_drawer(style_key.split(' (')[0].lower().replace(' & ', '_and_'), force_pbr=True))

    # Styles are now automatically added to the menu by the main application via register_3d_style.

# --- MAIN WIDGET ---

class AdvancedGraphicsWidget(QWidget):
    def __init__(self, parent_window):
        super().__init__(parent_window)
        self.hide()  # Prevent gray bar from appearing before widget is added to layout
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
        #self.show_toolbar_button = True
        self._init_default_presets()
        
        self.init_ui()
        self.update_preset_combo()
        
        # Delay load to ensure UI is ready
        QTimer.singleShot(100, self.load_settings)

        # Force Lighting Update at startup (Fixes lighting not applied bug)
        # Force Lighting Update at startup (Fixes lighting not applied bug)
        # Use robust retry logic to wait for Plotter/Renderer
        self._init_attempts = 0
        self._max_init_attempts = 20 # 20 * 250ms = 5 seconds
        QTimer.singleShot(250, self._retry_init_lighting)
        
        # Save on exit
        QCoreApplication.instance().aboutToQuit.connect(self.save_settings)
        
        # Polling for external style changes
        self._last_polled_style = ""
        self._sync_timer = QTimer(self)
        self._sync_timer.timeout.connect(self._poll_style_change)
        self._sync_timer.start(500)
        
    def _retry_init_lighting(self):
        """Robustly retries lighting initialization until renderer is ready."""
        # Check if we can access the renderer
        ready = False
        try:
             plotter = self.safe_plotter
             if plotter and hasattr(plotter, 'renderer') and plotter.renderer:
                 ready = True
        except:
             pass
        
        if ready:
             # Apply everything
             self.update_lights()
             self.update_atoms_pbr()
             # Enforce scene effects for persistence
             self.sync_style_ui(getattr(self.mw, 'current_3d_style', ''))
             
             logging.info("Advanced Rendering: Lighting & Effects Initialized Successfully.")
        else:
             self._init_attempts += 1
             if self._init_attempts < self._max_init_attempts:
                 QTimer.singleShot(250, self._retry_init_lighting)
             else:
                 logging.warning("Advanced Rendering: Timed out waiting for renderer initialization.")

        # Apply persisted style on startup -> REMOVED as per request "default shou be normal ball and stick"
        # QTimer.singleShot(200, self._apply_persisted_style)

    # def _apply_persisted_style(self):
    #     """Loads and applies the last used 3D style from plugin settings."""
    #     try:
    #         saved_style = self.presets.get("last_active_style", "")
    #         current = getattr(self.mw, 'current_3d_style', '')
    #         if saved_style and saved_style != current:
    #             if hasattr(self.mw, 'set_3d_style'):
    #                 self.mw.set_3d_style(saved_style)
    #     except Exception as e:
    #         logging.warning(f"Failed to apply persisted style: {e}")

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

        light_group.setLayout(light_group_layout)
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
        
        # Apply/Clear Removed as per request
        # button_layout = QHBoxLayout()
        # btn_apply = QPushButton("Force Apply")
        # btn_apply.setToolTip("Switch to Advanced Rendering mode for current style")
        # btn_apply.clicked.connect(self.force_enable_advanced_style)
        # button_layout.addWidget(btn_apply)
        
        # btn_clear = QPushButton("Reset Actors")
        # btn_clear.clicked.connect(self.clear_atom_settings)
        # button_layout.addWidget(btn_clear)
        # atoms_layout.addLayout(button_layout)
        
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
        
        atoms_layout.addStretch()
        tab_atoms.setLayout(atoms_layout)
        tabs.addTab(tab_atoms, "Atoms")
        
        layout.addWidget(tabs)
        
        # Toolbar Preference
        #self.check_toolbar_btn = QCheckBox("Show Button in Plugin Toolbar (Requires Restart)")
        #self.check_toolbar_btn.setChecked(self.show_toolbar_button)
        #self.check_toolbar_btn.toggled.connect(self.on_toolbar_btn_toggled)
        #layout.addWidget(self.check_toolbar_btn)
        
        self.setLayout(layout)

    #def on_toolbar_btn_toggled(self, checked):
    #    self.show_toolbar_button = checked

    # --- CORE LOGIC ---

    # def force_enable_advanced_style(self):
    #     """Switches the current style to its Advanced Rendering counterpart."""
    #     curr = getattr(self.mw, 'current_3d_style', 'ball_and_stick')
    #     # If already advanced, just update settings
    #     if "(Advanced Rendering)" in curr:
    #         self.update_atoms_pbr()
    #         return
            
    #     # Map standard styles to advanced styles
    #     style_map = {
    #         'ball_and_stick': "Ball & Stick (Advanced Rendering)",
    #         'cpk': "CPK (Advanced Rendering)",
    #         'wireframe': "Wireframe (Advanced Rendering)",
    #         'stick': "Stick (Advanced Rendering)"
    #     }
        
    #     # Try to find a match, defaulting to Ball & Stick Advanced if unknown
    #     target_style = style_map.get(curr, "Ball & Stick (Advanced Rendering)")
        
    #     # Set the style on the main window
    #     if hasattr(self.mw, 'set_3d_style'):
    #         self.mw.set_3d_style(target_style)
    #         # Ensure the checkbox is checked essentially by the style change, 
    #         # but we can also sync immediately just in case.
    #         self.sync_style_ui(target_style)
            
    def _poll_style_change(self):
        """Monitors current_3d_style and updates UI if changed."""
        if not self.mw: return
        curr = getattr(self.mw, 'current_3d_style', '')
        if curr != self._last_polled_style:
            self._last_polled_style = curr
            self.sync_style_ui(curr)
            # Persist the current style
            self.presets["last_active_style"] = curr
            self.save_settings()

    def apply_pbr_forced(self):
        """Used by PBR-specific styles to force PBR on without changing UI check state permanently."""
        self.update_atoms_pbr(force_pbr=True)

    # --- 修正版 AdvancedGraphicsWidget メソッド ---

    @property
    def safe_plotter(self):
        """
        Plotterを安全に取得するヘルパー。
        初期化時にPlotterがまだない場合でも、後から参照できるようにする。
        """
        if self.plotter:
            return self.plotter
        # まだ取得できていない場合は親から再取得を試みる
        if hasattr(self.mw, 'plotter'):
            self.plotter = self.mw.plotter
        return self.plotter

    def update_atoms_pbr(self, force_pbr=False):
        """Applies PBR settings to relevant actors. (Fixed Version)"""
        plotter = self.safe_plotter
        if not plotter or not hasattr(plotter, 'renderer'): return

        try:
            curr_style = getattr(self.mw, 'current_3d_style', '')
            is_advanced = "(Advanced Rendering)" in curr_style
            
            # --- FIX: Strict Policy for PBR ---
            if not is_advanced:
                 # Strict requirement: If not Advanced, PBR is Force Disable.
                 should_use_pbr = False
            else:
                 # Standard logic
                 should_use_pbr = self.use_atom_pbr or force_pbr
            
            # --- FIX 1: シグナルループ防止 ---
            if hasattr(self, 'check_atom_pbr'):
                self.check_atom_pbr.blockSignals(True) # シグナルを遮断
                self.check_atom_pbr.blockSignals(False) # 解除

            # --- FIX 2: アクター取得の互換性確保 ---
            actors_list = []
            if hasattr(plotter.renderer, 'actors'):
                # PyVistaのバージョンによって dict だったり list だったりするため吸収
                acts = plotter.renderer.actors
                if isinstance(acts, dict):
                    actors_list = list(acts.values())
                else:
                    actors_list = list(acts)
            elif hasattr(plotter.renderer, 'GetActors'):
                actors_list = plotter.renderer.GetActors()
            
            for actor in actors_list:
                if actor is None: continue
                
                # --- FIX 3: 厳密な型チェック ---
                # テキストやスカラバーなどの非ジオメトリアクターを除外
                if vtk and not actor.IsA("vtkActor"):
                    continue

                prop = actor.GetProperty()
                if not prop: continue

                # プロパティ変更の適用
                if should_use_pbr:
                    if hasattr(prop, 'SetInterpolationToPBR'):
                        prop.SetInterpolationToPBR()
                        prop.SetMetallic(self.atom_metallic)
                        prop.SetRoughness(self.atom_roughness)
                else:
                    # Phongに戻す (PBR非対応環境への配慮)
                    if hasattr(prop, 'SetInterpolationToPhong'):
                        prop.SetInterpolationToPhong()
                        # Phong用デフォルトリセット
                        prop.SetMetallic(0.0) 
                        prop.SetRoughness(0.5) 
                        # Specularなどもリセットした方が安全だが一旦最小限に
            
            plotter.render()
        except Exception as e:
            logging.error(f"PBR Update Error: {e}")


    # --- EVENT HANDLERS ---

    def on_atom_pbr_toggled(self, checked):
        self.use_atom_pbr = checked
        self.slider_atom_metallic.setEnabled(checked)
        self.slider_atom_roughness.setEnabled(checked)
        self.update_atoms_pbr()
        self.save_settings()
        
    def on_atom_metallic_changed(self, val):
        self.atom_metallic = val / 100.0
        if self.use_atom_pbr:
            self.update_atoms_pbr()
            self.save_settings()
        
    def on_atom_roughness_changed(self, val):
        self.atom_roughness = val / 100.0
        if self.use_atom_pbr:
             self.update_atoms_pbr()
             self.save_settings()

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
    
    def _clean_render_pipeline(self):
        """
        重要: パス切り替え時に発生するVTKエラーを防ぐため、
        一度レンダリングパイプラインを完全にリセットする。
        """
        if not self.plotter or not hasattr(self.plotter, 'renderer'): return
        
        try:
            # 既存のパスを強制解除 (これがReleaseGraphicsResourcesエラーを防ぐ鍵)
            if hasattr(self.plotter.renderer, 'SetPasses'):
                self.plotter.renderer.SetPasses(None)
            
            # Note: Do NOT call render() here. Rendering with None passes can cause
            # VTK to complain about missing resources or invalid state.
            # self.plotter.render() 
        except Exception as e:
            logging.warning(f"Pipeline clean error: {e}")

    def on_shadows_toggled(self, checked):
        if not self.plotter: return
        
        # PBRスタイル以外で強制的にチェックされた場合のガード
        curr_style = getattr(self.mw, 'current_3d_style', '')
        if checked and "(Advanced Rendering)" not in curr_style and not self.use_atom_pbr:
             self.check_shadows.blockSignals(True)
             self.check_shadows.setChecked(False)
             self.check_shadows.blockSignals(False)
             return

        self.use_shadows = checked
        
        # 排他制御
        if checked:
            self._disable_conflicting_effects(exclude="shadows")
            # ★修正: ここでパイプラインを洗浄して再有効化を可能にする
            self._clean_render_pipeline()

        try:
            if checked:
                self.plotter.enable_shadows()
                self.update_lights() # 影有効化後はライト更新が必須
            else:
                self.plotter.disable_shadows()
            
            # AA等の復帰
            if self.use_aa:
                self.plotter.enable_anti_aliasing()

            self.plotter.render()
        except Exception as e:
            logging.warning(f"Shadow error: {e}")
        self.save_settings()

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
        # --- FIX: Strict Light Control ---
        # Default styles (SRC) must keep their original lighting.
        # Only override lighting if we are definitely in an Advanced Requesting style.
        curr_style = getattr(self.mw, 'current_3d_style', '')
        if "(Advanced Rendering)" not in curr_style:
            return

        plotter = self.safe_plotter
        if not plotter or not hasattr(plotter, 'renderer'): return
        
        try:
            # Calculate new position (Spherical -> Cartesian) relative to Camera
            # r = distance
            r = 100.0 
            
            # Convert degrees to radians
            rad_azi = math.radians(self.light_azimuth)
            rad_ele = math.radians(self.light_elevation)
            
            # Coordinate System for Camera Light:
            # Y is UP (Vertical on screen)
            # X is RIGHT (Horizontal on screen)
            # Z is Depth (Towards user or into screen depending on handedness, usually +Z is towards viewer in Camera space)
            
            # Spherical conversion:
            # y = r * sin(ele)
            # x = r * cos(ele) * sin(azi)
            # z = r * cos(ele) * cos(azi)
            
            # Adjust Z to ensure light comes from "front" (viewers side) when azimuth is 0
            # If Z is "towards viewer", strictly positive Z is front.
            
            y = r * math.sin(rad_ele)
            h_proj = r * math.cos(rad_ele)
            x = h_proj * math.sin(rad_azi)
            z = h_proj * math.cos(rad_azi)

            lights = plotter.renderer.GetLights()
            
            # Ensure at least one light exists
            if lights.GetNumberOfItems() == 0:
                new_light = vtk.vtkLight()
                new_light.SetLightTypeToCameraLight()
                new_light.SetSwitch(True)
                plotter.renderer.AddLight(new_light)
                lights = plotter.renderer.GetLights() # Refresh collection

            # Modify the first light to be our Main Light
            # Disable others to prevent conflict if "Override" is implied
            
            main_light = None
            
            is_first = True
            for light in lights:
                if is_first:
                    main_light = light
                    is_first = False
                    
                    # Force Type to CameraLight so it is fixed to Camera Angle
                    if main_light.GetLightType() != vtk.VTK_LIGHT_TYPE_CAMERA_LIGHT:
                        main_light.SetLightTypeToCameraLight()
                    
                    # Position is relative to Camera (0,0,0 is Focus/Eye depending on impl, usually Eye)
                    # For CameraLight, FocalPoint is also relative.
                    # We want it pointing at the center of view (0,0,0) from (x,y,z)
                    main_light.SetPosition(x, y, z)
                    main_light.SetFocalPoint(0, 0, 0)
                    
                    main_light.SetIntensity(self.light_intensity)
                    main_light.SetSwitch(True)
                    main_light.SetPositional(True) 
                else:
                    # Disable other lights to ensure "Override"
                    light.SetSwitch(False)
            
            plotter.render()
        except Exception as e:
            logging.error(f"Lighting Update Error: {e}")

    def on_light_changed(self, val):
        self.light_intensity = val / 100.0
        self.update_lights()

    def on_ssao_toggled(self, checked):
        # Strict Policy
        curr_style = getattr(self.mw, 'current_3d_style', '')
        if checked and "(Advanced Rendering)" not in curr_style:
             self.check_ssao.blockSignals(True)
             self.check_ssao.setChecked(False)
             self.check_ssao.blockSignals(False)
             return

        self.use_ssao = checked
        if checked:
            self._disable_conflicting_effects(exclude="ssao")
        try:
            if checked: self.plotter.enable_ssao()
            else: self.plotter.disable_ssao()
            self.plotter.render()
        except Exception: pass
    
    def on_depth_peeling_toggled(self, checked):
        # Strict Policy
        curr_style = getattr(self.mw, 'current_3d_style', '')
        if checked and "(Advanced Rendering)" not in curr_style:
             self.check_depth.blockSignals(True)
             self.check_depth.setChecked(False)
             self.check_depth.blockSignals(False)
             return

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
        # Strict Policy
        curr_style = getattr(self.mw, 'current_3d_style', '')
        if checked and "(Advanced Rendering)" not in curr_style:
             self.check_edl.blockSignals(True)
             self.check_edl.setChecked(False)
             self.check_edl.blockSignals(False)
             return

        self.use_edl = checked
        if checked:
            self._disable_conflicting_effects(exclude="edl")
            # ★修正: パイプライン洗浄
            self._clean_render_pipeline()
        
        self.slider_edl.setEnabled(checked)
        self.apply_edl()
        
    def on_edl_strength_changed(self, val):
        self.edl_strength = val / 100.0
        if self.use_edl:
             self.update_edl_strength()
             self.save_settings()

    def update_edl_strength(self):
        """スライダーの値をEDLパスに適用"""
        if not self.plotter or not self.use_edl: return
        try:
            # PyVistaが保持しているEDLパスインスタンスにアクセス
            # バージョンによって _edl_pass または edl_pass の可能性があるため両方チェック
            edl_pass = getattr(self.plotter.renderer, '_edl_pass', None)
            if not edl_pass:
                edl_pass = getattr(self.plotter.renderer, 'edl_pass', None)

            if edl_pass:
                # ★修正: Radius(半径)だけでなくDist(距離)も変更しないと見た目が変わらない
                # slider: 1-100
                # Radius: ぼかし範囲 (通常 1.0 - 10.0)
                # Dist: 深度の強調度合い (通常 0.0001 - 0.01 程度だがシーンスケールによる)
                
                # 半径の設定
                if hasattr(edl_pass, 'SetRadius'):
                    edl_pass.SetRadius(self.edl_strength * 10.0)
                
                # 距離の設定 (これを追加することで凹凸がはっきりする)
                if hasattr(edl_pass, 'SetDist'):
                    # シーンによって適切な値は異なるが、スライダーで調整できるようにする
                    # 1.0は大きすぎる場合が多いのでスケーリングする
                    edl_pass.SetDist(self.edl_strength * 5.0) 

            self.plotter.render()
        except Exception as e:
            logging.warning(f"EDL Update Error: {e}")
        
    def apply_edl(self):
        if not self.plotter: return
        try:
            if self.use_edl:
                self.plotter.enable_eye_dome_lighting()
                self.update_edl_strength() # 値を適用
            else:
                self.plotter.disable_eye_dome_lighting()
            
            if self.use_aa:
                self.plotter.enable_anti_aliasing()

            self.plotter.render()
        except Exception as e:
            logging.warning(f"EDL Apply error: {e}")

    def closeEvent(self, event):
        """終了時に特殊効果をOFFにしてリソースを解放する"""
        try:
            if self.plotter:
                self.plotter.disable_eye_dome_lighting()
                self.plotter.disable_shadows()
                self.plotter.disable_ssao()
        except:
            pass
        super().closeEvent(event)

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
            "atom_roughness": self.atom_roughness,
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
        """Syncs widget state and scene effects with active style."""
        is_active_pbr = "(Advanced Rendering)" in active_style_name
        
        # --- FIX: Strictly enforce PBR Checkbox/State based on Style ---
        # This prevents PBR from "leaking" into Standard styles
        if hasattr(self, 'check_atom_pbr'):
            self.check_atom_pbr.blockSignals(True)
            if is_active_pbr:
                self.check_atom_pbr.setChecked(True)
                self.use_atom_pbr = True
                self.slider_atom_metallic.setEnabled(True)
                self.slider_atom_roughness.setEnabled(True)
            else:
                self.check_atom_pbr.setChecked(False)
                self.use_atom_pbr = False
                self.slider_atom_metallic.setEnabled(False)
                self.slider_atom_roughness.setEnabled(False)
            self.check_atom_pbr.blockSignals(False)
        
        # Apply the PBR state immediately
        self.update_atoms_pbr()

        # Enforce Scene Effects based on whether we are in "Advanced" mode or not
        self._enforce_scene_state(enable=is_active_pbr)


    def _enforce_scene_state(self, enable=True):
        """
        If enable=True: Apply settings based on internal state (self.use_shadows, etc.)
        If enable=False: Force disable all advanced effects (Revert to System Defaults)
        """
        if not self.plotter: return

        try:
            # Shadows
            if enable and self.use_shadows:
                self.plotter.enable_shadows()
            else:
                self.plotter.disable_shadows()

            # SSAO
            if enable and self.use_ssao:
                self.plotter.enable_ssao()
            else:
                self.plotter.disable_ssao()
            
            # EDL
            if enable and self.use_edl:
                self.plotter.enable_eye_dome_lighting()
                # Restore strength
                if hasattr(self.plotter.renderer, '_edl_pass'):
                     edl_pass = self.plotter.renderer._edl_pass
                     if edl_pass and hasattr(edl_pass, 'SetEDLStrength'):
                         edl_pass.SetEDLStrength(self.edl_strength)
            else:
                self.plotter.disable_eye_dome_lighting()

            # Anti-Aliasing (Maybe safer to leave if system uses it? But user said "system setting")
            # Usually PyVista enables AA by default if asked. Let's assume plugin controls explicit overrides.
            if enable and self.use_aa:
                 self.plotter.enable_anti_aliasing()
            elif not enable:
                 # Don't strictly disable if system might want it, but usually safe to disable explicit override
                 self.plotter.disable_anti_aliasing()

            # Depth Peeling (Transparency)
            if enable and self.use_depth_peeling:
                self.plotter.enable_depth_peeling()
            elif not enable:
                self.plotter.disable_depth_peeling()
                
            self.plotter.render()
        except Exception as e:
            logging.warning(f"Scene Sync Error: {e}")

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