import os
import json
import logging
from functools import partial
from PyQt6.QtWidgets import QMenu, QDialog, QVBoxLayout, QCheckBox, QDialogButtonBox, QLabel, QApplication
from PyQt6.QtGui import QAction
from PyQt6.QtCore import QTimer
from rdkit import Chem
from rdkit.Chem import AllChem

PLUGIN_NAME = "Structural Updater"
PLUGIN_VERSION = "20260128"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Applies 2D structural changes to 3D conformation without full re-embedding."

# Global reference to hold the original methods to prevent GC and allow restoration
_ORIGINAL_METHODS = {}
_PLUGIN_INSTANCE = None

class SettingsDialog(QDialog):
    def __init__(self, parent=None, current_enabled=True):
        super().__init__(parent)
        self.setWindowTitle("Structural Updater Settings")
        self.enabled = current_enabled
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        
        self.chk_enable = QCheckBox("Enable Plugin")
        self.chk_enable.setChecked(self.enabled)
        layout.addWidget(self.chk_enable)
        
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)
        
        self.setLayout(layout)
        
    def accept(self):
        self.enabled = self.chk_enable.isChecked()
        super().accept()

class StructuralUpdaterPlugin:
    def __init__(self, context):
        self.context = context
        self.mw = context.get_main_window()
        # Single script: settings file is named after the script (e.g. script.json)
        # self.plugin_dir is not really needed if we just use __file__
        base_path = os.path.splitext(os.path.abspath(__file__))[0]
        self.settings_file = base_path + ".json"
        
        self.enabled = True
        
        # State to track if we are in "Apply Mode" vs "Convert Mode"
        self.apply_mode_active = False
        
        # Load settings
        self.load_settings()
        
        # Monkey patch MainWindow
        self.patch_mainwindow()
        
        # Add Menu Actions
        self.context.add_menu_action("Settings/Structural Updater...", self.open_settings)
        
        # Auto-detection timer
        self.timer = QTimer(self.mw)
        self.timer.timeout.connect(self.check_state)
        self.timer.start(1000) # Check every 1 second

    def check_state(self):
        if not self.enabled: return
        
        # Do not interfere if calculation is running (Halt button active)
        if self.mw.convert_button.text() == "Halt conversion":
            return
            
        # Check if 3D molecule exists and has original_atom_id props
        has_3d = False
        if self.mw.current_mol and self.mw.current_mol.GetNumAtoms() > 0:
            # Check for at least one original_atom_id
            for atom in self.mw.current_mol.GetAtoms():
                if atom.HasProp("_original_atom_id"):
                    has_3d = True
                    break
        
        # Update mode logic
        if has_3d:
            if not self.apply_mode_active:
                self.apply_mode_active = True
                self.mw.convert_button.setText("Apply 2D Changes to 3D")
        else:
            if self.apply_mode_active:
                self.apply_mode_active = False
                self.mw.convert_button.setText("Convert 2D to 3D")
        
    def load_settings(self):
        try:
            if os.path.exists(self.settings_file):
                with open(self.settings_file, 'r') as f:
                    data = json.load(f)
                    self.enabled = data.get('enabled', True)
            else:
                # Create default settings if file doesn't exist
                self.save_settings()
        except Exception as e:
            print(f"[{PLUGIN_NAME}] Error loading settings: {e}")

    def save_settings(self):
        try:
            data = {'enabled': self.enabled}
            with open(self.settings_file, 'w') as f:
                json.dump(data, f)
        except Exception as e:
            print(f"[{PLUGIN_NAME}] Error saving settings: {e}")

    def open_settings(self):
        dlg = SettingsDialog(self.mw, self.enabled)
        if dlg.exec():
            self.enabled = dlg.enabled
            self.save_settings()
            
            # If disabled, restore original button text immediately
            if not self.enabled:
                self.mw.convert_button.setText("Convert 2D to 3D")
                
            status = "Enabled" if self.enabled else "Disabled"
            self.mw.statusBar().showMessage(f"Structural Updater: {status}")
            print(f"[{PLUGIN_NAME}] Enabled: {self.enabled}")

    def patch_mainwindow(self):
        # Store original methods
        if 'trigger_conversion' not in _ORIGINAL_METHODS:
            _ORIGINAL_METHODS['trigger_conversion'] = self.mw.trigger_conversion
        if 'on_calculation_finished' not in _ORIGINAL_METHODS:
            _ORIGINAL_METHODS['on_calculation_finished'] = self.mw.on_calculation_finished
        # Do not override show_convert_menu as per user request

        # Replace methods
        self.mw.trigger_conversion = self.new_trigger_conversion
        self.mw.on_calculation_finished = self.new_on_calculation_finished
        
        # CRITICAL: Reconnect the button signal!
        # The existing connection points to the OLD trigger_conversion function object.
        try:
            self.mw.convert_button.clicked.disconnect()
        except Exception:
            pass 
        self.mw.convert_button.clicked.connect(self.new_trigger_conversion)

    def new_trigger_conversion(self):
        """Replacement for MainWindow.trigger_conversion"""
        # Check if a temporary conversion mode is set (via right-click menu)
        # If so, we must respect it and bypass our "Apply" logic.
        temp_mode = getattr(self.mw, "_temp_conv_mode", None)
        
        # 1. Check if we should use our logic
        # Only use our logic if Plugin is Enabled AND Apply Mode is active AND No temp mode override
        if self.enabled and self.apply_mode_active and not temp_mode:
            # Run our custom "Apply 2D to 3D" logic
            self.apply_changes_to_3d()
        else:
            # Run original conversion logic
            _ORIGINAL_METHODS['trigger_conversion']()

    def new_on_calculation_finished(self, result):
        """Replacement for MainWindow.on_calculation_finished"""
        # 1. Run original (handles UI reset, loading mol into mw.current_mol, etc.)
        _ORIGINAL_METHODS['on_calculation_finished'](result)
        
        # Force reconnection to our method because original method resets it to self.trigger_conversion
        # ensuring it points to our wrapper
        try:
            self.mw.convert_button.clicked.disconnect()
        except Exception:
            pass
        self.mw.convert_button.clicked.connect(self.new_trigger_conversion)
        
        # 2. If plugin is enabled and conversion was successful, switch to "Apply Mode"
        if self.enabled:
            # Check if we have a valid molecule now
            if self.mw.current_mol and self.mw.current_mol.GetNumAtoms() > 0:
                self.apply_mode_active = True
                self.mw.convert_button.setText("Apply 2D Changes to 3D")
            else:
                self.apply_mode_active = False
                self.mw.convert_button.setText("Convert 2D to 3D")

    def force_full_conversion(self):
        """Forces a standard conversion, resetting the Apply Mode."""
        self.apply_mode_active = False 
        # Button text will be updated by on_calculation_finished later, 
        # but good to reset it visually now or let the conversion process handle it.
        # Original trigger_conversion clears plotter etc.
        _ORIGINAL_METHODS['trigger_conversion']()

    def apply_changes_to_3d(self):
        """
        The core logic: 
        1. Get current 2D structure (New)
        2. Get current 3D structure (Old)
        3. Match atoms via _original_atom_id
        4. "Direct" Construction: Manually build conformer from old coords + random new atoms
        5. Constrained Optimization
        """
        self.mw.statusBar().showMessage("Applying 2D changes to 3D structure...")
        QApplication.processEvents() # Ensure message is shown
        
        # 1. Get New 2D Mol
        new_mol = self.mw.data.to_rdkit_mol()
        if not new_mol:
            self.mw.statusBar().showMessage("Error: Invalid 2D structure.")
            return

        # 2. Get Old 3D Mol
        old_mol = self.mw.current_mol
        if not old_mol:
            # If no 3D mol exists, fall back to full conversion
            self.force_full_conversion()
            return
            
        # 3. Create Coordinate Map
        # coordMap is {atomIdx: Point3D} for the NEW molecule.
        coord_map = {}
        
        # We need to map New Mol Atom Idx -> Original ID -> Old Mol Atom Idx -> Coordinates
        
        # Helper to safely get original atom ID
        def get_original_id(at):
            try:
                if at.HasProp("_original_atom_id"):
                    return at.GetIntProp("_original_atom_id")
            except Exception:
                pass
            return None

        # Build map of Old Mol: Original ID -> Atom Idx
        old_id_to_idx = {}
        for atom in old_mol.GetAtoms():
            oid = get_original_id(atom)
            if oid is not None:
                old_id_to_idx[oid] = atom.GetIdx()
        
        # Get Conformer from Old Mol
        try:
            old_conf = old_mol.GetConformer()
        except:
            self.force_full_conversion()
            return

        # Prepare New Mol
        new_mol = Chem.AddHs(new_mol) # Important: Add Hydrogens for 3D embedding
        
        # Helper to get heavy neighbors info: list of (neighbor_id, bond_type, atomic_num)
        def get_env_signature(at):
            sig = []
            for bond in at.GetBonds():
                nbr = bond.GetOtherAtom(at)
                if nbr.GetAtomicNum() == 1: continue # Skip H
                
                nid = -1
                if nbr.HasProp("_original_atom_id"):
                    nid = nbr.GetIntProp("_original_atom_id")
                    
                sig.append((nid, bond.GetBondType(), nbr.GetAtomicNum()))
            sig.sort()
            return sig

        # Build Coord Map for New Mol
        matches_count = 0
        for atom in new_mol.GetAtoms():
            oid = get_original_id(atom)
            if oid is not None and oid in old_id_to_idx:
                old_idx = old_id_to_idx[oid]
                old_atom = old_mol.GetAtomWithIdx(old_idx)
                
                # Check if central atom identity changed (e.g. C -> N)
                atom_conserved = (atom.GetAtomicNum() == old_atom.GetAtomicNum())
                
                if atom_conserved:
                    # Heavy Atom: Always try to fix it if it's the same element.
                    # "Minimum matching using atom id" - we prioritize keeping the scaffold matching.
                    pos = old_conf.GetAtomPosition(old_idx)
                    coord_map[atom.GetIdx()] = pos
                    matches_count += 1
                    
                    # Restore Hydrogen Fixation (User Request)
                    # Only fix Hydrogens if the local environment (connectivity/bonds) is conserved.
                    
                    new_env = get_env_signature(atom)
                    old_env = get_env_signature(old_atom)
                    
                    if new_env == old_env:
                        # Environment conserved: Fix Hydrogens as well
                        new_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
                        old_neighbors = [nbr for nbr in old_atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
                        
                        count_h = min(len(new_neighbors), len(old_neighbors))
                        # Naive mapping: new[i] -> old[i]. 
                        # Since Hs are indistinguishable, this preserves the geometric "slots".
                        for i in range(count_h):
                            h_new = new_neighbors[i]
                            h_old = old_neighbors[i]
                            h_pos = old_conf.GetAtomPosition(h_old.GetIdx())
                            coord_map[h_new.GetIdx()] = h_pos
        
        # If no matching data or too few matches to define orientation, generate from scratch
        # We generally need at least 3 points to define 3D orientation.
        min_matches = 3
        if new_mol.GetNumAtoms() < 3:
            min_matches = 1 # Very small molecule
            
        if matches_count < min_matches:
            self.mw.statusBar().showMessage(f"Notice: Too few matching atoms ({matches_count}). Performing full conversion.")
            self.force_full_conversion()
            return
            
        try:
            # 4. Constrained Embedding with RDKit
            # We use EmbedMolecule with coordMap to guide the embedding.
            # useRandomCoords=True is critical to avoid O(N^3) freezing on large rings.
            
            # Primary attempt: Strict constraints, good box size, limit attempts to prevent freeze
            res = AllChem.EmbedMolecule(new_mol, 
                                      coordMap=coord_map, 
                                      useRandomCoords=True, 
                                      randomSeed=0xf00d,
                                      boxSizeMult=2.0,
                                      maxAttempts=5)
            
            if res != 0:
                 # Fallback 1: Less strict (relaxed box size, clear confs)
                 res = AllChem.EmbedMolecule(new_mol, 
                                           coordMap=coord_map, 
                                           useRandomCoords=True, 
                                           randomSeed=0xf00d,
                                           clearConfs=True,
                                           useBasicKnowledge=False,
                                           boxSizeMult=4.0,
                                           maxAttempts=5)
            
            if res != 0:
                 # Fallback 2: Free embedding (NO coordMap), then Snap
                 # Guarantees some structure is generated
                 res = AllChem.EmbedMolecule(new_mol, 
                                           useRandomCoords=True, 
                                           randomSeed=0xf00d,
                                           boxSizeMult=2.0,
                                           maxAttempts=5)
                 
                 if res == 0:
                     # Free embedding succeeded. We rely on the Snap step below to fix positions.
                     pass 
                 else:
                     raise ValueError("Failed to generate any 3D conformer.")

            # C. Snap Matches to Exact Old Positions (Hard Constraint Enforcement)
            # Embedding with coordMap is a "guide" (soft constraint).
            # We now Enforce hard constraints (overwrite positions).
            if new_mol.GetNumConformers() > 0:
                pass_conf = new_mol.GetConformer()
                for new_idx, pos in coord_map.items():
                    pass_conf.SetAtomPosition(new_idx, pos)
                
            # D. Optimization (Constrained)
            # Now we run the same constrained optimization as before to clean up the new connections.
            # (Logic continues below...)
            # Optimization (constrained) via UFF
            try:
                # Use UFF ForceField to optimize new atoms while keeping old atoms fixed
                if AllChem.MMFFHasAllMoleculeParams(new_mol):
                     ff = AllChem.MMFFGetMoleculeForceField(new_mol, AllChem.MMFFGetMoleculeProperties(new_mol))
                else:
                     ff = AllChem.UFFGetMoleculeForceField(new_mol)
                
                if ff:
                    # Add fixed point constraints for all matched atoms
                    for atom_idx in coord_map.keys():
                        ff.AddFixedPoint(atom_idx)
                    
                    # Minimize with iteration limit to prevent hanging on complex/strained systems
                    ff.Minimize(maxIts=200)
            except Exception as e:
                # If optimization explodes (e.g. Invariant Violation), fallback is safer
                # print(f"[{PLUGIN_NAME}] Optimization failed: {e}") # Suppress scary error
                self.mw.statusBar().showMessage("Notice: Structure optimization unstable. Re-generating full 3D structure...")
                self.force_full_conversion()
                return

        except Exception as e: # Error during embedding
            print(f"[{PLUGIN_NAME}] Embedding failed: {e}")
            self.force_full_conversion()
            return
            


        # 5. Success - Update UI
        self.mw.on_calculation_finished(new_mol)
        
        self.mw.statusBar().showMessage("Applied 2D changes to 3D structure.")

def initialize(context):
    global _PLUGIN_INSTANCE
    _PLUGIN_INSTANCE = StructuralUpdaterPlugin(context)

def finalize():
    # Restore original methods
    mw = _PLUGIN_INSTANCE.mw
    for name, method in _ORIGINAL_METHODS.items():
        setattr(mw, name, method)
    
    # Restore button signal?
    # It takes effort, but generally resetting methods is enough for cleanup.
    # The button connection might persist until restart, usually acceptable for python plugins.
