#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
POV-Ray Export Plugin for MoleditPy
Exports molecular structures as POV-Ray scene files for high-quality ray-traced rendering
"""

PLUGIN_NAME = "POV-Ray Export"
PLUGIN_VERSION = "2026.02.05"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Export molecular structures as POV-Ray scene files for professional ray-traced rendering"
PLUGIN_DEPENDENCIES = ["rdkit", "numpy", "PyQt6"]


def initialize(context):
    """Initialize the POV-Ray export plugin."""
    context.add_export_action("Export to POV-Ray Scene...", lambda: export_to_povray(context))


def export_to_povray(context):
    """Main export function that generates a POV-Ray scene file."""
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
    import os
    
    # Get the main window and current molecule
    mw = context.get_main_window()
    mol = context.current_molecule
    
    # Check if there's a molecule to export
    if not mol or mol.GetNumAtoms() == 0:
        QMessageBox.warning(mw, "No Molecule", "No molecule loaded. Please create or load a molecule first.")
        return
    
    # Check if the molecule has 3D coordinates
    if mol.GetNumConformers() == 0:
        QMessageBox.warning(mw, "No 3D Structure", 
                          "The molecule has no 3D coordinates. Please generate a 3D structure first.")
        return
    
    # Get default directory from current file path
    default_dir = ""
    default_name = "molecule"
    try:
        if hasattr(mw, 'current_file_path') and mw.current_file_path:
            default_dir = os.path.dirname(mw.current_file_path)
            base_name = os.path.splitext(os.path.basename(mw.current_file_path))[0]
            default_name = base_name
    except Exception:
        pass
    
    default_path = os.path.join(default_dir, default_name) if default_dir else default_name
    
    # Ask user where to save the file
    file_path, _ = QFileDialog.getSaveFileName(
        mw, 
        "Export to POV-Ray Scene", 
        default_path, 
        "POV-Ray Scene Files (*.pov);;All Files (*)"
    )
    
    if not file_path:
        return
    
    # Ensure .pov extension
    if not file_path.lower().endswith('.pov'):
        file_path += '.pov'
    
    try:
        # Generate and save the POV-Ray scene
        scene_content, exp_w, exp_h = generate_povray_scene(mol, mw)
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(scene_content)
        
        mw.statusBar().showMessage(f"POV-Ray scene exported to {file_path}")
        
        # Show success message with instructions
        QMessageBox.information(
            mw, 
            "Export Successful",
            f"POV-Ray scene saved to:\n{file_path}\n\n"
            "To render high quality:\n"
            "1. Install POV-Ray (www.povray.org)\n"
            "2. Open the .pov file in POV-Ray\n"
            "3. IMPORTANT: Select your resolution (e.g., 1920x1080 or better) from the dropdown menu in the top-left area of the POV-Ray window.\n\n"
            "Or use command line for recommended high quality:\n"
            f"povray +W{exp_w} +H{exp_h} +A0.01 +AM2 +R3 {os.path.basename(file_path)}"
        )
        
    except Exception as e:
        QMessageBox.critical(mw, "Export Failed", f"Failed to export POV-Ray scene:\n{str(e)}")
        mw.statusBar().showMessage(f"Export failed: {e}")


def _calculate_double_bond_offset(mol, bond, conf):
    """Calculate the offset direction for a double bond to make it planar."""
    import numpy as np
    begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
    end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
    
    begin_pos = np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx()))
    end_pos = np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
    
    bond_vec = end_pos - begin_pos
    bond_length = np.linalg.norm(bond_vec)
    if bond_length == 0:
        return np.array([0, 0, 1])
    
    bond_unit = bond_vec / bond_length
    
    # Check neighbors
    begin_neighbors = [np.array(conf.GetAtomPosition(n.GetIdx())) for n in begin_atom.GetNeighbors() if n.GetIdx() != bond.GetEndAtomIdx()]
    end_neighbors = [np.array(conf.GetAtomPosition(n.GetIdx())) for n in end_atom.GetNeighbors() if n.GetIdx() != bond.GetBeginAtomIdx()]
    
    normal_candidates = []
    
    for neighbors, pos in [(begin_neighbors, begin_pos), (end_neighbors, end_pos)]:
        for n_pos in neighbors:
            vec = n_pos - pos
            if np.linalg.norm(vec) > 1e-6:
                normal = np.cross(bond_vec, vec)
                norm_len = np.linalg.norm(normal)
                if norm_len > 1e-6:
                    normal_candidates.append(normal / norm_len)
    
    if normal_candidates:
        ref = normal_candidates[0]
        aligned = [n if np.dot(n, ref) >= 0 else -n for n in normal_candidates]
        avg_normal = np.mean(aligned, axis=0)
        avg_norm_len = np.linalg.norm(avg_normal)
        if avg_norm_len > 1e-6:
            avg_normal /= avg_norm_len
            offset_dir = np.cross(bond_unit, avg_normal)
            off_len = np.linalg.norm(offset_dir)
            if off_len > 1e-6:
                return offset_dir / off_len
    
    # Fallback
    v_arb = np.array([0, 0, 1])
    if np.allclose(np.abs(np.dot(bond_unit, v_arb)), 1.0):
        v_arb = np.array([0, 1, 0])
    off_dir = np.cross(bond_unit, v_arb)
    return off_dir / np.linalg.norm(off_dir)


def generate_povray_scene(mol, mw):
    """Generate the POV-Ray scene file content."""
    import math
    import numpy as np
    from rdkit import Chem
    pt = Chem.GetPeriodicTable()
    
    # Get 3D coordinates
    conf = mol.GetConformer()
    
    # Extract actual rendering data from the 3D view
    color_map = getattr(mw, '_3d_color_map', {})
    current_style = getattr(mw, 'current_3d_style', 'ball_and_stick')
    settings = getattr(mw, 'settings', {})
    
    # Get actual atom radii from the glyph source if available
    actual_radii = None
    if hasattr(mw, 'glyph_source') and mw.glyph_source is not None:
        try:
            actual_radii = mw.glyph_source['radii']
        except Exception:
            pass
            
    # Bond parameters from settings
    if current_style == 'wireframe':
        cyl_radius = settings.get('wireframe_bond_radius', 0.01)
        double_radius_factor = settings.get('wireframe_double_bond_radius_factor', 0.8)
        triple_radius_factor = settings.get('wireframe_triple_bond_radius_factor', 0.75)
        double_offset_factor = settings.get('wireframe_double_bond_offset_factor', 3.0)
        triple_offset_factor = settings.get('wireframe_triple_bond_offset_factor', 3.0)
    elif current_style == 'stick':
        cyl_radius = settings.get('stick_bond_radius', 0.15)
        double_radius_factor = settings.get('stick_double_bond_radius_factor', 0.60)
        triple_radius_factor = settings.get('stick_triple_bond_radius_factor', 0.40)
        double_offset_factor = settings.get('stick_double_bond_offset_factor', 1.5)
        triple_offset_factor = settings.get('stick_triple_bond_offset_factor', 1.0)
    else:  # ball_and_stick
        cyl_radius = settings.get('ball_stick_bond_radius', 0.1)
        double_radius_factor = settings.get('ball_stick_double_bond_radius_factor', 0.8)
        triple_radius_factor = settings.get('ball_stick_triple_bond_radius_factor', 0.75)
        double_offset_factor = settings.get('ball_stick_double_bond_offset_factor', 2.0)
        triple_offset_factor = settings.get('ball_stick_triple_bond_offset_factor', 2.0)

    # Background color from settings
    bg_hex = settings.get('background_color', '#4f4f4f')
    try:
        from PyQt6.QtGui import QColor
        q_bg = QColor(bg_hex)
        bg_rgb = (q_bg.redF(), q_bg.greenF(), q_bg.blueF())
    except Exception:
        bg_rgb = (0.31, 0.31, 0.31) # Default gray
    
    # Calculate molecule center and bounds for camera positioning (fallback only)
    atom_positions = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atom_positions.append(np.array([pos.x, pos.y, pos.z]))
    
    if atom_positions:
        center_x = sum(p[0] for p in atom_positions) / len(atom_positions)
        center_y = sum(p[1] for p in atom_positions) / len(atom_positions)
        center_z = sum(p[2] for p in atom_positions) / len(atom_positions)
        
        # Calculate bounds
        max_dist = 0
        for p in atom_positions:
            dist = math.sqrt((p[0] - center_x)**2 + (p[1] - center_y)**2 + (p[2] - center_z)**2)
            max_dist = max(max_dist, dist)
    else:
        center_x = center_y = center_z = 0
        max_dist = 5
    
    # Get current camera and resolution from plotter if available
    camera_data = {
        'location': (center_x, center_y, center_z + max_dist * 3),
        'look_at': (center_x, center_y, center_z),
        'sky': (0, 1, 0),
        'angle': 45,
        'orthographic': False,
        'width': 1920, # Initial fallback
        'height': 1080
    }
    
    if hasattr(mw, 'plotter'):
        # Get resolution
        if hasattr(mw.plotter, 'window_size'):
            camera_data['width'], camera_data['height'] = mw.plotter.window_size
        
        # Get camera
        if hasattr(mw.plotter, 'camera'):
            vtk_camera = mw.plotter.camera
            pos = vtk_camera.GetPosition()
            focal = vtk_camera.GetFocalPoint()
            up = vtk_camera.GetViewUp()
            
            camera_data['location'] = pos
            camera_data['look_at'] = focal
            camera_data['sky'] = up
            
            # Determine if orthographic
            if vtk_camera.GetParallelProjection():
                camera_data['orthographic'] = True
                # For orthographic, we need to scale the "right" and "up" vectors by ParallelScale
                camera_data['parallel_scale'] = vtk_camera.GetParallelScale()
                
            # Get FOV (VTK ViewAngle is vertical FOV)
            camera_data['angle'] = vtk_camera.GetViewAngle()

    # CPK colors for common elements (RGB 0-1 range) - fallback only
    cpk_colors = {
        'H': (1.0, 1.0, 1.0), 'He': (0.85, 1.0, 1.0), 'Li': (0.8, 0.5, 1.0), 'Be': (0.76, 1.0, 0.0),
        'B': (1.0, 0.7, 0.7), 'C': (0.56, 0.56, 0.56), 'N': (0.19, 0.31, 0.97), 'O': (1.0, 0.05, 0.05),
        'F': (0.56, 0.88, 0.31), 'Ne': (0.7, 0.89, 0.96), 'Na': (0.67, 0.36, 0.95), 'Mg': (0.54, 1.0, 0.0),
        'Al': (0.75, 0.65, 0.65), 'Si': (0.94, 0.78, 0.63), 'P': (1.0, 0.65, 0.0), 'S': (1.0, 1.0, 0.19),
        'Cl': (0.12, 0.94, 0.12), 'Ar': (0.5, 0.82, 0.89), 'K': (0.56, 0.25, 0.83), 'Ca': (0.24, 1.0, 0.0),
        'Fe': (0.88, 0.4, 0.2), 'Cu': (0.78, 0.5, 0.2), 'Zn': (0.49, 0.5, 0.69), 'Br': (0.65, 0.16, 0.16),
        'I': (0.58, 0.0, 0.58), 'Au': (1.0, 0.82, 0.14),
    }
    
    # Smart Resolution Scaling
    res_scale = 1.0
    if camera_data['width'] < 1280:
        res_scale = 2.0
    elif camera_data['width'] < 2560:
        res_scale = 1.5
    else:
        res_scale = 1.0
        
    exp_w = int(camera_data['width'] * res_scale)
    exp_h = int(camera_data['height'] * res_scale)
    
    # Extract specular/roughness from settings
    specular_val = settings.get('specular', 0.2)
    specular_power = settings.get('specular_power', 20)
    # Convert specular power to roughness (approximate)
    roughness_val = max(0.01, min(1.0, 1.0 - (specular_power / 100.0)))
    
    # Van der Waals radii (in Angstroms) - fallback only
    vdw_radii = {
        'C': 1.70, 'H': 1.20, 'O': 1.52, 'N': 1.55, 'S': 1.80,
        'P': 1.80, 'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98,
    }
    
    # Start building the scene
    lines = []
    
    # Header
    lines.append('// POV-Ray Scene File')
    lines.append(f'// Generated from MoleditPy - {PLUGIN_NAME} v{PLUGIN_VERSION}')
    lines.append('// Molecular structure with atoms and bonds')
    lines.append('//')
    lines.append('// NOTE: To get high-resolution output, select your desired resolution')
    lines.append('// (e.g., 1280x1024) from the dropdown menu in the POV-Ray GUI.')
    lines.append('//')
    lines.append(f'// Recommended Command Line (Smart Scaling):')
    lines.append(f'// povray +W{exp_w} +H{exp_h} +A0.01 +AM2 +R3 <filename>.pov')
    lines.append('')
    
    # Global settings
    lines.append('#include "colors.inc"')
    lines.append('#include "textures.inc"')
    lines.append('#include "finish.inc"')
    lines.append('#include "rad_def.inc"')
    lines.append('')
    
    lines.append('global_settings {')
    lines.append('  assumed_gamma 1.0')
    lines.append('  max_trace_level 8')
    lines.append('  // Radiosity for realistic lighting')
    lines.append('  radiosity {')
    lines.append('    Rad_Settings(Radiosity_Final, off, off)')
    lines.append('  }')
    lines.append('}')
    lines.append('')
    
    # Background and Environment
    lines.append(f'background {{ color rgb <{bg_rgb[0]:.3f}, {bg_rgb[1]:.3f}, {bg_rgb[2]:.3f}> }}')
    lines.append('')
    
    lines.append('// Sky sphere for reflections - using bozo for smooth, edge-free background')
    lines.append('sky_sphere {')
    lines.append('  pigment {')
    lines.append('    bozo')
    lines.append('    turbulence 0.6')
    lines.append('    color_map {')
    lines.append(f'      [0.0 color rgb <{bg_rgb[0]*0.95:.3f}, {bg_rgb[1]*0.95:.3f}, {bg_rgb[2]*0.95:.3f}>]')
    lines.append(f'      [0.5 color rgb <{bg_rgb[0]:.3f}, {bg_rgb[1]:.3f}, {bg_rgb[2]:.3f}>]')
    lines.append(f'      [1.0 color rgb <{min(1.0, bg_rgb[0]*1.05):.3f}, {min(1.0, bg_rgb[1]*1.05):.3f}, {min(1.0, bg_rgb[2]*1.05):.3f}>]')
    lines.append('    }')
    lines.append('    scale 0.5')
    lines.append('  }')
    lines.append('}')
    lines.append('')
    
    # Camera
    loc = camera_data['location']
    look = camera_data['look_at']
    sky = camera_data['sky']
    
    lines.append('camera {')
    if camera_data['orthographic']:
        lines.append('  orthographic')
        # Scale the view area
        s = camera_data['parallel_scale']
        lines.append(f'  up y*{s*2}')
        lines.append(f'  right -x*{s*2}*image_width/image_height')
    else:
        lines.append('  perspective')
        lines.append(f'  angle {camera_data["angle"]}')
        # Use negative x for right vector to match VTK's right-handed system
        lines.append('  right -x*image_width/image_height')
        
    lines.append(f'  location <{loc[0]:.4f}, {loc[1]:.4f}, {loc[2]:.4f}>')
    lines.append(f'  look_at <{look[0]:.4f}, {look[1]:.4f}, {look[2]:.4f}>')
    lines.append(f'  sky <{sky[0]:.4f}, {sky[1]:.4f}, {sky[2]:.4f}>')
    lines.append('}')
    lines.append('')
    
    # Lighting
    # We'll use a main area light for soft shadows and some fill lights
    light_dist = max_dist * 5
    lines.append('// High-quality Area Light')
    lines.append('light_source {')
    lines.append(f'  <{loc[0] + light_dist}, {loc[1] + light_dist}, {loc[2] + light_dist}>')
    lines.append('  color rgb <1.0, 1.0, 1.0>')
    lines.append('  area_light <5, 0, 0>, <0, 5, 0>, 5, 5')
    lines.append('  adaptive 1')
    lines.append('  jitter')
    lines.append('  circular orient')
    lines.append('}')
    lines.append('')
    
    lines.append('// Fill light')
    lines.append('light_source {')
    lines.append(f'  <{loc[0] - light_dist}, {loc[1]}, {loc[2] - light_dist}>')
    lines.append('  color rgb <0.3, 0.3, 0.3>')
    lines.append('  shadowless')
    lines.append('}')
    lines.append('')
    
    lines.append('// Back light for rim highlights')
    lines.append('light_source {')
    lines.append(f'  <{look[0]}, {look[1] - light_dist}, {look[2] - light_dist}>')
    lines.append('  color rgb <0.2, 0.2, 0.2>')
    lines.append('  shadowless')
    lines.append('}')
    lines.append('')
    
    # Define common finish for atoms
    lines.append('/*')
    lines.append('// Alternative Finish for atoms - Shiny Plastic')
    lines.append('#declare AtomFinish = finish {')
    lines.append('  ambient 0.1')
    lines.append('  diffuse 0.7')
    lines.append(f'  specular {specular_val:.2f}')
    lines.append(f'  roughness {roughness_val:.3f}')
    lines.append('  reflection { 0.05 }')
    lines.append('}')
    lines.append('// Alternative Finish for bonds - Matte')
    lines.append('#declare BondFinish = finish {')
    lines.append('  ambient 0.1')
    lines.append('  diffuse 0.8')
    lines.append(f'  specular {specular_val/2:.2f}')
    lines.append(f'  roughness {roughness_val*2 if roughness_val*2 < 1.0 else 1.0:.3f}')
    lines.append('}')
    lines.append('*/')
    lines.append('')
    
    lines.append('// Finish for atoms - ultra-metallic finish')
    lines.append('#declare AtomFinish = finish {')
    lines.append('  ambient 0.05')
    lines.append('  diffuse 0.4')
    lines.append('  brilliance 2.0')
    lines.append(f'  specular {min(1.0, specular_val * 2):.2f}')
    lines.append(f'  roughness {roughness_val * 0.5:.3f}')
    lines.append('  metallic 0.8')
    lines.append('  reflection {')
    lines.append('    0.3, 0.7')
    lines.append('    fresnel on')
    lines.append('    falloff 1.0')
    lines.append('    exponent 1.0')
    lines.append('    metallic 0.8')
    lines.append('  }')
    lines.append('  conserve_energy')
    lines.append('}')
    lines.append('')
    
    # Define common finish for bonds
    lines.append('// Finish for bonds - sleek chrome-like metal')
    lines.append('#declare BondFinish = finish {')
    lines.append('  ambient 0.1')
    lines.append('  diffuse 0.5')
    lines.append(f'  specular {min(1.0, specular_val):.2f}')
    lines.append(f'  roughness {max(0.001, roughness_val * 0.2):.3f}')
    lines.append('  metallic 0.6')
    lines.append('  reflection { 0.25, 0.5 metallic 0.6 }')
    lines.append('}')
    lines.append('')
    
    # Create atoms
    if current_style != 'wireframe':
        lines.append('// ===== ATOMS =====')
        lines.append('')
        
        # In stick style, terminal atoms with double/triple bonds are split
        atoms_to_draw = [] # list of (pos, radius, color, symbol, idx)
        
        if current_style == 'stick':
            skip_indices = set()
            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                if atom.GetDegree() == 1:
                    bond = atom.GetBonds()[0]
                    bt = bond.GetBondType()
                    if bt in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
                        other_idx = bond.GetBeginAtomIdx() if bond.GetEndAtomIdx() == idx else bond.GetEndAtomIdx()
                        pos_i = np.array(conf.GetAtomPosition(idx))
                        pos_other = np.array(conf.GetAtomPosition(other_idx))
                        bond_vec = pos_i - pos_other
                        bond_unit = bond_vec / np.linalg.norm(bond_vec)
                        
                        # Get offset parameters
                        if bt == Chem.BondType.DOUBLE:
                            off_dir = _calculate_double_bond_offset(mol, bond, conf)
                            off_dist = cyl_radius * double_offset_factor / 2
                            sphere_rad = cyl_radius * double_radius_factor
                            offsets = [off_dir * off_dist, -off_dir * off_dist]
                        else: # TRIPLE
                            v_arb = np.array([0, 0, 1])
                            if np.allclose(np.abs(np.dot(bond_unit, v_arb)), 1.0): v_arb = np.array([0, 1, 0])
                            off_dir = np.cross(bond_unit, v_arb)
                            off_dir /= np.linalg.norm(off_dir)
                            off_dist = cyl_radius * triple_offset_factor
                            sphere_rad = cyl_radius * triple_radius_factor
                            offsets = [np.zeros(3), off_dir * off_dist, -off_dir * off_dist]
                            
                        # Get color
                        color_key = f'atom_{idx}'
                        if color_key in color_map:
                            color = tuple(c / 255.0 for c in color_map[color_key])
                        else:
                            color = cpk_colors.get(atom.GetSymbol(), (0.8, 0.8, 0.8))
                            
                        for off in offsets:
                            atoms_to_draw.append((pos_i + off, sphere_rad, color, atom.GetSymbol(), idx))
                        skip_indices.add(idx)
            
            # Add remaining atoms
            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                if idx in skip_indices: continue
                pos = conf.GetAtomPosition(idx)
                symbol = atom.GetSymbol()
                color_key = f'atom_{idx}'
                color = tuple(c / 255.0 for c in color_map[color_key]) if color_key in color_map else cpk_colors.get(symbol, (0.8, 0.8, 0.8))
                radius = float(actual_radii[idx]) if actual_radii is not None and idx < len(actual_radii) else cyl_radius
                atoms_to_draw.append((np.array([pos.x, pos.y, pos.z]), radius, color, symbol, idx))
        else:
            # Standard ball_and_stick or cpk
            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                pos = conf.GetAtomPosition(idx)
                symbol = atom.GetSymbol()
                color_key = f'atom_{idx}'
                color = tuple(c / 255.0 for c in color_map[color_key]) if color_key in color_map else cpk_colors.get(symbol, (0.8, 0.8, 0.8))
                
                if actual_radii is not None and idx < len(actual_radii):
                    radius = float(actual_radii[idx])
                else:
                    if current_style == 'cpk':
                        scale = settings.get('cpk_atom_scale', 1.0)
                        try:
                            r = pt.GetRvdw(pt.GetAtomicNumber(symbol))
                            radius = (r if r > 0.1 else 1.5) * scale
                        except:
                            radius = 1.5 * scale
                    else: # ball_and_stick
                        scale = settings.get('ball_stick_atom_scale', 1.0)
                        try:
                            # Use 0.3 factor matching app's constants.py VDW_RADII
                            r = pt.GetRvdw(pt.GetAtomicNumber(symbol))
                            radius = (r if r > 0.1 else 1.5) * 0.3 * scale
                        except:
                            radius = 0.4 * scale
                atoms_to_draw.append((np.array([pos.x, pos.y, pos.z]), radius, color, symbol, idx))

        for pos, radius, color, symbol, idx in atoms_to_draw:
            lines.append(f'// Atom {idx}: {symbol}')
            lines.append('sphere {')
            lines.append(f'  <{pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f}>, {radius:.4f}')
            lines.append('  texture {')
            lines.append(f'    pigment {{ color rgb <{color[0]:.3f}, {color[1]:.3f}, {color[2]:.3f}> }}')
            lines.append('    finish { AtomFinish }')
            lines.append('  }')
            lines.append('}')
            lines.append('')
    # Always populate atom_positions for bonds
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        atom_positions[idx] = np.array([pos.x, pos.y, pos.z])
    
    # Create bonds
    if current_style != 'cpk':
        lines.append('// ===== BONDS =====')
        lines.append('')
        
        bond_counter = 0
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
            
            pos1 = atom_positions[begin_idx]
            pos2 = atom_positions[end_idx]
            d = pos2 - pos1
            h = np.linalg.norm(d)
            if h == 0: continue
            
            # Helper to add a cylinder segment to POV-Ray lines
            def add_pov_segment(p1, p2, radius, rgb_255, comment):
                rgb = tuple(c / 255.0 for c in rgb_255)
                lines.append(f'// {comment}')
                lines.append('cylinder {')
                lines.append(f'  <{p1[0]:.4f}, {p1[1]:.4f}, {p1[2]:.4f}>,')
                lines.append(f'  <{p2[0]:.4f}, {p2[1]:.4f}, {p2[2]:.4f}>,')
                lines.append(f'  {radius:.4f}')
                lines.append('  texture {')
                lines.append(f'    pigment {{ color rgb <{rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f}> }}')
                lines.append('    finish { BondFinish }')
                lines.append('  }')
                lines.append('}')
    
            if bond_type == Chem.rdchem.BondType.SINGLE or bond_type == Chem.rdchem.BondType.AROMATIC:
                if f'bond_{bond_counter}' in color_map:
                    # Uniform color
                    add_pov_segment(pos1, pos2, cyl_radius, color_map[f'bond_{bond_counter}'], f'Bond {bond_counter}')
                elif f'bond_{bond_counter}_start' in color_map:
                    # CPK Split
                    mid = (pos1 + pos2) / 2
                    add_pov_segment(pos1, mid, cyl_radius, color_map[f'bond_{bond_counter}_start'], f'Bond {bond_counter} start')
                    add_pov_segment(mid, pos2, cyl_radius, color_map[f'bond_{bond_counter}_end'], f'Bond {bond_counter} end')
                else:
                    # Fallback to gray
                    add_pov_segment(pos1, pos2, cyl_radius, (153, 153, 153), f'Bond {bond_counter} (fallback)')
            else:
                if bond_type == Chem.rdchem.BondType.DOUBLE:
                    r = cyl_radius * double_radius_factor
                    off_dir = _calculate_double_bond_offset(mol, bond, conf)
                    s_double = cyl_radius * double_offset_factor
                    p1_start, p1_end = pos1 + off_dir * (s_double / 2), pos2 + off_dir * (s_double / 2)
                    p2_start, p2_end = pos1 - off_dir * (s_double / 2), pos2 - off_dir * (s_double / 2)
                    
                    # Colors for first segment
                    if f'bond_{bond_counter}_1' in color_map:
                        add_pov_segment(p1_start, p1_end, r, color_map[f'bond_{bond_counter}_1'], f'Bond {bond_counter} segment 1')
                    else:
                        mid = (p1_start + p1_end) / 2
                        add_pov_segment(p1_start, mid, r, color_map.get(f'bond_{bond_counter}_1_start', (153,153,153)), f'Bond {bond_counter} seg 1 start')
                        add_pov_segment(mid, p1_end, r, color_map.get(f'bond_{bond_counter}_1_end', (153,153,153)), f'Bond {bond_counter} seg 1 end')
                    
                    # Colors for second segment
                    if f'bond_{bond_counter}_2' in color_map:
                        add_pov_segment(p2_start, p2_end, r, color_map[f'bond_{bond_counter}_2'], f'Bond {bond_counter} segment 2')
                    else:
                        mid = (p2_start + p2_end) / 2
                        add_pov_segment(p2_start, mid, r, color_map.get(f'bond_{bond_counter}_2_start', (153,153,153)), f'Bond {bond_counter} seg 2 start')
                        add_pov_segment(mid, p2_end, r, color_map.get(f'bond_{bond_counter}_2_end', (153,153,153)), f'Bond {bond_counter} seg 2 end')
                        
                elif bond_type == Chem.rdchem.BondType.TRIPLE:
                    r = cyl_radius * triple_radius_factor
                    s_triple = cyl_radius * triple_offset_factor
                    
                    if f'bond_{bond_counter}_1' in color_map:
                        add_pov_segment(pos1, pos2, r, color_map[f'bond_{bond_counter}_1'], f'Bond {bond_counter} center')
                    else:
                        mid = (pos1 + pos2) / 2
                        ccol1 = color_map.get(f'bond_{bond_counter}_1_start', (153, 153, 153))
                        ccol2 = color_map.get(f'bond_{bond_counter}_1_end', (153, 153, 153))
                        add_pov_segment(pos1, mid, r, ccol1, f'Bond {bond_counter} center start')
                        add_pov_segment(mid, pos2, r, ccol2, f'Bond {bond_counter} center end')
                    
                    # Sides
                    for i, sign in enumerate([1, -1], 2):
                        offset = off_dir * s_triple * sign
                        p_start, p_end = pos1 + offset, pos2 + offset
                        suffix = f'_{i}'
                        if f'bond_{bond_counter}{suffix}' in color_map:
                            add_pov_segment(p_start, p_end, r, color_map[f'bond_{bond_counter}{suffix}'], f'Bond {bond_counter} side {i-1}')
                        else:
                            mid = (p_start + p_end) / 2
                            add_pov_segment(p_start, mid, r, color_map.get(f'bond_{bond_counter}{suffix}_start', (153,153,153)), f'Bond {bond_counter} side {i-1} start')
                            add_pov_segment(mid, p_end, r, color_map.get(f'bond_{bond_counter}{suffix}_end', (153,153,153)), f'Bond {bond_counter} side {i-1} end')

            bond_counter += 1
        lines.append('')
    
    # Optional reflecting floor (commented out by default to avoid "molecular shadow")
    lines.append('// Optional reflecting floor (uncomment to enable)')
    lines.append('/*')
    lines.append('plane {')
    lines.append(f'  y, {center_y - max_dist - 2}')
    lines.append('  texture {')
    lines.append(f'    pigment {{ color rgb <{bg_rgb[0]*0.8:.3f}, {bg_rgb[1]*0.8:.3f}, {bg_rgb[2]*0.8:.3f}> }}')
    lines.append('    finish { ')
    lines.append('      ambient 0.1 diffuse 0.6 ')
    lines.append('      reflection { 0.1, 0.3 fresnel on } ')
    lines.append('      conserve_energy')
    lines.append('    }')
    lines.append('  }')
    lines.append('}')
    lines.append('*/')
    lines.append('')
    
    lines.append(f'// End of scene - {mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds')
    
    return '\n'.join(lines), exp_w, exp_h
