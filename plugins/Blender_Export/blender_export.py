#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Blender Export Plugin for MoleditPy
Exports molecular structures as Blender Python scripts
"""

PLUGIN_NAME = "Blender Export"
PLUGIN_VERSION = "2026.02.05"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Export molecular structures as Blender Python scripts that create 3D visualizations"


def initialize(context):
    """Initialize the Blender export plugin."""
    context.add_export_action("Export to Blender Script...", lambda: export_to_blender(context))


def export_to_blender(context):
    """Main export function that generates a Blender Python script."""
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
    default_name = "molecule_blender"
    try:
        if hasattr(mw, 'current_file_path') and mw.current_file_path:
            default_dir = os.path.dirname(mw.current_file_path)
            base_name = os.path.splitext(os.path.basename(mw.current_file_path))[0]
            default_name = f"{base_name}_blender"
    except Exception:
        pass
    
    default_path = os.path.join(default_dir, default_name) if default_dir else default_name
    
    # Ask user where to save the file
    file_path, _ = QFileDialog.getSaveFileName(
        mw, 
        "Export to Blender Script", 
        default_path, 
        "Python Scripts (*.py);;All Files (*)"
    )
    
    if not file_path:
        return
    
    # Ensure .py extension
    if not file_path.lower().endswith('.py'):
        file_path += '.py'
    
    try:
        # Generate and save the Blender script
        script_content = generate_blender_script(mol, mw)
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(script_content)
        
        mw.statusBar().showMessage(f"Blender script exported to {file_path}")
        
        # Show success message with instructions
        QMessageBox.information(
            mw, 
            "Export Successful",
            f"Blender script saved to:\n{file_path}\n\n"
            "To use in Blender:\n"
            "1. Open Blender\n"
            "2. Switch to the Scripting workspace\n"
            "3. Open the exported .py file\n"
            "4. Click 'Run Script' or press Alt+P"
        )
        
    except Exception as e:
        QMessageBox.critical(mw, "Export Failed", f"Failed to export Blender script:\n{str(e)}")
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


def generate_blender_script(mol, mw):
    """Generate the Blender Python script content."""
    from rdkit import Chem
    import math
    import numpy as np
    
    # Get 3D coordinates
    conf = mol.GetConformer()
    
    # Extract actual rendering data from the 3D view
    mw = mw # Rename for clarity if needed
    settings = getattr(mw, 'settings', {})
    
    # Kekulize if enabled in settings to match 3D view indices and bond types
    mol_to_draw = mol
    if settings.get('display_kekule_3d', False):
        try:
            mol_to_draw = Chem.Mol(mol)
            Chem.Kekulize(mol_to_draw, clearAromaticFlags=True)
        except Exception:
            mol_to_draw = mol
            
    # Get atom and bond data
    atoms_data = []
    bonds_data = []
    
    color_map = getattr(mw, '_3d_color_map', {})
    raw_style = getattr(mw, 'current_3d_style', 'ball_and_stick')
    current_style = str(raw_style).lower().replace(' ', '_').replace('&', '_and_')
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
        bond_radius = settings.get('wireframe_bond_radius', 0.01)
        double_radius_factor = settings.get('wireframe_double_bond_radius_factor', 0.8)
        triple_radius_factor = settings.get('wireframe_triple_bond_radius_factor', 0.75)
        double_offset_factor = settings.get('wireframe_double_bond_offset_factor', 3.0)
        triple_offset_factor = settings.get('wireframe_triple_bond_offset_factor', 3.0)
    elif current_style == 'stick':
        bond_radius = settings.get('stick_bond_radius', 0.15)
        double_radius_factor = settings.get('stick_double_bond_radius_factor', 0.60)
        triple_radius_factor = settings.get('stick_triple_bond_radius_factor', 0.40)
        double_offset_factor = settings.get('stick_double_bond_offset_factor', 1.5)
        triple_offset_factor = settings.get('stick_triple_bond_offset_factor', 1.0)
    else:  # ball_and_stick
        bond_radius = settings.get('ball_stick_bond_radius', 0.1)
        double_radius_factor = settings.get('ball_stick_double_bond_radius_factor', 0.8)
        triple_radius_factor = settings.get('ball_stick_triple_bond_radius_factor', 0.75)
        double_offset_factor = settings.get('ball_stick_double_bond_offset_factor', 2.0)
        triple_offset_factor = settings.get('ball_stick_triple_bond_offset_factor', 2.0)

    # Periodic Table for accurate radii
    pt = Chem.GetPeriodicTable()
    
    # Atom radius calculation logic (matches main_window_view_3d.py)
    def determine_atom_radius(atom_idx, symbol):
        if current_style == 'cpk':
            scale = settings.get('cpk_atom_scale', 1.0)
            try:
                r = pt.GetRvdw(pt.GetAtomicNumber(symbol))
                return (r if r > 0.1 else 1.5) * scale
            except:
                return 1.5 * scale
        elif current_style == 'stick':
            return settings.get('stick_bond_radius', 0.15)
        elif current_style == 'wireframe':
            return 0.01 # Matches app setting
        else: # ball_and_stick
            scale = settings.get('ball_stick_atom_scale', 1.0)
            try:
                # App VDW_RADII is already scaled by 0.3 in constants.py
                r = pt.GetRvdw(pt.GetAtomicNumber(symbol))
                return (r if r > 0.1 else 1.5) * 0.3 * scale
            except:
                return 0.4 * scale # App default fallback

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
    
    # Extract specular/roughness from settings
    specular_val = settings.get('specular', 0.2)
    specular_power = settings.get('specular_power', 20)
    # Convert specular power to roughness (approximate)
    roughness_val = max(0.05, min(1.0, 1.0 - (specular_power / 100.0)))
    
    # Background color from settings
    bg_rgb = (0.31, 0.31, 0.31) # Default gray
    
    # Calculate molecule center and bounds for camera positioning (fallback)
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
    
    # Get current camera from plotter if available
    camera_data = {
        'location': (float(center_x), float(center_y), float(center_z + max_dist * 3)),
        'look_at': (float(center_x), float(center_y), float(center_z)),
        'sky': (0.0, 1.0, 0.0),
        'angle': 45.0,
        'orthographic': False
    }
    
    if hasattr(mw, 'plotter') and mw.plotter is not None:
        if hasattr(mw.plotter, 'camera'):
            vtk_camera = mw.plotter.camera
            pos = vtk_camera.GetPosition()
            focal = vtk_camera.GetFocalPoint()
            up = vtk_camera.GetViewUp()
            
            camera_data['location'] = tuple(float(c) for c in pos)
            camera_data['look_at'] = tuple(float(c) for c in focal)
            camera_data['sky'] = tuple(float(c) for c in up)
            
            if vtk_camera.GetParallelProjection():
                camera_data['orthographic'] = True
                camera_data['parallel_scale'] = float(vtk_camera.GetParallelScale())
                
            camera_data['angle'] = float(vtk_camera.GetViewAngle())
    
    # Extract atom data
    if current_style != 'wireframe':
        for atom in mol_to_draw.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            symbol = atom.GetSymbol()
            
            # Get color from actual 3D view if available, otherwise use CPK
            color_key = f'atom_{idx}'
            if color_key in color_map:
                # Color from 3D view (RGB 0-255)
                rgb_255 = color_map[color_key]
                color = tuple(c / 255.0 for c in rgb_255)
            else:
                # Fallback to CPK color
                color = cpk_colors.get(symbol, (144/255.0, 144/255.0, 144/255.0))
            
            radius = determine_atom_radius(idx, symbol)
            
            atoms_data.append({
                'idx': idx,
                'symbol': symbol,
                'pos': (float(pos.x), float(pos.y), float(pos.z)),
                'color': tuple(float(c) for c in color),
                'radius': float(radius),
                'degree': atom.GetDegree()
            })
    
    # Extract bond data
    if current_style != 'cpk':
        bond_counter = 0
        for bond in mol_to_draw.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
            
            # Get positions
            pos1 = conf.GetAtomPosition(begin_idx)
            pos2 = conf.GetAtomPosition(end_idx)
            p1 = np.array([pos1.x, pos1.y, pos1.z])
            p2 = np.array([pos2.x, pos2.y, pos2.z])
            d = p2 - p1
            h = np.linalg.norm(d)
            if h == 0: continue
            
            # Multi-bond variables
            v1 = d / h
            if bond_type == Chem.rdchem.BondType.DOUBLE:
                off_dir = _calculate_double_bond_offset(mol_to_draw, bond, conf)
            else:
                v_arb = np.array([0.0, 0.0, 1.0])
                if np.allclose(np.abs(np.dot(v1, v_arb)), 1.0): v_arb = np.array([0.0, 1.0, 0.0])
                off_dir = np.cross(v1, v_arb)
                off_dir /= np.linalg.norm(off_dir)
    
            if bond_type == Chem.rdchem.BondType.SINGLE or bond_type == Chem.rdchem.BondType.AROMATIC:
                if f'bond_{bond_counter}' in color_map:
                    # Uniform color
                    color = tuple(float(c) / 255.0 for c in color_map[f'bond_{bond_counter}'])
                    bonds_data.append({'pos1': tuple(float(c) for c in p1), 'pos2': tuple(float(c) for c in p2), 'radius': float(bond_radius), 'color': color})
                elif f'bond_{bond_counter}_start' in color_map:
                    # CPK Split
                    mid = (p1 + p2) / 2
                    c_start = tuple(float(c) / 255.0 for c in color_map[f'bond_{bond_counter}_start'])
                    c_end = tuple(float(c) / 255.0 for c in color_map[f'bond_{bond_counter}_end'])
                    bonds_data.append({'pos1': tuple(float(c) for c in p1), 'pos2': tuple(float(c) for c in mid), 'radius': float(bond_radius), 'color': c_start})
                    bonds_data.append({'pos1': tuple(float(c) for c in mid), 'pos2': tuple(float(c) for c in p2), 'radius': float(bond_radius), 'color': c_end})
                else:
                    bonds_data.append({'pos1': tuple(float(c) for c in p1), 'pos2': tuple(float(c) for c in p2), 'radius': float(bond_radius), 'color': (0.6, 0.6, 0.6)})
            else:
                if bond_type == Chem.rdchem.BondType.DOUBLE:
                    r = float(bond_radius * double_radius_factor)
                    s_double = float(bond_radius * double_offset_factor)
                    s1_start, s1_end = p1 + off_dir * (s_double / 2), p2 + off_dir * (s_double / 2)
                    s2_start, s2_end = p1 - off_dir * (s_double / 2), p2 - off_dir * (s_double / 2)
                    
                    # Segment 1
                    if f'bond_{bond_counter}_1' in color_map:
                        color = tuple(float(c) / 255.0 for c in color_map[f'bond_{bond_counter}_1'])
                        bonds_data.append({'pos1': tuple(float(c) for c in s1_start), 'pos2': tuple(float(c) for c in s1_end), 'radius': r, 'color': color})
                    else:
                        mid = (s1_start + s1_end) / 2
                        c_start = tuple(float(c) / 255.0 for c in color_map.get(f'bond_{bond_counter}_1_start', (153,153,153)))
                        c_end = tuple(float(c) / 255.0 for c in color_map.get(f'bond_{bond_counter}_1_end', (153,153,153)))
                        bonds_data.append({'pos1': tuple(float(c) for c in s1_start), 'pos2': tuple(float(c) for c in mid), 'radius': r, 'color': c_start})
                        bonds_data.append({'pos1': tuple(float(c) for c in mid), 'pos2': tuple(float(c) for c in s1_end), 'radius': r, 'color': c_end})
                    # Segment 2
                    if f'bond_{bond_counter}_2' in color_map:
                        color = tuple(float(c) / 255.0 for c in color_map[f'bond_{bond_counter}_2'])
                        bonds_data.append({'pos1': tuple(float(c) for c in s2_start), 'pos2': tuple(float(c) for c in s2_end), 'radius': r, 'color': color})
                    else:
                        mid = (s2_start + s2_end) / 2
                        c_start = tuple(float(c) / 255.0 for c in color_map.get(f'bond_{bond_counter}_2_start', (153,153,153)))
                        c_end = tuple(float(c) / 255.0 for c in color_map.get(f'bond_{bond_counter}_2_end', (153,153,153)))
                        bonds_data.append({'pos1': tuple(float(c) for c in s2_start), 'pos2': tuple(float(c) for c in mid), 'radius': r, 'color': c_start})
                        bonds_data.append({'pos1': tuple(float(c) for c in mid), 'pos2': tuple(float(c) for c in s2_end), 'radius': r, 'color': c_end})
                
                elif bond_type == Chem.rdchem.BondType.TRIPLE:
                    r = float(bond_radius * triple_radius_factor)
                    s_triple = float(bond_radius * triple_offset_factor)
                    # Center
                    if f'bond_{bond_counter}_1' in color_map:
                        color = tuple(float(c) / 255.0 for c in color_map[f'bond_{bond_counter}_1'])
                        bonds_data.append({'pos1': tuple(float(c) for c in p1), 'pos2': tuple(float(c) for c in p2), 'radius': r, 'color': color})
                    else:
                        mid = (p1 + p2) / 2
                        c_start = tuple(float(c) / 255.0 for c in color_map.get(f'bond_{bond_counter}_1_start', (153,153,153)))
                        c_end = tuple(float(c) / 255.0 for c in color_map.get(f'bond_{bond_counter}_1_end', (153,153,153)))
                        bonds_data.append({'pos1': tuple(float(c) for c in p1), 'pos2': tuple(float(c) for c in mid), 'radius': r, 'color': c_start})
                        bonds_data.append({'pos1': tuple(float(c) for c in mid), 'pos2': tuple(float(c) for c in p2), 'radius': r, 'color': c_end})
                    # Sides
                    for i, sign in enumerate([1, -1], 2):
                        offset = off_dir * s_triple * sign
                        s_start, s_end = p1 + offset, p2 + offset
                        suffix = f'_{i}'
                        if f'bond_{bond_counter}{suffix}' in color_map:
                            color = tuple(float(c) / 255.0 for c in color_map[f'bond_{bond_counter}{suffix}'])
                            bonds_data.append({'pos1': tuple(float(c) for c in s_start), 'pos2': tuple(float(c) for c in s_end), 'radius': r, 'color': color})
                        else:
                            mid = (s_start + s_end) / 2
                            c_start = tuple(float(c) / 255.0 for c in color_map.get(f'bond_{bond_counter}{suffix}_start', (153,153,153)))
                            c_end = tuple(float(c) / 255.0 for c in color_map.get(f'bond_{bond_counter}{suffix}_end', (153,153,153)))
                            bonds_data.append({'pos1': tuple(float(c) for c in s_start), 'pos2': tuple(float(c) for c in mid), 'radius': r, 'color': c_start})
                            bonds_data.append({'pos1': tuple(float(c) for c in mid), 'pos2': tuple(float(c) for c in s_end), 'radius': r, 'color': c_end})
            
            bond_counter += 1

    
    # Generate the script
    script_lines = []
    
    # Header
    script_lines.append('"""')
    script_lines.append('Molecular Structure for Blender')
    script_lines.append(f'Generated from MoleditPy - {PLUGIN_NAME} v{PLUGIN_VERSION}')
    script_lines.append('')
    script_lines.append('Instructions:')
    script_lines.append('1. Open this file in Blender\'s Text Editor')
    script_lines.append('2. Press Alt+P or click "Run Script"')
    script_lines.append('3. The molecular structure will be created in the 3D viewport')
    script_lines.append('"""')
    script_lines.append('')
    
    # Imports
    script_lines.append('import bpy')
    script_lines.append('import math')
    script_lines.append('from mathutils import Vector')
    # Settings at the top
    script_lines.append('# Configuration Settings')
    script_lines.append('CLEAR_SCENE = True  # Set to True to clear the entire scene before creating the molecule')
    script_lines.append('')
    
    # Clear existing objects function
    script_lines.append('def clear_scene():')
    script_lines.append('    """Clear existing objects if enabled"""')
    script_lines.append('    if CLEAR_SCENE:')
    script_lines.append('        print("Clearing scene...")')
    script_lines.append('        bpy.ops.object.select_all(action=\'SELECT\')')
    script_lines.append('        bpy.ops.object.delete()')
    script_lines.append('')
    script_lines.append('')
    
    # Create material function
    script_lines.append('def create_material(name, color):')
    script_lines.append('    """Create a material with the given color"""')
    script_lines.append('    if name in bpy.data.materials:')
    script_lines.append('        return bpy.data.materials[name]')
    script_lines.append('        ')
    script_lines.append('    mat = bpy.data.materials.new(name=name)')
    script_lines.append('    mat.use_nodes = True')
    script_lines.append('    mat.diffuse_color = (*color, 1.0)')
    script_lines.append('    nodes = mat.node_tree.nodes')
    script_lines.append('    nodes.clear()')
    script_lines.append('    ')
    script_lines.append('    # Create Principled BSDF node')
    script_lines.append('    bsdf = nodes.new(type=\'ShaderNodeBsdfPrincipled\')')
    script_lines.append('    try:')
    script_lines.append('        # Premium Metallic Look')
    script_lines.append('        bsdf.inputs[\'Base Color\'].default_value = (*color, 1.0)')
    script_lines.append('        bsdf.inputs[\'Metallic\'].default_value = 0.9')
    script_lines.append('        bsdf.inputs[\'Roughness\'].default_value = 0.2')
    script_lines.append('    except:')
    script_lines.append('        # Fallback for different Blender versions')
    script_lines.append('        bsdf.inputs[0].default_value = (*color, 1.0)')
    script_lines.append('        if len(bsdf.inputs) > 6:')
    script_lines.append('            bsdf.inputs[6].default_value = 0.9 # Metallic')
    script_lines.append('            bsdf.inputs[9].default_value = 0.2 # Roughness')
    script_lines.append('    ')
    script_lines.append('    # Create output node')
    script_lines.append('    output = nodes.new(type=\'ShaderNodeOutputMaterial\')')
    script_lines.append('    ')
    script_lines.append('    # Link nodes')
    script_lines.append('    mat.node_tree.links.new(bsdf.outputs[\'BSDF\'], output.inputs[\'Surface\'])')
    script_lines.append('    ')
    script_lines.append('    return mat')
    script_lines.append('')
    
    # Create point light function
    script_lines.append('def create_light(location, energy=2.0, name="Sun_Light"):')
    script_lines.append('    """Create a sun light at the given location pointing at the scene center"""')
    script_lines.append('    light_data = bpy.data.lights.new(name=name, type=\'SUN\')')
    script_lines.append('    light_data.energy = energy')
    script_lines.append('    try:')
    script_lines.append('        light_data.use_contact_shadow = True')
    script_lines.append('    except:')
    script_lines.append('        pass')
    script_lines.append('    light_object = bpy.data.objects.new(name=name, object_data=light_data)')
    script_lines.append('    bpy.context.collection.objects.link(light_object)')
    script_lines.append('    light_object.location = location')
    script_lines.append('    ')
    script_lines.append('    # Point sun direction from camera to origin')
    script_lines.append('    from mathutils import Vector')
    script_lines.append('    dir_vec = -Vector(location).normalized()')
    script_lines.append('    light_object.rotation_euler = dir_vec.to_track_quat(\'-Z\', \'Y\').to_euler()')
    script_lines.append('    return light_object')
    script_lines.append('')
    script_lines.append('def create_rim_light(location, energy=100, name="Rim_Light"):')
    script_lines.append('    """Create a rim light to make objects pop"""')
    script_lines.append('    light_data = bpy.data.lights.new(name=name, type=\'POINT\')')
    script_lines.append('    light_data.energy = energy')
    script_lines.append('    light_data.color = (0.8, 0.9, 1.0) # Slightly cool rim')
    script_lines.append('    light_object = bpy.data.objects.new(name=name, object_data=light_data)')
    script_lines.append('    bpy.context.collection.objects.link(light_object)')
    script_lines.append('    light_object.location = location')
    script_lines.append('    return light_object')
    script_lines.append('')
    script_lines.append('')
    
    # Create atom function
    script_lines.append('def create_atom(location, radius, color, name):')
    script_lines.append('    """Create a sphere representing an atom"""')
    script_lines.append('    bpy.ops.mesh.primitive_uv_sphere_add(')
    script_lines.append('        radius=radius,')
    script_lines.append('        location=location,')
    script_lines.append('        segments=64,')
    script_lines.append('        ring_count=32')
    script_lines.append('    )')
    script_lines.append('    obj = bpy.context.active_object')
    script_lines.append('    obj.name = name')
    script_lines.append('    ')
    script_lines.append('    # Create and assign material')
    script_lines.append('    mat = create_material(f"Mat_Atom_{name}", color)')
    script_lines.append('    obj.data.materials.append(mat)')
    script_lines.append('    ')
    script_lines.append('    # Enable smooth shading')
    script_lines.append('    bpy.ops.object.shade_smooth()')
    script_lines.append('    # For newer Blender versions, ensure auto-smooth is off for spheres')
    script_lines.append('    # but for atoms we generally want pure smooth.')
    script_lines.append('    # If the user says \'not smooth\', it might be because of low resolution or auto-smooth limits.')
    script_lines.append('    if hasattr(obj.data, \'use_auto_smooth\'):')
    script_lines.append('        obj.data.use_auto_smooth = False')
    script_lines.append('    ')
    script_lines.append('    return obj')
    script_lines.append('')
    
    # Create camera function
    script_lines.append('def create_camera(location, look_at, sky, angle=45.0, orthographic=False, parallel_scale=5.0):')
    script_lines.append('    """Create and position the camera to match MoleditPy view exactly"""')
    script_lines.append('    bpy.ops.object.camera_add(location=location)')
    script_lines.append('    cam = bpy.context.active_object')
    script_lines.append('    ')
    script_lines.append('    if orthographic:')
    script_lines.append('        cam.data.type = \'ORTHO\'')
    script_lines.append('        cam.data.ortho_scale = parallel_scale * 2.0')
    script_lines.append('    else:')
    script_lines.append('        cam.data.type = \'PERSP\'')
    script_lines.append('        cam.data.lens_unit = \'FOV\'')
    script_lines.append('        # VTK angle is vertical FOV')
    script_lines.append('        cam.data.sensor_fit = \'VERTICAL\'')
    script_lines.append('        cam.data.angle = math.radians(angle)')
    script_lines.append('    ')
    script_lines.append('    cam.name = "MoleditPy_Camera"')
    script_lines.append('    ')
    script_lines.append('    # Ensure point of view matches target exactly')
    script_lines.append('    from mathutils import Matrix, Vector')
    script_lines.append('    cam_pos = Vector(location)')
    script_lines.append('    target = Vector(look_at)')
    script_lines.append('    up_vec = Vector(sky).normalized()')
    script_lines.append('    ')
    script_lines.append('    # Forward vector (what the camera is looking at)')
    script_lines.append('    forward = (target - cam_pos).normalized()')
    script_lines.append('    ')
    script_lines.append('    # In Blender, local Z+ is BACKWARDS. So local Z- is forward.')
    script_lines.append('    z_axis = -forward')
    script_lines.append('    # Local X is Right')
    script_lines.append('    x_axis = up_vec.cross(z_axis).normalized()')
    script_lines.append('    # Local Y is Up')
    script_lines.append('    y_axis = z_axis.cross(x_axis).normalized()')
    script_lines.append('    ')
    script_lines.append('    # Construct rotation matrix')
    script_lines.append('    rot_mat = Matrix((x_axis, y_axis, z_axis)).transposed()')
    script_lines.append('    cam.rotation_euler = rot_mat.to_euler()')
    script_lines.append('    ')
    script_lines.append('    bpy.context.scene.camera = cam')
    script_lines.append('    return cam')
    script_lines.append('')
    
    # Create bond function
    script_lines.append('def create_bond(pos1, pos2, radius=0.1, color=(0.5, 0.5, 0.5), name="Bond"):')
    script_lines.append('    """Create a cylinder representing a bond"""')
    script_lines.append('    # Calculate bond vector and properties')
    script_lines.append('    v1 = Vector(pos1)')
    script_lines.append('    v2 = Vector(pos2)')
    script_lines.append('    direction = v2 - v1')
    script_lines.append('    length = direction.length')
    script_lines.append('    center = (v1 + v2) / 2')
    script_lines.append('    ')
    script_lines.append('    # Create cylinder')
    script_lines.append('    bpy.ops.mesh.primitive_cylinder_add(')
    script_lines.append('        radius=radius,')
    script_lines.append('        depth=length,')
    script_lines.append('        location=center,')
    script_lines.append('        vertices=64')
    script_lines.append('    )')
    script_lines.append('    obj = bpy.context.active_object')
    script_lines.append('    obj.name = name')
    script_lines.append('    ')
    script_lines.append('    # Align cylinder with bond direction')
    script_lines.append('    # Use rotation_difference for robust orientation')
    script_lines.append('    z_axis = Vector((0, 0, 1))')
    script_lines.append('    direction.normalize()')
    script_lines.append('    quat = z_axis.rotation_difference(direction)')
    script_lines.append('    obj.rotation_mode = \'QUATERNION\'')
    script_lines.append('    obj.rotation_quaternion = quat')
    script_lines.append('    ')
    script_lines.append('    # Create and assign material')
    script_lines.append('    mat = create_material(f"Mat_Bond_{name}", color)')
    script_lines.append('    obj.data.materials.append(mat)')
    script_lines.append('    ')
    script_lines.append('    # Enable smooth shading by angle (prevents blobby caps)')
    script_lines.append('    try:')
    script_lines.append('        # Newer Blender API (4.1+)')
    script_lines.append('        bpy.ops.object.shade_smooth_by_angle(angle=math.radians(30))')
    script_lines.append('    except:')
    script_lines.append('        try:')
    script_lines.append('            # Older Blender API')
    script_lines.append('            bpy.ops.object.shade_smooth()')
    script_lines.append('            obj.data.use_auto_smooth = True')
    script_lines.append('            obj.data.auto_smooth_angle = math.radians(30)')
    script_lines.append('        except:')
    script_lines.append('            bpy.ops.object.shade_smooth()')
    script_lines.append('    ')
    script_lines.append('    return obj')
    script_lines.append('')
    
    # Main function
    script_lines.append('def create_molecule():')
    script_lines.append('    """Create the complete molecular structure"""')
    script_lines.append('    print("Creating molecular structure...")')
    script_lines.append('    ')
    script_lines.append('    # Clear existing objects')
    script_lines.append('    clear_scene()')
    script_lines.append('    ')
    script_lines.append('    # Create atoms')
    script_lines.append('    atoms = []')
    script_lines.append('    # Using a cache for materials to avoid creating thousands of them')
    script_lines.append('    material_cache = {}')
    script_lines.append('    ')
    script_lines.append('    def get_cached_material(color_tuple, prefix="Mat"):')
    script_lines.append('        color_key = tuple(round(c, 3) for c in color_tuple)')
    script_lines.append('        mat_key = (prefix, color_key)')
    script_lines.append('        if mat_key not in material_cache:')
    script_lines.append(f'            mat_name = f"{{prefix}}_{{color_key[0]}}_{{color_key[1]}}_{{color_key[2]}}"')
    script_lines.append('            material_cache[mat_key] = create_material(mat_name, color_tuple)')
    script_lines.append('        return material_cache[mat_key]')
    script_lines.append('')
    
    # Add cinematic lighting setup
    script_lines.append('    # Add high-end 3-point lighting setup')
    lp = np.array(camera_data['location'])
    la = np.array(camera_data['look_at'])
    sky = np.array(camera_data['sky'])
    
    fwd = la - lp
    fwd_norm = np.linalg.norm(fwd)
    if fwd_norm > 1e-6:
        fwd_u = fwd / fwd_norm
        sky_u = sky / np.linalg.norm(sky) if np.linalg.norm(sky) > 1e-6 else np.array([0, 1, 0])
        right_u = np.cross(fwd_u, sky_u)
        rn = np.linalg.norm(right_u)
        right_u = right_u / rn if rn > 1e-6 else np.array([1, 0, 0])
        
        # Main Sun (Key)
        lp_key = lp + (right_u + sky_u) * (fwd_norm * 0.2)
        script_lines.append(f'    create_light(({lp_key[0]:.3f}, {lp_key[1]:.3f}, {lp_key[2]:.3f}), energy=2.5, name="Key_Light")')
        
        # Rim Light (Back)
        lp_rim = la - fwd_u * (fwd_norm * 0.5) + (sky_u - right_u) * (fwd_norm * 0.4)
        script_lines.append(f'    create_rim_light(({lp_rim[0]:.3f}, {lp_rim[1]:.3f}, {lp_rim[2]:.3f}), energy={float(100 * max_dist**2):.1f})')
    else:
        script_lines.append(f'    create_light({camera_data["location"]}, energy=2.0)')
    script_lines.append('')
    
    # Enable Cinematic Viewport Effects
    script_lines.append('    # Enable Cinematic Viewport FX')
    script_lines.append('    scene = bpy.context.scene')
    script_lines.append('    if hasattr(scene, "eevee"):')
    script_lines.append('        try: scene.eevee.use_bloom = True')
    script_lines.append('        except: pass')
    script_lines.append('        try: scene.eevee.use_gtao = True')
    script_lines.append('        except: ')
    script_lines.append('            try: scene.eevee.use_ambient_occlusion = True')
    script_lines.append('            except: pass')
    script_lines.append('        try: scene.eevee.use_ssr = True')
    script_lines.append('        except: ')
    script_lines.append('            try: scene.eevee.use_raytracing = True')
    script_lines.append('            except: pass')
    script_lines.append('        try: scene.eevee.bloom_intensity = 0.05')
    script_lines.append('        except: pass')
    script_lines.append('')
    
    # Add atom creation code
    for atom in atoms_data:
        color = atom["color"]
        pos = atom["pos"]
        radius = atom["radius"]
        idx = atom["idx"]
        symbol = atom["symbol"]
        
        # Stick mode terminal split
        if current_style == 'stick' and atom["degree"] == 1:
            rd_atom = mol_to_draw.GetAtomWithIdx(idx)
            bond = rd_atom.GetBonds()[0]
            bt = bond.GetBondType()
            if bt in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
                other_idx = bond.GetBeginAtomIdx() if bond.GetEndAtomIdx() == idx else bond.GetEndAtomIdx()
                pos_i = np.array(pos)
                pos_other = np.array(conf.GetAtomPosition(other_idx))
                bv = pos_i - pos_other
                bu = bv / np.linalg.norm(bv)
                
                if bt == Chem.BondType.DOUBLE:
                    od = _calculate_double_bond_offset(mol_to_draw, bond, conf)
                    odist = bond_radius * double_offset_factor / 2
                    srad = bond_radius * double_radius_factor
                    offsets = [od * odist, -od * odist]
                else: # TRIPLE
                    v_arb = np.array([0, 0, 1])
                    if np.allclose(np.abs(np.dot(bu, v_arb)), 1.0): v_arb = np.array([0, 1, 0])
                    od = np.cross(bu, v_arb)
                    od /= np.linalg.norm(od)
                    odist = bond_radius * triple_offset_factor
                    srad = bond_radius * triple_radius_factor
                    offsets = [np.zeros(3), od * odist, -od * odist]
                
                for j, off in enumerate(offsets):
                    p_off = tuple(float(val) for val in (pos_i + off))
                    c_val = tuple(float(val) for val in color)
                    script_lines.append(f'    obj = create_atom(')
                    script_lines.append(f'        location={p_off},')
                    script_lines.append(f'        radius={srad:.3f},')
                    script_lines.append(f'        color={c_val},')
                    script_lines.append(f'        name="Atom_{idx}_{symbol}_{j}"')
                    script_lines.append(f'    )')
                    script_lines.append(f'    mat = get_cached_material({c_val}, "AtomMat")')
                    script_lines.append(f'    if obj.data.materials: obj.data.materials[0] = mat')
                    script_lines.append(f'    else: obj.data.materials.append(mat)')
                    script_lines.append(f'    atoms.append(obj)')
                continue

        p_val = tuple(float(c) for c in pos)
        c_val = tuple(float(c) for c in color)
        radius_val = float(radius)
        script_lines.append(f'    obj = create_atom(')
        script_lines.append(f'        location={p_val},')
        script_lines.append(f'        radius={radius_val:.3f},')
        script_lines.append(f'        color={c_val},')
        script_lines.append(f'        name="Atom_{idx}_{symbol}"')
        script_lines.append(f'    )')
        script_lines.append(f'    mat = get_cached_material({c_val}, "AtomMat")')
        script_lines.append(f'    if obj.data.materials: obj.data.materials[0] = mat')
        script_lines.append(f'    else: obj.data.materials.append(mat)')
        script_lines.append(f'    atoms.append(obj)')
    
    script_lines.append('    ')
    script_lines.append('    # Create bonds')
    script_lines.append('    bonds = []')
    
    # Add bond creation code
    for i, bond in enumerate(bonds_data):
        p1 = tuple(float(c) for c in bond['pos1'])
        p2 = tuple(float(c) for c in bond['pos2'])
        radius = float(bond.get('radius', 0.1))
        color = tuple(float(c) for c in bond.get('color', (0.6, 0.6, 0.6)))
        
        script_lines.append(f'    obj = create_bond(')
        script_lines.append(f'        pos1={p1},')
        script_lines.append(f'        pos2={p2},')
        script_lines.append(f'        radius={radius:.3f},')
        script_lines.append(f'        color={color},')
        script_lines.append(f'        name="Bond_{i}"')
        script_lines.append(f'    )')
        script_lines.append(f'    mat = get_cached_material({color}, "BondMat")')
        script_lines.append(f'    if obj.data.materials: obj.data.materials[0] = mat')
        script_lines.append(f'    else: obj.data.materials.append(mat)')
        script_lines.append(f'    bonds.append(obj)')
    
    script_lines.append('    ')
    script_lines.append('    # Join all objects into a collection for organization')
    script_lines.append('    collection_name = "Molecule"')
    script_lines.append('    if collection_name not in bpy.data.collections:')
    script_lines.append('        mol_collection = bpy.data.collections.new(collection_name)')
    script_lines.append('        bpy.context.scene.collection.children.link(mol_collection)')
    script_lines.append('    else:')
    script_lines.append('        mol_collection = bpy.data.collections[collection_name]')
    script_lines.append('    ')
    script_lines.append('    # Move all objects to molecule collection')
    script_lines.append('    for obj in atoms + bonds:')
    script_lines.append('        if obj.name in bpy.context.scene.collection.objects:')
    script_lines.append('            bpy.context.scene.collection.objects.unlink(obj)')
    script_lines.append('        if obj.name not in mol_collection.objects:')
    script_lines.append('            mol_collection.objects.link(obj)')
    script_lines.append('    ')
    script_lines.append('    # Reset 3D view and set shading')
    script_lines.append('    for area in bpy.context.screen.areas:')
    script_lines.append('        if area.type == \'VIEW_3D\':')
    script_lines.append('            # Skip view_all as it overrides the specific camera placement')
    script_lines.append('            pass')
    script_lines.append('            ')
    script_lines.append('            # Set shading to Material Preview (to see colors)')
    script_lines.append('            area.spaces[0].shading.type = \'MATERIAL\'')
    script_lines.append('            # Also set background to \'WORLD\' for consistency')
    script_lines.append('            area.spaces[0].shading.background_type = \'WORLD\'')
    
    # Create camera
    if 'camera_data' in locals():
        script_lines.append('    # Create camera')
        script_lines.append(f'    create_camera(')
        script_lines.append(f'        location=({camera_data["location"][0]:.6f}, {camera_data["location"][1]:.6f}, {camera_data["location"][2]:.6f}),')
        script_lines.append(f'        look_at=({camera_data["look_at"][0]:.6f}, {camera_data["look_at"][1]:.6f}, {camera_data["look_at"][2]:.6f}),')
        script_lines.append(f'        sky={camera_data["sky"]},')
        script_lines.append(f'        angle={camera_data["angle"]:.6f},')
        script_lines.append(f'        orthographic={camera_data["orthographic"]},')
        script_lines.append(f'        parallel_scale={camera_data.get("parallel_scale", 5.0):.6f}')
        script_lines.append(f'    )')
        script_lines.append('    ')
        script_lines.append('    # Switch to camera view')
        script_lines.append('    for area in bpy.context.screen.areas:')
        script_lines.append('        if area.type == \'VIEW_3D\':')
        script_lines.append('            area.spaces[0].region_3d.view_perspective = \'CAMERA\'')
    
    script_lines.append(f'    print("Molecular structure created: {len(atoms_data)} atoms, {len(bonds_data)} bonds")')
    script_lines.append('')
    
    # Run the script
    script_lines.append('# Run the molecule creation')
    script_lines.append('if __name__ == "__main__":')
    script_lines.append('    create_molecule()')
    
    return '\n'.join(script_lines)
