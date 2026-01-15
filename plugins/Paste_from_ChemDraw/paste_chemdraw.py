from PyQt6.QtWidgets import QApplication, QMessageBox
from PyQt6.QtCore import QPointF, QTimer
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

# --- Plugin Basic Information ---
PLUGIN_NAME = "Paste from ChemDraw"
PLUGIN_VERSION = "2026.01.15"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Paste chemical structures from ChemDraw clipboard data. Developed on ChemDraw version 25.5."

def initialize(context):
    """
    Modern plugin initialization hook.
    Registers the paste action in the Edit menu with Ctrl+Shift+V shortcut.
    """
    context.add_menu_action(
        "Edit/Paste from ChemDraw",
        lambda: run(context.get_main_window()),
        shortcut="Ctrl+Shift+V"
    )

def run(main_window):
    """
    ChemDrawからコピーされたMDLCTデータ（MOL形式テキスト）を
    クリップボードから取得し、RDKit経由でキャンバスに挿入します。
    """
    
    # 1. クリップボード準備
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
        
    clipboard = app.clipboard()
    mime_data = clipboard.mimeData()
    
    # Check multiple formats
    TARGET_FORMATS = [
        'application/x-qt-windows-mime;value="MDLCT"',
        'application/x-qt-windows-mime;value="MDLSK"'
    ]
    
    mol = None
    mol_text = None
    fmt_found = None

    # 2. データの取得
    for fmt in TARGET_FORMATS:
        if mime_data.hasFormat(fmt):
            try:
                byte_data = mime_data.data(fmt)
                # Use cp932 (Shift-JIS) for Windows clipboard compatibility
                text = byte_data.data().decode('cp932', errors='ignore')
                
                # If V2000 exists, it's a strong candidate
                if text and "V2000" in text:
                    mol_text = text
                    fmt_found = fmt
                    break
                elif text and mol_text is None:
                    # Weak candidate
                    mol_text = text
                    fmt_found = fmt
            except Exception as e:
                print(f"[{PLUGIN_NAME}] Decode Error for {fmt}: {e}")

    # Log debug info
    # try:
    #     with open("chemdraw_paste_debug.txt", "w", encoding="utf-8") as f:
    #         if mol_text:
    #             f.write(f"Format Found: {fmt_found}\n")
    #             f.write(f"Raw Length: {len(mol_text)}\n")
    #             f.write(f"Content Sample (First 500 chars):\n{mol_text[:500]}\n")
    #         else:
    #             f.write("No matching format found.\n")
    # except:
    #    pass

    # Parsing Attempt
    if mol_text:
        try:
            lines = mol_text.splitlines()
            found_v2000 = False
            
            # Case 1: Flat text (suspicious line count)
            if len(lines) < 3:
                 # print(f"[{PLUGIN_NAME}] Suspicious line count ({len(lines)}). Attempting reconstruction...")
                 mol = reconstruct_from_flat_text(mol_text)
                 if mol: found_v2000 = True
            
            # Case 2: Standard parse if reconstruction didn't run or fail
            if not found_v2000:
                for i, line in enumerate(lines):
                     if "V2000" in line:
                         found_v2000 = True
                         start = max(0, i-3)
                         mol_text_clean = "\n".join(lines[start:])
                         mol = Chem.MolFromMolBlock(mol_text_clean)
                         break

            # Case 3: Fallback direct parse
            if not found_v2000 and mol is None:
                 mol = Chem.MolFromMolBlock(mol_text)
                 if mol is None:
                     mol = reconstruct_from_flat_text(mol_text)

        except Exception as e:
            QMessageBox.critical(main_window, PLUGIN_NAME, f"Paste Error: {e}")
            # print(f"[{PLUGIN_NAME}] Error: {e}")

    # 3. Text fallback
    if mol is None and mime_data.hasText():
        try:
            text = mime_data.text().strip()
            if "V2000" in text or "M  END" in text:
                mol = Chem.MolFromMolBlock(text)
                if mol is None:
                     mol = reconstruct_from_flat_text(text)
        except:
            pass

    # 4. Drawing Logic
    if mol is not None:
        try:
            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)
            
            try:
                Chem.Kekulize(mol, clearAromaticFlags=True)
                AllChem.AssignStereochemistry(mol, force=True)
            except:
                pass

            SCALE_FACTOR = 40.0
            view_center = QPointF(0, 0)
            if hasattr(main_window, 'view_2d') and main_window.view_2d:
                 viewport_rect = main_window.view_2d.viewport().rect()
                 view_center = main_window.view_2d.mapToScene(viewport_rect.center())
            
            conf = mol.GetConformer()
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            cx = sum(p.x for p in positions) / len(positions) if positions else 0
            cy = sum(p.y for p in positions) / len(positions) if positions else 0

            rdkit_idx_to_item = {}

            if hasattr(main_window, 'scene') and main_window.scene:
                # Create Atoms
                for i in range(mol.GetNumAtoms()):
                    atom = mol.GetAtomWithIdx(i)
                    pos = conf.GetAtomPosition(i)
                    sx = (pos.x - cx) * SCALE_FACTOR + view_center.x()
                    sy = -(pos.y - cy) * SCALE_FACTOR + view_center.y() 
                    
                    if hasattr(main_window.scene, 'create_atom'):
                        a_id = main_window.scene.create_atom(
                            atom.GetSymbol(), 
                            QPointF(sx, sy), 
                            charge=atom.GetFormalCharge()
                        )
                        if hasattr(main_window, 'data') and a_id in main_window.data.atoms:
                            rdkit_idx_to_item[i] = main_window.data.atoms[a_id]['item']

                # Create Bonds
                for bond in mol.GetBonds():
                    idx1 = bond.GetBeginAtomIdx()
                    idx2 = bond.GetEndAtomIdx()
                    
                    if idx1 in rdkit_idx_to_item and idx2 in rdkit_idx_to_item:
                        item1 = rdkit_idx_to_item[idx1]
                        item2 = rdkit_idx_to_item[idx2]
                        order = int(bond.GetBondTypeAsDouble())
                        
                        stereo = 0
                        dir_ = bond.GetBondDir()
                        if dir_ == Chem.BondDir.BEGINWEDGE: stereo = 1
                        elif dir_ == Chem.BondDir.BEGINDASH: stereo = 2
                        
                        if hasattr(main_window.scene, 'create_bond'):
                            main_window.scene.create_bond(item1, item2, bond_order=order, bond_stereo=stereo)
            
                if hasattr(main_window, 'has_unsaved_changes'):
                    main_window.has_unsaved_changes = True
                if hasattr(main_window, 'update_window_title'):
                    main_window.update_window_title()
                
                # Redraw 2D scene to ensure implicit hydrogens are counted
                # Use QTimer to defer the update so Qt can process the new items first
                def update_scene():
                    if hasattr(main_window, 'scene') and main_window.scene:
                        # Trigger geometry recalculation (main code handles implicit H)
                        for atom_item in rdkit_idx_to_item.values():
                            if hasattr(atom_item, 'update'):
                                atom_item.update()
                        main_window.scene.update()
                        
                QTimer.singleShot(100, update_scene)  # 100ms delay

        except Exception as e:
            QMessageBox.critical(main_window, PLUGIN_NAME, f"Paste Error: {e}")
            # print(f"[{PLUGIN_NAME}] Error: {e}")
    else:
        # Failure Message
        msg = "No valid MDLCT data found."
        # if mol_text:
        #      msg += f"\n(Data found in {fmt_found} but parsing failed.\nSee chemdraw_paste_debug.txt)"
        QMessageBox.warning(main_window, PLUGIN_NAME, msg)


def reconstruct_from_flat_text(text):
    """
    ChemDraw clipboard text sometimes lacks newlines and has phantom tokens.
    Reconstructs a valid MolBlock from such flat text.
    """
    import re
    # debug_log = []
    # debug_log.append(f"[{PLUGIN_NAME}] Attempting reconstruction...")
    
    # Sanitize: Remove non-printable characters (control chars like \x01, \x16 etc.)
    # Keep only standard printable ASCII (0x20-0x7E) and whitespace (\r, \n, \t)
    original_len = len(text)
    text = re.sub(r'[^\x20-\x7E\s]', '', text)
    # if len(text) != original_len:
    #      debug_log.append(f"Sanitized input: removed {original_len - len(text)} non-printable characters.")

    try:
        idx_v2000 = text.find("V2000")
        if idx_v2000 == -1: 
            # debug_log.append("V2000 not found")
            # _write_debug_log(debug_log)
            return None
        
        header_part = text[:idx_v2000]
        rest_part = text[idx_v2000+5:]
        
        header_tokens = header_part.split()
        # debug_log.append(f"Header tokens: {header_tokens}")
        
        if len(header_tokens) < 10: 
            # debug_log.append("Not enough header tokens")
            # _write_debug_log(debug_log)
            return None
            
        try:
            n_atoms = int(header_tokens[-10])
            n_bonds = int(header_tokens[-9])
            # debug_log.append(f"Atoms: {n_atoms}, Bonds: {n_bonds}")
        except ValueError as e:
            # debug_log.append(f"Failed to parse counts: {e}")
            # _write_debug_log(debug_log)
            return None
        
        new_mol = ["Header", "  ChemDraw", "", f"{n_atoms:3d}{n_bonds:3d}  0  0  0  0  0  0  0  0999 V2000"]
        
        raw_tokens = rest_part.split()
        # debug_log.append(f"Raw tokens (first 20): {raw_tokens[:20]}")
        
        current_idx = 0
        
        for i in range(n_atoms):
            if current_idx >= len(raw_tokens): 
                 # debug_log.append(f"Ran out of tokens at atom {i}")
                 break
            
            # Skip 'F' phantom token logic
            if raw_tokens[current_idx] == 'F':
                current_idx += 1
            # Try to skip obvious garbage (non-numeric that is not a symbol)?
            # But be careful not to skip symbols like 'C', 'N'.
            # For now, stick to 'F' skip.
                
            x = raw_tokens[current_idx]; current_idx += 1
            y = raw_tokens[current_idx]; current_idx += 1
            z = raw_tokens[current_idx]; current_idx += 1
            sym = raw_tokens[current_idx]; current_idx += 1
            
            # Consume 12 zeros (V2000 standard fields)
            rest = []
            for _ in range(12):
                if current_idx < len(raw_tokens):
                    rest.append(raw_tokens[current_idx])
                    current_idx += 1
            
            line = f"{float(x):10.4f}{float(y):10.4f}{float(z):10.4f} {sym:<3}" + \
                   "".join([f"{int(val):3d}" for val in rest])
            new_mol.append(line)
            
        for i in range(n_bonds):
            if current_idx >= len(raw_tokens): break
            a1 = raw_tokens[current_idx]; current_idx += 1
            a2 = raw_tokens[current_idx]; current_idx += 1
            type = raw_tokens[current_idx]; current_idx += 1
            stereo = raw_tokens[current_idx]; current_idx += 1
            
            line = f"{int(a1):3d}{int(a2):3d}{int(type):3d}{int(stereo):3d}  0  0  0"
            new_mol.append(line)
            
            # Consume remaining zeros until next non-zero or M tag
            while current_idx < len(raw_tokens) and raw_tokens[current_idx] == '0':
                 current_idx += 1
        
        # Preserve M lines (M  CHG, M  ISO, M  RAD, etc.) from original text
        # These lines contain formal charges and other properties
        m_lines = []
        while current_idx < len(raw_tokens):
            if raw_tokens[current_idx] == 'M':
                # Collect tokens until we hit another M or END
                line_tokens = [raw_tokens[current_idx]]
                current_idx += 1
                
                # Special handling for M  END
                if current_idx < len(raw_tokens) and raw_tokens[current_idx] == 'END':
                    break  # Don't add M END to m_lines, we'll add it separately
                    
                # Collect the rest of this M line
                while current_idx < len(raw_tokens):
                    token = raw_tokens[current_idx]
                    if token == 'M':  # Start of next M line
                        break
                    line_tokens.append(token)
                    current_idx += 1
                    
                # Reconstruct this M line
                m_lines.append('  '.join(line_tokens))
            else:
                current_idx += 1
        
        # Add M lines before M  END
        for m_line in m_lines:
            new_mol.append(m_line)
        
        new_mol.append("M  END")
        
        full_block = "\n".join(new_mol)
        # debug_log.append("Reconstruction successful (raw block created).")
        # debug_log.append(f"Block Preview:\n{full_block}")
        
        # _write_debug_log(debug_log)
        return Chem.MolFromMolBlock(full_block)
        
    except Exception as e:
        # debug_log.append(f"Reconstruction Exception: {e}")
        import traceback
        # debug_log.append(traceback.format_exc())
        # _write_debug_log(debug_log)
        return None

# def _write_debug_log(lines):
#     try:
#         with open("chemdraw_paste_debug.txt", "a", encoding="utf-8") as f:
#             f.write("\n--- Reconstruction Log ---\n")
#             f.write("\n".join(lines))
#             f.write("\n--------------------------\n")
#     except:
#         pass
