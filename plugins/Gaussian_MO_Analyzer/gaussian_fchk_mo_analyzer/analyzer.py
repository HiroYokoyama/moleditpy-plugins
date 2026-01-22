import numpy as np
import re
import math

class FCHKReader:
    """
    Parses Gaussian FCHK files and holds data in a dictionary.
    """
    def __init__(self, filepath):
        self.filepath = filepath
        self.data = {}
        self.parse()

    def parse(self):
        with open(self.filepath, 'r') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            # Match Array Header: "Atomic numbers   I   N=   12"
            match_array = re.match(r'^(.*)\s+([IRC])\s+N=\s+(\d+)$', line)
            
            # Match Scalar: "Number of electrons   I   42" or "Charge   R   0.0"
            # Note: Scalars usually don't have N=, just the value at the end.
            match_scalar = re.match(r'^(.*)\s+([IRC])\s+([-\d\.E\+]+)$', line)

            if match_array:
                tag = match_array.group(1).strip()
                dtype = match_array.group(2)
                count = int(match_array.group(3))
                
                raw_values = []
                i += 1
                
                # Read data lines until we have enough values
                while len(raw_values) < count and i < len(lines):
                    # Robust splitting for packed Fortran numbers
                    l_str = lines[i].replace('D+', 'E+').replace('D-', 'E-').replace('\n', '').strip()
                    l_str = l_str.replace('-', ' -').replace('E -', 'E-') 
                    vals = l_str.split()
                    raw_values.extend(vals)
                    i += 1
                
                # Convert to numpy array
                if dtype == 'I':
                    self.data[tag] = np.array(raw_values[:count], dtype=int)
                elif dtype == 'R':
                    self.data[tag] = np.array(raw_values[:count], dtype=float)
                else: 
                    self.data[tag] = raw_values[:count]
                    
            elif match_scalar:
                tag = match_scalar.group(1).strip()
                dtype = match_scalar.group(2)
                val_str = match_scalar.group(3)
                
                # Store as single-element array for consistency/easy access
                if dtype == 'I':
                    self.data[tag] = np.array([int(val_str)], dtype=int)
                elif dtype == 'R':
                    try:
                        v = float(val_str.replace('D', 'E'))
                        self.data[tag] = np.array([v], dtype=float)
                    except:
                        pass # Ignore if parse fail
                        
                i += 1
            else:
                i += 1

    def get(self, key, default=None):
        return self.data.get(key, default)

    def get_xyz_block(self):
        """
        Reconstructs an XYZ block string (in Angstroms) for RDKit.
        """
        atoms = self.get("Atomic numbers")
        coords = self.get("Current cartesian coordinates")
        
        if atoms is None or coords is None:
            return None
            
        coords = coords.reshape(-1, 3)
        n_atoms = len(atoms)
        
        # Periodic Table map (simplified or import from rdkit)
        # Since we are in the plugin env, we can try importing rdkit, 
        # or just use a small lookup if we want to keep analyzer pure?
        # Let's import rdkit here as it is guaranteed in MoleditPy environment
        from rdkit import Chem
        pt = Chem.GetPeriodicTable()
        
        BOHR_TO_ANG = 0.529177249
        
        out = [f"{n_atoms}", "Generated from FCHK"]
        for i in range(n_atoms):
            sym = pt.GetElementSymbol(int(atoms[i]))
            x, y, z = coords[i] * BOHR_TO_ANG
            out.append(f"{sym} {x:.6f} {y:.6f} {z:.6f}")
            
        return "\n".join(out)

class BasisSetEngine:
    """
    Constructs basis set from FCHK data and evaluates it on a grid.
    """
    def __init__(self, fchk_data: FCHKReader):
        self.fchk = fchk_data
        self.shells = [] 
        self.n_basis = 0
        self._prepare_basis_set()

    def _normalization_prefactor(self, alpha, l, m, n):
        L = l + m + n
        # Simple implementation for S, P, D
        # Simple implementation for S, P, D, F, G
        fact = {0: 1, 1: 1, 2: 2, 3: 6, 4: 24}
        fact2 = {0: 1, 1: 2, 2: 24, 3: 720, 4: 40320} # (2n)! mapped by n
        
        # FIX: Access fact2 with n (l), not 2*n, because keys recall (2n)!
        num = (8 * alpha)**L * fact[l] * fact[m] * fact[n]
        den = fact2[l] * fact2[m] * fact2[n]
        
        return ((2 * alpha / np.pi)**0.75) * np.sqrt(num / den)

    def _prepare_basis_set(self):
        shell_types = self.fchk.get("Shell types")
        prim_per_shell = self.fchk.get("Number of primitives per shell")
        shell_to_atom = self.fchk.get("Shell to atom map")
        prim_exps = self.fchk.get("Primitive exponents")
        cont_coeffs = self.fchk.get("Contraction coefficients")
        
        # Check for P-modifiers (used for SP shells P-component)
        # Try both tag variations
        p_modifiers = self.fchk.get("P(S=P) Contraction coefficients")
        if p_modifiers is None:
            p_modifiers = self.fchk.get("P(S=P) modifiers")
        
        coords = self.fchk.get("Current cartesian coordinates").reshape(-1, 3) # Bohr
        
        exp_ptr = 0
        coeff_ptr = 0
        
        self.angular_paths = {
            0: [(0,0,0)],
            1: [(1,0,0), (0,1,0), (0,0,1)],
            2: [(2,0,0), (0,2,0), (0,0,2), (1,1,0), (1,0,1), (0,1,1)],
            3: [ # F (10)
                (3,0,0), (0,3,0), (0,0,3), 
                (1,2,0), (2,1,0), (2,0,1), (1,0,2), 
                (0,1,2), (0,2,1), (1,1,1)
            ],
            4: [ # G (15) - Gaussian ordering: ZZZZ, YZZZ, YYZZ, YYYZ, YYYY, XZZZ, XYZZ, XYYZ, XYYY, XXZZ, XXYZ, XXYY, XXXZ, XXXY, XXXX
                 (0,0,4), (0,1,3), (0,2,2), (0,3,1), (0,4,0),
                 (1,0,3), (1,1,2), (1,2,1), (1,3,0),
                 (2,0,2), (2,1,1), (2,2,0),
                 (3,0,1), (3,1,0),
                 (4,0,0)
            ]
        }

        basis_idx_counter = 0

        for i, stype in enumerate(shell_types):
            n_prim = prim_per_shell[i]
            atom_idx = shell_to_atom[i] - 1
            atom_center = coords[atom_idx]

            exps = prim_exps[exp_ptr : exp_ptr + n_prim]
            # S coefficients (always from main array)
            # Note: cont_coeffs is also aligned with primitives?
            # Standard FCHK: "Contraction coefficients" N=N_primitives. 1:1 map.
            coeffs_s = cont_coeffs[coeff_ptr : coeff_ptr + n_prim]
            
            current_shells_to_add = []

            if stype == -1: # SP Shell
                # Check where P coeffs are
                if p_modifiers is not None and len(p_modifiers) > 0:
                    # P coeffs in modifiers array (aligned with exp_ptr)
                    coeffs_p = p_modifiers[exp_ptr : exp_ptr + n_prim]
                    # coeff_ptr consumes S only from main array
                    coeff_ptr += n_prim
                else:
                    # Legacy/Standard fallback: P coeffs follow S coeffs in cont_coeffs
                    # This implies cont_coeffs has 2*N entries for these?
                    # or that they are interleaved?
                    # Standard fallback logic usually implies packed.
                    coeffs_p = cont_coeffs[coeff_ptr + n_prim : coeff_ptr + 2 * n_prim]
                    coeff_ptr += 2 * n_prim
                
                current_shells_to_add.append({
                    'type': 0, 'center': atom_center, 'exps': exps, 'coeffs': coeffs_s, 'lmn': self.angular_paths[0]
                })
                current_shells_to_add.append({
                    'type': 1, 'center': atom_center, 'exps': exps, 'coeffs': coeffs_p, 'lmn': self.angular_paths[1]
                })
                
                exp_ptr += n_prim

            elif stype >= 0: # S, P, D...
                lmn_list = self.angular_paths.get(stype)
                if lmn_list is None:
                    print(f"Warning: Shell type {stype} not supported yet. Skipping.")
                    exp_ptr += n_prim
                    coeff_ptr += n_prim
                    continue

                current_shells_to_add.append({
                    'type': stype, 'center': atom_center, 'exps': exps, 'coeffs': coeffs_s, 'lmn': lmn_list
                })
                
                exp_ptr += n_prim
                coeff_ptr += n_prim

            for sh in current_shells_to_add:
                # 1. Normalize primitives (compute N_i)
                # 2. Normalize contraction (compute N_cont)
                
                # Pre-compute primitive norms for this angular momentum
                # Assuming all components in the shell (e.g. Px, Py, Pz) share the same radial part normalization?
                # Yes, L is constant for the shell.
                l_sample, m_sample, n_sample = sh['lmn'][0]
                total_L = l_sample + m_sample + n_sample
                
                prim_norms = [self._normalization_prefactor(a, l_sample, m_sample, n_sample) for a in sh['exps']]
                prim_norms = np.array(prim_norms)
                
                # Check Contraction Normalization
                # <Chi|Chi> = sum_ij c_i c_j <phi_i|phi_j>
                # <phi_i|phi_j> = S_ij (overlap of normalized primitives)
                # Overlap of normalized primitives with exponents a, b:
                # S_ab = (2*sqrt(a*b)/(a+b))**(L+1.5) ?
                # Let's check formula.
                # For S-orbitals: (2*sqrt(ab)/(a+b))^3/2
                # For general L: (2*sqrt(ab)/(a+b))^(L+1.5)
                # Provided they are centered at same origin (which they are).
                
                norm_sq = 0.0
                n_prims = len(sh['exps'])
                c = sh['coeffs']
                alpha = sh['exps']
                
                # Calculate Contraction Norm Square
                for idx_i in range(n_prims):
                    for idx_j in range(n_prims):
                        ai = alpha[idx_i]
                        aj = alpha[idx_j]
                        overlap = (2.0 * np.sqrt(ai * aj) / (ai + aj)) ** (total_L + 1.5)
                        norm_sq += c[idx_i] * c[idx_j] * prim_norms[idx_i] * prim_norms[idx_j] * overlap
                
                contraction_norm = np.sqrt(norm_sq)
                
                # Scale coefficients so that the final function is normalized
                # We need final coeff C_final_i = C_input_i * PrimNorm_i / ContractionNorm
                
                scale = 1.0 / contraction_norm if contraction_norm > 1e-12 else 1.0
                
                norm_coeffs = []
                for (l, m, n) in sh['lmn']:
                    # We re-calculate primitive norm for specific L,M,N just in case factor depends on it
                    # (Actually _normalization_prefactor DOES depend on l,m,n factorial distribution!)
                    # Wait, do we normalize Px, Py, Pz differently?
                    # The formula relies on L+m+n. All P have same total L.
                    # Factorials: Px (1,0,0) -> 1!0!0! / (2!0!0!) = 1/2. 
                    # Py (0,1,0) -> same. 
                    # Dxy (1,1,0) vs Dxx (2,0,0)?
                    # Dxy: 1!1!0! / (2!2!0!) = 1/4
                    # Dxx: 2!0!0! / (4!0!0!) = 2/24 = 1/12
                    # So Cartesian D norms ARE different.
                    
                    # However, the overlap formula used above assumed generic radial overlap.
                    # The contraction coefficients in FCHK are usually shared for the whole shell.
                    # This implies the contraction applies to the RADIAL part primarily?
                    # If I normalize the Radial part times Spherical harmonic... 
                    # But here I am using Cartesian.
                    
                    # Standard Strategy:
                    # Normalize each Cartesian component independently. 
                    # Because they are orthogonal (mostly).
                    # Actually Px, Py, Pz are orthogonal.
                    # So we can normalize "Px" contraction.
                    
                    real_prim_norms = [self._normalization_prefactor(a, l, m, n) for a in sh['exps']]
                    real_prim_norms = np.array(real_prim_norms)
                    
                    # Recalculate contraction norm for THIS component (l,m,n)
                    # Because prim_norm changes, the 'N*N' term in the sum changes.
                    
                    
                    # comp_norm_sq = 0.0
                    # for idx_i in range(n_prims):
                    #     for idx_j in range(n_prims):
                    #         ai = alpha[idx_i]
                    #         aj = alpha[idx_j]
                    #         overlap = (2.0 * np.sqrt(ai * aj) / (ai + aj)) ** ((l+m+n) + 1.5)
                    #         comp_norm_sq += c[idx_i] * c[idx_j] * real_prim_norms[idx_i] * real_prim_norms[idx_j] * overlap
                    
                    # comp_contraction_norm = np.sqrt(comp_norm_sq)
                    # comp_scale = 1.0 / comp_contraction_norm if comp_contraction_norm > 1e-12 else 1.0

                    # Gaussian FCHK coefficients are usually already normalized.
                    comp_scale = 1.0
                    
                    # Store final coefficients: c_input * prim_norm * scale
                    c_final = sh['coeffs'] * real_prim_norms * comp_scale
                    norm_coeffs.append(c_final)
                
                sh['norm_coeffs'] = np.array(norm_coeffs)
                sh['start_idx'] = basis_idx_counter
                basis_idx_counter += len(sh['lmn'])
                
                self.shells.append(sh)
        
        self.n_basis = basis_idx_counter

    def evaluate_mo_on_grid(self, mo_idx: int, grid_coords: np.ndarray, mo_coeffs_all: np.ndarray, progress_callback=None):
        start = mo_idx * self.n_basis
        end = start + self.n_basis
        if end > len(mo_coeffs_all):
             raise ValueError("MO index out of range or Basis count mismatch")
        
        my_coeffs = mo_coeffs_all[start:end]
        phi_mo = np.zeros(grid_coords.shape[0])
        
        for i, sh in enumerate(self.shells):
            r_vec = grid_coords - sh['center']
            r2 = np.sum(r_vec**2, axis=1)
            E = np.exp(-np.outer(sh['exps'], r2))
            
            current_basis_idx = sh['start_idx']
            
            for ang_idx, (l, m, n) in enumerate(sh['lmn']):
                c_mo = my_coeffs[current_basis_idx]
                current_basis_idx += 1
                
                if abs(c_mo) < 1e-8:
                    continue

                ang_val = np.ones(grid_coords.shape[0])
                if l > 0: ang_val *= r_vec[:, 0]**l
                if m > 0: ang_val *= r_vec[:, 1]**m
                if n > 0: ang_val *= r_vec[:, 2]**n
                
                contracted_radial = np.dot(sh['norm_coeffs'][ang_idx], E)
                phi_mo += c_mo * ang_val * contracted_radial
                
        return phi_mo

    def get_basis_labels(self):
        """
        Returns a list of labels for the basis functions, e.g. "1 C 1s", "1 C 2px"...
        """
        labels = []
        
        # Shell path labels
        # 0: S
        # 1: Px, Py, Pz
        # 2: XX, YY, ZZ, XY, XZ, YZ (6D)
        
        type_map = {
            0: ["S"],
            1: ["Px", "Py", "Pz"],
            2: ["dXX", "dYY", "dZZ", "dXY", "dXZ", "dYZ"],
            3: ["fXXX", "fYYY", "fZZZ", "fXYY", "fXXY", "fXXZ", "fXZZ", "fYZZ", "fYYZ", "fXYZ"],
            4: ["gZZZZ", "gYZZZ", "gYYZZ", "gYYYZ", "gYYYY", "gXZZZ", "gXYZZ", "gXYYZ", "gXYYY", "gXXZZ", "gXXYZ", "gXXYY", "gXXXZ", "gXXXY", "gXXXX"]
        }
        
        # We need atom symbol map
        atom_list = self.fchk.get("Atomic numbers")
        
        # Need to reconstruct which shell belongs to which atom
        # We can replay the shell loop
        shell_to_atom = self.fchk.get("Shell to atom map")
        shell_types = self.fchk.get("Shell types")
        
        # Periodic Table 
        # (Simple lookup or import. Since we used rdkit in FCHKReader, we assume availability or fallback)
        try:
            from rdkit import Chem
            pt = Chem.GetPeriodicTable()
            to_sym = lambda z: pt.GetElementSymbol(int(z))
        except:
            to_sym = lambda z: f"Atom{z}"

        for i, stype in enumerate(shell_types):
            atom_idx = shell_to_atom[i] - 1
            atom_num = atom_list[atom_idx]
            sym = to_sym(atom_num)
            
            # Determine n (principal quantum number) guess? 
            # Hard to parse exact n from just S/P types without electron config map.
            # But we can just say "C1 S", "C1 Px"...
            
            comps = []
            if stype == -1: # SP
                comps = ["S", "Px", "Py", "Pz"]
            else:
                comps = type_map.get(stype, [f"?{stype}?"])
                
            for c in comps:
                labels.append(f"{atom_idx+1} {sym} {c}")
                
        return labels

    def evaluate_basis_on_grid(self, basis_idx: int, grid_coords: np.ndarray, progress_callback=None):
        """
        Evaluates a single basis function (AO) on the grid.
        """
        phi = np.zeros(grid_coords.shape[0])
        
        # Find which shell this basis_idx belongs to
        # We need to iterate to find the range
        
        current_idx = 0
        target_shell = None
        target_ang_idx = 0
        
        for sh in self.shells:
            n_funcs = len(sh['lmn'])
            if current_idx + n_funcs > basis_idx:
                target_shell = sh
                target_ang_idx = basis_idx - current_idx
                break
            current_idx += n_funcs
            
        if target_shell is None:
            return phi
            
        sh = target_shell
        l, m, n = sh['lmn'][target_ang_idx]
        
        r_vec = grid_coords - sh['center']
        r2 = np.sum(r_vec**2, axis=1)
        E = np.exp(-np.outer(sh['exps'], r2))
        
        # Contraction
        # AO = sum_k ( d_k * N_k * x^l... * exp )
        #    = (x^l...) * sum_k ( coeff_k * E_k )
        
        ang_val = np.ones(grid_coords.shape[0])
        if l > 0: ang_val *= r_vec[:, 0]**l
        if m > 0: ang_val *= r_vec[:, 1]**m
        if n > 0: ang_val *= r_vec[:, 2]**n
        
        contracted_radial = np.dot(sh['norm_coeffs'][target_ang_idx], E)
        phi = ang_val * contracted_radial
        
        return phi

class CubeWriter:
    @staticmethod
    def write(filepath, atoms, atom_nos, origin, vectors, data, comment=""):
        nx, ny, nz = data.shape
        n_atoms = len(atoms)
        
        with open(filepath, 'w') as f:
            f.write(f"Moleditpy Cube File: {comment}\n")
            f.write("Generated by FCHK-to-CUBE Plugin\n")
            f.write(f"{n_atoms:5d} {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
            f.write(f"{nx:5d} {vectors[0,0]:12.6f} {vectors[0,1]:12.6f} {vectors[0,2]:12.6f}\n")
            f.write(f"{ny:5d} {vectors[1,0]:12.6f} {vectors[1,1]:12.6f} {vectors[1,2]:12.6f}\n")
            f.write(f"{nz:5d} {vectors[2,0]:12.6f} {vectors[2,1]:12.6f} {vectors[2,2]:12.6f}\n")
            
            for i in range(n_atoms):
                z_no = atom_nos[i]
                f.write(f"{z_no:5d} {float(z_no):12.6f} {atoms[i][0]:12.6f} {atoms[i][1]:12.6f} {atoms[i][2]:12.6f}\n")
            
            vals = data.flatten()
            for i, v in enumerate(vals):
                f.write(f"{v:13.5E}")
                if (i + 1) % 6 == 0:
                    f.write("\n")
            if len(vals) % 6 != 0:
                f.write("\n")
