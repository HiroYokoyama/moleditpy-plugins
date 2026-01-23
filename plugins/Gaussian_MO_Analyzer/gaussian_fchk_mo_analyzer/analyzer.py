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
    Supports both Cartesian and Spherical basis functions.
    """
    def __init__(self, fchk_data: FCHKReader):
        self.fchk = fchk_data
        self.shells = [] 
        self.n_basis = 0
        self._prepare_basis_set()

    def _normalization_prefactor(self, alpha, l, m, n):
        L = l + m + n
        # Simple implementation for S, P, D, F, G
        fact = {0: 1, 1: 1, 2: 2, 3: 6, 4: 24, 5: 120}
        fact2 = {0: 1, 1: 2, 2: 24, 3: 720, 4: 40320, 5: 3628800} # (2n)! mapped by n
        
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
        
        # Pure/Cartesian Flags (0=Pure/Spherical, 1=Cartesian)
        # Default to 1 (Cartesian) if not found, to match old behavior
        d_is_cart = True
        f_is_cart = True
        
        flag_d = self.fchk.get("Pure/Cartesian d shells")
        if flag_d is not None and len(flag_d) > 0 and flag_d[0] == 0:
            d_is_cart = False
            
        flag_f = self.fchk.get("Pure/Cartesian f shells")
        if flag_f is not None and len(flag_f) > 0 and flag_f[0] == 0:
            f_is_cart = False

        # Build Basis Definitions (Linear Combinations)
        # Structure: Type -> [  List of Basis Functions in this Shell ]
        # Each Basis Function -> List of Components (coeff, (l,m,n))
        
        # Cartesian Definitions (Simple 1:1 map)
        cart_s = [ [(1.0, (0,0,0))] ] # S
        cart_p = [ [(1.0, (1,0,0))], [(1.0, (0,1,0))], [(1.0, (0,0,1))] ] # Px, Py, Pz
        cart_d = [ 
            [(1.0, (2,0,0))], [(1.0, (0,2,0))], [(1.0, (0,0,2))], 
            [(1.0, (1,1,0))], [(1.0, (1,0,1))], [(1.0, (0,1,1))]
        ] # XX, YY, ZZ, XY, XZ, YZ
        cart_f = [
            [(1.0, (3,0,0))], [(1.0, (0,3,0))], [(1.0, (0,0,3))],
            [(1.0, (1,2,0))], [(1.0, (2,1,0))], [(1.0, (2,0,1))],
            [(1.0, (1,0,2))], [(1.0, (0,1,2))], [(1.0, (0,2,1))],
            [(1.0, (1,1,1))]
        ] # XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ
        
        self.basis_definitions = {
            0: cart_s,
            1: cart_p,
            -1: cart_s + cart_p # SP
        }
        
        # Handle D Shells
        if d_is_cart:
            self.basis_definitions[2] = cart_d
        else:
            # Spherical D (5 components)
            # Gaussian Spherical D ordering: d(0), d(+1), d(-1), d(+2), d(-2)
            # Mapped as:
            # 1. d(3z^2-r^2) -> 2z^2 - x^2 - y^2 (Need -0.5, -0.5, 1.0)
            # 2. d(xz)
            # 3. d(yz)
            # 4. d(x^2-y^2) -> sqrt(3)/2 ?? Wait, standard real harmonic is X^2-Y^2.
            #    However, normalization implies specific coeffs.
            #    Standard un-normalized Cartesian form:
            #    d(0) = 2z^2 - x^2 - y^2
            #    d(x^2-y^2) = x^2 - y^2
            #    BUT Gaussian internal might require normalized angular parts.
            #    Let's use the coefficients that construct the SHAPE correctly.
            #    The radial parts are normalized separately.
            
            # Using common conversion (approximate for shape):
            # Spherical D (5 components)
            # Normalized to sum-of-squares of weights = 1 (assuming normalized Cartesians)
            # Calculated via calculate_norm.py accounting for overlaps.
            
            # D0: 2z^2 - x^2 - y^2. Weights (2, -1, -1). Norm factor 0.5.
            # D2: x^2 - y^2. Norm factor sqrt(3)/2 = 0.8660254
            
            sph_d = [
                [ (-0.5, (2,0,0)), (-0.5, (0,2,0)), (1.0, (0,0,2)) ], # d(3z^2-r^2) -> 0.5*(2z^2-x^2-y^2) gives -0.5x^2... Correct.
                [ (1.0, (1,0,1)) ], # d(xz)
                [ (1.0, (0,1,1)) ], # d(yz)
                [ (0.8660254037844387, (2,0,0)), (-0.8660254037844387, (0,2,0)) ], # d(x^2-y^2)
                [ (1.0, (1,1,0)) ]  # d(xy)
            ]
            self.basis_definitions[2] = sph_d

        # Handle F Shells
        if f_is_cart:
            self.basis_definitions[3] = cart_f
        else:
            # Spherical F (7 components) - Normalized
            # Calculated via calculate_norm.py accounting for overlaps.
            
            # F0 (5z^2-3r^2)z -> 2z^3 - 3x^2z - 3y^2z. Norm Factor: 0.24065403
            f0_n = 0.24065403
            # F1 (5z^2-r^2)x -> 4xz^2 - x^3 - xy^2. Norm Factor: 0.28116020
            f1_n = 0.28116020
            # F2 (x^2-y^2)z -> x^2z - y^2z. Norm Factor: 0.86602540
            f2_n = 0.86602540
            # F3 (x^3-3xy^2). Norm Factor: 0.36969351
            f3_n = 0.36969351
            
            sph_f = [
                # f(0): 2ZZZ - 3XXZ - 3YYZ.
                [ (2.0*f0_n, (0,0,3)), (-3.0*f0_n, (2,0,1)), (-3.0*f0_n, (0,2,1)) ], 
                # 1: x(5z^2-r^2) -> 4XZZ - XXX - XYY
                [ (4.0*f1_n, (1,0,2)), (-1.0*f1_n, (3,0,0)), (-1.0*f1_n, (1,2,0)) ], 
                # 2: y(5z^2-r^2) -> 4YZZ - XXY - YYY
                [ (4.0*f1_n, (0,1,2)), (-1.0*f1_n, (2,1,0)), (-1.0*f1_n, (0,3,0)) ],
                # 3: z(x^2-y^2) -> X2Z - Y2Z
                [ (1.0*f2_n, (2,0,1)), (-1.0*f2_n, (0,2,1)) ],
                # 4: xyz
                [ (1.0, (1,1,1)) ],
                # 5: x(x^2-3y^2) -> XXX - 3XYY
                [ (1.0*f3_n, (3,0,0)), (-3.0*f3_n, (1,2,0)) ],
                # 6: y(3x^2-y^2) -> 3XXY - YYY
                [ (3.0*f3_n, (2,1,0)), (-1.0*f3_n, (0,3,0)) ]
            ]
            self.basis_definitions[3] = sph_f

        # Check for P-modifiers (used for SP shells P-component)
        # Try both tag variations
        p_modifiers = self.fchk.get("P(S=P) Contraction coefficients")
        if p_modifiers is None:
            p_modifiers = self.fchk.get("P(S=P) modifiers")
        
        coords = self.fchk.get("Current cartesian coordinates").reshape(-1, 3) # Bohr
        
        exp_ptr = 0
        coeff_ptr = 0
        basis_idx_counter = 0

        for i, stype in enumerate(shell_types):
            n_prim = prim_per_shell[i]
            atom_idx = shell_to_atom[i] - 1
            atom_center = coords[atom_idx]

            exps = prim_exps[exp_ptr : exp_ptr + n_prim]
            coeffs_s = cont_coeffs[coeff_ptr : coeff_ptr + n_prim]
            
            current_shells_to_add = []

            if stype == -1: # SP Shell
                if p_modifiers is not None and len(p_modifiers) > 0:
                    coeffs_p = p_modifiers[exp_ptr : exp_ptr + n_prim]
                    coeff_ptr += n_prim
                else:
                    coeffs_p = cont_coeffs[coeff_ptr + n_prim : coeff_ptr + 2 * n_prim]
                    coeff_ptr += 2 * n_prim
                
                current_shells_to_add.append({
                    'type': 0, 'center': atom_center, 'exps': exps, 'coeffs': coeffs_s, 'defs': self.basis_definitions[0]
                })
                current_shells_to_add.append({
                    'type': 1, 'center': atom_center, 'exps': exps, 'coeffs': coeffs_p, 'defs': self.basis_definitions[1]
                })
                
                exp_ptr += n_prim

            else: # S, P, D, F... AND -2(D), -3(F) handling
                # Handle negative types -2, -3 etc which are just D, F...
                effective_type = stype
                if effective_type < -1:
                    effective_type = abs(effective_type)
                
                defs_list = self.basis_definitions.get(effective_type)
                if defs_list is None:
                    print(f"Warning: Shell type {stype} (Effective {effective_type}) not supported. Skipping.")
                    exp_ptr += n_prim
                    coeff_ptr += n_prim
                    continue

                current_shells_to_add.append({
                    'type': effective_type, 'center': atom_center, 'exps': exps, 'coeffs': coeffs_s, 'defs': defs_list
                })
                
                exp_ptr += n_prim
                coeff_ptr += n_prim

            for sh in current_shells_to_add:
                func_norm_coeffs = []
                
                # Precompute coefficients for efficient evaluation
                for basis_func_def in sh['defs']:
                    comps_data = []
                    for weight, (l, m, n) in basis_func_def:
                        # Calculate prefactors for this l,m,n for all exponents
                        prim_norms = np.array([self._normalization_prefactor(a, l, m, n) for a in sh['exps']])
                        
                        # Apply Contraction Coeffs & Weights
                        comp_coeffs = sh['coeffs'] * prim_norms * weight
                        
                        comps_data.append( {
                            'l': l, 'm': m, 'n': n,
                            'coeffs': comp_coeffs 
                        })
                    
                    func_norm_coeffs.append(comps_data)
                
                sh['basis_data'] = func_norm_coeffs
                sh['start_idx'] = basis_idx_counter
                basis_idx_counter += len(sh['defs'])
                
                self.shells.append(sh)
        
        self.n_basis = basis_idx_counter

    def evaluate_mo_on_grid(self, mo_idx: int, grid_coords: np.ndarray, mo_coeffs_all: np.ndarray, progress_callback=None):
        start = mo_idx * self.n_basis
        end = start + self.n_basis
        if end > len(mo_coeffs_all):
             raise ValueError(f"MO index {mo_idx} out of range (need {end} coeffs, have {len(mo_coeffs_all)})")
        
        my_coeffs = mo_coeffs_all[start:end]
        phi_mo = np.zeros(grid_coords.shape[0])
        
        for i, sh in enumerate(self.shells):
            r_vec = grid_coords - sh['center']
            r2 = np.sum(r_vec**2, axis=1)
            E = np.exp(-np.outer(sh['exps'], r2))
            
            current_basis_idx = sh['start_idx']
            
            # Iterate over Basis Functions in this Shell
            for basis_func_data in sh['basis_data']:
                c_mo = my_coeffs[current_basis_idx]
                current_basis_idx += 1
                
                if abs(c_mo) < 1e-8:
                    continue
                
                val_accum = np.zeros(grid_coords.shape[0])
                
                for comp in basis_func_data:
                    l, m, n = comp['l'], comp['m'], comp['n']
                    coeffs_prim = comp['coeffs']
                    
                    ang_val = np.ones(grid_coords.shape[0])
                    if l > 0: ang_val *= r_vec[:, 0]**l
                    if m > 0: ang_val *= r_vec[:, 1]**m
                    if n > 0: ang_val *= r_vec[:, 2]**n
                    
                    contracted_radial = np.dot(coeffs_prim, E)
                    val_accum += ang_val * contracted_radial
                
                phi_mo += c_mo * val_accum
                
        return phi_mo

    def get_basis_labels(self):
        """
        Returns a list of labels for the basis functions.
        """
        labels = []
        atom_list = self.fchk.get("Atomic numbers")
        shell_to_atom = self.fchk.get("Shell to atom map")
        shell_types = self.fchk.get("Shell types")
        
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
            
            comps = []
            if stype == -1: # SP
                comps = ["S", "Px", "Py", "Pz"]
            elif stype >= 0:
                count = len(self.basis_definitions.get(stype, []))
                comps = [f"L{stype}_{k}" for k in range(count)]
                
            for c in comps:
                labels.append(f"{atom_idx+1} {sym} {c}")
                
        return labels

    def evaluate_basis_on_grid(self, basis_idx: int, grid_coords: np.ndarray, progress_callback=None):
        """
        Evaluates a single basis function (AO) on the grid.
        """
        phi = np.zeros(grid_coords.shape[0])
        
        current_idx = 0
        target_shell = None
        target_func_idx = 0
        
        for sh in self.shells:
            n_funcs = len(sh['basis_data'])
            if current_idx + n_funcs > basis_idx:
                target_shell = sh
                target_func_idx = basis_idx - current_idx
                break
            current_idx += n_funcs
            
        if target_shell is None:
            return phi
            
        sh = target_shell
        basis_func_data = sh['basis_data'][target_func_idx]
        
        r_vec = grid_coords - sh['center']
        r2 = np.sum(r_vec**2, axis=1)
        E = np.exp(-np.outer(sh['exps'], r2))
        
        for comp in basis_func_data:
            l, m, n = comp['l'], comp['m'], comp['n']
            coeffs_prim = comp['coeffs']
            
            ang_val = np.ones(grid_coords.shape[0])
            if l > 0: ang_val *= r_vec[:, 0]**l
            if m > 0: ang_val *= r_vec[:, 1]**m
            if n > 0: ang_val *= r_vec[:, 2]**n
            
            contracted_radial = np.dot(coeffs_prim, E)
            phi += ang_val * contracted_radial
        
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
