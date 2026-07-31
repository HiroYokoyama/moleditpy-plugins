import numpy as np
import re


class UnsupportedBasisError(ValueError):
    """The FCHK uses basis functions this plugin cannot render exactly.

    Raised instead of rendering a partial or mis-indexed orbital — the
    message is written for the user and is shown verbatim in the dialog.
    """


#: Factorials indexed by n: _FACT[n] = n!, _FACT2[n] = (2n)!
_FACT = {0: 1, 1: 1, 2: 2, 3: 6, 4: 24, 5: 120, 6: 720, 7: 5040, 8: 40320}
_FACT2 = {
    0: 1,
    1: 2,
    2: 24,
    3: 720,
    4: 40320,
    5: 3628800,
    6: 479001600,
    7: 87178291200,
    8: 20922789888000,
}


class FCHKReader:
    """
    Parses Gaussian FCHK files and holds data in a dictionary.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.data = {}
        self.parse()

    def parse(self):
        with open(self.filepath, "r") as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            # Match Array Header: "Atomic numbers   I   N=   12"
            match_array = re.match(r"^(.*)\s+([IRC])\s+N=\s+(\d+)$", line)

            # Match Scalar: "Number of electrons   I   42" or "Charge   R   0.0"
            # Note: Scalars usually don't have N=, just the value at the end.
            match_scalar = re.match(r"^(.*)\s+([IRC])\s+([-\d\.E\+]+)$", line)

            if match_array:
                tag = match_array.group(1).strip()
                dtype = match_array.group(2)
                count = int(match_array.group(3))

                raw_values = []
                i += 1

                # Read data lines until we have enough values
                while len(raw_values) < count and i < len(lines):
                    # Robust splitting for packed Fortran numbers
                    l_str = (
                        lines[i]
                        .replace("D+", "E+")
                        .replace("D-", "E-")
                        .replace("\n", "")
                        .strip()
                    )
                    l_str = l_str.replace("-", " -").replace("E -", "E-")
                    vals = l_str.split()
                    raw_values.extend(vals)
                    i += 1

                # Convert to numpy array
                if dtype == "I":
                    self.data[tag] = np.array(raw_values[:count], dtype=int)
                elif dtype == "R":
                    self.data[tag] = np.array(raw_values[:count], dtype=float)
                else:
                    self.data[tag] = raw_values[:count]

            elif match_scalar:
                tag = match_scalar.group(1).strip()
                dtype = match_scalar.group(2)
                val_str = match_scalar.group(3)

                # Store as single-element array for consistency/easy access
                if dtype == "I":
                    self.data[tag] = np.array([int(val_str)], dtype=int)
                elif dtype == "R":
                    try:
                        v = float(val_str.replace("D", "E"))
                        self.data[tag] = np.array([v], dtype=float)
                    except:
                        pass  # Ignore if parse fail

                i += 1
            else:
                i += 1

    def get(self, key, default=None):
        return self.data.get(key, default)

    def get_xyz_block(self):
        """
        Reconstructs an XYZ block string (in Angstroms) for RDKit.
        """
        atoms = self.get("Atomic numbers", None)
        coords = self.get("Current cartesian coordinates", None)

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
        num = (8 * alpha) ** L * _FACT[l] * _FACT[m] * _FACT[n]
        den = _FACT2[l] * _FACT2[m] * _FACT2[n]

        return ((2 * alpha / np.pi) ** 0.75) * np.sqrt(num / den)

    def _n_cart_part(self, l, m, n):
        """Angular part of the Cartesian Gaussian normalization."""
        return np.sqrt(
            (_FACT[l] * _FACT[m] * _FACT[n]) / (_FACT2[l] * _FACT2[m] * _FACT2[n])
        )

    @staticmethod
    def _odd_double_factorial(k):
        """(k-1)!! for even k — the 1-D Gaussian moment integral factor."""
        result = 1.0
        while k > 1:
            result *= k - 1
            k -= 2
        return result

    def _angular_overlap(self, a, b):
        """Overlap of two normalized Cartesian Gaussians in the same shell.

        Scaled so a single component overlaps itself as 1.0. The radial part
        is shared by every component of a shell, so it cancels and the result
        does not depend on the exponent.
        """
        moment = 1.0
        for ka, kb in zip(a, b):
            k = ka + kb
            if k % 2:
                return 0.0  # odd power integrates to zero over all space
            moment *= self._odd_double_factorial(k)
        return (2.0 ** sum(a)) * moment * self._n_cart_part(*a) * self._n_cart_part(*b)

    def _build_sph_defs(self, polys):
        """Turn real-solid-harmonic polynomials into component weights.

        Each component is evaluated as N_cart(l,m,n) * x^l y^m z^n, and
        N_cart differs between topologies (N(0,0,3) is N(2,0,1)/sqrt(5)).
        The previous hand-typed table applied one scalar per component, which
        squashed any component mixing topologies — it put the f(z^3) nodal
        cone at 28.6 degrees instead of 39.2. Dividing out N_cart and then
        normalizing against the exact angular overlap fixes the shape and the
        scale together, and matches the s/p/d shells.
        """
        defs = []
        for poly in polys:
            weights = [(c_poly / self._n_cart_part(*lmn), lmn) for c_poly, lmn in poly]
            norm_sq = sum(
                wa * wb * self._angular_overlap(a, b)
                for wa, a in weights
                for wb, b in weights
            )
            scale = 1.0 / np.sqrt(norm_sq) if norm_sq > 0 else 1.0
            defs.append([(w * scale, lmn) for w, lmn in weights])
        return defs

    #: Human-readable names for the shell types Gaussian writes.
    _SHELL_NAMES = {
        0: "S",
        1: "P",
        -1: "SP",
        2: "6D (Cartesian)",
        -2: "5D (pure)",
        3: "10F (Cartesian)",
        -3: "7F (pure)",
        4: "15G (Cartesian)",
        -4: "9G (pure)",
        5: "21H (Cartesian)",
        -5: "11H (pure)",
    }

    def _reject_unsupported_shells(self, shell_types, d_is_cart, f_is_cart):
        """Refuse the whole file if any shell type cannot be rendered exactly.

        Rendering only the shells we understand is not an option: the FCHK MO
        coefficient block covers every basis function, so an unrenderable
        shell either desyncs the coefficient indexing or silently drops part
        of the orbital. Either way the picture is wrong with nothing on screen
        to say so. Better to refuse the file and name the shell type.
        """
        offenders = []
        for stype in shell_types:
            stype = int(stype)
            expected = self._shell_n_functions(stype, d_is_cart, f_is_cart)
            defs = self.basis_definitions.get(abs(stype) if stype < -1 else stype, None)
            if defs is None or len(defs) != expected:
                name = self._SHELL_NAMES.get(stype, f"type {stype}")
                if name not in offenders:
                    offenders.append(name)

        if offenders:
            raise UnsupportedBasisError(
                "This FCHK uses basis functions the analyzer cannot render "
                "exactly: " + ", ".join(offenders) + ".\n\n"
                "Visualization is blocked rather than showing an orbital that "
                "would be silently incorrect.\n\n"
                "Re-run the Gaussian job with spherical harmonics (5D 7F, the "
                "default for most modern basis sets) to produce a file this "
                "plugin can render."
            )

    @staticmethod
    def _shell_n_functions(stype, d_is_cart=True, f_is_cart=True):
        """Number of basis functions a Gaussian shell type contributes.

        Purity comes from the file's "Pure/Cartesian" flags rather than the
        sign of the shell type, matching how the definitions above are
        chosen. Gaussian writes the two consistently, but the flags are what
        this code acts on, so the count has to agree with them or a
        perfectly renderable file would be refused.

        Gaussian's "Pure/Cartesian f shells" flag governs f and every higher
        shell; S and P are the same either way. -1 is the SP special case.
        """
        if stype == -1:
            return 4
        total_l = abs(stype)
        if total_l <= 1:
            return 2 * total_l + 1
        is_cart = d_is_cart if total_l == 2 else f_is_cart
        if is_cart:
            return (total_l + 1) * (total_l + 2) // 2
        return 2 * total_l + 1

    def _prepare_basis_set(self):
        shell_types = self.fchk.get("Shell types", None)
        prim_per_shell = self.fchk.get("Number of primitives per shell", None)
        shell_to_atom = self.fchk.get("Shell to atom map", None)
        prim_exps = self.fchk.get("Primitive exponents", None)
        cont_coeffs = self.fchk.get("Contraction coefficients", None)

        # Pure/Cartesian Flags (0=Pure/Spherical, 1=Cartesian)
        # Default to 1 (Cartesian) if not found, to match old behavior
        d_is_cart = True
        f_is_cart = True

        flag_d = self.fchk.get("Pure/Cartesian d shells", None)
        if flag_d is not None and len(flag_d) > 0 and flag_d[0] == 0:
            d_is_cart = False

        flag_f = self.fchk.get("Pure/Cartesian f shells", None)
        if flag_f is not None and len(flag_f) > 0 and flag_f[0] == 0:
            f_is_cart = False

        # Build Basis Definitions (Linear Combinations)
        # Structure: Type -> [  List of Basis Functions in this Shell ]
        # Each Basis Function -> List of Components (coeff, (l,m,n))

        # Cartesian Definitions (Simple 1:1 map)
        cart_s = [[(1.0, (0, 0, 0))]]  # S
        cart_p = [
            [(1.0, (1, 0, 0))],
            [(1.0, (0, 1, 0))],
            [(1.0, (0, 0, 1))],
        ]  # Px, Py, Pz
        cart_d = [
            [(1.0, (2, 0, 0))],
            [(1.0, (0, 2, 0))],
            [(1.0, (0, 0, 2))],
            [(1.0, (1, 1, 0))],
            [(1.0, (1, 0, 1))],
            [(1.0, (0, 1, 1))],
        ]  # XX, YY, ZZ, XY, XZ, YZ
        cart_f = [
            [(1.0, (3, 0, 0))],
            [(1.0, (0, 3, 0))],
            [(1.0, (0, 0, 3))],
            [(1.0, (1, 2, 0))],
            [(1.0, (2, 1, 0))],
            [(1.0, (2, 0, 1))],
            [(1.0, (1, 0, 2))],
            [(1.0, (0, 1, 2))],
            [(1.0, (0, 2, 1))],
            [(1.0, (1, 1, 1))],
        ]  # XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ

        self.basis_definitions = {
            0: cart_s,
            1: cart_p,
            -1: cart_s + cart_p,  # SP
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
                [
                    (-0.5, (2, 0, 0)),
                    (-0.5, (0, 2, 0)),
                    (1.0, (0, 0, 2)),
                ],  # d(3z^2-r^2) -> 0.5*(2z^2-x^2-y^2) gives -0.5x^2... Correct.
                [(1.0, (1, 0, 1))],  # d(xz)
                [(1.0, (0, 1, 1))],  # d(yz)
                [
                    (0.8660254037844387, (2, 0, 0)),
                    (-0.8660254037844387, (0, 2, 0)),
                ],  # d(x^2-y^2)
                [(1.0, (1, 1, 0))],  # d(xy)
            ]
            self.basis_definitions[2] = sph_d

        # Handle F Shells
        if f_is_cart:
            self.basis_definitions[3] = cart_f
        else:
            # Spherical F (7 components), Gaussian m order 0,+1,-1,+2,-2,+3,-3.
            # Polynomials are the real solid harmonics r^3 * Y_3m in Cartesian
            # monomials; _build_sph_defs supplies the normalization.
            f_polys = [
                # f(0): z(5z^2-3r^2) -> 2ZZZ - 3XXZ - 3YYZ
                [(2.0, (0, 0, 3)), (-3.0, (2, 0, 1)), (-3.0, (0, 2, 1))],
                # f(+1): x(5z^2-r^2) -> 4XZZ - XXX - XYY
                [(4.0, (1, 0, 2)), (-1.0, (3, 0, 0)), (-1.0, (1, 2, 0))],
                # f(-1): y(5z^2-r^2) -> 4YZZ - XXY - YYY
                [(4.0, (0, 1, 2)), (-1.0, (2, 1, 0)), (-1.0, (0, 3, 0))],
                # f(+2): z(x^2-y^2)
                [(1.0, (2, 0, 1)), (-1.0, (0, 2, 1))],
                # f(-2): xyz
                [(1.0, (1, 1, 1))],
                # f(+3): x(x^2-3y^2)
                [(1.0, (3, 0, 0)), (-3.0, (1, 2, 0))],
                # f(-3): y(3x^2-y^2)
                [(3.0, (2, 1, 0)), (-1.0, (0, 3, 0))],
            ]
            self.basis_definitions[3] = self._build_sph_defs(f_polys)

        # Spherical G (9 components), same m ordering. Gaussian's
        # "Pure/Cartesian f shells" flag governs f and every higher shell, so
        # g purity follows f_is_cart. Cartesian 15G is NOT defined here: its
        # 15-monomial ordering cannot be confirmed without a sample file, and
        # a guessed ordering would render a confidently wrong orbital. Such
        # files are rejected outright by _reject_unsupported_shells.
        #
        # Before this, a g shell was skipped without reserving its
        # coefficient slots, which shifted every later shell's start_idx AND
        # shrank n_basis -- so mo_coeffs_all[mo_idx * n_basis : ...] strided
        # into the wrong place and every orbital past the first was garbage.
        g_polys = [
            # g(0): (35z^4 - 30z^2 r^2 + 3r^4)/8
            [
                (1.0, (0, 0, 4)),
                (-3.0, (2, 0, 2)),
                (-3.0, (0, 2, 2)),
                (0.375, (4, 0, 0)),
                (0.375, (0, 4, 0)),
                (0.75, (2, 2, 0)),
            ],
            # g(+1): xz(7z^2 - 3r^2)
            [(4.0, (1, 0, 3)), (-3.0, (3, 0, 1)), (-3.0, (1, 2, 1))],
            # g(-1): yz(7z^2 - 3r^2)
            [(4.0, (0, 1, 3)), (-3.0, (0, 3, 1)), (-3.0, (2, 1, 1))],
            # g(+2): (x^2-y^2)(7z^2 - r^2)
            [(6.0, (2, 0, 2)), (-6.0, (0, 2, 2)), (-1.0, (4, 0, 0)), (1.0, (0, 4, 0))],
            # g(-2): xy(7z^2 - r^2)
            [(6.0, (1, 1, 2)), (-1.0, (3, 1, 0)), (-1.0, (1, 3, 0))],
            # g(+3): xz(x^2 - 3y^2)
            [(1.0, (3, 0, 1)), (-3.0, (1, 2, 1))],
            # g(-3): yz(3x^2 - y^2)
            [(3.0, (2, 1, 1)), (-1.0, (0, 3, 1))],
            # g(+4): x^4 - 6x^2y^2 + y^4
            [(1.0, (4, 0, 0)), (-6.0, (2, 2, 0)), (1.0, (0, 4, 0))],
            # g(-4): xy(x^2 - y^2)
            [(4.0, (3, 1, 0)), (-4.0, (1, 3, 0))],
        ]
        if not f_is_cart:
            self.basis_definitions[4] = self._build_sph_defs(g_polys)

        self._reject_unsupported_shells(shell_types, d_is_cart, f_is_cart)

        # Check for P-modifiers (used for SP shells P-component)
        # Try both tag variations
        p_modifiers = self.fchk.get("P(S=P) Contraction coefficients", None)
        if p_modifiers is None:
            p_modifiers = self.fchk.get("P(S=P) modifiers", None)

        coords = self.fchk.get("Current cartesian coordinates").reshape(-1, 3)  # Bohr

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

            if stype == -1:  # SP Shell
                if p_modifiers is not None and len(p_modifiers) > 0:
                    coeffs_p = p_modifiers[exp_ptr : exp_ptr + n_prim]
                    coeff_ptr += n_prim
                else:
                    coeffs_p = cont_coeffs[coeff_ptr + n_prim : coeff_ptr + 2 * n_prim]
                    coeff_ptr += 2 * n_prim

                current_shells_to_add.append(
                    {
                        "type": 0,
                        "center": atom_center,
                        "exps": exps,
                        "coeffs": coeffs_s,
                        "defs": self.basis_definitions[0],
                    }
                )
                current_shells_to_add.append(
                    {
                        "type": 1,
                        "center": atom_center,
                        "exps": exps,
                        "coeffs": coeffs_p,
                        "defs": self.basis_definitions[1],
                    }
                )

                exp_ptr += n_prim

            else:  # S, P, D, F... AND -2(D), -3(F) handling
                # Handle negative types -2, -3 etc which are just D, F...
                effective_type = stype
                if effective_type < -1:
                    effective_type = abs(effective_type)

                # _reject_unsupported_shells has already refused any type
                # without an exactly-sized definition, so this lookup cannot
                # come back None. Skipping a shell here would desync every
                # later shell's coefficient indexing.
                defs_list = self.basis_definitions[effective_type]

                current_shells_to_add.append(
                    {
                        "type": effective_type,
                        "center": atom_center,
                        "exps": exps,
                        "coeffs": coeffs_s,
                        "defs": defs_list,
                    }
                )

                exp_ptr += n_prim
                coeff_ptr += n_prim

            for sh in current_shells_to_add:
                func_norm_coeffs = []

                # Precompute coefficients for efficient evaluation
                for basis_func_def in sh["defs"]:
                    comps_data = []
                    for weight, (l, m, n) in basis_func_def:
                        # Calculate prefactors for this l,m,n for all exponents
                        prim_norms = np.array(
                            [
                                self._normalization_prefactor(a, l, m, n)
                                for a in sh["exps"]
                            ]
                        )

                        # Apply Contraction Coeffs & Weights
                        comp_coeffs = sh["coeffs"] * prim_norms * weight

                        comps_data.append(
                            {"l": l, "m": m, "n": n, "coeffs": comp_coeffs}
                        )

                    func_norm_coeffs.append(comps_data)

                sh["basis_data"] = func_norm_coeffs
                sh["start_idx"] = basis_idx_counter
                basis_idx_counter += len(sh["defs"])

                self.shells.append(sh)

        self.n_basis = basis_idx_counter

    def evaluate_mo_on_grid(
        self,
        mo_idx: int,
        grid_coords: np.ndarray,
        mo_coeffs_all: np.ndarray,
        progress_callback=None,
    ):
        start = mo_idx * self.n_basis
        end = start + self.n_basis
        if end > len(mo_coeffs_all):
            raise ValueError(
                f"MO index {mo_idx} out of range (need {end} coeffs, have {len(mo_coeffs_all)})"
            )

        my_coeffs = mo_coeffs_all[start:end]
        phi_mo = np.zeros(grid_coords.shape[0])

        for i, sh in enumerate(self.shells):
            r_vec = grid_coords - sh["center"]
            r2 = np.sum(r_vec**2, axis=1)
            E = np.exp(-np.outer(sh["exps"], r2))

            current_basis_idx = sh["start_idx"]

            # Iterate over Basis Functions in this Shell
            for basis_func_data in sh["basis_data"]:
                c_mo = my_coeffs[current_basis_idx]
                current_basis_idx += 1

                if abs(c_mo) < 1e-8:
                    continue

                val_accum = np.zeros(grid_coords.shape[0])

                for comp in basis_func_data:
                    l, m, n = comp["l"], comp["m"], comp["n"]
                    coeffs_prim = comp["coeffs"]

                    ang_val = np.ones(grid_coords.shape[0])
                    if l > 0:
                        ang_val *= r_vec[:, 0] ** l
                    if m > 0:
                        ang_val *= r_vec[:, 1] ** m
                    if n > 0:
                        ang_val *= r_vec[:, 2] ** n

                    contracted_radial = np.dot(coeffs_prim, E)
                    val_accum += ang_val * contracted_radial

                phi_mo += c_mo * val_accum

        return phi_mo

    def get_basis_labels(self):
        """
        Returns a list of labels for the basis functions.
        """
        labels = []
        atom_list = self.fchk.get("Atomic numbers", None)
        shell_to_atom = self.fchk.get("Shell to atom map", None)
        shell_types = self.fchk.get("Shell types", None)

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
            if stype == -1:  # SP
                comps = ["S", "Px", "Py", "Pz"]
            else:
                # Negative types below -1 (e.g. -2, -3) denote pure/spherical
                # D, F, ... shells; mirror the abs() handling used in
                # _prepare_basis_set so label counts match n_basis.
                effective_type = abs(stype) if stype < -1 else stype
                count = len(self.basis_definitions.get(effective_type, []))
                comps = [f"L{effective_type}_{k}" for k in range(count)]

            for c in comps:
                labels.append(f"{atom_idx + 1} {sym} {c}")

        return labels

    def evaluate_basis_on_grid(
        self, basis_idx: int, grid_coords: np.ndarray, progress_callback=None
    ):
        """
        Evaluates a single basis function (AO) on the grid.
        """
        phi = np.zeros(grid_coords.shape[0])

        current_idx = 0
        target_shell = None
        target_func_idx = 0

        for sh in self.shells:
            n_funcs = len(sh["basis_data"])
            if current_idx + n_funcs > basis_idx:
                target_shell = sh
                target_func_idx = basis_idx - current_idx
                break
            current_idx += n_funcs

        if target_shell is None:
            return phi

        sh = target_shell
        basis_func_data = sh["basis_data"][target_func_idx]

        r_vec = grid_coords - sh["center"]
        r2 = np.sum(r_vec**2, axis=1)
        E = np.exp(-np.outer(sh["exps"], r2))

        for comp in basis_func_data:
            l, m, n = comp["l"], comp["m"], comp["n"]
            coeffs_prim = comp["coeffs"]

            ang_val = np.ones(grid_coords.shape[0])
            if l > 0:
                ang_val *= r_vec[:, 0] ** l
            if m > 0:
                ang_val *= r_vec[:, 1] ** m
            if n > 0:
                ang_val *= r_vec[:, 2] ** n

            contracted_radial = np.dot(coeffs_prim, E)
            phi += ang_val * contracted_radial

        return phi


class CubeWriter:
    @staticmethod
    def write(filepath, atoms, atom_nos, origin, vectors, data, comment=""):
        nx, ny, nz = data.shape
        n_atoms = len(atoms)

        with open(filepath, "w") as f:
            f.write(f"Moleditpy Cube File: {comment}\n")
            f.write("Generated by FCHK-to-CUBE Plugin\n")
            f.write(
                f"{n_atoms:5d} {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n"
            )
            f.write(
                f"{nx:5d} {vectors[0, 0]:12.6f} {vectors[0, 1]:12.6f} {vectors[0, 2]:12.6f}\n"
            )
            f.write(
                f"{ny:5d} {vectors[1, 0]:12.6f} {vectors[1, 1]:12.6f} {vectors[1, 2]:12.6f}\n"
            )
            f.write(
                f"{nz:5d} {vectors[2, 0]:12.6f} {vectors[2, 1]:12.6f} {vectors[2, 2]:12.6f}\n"
            )

            for i in range(n_atoms):
                z_no = atom_nos[i]
                f.write(
                    f"{z_no:5d} {float(z_no):12.6f} {atoms[i][0]:12.6f} {atoms[i][1]:12.6f} {atoms[i][2]:12.6f}\n"
                )

            vals = data.flatten()
            for i, v in enumerate(vals):
                f.write(f"{v:13.5E}")
                if (i + 1) % 6 == 0:
                    f.write("\n")
            if len(vals) % 6 != 0:
                f.write("\n")
