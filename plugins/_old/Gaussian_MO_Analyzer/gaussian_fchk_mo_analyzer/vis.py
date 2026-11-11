from itertools import islice
import os
import numpy as np
import pyvista as pv

# --- Constants ---
BOHR_TO_ANG = 0.529177249


class CubeVisualizer:
    def __init__(self, mw):
        self.mw = mw
        self.plotter = mw.plotter
        self.current_grid = None
        self.actors = {}

    def load_file(self, filename):
        try:
            # Parse using simple internal parser or pyvista if robust
            # We use internal to ensure consistency with our writer
            meta = self._parse_cube(filename)
            self.current_grid = self._build_grid(meta)
            return True
        except Exception as e:
            print(f"Error loading cube: {e}")
            return False

    def _parse_cube(self, filename):
        """Parse Cube records using the PySCF Calculator's strict parser logic.

        Kept self-contained because each viewer is installed independently.
        """
        def number(token):
            return float(token.replace("D", "E").replace("d", "e"))
        if not filename or not os.path.exists(filename):
            raise FileNotFoundError(f"File not found: {filename}")
        if os.path.getsize(filename) > 256 * 1024 * 1024:
            raise ValueError("Cube file exceeds the 256 MiB loading limit.")

        with open(filename, "r", encoding="utf-8", errors="strict") as f:
            lines = f.readlines()

        if len(lines) < 6:
            raise ValueError("File too short to be a Cube file.")

        # --- Header Parsing ---
        try:
            # Line 3: Natoms, Origin
            tokens = lines[2].split()
            if len(tokens) < 4:
                raise ValueError("Invalid Origin line")
            n_atoms_raw = int(tokens[0])
            n_atoms = abs(n_atoms_raw)
            origin_raw = np.array([number(tokens[1]), number(tokens[2]), number(tokens[3])])

            def parse_vec(line):
                t = line.split()
                if len(t) < 4:
                    raise ValueError("Invalid Axis line")
                return int(t[0]), np.array([number(t[1]), number(t[2]), number(t[3])])

            nx, x_vec_raw = parse_vec(lines[3])
            ny, y_vec_raw = parse_vec(lines[4])
            nz, z_vec_raw = parse_vec(lines[5])

            is_angstrom_header = nx < 0 or ny < 0 or nz < 0
            nx, ny, nz = abs(nx), abs(ny), abs(nz)
            if not all(np.all(np.isfinite(v)) for v in (origin_raw, x_vec_raw, y_vec_raw, z_vec_raw)):
                raise ValueError("Cube header contains nonfinite coordinates.")
            if np.linalg.matrix_rank(np.array([x_vec_raw, y_vec_raw, z_vec_raw])) != 3:
                raise ValueError("Cube axis vectors must span three dimensions.")
            if min(nx, ny, nz) == 0 or nx * ny * nz > 2_000_000:
                raise ValueError("Cube grid must contain 1 to 2,000,000 points.")

        except (ValueError, IndexError) as e:
            raise ValueError(f"Header parsing failed: {e}") from e

        # --- Atoms Parsing ---
        atoms = []
        current_line = 6

        for _ in range(n_atoms):
            if current_line >= len(lines):
                raise ValueError("Incomplete cube atom block")
            line = lines[current_line].split()
            current_line += 1

            try:
                if len(line) < 5:
                    # Potentially empty line or malformed
                    raise ValueError("Malformed cube atom line")
                atomic_num = int(line[0])
                if not 0 <= atomic_num <= 118:
                    raise ValueError("Invalid cube atomic number")
                x, y, z = number(line[2]), number(line[3]), number(line[4])
                if not np.all(np.isfinite([number(line[1]), x, y, z])):
                    raise ValueError("Nonfinite cube atom coordinates")
                atoms.append((atomic_num, np.array([x, y, z])))
            except (ValueError, IndexError) as exc:
                raise ValueError("Malformed cube atom line") from exc

        # A negative atom count means the cube holds several data sets. The
        # DSET_IDS block (count + that many ids, possibly wrapped) sits *after*
        # the atom lines -- skipping a line before them consumed the first atom.
        n_datasets = int(tokens[4]) if len(tokens) > 4 else 1
        if n_atoms_raw < 0:
            if n_datasets != 1:
                raise ValueError("Negative NATOMS requires NVAL absent or equal to one.")
            try:
                parts = lines[current_line].split()
                n_datasets = int(parts[0])
                if n_datasets < 1:
                    raise ValueError("Cube dataset count must be positive")
                identifiers = [int(value) for value in parts[1:]]
                current_line += 1
                while len(identifiers) < n_datasets and current_line < len(lines):
                    identifiers.extend(int(value) for value in lines[current_line].split())
                    current_line += 1
                if len(identifiers) != n_datasets:
                    raise ValueError("Incomplete cube dataset identifiers")
            except (ValueError, IndexError) as exc:
                raise ValueError("Malformed cube dataset identifiers") from exc
        if n_datasets < 1 or nx * ny * nz * n_datasets > 2_000_000:
            raise ValueError("Cube dataset size exceeds the loading limit.")

        # --- Volumetric Data Parsing ---
        # Find start of data
        while current_line < len(lines):
            line_content = lines[current_line].strip()
            parts = line_content.split()
            if not parts:
                current_line += 1
                continue

            break

        if current_line >= len(lines):
            raise ValueError("Cube contains no volumetric data")
        else:
            # Validate every token; whitespace wrapping is independent of grid shape.
            n_points = nx * ny * nz
            expected_size = n_points * n_datasets
            tokens = (token for line in lines[current_line:] for token in line.split())
            try:
                data_values = np.fromiter((number(token) for token in islice(tokens, expected_size + 1)), dtype=float)
            except ValueError as exc:
                raise ValueError("Malformed cube volumetric data") from exc

        actual_size = len(data_values)

        if actual_size != expected_size:
            raise ValueError(f"Cube data size mismatch: expected {expected_size}, found {actual_size}" + (" or more" if actual_size > expected_size else ""))
        if not np.all(np.isfinite(data_values)):
            raise ValueError("Cube contains nonfinite volumetric values")

        if n_datasets > 1:
            # Values are interleaved point by point; the first n_points of the
            # mixed stream belong to no single orbital.
            data_values = data_values[0 : n_points * n_datasets : n_datasets]

        return {
            "atoms": atoms,
            "origin": origin_raw,
            "vectors": (x_vec_raw, y_vec_raw, z_vec_raw),
            "dims": (nx, ny, nz),
            "data": data_values,
            "is_angstrom": is_angstrom_header,
        }

    def _build_grid(self, meta):
        nx, ny, nz = meta["dims"]
        origin = meta["origin"]
        xv, yv, zv = meta["vectors"]

        grid = pv.StructuredGrid()

        # Create coordinates
        x = np.arange(nx)
        y = np.arange(ny)
        z = np.arange(nz)
        gx, gy, gz = np.meshgrid(x, y, z, indexing="ij")

        # Flatten (X-fastest for PyVista StructuredGrid points)
        gx = gx.flatten(order="F")
        gy = gy.flatten(order="F")
        gz = gz.flatten(order="F")

        # A negative voxel count on the NX line means the grid is already in
        # Angstrom (Gaussian cube convention). Converting those as if they
        # were Bohr inflated the whole isosurface by 1/0.529 = 1.89x.
        scale = 1.0 if meta.get("is_angstrom") else BOHR_TO_ANG

        points = origin + np.outer(gx, xv) + np.outer(gy, yv) + np.outer(gz, zv)
        points *= scale

        grid.points = points
        grid.dimensions = [nx, ny, nz]

        # Add data
        # Cube file data is Z-fastest (C-order flatten of (nx,ny,nz))
        # Points are X-fastest (F-order)
        raw_data = meta["data"]
        # Reshape to 3D using Input order (C)
        vol_3d = raw_data.reshape((nx, ny, nz), order="C")
        # Flatten using Grid Point order (F)
        grid.point_data["values"] = vol_3d.flatten(order="F")

        return grid

    def show_iso(
        self,
        isovalue=0.02,
        color_p="red",
        color_n="blue",
        opacity=0.5,
        smooth_shading=True,
    ):
        if not self.current_grid:
            return

        self.plotter.remove_actor("fchk_iso_p")
        self.plotter.remove_actor("fchk_iso_n")

        try:
            iso_p = self.current_grid.contour([isovalue], scalars="values")
            if iso_p.n_points > 0:
                self.plotter.add_mesh(
                    iso_p,
                    color=color_p,
                    opacity=opacity,
                    name="fchk_iso_p",
                    smooth_shading=smooth_shading,
                )

            iso_n = self.current_grid.contour([-isovalue], scalars="values")
            if iso_n.n_points > 0:
                self.plotter.add_mesh(
                    iso_n,
                    color=color_n,
                    opacity=opacity,
                    name="fchk_iso_n",
                    smooth_shading=smooth_shading,
                )

            self.plotter.render()
        except Exception as e:
            print(f"Iso error: {e}")
