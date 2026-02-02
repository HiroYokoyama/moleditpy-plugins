import numpy as np
import pyvista as pv
from PyQt6.QtGui import QColor

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
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Skip comments
        tokens = lines[2].split()
        n_atoms = abs(int(tokens[0]))
        origin = np.array([float(x) for x in tokens[1:4]])
        
        nx, *x_vec = lines[3].split()
        nx = int(nx)
        x_vec = np.array([float(x) for x in x_vec])
        
        ny, *y_vec = lines[4].split()
        ny = int(ny)
        y_vec = np.array([float(x) for x in y_vec])
        
        nz, *z_vec = lines[5].split()
        nz = int(nz)
        z_vec = np.array([float(x) for x in z_vec])
        
        # Data starts after atoms
        start_line = 6 + n_atoms
        full_str = " ".join(lines[start_line:])
        data = np.fromstring(full_str, sep=' ')
        
        return {
            "dims": (abs(nx), abs(ny), abs(nz)),
            "origin": origin,
            "vectors": (x_vec, y_vec, z_vec),
            "data": data,
            "is_angstrom": (nx < 0) 
        }

    def _build_grid(self, meta):
        nx, ny, nz = meta['dims']
        origin = meta['origin']
        xv, yv, zv = meta['vectors']
        
        grid = pv.StructuredGrid()
        
        # Create coordinates
        x = np.arange(nx)
        y = np.arange(ny)
        z = np.arange(nz)
        gx, gy, gz = np.meshgrid(x, y, z, indexing='ij')
        
        # Flatten (X-fastest for PyVista StructuredGrid points)
        gx = gx.flatten(order='F')
        gy = gy.flatten(order='F')
        gz = gz.flatten(order='F')
        
        scale = BOHR_TO_ANG
        
        points = origin + np.outer(gx, xv) + np.outer(gy, yv) + np.outer(gz, zv)
        points *= scale # Convert Bohr -> Angstrom
        
        grid.points = points
        grid.dimensions = [nx, ny, nz]
        
        # Add data
        # Cube file data is Z-fastest (C-order flatten of (nx,ny,nz))
        # Points are X-fastest (F-order)
        raw_data = meta['data']
        # Reshape to 3D using Input order (C)
        vol_3d = raw_data.reshape((nx, ny, nz), order='C')
        # Flatten using Grid Point order (F)
        grid.point_data["values"] = vol_3d.flatten(order='F')
        
        return grid

    def show_iso(self, isovalue=0.02, color_p='red', color_n='blue', opacity=0.5, smooth_shading=True):
        if not self.current_grid: return
        
        self.plotter.remove_actor("fchk_iso_p")
        self.plotter.remove_actor("fchk_iso_n")
        
        try:
            iso_p = self.current_grid.contour([isovalue], scalars="values")
            if iso_p.n_points > 0:
                self.plotter.add_mesh(iso_p, color=color_p, opacity=opacity, name="fchk_iso_p", smooth_shading=smooth_shading)
                
            iso_n = self.current_grid.contour([-isovalue], scalars="values")
            if iso_n.n_points > 0:
                self.plotter.add_mesh(iso_n, color=color_n, opacity=opacity, name="fchk_iso_n", smooth_shading=smooth_shading)
                
            self.plotter.render()
        except Exception as e:
            print(f"Iso error: {e}")
