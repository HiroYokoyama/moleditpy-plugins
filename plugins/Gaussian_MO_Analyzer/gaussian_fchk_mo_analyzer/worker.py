import numpy as np
from PyQt6.QtCore import QThread, pyqtSignal
from .analyzer import CubeWriter

class CalcWorker(QThread):
    progress_sig = pyqtSignal(int)
    finished_sig = pyqtSignal(bool, str) # success, message/path

    def __init__(self, engine, mo_idx, n_points, margin, atoms, mo_coeffs, output_path, mode="MO"):
        super().__init__()
        self.engine = engine
        self.mo_idx = mo_idx
        self.n_points = n_points
        self.margin = margin
        self.atoms = atoms
        self.mo_coeffs = mo_coeffs
        self.output_path = output_path
        self.mode = mode
        self._is_cancelled = False

    def run(self):
        try:
            # Grid Definition
            min_c = np.min(self.atoms, axis=0) - self.margin
            max_c = np.max(self.atoms, axis=0) + self.margin
            span = max_c - min_c
            
            nx = ny = nz = self.n_points
            
            # Calculate step sizes (dx, dy, dz)
            dx = span[0] / (nx - 1)
            dy = span[1] / (ny - 1)
            dz = span[2] / (nz - 1)
            
            # Generate points
            x = np.linspace(min_c[0], max_c[0], nx)
            y = np.linspace(min_c[1], max_c[1], ny)
            z = np.linspace(min_c[2], max_c[2], nz)
            
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
            grid_points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
            
            n_total = grid_points.shape[0]
            chunk_size = 50000
            result_flat = np.zeros(n_total)
            
            for i, start in enumerate(range(0, n_total, chunk_size)):
                if self._is_cancelled: return
                
                end = min(start + chunk_size, n_total)
                chunk_pts = grid_points[start:end]
                
                # Evaluate
                if self.mode == "MO":
                    val = self.engine.evaluate_mo_on_grid(self.mo_idx, chunk_pts, self.mo_coeffs)
                else:
                    # self.mo_idx is basis_idx in this context
                    val = self.engine.evaluate_basis_on_grid(self.mo_idx, chunk_pts)
                    
                result_flat[start:end] = val
                
                # Progress
                pct = int((end / n_total) * 100)
                self.progress_sig.emit(pct)
            
            # Write
            vol_data = result_flat.reshape(nx, ny, nz)
            vectors = np.diag([dx, dy, dz])
            
            atom_nos = self.engine.fchk.get("Atomic numbers")
            
            label = f"MO {self.mo_idx+1}" if self.mode == "MO" else f"Basis {self.mo_idx+1}"
            CubeWriter.write(self.output_path, self.atoms, atom_nos, min_c, vectors, vol_data, comment=label)
            
            self.finished_sig.emit(True, self.output_path)
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.finished_sig.emit(False, str(e))
