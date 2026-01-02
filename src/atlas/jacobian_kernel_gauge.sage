import sys
import os
import random as rnd
from sage.all import *

# Path setup
sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.intrinsic_lattice import IntrinsicLattice
from src.posgeom.moment_map_laplacian import MomentMapLaplacian

# --- Utilities ---
def bracket(l1, l2):
    return l1[0]*l2[1] - l1[1]*l2[0]

def compute_s_ij(lambdas, tildes, n=6):
    s = {}
    for i in range(n):
        for j in range(i+1, n):
            li = lambdas[i]
            lj = lambdas[j]
            ti = tildes[i]
            tj = tildes[j]
            ang = bracket(li, lj)
            sq = bracket(tj, ti) 
            val = ang * sq
            s[(i,j)] = val
            s[(j,i)] = val
    return s

def solve_conservation(lambdas, tildes_free, n):
    rhs_0 = 0
    rhs_1 = 0
    for i in range(n-2):
        rhs_0 -= lambdas[i] * tildes_free[i][0]
        rhs_1 -= lambdas[i] * tildes_free[i][1]
        
    M = matrix(RR, [[lambdas[n-2][0], lambdas[n-1][0]], 
                    [lambdas[n-2][1], lambdas[n-1][1]]])
    try:
        sol_x = M.solve_right(rhs_0)
        sol_y = M.solve_right(rhs_1)
        tildes = {}
        for i in range(n-2): tildes[i] = tildes_free[i]
        tildes[n-2] = vector(RR, [sol_x[0], sol_y[0]])
        tildes[n-1] = vector(RR, [sol_x[1], sol_y[1]])
        return tildes
    except:
        return None

class CanonicalJacobianEvaluator:
    def __init__(self, n=6):
        self.n = n
        self.basis_edges = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5)]
        self.charts = {} # Cache for chart geometry
        self.kernel_refs = {} # Cache for reference kernels k1, k2
        
    def get_chart_geom(self, roots):
        roots = tuple(sorted(list(roots)))
        if roots not in self.charts:
            exponents, edge_order = get_forest_exponents(self.n, roots)
            lattice = IntrinsicLattice(exponents)
            mml = MomentMapLaplacian(self.n, roots, edge_order)
            self.charts[roots] = {
                "lattice": lattice,
                "mml": mml,
                "edge_order": edge_order
            }
        return self.charts[roots]

    def _compute_point_data(self, ts, ts_tilde_free, x_spinor, y_spinor):
        n = self.n
        lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
        tildes = solve_conservation(lambdas, ts_tilde_free, n)
        if tildes is None: return None
        
        C = {}
        for i in range(n): C[i] = bracket(lambdas[i], x_spinor) * bracket(lambdas[i], y_spinor)
        
        return lambdas, tildes, C

    def _compute_z_vals(self, lambdas, tildes, C, edge_order):
        z_vals = []
        for (u, v) in edge_order:
            ang = bracket(lambdas[u], lambdas[v])
            sq = bracket(tildes[u], tildes[v]) 
            if abs(ang) < 1e-15: ang = 1e-15
            val = (sq / ang) * C[u] * C[v]
            z_vals.append(val)
        return z_vals

    def compute_jacobian_matrices(self, roots, ts, ts_tilde_free, x_spinor, y_spinor, geom=None):
        if geom is None: geom = self.get_chart_geom(roots)
        lattice = geom['lattice']
        mml = geom['mml']
        edge_order = geom['edge_order']
        B_mat = lattice.B
        
        point_data = self._compute_point_data(ts, ts_tilde_free, x_spinor, y_spinor)
        if point_data is None: return None, None
        lambdas, tildes, C = point_data
        
        z_vals = self._compute_z_vals(lambdas, tildes, C, edge_order)
        X, H = mml.compute_X_H(z_vals)
        
        # Central t
        diff_X = X - lattice.a0
        try:
             # Try solve_right first
            t_center = B_mat.solve_right(diff_X)
        except:
            # LS fallback
            try:
                BT = B_mat.transpose()
                t_center = (BT * B_mat).solve_right(BT * diff_X)
            except:
                return None, None
                
        # Central s
        s_dict = compute_s_ij(lambdas, tildes, self.n)
        s_center = vector(RR, [s_dict[edge] for edge in self.basis_edges])
        
        # Perturbations
        delta = 1e-5
        dim_t = lattice.dim # 11
        
        grads_t = []
        grads_s = []
        
        # We need 11 independent directions to form full A and B matrices
        # Try a bit more than 11 times
        for _ in range(dim_t + 10):
            d_ts = [rnd.gauss(0, 1) * delta for _ in range(self.n)]
            d_tildes = [vector(RR, [rnd.gauss(0, 1) * delta, rnd.gauss(0, 1) * delta]) for _ in range(self.n-2)]
            
            ts_p = [ts[i] + d_ts[i] for i in range(self.n)]
            ts_tilde_free_p = {i: ts_tilde_free[i] + d_tildes[i] for i in range(self.n-2)}
            
            p_data = self._compute_point_data(ts_p, ts_tilde_free_p, x_spinor, y_spinor)
            if p_data is None: continue
            l_p, t_p_vecs, C_p = p_data
            
            z_p = self._compute_z_vals(l_p, t_p_vecs, C_p, edge_order)
            X_p, _ = mml.compute_X_H(z_p)
            
            try:
                t_p_intrinsic = B_mat.solve_right(X_p - lattice.a0)
            except:
                try:
                    diff_p = X_p - lattice.a0
                    BT = B_mat.transpose()
                    t_p_intrinsic = (BT * B_mat).solve_right(BT * diff_p)
                except: continue
                
            s_p_dict = compute_s_ij(l_p, t_p_vecs, self.n)
            s_p = vector(RR, [s_p_dict[edge] for edge in self.basis_edges])
            
            grads_t.append((t_p_intrinsic - t_center)/delta)
            grads_s.append((s_p - s_center)/delta)
            
            if len(grads_t) >= dim_t: break
            
        if len(grads_t) < dim_t: return None, None
        
        # A matrix: 11x11 (dt/du)
        A = matrix(RR, grads_t[:dim_t]).transpose()
        # B matrix: 9x11 (ds/du)
        B = matrix(RR, grads_s[:dim_t]).transpose()
        
        if abs(A.det()) < 1e-20: return None, None
        
        # J_s = B * A^-1
        J_s = B * A.inverse()
        
        return J_s, (X, H)

    def get_reference_kernels(self, roots):
        roots = tuple(sorted(list(roots)))
        if roots in self.kernel_refs:
            return self.kernel_refs[roots]
            
        # Compute reference kernels at a generic point
        # Use a fixed seed or just a generic point construction
        # To be consistent across runs, maybe hardcode or ensure stability?
        # For now, generate on fly if missing.
        
        # Generic point
        ts = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        ts_tilde_free = {
            0: vector(RR, [0.1, 0.2]),
            1: vector(RR, [0.3, 0.4]),
            2: vector(RR, [0.5, 0.6]),
            3: vector(RR, [0.7, 0.8])
        }
        x_s = vector(RR, [1, -2.0])
        y_s = vector(RR, [1, 12.0])
        
        # Adjust per chart to ensure we are inside?
        # Forest polytope is defined for positive z.
        # This generic kinematic point might not be positive for all charts.
        # But derivatives are defined locally. 
        # We need a point where z > 0 for THIS chart.
        # This is tricky. Finding an interior point is hard.
        # BUT: For the kernel definition, we just need ANY valid Jacobian J_s at some point.
        # Let's try to find a valid point for this chart.
        
        # Strategy: Sample until valid
        for _ in range(50):
            ts_rnd = [rnd.uniform(0, 10) for _ in range(self.n)]
            ts_tilde_free_rnd = {i: vector(RR, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for i in range(self.n-2)}
            
            J_s, _ = self.compute_jacobian_matrices(roots, ts_rnd, ts_tilde_free_rnd, x_s, y_s)
            if J_s is not None:
                # Compute kernel
                # J_s is 9x11
                K = J_s.right_kernel() # Returns subspace
                if K.dimension() == 2:
                    basis = K.basis()
                    self.kernel_refs[roots] = (basis[0], basis[1])
                    return (basis[0], basis[1])
                    
        raise ValueError(f"Could not find reference kernel for chart {roots}")

    def evaluate(self, roots, ts, ts_tilde_free, x_spinor, y_spinor):
        # 1. Compute Geometry
        geom = self.get_chart_geom(roots)
        
        # 2. Compute Jacobian J_s = ds/dt (9x11)
        result = self.compute_jacobian_matrices(roots, ts, ts_tilde_free, x_spinor, y_spinor, geom)
        J_s, matrices = result
        if J_s is None: return 0
        X, H = matrices
        
        # 3. Compute Omega_11
        lattice = geom['lattice']
        B_mat = lattice.B
        H_int = B_mat.transpose() * H * B_mat
        det_H = H_int.det()
        if abs(det_H) < 1e-20: return 0
        Omega_11 = 1.0 / det_H
        
        # 4. Get Reference Kernels k1, k2
        k1, k2 = self.get_reference_kernels(roots)
        
        # 5. Form J_full (11x11)
        # Top 9 rows: J_s
        # Bottom 2 rows: k1^T, k2^T
        
        # Convert J_s to list of lists
        rows = [list(row) for row in J_s.rows()]
        rows.append(list(k1))
        rows.append(list(k2))
        
        J_full = matrix(RR, rows)
        det_J = J_full.det()
        
        if abs(det_J) < 1e-20: return 0
        
        # 6. Canonical Form
        # Omega_9 = Omega_11 / det(J_full)
        return Omega_11 / det_J

def run_stability_check():
    print("Running Stability Check for Canonical Jacobian...")
    
    n = 6
    roots = (0, 1, 2)
    evaluator = CanonicalJacobianEvaluator(n)
    
    # 1. Warmup / Init Reference
    print(f"Initializing reference kernel for {roots}...")
    k1, k2 = evaluator.get_reference_kernels(roots)
    print(f"Reference Kernel Norms: {k1.norm():.4f}, {k2.norm():.4f}")
    
    # 2. Sample points and check values
    vals = []
    
    print("Sampling 20 points...")
    for i in range(20):
        ts = [rnd.uniform(0, 10) for _ in range(n)]
        ts_tilde_free = {i: vector(RR, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for i in range(n-2)}
        x_s = vector(RR, [1, -2.0])
        y_s = vector(RR, [1, 12.0])
        
        val = evaluator.evaluate(roots, ts, ts_tilde_free, x_s, y_s)
        
        # Also compute M_MHV for comparison (ratio should be somewhat stable if single chart?)
        # No, single chart doesn't match M_MHV. 
        # Just check if val is non-zero and not exploding randomly (unless near boundary)
        
        if val != 0:
            vals.append(val)
            
    print(f"Computed {len(vals)} valid values.")
    print("First 5 values:", ["%.4e" % v for v in vals[:5]])
    
    # Check sign stability
    signs = [1 if v > 0 else -1 for v in vals]
    print(f"Signs: {signs}")
    
    # With fixed kernel gauge, the sign should be relatively consistent 
    # (unless crossing a geometric boundary of the chart parametrization, 
    # but we are in a convex polytope? No, the parametrization map might flip orientation?)
    # Actually, the map t -> s is non-linear.
    
    # Key check: Do we see massive numerical noise?
    # Values should be smooth.
    
    return True

if __name__ == "__main__":
    run_stability_check()

