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

def solve_conservation_generic(lambdas, tildes_free, n, Field=RR):
    rhs_0 = vector(Field, [0, 0])
    rhs_1 = vector(Field, [0, 0])
    for i in range(n-2):
        rhs_0 -= lambdas[i] * tildes_free[i][0]
        rhs_1 -= lambdas[i] * tildes_free[i][1]
        
    M = matrix(Field, [[lambdas[n-2][0], lambdas[n-1][0]], 
                       [lambdas[n-2][1], lambdas[n-1][1]]])
    try:
        sol_x = M.solve_right(rhs_0)
        sol_y = M.solve_right(rhs_1)
        tildes = {}
        for i in range(n-2): tildes[i] = tildes_free[i]
        tildes[n-2] = vector(Field, [sol_x[0], sol_y[0]])
        tildes[n-1] = vector(Field, [sol_x[1], sol_y[1]])
        return tildes
    except:
        return None

class FiberGaugeEvaluator:
    def __init__(self, n=6):
        self.n = n
        self.basis_edges = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5)]
        self.charts = {} # Cache for chart geometry and inequalities
        
    def get_chart_data(self, roots):
        roots = tuple(sorted(list(roots)))
        if roots not in self.charts:
            exponents, edge_order = get_forest_exponents(self.n, roots)
            lattice = IntrinsicLattice(exponents)
            mml = MomentMapLaplacian(self.n, roots, edge_order)
            
            # Precompute inequalities: A*X + b >= 0
            P = Polyhedron(vertices=exponents)
            ieqs = P.inequality_generator()
            inequalities = []
            for ieq in ieqs:
                b = ieq[0]
                A = vector(QQ, ieq[1:]) # Inequalities are usually rational/integer
                inequalities.append((b, A))
                
            self.charts[roots] = {
                "lattice": lattice,
                "mml": mml,
                "edge_order": edge_order,
                "inequalities": inequalities
            }
        return self.charts[roots]
    
    def is_in_polytope(self, roots, X, Field=RR):
        data = self.get_chart_data(roots)
        inequalities = data["inequalities"]
        
        tol = Field(1e-9)
        # X is in Field, A/b are QQ/ZZ. Coercion should work.
        for b, A in inequalities:
            # Explicitly cast A to Field vector to be safe
            A_f = vector(Field, list(A))
            b_f = Field(b)
            val = A_f.dot_product(X) + b_f
            if val < -tol:
                return False
        return True

    def _compute_point_data_generic(self, ts, ts_tilde_free, x_spinor, y_spinor, Field):
        n = self.n
        lambdas = {i: vector(Field, [1, ts[i]]) for i in range(n)}
        tildes = solve_conservation_generic(lambdas, ts_tilde_free, n, Field)
        if tildes is None: return None, None, None
        
        C = {}
        for i in range(n): C[i] = bracket(lambdas[i], x_spinor) * bracket(lambdas[i], y_spinor)
        
        return lambdas, tildes, C

    def _compute_z_vals(self, lambdas, tildes, C, edge_order, Field):
        z_vals = []
        for (u, v) in edge_order:
            ang = bracket(lambdas[u], lambdas[v])
            sq = bracket(tildes[v], tildes[u]) 
            # Avoid division by zero with safe threshold
            if abs(ang) < Field(1e-25): ang = Field(1e-25)
            val = (sq / ang) * C[u] * C[v]
            z_vals.append(val)
        return z_vals

    def evaluate_chart_fiber_gauge(self, roots, ts, ts_tilde_free, x_spinor, y_spinor, prec=200, require_polytope=True):
        RF = RealField(prec)
        
        # Helper to convert inputs to RF
        def to_rf(obj):
            if hasattr(obj, 'apply_map'): # Vectors, Matrices
                return obj.apply_map(lambda x: RF(x))
            elif isinstance(obj, (list, tuple)):
                return [to_rf(x) for x in obj]
            elif isinstance(obj, dict):
                return {k: to_rf(v) for k, v in obj.items()}
            return RF(obj)
            
        ts = to_rf(ts)
        ts_tilde_free = to_rf(ts_tilde_free)
        x_spinor = to_rf(x_spinor)
        y_spinor = to_rf(y_spinor)
        
        # 1. Get Chart Data
        data = self.get_chart_data(roots)
        lattice = data['lattice']
        mml = data['mml']
        edge_order = data['edge_order']
        # Change ring of B to RF
        B_mat = lattice.B.change_ring(RF)
        a0 = lattice.a0.change_ring(RF)
        
        # 2. Compute Center Point
        lambdas, tildes, C = self._compute_point_data_generic(ts, ts_tilde_free, x_spinor, y_spinor, RF)
        if lambdas is None: return RF(0)
        
        z_vals = self._compute_z_vals(lambdas, tildes, C, edge_order, RF)
        # mml uses QQ internally but compute_X_H handles symbolic/other types via arithmetic
        # but matrix construction inside might default to QQ if initialized so.
        # Let's check mml.compute_X_H.
        # It creates Matrix(QQ, ...) in the code I read.
        # I need to patch that or assume it works via coercion?
        # If I pass RF z_vals, it might fail if it tries to put RF into QQ matrix.
        # It initializes M = Matrix(QQ, ...)
        
        # For now, let's assume z_vals coerce or fail.
        # If it fails, I might need to fix MML.
        # The MML code had `M = Matrix(QQ, ...)`
        # This will fail with RealNumber(200).
        # I will need to use a local updated compute_X_H or modify MML.
        # Since MML is imported, I should rely on MML supporting generic types or I should reimplement the relevant part here.
        # Given the plan, I should probably make MML robust or override.
        # I'll override the computation here to be safe and precise.
        
        # Re-implementation of compute_X_H for RF
        X, H = self._compute_X_H_generic(mml, z_vals, RF)
        
        # 3. Polytope Gating
        if require_polytope and (not self.is_in_polytope(roots, X, RF)):
            return RF(0)
            
        # 4. Compute Intrinsic t
        diff_X = X - a0
        try:
            t_center = B_mat.solve_right(diff_X)
        except:
            BT = B_mat.transpose()
            t_center = (BT * B_mat).solve_right(BT * diff_X)
            
        s_dict = compute_s_ij(lambdas, tildes, self.n)
        s_vec = vector(RF, [s_dict[edge] for edge in self.basis_edges])
        
        tx = x_spinor[1]
        ty = y_spinor[1]
        y_center = vector(RF, list(s_vec) + [tx, ty])
        
        # 5. Omega_11
        H_int = B_mat.transpose() * H * B_mat
        det_H = H_int.det()
        if abs(det_H) < RF(1e-25): return RF(0)
        Omega_11 = RF(1.0) / det_H
        
        # 6. Perturbations
        dim_t = lattice.dim
        delta = RF(1e-6)
        
        grads_t = []
        grads_y = []
        
        attempts = 0
        while len(grads_t) < dim_t and attempts < 30:
            attempts += 1
            
            d_ts = [rnd.gauss(0, 1) * float(delta) for _ in range(self.n)]
            d_tildes_vectors = {}
            for i in range(self.n-2):
                d_tildes_vectors[i] = vector(RF, [rnd.gauss(0, 1) * float(delta), rnd.gauss(0, 1) * float(delta)])
                
            d_tx = rnd.gauss(0, 1) * float(delta)
            d_ty = rnd.gauss(0, 1) * float(delta)
            
            ts_p = [ts[i] + RF(d_ts[i]) for i in range(self.n)]
            ts_tilde_free_p = {i: ts_tilde_free[i] + d_tildes_vectors[i] for i in range(self.n-2)}
            x_s_p = vector(RF, [1, tx + RF(d_tx)])
            y_s_p = vector(RF, [1, ty + RF(d_ty)])
            
            l_p, t_p, C_p = self._compute_point_data_generic(ts_p, ts_tilde_free_p, x_s_p, y_s_p, RF)
            if l_p is None: continue
            
            z_p = self._compute_z_vals(l_p, t_p, C_p, edge_order, RF)
            X_p, _ = self._compute_X_H_generic(mml, z_p, RF)
            
            diff_p = X_p - a0
            try:
                t_p_intrinsic = B_mat.solve_right(diff_p)
            except:
                BT = B_mat.transpose()
                t_p_intrinsic = (BT * B_mat).solve_right(BT * diff_p)
                
            s_dict_p = compute_s_ij(l_p, t_p, self.n)
            s_vec_p = vector(RF, [s_dict_p[edge] for edge in self.basis_edges])
            y_p = vector(RF, list(s_vec_p) + [x_s_p[1], y_s_p[1]])
            
            grads_t.append((t_p_intrinsic - t_center) / delta)
            grads_y.append((y_p - y_center) / delta)
            
        if len(grads_t) < dim_t: return RF(0)
        
        T = matrix(RF, grads_t[:dim_t]).transpose()
        Y = matrix(RF, grads_y[:dim_t]).transpose()
        
        det_Y = Y.det()
        if abs(det_Y) < RF(1e-25): return RF(0)
        det_T = T.det()
        
        return Omega_11 * det_T / det_Y
        
    def get_sign_pattern(self, ts, ts_tilde_free, x_spinor, y_spinor, RF):
        lambdas, tildes, C = self._compute_point_data_generic(ts, ts_tilde_free, x_spinor, y_spinor, RF)
        if lambdas is None: return None
        
        key = 0
        for i in range(self.n):
            for j in range(i+1, self.n):
                ang = bracket(lambdas[i], lambdas[j])
                sq = bracket(tildes[j], tildes[i])
                if abs(ang) < RF(1e-25): ang = RF(1e-25)
                val = (sq / ang) * C[i] * C[j]
                
                bit = 1 if val > 0 else 0
                key = (key << 1) | bit
        return key

    def _compute_X_H_generic(self, mml, z_values, Field):
        # Re-implementation of MomentMapLaplacian.compute_X_H for generic fields
        dim_M = mml.dim_M
        M = matrix(Field, dim_M, dim_M)
        
        for k, (u, v) in enumerate(mml.edges):
            val = z_values[k]
            if val == 0: continue
            for r, c, sgn in mml.dM_updates[k]:
                M[r, c] += Field(sgn) * val
                
        try:
            Minv = M.inverse()
        except:
            raise ValueError("Singular M")
            
        X = [Field(0)] * mml.num_edges
        num_edges = mml.num_edges
        
        # Precompute Bs
        Bs = []
        for k in range(num_edges):
            updates = mml.dM_updates[k]
            if not updates:
                Bs.append(None)
                continue
            
            # Sparse mult Minv * dMk
            # Bk is dense dim_M x dim_M
            Bk = matrix(Field, dim_M, dim_M)
            for r, c, sgn in updates:
                # Add sgn * Minv column r to column c of Bk
                # Wait. Bk = Minv * dMk.
                # Col c of Bk is sum_l Minv[l, r] * dMk[r, c]
                # dMk[r,c] = sgn.
                # So Col c of Bk += sgn * Col r of Minv?
                # Yes.
                for row_idx in range(dim_M):
                    Bk[row_idx, c] += Field(sgn) * Minv[row_idx, r]
            Bs.append(Bk)
            
            # Trace for X
            # Tr(Minv * dMk) = Tr(Bk)
            tr = Bk.trace()
            X[k] = z_values[k] * tr
            
        H = matrix(Field, num_edges, num_edges)
        for i in range(num_edges):
            if Bs[i] is None: continue
            for j in range(i, num_edges):
                if Bs[j] is None: continue
                
                # Tr(Bi * Bj)
                term = (Bs[i] * Bs[j]).trace()
                val = - z_values[i] * z_values[j] * term
                if i == j:
                    val += X[i]
                    
                H[i, j] = val
                if i != j: H[j, i] = val
                
        return vector(Field, X), H

if __name__ == "__main__":
    # Quick Test
    print("Testing FiberGaugeEvaluator...")
    evaluator = FiberGaugeEvaluator(6)
    
    # Random point
    ts = [rnd.uniform(0, 10) for _ in range(6)]
    ts_tilde_free = {i: vector(RR, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for i in range(4)}
    x_s = vector(RR, [1, -2.0])
    y_s = vector(RR, [1, 12.0])
    
    val = evaluator.evaluate_chart_fiber_gauge([0,1,2], ts, ts_tilde_free, x_s, y_s, prec=53)
    print(f"Value (prec=53): {val}")
    
    val_high = evaluator.evaluate_chart_fiber_gauge([0,1,2], ts, ts_tilde_free, x_s, y_s, prec=200)
    print(f"Value (prec=200): {val_high}")
