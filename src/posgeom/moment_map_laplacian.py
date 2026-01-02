import sys
import os
import time

# Ensure we can import from src
sys.path.append(os.getcwd())

try:
    from sage.all import Matrix, vector, QQ, ZZ
except ImportError:
    pass # Will fail if not run with sage, handled by checks

class MomentMapLaplacian:
    def __init__(self, n, roots, edge_ordering=None):
        self.n = n
        self.roots = sorted(roots)
        self.non_roots = sorted([i for i in range(n) if i not in roots])
        
        # M is indexed by non_roots
        # Map non_root vertex index to matrix index 0..k-1
        self.v_to_idx = {v: i for i, v in enumerate(self.non_roots)}
        self.dim_M = len(self.non_roots)
        
        # Edge ordering
        if edge_ordering is None:
            # Canonical ordering
            self.edges = []
            for i in range(n):
                for j in range(i + 1, n):
                    self.edges.append((i, j))
        else:
            self.edges = edge_ordering
            
        self.num_edges = len(self.edges)
        
        # Precompute derivatives of M
        # dM[edge_idx] = sparse matrix or list of (row, col, val)
        self.dM_updates = []
        
        for k, (u, v) in enumerate(self.edges):
            # Edge (u, v) with u < v
            updates = []
            
            # Diagonal terms in Laplacian: sum_{k!=i} z_{ik}
            # z_{uv} appears in L_{uu} and L_{vv} with +1
            
            if u in self.v_to_idx:
                idx_u = self.v_to_idx[u]
                updates.append((idx_u, idx_u, 1))
                
            if v in self.v_to_idx:
                idx_v = self.v_to_idx[v]
                updates.append((idx_v, idx_v, 1))
                
            # Off-diagonal terms: L_{uv} = -z_{uv}
            if u in self.v_to_idx and v in self.v_to_idx:
                idx_u = self.v_to_idx[u]
                idx_v = self.v_to_idx[v]
                updates.append((idx_u, idx_v, -1))
                updates.append((idx_v, idx_u, -1))
                
            self.dM_updates.append(updates)
            
    def compute_X_H(self, z_values):
        """
        Computes the moment map X and Hessian H in ambient log-coordinates.
        
        Args:
            z_values: list or vector of edge variables z_{ij} corresponding to self.edges
            
        Returns:
            X: vector of length num_edges
            H: matrix of size num_edges x num_edges
        """
        # Construct M
        M = Matrix(QQ, self.dim_M, self.dim_M)
        
        # We can construct M by summing dM * z
        # But iterating is faster
        # M_{ab} = L_{non_roots[a], non_roots[b]}
        
        # Initialize diagonals to 0
        
        for k, (u, v) in enumerate(self.edges):
            val = z_values[k]
            if val == 0: continue
            
            # Apply updates
            for r, c, sgn in self.dM_updates[k]:
                M[r, c] += sgn * val
                
        # Invert M
        try:
            Minv = M.inverse()
        except ZeroDivisionError:
            raise ValueError("Matrix M is singular")
            
        # Compute Tr(Minv * dM/dz_e) for each e
        # This is equivalent to sum_{i,j} Minv_{ji} * (dM/dz_e)_{ij}
        # Since dM is sparse (at most 4 entries), this is fast.
        
        X = [0] * self.num_edges
        tr_Minv_dM = [0] * self.num_edges
        
        for k in range(self.num_edges):
            updates = self.dM_updates[k]
            if not updates:
                continue
                
            trace_val = 0
            for r, c, sgn in updates:
                trace_val += sgn * Minv[c, r]
            
            tr_Minv_dM[k] = trace_val
            X[k] = z_values[k] * trace_val
            
        # Compute Hessian
        # H_{ef} = z_e z_f [ -Tr(Minv dM_f Minv dM_e) ] + delta_{ef} X_e
        
        H = Matrix(QQ, self.num_edges, self.num_edges)
        
        # We need A_e = Minv * dM_e
        # Computing A_e explicitly is expensive (dense matrix mult).
        # But dM_e is sparse. A_e col j = sum_k Minv_{jk} (dM_e)_{kj}
        # A_e is Minv * (sparse).
        # We need Tr(A_f A_e).
        
        # Optimization:
        # H_{ef} = - z_e z_f * sum_{r,c} (Minv * dM_f)_{rc} * (Minv * dM_e)_{cr}
        # Actually Tr(M^{-1} dM_f M^{-1} dM_e)
        # Let B_f = M^{-1} dM_f. We want Tr(B_f B_e).
        
        # B_e is sparse-ish? No, Minv is dense.
        # But dM_e has few entries.
        # B_e[i, j] = sum_k Minv[i, k] dM_e[k, j]
        # dM_e only has entries at (u,u), (v,v), (u,v), (v,u) (mapped to indices)
        
        # Let's precompute B_e for all e where z_e != 0 (or relevant)
        # Actually we need it for all e?
        # If z_e is 0, then X_e is 0 and H_{ef} = 0 (due to z_e factor), EXCEPT maybe if limit?
        # Assuming z_e != 0 for now (interior of polytope).
        
        Bs = [] # List of matrices or better, just needed components?
        # We need Tr(B_f B_e).
        # Let's store B_e as a dense matrix for now. dimension is small (3x3 for n=6).
        
        for k in range(self.num_edges):
            updates = self.dM_updates[k]
            if not updates:
                Bs.append(None)
                continue
                
            # Construct B_k = Minv * dM_k
            # column j of B_k = Minv * (col j of dM_k)
            # col j of dM_k has non-zeros only at specific rows
            
            # Since dim is small (3), standard matrix mult is fine.
            # dM_k matrix
            dMk = Matrix(QQ, self.dim_M, self.dim_M)
            for r, c, sgn in updates:
                dMk[r, c] = sgn
            
            Bk = Minv * dMk
            Bs.append(Bk)
            
        for i in range(self.num_edges):
            if Bs[i] is None: continue
            
            for j in range(i, self.num_edges):
                if Bs[j] is None: continue
                
                # Tr(B_i * B_j)
                term = (Bs[i] * Bs[j]).trace()
                
                val = - z_values[i] * z_values[j] * term
                
                if i == j:
                    val += X[i]
                
                H[i, j] = val
                if i != j:
                    H[j, i] = val
                    
        return vector(QQ, X), H

def main():
    import argparse
    from src.posgeom.forest_polytope import get_forest_exponents
    from src.posgeom.intrinsic_lattice import IntrinsicLattice
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--samples", type=int, default=1)
    args = parser.parse_args()
    
    n = args.n
    roots = [0, 1, 2]
    
    print(f"Testing MomentMapLaplacian for n={n}, roots={roots}")
    
    # Get edge ordering
    _, edges = get_forest_exponents(n, roots)
    
    mml = MomentMapLaplacian(n, roots, edges)
    
    # Test on random z
    # Make sure to use positive z
    from sage.all import RR
    import random
    
    for s in range(args.samples):
        z_vals = [random.randint(1, 10) for _ in range(mml.num_edges)]
        # Use QQ
        z_vals = [QQ(z) for z in z_vals]
        
        start = time.time()
        X, H = mml.compute_X_H(z_vals)
        dt = time.time() - start
        
        print(f"Sample {s}: Computed X (norm={float(X.norm(2)):.4f}) and H (det={H.det()})")
        print(f"Time: {dt*1000:.2f} ms")
        
        # Check projection
        # We need lattice to project
        exponents, _ = get_forest_exponents(n, roots)
        lattice = IntrinsicLattice(exponents)
        
        # Projected Hessian
        # H_int = B.T * H * B
        H_int = lattice.B.transpose() * H * lattice.B
        print(f"Projected H det: {H_int.det()}")

if __name__ == "__main__":
    main()

