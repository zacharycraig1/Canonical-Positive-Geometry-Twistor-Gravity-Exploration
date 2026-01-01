
import sys
from sage.all import *
load("correct_klt_proof.sage")

def test_hypothesis():
    print("Testing Single Monomial Hypothesis...")
    
    # Hypothesis: N ~ <0 1>^8 <2 5>^2 <3 4>
    # Weights: 8, 8, 2, 1, 1, 2
    
    def get_denom(twistor):
        val = QQ(1)
        for i in range(6):
            val *= twistor.get_angle(i, (i+1)%6)
        return val**2
        
    def get_monomial(tw):
        # <0 1>^8 <2 5>^2 <3 4>
        t1 = tw.get_angle(0, 1)**8
        t2 = tw.get_angle(2, 5)**2
        t3 = tw.get_angle(3, 4)
        return t1 * t2 * t3
        
    Z = sample_positive_Z_moment_curve(n=6, seed=777)
    tw = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    H = hodges_6pt_mhv(tw)[0]
    D = get_denom(tw)
    N_true = H * D
    
    N_hyp = get_monomial(tw)
    
    print(f"N_true: {float(N_true):.4e}")
    print(f"N_hyp:  {float(N_hyp):.4e}")
    
    if N_hyp != 0:
        ratio = N_true / N_hyp
        print(f"Ratio:  {float(ratio):.6f}")
        
        # Check if constant
        Z2 = sample_positive_Z_moment_curve(n=6, seed=888)
        tw2 = MomentumTwistor(n=6, Z=Z2, check_domain=True)
        N_true2 = hodges_6pt_mhv(tw2)[0] * get_denom(tw2)
        N_hyp2 = get_monomial(tw2)
        ratio2 = N_true2 / N_hyp2
        print(f"Ratio2: {float(ratio2):.6f}")
        
        if abs(ratio - ratio2) < 1e-5:
            print("[SUCCESS] Numerator IS the monomial (up to constant)!")
            print(f"Constant factor: {ratio}")
        else:
            print("[FAILURE] Not a single monomial. Searching basis...")
            run_solver()
            
def run_solver():
    print("\nRunning Exhaustive Solver for Basis...")
    target_degrees = [8, 8, 2, 1, 1, 2]
    
    # Simple iterator for this specific case
    # High degree at 0,1 means we almost certainly have <0 1>^k
    # Try k=8, k=7, k=6...
    
    import itertools
    
    edges_pool = list(combinations(range(6), 2))
    
    found_basis = []
    
    # We need to select 11 edges (with replacement) to match degrees
    # This is hard to iterate directly.
    # Better to iterate edge counts.
    
    # Randomized search to find ANY valid graph?
    # No, we need ALL valid graphs to form basis.
    
    # Improved recursive solver
    def solve(node_idx, current_degrees):
        if node_idx == 5:
            # Check if last node degree matches
            needed = target_degrees[5] - current_degrees[5]
            if needed == 0:
                return [[]]
            return []
            
        needed = target_degrees[node_idx] - current_degrees[node_idx]
        if needed < 0: return []
        
        neighbors = range(node_idx + 1, 6)
        
        res = []
        
        # Distribute 'needed' edges to neighbors
        # We can put at most 'target[neigh] - current[neigh]' edges to 'neigh'
        
        # Partition 'needed' into len(neighbors) parts
        # This is composition of integer 'needed' into k parts.
        
        def partitions(n, k, limits):
            if k == 0:
                if n == 0: yield []
                return
            
            # limits[0] is max we can give to first bin
            max_val = min(n, limits[0])
            for i in range(max_val + 1):
                for p in partitions(n - i, k - 1, limits[1:]):
                    yield [i] + p
                    
        limits = [target_degrees[j] - current_degrees[j] for j in neighbors]
        
        for p in partitions(needed, len(neighbors), limits):
            # p is list of counts for edges (node_idx, neighbor)
            added_edges = []
            next_degrees = list(current_degrees)
            
            for i, count in enumerate(p):
                neigh = neighbors[i]
                next_degrees[neigh] += count
                if count > 0:
                    added_edges.append(((node_idx, neigh), count))
            
            # Recurse
            sub = solve(node_idx + 1, next_degrees)
            for s in sub:
                res.append(added_edges + s)
                
        return res

    basis_structs = solve(0, [0]*6)
    print(f"Found {len(basis_structs)} basis monomials.")
    
    # If found, Fit
    if basis_structs:
        fit_basis(basis_structs)

def fit_basis(basis_structs):
    # Collect data
    num_basis = len(basis_structs)
    print(f"Fitting {num_basis} monomials...")
    
    def get_denom(twistor):
        val = QQ(1)
        for i in range(6):
            val *= twistor.get_angle(i, (i+1)%6)
        return val**2

    X = []
    y = []
    
    import numpy as np
    
    for i in range(num_basis + 10):
        Z = sample_positive_Z_moment_curve(n=6, seed=2000+i)
        tw = MomentumTwistor(n=6, Z=Z, check_domain=True)
        if not tw.domain_ok: continue
            
        H = hodges_6pt_mhv(tw)[0]
        D = get_denom(tw)
        N_val = H * D
        
        row = []
        for edges in basis_structs:
            val = QQ(1)
            for (u,v), count in edges:
                val *= (tw.get_angle(u, v) ** count)
            row.append(val)
        X.append(row)
        y.append(N_val)
        
    X_mat = np.array(X, dtype=float)
    y_vec = np.array(y, dtype=float)
    
    coeffs, resid, rank, s = np.linalg.lstsq(X_mat, y_vec, rcond=None)
    
    print(f"Residual: {resid}")
    
    # Print formula
    print("\nFormula:")
    for i, c in enumerate(coeffs):
        if abs(c) > 1e-4:
            monomial = ""
            for (u,v), count in basis_structs[i]:
                monomial += f"<{u}{v}>^{count}"
            print(f"  {c:+.4f} * {monomial}")

if __name__ == "__main__":
    test_hypothesis()







