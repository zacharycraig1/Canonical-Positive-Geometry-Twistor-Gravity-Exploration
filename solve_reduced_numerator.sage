
import sys
from sage.all import *
load("correct_klt_proof.sage")
from itertools import combinations

def solve_reduced_numerator():
    print("Solving Reduced Numerator P(2,3,4,5)...")
    
    # Target Degrees on nodes 2,3,4,5
    # 2: 2
    # 3: 1
    # 4: 1
    # 5: 2
    
    # Possible edges between {2,3,4,5}
    nodes = [2,3,4,5]
    edges = list(combinations(nodes, 2))
    # Edges: (2,3), (2,4), (2,5), (3,4), (3,5), (4,5)
    
    # Total edges needed: sum(degrees)/2 = 6/2 = 3.
    # We need to choose 3 edges (with replacement) to match degrees.
    
    basis_graphs = []
    
    import itertools
    # Iterate all combinations of 3 edges
    for p in itertools.combinations_with_replacement(edges, 3):
        # Check degrees
        degs = {n: 0 for n in nodes}
        for (u, v) in p:
            degs[u] += 1
            degs[v] += 1
            
        if degs[2]==2 and degs[3]==1 and degs[4]==1 and degs[5]==2:
            basis_graphs.append(p)
            
    print(f"Found {len(basis_graphs)} basis graphs.")
    for bg in basis_graphs:
        print(f"  {bg}")
        
    # Fit coefficients
    print("Fitting coefficients...")
    
    def get_denom(twistor):
        val = QQ(1)
        for i in range(6):
            val *= twistor.get_angle(i, (i+1)%6)
        return val**2
        
    X = []
    y = []
    
    # Need > len(basis) points
    for i in range(20):
        Z = sample_positive_Z_moment_curve(n=6, seed=3000+i)
        tw = MomentumTwistor(n=6, Z=Z, check_domain=True)
        if not tw.domain_ok: continue
            
        H = hodges_6pt_mhv(tw)[0]
        D = get_denom(tw)
        N = H * D
        
        # Reduced Numerator P = N / <0 1>^8
        ang01 = tw.get_angle(0, 1)
        P_val = N / (ang01**8)
        
        row = []
        for bg in basis_graphs:
            val = QQ(1)
            for (u, v) in bg:
                val *= tw.get_angle(u, v)
            row.append(val)
            
        X.append(row)
        y.append(P_val)
        
    import numpy as np
    X_mat = np.array(X, dtype=float)
    y_vec = np.array(y, dtype=float)
    
    coeffs, resid, rank, s = np.linalg.lstsq(X_mat, y_vec, rcond=None)
    
    print(f"Residual: {resid}")
    print("\nResult:")
    
    for i, c in enumerate(coeffs):
        if abs(c) > 1e-4:
            term = "".join([f"<{u}{v}>" for (u,v) in basis_graphs[i]])
            print(f"  {c:+.4f} * {term}")
            
    # Verify integer exactness using Sage
    print("\nVerifying Exact Rational Coefficients...")
    # Solve exactly: X * c = y
    # Just take first few rows
    M = matrix(QQ, X[:len(basis_graphs)], sparse=False)
    Y = vector(QQ, y[:len(basis_graphs)])
    try:
        sol = M.solve_right(Y)
        print("Exact Solution found:")
        for i, c in enumerate(sol):
            if c != 0:
                term = "".join([f"<{u}{v}>" for (u,v) in basis_graphs[i]])
                print(f"  {c} * {term}")
    except Exception as e:
        print(f"Exact solve failed: {e}")

if __name__ == "__main__":
    solve_reduced_numerator()








