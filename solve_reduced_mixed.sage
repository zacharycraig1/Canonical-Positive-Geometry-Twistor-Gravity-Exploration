
import sys
from sage.all import *
load("correct_klt_proof.sage")
from itertools import combinations, combinations_with_replacement

def solve_reduced_numerator_mixed():
    print("Solving Reduced Numerator P(2,3,4,5) with 4-brackets...")
    
    # Target Degrees on nodes 2,3,4,5: [2:2, 3:1, 4:1, 5:2]
    
    # 1. 2-bracket basis (same as before)
    nodes = [2,3,4,5]
    edges = list(combinations(nodes, 2))
    basis_2 = []
    for p in combinations_with_replacement(edges, 3):
        degs = {n: 0 for n in nodes}
        for (u, v) in p:
            degs[u] += 1
            degs[v] += 1
        if degs[2]==2 and degs[3]==1 and degs[4]==1 and degs[5]==2:
            basis_2.append(('2b', p))
            
    # 2. 4-bracket basis
    # We need a 4-bracket involving 2,3,4,5?
    # <2 3 4 5> uses 1 degree from each.
    # Remaining needed: 2:1, 3:0, 4:0, 5:1.
    # Need 2-bracket <2 5>.
    # So Term: <2 3 4 5> <2 5>.
    # Are there other 4-brackets?
    # No, only 4 nodes available.
    
    term_4b = ('4b', [(2,3,4,5), (2,5)]) # Format: 4-bracket tuple, then 2-bracket tuple
    
    basis = basis_2 + [term_4b]
    
    print(f"Basis size: {len(basis)}")
    for b in basis:
        print(f"  {b}")
        
    # Fit
    print("Fitting...")
    
    X = []
    y = []
    
    # Use higher precision or exact QQ to check
    # But for finding the relation, float is enough initially
    
    for i in range(20):
        Z = sample_positive_Z_moment_curve(n=6, seed=6000+i)
        tw = MomentumTwistor(n=6, Z=Z, check_domain=True)
        if not tw.domain_ok: continue
            
        H = hodges_6pt_mhv(tw)[0]
        # Denom
        D = QQ(1)
        for k in range(6): D *= tw.get_angle(k, (k+1)%6)
        D = D**2
        
        N = H * D
        ang01 = tw.get_angle(0, 1)
        P_val = N / (ang01**8)
        
        row = []
        for type, parts in basis:
            val = QQ(1)
            if type == '2b':
                for (u,v) in parts:
                    val *= tw.get_angle(u, v)
            elif type == '4b':
                # parts[0] is 4-bracket indices, parts[1] is 2-bracket
                # <i j k l>
                idx = parts[0]
                b4 = tw.get_four_bracket(*idx)
                b2 = tw.get_angle(parts[1][0], parts[1][1])
                val = b4 * b2
            row.append(val)
            
        X.append(row)
        y.append(P_val)
        
    import numpy as np
    X_mat = np.array(X, dtype=float)
    y_vec = np.array(y, dtype=float)
    
    coeffs, resid, rank, s = np.linalg.lstsq(X_mat, y_vec, rcond=None)
    
    print(f"Residual: {resid}")
    
    print("\nApprox Solution:")
    for i, c in enumerate(coeffs):
        if abs(c) > 1e-4:
            if basis[i][0] == '2b':
                term = "".join([f"<{u}{v}>" for u,v in basis[i][1]])
            else:
                term = f"<{basis[i][1][0]}> <{basis[i][1][1]}>"
            print(f"  {c:+.4f} * {term}")
            
    # Try EXACT solve
    print("\nAttempting Exact Rational Solve...")
    try:
        M = matrix(QQ, X[:len(basis)])
        Y = vector(QQ, y[:len(basis)])
        sol = M.solve_right(Y)
        print("EXACT SOLUTION FOUND!")
        for i, c in enumerate(sol):
            if c != 0:
                if basis[i][0] == '2b':
                    term = "".join([f"<{u}{v}>" for u,v in basis[i][1]])
                else:
                    idx4 = basis[i][1][0]
                    idx2 = basis[i][1][1]
                    term = f"<{idx4[0]}{idx4[1]}{idx4[2]}{idx4[3]}><{idx2[0]}{idx2[1]}>"
                print(f"  {c} * {term}")
    except Exception as e:
        print(f"Exact solve failed: {e}")

if __name__ == "__main__":
    solve_reduced_numerator_mixed()




