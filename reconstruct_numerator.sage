
import sys
from sage.all import *
load("correct_klt_proof.sage")
from itertools import combinations_with_replacement

def solve_diophantine(target_degrees):
    # Find all sets of edges (i,j) that satisfy the degrees
    # edges can be repeated
    
    n = 6
    num_edges = sum(target_degrees) // 2
    print(f"Target Degrees: {target_degrees}")
    print(f"Num Edges: {num_edges}")
    
    # Possible edges
    possible_edges = list(combinations(range(n), 2))
    
    solutions = []
    
    # We need to choose 11 edges.
    # The search space is roughly (15 choose 11) with replacement.
    # 15+11-1 choose 11 = 25 choose 11 = 4,457,400. Too big?
    # No, many won't satisfy degrees.
    # Better: Recursive search by node.
    
    # Or simplified constraint solver.
    # Use IntegerVectors? No.
    
    # Backtracking search
    # State: current_degrees (initially 0), current_edges
    
    # Optimization: Prioritize edges connected to high degree nodes (0, 1).
    # Node 0 needs 8 edges.
    # Node 1 needs 8 edges.
    # Node 2 needs 2.
    # Node 3 needs 1.
    # Node 4 needs 1.
    # Node 5 needs 2.
    
    # Deterministic construction
    # We can iterate over possible edges connected to Node 3 (only 1 edge).
    # It must connect to something.
    # Say (3, k).
    
    valid_monomials = []
    
    # Helper to format solution
    def to_monomial(edge_counts):
        # edge_counts: dict (i,j) -> count
        return edge_counts

    # Recursive solver
    def solve(node_idx, current_degrees):
        if node_idx == n:
            if all(current_degrees[i] == target_degrees[i] for i in range(n)):
                return [[]]
            return []
        
        # We need to add edges (node_idx, k) for k > node_idx
        # such that degree of node_idx becomes target[node_idx]
        needed = target_degrees[node_idx] - current_degrees[node_idx]
        if needed < 0: return []
        
        # Possible neighbors k > node_idx
        neighbors = range(node_idx + 1, n)
        
        # Distribute 'needed' edges among neighbors
        # x_k = count of edge (node_idx, k)
        # sum x_k = needed
        # And we must not exceed target[k] for the neighbor
        
        # Generate partitions of 'needed' into len(neighbors) bins
        import itertools
        
        res = []
        
        # simple recursion for distribution
        def distribute(k_idx, remain):
            if k_idx == len(neighbors):
                if remain == 0:
                    yield []
                return
            
            neigh = neighbors[k_idx]
            # Max edges we can add to this neighbor
            max_add = target_degrees[neigh] - current_degrees[neigh]
            
            # Try adding c edges
            for c in range(min(remain, max_add) + 1):
                for rest in distribute(k_idx + 1, remain - c):
                    yield [(neigh, c)] + rest

        for dist in distribute(0, needed):
            # Apply distribution
            next_degrees = list(current_degrees)
            edges_added = []
            
            for neigh, count in dist:
                next_degrees[neigh] += count
                if count > 0:
                    edges_added.append(((node_idx, neigh), count))
            
            # Recurse
            sub_sols = solve(node_idx + 1, next_degrees)
            for sol in sub_sols:
                res.append(edges_added + sol)
                
        return res

    raw_sols = solve(0, [0]*n)
    print(f"Found {len(raw_sols)} basis monomials.")
    return raw_sols

def reconstruct_numerator():
    print("Reconstructing N(Z) from Angle Bracket Basis...")
    
    target_weights = [8, 8, 2, 1, 1, 2]
    basis_structs = solve_diophantine(target_weights)
    
    # Evaluate basis on sample points
    # Need len(basis) points
    num_basis = len(basis_structs)
    
    # Gather data
    points = [] # list of (twistor, N_val)
    
    def get_denom(twistor):
        val = QQ(1)
        for i in range(6):
            val *= twistor.get_angle(i, (i+1)%6)
        return val**2
        
    print(f"Collecting {num_basis + 10} points...")
    
    X = []
    y = []
    
    for i in range(num_basis + 10):
        Z = sample_positive_Z_moment_curve(n=6, seed=1000+i)
        tw = MomentumTwistor(n=6, Z=Z, check_domain=True)
        if not tw.domain_ok: continue
            
        H = hodges_6pt_mhv(tw)[0]
        D = get_denom(tw)
        N_val = H * D
        
        # Eval basis
        row = []
        for edges in basis_structs:
            # edges is list of ((u,v), count)
            val = QQ(1)
            for (u,v), count in edges:
                ang = tw.get_angle(u, v)
                val *= (ang ** count)
            row.append(val)
            
        X.append(row)
        y.append(N_val)
        
    # Solve linear system
    import numpy as np
    X_mat = np.array(X, dtype=float) # Use float for solve? Or Sage?
    # Use Sage matrix for exact arithmetic if possible, but might be slow/huge coefficients.
    # Let's use float first to check residuals.
    y_vec = np.array(y, dtype=float)
    
    try:
        # Least squares
        coeffs, resid, rank, s = np.linalg.lstsq(X_mat, y_vec, rcond=None)
        
        print(f"Residual: {resid}")
        print(f"Rank: {rank} / {num_basis}")
        
        # Identify non-zero coefficients
        significant = []
        for i, c in enumerate(coeffs):
            if abs(c) > 1e-4:
                monomial_str = ""
                for (u,v), count in basis_structs[i]:
                    monomial_str += f"<{u}{v}>^{count} "
                significant.append((c, monomial_str))
                
        print("\nReconstructed Formula:")
        significant.sort(key=lambda x: abs(x[0]), reverse=True)
        for c, m in significant:
            print(f"  {c:+.4f} * {m}")
            
        # Check integer coefficients
        print("\nChecking for integer coefficients...")
        is_integer = True
        for c, m in significant:
            if abs(c - round(c)) > 1e-3:
                is_integer = False
                break
        
        if is_integer:
            print("[SUCCESS] Coefficients appear to be integers!")
        else:
            print("[WARNING] Coefficients are not integers.")
            
    except Exception as e:
        print(f"Solver failed: {e}")

if __name__ == "__main__":
    reconstruct_numerator()







