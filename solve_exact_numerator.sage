
import sys
from sage.all import *
load("correct_klt_proof.sage")
import itertools

# =============================================================================
# 1. BASIS GENERATION
# =============================================================================

def generate_basis(target_degrees):
    """
    Generate all multigraphs (monomials) with the given degree sequence.
    Returns list of edge lists: [((u,v), count), ...]
    """
    n = len(target_degrees)
    
    def solve(node_idx, current_degrees):
        if node_idx == n - 1:
            # Check last node
            needed = target_degrees[node_idx] - current_degrees[node_idx]
            if needed == 0:
                return [[]]
            return []
            
        needed = target_degrees[node_idx] - current_degrees[node_idx]
        if needed < 0: return []
        
        neighbors = range(node_idx + 1, n)
        
        # Partition 'needed' edges among neighbors
        # Limit for each neighbor is target[neigh] - current[neigh]
        limits = [target_degrees[j] - current_degrees[j] for j in neighbors]
        
        res = []
        
        def partitions(total, k, limits):
            if k == 0:
                if total == 0: yield []
                return
            
            # Max we can put in first bin
            limit = limits[0]
            # Also constrained by total
            max_val = min(total, limit)
            
            for i in range(max_val + 1):
                # Recurse
                for rest in partitions(total - i, k - 1, limits[1:]):
                    yield [i] + rest

        for p in partitions(needed, len(neighbors), limits):
            # p is counts for neighbors
            added_edges = []
            next_degrees = list(current_degrees)
            
            for i, count in enumerate(p):
                if count > 0:
                    neigh = neighbors[i]
                    next_degrees[neigh] += count
                    added_edges.append(((node_idx, neigh), count))
            
            # Recurse
            sub_sols = solve(node_idx + 1, next_degrees)
            for s in sub_sols:
                res.append(added_edges + s)
                
        return res

    return solve(0, [0]*n)

# =============================================================================
# 2. DATA COLLECTION
# =============================================================================

def collect_exact_data(basis_structs, num_samples):
    """
    Collect samples of (Monomial_Values, N_Exact) over QQ.
    """
    X_rows = []
    y_vals = []
    
    count = 0
    seed_offset = 5000
    
    print(f"Collecting {num_samples} samples...")
    
    while count < num_samples:
        # Generate point
        Z = sample_positive_Z_moment_curve(n=6, seed=seed_offset + count)
        tw = MomentumTwistor(n=6, Z=Z, check_domain=True)
        
        if not tw.domain_ok:
            seed_offset += 1
            continue
            
        # Compute Exact Hodges
        H_res = hodges_6pt_mhv(tw) # Returns (val, reason)
        if isinstance(H_res, tuple):
            H = H_res[0]
        else:
            H = H_res
            
        if H is None:
            seed_offset += 1
            continue
            
        # Compute Cyclic Denom
        D = QQ(1)
        for i in range(6):
            D *= tw.get_angle(i, (i+1)%6)
        D = D**2
        
        # Numerator Value
        N_val = H * D
        
        # Evaluate Basis
        row = []
        for edges in basis_structs:
            # edges: [((u,v), count), ...]
            val = QQ(1)
            for (u,v), count in edges:
                ang = tw.get_angle(u, v)
                val *= (ang ** count)
            row.append(val)
            
        X_rows.append(row)
        y_vals.append(N_val)
        count += 1
        
        if count % 10 == 0:
            print(f"  Collected {count}/{num_samples}")
            
    return X_rows, y_vals

# =============================================================================
# 3. EXACT SOLVER
# =============================================================================

def solve_exact():
    print("Generating Basis for Weights [8, 8, 2, 1, 1, 2]...")
    target_weights = [8, 8, 2, 1, 1, 2]
    basis = generate_basis(target_weights)
    print(f"Basis size: {len(basis)}")
    
    # We need at least len(basis) samples
    num_samples = len(basis) + 5
    X, y = collect_exact_data(basis, num_samples)
    
    print("Building Matrix...")
    M = matrix(QQ, X)
    Y = vector(QQ, y)
    
    print("Solving linear system exactly...")
    try:
        # Solve M * c = Y
        # Use simple backslash or solve_right
        # Check rank first
        rk = M.rank()
        print(f"Matrix Rank: {rk} / {len(basis)}")
        
        if rk < len(basis):
            print("Warning: System is underdetermined. Solution may not be unique.")
            # We can still find a particular solution
            
        coeffs = M.solve_right(Y)
        
        print("\n=== EXACT SOLUTION FOUND ===")
        terms = []
        for i, c in enumerate(coeffs):
            if c != 0:
                # Format monomial
                monomial_parts = []
                for (u,v), count in basis[i]:
                    if count == 1:
                        monomial_parts.append(f"<{u}{v}>")
                    else:
                        monomial_parts.append(f"<{u}{v}>^{count}")
                monomial_str = "".join(monomial_parts)
                terms.append((c, monomial_str))
                
        # Sort by magnitude of coefficient
        terms.sort(key=lambda x: abs(x[0]), reverse=True)
        
        for c, m in terms:
            print(f"{c} * {m}")
            
        # Verify solution on one extra point
        print("\nVerifying on test point...")
        Z_test = sample_positive_Z_moment_curve(n=6, seed=9999)
        tw_test = MomentumTwistor(n=6, Z=Z_test, check_domain=True)
        H_test = hodges_6pt_mhv(tw_test)[0]
        
        D_test = QQ(1)
        for i in range(6): D_test *= tw_test.get_angle(i, (i+1)%6)
        D_test = D_test**2
        
        N_actual = H_test * D_test
        N_pred = QQ(0)
        
        # Re-evaluate polynomial
        # Need to map 'm' string back to calculation or just reuse basis[i] with coeffs[i]
        
        for i, c in enumerate(coeffs):
            if c == 0: continue
            val = QQ(1)
            for (u,v), count in basis[i]:
                val *= (tw_test.get_angle(u, v) ** count)
            N_pred += c * val
            
        print(f"Actual: {N_actual}")
        print(f"Pred:   {N_pred}")
        print(f"Diff:   {N_actual - N_pred}")
        
        if N_actual == N_pred:
            print("[SUCCESS] Exact Match Verified!")
        else:
            print("[FAILURE] Verification Failed.")
            
    except Exception as e:
        print(f"Solve failed: {e}")

if __name__ == "__main__":
    solve_exact()








