import sys
import os
import json
import itertools
import random as rnd
from sage.all import *

# Path setup
sys.path.append(os.getcwd())

# Load dependency
load("src/atlas/jacobian_kernel_gauge.sage")

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

def bracket(l1, l2):
    return l1[0]*l2[1] - l1[1]*l2[0]

def compute_M_MHV_value(n, lambdas, tildes, x_spinor, y_spinor):
    # Standard MHV formula or numerical evaluation
    # Using the same one as before for consistency
    roots = [0, 1, 2]
    # We can use the forest sum directly or just copy the logic
    # For speed, let's copy the logic from previous files
    
    from src.posgeom.forest_polytope import get_forest_exponents
    exponents, edge_order = get_forest_exponents(n, roots)
    
    C = {}
    for i in range(n):
        C[i] = bracket(lambdas[i], x_spinor) * bracket(lambdas[i], y_spinor)
        
    z_vals = []
    for (u, v) in edge_order:
        ang = bracket(lambdas[u], lambdas[v])
        sq = bracket(tildes[u], tildes[v]) 
        if abs(ang) < 1e-15: ang = 1e-15
        val = (sq / ang) * C[u] * C[v]
        z_vals.append(val)
        
    M_dim = n - len(roots)
    M_mat = matrix(RR, M_dim, M_dim)
    non_roots = sorted([i for i in range(n) if i not in roots])
    v_map = {v: i for i, v in enumerate(non_roots)}
    diags = {v: 0 for v in non_roots}
    edge_dict = {edge: val for edge, val in zip(edge_order, z_vals)}
    
    for u in range(n):
        for v in range(u+1, n):
            if (u, v) in edge_dict:
                val = edge_dict[(u, v)]
                if u in diags: diags[u] += val
                if v in diags: diags[v] += val
                if u in v_map and v in v_map:
                    M_mat[v_map[u], v_map[v]] = -val
                    M_mat[v_map[v], v_map[u]] = -val
                    
    for v in non_roots:
        M_mat[v_map[v], v_map[v]] = diags[v]
        
    F_z = M_mat.det()
    
    prod_C_sq = 1
    for k in non_roots: prod_C_sq *= (C[k]**2)
    prod_roots_sq = 1
    for i in range(len(roots)):
        r1 = roots[i]
        r2 = roots[(i+1) % len(roots)]
        prod_roots_sq *= (bracket(lambdas[r1], lambdas[r2])**2)
        
    xy_bracket = bracket(x_spinor, y_spinor)
    if abs(prod_C_sq) < 1e-20 or abs(prod_roots_sq) < 1e-20: return 0
    
    M_MHV = -(xy_bracket**8) * F_z / (prod_C_sq * prod_roots_sq)
    return M_MHV

def validate_atlas_sum():
    print("Validating Atlas Sum with Canonical Jacobian...")
    n = 6
    
    # Load Signs
    signs_path = "RESULTS/atlas_signs.json"
    atlas_signs = {}
    if os.path.exists(signs_path):
        with open(signs_path, 'r') as f:
            data = json.load(f)
            for item in data:
                roots = tuple(sorted(item['roots']))
                atlas_signs[roots] = item['sign']
    else:
        print("Warning: Signs file not found. Using all +1.")
        
    # Initialize Evaluator
    evaluator = CanonicalJacobianEvaluator(n)
    
    # 1. Pointwise Check
    print("\n--- Pointwise Equality Check ---")
    
    # We will also attempt to FIT the coefficients here to see what they SHOULD be.
    # This acts as a check on our topological signs.
    
    num_points = 15
    rows = []
    b_vals = []
    
    all_roots_list = list(itertools.combinations(range(n), 3))
    
    for k in range(num_points):
        ts = [rnd.uniform(0, 10) for _ in range(n)]
        ts_tilde_free = {i: vector(RR, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for i in range(n-2)}
        x_s = vector(RR, [1, -2.0])
        y_s = vector(RR, [1, 12.0])
        
        lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
        tildes = solve_conservation(lambdas, ts_tilde_free, n)
        if tildes is None: continue
        
        M_mhv = compute_M_MHV_value(n, lambdas, tildes, x_s, y_s)
        
        row = []
        for roots in all_roots_list:
            roots_t = tuple(sorted(list(roots)))
            val = evaluator.evaluate(roots_t, ts, ts_tilde_free, x_s, y_s)
            row.append(val)
            
        # Check current sum
        current_sum = 0
        for i, roots in enumerate(all_roots_list):
            roots_t = tuple(sorted(list(roots)))
            s = atlas_signs.get(roots_t, 0)
            current_sum += s * row[i]
            
        ratio = current_sum / M_mhv if abs(M_mhv) > 1e-10 else 0
        print(f"Pt {k}: M_MHV={M_mhv:.2e} Sum={current_sum:.2e} Ratio={ratio:.2f}")
        
        rows.append(row)
        b_vals.append(M_mhv)
        
    # Least Squares Fit to find TRUE coefficients
    print("\n--- Fitting True Coefficients ---")
    A = matrix(RR, rows)
    b = vector(RR, b_vals)
    
    try:
        # coefficients x: A*x = b
        # x = (A^T A)^-1 A^T b
        AT = A.transpose()
        x = (AT * A).solve_right(AT * b)
        
        print(f"{'Chart':<15} | {'Sign (Topo)':<10} | {'Fitted Coeff':<15}")
        print("-" * 50)
        
        discrepancies = 0
        new_signs = []
        
        for i, roots in enumerate(all_roots_list):
            roots_t = tuple(sorted(list(roots)))
            topo_sign = atlas_signs.get(roots_t, 0)
            fitted = float(x[i])
            
            # Check agreement
            # If fitted is close to integer
            fitted_int = int(round(fitted))
            match = "OK" if abs(fitted - topo_sign) < 0.1 else "MISMATCH"
            if match == "MISMATCH": discrepancies += 1
            
            print(f"{str(roots):<15} | {topo_sign:<10} | {fitted:<15.4f} [{match}]")
            
            new_signs.append({
                "roots": list(roots),
                "sign": fitted_int, # Assuming it should be integer
                "fitted": fitted
            })
            
        print(f"\nTotal Discrepancies: {discrepancies}")
        
        if discrepancies > 0:
            print("Updating signs with fitted values (assuming they are correct)...")
            out_path = "RESULTS/atlas_signs_fitted.json"
            with open(out_path, 'w') as f:
                json.dump(new_signs, f, indent=2)
            print(f"Saved fitted signs to {out_path}")
            
    except Exception as e:
        print(f"Fitting failed: {e}")

if __name__ == "__main__":
    validate_atlas_sum()

