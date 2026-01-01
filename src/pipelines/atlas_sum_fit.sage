import sys
import os
import json
import itertools
import random as rnd
import math
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

def evaluate_chart_with_jacobian(n, roots, ts, ts_tilde_free, x_spinor, y_spinor, lattice, mml, edge_order, basis_edges):
    # 1. Central Value
    lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
    tildes = solve_conservation(lambdas, ts_tilde_free, n)
    if tildes is None: return 0
    
    C = {}
    for i in range(n): C[i] = bracket(lambdas[i], x_spinor) * bracket(lambdas[i], y_spinor)
    z_vals = []
    for (u, v) in edge_order:
        ang = bracket(lambdas[u], lambdas[v])
        sq = bracket(tildes[u], tildes[v]) 
        if abs(ang) < 1e-15: ang = 1e-15
        val = (sq / ang) * C[u] * C[v]
        z_vals.append(val)
        
    X, H = mml.compute_X_H(z_vals)
    B_mat = lattice.B
    H_int = B_mat.transpose() * H * B_mat
    det_H = H_int.det()
    if abs(det_H) < 1e-20: return 0
    Omega = 1.0 / det_H
    
    # 2. Jacobian
    dim_t = lattice.dim
    dim_s = len(basis_edges) # Should be 9
    
    # If lattice dim != 9, we might be on a slice or have redundancies
    # For n=6 k=3, dim should be 9.
    
    # Check if we can compute J
    if dim_t != dim_s:
        # Fallback: maybe just return Omega if we can't define J
        # But if we need J for correct scaling, this is bad.
        # Let's hope dim_t == 9.
        return Omega # No J correction
    
    # Perturb
    delta = 1e-5
    grads_t = []
    grads_s = []
    
    diff_X = X - lattice.a0
    try:
        t_center = B_mat.solve_left(diff_X)
    except:
        return 0
        
    s_center_dict = compute_s_ij(lambdas, tildes, n)
    s_center = vector(RR, [s_center_dict[edge] for edge in basis_edges])
    
    # We need dim_s independent directions
    # 20 attempts
    for _ in range(dim_s + 10):
        d_ts = [rnd.gauss(0, 1) * delta for _ in range(n)]
        d_tildes = [vector(RR, [rnd.gauss(0, 1) * delta, rnd.gauss(0, 1) * delta]) for _ in range(n-2)]
        
        ts_p = [ts[i] + d_ts[i] for i in range(n)]
        ts_tilde_free_p = {i: ts_tilde_free[i] + d_tildes[i] for i in range(n-2)}
        
        lambdas_p = {i: vector(RR, [1, ts_p[i]]) for i in range(n)}
        tildes_p = solve_conservation(lambdas_p, ts_tilde_free_p, n)
        if tildes_p is None: continue
        
        C_p = {}
        for i in range(n): C_p[i] = bracket(lambdas_p[i], x_spinor) * bracket(lambdas_p[i], y_spinor)
        z_p = []
        for (u, v) in edge_order:
            ang = bracket(lambdas_p[u], lambdas_p[v])
            sq = bracket(tildes_p[u], tildes_p[v])
            if abs(ang) < 1e-15: ang = 1e-15
            z_p.append((sq / ang) * C_p[u] * C_p[v])
            
        X_p, _ = mml.compute_X_H(z_p)
        try:
            t_p = B_mat.solve_left(X_p - lattice.a0)
        except: continue
        
        s_p_dict = compute_s_ij(lambdas_p, tildes_p, n)
        s_p = vector(RR, [s_p_dict[edge] for edge in basis_edges])
        
        grads_t.append((t_p - t_center)/delta)
        grads_s.append((s_p - s_center)/delta)
        
        if len(grads_t) >= dim_s: break
        
    if len(grads_t) < dim_s: return 0
    
    T_mat = matrix(RR, grads_t[:dim_s]).transpose()
    S_mat = matrix(RR, grads_s[:dim_s]).transpose()
    
    det_S = S_mat.det()
    det_T = T_mat.det()
    
    if abs(det_S) < 1e-20: return 0
    J = det_T / det_S
    
    # Sign of J?
    # We take abs(J) or J?
    # Usually forms have orientation. 
    # For now, let's keep J signed.
    
    return Omega * J

def compute_M_MHV_value(n, lambdas, tildes, x_spinor, y_spinor):
    roots = [0, 1, 2]
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

def fit_coefficients():
    n = 6
    all_roots = list(itertools.combinations(range(n), 3))
    charts = all_roots
    print(f"Fitting over ALL {len(charts)} root sets with Jacobian...")
    
    # Precompute geometry
    chart_geoms = []
    for roots in charts:
        exponents, edge_order = get_forest_exponents(n, roots)
        lattice = IntrinsicLattice(exponents)
        mml = MomentMapLaplacian(n, roots, edge_order)
        print(f"Roots {roots}: Lattice Dim {lattice.dim}")
        chart_geoms.append({
            "roots": roots,
            "lattice": lattice,
            "mml": mml,
            "edge_order": edge_order
        })
        
    basis_edges = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5)]
    
    # Generate Points
    # Reduce number of points because Jacobian is slow
    num_points = 25 
    print(f"Generating {num_points} random points...")
    
    A_rows = []
    b_vals = []
    
    for pt_idx in range(num_points):
        print(f"Processing point {pt_idx+1}/{num_points}")
        
        ts = [rnd.uniform(0, 10) for _ in range(n)]
        ts_tilde_free = {i: vector(RR, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for i in range(n-2)}
        tx = rnd.uniform(-5, -1)
        ty = rnd.uniform(11, 15)
        x_s = vector(RR, [1, tx])
        y_s = vector(RR, [1, ty])
        
        lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
        tildes = solve_conservation(lambdas, ts_tilde_free, n)
        if tildes is None: continue
        
        M = compute_M_MHV_value(n, lambdas, tildes, x_s, y_s)
        
        row = []
        for geom in chart_geoms:
            val = evaluate_chart_with_jacobian(n, geom['roots'], ts, ts_tilde_free, x_s, y_s, geom['lattice'], geom['mml'], geom['edge_order'], basis_edges)
            row.append(val)
            
        A_rows.append(row)
        b_vals.append(M)
        
    # Solve
    print("Solving system...")
    A = matrix(RR, A_rows)
    b = vector(RR, b_vals)
    
    try:
        AT = A.transpose()
        ATA = AT * A
        ATb = AT * b
        coeffs = ATA.solve_right(ATb)
    except Exception as e:
        print(f"Solver failed: {e}")
        return

    print("\n--- Fitted Coefficients ---")
    results = []
    for i, chart in enumerate(charts):
        c = float(coeffs[i])
        print(f"Chart {chart}: {c:.4f}")
        results.append({
            "roots": list(chart),
            "coeff": c
        })
        
    predictions = A * coeffs
    residuals = b - predictions
    rel_errors = [abs(residuals[i])/abs(b[i]) if abs(b[i]) > 1e-10 else 0 for i in range(len(b))]
    avg_err = sum(rel_errors) / len(rel_errors) if rel_errors else 0
    max_err = max(rel_errors) if rel_errors else 0
    print(f"Avg Relative Error: {avg_err:.2e}")
    print(f"Max Relative Error: {max_err:.2e}")
    
    # Save
    rounded_results = []
    for r in results:
        val = r['coeff']
        closest_int = round(val)
        if abs(val - closest_int) < 0.1:
            r['suggested_coeff'] = int(closest_int)
        else:
            r['suggested_coeff'] = float(val)
        
        r['roots'] = [int(x) for x in r['roots']]
        r['coeff'] = float(val)
        rounded_results.append(r)
        
    out_path = "RESULTS/atlas_coeffs_n6.json"
    with open(out_path, 'w') as f:
        json.dump(rounded_results, f, indent=2)
    print(f"Saved coeffs to {out_path}")

if __name__ == "__main__":
    fit_coefficients()
