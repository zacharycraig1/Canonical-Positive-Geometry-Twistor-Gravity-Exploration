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

# --- Utilities for Gate B ---
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

def solve_conservation(lambdas, tildes_free, n, fixed_indices=None):
    if fixed_indices is None:
        fixed_indices = [n-2, n-1]
        
    idx1, idx2 = fixed_indices
    free_indices = [k for k in range(n) if k not in fixed_indices]
    
    rhs_0 = 0
    rhs_1 = 0
    for i in free_indices:
        # tildes_free should be keyed by index
        if i in tildes_free:
             t_vec = tildes_free[i]
             rhs_0 -= lambdas[i] * t_vec[0]
             rhs_1 -= lambdas[i] * t_vec[1]
        else:
             # Should not happen if tildes_free covers all free indices
             return None
        
    M = matrix(RR, [[lambdas[idx1][0], lambdas[idx2][0]], 
                    [lambdas[idx1][1], lambdas[idx2][1]]])
    try:
        sol_x = M.solve_right(rhs_0)
        sol_y = M.solve_right(rhs_1)
        tildes = {}
        for i in free_indices: tildes[i] = tildes_free[i]
        tildes[idx1] = vector(RR, [sol_x[0], sol_y[0]])
        tildes[idx2] = vector(RR, [sol_x[1], sol_y[1]])
        return tildes
    except:
        return None

def compute_physics_and_geometry(n, roots, lambdas, tildes, x_spinor, y_spinor, lattice, mml, edge_order):
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
        
    # M_MHV
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
    if abs(prod_C_sq) < 1e-20 or abs(prod_roots_sq) < 1e-20: M_MHV = 0
    else: M_MHV = -(xy_bracket**8) * F_z / (prod_C_sq * prod_roots_sq)
        
    # Omega
    X, H = mml.compute_X_H(z_vals)
    B_mat = lattice.B
    H_int = B_mat.transpose() * H * B_mat
    det_H = H_int.det()
    if abs(det_H) < 1e-20: Omega = float('inf')
    else: Omega = 1.0 / det_H
        
    return M_MHV, Omega

def run_gate_b_pair(n, roots, pair):
    # Setup Geometry
    try:
        exponents, edge_order = get_forest_exponents(n, roots)
        lattice = IntrinsicLattice(exponents)
        mml = MomentMapLaplacian(n, roots, edge_order)
    except Exception as e:
        # Some root sets might be invalid (disconnected etc), though for n=6 k=3 usually fine
        # Assuming provided valid roots
        return "SetupError"
    
    # Probe s_ij -> 0
    i_idx, j_idx = pair
    
    # Choose fixed indices for conservation solver that are NOT i or j
    # This avoids singular matrix when i,j become collinear
    candidates_fixed = [k for k in range(n) if k != i_idx and k != j_idx]
    fixed_indices = candidates_fixed[-2:] # Pick last two available
    free_indices = [k for k in range(n) if k not in fixed_indices]
    
    # Random base config
    ts_base = [rnd.uniform(0, 10) for _ in range(n)]
    # tildes_free now needs to map free_indices to vectors
    ts_tilde_base_list = [vector(RR, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for _ in range(n-2)]
    ts_tilde_free = {k: v for k, v in zip(free_indices, ts_tilde_base_list)}
    
    tx = rnd.uniform(-5, -1)
    ty = rnd.uniform(11, 15)
    x_s = vector(RR, [1, tx])
    y_s = vector(RR, [1, ty])
    
    epsilons = [1e-4, 1e-6]
    results = []
    
    for eps in epsilons:
        ts = list(ts_base)
        # Make j collinear with i
        ts[j_idx] = ts[i_idx] + eps
        
        lambdas = {k: vector(RR, [1, ts[k]]) for k in range(n)}
        tildes = solve_conservation(lambdas, ts_tilde_free, n, fixed_indices=fixed_indices)
        if tildes is None: return "SolverFail"
        
        M, Om = compute_physics_and_geometry(n, roots, lambdas, tildes, x_s, y_s, lattice, mml, edge_order)
        results.append((M, Om))
        
    # Check Slope
    m1, o1 = results[0]
    m2, o2 = results[1]
    e1, e2 = epsilons[0], epsilons[1]
    
    if abs(m1) < 1e-20 or abs(m2) < 1e-20: slope_M = 0
    else: slope_M = (math.log(abs(m2)) - math.log(abs(m1))) / (math.log(e2) - math.log(e1))
    
    if abs(o1) < 1e-20 or abs(o2) < 1e-20: slope_Om = 0
    else: slope_Om = (math.log(abs(o2)) - math.log(abs(o1))) / (math.log(e2) - math.log(e1))
    
    return slope_M, slope_Om

# --- Main Scan ---
def scan_all_pairs():
    n = 6
    all_roots = list(itertools.combinations(range(n), 3))
    pairs = list(itertools.combinations(range(n), 2))
    
    print(f"Scanning {len(all_roots)} root sets across {len(pairs)} pairs for Gate B...")
    
    final_data = []
    
    # We want to know for each pair, which roots work.
    # Structure: List of objects
    # { "pair": [i,j], "roots": [r1,r2,r3], "slopes": [m, o], "match": bool }
    
    # To save time, maybe we group by pair first?
    # Or loop pairs then roots.
    
    total_checks = len(pairs) * len(all_roots)
    count = 0
    
    for pair in pairs:
        print(f"\n--- Channel s_{pair[0]}{pair[1]} ---")
        
        for roots in all_roots:
            roots = list(roots)
            
            try:
                s_m, s_o = run_gate_b_pair(n, roots, pair)
                if isinstance(s_m, str): # Error string
                    print(f"Roots {roots}: {s_m}")
                    continue
                    
                match = (abs(s_m - s_o) < 0.2)
                
                # Only print matches or significant failures to keep log clean?
                # Actually, print summary line
                # print(f"Roots {roots}: M~{s_m:.1f} Om~{s_o:.1f} Match={match}")
                
                final_data.append({
                    "channel": f"s_{int(pair[0])}{int(pair[1])}",
                    "pair": [int(p) for p in pair],
                    "roots": [int(r) for r in roots],
                    "gravity_slope": float(s_m),
                    "omega_slope": float(s_o),
                    "match": bool(match)
                })
                
            except Exception as e:
                print(f"Roots {roots} Pair {pair} Error: {e}")
                
            count += 1
            if count % 20 == 0:
                print(f"Progress: {count}/{total_checks}")

    # Save
    out_path = "RESULTS/atlas_sweep_all_pairs_n6.json"
    with open(out_path, "w") as f:
        json.dump(final_data, f, indent=2)
    print(f"Saved complete scan to {out_path}")

if __name__ == "__main__":
    scan_all_pairs()

