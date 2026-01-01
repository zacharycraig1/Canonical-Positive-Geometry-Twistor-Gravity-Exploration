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

def run_gate_b(n, roots):
    # Setup Geometry
    exponents, edge_order = get_forest_exponents(n, roots)
    lattice = IntrinsicLattice(exponents)
    mml = MomentMapLaplacian(n, roots, edge_order)
    
    # Probe s_01 -> 0
    ts_base = [rnd.uniform(0, 10) for _ in range(n)]
    ts_tilde_base = [vector(RR, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for _ in range(n-2)]
    tx = rnd.uniform(-5, -1)
    ty = rnd.uniform(11, 15)
    x_s = vector(RR, [1, tx])
    y_s = vector(RR, [1, ty])
    
    epsilons = [1e-4, 1e-6]
    results = []
    
    for eps in epsilons:
        ts = list(ts_base)
        ts[1] = ts[0] + eps
        lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
        tildes = solve_conservation(lambdas, {i: ts_tilde_base[i] for i in range(n-2)}, n)
        if tildes is None: return "Fail"
        
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
    
    return f"M~{slope_M:.1f}, Om~{slope_Om:.1f}"

# --- Main Scan ---
def scan_roots():
    n = 6
    all_roots = list(itertools.combinations(range(n), 3))
    
    print(f"Scanning {len(all_roots)} root sets for Gate B (s_01 -> 0 scaling)...")
    print(f"{'Roots':<15} | {'Gate B Slope':<20} | {'Match?'}")
    print("-" * 50)
    
    results = []
    
    for roots in all_roots:
        roots = list(roots)
        
        try:
            res_str = run_gate_b(n, roots)
            # Parse slopes
            parts = res_str.split(',')
            s_m = float(parts[0].split('~')[1])
            s_o = float(parts[1].split('~')[1])
            
            match = "YES" if abs(s_m - s_o) < 0.2 else "NO"
            print(f"{str(roots):<15} | {res_str:<20} | {match}")
            
            results.append({
                "roots": roots,
                "gate_b": res_str,
                "match": match
            })
        except Exception as e:
            print(f"{str(roots):<15} | Error: {e}")
            
    # Save
    with open("RESULTS/atlas_sweep_n6.json", "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    scan_roots()
