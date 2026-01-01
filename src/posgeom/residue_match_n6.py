import sys
import os
import argparse
import random
import statistics
from sage.all import RR, vector, matrix

sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.intrinsic_lattice import IntrinsicLattice
from src.posgeom.moment_map_laplacian import MomentMapLaplacian

def bracket(l1, l2):
    return l1[0]*l2[1] - l1[1]*l2[0]

def compute_amplitudes_at_kinematics(n, roots, lambdas, tildes, x, y, lattice, mml, edge_order):
    # Compute z_ij
    C = {}
    for i in range(n):
        C[i] = bracket(lambdas[i], x) * bracket(lambdas[i], y)
        
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
            else: continue
            if u in diags: diags[u] += val
            if v in diags: diags[v] += val
            if u in v_map and v in v_map:
                idx_u = v_map[u]
                idx_v = v_map[v]
                M_mat[idx_u, idx_v] = -val
                M_mat[idx_v, idx_u] = -val
    for v in non_roots:
        M_mat[v_map[v], v_map[v]] = diags[v]
        
    F_z = M_mat.det()
    
    prod_C_sq = 1
    for k in non_roots:
        prod_C_sq *= (C[k]**2)
    prod_roots_sq = 1
    num_roots = len(roots)
    for i in range(num_roots):
        r1 = roots[i]
        r2 = roots[(i+1) % num_roots]
        prod_roots_sq *= (bracket(lambdas[r1], lambdas[r2])**2)
        
    xy_bracket = bracket(x, y)
    if prod_C_sq == 0 or prod_roots_sq == 0:
        M_MHV = 0 # Singular
    else:
        M_MHV = -(xy_bracket**8) * F_z / (prod_C_sq * prod_roots_sq)
        
    # Omega
    X, H = mml.compute_X_H(z_vals)
    B_mat = lattice.B
    H_int = B_mat.transpose() * H * B_mat
    det_H = H_int.det()
    
    if abs(det_H) < 1e-20:
        Omega = float('inf')
    else:
        Omega = 1.0 / det_H
        
    return M_MHV, Omega

def residue_match_n6():
    n = 6
    roots = [0, 1, 2]
    
    print("Initializing geometry...")
    exponents, edge_order = get_forest_exponents(n, roots)
    lattice = IntrinsicLattice(exponents)
    mml = MomentMapLaplacian(n, roots, edge_order)
    
    # Define limit: s_01 -> 0
    # Parametrize by epsilon
    # t_1 = t_0 + epsilon
    
    print("\nChecking limit s_01 -> 0 (Collinear 0,1)...")
    
    epsilons = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    
    # Fix base kinematics
    ts_base = sorted([random.uniform(0, 10) for _ in range(n)])
    ts_tilde = sorted([random.uniform(0, 10) for _ in range(n)])
    tx = random.uniform(-5, -1)
    ty = random.uniform(11, 15)
    x = vector(RR, [1, tx])
    y = vector(RR, [1, ty])
    
    log_eps = []
    log_M = []
    log_Om = []
    
    for eps in epsilons:
        # Modify t_1
        ts = list(ts_base)
        ts[1] = ts[0] + eps
        
        lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
        tildes = {i: vector(RR, [1, ts_tilde[i]]) for i in range(n)}
        
        M, Om = compute_amplitudes_at_kinematics(n, roots, lambdas, tildes, x, y, lattice, mml, edge_order)
        
        from math import log
        if M != 0 and Om != float('inf') and Om != 0:
            log_eps.append(log(eps))
            log_M.append(log(abs(M)))
            log_Om.append(log(abs(Om)))
            
            print(f"eps={eps:.1e}: M={M:.2e}, Om={Om:.2e}, Ratio={M/Om:.2e}")
        else:
            print(f"eps={eps:.1e}: Singular values")
            
    # Compute slopes
    if len(log_eps) > 1:
        slope_M = (log_M[-1] - log_M[0]) / (log_eps[-1] - log_eps[0])
        slope_Om = (log_Om[-1] - log_Om[0]) / (log_eps[-1] - log_eps[0])
        
        print("\nScaling Analysis (Power Law):")
        print(f"Slope M (Gravity): {slope_M:.4f} (Expected ~ -2 for 1/s^2? or -1?)")
        print(f"Slope Omega (Geom): {slope_Om:.4f}")
        
        if abs(slope_M - slope_Om) < 0.1:
            print("SUCCESS: Scaling dimensions match!")
        else:
            print("FAILURE: Scaling dimensions mismatch.")
            
    # Check s_012 -> 0
    print("\nChecking limit s_012 -> 0 (Collinear 0,1,2)...")
    log_eps = []
    log_M = []
    log_Om = []
    
    for eps in epsilons:
        ts = list(ts_base)
        ts[1] = ts[0] + eps
        ts[2] = ts[0] + 2*eps
        
        lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
        tildes = {i: vector(RR, [1, ts_tilde[i]]) for i in range(n)}
        
        M, Om = compute_amplitudes_at_kinematics(n, roots, lambdas, tildes, x, y, lattice, mml, edge_order)
        
        from math import log
        if M != 0 and Om != float('inf') and Om != 0:
            log_eps.append(log(eps))
            log_M.append(log(abs(M)))
            log_Om.append(log(abs(Om)))
            print(f"eps={eps:.1e}: M={M:.2e}, Om={Om:.2e}")

    if len(log_eps) > 1:
        slope_M = (log_M[-1] - log_M[0]) / (log_eps[-1] - log_eps[0])
        slope_Om = (log_Om[-1] - log_Om[0]) / (log_eps[-1] - log_eps[0])
        
        print("\nScaling Analysis:")
        print(f"Slope M: {slope_M:.4f}")
        print(f"Slope Omega: {slope_Om:.4f}")
        
        if abs(slope_M - slope_Om) < 0.1:
            print("SUCCESS: Scaling dimensions match!")
        else:
            print("FAILURE: Scaling dimensions mismatch.")

if __name__ == "__main__":
    residue_match_n6()



