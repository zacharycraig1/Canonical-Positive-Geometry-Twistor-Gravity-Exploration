import sys
import os
import time
import argparse
from math import isclose

# Ensure we can import from src
sys.path.append(os.getcwd())

from sage.all import QQ, RR, matrix, vector, ComplexField, RealField

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.intrinsic_lattice import IntrinsicLattice
from src.posgeom.moment_map_laplacian import MomentMapLaplacian

def bracket(l1, l2):
    return l1[0]*l2[1] - l1[1]*l2[0]

def generate_positive_kinematics_moment_curve(n):
    import random
    ts = sorted([random.uniform(0, 10) for _ in range(n)])
    
    lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
    ts_tilde = sorted([random.uniform(0, 10) for _ in range(n)])
    tildes = {i: vector(RR, [1, ts_tilde[i]]) for i in range(n)}
    
    # x, y
    tx = random.uniform(-5, -1) # Outside range
    ty = random.uniform(11, 15) # Outside range
    x = vector(RR, [1, tx])
    y = vector(RR, [1, ty])
    
    return lambdas, tildes, x, y

def map_sweep_n6(samples=20):
    n = 6
    roots = [0, 1, 2]
    
    print(f"Running map sweep for n={n}, roots={roots}, samples={samples}")
    
    # 1. Setup Geometry (once)
    exponents, edge_order = get_forest_exponents(n, roots)
    lattice = IntrinsicLattice(exponents)
    mml = MomentMapLaplacian(n, roots, edge_order)
    
    print(f"Lattice dim: {lattice.dim}")
    print(f"Covolume: {lattice.covolume}")
    
    norm_ratios = []
    
    for s in range(samples):
        # 2. Kinematics
        lambdas, tildes, x, y = generate_positive_kinematics_moment_curve(n)
        
        # 3. Compute z_ij
        C = {}
        for i in range(n):
            C[i] = bracket(lambdas[i], x) * bracket(lambdas[i], y)
            
        z_vals = []
        for (u, v) in edge_order:
            ang = bracket(lambdas[u], lambdas[v])
            sq = bracket(tildes[u], tildes[v])
            if abs(ang) < 1e-10: ang = 1e-10
            val = (sq / ang) * C[u] * C[v]
            z_vals.append(val)
            
        # 4. Compute Amplitude M_MHV
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
        # Note: M_MHV formula with explicit prefactors
        M_MHV = -(xy_bracket**8) * F_z / (prod_C_sq * prod_roots_sq)
        
        # 5. Compute Omega
        X, H = mml.compute_X_H(z_vals)
        B_mat = lattice.B
        H_int = B_mat.transpose() * H * B_mat
        det_H = H_int.det()
        
        if abs(det_H) < 1e-15: continue
        Omega = 1.0 / det_H
        
        ratio = M_MHV / Omega
        ratio_norm = ratio / (xy_bracket**8)
        
        norm_ratios.append(ratio_norm)
        
        if s % 5 == 0:
            print(f"Sample {s}: Ratio/xy^8 = {ratio_norm:.6e}")

    if not norm_ratios:
        print("No valid samples.")
        return

    n_samples = len(norm_ratios)
    sum_val = sum(norm_ratios)
    mean_val = sum_val / n_samples
    if n_samples > 1:
        variance = sum((x - mean_val)**2 for x in norm_ratios) / (n_samples - 1)
        stdev_val = variance.sqrt()
    else:
        stdev_val = 0
    
    print("\nResults (Normalized by <xy>^8):")
    print(f"Mean:    {float(mean_val):.6e}")
    print(f"Std Dev: {float(stdev_val):.6e}")
    if mean_val != 0:
        rel_std = abs(stdev_val/mean_val)
        print(f"Rel Std Dev: {float(rel_std):.6e}")
        if rel_std < 1e-5:
            print("SUCCESS: Ratio is stable!")
        else:
            print("FAILURE: Ratio is not stable.")
    else:
        print("Mean is zero.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=20)
    args = parser.parse_args()
    map_sweep_n6(args.samples)
