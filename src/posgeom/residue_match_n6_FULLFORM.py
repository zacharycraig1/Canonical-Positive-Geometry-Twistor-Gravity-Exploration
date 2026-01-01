import sys
import os
import random
import math
import numpy as np
from sage.all import RR, vector, matrix

sys.path.append(os.getcwd())

from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.intrinsic_lattice import IntrinsicLattice
from src.posgeom.moment_map_laplacian import MomentMapLaplacian

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
            
            # <i j> [j i]
            # [j i] = det(tj, ti) = -(ti0*tj1 - ti1*tj0) = -bracket(ti, tj)
            # Standard conventions vary. Let's use <ij>[ji].
            ang = bracket(li, lj)
            sq = bracket(tj, ti) # Note order for square bracket
            
            val = ang * sq
            s[(i,j)] = val
            s[(j,i)] = val
    return s

def solve_conservation(lambdas, tildes_free, n):
    """
    Given all lambdas and n-2 tildes, solve for last 2 tildes
    to satisfy sum |i> [i| = 0.
    """
    # L_mat * [t_{n-2}, t_{n-1}] = - sum_{0..n-3} |i> [i|
    # Each term |i>[i| is a 2x2 matrix.
    # We want sum_{i=0}^{n-1} lambda_i * tilde_i^T = 0?
    # Usually sum |i> [i| = 0 means sum_i lambda_i^a tilde_i^{\dot{a}} = 0.
    
    # RHS = - sum_{i=0}^{n-3} lambda_i * tilde_i
    rhs_0 = 0
    rhs_1 = 0
    for i in range(n-2):
        rhs_0 -= lambdas[i] * tildes_free[i][0] # x component
        rhs_1 -= lambdas[i] * tildes_free[i][1] # y component
        
    # We need tilde_{n-2}, tilde_{n-1} such that:
    # lambda_{n-2} * tilde_{n-2}^x + lambda_{n-1} * tilde_{n-1}^x = rhs_0
    # lambda_{n-2} * tilde_{n-2}^y + lambda_{n-1} * tilde_{n-1}^y = rhs_1
    
    # This is a linear system for scalar components?
    # No, vectors.
    # lambda_{n-2} (2-vec) * scalar + ... = vector
    # Matrix M = [lambda_{n-2}, lambda_{n-1}] (2x2)
    # M * [tilde_{n-2}^x, tilde_{n-1}^x]^T = rhs_0
    
    M = matrix(RR, [[lambdas[n-2][0], lambdas[n-1][0]], 
                    [lambdas[n-2][1], lambdas[n-1][1]]])
    
    try:
        sol_x = M.solve_right(rhs_0)
        sol_y = M.solve_right(rhs_1)
        
        tildes = {}
        for i in range(n-2):
            tildes[i] = tildes_free[i]
        tildes[n-2] = vector(RR, [sol_x[0], sol_y[0]])
        tildes[n-1] = vector(RR, [sol_x[1], sol_y[1]])
        return tildes
    except:
        return None

def compute_physics_and_geometry(n, roots, lambdas, tildes, x_spinor, y_spinor, lattice, mml, edge_order):
    # 1. Physics: M_MHV (Hodges)
    # Compute z_ij for map
    C = {}
    for i in range(n):
        C[i] = bracket(lambdas[i], x_spinor) * bracket(lambdas[i], y_spinor)
        
    z_vals = []
    for (u, v) in edge_order:
        ang = bracket(lambdas[u], lambdas[v])
        sq = bracket(tildes[u], tildes[v]) # check sign? [uv]
        if abs(ang) < 1e-15: ang = 1e-15
        # Phase V Map: z_ij = [ij] / <ij> * C_i * C_j
        # Note: sq/ang is [uv]/<uv>.
        val = (sq / ang) * C[u] * C[v]
        z_vals.append(val)
        
    # M_MHV using Matrix Tree
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
    for k in non_roots:
        prod_C_sq *= (C[k]**2)
    prod_roots_sq = 1
    for i in range(len(roots)):
        r1 = roots[i]
        r2 = roots[(i+1) % len(roots)]
        prod_roots_sq *= (bracket(lambdas[r1], lambdas[r2])**2)
        
    xy_bracket = bracket(x_spinor, y_spinor)
    if abs(prod_C_sq) < 1e-20 or abs(prod_roots_sq) < 1e-20:
        M_MHV = 0
    else:
        M_MHV = -(xy_bracket**8) * F_z / (prod_C_sq * prod_roots_sq)
        
    # 2. Geometry: Canonical Function Omega
    X, H = mml.compute_X_H(z_vals)
    B_mat = lattice.B
    H_int = B_mat.transpose() * H * B_mat
    det_H = H_int.det()
    
    if abs(det_H) < 1e-20:
        Omega = float('inf')
    else:
        Omega = 1.0 / det_H
        
    # 3. Intrinsic coords t
    # X = a0 + B*t => t = (B^T B)^-1 B^T (X - a0)
    # Use least squares or solve
    diff = X - lattice.a0
    try:
        t_vec = B_mat.solve_left(diff) # B * t = diff
    except:
        t_vec = vector(RR, [0]*lattice.dim) # Fallback
        
    return M_MHV, Omega, t_vec

def residue_match_fullform():
    n = 6
    roots = [0, 1, 2]
    
    print("Initializing geometry...")
    exponents, edge_order = get_forest_exponents(n, roots)
    lattice = IntrinsicLattice(exponents)
    mml = MomentMapLaplacian(n, roots, edge_order)
    
    # Setup Probe: s_01 -> 0
    # Parametrize by eps
    
    ts_base = [random.uniform(0, 10) for _ in range(n)]
    ts_tilde_base = [vector(RR, [random.uniform(-1,1), random.uniform(-1,1)]) for _ in range(n-2)]
    
    tx = random.uniform(-5, -1)
    ty = random.uniform(11, 15)
    x_spinor = vector(RR, [1, tx])
    y_spinor = vector(RR, [1, ty])
    
    # Basis of s_ij to compute Jacobian against
    # We select 9 edges to form a basis for K_6
    # E.g. (0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5)
    basis_edges = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5)]
    
    print("\n--- Probe: Collinear 0,1 (s_01 -> 0) ---")
    
    epsilons = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    
    results = []
    
    for eps in epsilons:
        # Base Point
        ts = list(ts_base)
        ts[1] = ts[0] + eps # Make 0,1 collinear
        
        lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
        tildes_free = {i: ts_tilde_base[i] for i in range(n-2)}
        tildes = solve_conservation(lambdas, tildes_free, n)
        
        if tildes is None:
            print(f"Skipping eps={eps}: Solver failed")
            continue
            
        # Central Value
        M, Om, t_center = compute_physics_and_geometry(n, roots, lambdas, tildes, x_spinor, y_spinor, lattice, mml, edge_order)
        
        # Jacobian Computation
        # Perturb input params (ts and tildes_free) to get 9 gradients
        # We need 9 linearly independent directions in s-space.
        # Inputs: 6 ts + 8 tildes coords = 14 vars.
        # Random perturbations
        
        delta = eps * 1e-4
        num_vars = 14
        
        grads_t = []
        grads_s = []
        
        for _ in range(12): # Generate more than 9 to be safe
            # Random perturbation vector
            d_ts = [random.gauss(0, 1) * delta for _ in range(n)]
            d_tildes = [vector(RR, [random.gauss(0, 1) * delta, random.gauss(0, 1) * delta]) for _ in range(n-2)]
            
            # Perturbed Point
            ts_p = [ts[i] + d_ts[i] for i in range(n)]
            lambdas_p = {i: vector(RR, [1, ts_p[i]]) for i in range(n)}
            tildes_free_p = {i: tildes_free[i] + d_tildes[i] for i in range(n-2)}
            tildes_p = solve_conservation(lambdas_p, tildes_free_p, n)
            
            if tildes_p is None: continue
            
            # Compute new t and s
            _, _, t_p = compute_physics_and_geometry(n, roots, lambdas_p, tildes_p, x_spinor, y_spinor, lattice, mml, edge_order)
            s_p = compute_s_ij(lambdas_p, tildes_p, n)
            
            # Compute original s
            s_center = compute_s_ij(lambdas, tildes, n)
            
            diff_t = (t_p - t_center) / delta # Scaling doesn't strictly matter for ratio
            diff_s = vector(RR, [s_p[edge] - s_center[edge] for edge in basis_edges]) / delta
            
            grads_t.append(diff_t)
            grads_s.append(diff_s)
            
            if len(grads_t) >= 9: break
            
        # Form Matrices
        if len(grads_t) < 9:
            print(f"eps={eps}: Failed to generate gradients")
            continue
            
        T_mat = matrix(RR, grads_t[:9]).transpose() # 11 x 9
        S_mat = matrix(RR, grads_s[:9]).transpose() # 9 x 9
        
        det_S = S_mat.det()
        if abs(det_S) < 1e-20:
            J_minor = 0
        else:
            # Choose 9 rows of T (first 9 for now - intrinsic coords are arbitrary)
            # Actually, try to find max minor
            # For consistency, let's just pick first 9 rows (indices 0..8)
            # Intrinsic coords usually well-conditioned
            T_sub = T_mat[:9, :] 
            det_T = T_sub.det()
            J_minor = det_T / det_S
            
        results.append({
            "eps": eps,
            "M": M,
            "Om": Om,
            "J": J_minor
        })
        
        print(f"eps={eps:.1e} | M={M:.2e} | Om={Om:.2e} | J={J_minor:.2e} | Om*J={Om*J_minor:.2e}")
        
    # Analysis
    if len(results) >= 2:
        r1, r2 = results[-2], results[-1]
        e1, e2 = r1["eps"], r2["eps"]
        
        def slope(val1, val2):
            if abs(val1) < 1e-20 or abs(val2) < 1e-20: return 0
            return (math.log(abs(val2)) - math.log(abs(val1))) / (math.log(e2) - math.log(e1))
            
        s_M = slope(r1["M"], r2["M"])
        s_Om = slope(r1["Om"], r2["Om"])
        s_Full = slope(r1["Om"] * r1["J"], r2["Om"] * r2["J"])
        
        print("\n--- Scaling Analysis (Power Law) ---")
        print(f"M (Gravity):    ~ eps^{s_M:.2f}")
        print(f"Omega (Scalar): ~ eps^{s_Om:.2f}")
        print(f"Omega * J:      ~ eps^{s_Full:.2f}")
        
        print("\nConclusion:")
        if abs(s_M - s_Om) < 0.2:
            print("[D1] SCALAR MATCH: The canonical function matches gravity scaling!")
        elif abs(s_M - s_Full) < 0.2:
            print("[D3] FULL FORM MATCH: The form (with Jacobian) matches gravity scaling!")
        else:
            print("NO MATCH: Neither scalar nor form matches gravity.")

if __name__ == "__main__":
    residue_match_fullform()


