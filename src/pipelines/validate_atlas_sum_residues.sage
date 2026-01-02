import sys
import os
import json
import itertools
import random as rnd
import math
from sage.all import *

DEBUG_PRINT = True

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

def solve_conservation(lambdas, tildes_free, n, fixed_indices=None):
    if fixed_indices is None: fixed_indices = [n-2, n-1]
    idx1, idx2 = fixed_indices
    free_indices = [k for k in range(n) if k not in fixed_indices]
    
    rhs_0 = 0
    rhs_1 = 0
    for i in free_indices:
        if i in tildes_free:
             t_vec = tildes_free[i]
             rhs_0 -= lambdas[i] * t_vec[0]
             rhs_1 -= lambdas[i] * t_vec[1]
             
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

def compute_M_MHV_value(n, lambdas, tildes, x_spinor, y_spinor):
    roots = [0, 1, 2]
    exponents, edge_order = get_forest_exponents(n, roots)
    C = {}
    for i in range(n): C[i] = bracket(lambdas[i], x_spinor) * bracket(lambdas[i], y_spinor)
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
    for v in non_roots: M_mat[v_map[v], v_map[v]] = diags[v]
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
    return -(xy_bracket**8) * F_z / (prod_C_sq * prod_roots_sq)

def evaluate_chart_with_jacobian(n, roots, ts, ts_tilde_free, x_spinor, y_spinor, lattice, mml, edge_order, basis_edges, fixed_indices):
    lambdas = {i: vector(RR, [1, ts[i]]) for i in range(n)}
    tildes = solve_conservation(lambdas, ts_tilde_free, n, fixed_indices)
    if tildes is None: 
        if DEBUG_PRINT: print("Tildes None")
        return 0
    
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
    if abs(det_H) < 1e-20: 
        if DEBUG_PRINT: print(f"det_H=0 roots={roots}")
        return 0
    Omega = 1.0 / det_H
    
    dim_t = lattice.dim
    dim_s = len(basis_edges)
    
    delta = 1e-5
    grads_t = []
    grads_s = []
    diff_X = X - lattice.a0
    try:
        # Try solve_right (B t = diff)
        t_center = B_mat.solve_right(diff_X)
    except: 
        if DEBUG_PRINT: print("t_center solve fail (trying LS)")
        # Try least squares
        try:
            # (B^T B) t = B^T diff
            # t = (B^T B)^-1 B^T diff
            # Use solve_right on B^T B
            BT = B_mat.transpose()
            t_center = (BT * B_mat).solve_right(BT * diff_X)
        except:
             if DEBUG_PRINT: print("t_center LS fail")
             return 0
    s_center_dict = compute_s_ij(lambdas, tildes, n)
    s_center = vector(RR, [s_center_dict[edge] for edge in basis_edges])
    
    for _ in range(dim_s + 10):
        d_ts = [rnd.gauss(0, 1) * delta for _ in range(n)]
        d_tildes = [vector(RR, [rnd.gauss(0, 1) * delta, rnd.gauss(0, 1) * delta]) for _ in range(n-2)]
        ts_p = [ts[i] + d_ts[i] for i in range(n)]
        
        # tildes_free needs to handle indices
        free_indices = [k for k in range(n) if k not in fixed_indices]
        # map back
        ts_tilde_free_p = {}
        for idx, k in enumerate(free_indices):
             # d_tildes is list size n-2.
             ts_tilde_free_p[k] = ts_tilde_free[k] + d_tildes[idx]
             
        lambdas_p = {i: vector(RR, [1, ts_p[i]]) for i in range(n)}
        tildes_p = solve_conservation(lambdas_p, ts_tilde_free_p, n, fixed_indices)
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
            t_p = B_mat.solve_right(X_p - lattice.a0)
        except: 
            # Least squares fallback
            try:
                diff_p = X_p - lattice.a0
                BT = B_mat.transpose()
                t_p = (BT * B_mat).solve_right(BT * diff_p)
            except: continue
        s_p_dict = compute_s_ij(lambdas_p, tildes_p, n)
        s_p = vector(RR, [s_p_dict[edge] for edge in basis_edges])
        grads_t.append((t_p - t_center)/delta)
        grads_s.append((s_p - s_center)/delta)
        if len(grads_t) >= dim_s: break
        
    if len(grads_t) < dim_s: 
        # print("Not enough gradients")
        return 0
    T_mat = matrix(RR, grads_t[:dim_s]).transpose()
    S_mat = matrix(RR, grads_s[:dim_s]).transpose()
    det_S = S_mat.det()
    
    # T_mat is (dim_t x dim_s). If dim_t > dim_s, we need to pick a minor.
    # We pick the top square block (first dim_s intrinsic coords).
    if T_mat.nrows() > dim_s:
        det_T = T_mat[:dim_s, :].det()
    else:
        det_T = T_mat.det()
        
    if abs(det_S) < 1e-20:
        if DEBUG_PRINT: print("det_S=0")
        return 0
    
    J = det_T / det_S
    if DEBUG_PRINT: 
         # Only print if non-zero or interesting
         if abs(Omega*J) > 1e-10 or roots == (0,1,2):
             print(f"Roots={roots} Om={Omega:.2e} det_T={det_T:.2e} det_S={det_S:.2e} J={J:.2e} Val={Omega*J:.2e}")
             
    return Omega * J

def validate_atlas_sum():
    n = 6
    cover_path = "RESULTS/atlas_cover_n6.json"
    if not os.path.exists(cover_path): return
    with open(cover_path, 'r') as f: cover_data = json.load(f)
    charts = [tuple(item['roots']) for item in cover_data]
    print(f"Validating sum over {len(charts)} cover charts...")
    
    chart_geoms = []
    for roots in charts:
        exponents, edge_order = get_forest_exponents(n, roots)
        lattice = IntrinsicLattice(exponents)
        mml = MomentMapLaplacian(n, roots, edge_order)
        chart_geoms.append({
            "roots": roots, "lattice": lattice, "mml": mml, "edge_order": edge_order
        })
        
    basis_edges = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (1,3), (1,4), (1,5)]
    pairs = list(itertools.combinations(range(n), 2))
    
    print(f"{'Channel':<10} | {'M Slope':<10} | {'Sum Slope':<10} | {'Match?'}")
    print("-" * 50)
    
    for pair in pairs:
        i_idx, j_idx = pair
        candidates_fixed = [k for k in range(n) if k != i_idx and k != j_idx]
        fixed_indices = candidates_fixed[-2:]
        free_indices = [k for k in range(n) if k not in fixed_indices]
        
        ts_base = [rnd.uniform(0, 10) for _ in range(n)]
        ts_tilde_base_list = [vector(RR, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for _ in range(n-2)]
        ts_tilde_free = {k: v for k, v in zip(free_indices, ts_tilde_base_list)}
        tx = rnd.uniform(-5, -1)
        ty = rnd.uniform(11, 15)
        x_s = vector(RR, [1, tx])
        y_s = vector(RR, [1, ty])
        
        epsilons = [1e-4, 1e-5] # Closer
        results = []
        
        for eps in epsilons:
            ts = list(ts_base)
            ts[j_idx] = ts[i_idx] + eps
            
            lambdas = {k: vector(RR, [1, ts[k]]) for k in range(n)}
            tildes = solve_conservation(lambdas, ts_tilde_free, n, fixed_indices)
            if tildes is None: 
                results.append((0, 0))
                continue
            
            M = compute_M_MHV_value(n, lambdas, tildes, x_s, y_s)
            
            sum_Om = 0
            # Debug: check contribution of each chart
            contributions = []
            for geom in chart_geoms:
                val = evaluate_chart_with_jacobian(n, geom['roots'], ts, ts_tilde_free, x_s, y_s, geom['lattice'], geom['mml'], geom['edge_order'], basis_edges, fixed_indices)
                sum_Om += val
                contributions.append(val)
                
            # print(f"Eps {eps}: M={M:.2e} Sum={sum_Om:.2e} Components={['%.2e'%x for x in contributions]}")
            if i_idx == 0 and j_idx == 1:
                print(f"Eps {eps}: M={M:.2e} Sum={sum_Om:.2e} Components={['%.2e'%x for x in contributions]}")
            results.append((M, sum_Om))
            
        m1, s1 = results[0]
        m2, s2 = results[1]
        e1, e2 = epsilons[0], epsilons[1]
        
        if abs(m1) < 1e-20 or abs(m2) < 1e-20: sl_M = 0
        else: sl_M = (math.log(abs(m2)) - math.log(abs(m1))) / (math.log(e2) - math.log(e1))
        
        if abs(s1) < 1e-20 or abs(s2) < 1e-20: sl_S = 0
        else: sl_S = (math.log(abs(s2)) - math.log(abs(s1))) / (math.log(e2) - math.log(e1))
        
        match = "YES" if abs(sl_M - sl_S) < 0.2 and abs(sl_M + 1.0) < 0.5 else "NO"
        print(f"s_{i_idx}{j_idx:<5} | {sl_M:<10.2f} | {sl_S:<10.2f} | {match}")

if __name__ == "__main__":
    validate_atlas_sum()

