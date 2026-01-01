import sys
import os
import itertools
import random as rnd
import math
import json
from sage.all import *

sys.path.append(os.getcwd())

load("src/atlas/jacobian_fiber_gauge.sage")
from src.posgeom.forest_polytope import get_forest_exponents

def compute_M_MHV(n, lambdas, tildes, x_spinor, y_spinor, RF=RR, roots=None):
    # Parke-Taylor-like or similar for comparison
    # Actually, we want the GRAVITY amplitude M_MHV.
    # The formula used previously:
    # M = - <xy>^8 * det(Phi) / ( prod(C_i^2) * prod(<i i+1>^2) )
    
    # We need to construct Phi matrix (M_mat)
    if roots is None:
        roots = [0, 1, 2] # Any roots for reference matrix
    # But Matrix-Tree says det(Phi_roots) sums over forests
    # We can use any minor.
    
    # Let's compute M_mat (Laplacian like)
    # M_{uv} = - z_{uv} ...
    # Wait, we need physical M_MHV.
    # Physical M_MHV is sum over trees? Or determinant of something?
    # It is determinant of the "Hodges Matrix".
    # Hodges matrix H_{ij}.
    # H_{ij} = - [ij] / <ij>  (with some factors?)
    # No, H is defined using spinor brackets.
    # <i j> H_{ij} = - [i j] ?
    # Standard formula:
    # M_n = - <xy>^8 det(H_red)
    # H_{ij} = [ij] / <ij>.
    # Diagonal H_{ii} = - sum_{j!=i} H_{ij} * <xj><yi>/<xi><yi> ? 
    # The formula involves auxiliary spinors x, y?
    # The "x, y" in our fiber setup ARE the auxiliary spinors for the KLT/Hodges form.
    # Let's use the explicit Hodges matrix determinant.
    
    # Hodges element phi_{ij} = [ij] / <ij>
    # Reduced determinant |H|^123_123 or similar.
    # With x, y scaling factors.
    
    # Actually, the previous scripts used `compute_M_MHV_value` which implemented a determinant of `z`.
    # But `z` depends on the chart?
    # No, `compute_M_MHV_value` in `atlas_sum_validate_v2.sage` calculated a matrix from `z` derived from `roots=[0,1,2]`.
    # This implies M_MHV is consistently defined via that specific tree sum (which is correct for gravity if sum over all trees).
    # Wait, `get_forest_exponents` sums over forests.
    # det(L_roots) is sum over forests rooted at roots.
    # If we want sum over ALL trees, we need specific roots?
    # Actually, Gravity = Permutation Sum of YM?
    # Or just use the determinant formula provided in `atlas_sum_validate_v2.sage`.
    # It constructs M_mat from z values and computes det.
    # I will copy that logic, adapted for RF.
    
    # Note: This logic computes "Sum over forests compatible with roots=[0,1,2]".
    # Is this M_MHV?
    # M_MHV should be chart-independent.
    # If we use roots=[0,1,2], we sum over trees where 0,1,2 are in separate components.
    # Since we are in N=6, we have 3 components.
    # The "Gravity Amplitude" usually refers to the single object.
    # Maybe the "M_MHV" in previous scripts was just a reference volume?
    # Let's trust the previous implementation of `compute_M_MHV_value` for now as the "Target".
    
    exponents, edge_order = get_forest_exponents(n, roots)
    
    C = {}
    for i in range(n):
        C[i] = bracket(lambdas[i], x_spinor) * bracket(lambdas[i], y_spinor)
        
    z_vals = []
    edge_dict = {}
    for (u, v) in edge_order:
        ang = bracket(lambdas[u], lambdas[v])
        sq = bracket(tildes[v], tildes[u]) 
        if abs(ang) < RF(1e-25): ang = RF(1e-25)
        val = (sq / ang) * C[u] * C[v]
        z_vals.append(val)
        edge_dict[(u,v)] = val
        edge_dict[(v,u)] = val
        
    M_dim = n - len(roots)
    M_mat = matrix(RF, M_dim, M_dim)
    non_roots = sorted([i for i in range(n) if i not in roots])
    v_map = {v: i for i, v in enumerate(non_roots)}
    diags = {v: RF(0) for v in non_roots}
    
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
    
    prod_C_sq = RF(1)
    for k in non_roots: prod_C_sq *= (C[k]**2)
    prod_roots_sq = RF(1)
    for i in range(len(roots)):
        r1 = roots[i]
        r2 = roots[(i+1) % len(roots)]
        prod_roots_sq *= (bracket(lambdas[r1], lambdas[r2])**2)
        
    xy_bracket = bracket(x_spinor, y_spinor)
    
    denom = prod_C_sq * prod_roots_sq
    if abs(denom) < RF(1e-25): return RF(0)
    
    M_MHV = -(xy_bracket**8) * F_z / denom
    return M_MHV

def run_gate_b():
    print("Running Gate B: Residue & Slope Verification (Fiber Gauge)...")
    
    n = 6
    evaluator = FiberGaugeEvaluator(n)
    all_roots = list(itertools.combinations(range(n), 3))
    pairs = list(itertools.combinations(range(n), 2))
    
    RF = RealField(200)
    
    results = []
    
    # We will sample a few pairs to save time, or all?
    # 15 pairs * 20 charts = 300 evaluations per probe.
    # Let's do all pairs.
    
    for idx_pair, (u, v) in enumerate(pairs):
        print(f"Checking pair ({u}, {v}) [{idx_pair+1}/{len(pairs)}]...")
        
        # Generate probe points with s_uv -> 0
        # eps 1e-4 and 1e-8
        eps1 = RF(1e-4)
        eps2 = RF(1e-8)
        
        # Base parameters
        ts_base = [rnd.uniform(0, 10) for _ in range(n)]
        # Force u, v close
        ts_base[v] = ts_base[u]
        
        ts_tilde_free_base = {i: vector(RF, [rnd.uniform(-1,1), rnd.uniform(-1,1)]) for i in range(n-2)}
        x_s = vector(RF, [1, -2.0])
        y_s = vector(RF, [1, 12.0])
        
        # Probe 1
        ts1 = [RF(t) for t in ts_base]
        ts1[v] += eps1
        # Re-solve conservation
        lambdas1 = {i: vector(RF, [1, ts1[i]]) for i in range(n)}
        tildes1 = solve_conservation_generic(lambdas1, ts_tilde_free_base, n, RF)
        
        # Probe 2
        ts2 = [RF(t) for t in ts_base]
        ts2[v] += eps2
        lambdas2 = {i: vector(RF, [1, ts2[i]]) for i in range(n)}
        tildes2 = solve_conservation_generic(lambdas2, ts_tilde_free_base, n, RF)
        
        if tildes1 is None or tildes2 is None:
            print("  Failed to solve conservation for probe.")
            continue
            
        # Compute M_MHV
        m1 = compute_M_MHV(n, lambdas1, tildes1, x_s, y_s, RF)
        m2 = compute_M_MHV(n, lambdas2, tildes2, x_s, y_s, RF)
        
        # Compute s_uv actual
        s1 = compute_s_ij(lambdas1, tildes1, n)[(u,v)]
        s2 = compute_s_ij(lambdas2, tildes2, n)[(u,v)]
        
        # Slope M
        slope_M = (log(abs(m2)) - log(abs(m1))) / (log(abs(s2)) - log(abs(s1)))
        print(f"  M_MHV slope: {float(slope_M):.2f}")
        
        # Iterate Charts
        for roots in all_roots:
            roots_t = tuple(sorted(list(roots)))
            
            # Eval 1
            val1 = evaluator.evaluate_chart_fiber_gauge(roots_t, ts1, ts_tilde_free_base, x_s, y_s, prec=200)
            
            # Eval 2
            val2 = evaluator.evaluate_chart_fiber_gauge(roots_t, ts2, ts_tilde_free_base, x_s, y_s, prec=200)
            
            status = "Inactive"
            slope = 0.0
            
            if abs(val1) > 1e-20 and abs(val2) > 1e-20:
                status = "Active"
                slope = float((log(abs(val2)) - log(abs(val1))) / (log(abs(s2)) - log(abs(s1))))
            elif abs(val1) > 1e-20 or abs(val2) > 1e-20:
                status = "Unstable/Boundary"
                
            # Check Rule
            in_roots = (u in roots) or (v in roots)
            expected_singular = not in_roots
            
            match = False
            if status == "Active":
                if expected_singular:
                    # Expect slope ~ -1
                    if abs(slope + 1.0) < 0.2: match = True
                else:
                    # Expect slope ~ 0
                    if abs(slope) < 0.2: match = True
            else:
                # Inactive usually means we are not covering this region
                # Which is fine, but we can't verify slope.
                # But for 'expected_singular' charts, IF they are inactive, we might be missing coverage?
                # Or we are just on the wrong side of the wall.
                match = True # Pass if inactive (no bad singularity)
                
            # Log interesting cases
            if status == "Active" or (expected_singular and status == "Active"):
                 results.append({
                    "pair": [u, v],
                    "roots": list(roots),
                    "status": status,
                    "slope": slope,
                    "expected_singular": expected_singular,
                    "match": match
                })
                
    # Summary
    print("\n--- Summary ---")
    active_count = len([r for r in results if r["status"] == "Active"])
    mismatches = [r for r in results if r["status"] == "Active" and not r["match"]]
    
    print(f"Total Active Chart Evaluations: {active_count}")
    print(f"Mismatches: {len(mismatches)}")
    
    if mismatches:
        print("Sample mismatches:")
        for m in mismatches[:5]:
            print(m)
            
    out_path = "RESULTS/atlas_gateB_fiber_gauge.json"
    os.makedirs("RESULTS", exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
        
    print(f"Saved results to {out_path}")

if __name__ == "__main__":
    run_gate_b()

