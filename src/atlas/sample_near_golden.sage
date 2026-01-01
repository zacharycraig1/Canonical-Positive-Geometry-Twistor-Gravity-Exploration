# src/atlas/sample_near_golden.sage
#
# Generates a cloud of valid kinematic points near a "golden" seed.
# Used for coefficient fitting (Gate C).
#
# Logic:
#   1. Load RESULTS/golden_kinematics.json
#   2. Perturb continuous variables (u_free, tx, ty) slightly.
#   3. Keep s_free (signs) fixed.
#   4. Reconstruct and verify pos=15 and active_charts > 0.
#   5. Save valid points to RESULTS/valid_points_for_fit.json

import json
import os
import numpy as np
from sage.all import RealField, vector, load

# Load helpers
load("src/atlas/jacobian_fiber_gauge.sage")
load("src/atlas/solve_conservation_pair.sage")

N = 6
TARGET_POINTS = 80
SIGMA_U = 1e-3
SIGMA_T = 1e-3

def main():
    if not os.path.exists("RESULTS/golden_kinematics.json"):
        print("Error: RESULTS/golden_kinematics.json not found.")
        print("Run the search script first to find a golden point.")
        return

    with open("RESULTS/golden_kinematics.json", "r") as f:
        golden = json.load(f)

    print("Loaded golden point.")
    print(f"  ts: {golden['ts']}")
    print(f"  pair: {golden['pair']}")
    
    RF = RealField(200)
    evaluator = FiberGaugeEvaluator(N)
    
    # Base values
    ts = np.array(golden["ts"], dtype=np.float64)
    u_base = np.array(golden["u_free"], dtype=np.float64)
    tx_base = float(golden["tx"])
    ty_base = float(golden["ty"])
    s_free = golden["s_free"] # Fixed discrete signs
    pair = tuple(golden["pair"])
    free_idx = golden["free_idx"]
    
    # We also keep ts fixed for this sampling to keep coefficients constant?
    # Usually we want to sample ts as well for general fitting, but for "nearby" 
    # sampling often ts is fixed or also perturbed. 
    # If we want to fit coefficients *functions of t*, we need to vary t.
    # If we want to fit coefficients at *this* t, we fix t.
    # The instructions say "sample near golden", implying local exploration. 
    # Assuming we want to populate the region for this specific chart/topology.
    # We will perturb ts slightly as well if we want to test robustness, 
    # but for Gate C usually we want *many* points.
    # Let's perturb ts slightly too, to avoid degeneracy.
    
    valid_points = []
    
    # Add the golden point itself first
    valid_points.append(golden)
    
    rng = np.random.default_rng(42)
    attempts = 0
    
    while len(valid_points) < TARGET_POINTS and attempts < TARGET_POINTS * 20:
        attempts += 1
        
        # Perturb
        ts_new = ts + rng.normal(0, SIGMA_T, size=N)
        ts_new.sort() # keep order
        
        u_new = u_base + rng.normal(0, SIGMA_U, size=len(u_base))
        tx_new = tx_base + rng.normal(0, SIGMA_T)
        ty_new = ty_base + rng.normal(0, SIGMA_T)
        
        # Verify order of tx, ty (re-sort if crossed, or reject)
        # Actually tx, ty are just auxiliary vars, their order wrt t defines the interval.
        # We should ensure they stay in the same interval or just check positivity.
        
        # Reconstruct in Sage
        ts_sage = [float(x) for x in ts_new]
        lambdas_sage = {i: vector(RF, [1, ts_sage[i]]) for i in range(N)}
        
        tildes_free_sage = {}
        for k, idx in enumerate(free_idx):
            tildes_free_sage[idx] = vector(RF, [float(s_free[k]), float(u_new[k])])
            
        tildes_recon = solve_conservation_pair(lambdas_sage, tildes_free_sage, pair, N, RF)
        if tildes_recon is None:
            continue
            
        x_sage = vector(RF, [1, float(tx_new)])
        y_sage = vector(RF, [1, float(ty_new)])
        
        # Check Positivity
        def br(u, v): return u[0]*v[1] - u[1]*v[0]
        C_poly = {i: br(lambdas_sage[i], x_sage) * br(lambdas_sage[i], y_sage) for i in range(N)}
        
        pos_count = 0
        z_vals = []
        for i in range(N):
            for j in range(i+1, N):
                ang = br(lambdas_sage[i], lambdas_sage[j])
                sq  = br(tildes_recon[j], tildes_recon[i])
                val = (sq / ang) * C_poly[i] * C_poly[j]
                zf = float(val)
                z_vals.append(zf)
                if zf > 0: pos_count += 1
        
        if pos_count != 15:
            continue
            
        # Check Charts
        tf_final = {int(i): [float(s_free[k]), float(u_new[k])] for k, i in enumerate(free_idx)}
        xs_final = [1.0, float(tx_new)]
        ys_final = [1.0, float(ty_new)]
        
        # We only need >0 charts, not full enumeration for valid points
        active_charts = 0
        hit_roots = []
        
        # Check a few random roots first for speed? Or just check all since N=6 is small (20 roots)
        # 6 choose 3 = 20. Fast enough.
        
        for roots in itertools.combinations(range(N), 3):
            val = evaluator.evaluate_chart_fiber_gauge(roots, ts_sage, tf_final, xs_final, ys_final, require_polytope=True)
            if abs(val) > 1e-20:
                active_charts += 1
                hit_roots.append(list(roots))
                
        if active_charts == 0:
            continue
            
        # Valid!
        pt = {
            "ts": ts_new.tolist(),
            "pair": list(pair),
            "free_idx": list(free_idx),
            "s_free": s_free,
            "u_free": u_new.tolist(),
            "tx": float(tx_new),
            "ty": float(ty_new),
            "z_values": z_vals,
            "active_charts": active_charts,
            "active_roots": hit_roots
        }
        valid_points.append(pt)
        if len(valid_points) % 10 == 0:
            print(f"  Collected {len(valid_points)} points...")

    print(f"Finished. Collected {len(valid_points)} valid points.")
    with open("RESULTS/valid_points_for_fit.json", "w") as f:
        json.dump(valid_points, f, indent=2)
    print("Saved RESULTS/valid_points_for_fit.json")

if __name__ == "__main__":
    main()


