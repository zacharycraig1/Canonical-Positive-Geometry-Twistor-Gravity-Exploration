
import numpy as np
from sage.all import *
import sys

# Load the reference implementation
load("src/atlas/solve_conservation_pair.sage")

NUM_POINTS = 6

# =========================================================
# Re-implementation of Closed-Form Logic from V10 Search
# =========================================================

def build_base_struct_numpy(ts, pair):
    free_idx = [i for i in range(NUM_POINTS) if i not in pair]
    a, b = int(pair[0]), int(pair[1])
    t = ts.astype(np.float64)
    ta, tb = float(t[a]), float(t[b])
    denom = tb - ta
    if abs(denom) < 1e-12: denom = 1e-12

    t_free = t[free_idx]
    B = (ta - t_free) / denom
    A = -np.ones_like(B) - B

    M1 = np.zeros((NUM_POINTS, 4), dtype=np.float64)
    for k, idx in enumerate(free_idx):
        M1[idx, k] = 1.0
    M1[a, :] = A
    M1[b, :] = B

    return {
        "ts": ts, "pair": pair, "free_idx": free_idx,
        "a": a, "b": b, "t": t, "M1": M1, 
        "ta": ta, "tb": tb, "denom": denom, "t_free": t_free
    }

def compute_t0_from_free_numpy(base, t0_free_vals):
    free_idx = base["free_idx"]
    a, b = base["a"], base["b"]
    ta, tb = base["ta"], base["tb"]
    denom = base["denom"]
    t_free = base["t_free"]
    
    vals_arr = np.array(t0_free_vals, dtype=np.float64)
    S00 = float(np.sum(vals_arr))
    S10 = float(np.sum(t_free * vals_arr))
    
    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0

    t0 = np.zeros(NUM_POINTS, dtype=np.float64)
    t0[free_idx] = vals_arr
    t0[a] = a0
    t0[b] = b0
    return t0

def run_consistency_check(trials=100):
    print(f"Running consistency check ({trials} trials)...")
    
    max_err_t0 = 0.0
    max_err_t1 = 0.0
    
    rng = np.random.default_rng(int(42))
    
    for k in range(trials):
        # 1. Random ts (sorted)
        ts_raw = np.sort(rng.uniform(1.0, 10.0, size=NUM_POINTS))
        # Ensure gap
        for i in range(1, NUM_POINTS):
            if ts_raw[i] <= ts_raw[i-1] + 1e-4:
                ts_raw[i] = ts_raw[i-1] + 1e-3
        
        # 2. Random Pair
        a = rng.integers(0, NUM_POINTS)
        b = rng.integers(0, NUM_POINTS)
        while a == b:
            b = rng.integers(0, NUM_POINTS)
        pair = tuple(sorted((a, b)))
        
        # 3. Build Base
        base = build_base_struct_numpy(ts_raw, pair)
        
        # 4. Random free components
        # 4 free indices
        t0_free_vals = rng.uniform(-10, 10, size=4)
        u_free_vals  = rng.uniform(-50, 50, size=4) # this is t1_free
        
        # ------------------------------------------------
        # A) Closed Form Calculation
        # ------------------------------------------------
        # t0
        t0_closed = compute_t0_from_free_numpy(base, t0_free_vals)
        
        # t1 = M1 @ u
        M1 = base["M1"]
        t1_closed = M1 @ u_free_vals
        
        # ------------------------------------------------
        # B) Conservation Solver (Sage/Exact-ish)
        # ------------------------------------------------
        lambdas = {i: vector(RR, [1, float(ts_raw[i])]) for i in range(NUM_POINTS)}
        
        tildes_free = {}
        free_idx = base["free_idx"]
        for idx_enum, real_idx in enumerate(free_idx):
            val_0 = float(t0_free_vals[idx_enum])
            val_1 = float(u_free_vals[idx_enum])
            tildes_free[real_idx] = vector(RR, [val_0, val_1])
            
        tildes_solved = solve_conservation_pair(lambdas, tildes_free, pair, NUM_POINTS, RR)
        
        if tildes_solved is None:
            print(f"Trial {k}: solve_conservation_pair returned None (singularity). Skipping.")
            continue
            
        # ------------------------------------------------
        # C) Compare
        # ------------------------------------------------
        t0_err = 0.0
        t1_err = 0.0
        
        for i in range(NUM_POINTS):
            # t0 check
            diff0 = abs(t0_closed[i] - tildes_solved[i][0])
            t0_err = max(t0_err, diff0)
            
            # t1 check
            diff1 = abs(t1_closed[i] - tildes_solved[i][1])
            t1_err = max(t1_err, diff1)
            
        max_err_t0 = max(max_err_t0, t0_err)
        max_err_t1 = max(max_err_t1, t1_err)
        
        if t0_err > 1e-9 or t1_err > 1e-9:
            print(f"FAIL at trial {k}")
            print(f"  Max t0 err: {t0_err}")
            print(f"  Max t1 err: {t1_err}")
            print(f"  Pair: {pair}")
            print(f"  ts: {ts_raw}")
            sys.exit(1)
            
    print("-" * 40)
    print("CONSISTENCY CHECK PASSED")
    print(f"Max t0 Error: {max_err_t0:.3e}")
    print(f"Max t1 Error: {max_err_t1:.3e}")
    print("-" * 40)

if __name__ == "__main__":
    run_consistency_check()

