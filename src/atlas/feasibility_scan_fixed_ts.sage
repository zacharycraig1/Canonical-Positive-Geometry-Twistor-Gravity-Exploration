# src/atlas/feasibility_scan_fixed_ts.sage
#
# deterministic feasibility scan for v9 model
# Checks if *any* parameters yield pos=15 for fixed ts vectors.

import os
import json
import itertools
import numpy as np

# Use SciPy if available
try:
    from scipy.optimize import linprog
    SCIPY_OK = True
except ImportError:
    linprog = None
    SCIPY_OK = False

NPTS = 6
U_CLIP = 200.0

# ==========================================
# Helpers (duplicated from v9 for standalone)
# ==========================================

EDGE_PAIRS = [(i, j) for i in range(NPTS) for j in range(i + 1, NPTS)]
EDGE_I = np.array([i for (i, j) in EDGE_PAIRS], dtype=np.int64)
EDGE_J = np.array([j for (i, j) in EDGE_PAIRS], dtype=np.int64)

def interval_mid(ts, idx):
    if idx == 0: return float(ts[0] - 2.0)
    if idx == NPTS: return float(ts[NPTS - 1] + 2.0)
    return float(0.5 * (ts[idx - 1] + ts[idx]))

def build_base_struct(ts, pair):
    free_idx = [i for i in range(NPTS) if i not in pair]
    a, b = int(pair[0]), int(pair[1])
    t = ts.astype(np.float64)
    ta, tb = float(t[a]), float(t[b])
    denom = tb - ta
    if abs(denom) < 1e-12: denom = 1e-12

    t_free = t[free_idx]
    B = (ta - t_free) / denom
    A = -np.ones_like(B) - B

    M1 = np.zeros((NPTS, 4), dtype=np.float64)
    for k, idx in enumerate(free_idx):
        M1[idx, k] = 1.0
    M1[a, :] = A
    M1[b, :] = B

    mids = np.array([interval_mid(ts, k) for k in range(NPTS + 1)], dtype=np.float64)
    combos = []
    for ia in range(NPTS + 1):
        tx = float(mids[ia])
        for ib in range(ia + 1, NPTS + 1):
            ty = float(mids[ib])
            if tx >= ty - 1e-6: continue
            Ci = (tx - t) * (ty - t)
            s = np.sign(Ci[EDGE_I] * Ci[EDGE_J]).astype(np.int8)
            s[s == 0] = 1
            combos.append((ia, ib, tx, ty, s))

    return {
        "ts": ts, "pair": pair, "free_idx": free_idx,
        "a": a, "b": b, "t": t, "M1": M1, "combos": combos,
        "ta": ta, "tb": tb, "denom": denom, "t_free": t_free
    }

def compute_t0_for_signs(base, s_free):
    free_idx = base["free_idx"]
    a, b = base["a"], base["b"]
    ta, tb = base["ta"], base["tb"]
    denom = base["denom"]
    t_free = base["t_free"]
    
    s_free_arr = np.array(s_free, dtype=np.float64)
    S00 = float(np.sum(s_free_arr))
    S10 = float(np.sum(t_free * s_free_arr))
    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0

    t0 = np.zeros(NPTS, dtype=np.float64)
    t0[free_idx] = s_free_arr
    t0[a] = a0
    t0[b] = b0
    return t0

def build_C_matrix(base, t0):
    M1 = base["M1"]
    C = np.zeros((len(EDGE_PAIRS), 4), dtype=np.float64)
    for e, (i, j) in enumerate(EDGE_PAIRS):
        C[e, :] = t0[j] * M1[i, :] - t0[i] * M1[j, :]
    return C

def solve_lp_max_margin(C, s_edge):
    A_ub = np.zeros((C.shape[0], 5), dtype=np.float64)
    b_ub = np.zeros((C.shape[0],), dtype=np.float64)
    A_ub[:, :4] = -(s_edge[:, None] * C)
    A_ub[:, 4] = 1.0
    c = np.array([0, 0, 0, 0, -1], dtype=np.float64)
    bounds = [(-U_CLIP, U_CLIP)] * 4 + [(0, None)]

    if SCIPY_OK:
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if not res.success: return False, None, -1.0
        x = res.x
        return (x[4] > 0), x[:4], float(x[4])
    
    return False, None, -1.0

# ==========================================
# Main Scan
# ==========================================

def run_scan():
    print("Starting Deterministic Feasibility Scan (v9 model)")
    print(f"  SCIPY_OK={SCIPY_OK}, U_CLIP={U_CLIP}")
    
    # Deterministic sets of ts
    ts_sets = [
        np.array([1, 2, 3, 4, 5, 6], dtype=np.float64),
        np.array([1, 2, 4, 8, 16, 32], dtype=np.float64),
        np.array([1.1, 2.3, 3.5, 5.7, 7.9, 11.0], dtype=np.float64), # Random-ish but fixed
    ]

    all_pairs = [(i, j) for i in range(NPTS) for j in range(i + 1, NPTS)]
    sign_patterns = list(itertools.product([-1, 1], repeat=4))
    
    results = []

    for idx, ts in enumerate(ts_sets):
        print(f"\n--- Checking ts set {idx+1}: {ts} ---")
        best_m_global = -1.0
        feasible_found = False
        
        pair_stats = {}

        for pair in all_pairs:
            base = build_base_struct(ts, pair)
            best_m_pair = -1.0
            
            # Loop signs
            for s_free in sign_patterns:
                t0 = compute_t0_for_signs(base, s_free)
                C_mat = build_C_matrix(base, t0)
                
                # Loop interval combos
                for (ia, ib, _, _, s_edge) in base["combos"]:
                    ok, u, m = solve_lp_max_margin(C_mat, s_edge)
                    if m > best_m_pair:
                        best_m_pair = m
            
            pair_stats[str(pair)] = best_m_pair
            if best_m_pair > best_m_global:
                best_m_global = best_m_pair
            
            if best_m_pair > 0:
                print(f"  Pair {pair}: FEASIBLE (m={best_m_pair:.3e})")
                feasible_found = True
            # else:
            #     print(f"  Pair {pair}: infeasible (best m={best_m_pair:.3e})")

        print(f"  Result for ts[{idx}]: {'FEASIBLE' if feasible_found else 'INFEASIBLE'}")
        print(f"  Best margin overall: {best_m_global:.3e}")
        
        results.append({
            "ts": ts.tolist(),
            "feasible": feasible_found,
            "best_margin": best_m_global,
            "pair_margins": pair_stats
        })

    print("\n=== SCAN SUMMARY ===")
    for res in results:
        print(f"ts={res['ts']} -> Feasible: {res['feasible']}, Max Margin: {res['best_margin']:.3e}")
        
    with open("RESULTS/feasibility_scan_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("Saved RESULTS/feasibility_scan_results.json")

if __name__ == "__main__":
    run_scan()


