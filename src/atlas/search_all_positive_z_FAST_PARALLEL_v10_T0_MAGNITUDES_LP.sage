import os
# Set thread limits to avoid oversubscription
for _k in [
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"
]:
    os.environ.setdefault(_k, "1")

import json, time, math, heapq
import itertools
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

# =========================
# CONFIG
# =========================
NUM_POINTS = 6

# Search Parameters
TOTAL_RESTARTS      = 12000
CHUNK_RESTARTS      = 150
SEED                = 42
MAX_WORKERS         = None   # auto
K_TOP_CANDIDATES    = 25     # Keep top K candidates per worker

# LP / Sampling Constants
U_CLIP_INIT         = 50.0   # Initial Bound for u variables
MARGIN_MIN          = 0.0    # Relaxed for Phase A
ZMIN_MIN            = 0.0    # Relaxed for Phase A
T0_MAGS             = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]

# Inner loop parameters
K_T0_TRIALS         = 64     # T0 patterns per (ts, pair)

# Pair bias (same as v8/v9)
PAIR_BIAS = {
    (0, 1): 4.0,
    (4, 5): 4.0,
    (0, 5): 2.0,
    (1, 4): 2.0,
}

EDGE_PAIRS = [(i, j) for i in range(NUM_POINTS) for j in range(i + 1, NUM_POINTS)]
EDGE_I = np.array([i for (i, j) in EDGE_PAIRS], dtype=np.int64)
EDGE_J = np.array([j for (i, j) in EDGE_PAIRS], dtype=np.int64)
ALL_ROOTS = list(itertools.combinations(range(NUM_POINTS), 3))

# SciPy LP
try:
    from scipy.optimize import linprog
    SCIPY_OK = True
except ImportError:
    linprog = None
    SCIPY_OK = False

def _pair_weights():
    pairs = [(i, j) for i in range(NUM_POINTS) for j in range(i + 1, NUM_POINTS)]
    w = np.array([PAIR_BIAS.get(p, 1.0) for p in pairs], dtype=np.float64)
    w /= np.sum(w)
    return pairs, w

PAIRS, PAIR_W = _pair_weights()

def pick_pair(rng):
    idx = rng.choice(len(PAIRS), p=PAIR_W)
    return PAIRS[idx]

def sample_ts(rng):
    # Enforce strict ordering with gap
    ts = np.sort(rng.uniform(1.0, 10.0, size=NUM_POINTS).astype(np.float64))
    for k in range(1, NUM_POINTS):
        if ts[k] <= ts[k - 1] + 1e-6:
            ts[k] = ts[k - 1] + 1e-5
    return ts

def interval_mid(ts, idx):
    if idx == 0: return float(ts[0] - 2.0)
    if idx == NUM_POINTS: return float(ts[NUM_POINTS - 1] + 2.0)
    return float(0.5 * (ts[idx - 1] + ts[idx]))

def build_base_struct(ts, pair, Lx=None, Ly=None):
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

    # XY Sampling Logic
    mids = np.array([interval_mid(ts, k) for k in range(NUM_POINTS + 1)], dtype=np.float64)
    
    # Override outer intervals if Lx/Ly provided
    if Lx is not None:
        mids[0] = float(ts[0] - Lx)
    if Ly is not None:
        mids[NUM_POINTS] = float(ts[NUM_POINTS-1] + Ly)

    combos = []
    for ia in range(NUM_POINTS + 1):
        tx = float(mids[ia])
        for ib in range(ia + 1, NUM_POINTS + 1):
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

def compute_t0_from_free(base, t0_free_vals):
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

def build_C_matrix(base, t0):
    M1 = base["M1"]
    C = np.zeros((len(EDGE_PAIRS), 4), dtype=np.float64)
    for e, (i, j) in enumerate(EDGE_PAIRS):
        C[e, :] = t0[j] * M1[i, :] - t0[i] * M1[j, :]
    return C

def solve_lp_max_margin(C, s_edge, u_clip_val):
    A_ub = np.zeros((C.shape[0], 5), dtype=np.float64)
    b_ub = np.zeros((C.shape[0],), dtype=np.float64)
    A_ub[:, :4] = -(s_edge[:, None] * C)
    A_ub[:, 4] = 1.0
    c = np.array([0, 0, 0, 0, -1], dtype=np.float64)
    bounds = [(-u_clip_val, u_clip_val)] * 4 + [(0, None)]

    if SCIPY_OK:
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if not res.success: return False, None, -1.0
        x = res.x
        return (x[4] > 0), x[:4], float(x[4])
    return False, None, -1.0

def compute_z_v10(base, t0, u, tx, ty):
    M1 = base["M1"]
    t = base["t"]
    t1 = M1 @ u

    ang = (t[EDGE_J] - t[EDGE_I])
    ang = np.where(np.abs(ang) < 1e-12, 1e-12, ang)
    ang_inv = 1.0 / ang

    C_poly = (tx - t) * (ty - t)
    Cprod = (C_poly[EDGE_I] * C_poly[EDGE_J])

    br = t0[EDGE_J] * t1[EDGE_I] - t1[EDGE_J] * t0[EDGE_I]
    z = (br * ang_inv) * Cprod
    return z, t0, t1

# =========================
# Worker Logic
# =========================

def worker_chunk(task_id, restarts, seed):
    rng = np.random.default_rng(int(seed))
    best_heap = [] # Store tuples: (pos, zmin, cand_dict)
    
    # Helper to push to heap
    def add_to_heap(cand):
        # Add a random tiebreaker to avoid comparing dicts
        tie = rng.random()
        item = (cand["pos"], cand["zmin"], tie, cand)
        if len(best_heap) < K_TOP_CANDIDATES:
            heapq.heappush(best_heap, item)
        else:
            heapq.heappushpop(best_heap, item)

    for _ in range(restarts):
        ts = sample_ts(rng)
        pair = pick_pair(rng)
        
        # Sample XY distances (LogUniform)
        # loguniform(a, b) -> exp(uniform(log(a), log(b)))
        Lx = float(np.exp(rng.uniform(np.log(2.0), np.log(500.0))))
        Ly = float(np.exp(rng.uniform(np.log(2.0), np.log(500.0))))
        
        base = build_base_struct(ts, pair, Lx=Lx, Ly=Ly)

        # Multi-T0 loop
        for _t0_iter in range(K_T0_TRIALS):
            # Sample free t0: sign * magnitude
            mags = rng.choice(T0_MAGS, size=4)
            signs = rng.choice([-1, 1], size=4)
            t0_free_vals = mags * signs

            t0 = compute_t0_from_free(base, t0_free_vals)
            C_mat = build_C_matrix(base, t0)
            
            # Iterate over all interval combos
            for (ia, ib, tx_mid, ty_mid, s_edge) in base["combos"]:
                
                # Adaptive U_CLIP loop
                u_clip = U_CLIP_INIT
                ok, u, m = solve_lp_max_margin(C_mat, s_edge, u_clip)
                
                if not ok:
                    # Retry once with larger U_CLIP
                    u_clip = 200.0
                    ok, u, m = solve_lp_max_margin(C_mat, s_edge, u_clip)
                
                if not ok or m is None:
                    continue

                if m <= 1e-6: # Margin tiny
                     u_clip = 800.0
                     ok_retry, u_retry, m_retry = solve_lp_max_margin(C_mat, s_edge, u_clip)
                     if ok_retry and m_retry > m:
                         u, m = u_retry, m_retry
                
                if m <= MARGIN_MIN:
                    continue
                
                # Verify solution
                z, _, _ = compute_z_v10(base, t0, u, tx_mid, ty_mid)
                pos = int(np.count_nonzero(z > 0))
                zmin = float(np.min(z))

                if zmin < ZMIN_MIN:
                    continue
                
                cand = {
                    "ts": ts.tolist(),
                    "pair": [int(pair[0]), int(pair[1])],
                    "free_idx": [int(i) for i in base["free_idx"]],
                    "t0_free_vals": t0_free_vals.tolist(),
                    "u_free": u.tolist(),
                    "tx": float(tx_mid),
                    "ty": float(ty_mid),
                    "margin": float(m),
                    "zmin": float(zmin),
                    "pos": pos,
                    "ia": ia,
                    "ib": ib,
                    "Lx": Lx, 
                    "Ly": Ly
                }
                
                # Capture in heap
                add_to_heap(cand)

                if pos == 15:
                     # Return immediately if found
                     return {"found": True, "candidate": cand, "heap": best_heap}

    return {"found": False, "candidate": None, "heap": best_heap}

# =========================
# Golden Validation
# =========================

def check_active_charts_and_save(candidate):
    """
    Heavy Sage checks: reconstruct, verify pos=15, check charts.
    Includes Rescue Strategies.
    """
    from sage.all import RealField, vector, load
    
    load("src/atlas/jacobian_fiber_gauge.sage")
    load("src/atlas/solve_conservation_pair.sage") 

    RF = RealField(200) 
    evaluator = FiberGaugeEvaluator(NUM_POINTS)

    ts = [float(x) for x in candidate["ts"]]
    pair = (int(candidate["pair"][0]), int(candidate["pair"][1]))
    free_idx = [int(i) for i in candidate["free_idx"]]
    t0_free_vals = candidate["t0_free_vals"]
    u_free_init = candidate["u_free"]
    
    lambdas_sage = {i: vector(RF, [1, ts[i]]) for i in range(NUM_POINTS)}
    tx_init = float(candidate["tx"])
    ty_init = float(candidate["ty"])
    
    def br(u, v):
        return u[0]*v[1] - u[1]*v[0]
    
    def check_point(u_curr, tx_curr, ty_curr):
        xs_c = vector(RF, [1, tx_curr])
        ys_c = vector(RF, [1, ty_curr])
        Cp = {i: br(lambdas_sage[i], xs_c) * br(lambdas_sage[i], ys_c) for i in range(NUM_POINTS)}

        tildes_free_sage = {}
        for k, idx in enumerate(free_idx):
            tildes_free_sage[idx] = vector(RF, [float(t0_free_vals[k]), float(u_curr[k])])
            
        tildes_recon = solve_conservation_pair(lambdas_sage, tildes_free_sage, pair, NUM_POINTS, RF)
        if tildes_recon is None:
            return 0, -1.0, [], None
            
        pos_count = 0
        zmin_val = None
        z_out = []
        for (i, j) in EDGE_PAIRS:
            ang = br(lambdas_sage[i], lambdas_sage[j])
            sq  = br(tildes_recon[j], tildes_recon[i])
            val = (sq / ang) * Cp[i] * Cp[j]
            zf = float(val)
            z_out.append(zf)
            if zf > 0:
                pos_count += 1
            zmin_val = zf if zmin_val is None else min(zmin_val, zf)
        return pos_count, zmin_val, z_out, tildes_recon

    # Initial check
    pos_count, zmin, z_vals, _ = check_point(u_free_init, tx_init, ty_init)
    
    print(f"  Reconstructed Positivity (hi-prec): {pos_count}/15  (zmin={zmin:.3e})")
    
    u_final = list(u_free_init)
    tx_final = tx_init
    ty_final = ty_init
    
    rescued = (pos_count == 15 and zmin > 1e-16) # use stricter zmin for golden
    
    if not rescued:
        print("  [rescue] Attempting rescue sequence...")
        
        # Strategy 0: LP Re-solve with Hi-Prec Signs
        # Reconstruct base to get C matrix
        try:
            Lx = candidate.get("Lx", None)
            Ly = candidate.get("Ly", None)
            base_recon = build_base_struct(np.array(ts), pair, Lx=Lx, Ly=Ly)
            t0_recon = compute_t0_from_free(base_recon, t0_free_vals)
            C_recon = build_C_matrix(base_recon, t0_recon)
            
            # Compute signs from hi-prec geometry
            # We need the signs of Ci*Cj at (tx, ty)
            # We already computed Cp in check_point, but didn't return it.
            # Let's recompute Cp quickly.
            xs_c = vector(RF, [1, tx_init])
            ys_c = vector(RF, [1, ty_init])
            Cp = {i: br(lambdas_sage[i], xs_c) * br(lambdas_sage[i], ys_c) for i in range(NUM_POINTS)}
            
            s_edge_new = np.zeros(len(EDGE_PAIRS), dtype=np.int8)
            for k, (i, j) in enumerate(EDGE_PAIRS):
                val = float(Cp[i] * Cp[j])
                s = 1 if val >= 0 else -1
                s_edge_new[k] = s
            
            # Resolve LP
            print("  [rescue] Re-solving LP with hi-prec signs...")
            ok_lp, u_lp, m_lp = solve_lp_max_margin(C_recon, s_edge_new, 800.0) # Use loose bound
            if ok_lp and m_lp > 0:
                 p_lp, zmin_lp, z_lp, _ = check_point(u_lp, tx_init, ty_init)
                 if p_lp > pos_count:
                     print(f"  [rescue] LP Re-solve improved pos: {pos_count} -> {p_lp}")
                     u_free_init = u_lp # Update starting point for grid sweep
                     pos_count = p_lp
                     zmin = zmin_lp
                     z_vals = z_lp
                     if p_lp == 15 and zmin_lp > 1e-16:
                         u_final = list(u_lp)
                         rescued = True
                         print(f"  [rescue] LP Re-solve SUCCESS (zmin={zmin_lp:.3e})")
        except Exception as e:
            print(f"  [rescue] LP Re-solve failed with error: {e}")

    if not rescued:
        # Strategy: Deterministic Grid Sweep
        # 8x8 Grid sweep around tx, ty
        steps = np.linspace(-0.05, 0.05, 9)
        best_r_pos = pos_count
        best_r_zmin = zmin
        
        for dx in steps:
            for dy in steps:
                if dx==0 and dy==0: continue
                p, zm, zv, _ = check_point(u_free_init, tx_init + dx, ty_init + dy)
                if p > best_r_pos or (p == best_r_pos and zm > best_r_zmin):
                    best_r_pos = p
                    best_r_zmin = zm
                    if p == 15 and zm > 1e-16:
                        tx_final = tx_init + dx
                        ty_final = ty_init + dy
                        pos_count = 15
                        zmin = zm
                        z_vals = zv
                        rescued = True
                        print(f"  [rescue] Grid success: dx={dx:.3f}, dy={dy:.3f}, zmin={zm:.3e}")
                        break
            if rescued: break
            
    if not rescued:
         print("  [reject] Rescue failed.")
         return False

    # Chart check
    tf_final = {int(i): [float(t0_free_vals[k]), float(u_final[k])] for k, i in enumerate(free_idx)}
    xs_final = [1.0, float(tx_final)]
    ys_final = [1.0, float(ty_final)]

    active_charts = 0
    hit_roots = []
    
    for roots in ALL_ROOTS:
        val = evaluator.evaluate_chart_fiber_gauge(roots, ts, tf_final, xs_final, ys_final, require_polytope=True)
        if abs(val) > 1e-20:
            active_charts += 1
            hit_roots.append(roots)

    if active_charts <= 0:
        print("  [reject] 15/15 but no active charts under polytope requirement.")
        return False

    os.makedirs("RESULTS", exist_ok=True)
    out = dict(candidate)
    out["u_free"] = u_final
    out["tx"] = tx_final
    out["ty"] = ty_final
    out["z_values"] = z_vals
    out["active_charts"] = int(active_charts)
    out["active_roots"] = [list(r) for r in hit_roots]
    out["meta"] = {
        "TOTAL_RESTARTS": int(TOTAL_RESTARTS),
        "CHUNK_RESTARTS": int(CHUNK_RESTARTS),
        "U_CLIP": float(U_CLIP_INIT),
        "timestamp": float(time.time())
    }

    with open("RESULTS/golden_kinematics.json", "w") as f:
        json.dump(out, f, indent=2)

    print(f"  ✨ GOLDEN FIND! Active Charts: {active_charts}")
    print("  Saved: RESULTS/golden_kinematics.json")
    return True

def main():
    print("Starting V10 Improved Search")
    print(f"  RESTARTS={TOTAL_RESTARTS}, CHUNK={CHUNK_RESTARTS}")
    print(f"  K_T0_TRIALS={K_T0_TRIALS}, T0_MAGS={T0_MAGS}")
    print(f"  U_CLIP_INIT={U_CLIP_INIT}")
    
    cores = os.cpu_count() or 1
    if MAX_WORKERS is None:
        workers = min(max(1, cores - 2), 12)
    else:
        workers = min(MAX_WORKERS, cores)

    num_chunks = int(math.ceil(TOTAL_RESTARTS / float(CHUNK_RESTARTS)))
    print(f"  workers={workers}, total_chunks={num_chunks}")

    started = time.time()
    
    # Global Heap: keep top K overall
    global_heap = []

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = []
        for k in range(num_chunks):
            restarts = CHUNK_RESTARTS if (k < num_chunks - 1) else (TOTAL_RESTARTS - CHUNK_RESTARTS * (num_chunks - 1))
            futures.append(ex.submit(worker_chunk, k, restarts, SEED + 1000 * k))

        done = 0
        found_golden = False
        
        for fut in as_completed(futures):
            done += 1
            res = fut.result()
            
            # Merge heaps
            w_heap = res.get("heap", [])
            for item in w_heap:
                if len(global_heap) < K_TOP_CANDIDATES:
                    heapq.heappush(global_heap, item)
                else:
                    heapq.heappushpop(global_heap, item)
            
            if res["found"] and not found_golden:
                cand = res["candidate"]
                print(f"[candidate] Found POS={cand['pos']} zmin={cand['zmin']:.3e}. Checking charts...")
                ok = check_active_charts_and_save(cand)
                if ok:
                    print("Golden point confirmed. Stopping search.")
                    found_golden = True
                    for f in futures: f.cancel()
                    break
            
            if done % 5 == 0:
                best_str = "None"
                if global_heap:
                    # heap is min-heap of top K. Max element is the "best" of the top K?
                    # No, heap items are (pos, zmin, cand). Smallest element is at [0].
                    # So global_heap[0] is the worst of the top K.
                    # The best of the top K is max(global_heap).
                    b = max(global_heap, key=lambda x: (x[0], x[1]))
                    best_str = f"POS={b[0]}/15 zmin={b[1]:.2e}"
                print(f"[status] chunks={done}/{num_chunks}, best_in_heap: {best_str}, elapsed={time.time()-started:.1f}s")
    
    print("Saving best candidates...")
    os.makedirs("RESULTS", exist_ok=True)
    
    sorted_cands = sorted(global_heap, key=lambda x: (x[0], x[1]), reverse=True)
    out_cands = [x[3] for x in sorted_cands]
    
    with open("RESULTS/v10_best_candidates.json", "w") as f:
        json.dump(out_cands, f, indent=2)
        
    meta = {
        "total_restarts": TOTAL_RESTARTS,
        "elapsed": time.time() - started,
        "top_pos": out_cands[0]["pos"] if out_cands else 0,
        "top_zmin": out_cands[0]["zmin"] if out_cands else 0
    }
    with open("RESULTS/v10_run_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    print("Done.")

if __name__ == "__main__":
    main()
