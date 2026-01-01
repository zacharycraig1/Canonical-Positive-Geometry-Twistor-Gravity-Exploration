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
TOTAL_RESTARTS      = 6000
CHUNK_RESTARTS      = 150
SEED                = 42
MAX_WORKERS         = 12   # As recommended in plan
K_TOP_CANDIDATES    = 10   # Keep top K per worker for Phase A
K_XY                = 6    # Random XY samples per interval

# LP / Sampling Constants
U_CLIP_INIT         = 50.0   # Initial Bound for u variables
Z_MARGIN_MIN        = 1e-8   # Robust threshold for m_z
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
    # Enforce strict ordering with larger gap (0.05)
    ts = np.sort(rng.uniform(1.0, 10.0, size=NUM_POINTS).astype(np.float64))
    for k in range(1, NUM_POINTS):
        if ts[k] <= ts[k - 1] + 0.05:
            ts[k] = ts[k - 1] + 0.05 + rng.uniform(0, 0.01)
    return ts

def interval_mid(ts, idx):
    if idx == 0: return float(ts[0] - 2.0)
    if idx == NUM_POINTS: return float(ts[NUM_POINTS - 1] + 2.0)
    return float(0.5 * (ts[idx - 1] + ts[idx]))

def get_interval_bounds(ts, idx, Lx=None, Ly=None):
    if idx == 0:
        # (-inf, ts[0])
        U = float(ts[0])
        # Use Lx to define effective L
        L = U - (Lx if Lx else 2.0)
        return L, U, True # True = "infinite" end (conceptually)
    elif idx == NUM_POINTS:
        # (ts[-1], +inf)
        L = float(ts[NUM_POINTS - 1])
        U = L + (Ly if Ly else 2.0)
        return L, U, True
    else:
        # (ts[idx-1], ts[idx])
        return float(ts[idx-1]), float(ts[idx]), False

def sample_in_interval(rng, L, U, is_outer):
    if is_outer:
        # Sample roughly Uniform(0.5, 3.0) away from boundary
        dist = rng.uniform(0.5, 3.0)
        # If it's the left outer (idx=0), we go L = U - dist, but here L, U are bounds.
        # Wait, get_interval_bounds returns effective L, U.
        # For idx=0, U=ts[0]. We want range (-inf, ts[0]).
        # For idx=N, L=ts[-1]. We want range (ts[-1], inf).
        pass
    
    # Actually, simpler logic:
    # If finite: L + 0.15*W ... U - 0.15*W
    width = U - L
    if not is_outer:
        low = L + 0.15 * width
        high = U - 0.15 * width
        if low >= high: return (L + U) * 0.5
        return rng.uniform(low, high)
    else:
        # For outer, we used Lx/Ly to define a range, but let's stick to the plan's logic:
        # "sample U - Uniform(0.5, 3.0) or L + Uniform(0.5, 3.0)"
        # My get_interval_bounds logic was mixing things up.
        # Let's handle it at call site or fix here.
        pass
    return (L + U) * 0.5 # fallback

def build_base_struct_v11(ts, pair, rng, Lx=None, Ly=None):
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

    combos = []
    
    # Pre-calculate interval bounds
    intervals = []
    for k in range(NUM_POINTS + 1):
        if k == 0:
            intervals.append((-float('inf'), float(ts[0])))
        elif k == NUM_POINTS:
            intervals.append((float(ts[NUM_POINTS-1]), float('inf')))
        else:
            intervals.append((float(ts[k-1]), float(ts[k])))

    for ia in range(NUM_POINTS + 1):
        for ib in range(ia + 1, NUM_POINTS + 1):
            # Sample K_XY points in this interval pair
            for _ in range(K_XY):
                # Sample tx
                La, Ua = intervals[ia]
                if ia == 0:
                    tx = Ua - rng.uniform(0.5, 3.0)
                elif ia == NUM_POINTS: # Should not happen if ib > ia
                    tx = La + rng.uniform(0.5, 3.0)
                else:
                    wa = Ua - La
                    tx = rng.uniform(La + 0.15*wa, Ua - 0.15*wa)
                
                # Sample ty
                Lb, Ub = intervals[ib]
                if ib == 0: # Should not happen
                    ty = Ub - rng.uniform(0.5, 3.0)
                elif ib == NUM_POINTS:
                    ty = Lb + rng.uniform(0.5, 3.0)
                else:
                    wb = Ub - Lb
                    ty = rng.uniform(Lb + 0.15*wb, Ub - 0.15*wb)
                
                if tx >= ty - 1e-6: continue

                Ci = (tx - t) * (ty - t)
                # Compute signs
                s = np.sign(Ci[EDGE_I] * Ci[EDGE_J]).astype(np.int8)
                s[s == 0] = 1
                
                # Compute scaling factors for Z-Margin
                # scale_e = abs(Cprod_e) * abs(ang_inv_e)
                # ang_inv_e = 1/abs(t[j]-t[i])
                ang_diff = np.abs(t[EDGE_J] - t[EDGE_I])
                ang_diff = np.where(ang_diff < 1e-12, 1e-12, ang_diff)
                ang_inv = 1.0 / ang_diff
                
                Cprod = np.abs(Ci[EDGE_I] * Ci[EDGE_J])
                scale = Cprod * ang_inv
                
                combos.append({
                    "ia": ia, "ib": ib,
                    "tx": tx, "ty": ty,
                    "s_edge": s,
                    "scale": scale
                })

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

def solve_lp_z_margin(C, s_edge, scale, u_clip_val):
    # Maximize m_z
    # Constraint: scale_e * (s_e * (C_e . u)) >= m_z
    # => -(scale_e * s_e * C_e) . u + m_z <= 0
    
    n_edges = C.shape[0]
    A_ub = np.zeros((n_edges, 5), dtype=np.float64)
    b_ub = np.zeros((n_edges,), dtype=np.float64)
    
    # A_ub[:, :4] = -(scale[:, None] * s_edge[:, None] * C)
    # Using broadcasting carefully
    factor = -(scale * s_edge)[:, None] # Shape (15, 1)
    A_ub[:, :4] = factor * C
    A_ub[:, 4] = 1.0 # coeff for m_z
    
    # Objective: maximize m_z => minimize -m_z
    c = np.array([0, 0, 0, 0, -1], dtype=np.float64)
    
    # Bounds: u in [-clip, clip], m_z in [0, None]
    bounds = [(-u_clip_val, u_clip_val)] * 4 + [(0, None)]

    if SCIPY_OK:
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if not res.success: return False, None, -1.0
        x = res.x
        return (x[4] > 0), x[:4], float(x[4])
    return False, None, -1.0

def compute_z_values(base, t0, u, tx, ty):
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
    best_heap = [] # Store tuples: (pos, zmin, tie, cand_dict)
    
    rejects_log = []
    
    def add_to_heap(cand):
        tie = rng.random()
        # Sort by zmin (larger is better) for top candidates? 
        # Actually heap is min-heap. If we want to keep TOP K candidates (largest zmin),
        # we push (zmin, tie, cand). Then heappushpop removes the smallest zmin.
        # Wait, if pos is not always 15, we prioritize pos.
        # So key is (pos, zmin).
        item = (cand["pos"], cand["zmin"], tie, cand)
        if len(best_heap) < K_TOP_CANDIDATES:
            heapq.heappush(best_heap, item)
        else:
            heapq.heappushpop(best_heap, item)

    for _ in range(restarts):
        ts = sample_ts(rng)
        pair = pick_pair(rng)
        
        # Lx/Ly are used for "outer" interval definition in logs, 
        # but in v11 sampling they are implicit or dynamic. 
        # We'll just generate them for record keeping if needed.
        Lx = float(np.exp(rng.uniform(np.log(2.0), np.log(500.0))))
        Ly = float(np.exp(rng.uniform(np.log(2.0), np.log(500.0))))
        
        base = build_base_struct_v11(ts, pair, rng, Lx, Ly)

        # Multi-T0 loop
        for _t0_iter in range(K_T0_TRIALS):
            # Sample free t0: sign * magnitude
            mags = rng.choice(T0_MAGS, size=4)
            signs = rng.choice([-1, 1], size=4)
            t0_free_vals = mags * signs

            t0 = compute_t0_from_free(base, t0_free_vals)
            C_mat = build_C_matrix(base, t0)
            
            # Iterate over all sampled (tx, ty) combos
            for combo in base["combos"]:
                tx = combo["tx"]
                ty = combo["ty"]
                s_edge = combo["s_edge"]
                scale = combo["scale"]
                
                # Z-Margin LP
                u_clip = U_CLIP_INIT
                ok, u, mz = solve_lp_z_margin(C_mat, s_edge, scale, u_clip)
                
                if not ok:
                    u_clip = 200.0
                    ok, u, mz = solve_lp_z_margin(C_mat, s_edge, scale, u_clip)
                
                if not ok or mz is None:
                    continue

                if mz < Z_MARGIN_MIN:
                    continue
                
                # Compute actual z
                z, _, _ = compute_z_values(base, t0, u, tx, ty)
                pos = int(np.count_nonzero(z > 0))
                zmin = float(np.min(z))
                
                # We want robust points. 
                # If pos=15, we check hi-prec.
                if pos == 15:
                    cand = {
                        "ts": ts.tolist(),
                        "pair": [int(pair[0]), int(pair[1])],
                        "free_idx": [int(i) for i in base["free_idx"]],
                        "t0_free_vals": t0_free_vals.tolist(),
                        "u_free": u.tolist(),
                        "tx": float(tx),
                        "ty": float(ty),
                        "mz_float": float(mz),
                        "zmin_float": float(zmin),
                        "pos_float": pos,
                        "ia": combo["ia"],
                        "ib": combo["ib"],
                        "Lx": Lx, "Ly": Ly
                    }
                    
                    # Hi-Prec Check
                    # We can't run full Sage here efficiently inside the worker if it spawns process.
                    # But we can assume verify_candidate_hiprec is available if we structure correctly.
                    # Actually, importing sage inside worker is slow.
                    # The plan says: "Run only high-precision reconstruction (check_point)."
                    # It's better to do this in the main process or have a sage-enabled worker.
                    # However, "check_point" relies on Sage's RealField.
                    # If this script is run with `sage -python` or `sage`, workers can use Sage.
                    # But multiprocessing with Sage can be tricky.
                    # Strategy: Use a helper function that imports sage locally.
                    
                    hp_res = check_point_hiprec(cand)
                    cand["pos_hi"] = hp_res["pos"]
                    cand["zmin_hi"] = hp_res["zmin"]
                    cand["neg_edges_hi"] = hp_res["neg_edges"]

                    if cand["pos_hi"] == 15:
                        cand["pos"] = 15
                        cand["zmin"] = cand["zmin_hi"]
                        add_to_heap(cand)
                    else:
                        # Reject log
                        rejects_log.append(cand)
                
    return {"heap": best_heap, "rejects": rejects_log}

def check_point_hiprec(cand):
    # Minimal Sage check
    try:
        from sage.all import RealField, vector
        RF = RealField(200)
    except ImportError:
        return {"pos": -1, "zmin": -1.0, "neg_edges": []}

    ts = [float(x) for x in cand["ts"]]
    pair = (int(cand["pair"][0]), int(cand["pair"][1]))
    free_idx = [int(i) for i in cand["free_idx"]]
    t0_free_vals = cand["t0_free_vals"]
    u_free = cand["u_free"]
    tx = float(cand["tx"])
    ty = float(cand["ty"])
    
    lambdas = {i: vector(RF, [1, ts[i]]) for i in range(NUM_POINTS)}
    
    # Reconstruct t0, t1 from free params
    # We need to solve conservation for tildes.
    # We can implement a pure-python or sage reconstruction of tildes.
    # Let's reuse the logic from `solve_conservation_pair.sage` effectively reimplemented here
    # or rely on `solve_conservation_pair` if we can import it.
    # To keep it self-contained and fast:
    
    ta, tb = ts[pair[0]], ts[pair[1]]
    denom = tb - ta
    
    # t0 reconstruction
    S00 = sum(t0_free_vals)
    S10 = sum(ts[k] * v for k, v in zip(free_idx, t0_free_vals))
    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0
    t0_vals = [0]*NUM_POINTS
    for k, idx in enumerate(free_idx): t0_vals[idx] = t0_free_vals[k]
    t0_vals[pair[0]] = a0
    t0_vals[pair[1]] = b0
    
    # t1 reconstruction
    S01 = sum(u_free)
    S11 = sum(ts[k] * v for k, v in zip(free_idx, u_free))
    b1 = (-S11 + ta * S01) / denom
    a1 = -S01 - b1
    t1_vals = [0]*NUM_POINTS
    for k, idx in enumerate(free_idx): t1_vals[idx] = u_free[k]
    t1_vals[pair[0]] = a1
    t1_vals[pair[1]] = b1
    
    # Calculate Z
    z_min = 1000.0
    pos_count = 0
    neg_edges = []
    
    xs = vector(RF, [1, tx])
    ys = vector(RF, [1, ty])
    
    def br(u, v): return u[0]*v[1] - u[1]*v[0]
    
    tildes = {}
    for i in range(NUM_POINTS):
        tildes[i] = vector(RF, [t0_vals[i], t1_vals[i]])
        
    Cp = {i: br(lambdas[i], xs) * br(lambdas[i], ys) for i in range(NUM_POINTS)}
    
    idx_list = list(range(NUM_POINTS))
    edge_idx = 0
    for i in range(NUM_POINTS):
        for j in range(i+1, NUM_POINTS):
            ang = br(lambdas[i], lambdas[j])
            sq = br(tildes[j], tildes[i])
            val = (sq / ang) * Cp[i] * Cp[j]
            zf = float(val)
            if zf > 0:
                pos_count += 1
            else:
                neg_edges.append(edge_idx)
            z_min = min(z_min, zf)
            edge_idx += 1
            
    return {"pos": pos_count, "zmin": z_min, "neg_edges": neg_edges}


def main():
    print("Starting V11 Z-Margin Robust Search")
    print(f"  RESTARTS={TOTAL_RESTARTS}, CHUNK={CHUNK_RESTARTS}")
    print(f"  Z_MARGIN_MIN={Z_MARGIN_MIN}, K_XY={K_XY}")
    
    cores = os.cpu_count() or 1
    if MAX_WORKERS is None:
        workers = min(max(1, cores - 2), 12)
    else:
        workers = min(MAX_WORKERS, cores)

    num_chunks = int(math.ceil(TOTAL_RESTARTS / float(CHUNK_RESTARTS)))
    print(f"  workers={workers}, total_chunks={num_chunks}")

    started = time.time()
    
    global_heap = []
    all_rejects = []

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = []
        for k in range(num_chunks):
            restarts = CHUNK_RESTARTS if (k < num_chunks - 1) else (TOTAL_RESTARTS - CHUNK_RESTARTS * (num_chunks - 1))
            futures.append(ex.submit(worker_chunk, k, restarts, SEED + 1000 * k))

        done = 0
        
        for fut in as_completed(futures):
            done += 1
            res = fut.result()
            
            # Merge heaps
            w_heap = res.get("heap", [])
            for item in w_heap:
                if len(global_heap) < K_TOP_CANDIDATES * 4: # Keep a larger pool globally
                    heapq.heappush(global_heap, item)
                else:
                    heapq.heappushpop(global_heap, item)
            
            # Collect rejects
            all_rejects.extend(res.get("rejects", []))
            
            if done % 5 == 0:
                best_str = "None"
                if global_heap:
                    # Best is max in heap (since we push (pos, zmin, ...))
                    # But heapq is a min-heap. So global_heap[0] is the WORST of the set.
                    # The best is the one with max (pos, zmin).
                    b = max(global_heap, key=lambda x: (x[0], x[1]))
                    best_str = f"POS={b[0]}/15 zmin={b[1]:.2e} mz={b[3]['mz_float']:.2e}"
                
                print(f"[status] chunks={done}/{num_chunks}, best: {best_str}, rejects={len(all_rejects)}, elapsed={time.time()-started:.1f}s")
    
    print("Search Complete.")
    
    # Save Rejects
    os.makedirs("RESULTS", exist_ok=True)
    with open("RESULTS/v11_rejects_pos15.jsonl", "w") as f:
        for r in all_rejects:
            f.write(json.dumps(r) + "\n")
    print(f"Saved {len(all_rejects)} rejected candidates to RESULTS/v11_rejects_pos15.jsonl")
    
    # Save Top Candidates
    # Sort descending
    sorted_cands = sorted(global_heap, key=lambda x: (x[0], x[1]), reverse=True)
    # Deduplicate based on zmin/tx/ty somewhat? Or just keep all.
    # Keep top 100
    top_cands = [x[3] for x in sorted_cands[:100]]
    
    with open("RESULTS/v11_top_candidates.json", "w") as f:
        json.dump(top_cands, f, indent=2)
    print(f"Saved {len(top_cands)} top candidates to RESULTS/v11_top_candidates.json")

if __name__ == "__main__":
    main()

