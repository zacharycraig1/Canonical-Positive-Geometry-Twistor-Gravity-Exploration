# src/atlas/search_all_positive_z_FAST_PARALLEL_v9_SIGNED_T0.sage
#
# v9: Extends v8 by allowing signed free tilde0 components.
#
# Core idea:
#   In v8, free tilde0 components were fixed to +1.
#   In v9, we iterate over all 16 sign patterns s_k \in {+1, -1} for the 4 free indices.
#   This retains linearity of the LP in u (since s_k are fixed per inner loop),
#   but significantly expands the gauge freedom to find feasible regions.
#
# Run:
#   sage src/atlas/search_all_positive_z_FAST_PARALLEL_v9_SIGNED_T0.sage

import os
for _k in [
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"
]:
    os.environ.setdefault(_k, "1")

import json, time, math
import itertools
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

# =========================
# CONFIG
# =========================
NUM_POINTS = 6

TOTAL_RESTARTS      = 500
CHUNK_RESTARTS      = 50
SEED                = 42
MAX_WORKERS         = None   # auto

U_CLIP              = 50.0   # Reduced to improve numerical stability
MARGIN_MIN          = 1e-8   # Require robust margin
ZMIN_MIN            = 1e-10  # Require robust zmin

# (tx,ty) intervals
REFINE_XY_TRIES     = 14

# Pair bias (same as v8)
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
except Exception:
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
    ts = np.sort(rng.uniform(1.0, 10.0, size=NUM_POINTS).astype(np.float64))
    for k in range(1, NUM_POINTS):
        if ts[k] <= ts[k - 1]:
            ts[k] = ts[k - 1] + 1e-6
    return ts


def interval_mid(ts, idx):
    if idx == 0:
        return float(ts[0] - 2.0)
    if idx == NUM_POINTS:
        return float(ts[NUM_POINTS - 1] + 2.0)
    return float(0.5 * (ts[idx - 1] + ts[idx]))


def interval_bounds(ts, idx):
    if idx == 0:
        return (-np.inf, float(ts[0]))
    if idx == NUM_POINTS:
        return (float(ts[NUM_POINTS - 1]), np.inf)
    return (float(ts[idx - 1]), float(ts[idx]))


def sample_tx_ty_in_intervals(ts, ia, ib, rng):
    lo_x, hi_x = interval_bounds(ts, ia)
    lo_y, hi_y = interval_bounds(ts, ib)

    def sample_in(lo, hi):
        if not np.isfinite(lo) and np.isfinite(hi):
            return float(hi - rng.uniform(0.5, 3.0))
        if np.isfinite(lo) and not np.isfinite(hi):
            return float(lo + rng.uniform(0.5, 3.0))
        if hi - lo < 1e-9:
            return float(lo)
        return float(rng.uniform(lo, hi))

    tx = sample_in(lo_x, hi_x)
    ty = sample_in(lo_y, hi_y)
    if tx >= ty - 1e-6:
        m = 0.5 * (tx + ty)
        tx = m - 1e-3
        ty = m + 1e-3
    return tx, ty


def build_base_struct(ts, pair):
    """
    Precompute parts of the structure that do NOT depend on the signs of free tilde0.
    - M1 matrix (linear map from u to t1)
    - Interval combos (sign(Ci Cj))
    """
    free_idx = [i for i in range(NUM_POINTS) if i not in pair]
    a, b = int(pair[0]), int(pair[1])

    t = ts.astype(np.float64)
    ta, tb = float(t[a]), float(t[b])
    denom = tb - ta
    if abs(denom) < 1e-12:
        denom = 1e-12

    # tilde1 is linear in free u variables (4 vars)
    # free: t1[i_k] = u_k
    # dependent:
    #   b1 = Σ_k ((ta - t_free[k]) / denom) u_k
    #   a1 = -Σ_k u_k - b1
    t_free = t[free_idx]
    B = (ta - t_free) / denom                     # coefficients for b1
    A = -np.ones_like(B) - B                      # coefficients for a1

    # Build matrix M1 of shape (NUM_POINTS,4): t1[i] = sum_k M1[i,k] u_k
    M1 = np.zeros((NUM_POINTS, 4), dtype=np.float64)
    for k, idx in enumerate(free_idx):
        M1[idx, k] = 1.0
    M1[a, :] = A
    M1[b, :] = B

    # interval combos and sign vectors for each combo
    mids = np.array([interval_mid(ts, k) for k in range(NUM_POINTS + 1)], dtype=np.float64)
    combos = []
    for ia in range(NUM_POINTS + 1):
        tx = float(mids[ia])
        for ib in range(ia + 1, NUM_POINTS + 1):
            ty = float(mids[ib])
            if tx >= ty - 1e-6:
                continue
            Ci = (tx - t) * (ty - t)
            # sign(C_i C_j) per edge:
            s = np.sign(Ci[EDGE_I] * Ci[EDGE_J]).astype(np.int8)
            # avoid zeros
            s[s == 0] = 1
            combos.append((ia, ib, tx, ty, s))

    return {
        "ts": ts,
        "pair": pair,
        "free_idx": free_idx,
        "a": a,
        "b": b,
        "t": t,
        "M1": M1,
        "combos": combos,
        "ta": ta,
        "tb": tb,
        "denom": denom,
        "t_free": t_free
    }


def compute_t0_for_signs(base, s_free):
    """
    Compute t0 vector given signs for free indices.
    s_free: array/list of 4 signs (+1 or -1)
    """
    free_idx = base["free_idx"]
    a, b = base["a"], base["b"]
    ta, tb = base["ta"], base["tb"]
    denom = base["denom"]
    t_free = base["t_free"]
    
    s_free_arr = np.array(s_free, dtype=np.float64)
    
    # tilde0 constants: free -> s_free, dependent chosen to solve Σ λ_i tilde0_i = 0
    S00 = float(np.sum(s_free_arr))
    S10 = float(np.sum(t_free * s_free_arr))

    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0

    t0 = np.zeros(NUM_POINTS, dtype=np.float64)
    t0[free_idx] = s_free_arr
    t0[a] = a0
    t0[b] = b0
    return t0


def build_C_matrix(base, t0):
    """
    Compute C matrix (15x4) for [j i] = t0[j]*t1[i] - t1[j]*t0[i]
    """
    M1 = base["M1"]
    C = np.zeros((len(EDGE_PAIRS), 4), dtype=np.float64)
    for e, (i, j) in enumerate(EDGE_PAIRS):
        # [j i] = t0[j] * (M1[i] . u) - t0[i] * (M1[j] . u)
        #       = (t0[j]*M1[i] - t0[i]*M1[j]) . u
        C[e, :] = t0[j] * M1[i, :] - t0[i] * M1[j, :]
    return C


def solve_lp_max_margin(C, s_edge):
    """
    Solve max-margin LP for u in R^4 and m>=0
    """
    # x = [u0,u1,u2,u3,m] length 5
    A_ub = np.zeros((C.shape[0], 5), dtype=np.float64)
    b_ub = np.zeros((C.shape[0],), dtype=np.float64)

    A_ub[:, :4] = -(s_edge[:, None] * C)
    A_ub[:, 4] = 1.0

    c = np.array([0, 0, 0, 0, -1], dtype=np.float64) # Maximize m
    bounds = [(-U_CLIP, U_CLIP)] * 4 + [(0, None)]

    if SCIPY_OK:
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if not res.success:
            return False, None, None
        x = res.x
        u = x[:4]
        m = x[4]
        return (m is not None and m > 0), u, float(m)
    
    # Fallback to Sage LP (omitted for brevity unless SciPy fails, assumption is SciPy is present or we rely on v8 fallback if needed)
    return False, None, None


def compute_z_v9(base, t0, u, tx, ty):
    """
    Reconstruct z_ij
    """
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


def worker_chunk(task_id, restarts, seed):
    rng = np.random.default_rng(int(seed))
    best_overall = None # Track best pos if no 15 found

    # Sign patterns (16)
    sign_patterns = list(itertools.product([-1, 1], repeat=4))

    for _ in range(restarts):
        ts = sample_ts(rng)
        pair = pick_pair(rng)
        
        base = build_base_struct(ts, pair)
        
        # Iterate over all sign patterns for free t0
        for s_free in sign_patterns:
            t0 = compute_t0_for_signs(base, s_free)
            C_mat = build_C_matrix(base, t0)
            
            # Iterate over all interval combos
            for (ia, ib, tx_mid, ty_mid, s_edge) in base["combos"]:
                ok, u, m = solve_lp_max_margin(C_mat, s_edge)
                
                if not ok or m is None or m <= MARGIN_MIN:
                    continue
                
                # Verify solution
                z, _, _ = compute_z_v9(base, t0, u, tx_mid, ty_mid)
                pos = int(np.count_nonzero(z > 0))
                zmin = float(np.min(z))

                if zmin < ZMIN_MIN:
                    continue
                
                cand = {
                    "ts": ts.tolist(),
                    "pair": [int(pair[0]), int(pair[1])],
                    "free_idx": [int(i) for i in base["free_idx"]],
                    "s_free": list(s_free),
                    "u_free": u.tolist(),
                    "tx": float(tx_mid),
                    "ty": float(ty_mid),
                    "margin": float(m),
                    "zmin": float(zmin),
                    "pos": pos,
                    "ia": ia,
                    "ib": ib
                }

                if best_overall is None or pos > best_overall["pos"] or (pos == best_overall["pos"] and zmin > best_overall["zmin"]):
                    best_overall = cand

                if pos == 15:
                    # Consistency check: Rebuild exactly and verify
                    base_chk = build_base_struct(ts, pair)
                    t0_chk = compute_t0_for_signs(base_chk, s_free)
                    # Use exact same u, tx, ty
                    z_chk, _, _ = compute_z_v9(base_chk, t0_chk, u, tx_mid, ty_mid)
                    pos_chk = int(np.count_nonzero(z_chk > 0))
                    zmin_chk = float(np.min(z_chk))
                    
                    if pos_chk != 15 or zmin_chk < ZMIN_MIN:
                         # Failed consistency check (likely boundary noise)
                         continue

                    # Found a candidate! Refine (tx, ty) slightly
                    # For v9, we return immediately to check charts
                    return {"found": True, "candidate": cand, "best": best_overall}

    return {"found": False, "candidate": None, "best": best_overall}


def check_active_charts_and_save(candidate):
    """
    Heavy Sage checks: reconstruct, verify pos=15, check charts.
    """
    from sage.all import RealField, vector, load
    
    # We reuse solve_conservation_pair but need to handle t0 carefully
    load("src/atlas/jacobian_fiber_gauge.sage")
    load("src/atlas/solve_conservation_pair.sage") # This is for generic reconstruction

    RF = RealField(200)
    evaluator = FiberGaugeEvaluator(NUM_POINTS)

    ts = [float(x) for x in candidate["ts"]]
    pair = (int(candidate["pair"][0]), int(candidate["pair"][1]))
    free_idx = [int(i) for i in candidate["free_idx"]]
    s_free = candidate["s_free"]
    u_free_init = candidate["u_free"]

    lambdas_sage = {i: vector(RF, [1, ts[i]]) for i in range(NUM_POINTS)}
    tx = float(candidate["tx"])
    ty = float(candidate["ty"])
    x_sage = vector(RF, [1, tx])
    y_sage = vector(RF, [1, ty])
    
    def br(u, v):
        return u[0]*v[1] - u[1]*v[0]
    
    C_poly = {i: br(lambdas_sage[i], x_sage) * br(lambdas_sage[i], y_sage) for i in range(NUM_POINTS)}

    def check_point(u_curr):
        tildes_free_sage = {}
        for k, idx in enumerate(free_idx):
            tildes_free_sage[idx] = vector(RF, [float(s_free[k]), float(u_curr[k])])
            
        tildes_recon = solve_conservation_pair(lambdas_sage, tildes_free_sage, pair, NUM_POINTS, RF)
        if tildes_recon is None:
            return 0, -1, None
            
        pos_count = 0
        zmin_val = None
        z_out = []
        for (i, j) in EDGE_PAIRS:
            ang = br(lambdas_sage[i], lambdas_sage[j])
            sq  = br(tildes_recon[j], tildes_recon[i])
            val = (sq / ang) * C_poly[i] * C_poly[j]
            zf = float(val)
            z_out.append(zf)
            if zf > 0:
                pos_count += 1
            zmin_val = zf if zmin_val is None else min(zmin_val, zf)
        return pos_count, zmin_val, z_out

    # Initial check
    pos_count, zmin, z_vals = check_point(u_free_init)
    
    print(f"  Reconstructed Positivity (hi-prec): {pos_count}/15  (zmin={zmin:.3e})")
    
    u_final = list(u_free_init)
    
    if pos_count != 15:
        print("  [reject] Point failed high-prec check.")
        print(f"   Debug: pair={pair}, free_idx={free_idx}, s_free={s_free}")
        # Identify flipped edges
        # We can't easily see which edges flipped without re-running compute_z_v9 equivalent or storing expectation
        # But z_vals has current values.
        neg_edges = [k for k, val in enumerate(z_vals) if val <= 0]
        print(f"   Negative edges indices: {neg_edges}")
        
        print("  [rescue] Attempting local rescue...")
        import random
        rng = random.Random(42)
        
        # Rescue loop
        best_rescue_pos = pos_count
        best_rescue_zmin = zmin
        
        for attempt in range(200):
            # Perturb u slightly
            u_perturb = [x + rng.gauss(0, 1e-6) for x in u_free_init]
            
            p, zm, zv = check_point(u_perturb)
            if p > best_rescue_pos or (p == best_rescue_pos and zm > best_rescue_zmin):
                best_rescue_pos = p
                best_rescue_zmin = zm
            
            if p == 15:
                print(f"  [rescue] Success! Found valid point at attempt {attempt} (zmin={zm:.3e})")
                u_final = u_perturb
                pos_count = 15
                zmin = zm
                z_vals = zv
                break
        
        if pos_count != 15:
             print(f"  [rescue] Failed. Best was {best_rescue_pos}/15 (zmin={best_rescue_zmin:.3e})")
             return False

    # Chart check
    tf_final = {int(i): [float(s_free[k]), float(u_final[k])] for k, i in enumerate(free_idx)}
    xs_final = [1.0, float(tx)]
    ys_final = [1.0, float(ty)]

    # Quick check
    quick_roots = ALL_ROOTS[:10]
    any_hit = False
    for roots in quick_roots:
        val = evaluator.evaluate_chart_fiber_gauge(roots, ts, tf_final, xs_final, ys_final, require_polytope=True)
        if abs(val) > 1e-20:
            any_hit = True
            break
    
    if not any_hit:
        print("  [reject] Quick check: 0 hits in first 10 charts.")
        # Try full check anyway if rescue succeeded? No, expensive.
        # But if we rescued, we really want to know.
        # Continue to full check.
    
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
    out["u_free"] = u_final # Save the rescued u
    out["z_values"] = z_vals
    out["active_charts"] = int(active_charts)
    out["active_roots"] = [list(r) for r in hit_roots]
    out["meta"] = {
        "TOTAL_RESTARTS": int(TOTAL_RESTARTS),
        "CHUNK_RESTARTS": int(CHUNK_RESTARTS),
        "U_CLIP": float(U_CLIP),
        "MARGIN_MIN": float(MARGIN_MIN),
        "timestamp": float(time.time()),
        "rescued": (u_final != u_free_init)
    }

    with open("RESULTS/golden_kinematics.json", "w") as f:
        json.dump(out, f, indent=2)

    print(f"  ✨ GOLDEN FIND! Active Charts: {active_charts}")
    print("  Saved: RESULTS/golden_kinematics.json")
    return True


def main():
    print("Starting v9 Signed-T0 LP-feasibility search.")
    print(f"  TOTAL_RESTARTS={TOTAL_RESTARTS}, CHUNK_RESTARTS={CHUNK_RESTARTS}")
    print(f"  U_CLIP={U_CLIP}, MARGIN_MIN={MARGIN_MIN}, ZMIN_MIN={ZMIN_MIN}")
    print(f"  SciPy available: {SCIPY_OK}")

    cores = os.cpu_count() or 1
    if MAX_WORKERS is None:
        workers = min(max(1, cores - 2), 12)
    else:
        workers = min(MAX_WORKERS, cores)

    num_chunks = int(math.ceil(TOTAL_RESTARTS / float(CHUNK_RESTARTS)))
    print(f"  cores={cores} -> workers={workers}, chunks={num_chunks}")

    started = time.time()
    best_global = None

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = []
        for k in range(num_chunks):
            restarts = CHUNK_RESTARTS if (k < num_chunks - 1) else (TOTAL_RESTARTS - CHUNK_RESTARTS * (num_chunks - 1))
            futures.append(ex.submit(worker_chunk, k, restarts, SEED + 1000 * k))

        done = 0
        for fut in as_completed(futures):
            done += 1
            res = fut.result()
            
            b = res.get("best", None)
            if b is not None:
                if best_global is None or b["pos"] > best_global["pos"] or (b["pos"] == best_global["pos"] and b["zmin"] > best_global["zmin"]):
                    best_global = b
                    print(f"[progress] Best POS={best_global['pos']}/15, zmin={best_global['zmin']:.3e} (t={time.time()-started:.1f}s)")

            if res["found"]:
                cand = res["candidate"]
                print(f"[candidate] Found POS={cand['pos']} zmin={cand['zmin']:.3e} margin={cand['margin']:.3e}. Checking charts...")
                
                ok = check_active_charts_and_save(cand)
                if ok:
                    print("Golden point confirmed. Stopping search.")
                    for f in futures:
                        f.cancel()
                    ex.shutdown(wait=False)
                    return
                else:
                    print("[continue] Candidate rejected.")

            if done % 1 == 0:
                print(f"[status] chunks={done}/{num_chunks}, elapsed={time.time()-started:.1f}s")

    print("Search exhausted.")
    if best_global:
        print(f"Best result: POS={best_global['pos']}/15, zmin={best_global['zmin']:.3e}")

if __name__ == "__main__":
    main()
