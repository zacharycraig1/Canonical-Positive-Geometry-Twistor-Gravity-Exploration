# src/atlas/search_all_positive_z_FAST_PARALLEL_v9_LP_SCALES_TOPK.sage
#
# v9: Major "should work" upgrades after v8 exhausted quickly with no golden:
#
# 1) Remove the biggest over-restriction of v8:
#    v8 fixed each free tilde vector to (1, u_k). That is only 1 DOF per free leg.
#    The true free space is 2 DOF per free leg (8 DOF total).
#
#    v9 keeps the LP *structure* by solving for slopes u via LP feasibility, BUT then
#    introduces per-leg positive scales s_k and passes full free vectors:
#         tilde_i = (s_i, s_i*u_i)
#    (s_i > 0)
#    This often fixes "active charts always 0" if the (1,u) restriction was killing charts.
#
# 2) Top-K candidate queue:
#    v8 only chart-checked "robust zmin" candidates. If all feasible points are near-boundary,
#    it may have checked none. v9 always returns the TOP_K best 15/15 candidates (by zmin, margin)
#    and the main process chart-checks them.
#
# 3) Fast dependent solve (float) in workers (NO Sage in workers):
#    Conservation Σ λ_i ⊗ tilde_i = 0 with λ_i=(1,t_i) gives 4 scalar equations:
#      Σ u0_i = 0, Σ u1_i = 0, Σ t_i u0_i = 0, Σ t_i u1_i = 0
#    For chosen dependent pair (a,b), solve two 2x2 systems to obtain tilde_a and tilde_b.
#    Sage high-precision reconstruction happens only at final verification.
#
# Output: RESULTS/golden_kinematics.json
#
# Run:
#   sage src/atlas/search_all_positive_z_FAST_PARALLEL_v9_LP_SCALES_TOPK.sage

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
N = 6

TOTAL_RESTARTS      = 8000         # cheap now; adjust up freely
CHUNK_RESTARTS      = 80
SEED                = 3
MAX_WORKERS         = None         # auto: min(cores-2, 12)

U_CLIP              = 30.0
MARGIN_MIN          = 0.0          # accept any m>0
TOP_K_PER_WORKER    = 6

# Scale sampling around LP slopes u:
SCALE_SAMPLES       = 18
S_LOGSIGMA          = 0.9          # lognormal spread; higher explores more
S_MIN               = 1e-3
S_MAX               = 1e3

# (tx,ty) intervals
REFINE_XY_TRIES     = 18           # refine inside winning interval pair to increase zmin

# Pair bias (empirical)
PAIR_BIAS = {
    (0, 1): 4.0,
    (4, 5): 4.0,
    (0, 5): 2.0,
    (1, 4): 2.0,
}

EDGE_PAIRS = [(i, j) for i in range(N) for j in range(i + 1, N)]
EDGE_I = np.array([i for (i, j) in EDGE_PAIRS], dtype=np.int64)
EDGE_J = np.array([j for (i, j) in EDGE_PAIRS], dtype=np.int64)
ALL_ROOTS = list(itertools.combinations(range(N), 3))

# SciPy LP
try:
    from scipy.optimize import linprog
    SCIPY_OK = True
except Exception:
    linprog = None
    SCIPY_OK = False


def _pair_weights():
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    w = np.array([PAIR_BIAS.get(p, 1.0) for p in pairs], dtype=np.float64)
    w /= np.sum(w)
    return pairs, w


PAIRS, PAIR_W = _pair_weights()


def pick_pair(rng):
    idx = rng.choice(len(PAIRS), p=PAIR_W)
    return PAIRS[idx]


def sample_ts(rng):
    ts = np.sort(rng.uniform(1.0, 10.0, size=N).astype(np.float64))
    for k in range(1, N):
        if ts[k] <= ts[k - 1]:
            ts[k] = ts[k - 1] + 1e-6
    return ts


def interval_mid(ts, idx):
    if idx == 0:
        return float(ts[0] - 2.0)
    if idx == N:
        return float(ts[N - 1] + 2.0)
    return float(0.5 * (ts[idx - 1] + ts[idx]))


def interval_bounds(ts, idx):
    if idx == 0:
        return (-np.inf, float(ts[0]))
    if idx == N:
        return (float(ts[N - 1]), np.inf)
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


def build_struct(ts, pair):
    """
    Build:
      - linear coefficients C_e in R^4 for bracket [j i] using slope vars u_k (free indices only)
      - interval combos and s_edge = sign(C_i C_j) per edge
    """
    free_idx = [i for i in range(N) if i not in pair]
    a, b = int(pair[0]), int(pair[1])

    t = ts.astype(np.float64)
    ta, tb = float(t[a]), float(t[b])
    denom = tb - ta
    if abs(denom) < 1e-12:
        denom = 1e-12

    # Choose tilde0 constants (same as prior scripts):
    S00 = float(len(free_idx))
    S10 = float(np.sum(t[free_idx]))

    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0

    t0 = np.zeros(N, dtype=np.float64)
    t0[free_idx] = 1.0
    t0[a] = a0
    t0[b] = b0

    # Slope model for tilde1: free tilde1[i_k]=u_k
    t_free = t[free_idx]
    B = (ta - t_free) / denom
    A = -np.ones_like(B) - B

    M1 = np.zeros((N, 4), dtype=np.float64)
    for k, idx in enumerate(free_idx):
        M1[idx, k] = 1.0
    M1[a, :] = A
    M1[b, :] = B

    # bracket [j i] = dot(C_e, u)
    Ccoef = np.zeros((len(EDGE_PAIRS), 4), dtype=np.float64)
    for e, (i, j) in enumerate(EDGE_PAIRS):
        Ccoef[e, :] = t0[j] * M1[i, :] - t0[i] * M1[j, :]

    mids = np.array([interval_mid(ts, k) for k in range(N + 1)], dtype=np.float64)
    combos = []
    for ia in range(N + 1):
        tx = float(mids[ia])
        for ib in range(ia + 1, N + 1):
            ty = float(mids[ib])
            if tx >= ty - 1e-6:
                continue
            Ci = (tx - t) * (ty - t)
            s_edge = np.sign(Ci[EDGE_I] * Ci[EDGE_J]).astype(np.int8)
            s_edge[s_edge == 0] = 1
            combos.append((ia, ib, tx, ty, s_edge))

    return {
        "ts": ts,
        "pair": pair,
        "free_idx": free_idx,
        "a": a,
        "b": b,
        "t": t,
        "t0": t0,
        "Ccoef": Ccoef,
        "combos": combos,
    }


def solve_lp_max_margin(Ccoef, s_edge):
    """
    Maximize m s.t. s_e*(C_e·u) >= m, |u_k|<=U_CLIP, m>=0.
    Returns (ok, u, m)
    """
    # x = [u0,u1,u2,u3,m]
    A_ub = np.zeros((Ccoef.shape[0], 5), dtype=np.float64)
    b_ub = np.zeros((Ccoef.shape[0],), dtype=np.float64)
    A_ub[:, :4] = -(s_edge[:, None] * Ccoef)
    A_ub[:, 4] = 1.0
    c = np.array([0, 0, 0, 0, -1], dtype=np.float64)
    bounds = [(-U_CLIP, U_CLIP)] * 4 + [(0, None)]

    if SCIPY_OK:
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if not res.success:
            return False, None, None
        x = res.x
        u = x[:4]
        m = x[4]
        return (m is not None and m > 0), u, float(m)

    # Fallback: no solver => skip
    return False, None, None


def solve_dependent_pair_float(ts, pair, tildes_free):
    """
    Given free tilde vectors for 4 indices, solve for dependent pair (a,b) tilde vectors
    using the 4 scalar equations:
        Σ u0_i = 0
        Σ u1_i = 0
        Σ t_i u0_i = 0
        Σ t_i u1_i = 0
    Returns full dict tildes[i] = np.array([u0,u1])
    """
    a, b = int(pair[0]), int(pair[1])
    t = np.array(ts, dtype=np.float64)
    ta, tb = float(t[a]), float(t[b])
    det = ta - tb
    if abs(det) < 1e-12:
        return None

    # sums over known indices
    S0 = 0.0
    S1 = 0.0
    T0 = 0.0
    T1 = 0.0
    for i, v in tildes_free.items():
        u0, u1 = float(v[0]), float(v[1])
        S0 += u0
        S1 += u1
        T0 += float(t[i]) * u0
        T1 += float(t[i]) * u1

    # equations:
    # u0a + u0b = -S0
    # ta*u0a + tb*u0b = -T0
    # Solve:
    # u0a = (tb*(-S0) - (-T0)) / (tb - ta) = (T0 - tb*S0)/(tb-ta)
    # u0b = -S0 - u0a
    denom = (tb - ta)
    if abs(denom) < 1e-12:
        denom = 1e-12

    u0a = (T0 - tb * S0) / denom
    u0b = -S0 - u0a

    u1a = (T1 - tb * S1) / denom
    u1b = -S1 - u1a

    out = {int(i): np.array(v, dtype=np.float64) for i, v in tildes_free.items()}
    out[a] = np.array([u0a, u1a], dtype=np.float64)
    out[b] = np.array([u0b, u1b], dtype=np.float64)
    return out


def compute_z_from_full(ts, tildes, tx, ty):
    t = np.array(ts, dtype=np.float64)
    # angles: <ij> = t_j - t_i (since λ=(1,t))
    ang = (t[EDGE_J] - t[EDGE_I])
    ang = np.where(np.abs(ang) < 1e-12, 1e-12, ang)
    ang_inv = 1.0 / ang

    # C_i = <i x><i y> = (tx - t_i)(ty - t_i)
    C = (tx - t) * (ty - t)
    Cprod = (C[EDGE_I] * C[EDGE_J])

    # [j i] = det(tilde_j, tilde_i)
    T = np.zeros((N, 2), dtype=np.float64)
    for i in range(N):
        T[i, :] = tildes[i]
    br = T[EDGE_J, 0] * T[EDGE_I, 1] - T[EDGE_J, 1] * T[EDGE_I, 0]

    z = (br * ang_inv) * Cprod
    return z


def sample_scales(rng, d):
    # lognormal positive scales, clipped
    s = np.exp(rng.normal(0.0, S_LOGSIGMA, size=(d,)))
    s = np.clip(s, S_MIN, S_MAX)
    return s


def candidate_search_one(struct, rng):
    """
    For one (ts,pair), try combos:
      - Solve LP for u (slopes)
      - For each feasible u, randomize scales s and compute full tildes, test 15/15
      - Keep best by zmin and margin.
    Returns best candidate dict or None.
    """
    ts = struct["ts"]
    pair = struct["pair"]
    free_idx = struct["free_idx"]
    t = struct["t"]
    combos = struct["combos"]
    Ccoef = struct["Ccoef"]

    best = None  # dict

    # Try combos in descending "promising" order: larger |C| tends to help; approximate via |prod(Ci)| median
    # We'll just iterate; it's cheap.
    for (ia, ib, tx_mid, ty_mid, s_edge) in combos:
        ok, u, m = solve_lp_max_margin(Ccoef, s_edge)
        if not ok:
            continue
        if m <= MARGIN_MIN:
            continue

        # Build base slopes for free indices
        u = np.asarray(u, dtype=np.float64)

        # refine tx,ty inside same intervals (optional)
        best_tx, best_ty = float(tx_mid), float(ty_mid)
        best_local = None  # (zmin, tx, ty, tildes_free)

        for rxy in range(max(1, REFINE_XY_TRIES)):
            if rxy == 0:
                tx, ty = float(tx_mid), float(ty_mid)
            else:
                tx, ty = sample_tx_ty_in_intervals(ts, ia, ib, rng)

            # Scale sampling
            for _ in range(SCALE_SAMPLES):
                s = sample_scales(rng, len(free_idx))
                tildes_free = {}
                for k, idx in enumerate(free_idx):
                    tildes_free[int(idx)] = np.array([s[k], s[k] * u[k]], dtype=np.float64)

                tildes_full = solve_dependent_pair_float(ts, pair, tildes_free)
                if tildes_full is None:
                    continue

                z = compute_z_from_full(ts, tildes_full, tx, ty)
                if np.any(~np.isfinite(z)):
                    continue
                pos = int(np.count_nonzero(z > 0))
                if pos != 15:
                    continue
                zmin = float(np.min(z))

                if best_local is None or zmin > best_local[0]:
                    best_local = (zmin, tx, ty, tildes_free, z)

        if best_local is None:
            continue

        zmin, txB, tyB, tildes_freeB, zB = best_local
        cand = {
            "ts": ts.tolist(),
            "pair": [int(pair[0]), int(pair[1])],
            "tildes_free": {str(k): [float(v[0]), float(v[1])] for k, v in tildes_freeB.items()},
            "tx": float(txB),
            "ty": float(tyB),
            "margin": float(m),
            "zmin": float(zmin),
        }

        if best is None or cand["zmin"] > best["zmin"]:
            best = cand

    return best


def worker_chunk(task_id, restarts, seed):
    rng = np.random.default_rng(int(seed))
    top = []  # list of candidates

    for _ in range(restarts):
        ts = sample_ts(rng)
        pair = pick_pair(rng)
        struct = build_struct(ts, pair)
        cand = candidate_search_one(struct, rng)
        if cand is None:
            continue

        top.append(cand)
        # keep top-K by zmin
        top.sort(key=lambda d: d["zmin"], reverse=True)
        if len(top) > TOP_K_PER_WORKER:
            top = top[:TOP_K_PER_WORKER]

    return {"top": top}


def check_active_charts_and_save(candidate):
    """
    Heavy Sage checks in main process:
      - reconstruct exact conservation via solve_conservation_pair (RealField 200)
      - verify 15/15 positivity
      - check active charts with require_polytope=True
      - additionally measure active charts w/o polytope gate for debugging
    """
    from sage.all import RealField, vector, load

    load("src/atlas/jacobian_fiber_gauge.sage")
    load("src/atlas/solve_conservation_pair.sage")

    RF = RealField(200)
    evaluator = FiberGaugeEvaluator(N)

    ts = [float(x) for x in candidate["ts"]]
    pair = (int(candidate["pair"][0]), int(candidate["pair"][1]))

    # Parse tildes_free
    tf = {int(k): vector(RF, [float(v[0]), float(v[1])]) for k, v in candidate["tildes_free"].items()}
    lambdas_sage = {i: vector(RF, [1, ts[i]]) for i in range(N)}

    tildes_recon = solve_conservation_pair(lambdas_sage, tf, pair, N, RF)
    if tildes_recon is None:
        print("  [reject] solve_conservation_pair returned None.")
        return False

    tx = float(candidate["tx"])
    ty = float(candidate["ty"])
    x_sage = vector(RF, [1, tx])
    y_sage = vector(RF, [1, ty])

    def br(u, v):
        return u[0]*v[1] - u[1]*v[0]

    C = {i: br(lambdas_sage[i], x_sage) * br(lambdas_sage[i], y_sage) for i in range(N)}

    z_vals = []
    pos_count = 0
    zmin = None
    for (i, j) in EDGE_PAIRS:
        ang = br(lambdas_sage[i], lambdas_sage[j])
        sq  = br(tildes_recon[j], tildes_recon[i])
        val = (sq / ang) * C[i] * C[j]
        zf = float(val)
        z_vals.append(zf)
        if zf > 0:
            pos_count += 1
        zmin = zf if zmin is None else min(zmin, zf)

    print(f"  Reconstructed Positivity (hi-prec): {pos_count}/15  (zmin={zmin:.3e})")
    if pos_count != 15:
        print("  [reject] Lost positivity after hi-prec reconstruction.")
        return False

    # Build evaluator args (expects free tildes only as floats)
    tf_final = {int(k): [float(v[0]), float(v[1])] for k, v in candidate["tildes_free"].items()}
    xs_final = [1.0, float(tx)]
    ys_final = [1.0, float(ty)]

    # Debug: charts without polytope requirement
    active_no_poly = 0
    for roots in ALL_ROOTS:
        val = evaluator.evaluate_chart_fiber_gauge(roots, ts, tf_final, xs_final, ys_final, require_polytope=False)
        if abs(val) > 1e-20:
            active_no_poly += 1

    # Actual gate: require_polytope=True
    active_charts = 0
    hit_roots = []
    for roots in ALL_ROOTS:
        val = evaluator.evaluate_chart_fiber_gauge(roots, ts, tf_final, xs_final, ys_final, require_polytope=True)
        if abs(val) > 1e-20:
            active_charts += 1
            hit_roots.append(roots)

    print(f"  Active charts (no polytope): {active_no_poly} / 20")
    print(f"  Active charts (require_polytope): {active_charts} / 20")

    if active_charts <= 0:
        print("  [reject] 15/15 but no active charts under polytope requirement.")
        return False

    os.makedirs("RESULTS", exist_ok=True)
    out = dict(candidate)
    out["z_values"] = z_vals
    out["active_charts_no_polytope"] = int(active_no_poly)
    out["active_charts"] = int(active_charts)
    out["active_roots"] = [list(r) for r in hit_roots]
    out["meta"] = {
        "TOTAL_RESTARTS": int(TOTAL_RESTARTS),
        "CHUNK_RESTARTS": int(CHUNK_RESTARTS),
        "U_CLIP": float(U_CLIP),
        "MARGIN_MIN": float(MARGIN_MIN),
        "SCALE_SAMPLES": int(SCALE_SAMPLES),
        "S_LOGSIGMA": float(S_LOGSIGMA),
        "TOP_K_PER_WORKER": int(TOP_K_PER_WORKER),
        "SCIPY_OK": bool(SCIPY_OK),
        "timestamp": float(time.time()),
    }

    with open("RESULTS/golden_kinematics.json", "w") as f:
        json.dump(out, f, indent=2)

    print(f"  ✨ GOLDEN FIND! Active Charts: {active_charts}")
    print("  Saved: RESULTS/golden_kinematics.json")
    return True


def main():
    print("Starting v9 LP-slopes + scale-lift + topK chart checking.")
    print(f"  TOTAL_RESTARTS={TOTAL_RESTARTS}, CHUNK_RESTARTS={CHUNK_RESTARTS}, workers auto")
    print(f"  SciPy linprog available: {SCIPY_OK}")

    import multiprocessing as mp
    cores = os.cpu_count() or 1
    if MAX_WORKERS is None:
        workers = min(max(1, cores - 2), 12)
    else:
        workers = min(MAX_WORKERS, cores)

    num_chunks = int(math.ceil(TOTAL_RESTARTS / float(CHUNK_RESTARTS)))

    started = time.time()
    print(f"  cores={cores} -> workers={workers}, chunks={num_chunks}")

    all_candidates = []

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = []
        for k in range(num_chunks):
            restarts = CHUNK_RESTARTS if (k < num_chunks - 1) else (TOTAL_RESTARTS - CHUNK_RESTARTS * (num_chunks - 1))
            futures.append(ex.submit(worker_chunk, k, restarts, SEED + 1000 * k))

        done = 0
        for fut in as_completed(futures):
            done += 1
            res = fut.result()
            top = res.get("top", [])
            all_candidates.extend(top)

            # periodically check best few candidates so far
            if done % 3 == 0 and len(all_candidates) > 0:
                all_candidates.sort(key=lambda d: d["zmin"], reverse=True)
                probe = all_candidates[:min(10, len(all_candidates))]
                print(f"[status] chunks_done={done}/{num_chunks}, candidates={len(all_candidates)}, best_zmin={probe[0]['zmin']:.3e}, elapsed={time.time()-started:.1f}s")

                # try chart check on the very best candidate so far (cheap early win)
                cand0 = probe[0]
                print("[check] Trying best candidate so far...")
                if check_active_charts_and_save(cand0):
                    for f in futures:
                        f.cancel()
                    try:
                        ex.shutdown(wait=False, cancel_futures=True)
                    except TypeError:
                        ex.shutdown(wait=False)
                    return

    print("Search complete. Final chart-check sweep on top candidates...")
    if len(all_candidates) == 0:
        print("No 15/15 candidates found at all.")
        return

    all_candidates.sort(key=lambda d: d["zmin"], reverse=True)
    for k, cand in enumerate(all_candidates[:50]):
        print(f"[final-check] k={k} zmin={cand['zmin']:.3e} margin={cand.get('margin',0):.3e}")
        if check_active_charts_and_save(cand):
            return

    print("No golden point found among top candidates.")
    best = all_candidates[0]
    print(f"Best zmin={best['zmin']:.3e}, margin={best.get('margin',0):.3e}")


if __name__ == "__main__":
    main()
