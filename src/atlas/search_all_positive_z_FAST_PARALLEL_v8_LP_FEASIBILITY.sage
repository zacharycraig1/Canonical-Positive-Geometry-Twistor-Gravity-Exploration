# src/atlas/search_all_positive_z_FAST_PARALLEL_v8_LP_FEASIBILITY.sage
#
# v8: Major algorithmic optimization: deterministic LP feasibility for 15/15 positivity
#
# Core idea:
#   For fixed (ts, dependence pair (a,b), interval-combo for (tx,ty)), the signs of C_i C_j are fixed.
#   Under the fiber-gauge parametrization used throughout this project:
#     - free tilde vectors are (tilde0=1, tilde1=u_k) for 4 free indices
#     - dependent tilde1 entries (a1,b1) are linear in u due to conservation
#   Therefore each Plücker bracket [j i] is a linear form in u, and:
#       z_{ij} > 0  <=>  sign(C_i C_j) * [j i](u) > 0
#   which is a system of strict linear inequalities in u (R^4).
#
# We solve a max-margin LP:
#   maximize m  subject to  s_e * (c_e · u) >= m  for all 15 edges e
#                          -U_CLIP <= u_k <= U_CLIP
#                           0 <= m
#
# If m > 0, we have a robust 15/15 positive point essentially instantly.
# Then we run the usual high-precision reconstruction + active chart/polytope gate.
#
# This can be orders-of-magnitude faster than stochastic SA when 15/15 is rare.
#
# Requires SciPy (linprog). If unavailable, the script falls back to a Sage LP solver.
#
# Output:
#   RESULTS/golden_kinematics.json
#
# Run:
#   sage src/atlas/search_all_positive_z_FAST_PARALLEL_v8_LP_FEASIBILITY.sage

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

TOTAL_RESTARTS      = 1200
CHUNK_RESTARTS      = 30
SEED                = 3
MAX_WORKERS         = None   # auto: min(cores-2, 12)

U_CLIP              = 20.0
MARGIN_MIN          = 1e-8    # accept LP solutions with margin >= this

# (tx,ty) intervals
REFINE_XY_TRIES     = 14      # refine inside chosen interval to increase zmin / robustness

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
      - t0_const (constants tilde0_i)
      - bracket linear coefficients c_e in R^4 for each edge e (for bracket [j i])
      - interval combos for (tx,ty): ia<ib in {0..N}
        and signs s_e = sign(C_i C_j) for each edge.
    """
    free_idx = [i for i in range(N) if i not in pair]
    a, b = int(pair[0]), int(pair[1])

    t = ts.astype(np.float64)
    ta, tb = float(t[a]), float(t[b])
    denom = tb - ta
    if abs(denom) < 1e-12:
        denom = 1e-12

    # tilde0 constants: free -> 1, dependent chosen to solve Σ λ_i tilde0_i = 0
    S00 = float(len(free_idx))
    S10 = float(np.sum(t[free_idx]))

    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0

    t0 = np.zeros(N, dtype=np.float64)
    t0[free_idx] = 1.0
    t0[a] = a0
    t0[b] = b0

    # tilde1 is linear in free u variables (4 vars)
    # free: t1[i_k] = u_k
    # dependent:
    #   b1 = Σ_k ((ta - t_free[k]) / denom) u_k
    #   a1 = -Σ_k u_k - b1
    t_free = t[free_idx]
    B = (ta - t_free) / denom                     # coefficients for b1
    A = -np.ones_like(B) - B                      # coefficients for a1

    # Build matrix M1 of shape (N,4): t1[i] = sum_k M1[i,k] u_k
    M1 = np.zeros((N, 4), dtype=np.float64)
    for k, idx in enumerate(free_idx):
        M1[idx, k] = 1.0
    M1[a, :] = A
    M1[b, :] = B

    # For each edge (i<j), bracket [j i] = t0[j]*t1[i] - t1[j]*t0[i] = dot(c_e, u)
    C = np.zeros((len(EDGE_PAIRS), 4), dtype=np.float64)
    for e, (i, j) in enumerate(EDGE_PAIRS):
        C[e, :] = t0[j] * M1[i, :] - t0[i] * M1[j, :]

    # interval combos and sign vectors for each combo
    mids = np.array([interval_mid(ts, k) for k in range(N + 1)], dtype=np.float64)
    combos = []
    for ia in range(N + 1):
        tx = float(mids[ia])
        for ib in range(ia + 1, N + 1):
            ty = float(mids[ib])
            if tx >= ty - 1e-6:
                continue
            Ci = (tx - t) * (ty - t)
            # sign(C_i C_j) per edge:
            s = np.sign(Ci[EDGE_I] * Ci[EDGE_J]).astype(np.int8)
            # avoid zeros: if any Ci hits 0 (extremely unlikely with midpoints), nudge
            s[s == 0] = 1
            combos.append((ia, ib, tx, ty, s))

    return {
        "ts": ts,
        "pair": pair,
        "free_idx": free_idx,
        "a": a,
        "b": b,
        "t": t,
        "t0": t0,
        "C": C,             # (15,4)
        "combos": combos,   # list of (ia,ib,tx,ty,s_edge[15])
    }


def solve_lp_max_margin(C, s_edge):
    """
    Solve max-margin LP for u in R^4 and m>=0:
      s_e * (C_e · u) >= m for all e
      -U_CLIP <= u_k <= U_CLIP
      m >= 0
    Returns (ok, u, m)
    """
    # x = [u0,u1,u2,u3,m] length 5
    # constraint: s*(C·u) - m >= 0  <=>  -s*C·u + m <= 0
    A_ub = np.zeros((C.shape[0], 5), dtype=np.float64)
    b_ub = np.zeros((C.shape[0],), dtype=np.float64)

    A_ub[:, :4] = -(s_edge[:, None] * C)
    A_ub[:, 4] = 1.0

    # objective: minimize -m => c = [0,0,0,0,-1]
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

    # Fallback to Sage LP solver
    try:
        from sage.all import MixedIntegerLinearProgram, QQ
        p = MixedIntegerLinearProgram(maximization=True, solver="GLPK")
        uvar = p.new_variable(real=True)
        mvar = p.new_variable(real=True)

        # maximize m
        p.set_objective(mvar[0])

        # bounds
        for k in range(4):
            p.add_constraint(uvar[k] >= -U_CLIP)
            p.add_constraint(uvar[k] <=  U_CLIP)
        p.add_constraint(mvar[0] >= 0)

        for e in range(C.shape[0]):
            expr = 0
            for k in range(4):
                expr += float(s_edge[e]) * float(C[e, k]) * uvar[k]
            p.add_constraint(expr - mvar[0] >= 0)

        p.solve()
        u = np.array([float(p.get_values(uvar[k])) for k in range(4)], dtype=np.float64)
        m = float(p.get_values(mvar[0]))
        return (m > 0), u, m
    except Exception:
        return False, None, None


def compute_z(ts, pair, free_idx, u, tx, ty):
    """
    Reconstruct the full tilde1 vector (including dependent) in the same linear form and compute z_ij.
    """
    t = np.array(ts, dtype=np.float64)
    a, b = int(pair[0]), int(pair[1])

    ta, tb = float(t[a]), float(t[b])
    denom = tb - ta
    if abs(denom) < 1e-12:
        denom = 1e-12

    t_free = t[free_idx]
    S01 = float(np.sum(u))
    S11 = float(np.sum(t_free * u))

    b1 = (-S11 + ta * S01) / denom
    a1 = -S01 - b1

    # tilde0 constants
    S00 = float(len(free_idx))
    S10 = float(np.sum(t[free_idx]))
    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0

    t0 = np.zeros(N, dtype=np.float64)
    t0[free_idx] = 1.0
    t0[a] = a0
    t0[b] = b0

    t1 = np.zeros(N, dtype=np.float64)
    for k, idx in enumerate(free_idx):
        t1[idx] = float(u[k])
    t1[a] = float(a1)
    t1[b] = float(b1)

    ang = (t[EDGE_J] - t[EDGE_I])
    ang = np.where(np.abs(ang) < 1e-12, 1e-12, ang)
    ang_inv = 1.0 / ang

    C = (tx - t) * (ty - t)
    Cprod = (C[EDGE_I] * C[EDGE_J])

    br = t0[EDGE_J] * t1[EDGE_I] - t1[EDGE_J] * t0[EDGE_I]
    z = (br * ang_inv) * Cprod
    return z, t0, t1


def find_positive_point_for_ts_pair(struct, rng):
    """
    For this (ts,pair), try all interval combos; solve LP; keep best by margin m.
    Also refine tx,ty within the winning interval pair to increase robustness.
    """
    best = None  # (m, pos, zmin, u, ia, ib, tx, ty, z)

    C = struct["C"]
    ts = struct["ts"]
    pair = struct["pair"]
    free_idx = struct["free_idx"]

    for (ia, ib, tx_mid, ty_mid, s_edge) in struct["combos"]:
        ok, u, m = solve_lp_max_margin(C, s_edge)
        if not ok or m is None or m < MARGIN_MIN:
            continue

        z, _, _ = compute_z(ts, pair, free_idx, u, tx_mid, ty_mid)
        pos = int(np.count_nonzero(z > 0))
        zmin = float(np.min(z))
        if pos != 15:
            # numerical mismatch (should be rare); skip
            continue

        cand = (float(m), pos, zmin, u.copy(), int(ia), int(ib), float(tx_mid), float(ty_mid), z.copy())
        if best is None or cand[0] > best[0]:
            best = cand

    if best is None:
        return None

    m, pos, zmin, u, ia, ib, tx0, ty0, z0 = best

    # refine tx,ty inside the same intervals to improve zmin (often helps downstream checks)
    best2 = (m, pos, zmin, u, ia, ib, tx0, ty0, z0)
    for _ in range(REFINE_XY_TRIES):
        tx, ty = sample_tx_ty_in_intervals(ts, ia, ib, rng)
        z, _, _ = compute_z(ts, pair, free_idx, u, tx, ty)
        pos2 = int(np.count_nonzero(z > 0))
        if pos2 != 15:
            continue
        zmin2 = float(np.min(z))
        # prioritize robustness zmin, then margin
        if zmin2 > best2[2]:
            best2 = (m, pos2, zmin2, u, ia, ib, tx, ty, z)

    return best2


def worker_chunk(task_id, restarts, seed):
    rng = np.random.default_rng(int(seed))
    best_overall = None

    for _ in range(restarts):
        ts = sample_ts(rng)
        pair = pick_pair(rng)
        struct = build_struct(ts, pair)

        out = find_positive_point_for_ts_pair(struct, rng)
        if out is None:
            continue

        m, pos, zmin, u, ia, ib, tx, ty, z = out
        cand = {
            "ts": ts.tolist(),
            "pair": [int(pair[0]), int(pair[1])],
            "free_idx": [int(i) for i in struct["free_idx"]],
            "u_free": u.tolist(),
            "tx": float(tx),
            "ty": float(ty),
            "margin": float(m),
            "zmin": float(zmin),
        }

        if best_overall is None or cand["zmin"] > best_overall["zmin"]:
            best_overall = cand

        # If zmin is robust, return early as a strong candidate for chart check.
        if cand["zmin"] >= 1e-8:
            return {"found": True, "candidate": cand, "best": best_overall}

    return {"found": False, "candidate": None, "best": best_overall}


def check_active_charts_and_save(candidate):
    """
    Heavy Sage checks only happen in main process:
      - reconstruct exact conservation via solve_conservation_pair (RealField 200)
      - verify 15/15 positivity
      - check active charts with require_polytope=True (early exit then full count)
    """
    from sage.all import RealField, vector, load

    load("src/atlas/jacobian_fiber_gauge.sage")
    load("src/atlas/solve_conservation_pair.sage")

    RF = RealField(200)
    evaluator = FiberGaugeEvaluator(N)

    ts = [float(x) for x in candidate["ts"]]
    pair = (int(candidate["pair"][0]), int(candidate["pair"][1]))
    free_idx = [i for i in range(N) if i not in pair]

    lambdas_sage = {i: vector(RF, [1, ts[i]]) for i in range(N)}
    u_free = candidate["u_free"]
    tildes_free_sage = {i: vector(RF, [1, float(u_free[k])]) for k, i in enumerate(free_idx)}

    tildes_recon = solve_conservation_pair(lambdas_sage, tildes_free_sage, pair, N, RF)
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

    # Chart check: early exit, then full enumeration for save
    tf_final = {int(i): [1.0, float(u_free[k])] for k, i in enumerate(free_idx)}
    xs_final = [1.0, float(tx)]
    ys_final = [1.0, float(ty)]

    # quick early exit
    quick_roots = ALL_ROOTS[:6]
    any_hit = False
    for roots in quick_roots:
        val = evaluator.evaluate_chart_fiber_gauge(roots, ts, tf_final, xs_final, ys_final, require_polytope=True)
        if abs(val) > 1e-20:
            any_hit = True
            break
    if not any_hit:
        print("  [reject] quick roots had 0 hits (likely no active charts).")
        return False

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
    out["z_values"] = z_vals
    out["active_charts"] = int(active_charts)
    out["active_roots"] = [list(r) for r in hit_roots]
    out["meta"] = {
        "TOTAL_RESTARTS": int(TOTAL_RESTARTS),
        "CHUNK_RESTARTS": int(CHUNK_RESTARTS),
        "U_CLIP": float(U_CLIP),
        "MARGIN_MIN": float(MARGIN_MIN),
        "SCIPY_OK": bool(SCIPY_OK),
        "timestamp": float(time.time()),
    }

    with open("RESULTS/golden_kinematics.json", "w") as f:
        json.dump(out, f, indent=2)

    print(f"  ✨ GOLDEN FIND! Active Charts: {active_charts}")
    print("  Saved: RESULTS/golden_kinematics.json")
    return True


def main():
    print("Starting v8 LP-feasibility search for 15/15 all-positive z_ij (then chart/polytope gate).")
    print(f"  TOTAL_RESTARTS={TOTAL_RESTARTS}, CHUNK_RESTARTS={CHUNK_RESTARTS}, U_CLIP={U_CLIP}")
    print(f"  SciPy linprog available: {SCIPY_OK}")

    import multiprocessing as mp
    cores = os.cpu_count() or 1
    if MAX_WORKERS is None:
        workers = min(max(1, cores - 2), 12)
    else:
        workers = min(MAX_WORKERS, cores)

    num_chunks = int(math.ceil(TOTAL_RESTARTS / float(CHUNK_RESTARTS)))

    best_global = None
    started = time.time()
    print(f"  cores={cores} -> workers={workers}, chunks={num_chunks}")

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
            if b is not None and (best_global is None or b["zmin"] > best_global["zmin"]):
                best_global = b
                print(f"[progress] best zmin={best_global['zmin']:.3e} margin={best_global.get('margin',0):.3e} t={time.time()-started:.1f}s")

            if res["found"]:
                cand = res["candidate"]
                print(f"[candidate] LP found 15/15 (zmin={cand['zmin']:.3e}, margin={cand['margin']:.3e}). Checking charts...")

                ok = check_active_charts_and_save(cand)
                if ok:
                    for f in futures:
                        f.cancel()
                    try:
                        ex.shutdown(wait=False, cancel_futures=True)
                    except TypeError:
                        ex.shutdown(wait=False)
                    return
                else:
                    print("[continue] Candidate rejected. Continuing search...")

            if done % 1 == 0:
                print(f"[status] chunks_done={done}/{num_chunks}, elapsed={time.time()-started:.1f}s")

    print("Search exhausted without golden point.")
    if best_global is not None:
        print(f"Best candidate zmin={best_global['zmin']:.3e} margin={best_global.get('margin',0):.3e}")


if __name__ == "__main__":
    main()
