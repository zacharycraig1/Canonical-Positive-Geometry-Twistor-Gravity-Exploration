# src/atlas/search_all_positive_z_FAST_PARALLEL.sage
#
# FAST + PARALLEL Simulated Annealing search for n=6 all-positive z_ij (15/15),
# with EXACT momentum conservation solved in closed form at every step.
#
# Designed to be "drop & run" on CoCalc (Linux) or any Sage install.
#
# Key features:
# - No SciPy minimize (huge speedup). Conservation enforced exactly (no penalty).
# - Parallel restarts via ProcessPoolExecutor (uses available cores).
# - Vectorized z-evaluation for all 15 edges.
# - Lazy-loads heavy Sage evaluator only for rare 15/15 candidates.
# - Early-stop via a shared stop_event (other workers wind down quickly once golden point is found).
#
# Output:
#   RESULTS/golden_kinematics.json
#
# Run:
#   sage src/atlas/search_all_positive_z_FAST_PARALLEL.sage
#
# Notes:
# - This script assumes your project has:
#     src/atlas/jacobian_fiber_gauge.sage
#     src/atlas/solve_conservation_pair.sage
#   providing FiberGaugeEvaluator and solve_conservation_pair as in your earlier pipeline.

import os, json, math, time
import itertools
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

# =========================
# CONFIG
# =========================
N = 6

TOTAL_RESTARTS     = 5000        # INCREASED for better coverage
STEPS_PER_RESTART  = 3000        # SA steps per restart
CHUNK_RESTARTS     = 50          # Batches for workers
SEED               = int(3)

# Parallelism
MAX_WORKERS        = 16          # Use all 16 cores

# SA temperature schedule
T_END              = 1e-4

# Softplus loss params
MARGIN             = 1e-6
BETA               = 0.5
SCALE_EPS          = 1e-12

# Proposal step sizes (auto-adapted mildly)
SIGMA_U_INIT       = 0.8         # for free tilde scalar params u_i (tildes free are [1,u])
SIGMA_XY_INIT      = 0.5         # for tx, ty

# Strict positivity threshold (use 0.0 for true >0; raise to demand margin)
POS_EPS            = 0.0

# Optional: require a safety margin on the minimal z_ij before invoking the heavy evaluator
ZMIN_TRIGGER       = 0.0         # e.g. 1e-10 if you want more robust positivity

# Logging
PRINT_EVERY_CHUNK  = 1

# Weighted dependent-pair restarts (bias toward empirically “good” pairs)
PAIR_BIAS = {
    (0, 1): 4.0,
    (4, 5): 4.0,
    (0, 5): 2.0,
    (1, 4): 2.0,
}

# =========================
# PRECOMPUTE INDICES
# =========================
EDGE_PAIRS = [(i, j) for i in range(N) for j in range(i + 1, N)]
EDGE_I = np.array([i for (i, j) in EDGE_PAIRS], dtype=np.int64)
EDGE_J = np.array([j for (i, j) in EDGE_PAIRS], dtype=np.int64)

ALL_ROOTS = list(itertools.combinations(range(N), 3))  # 20 roots


def _softplus_vec(x, scale):
    # softplus(x; scale) = scale*log(1+exp(x/scale)) stable form
    xs = x / scale
    return scale * (np.log1p(np.exp(-np.abs(xs))) + np.maximum(xs, 0.0))


def _pair_weights():
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    w = np.array([PAIR_BIAS.get(p, 1.0) for p in pairs], dtype=np.float64)
    w = w / np.sum(w)
    return pairs, w


PAIRS, PAIR_W = _pair_weights()


def pick_pair(rng):
    idx = rng.choice(len(PAIRS), p=PAIR_W)
    return PAIRS[idx]


def sample_ts(rng):
    # Moment-curve ordering (strictly increasing)
    ts = np.sort(rng.uniform(1.0, 10.0, size=N).astype(np.float64))
    for k in range(1, N):
        if ts[k] <= ts[k - 1]:
            ts[k] = ts[k - 1] + 1e-6
    return ts


def sample_in_interval(rng, lo, hi, base_scale=2.0):
    # Finite sampling for +/- inf endpoints
    if not np.isfinite(lo) and np.isfinite(hi):
        return float(hi - rng.uniform(0.5, base_scale))
    if np.isfinite(lo) and not np.isfinite(hi):
        return float(lo + rng.uniform(0.5, base_scale))
    if hi - lo < 1e-9:
        return float(lo)
    return float(rng.uniform(lo, hi))


def reconstruct_full_tildes(ts, pair, free_idx, u_free):
    """
    Gauge-fixed free tildes: tilde_i = [1, u_i] for i in free_idx.
    Solve conservation EXACTLY for pair indices (a,b).

    Conservation in spinor form:
      sum_i lambda_i ⊗ tilde_i = 0
    with lambda_i = [1, t_i].
    """
    T = np.empty((N, 2), dtype=np.float64)

    # Fill free indices
    for k, i in enumerate(free_idx):
        T[i, 0] = 1.0
        T[i, 1] = float(u_free[k])

    a, b = pair
    ta, tb = float(ts[a]), float(ts[b])
    denom = (tb - ta)
    if abs(denom) < 1e-12:
        denom = 1e-12

    # Compute partial sums from free indices only
    free0 = np.array([T[i, 0] for i in free_idx], dtype=np.float64)
    free1 = np.array([T[i, 1] for i in free_idx], dtype=np.float64)
    free_t = np.array([ts[i] for i in free_idx], dtype=np.float64)

    S00 = float(np.sum(free0))
    S01 = float(np.sum(free1))
    S10 = float(np.sum(free_t * free0))
    S11 = float(np.sum(free_t * free1))

    # Solve linear system for pair tildes:
    # a0 + b0 = -S00
    # ta*a0 + tb*b0 = -S10
    # similarly for component 1
    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0
    b1 = (-S11 + ta * S01) / denom
    a1 = -S01 - b1

    T[a, 0], T[a, 1] = a0, a1
    T[b, 0], T[b, 1] = b0, b1
    return T


def compute_z(ts, T, tx, ty):
    """
    z_ij = ([j i] / <i j>) * C_i*C_j
    where <i j> = t_j - t_i for lambda=[1,t],
    and [j i] = det(tilde_j, tilde_i).
    """
    t = ts
    C = (tx - t) * (ty - t)  # C_i

    br = T[EDGE_J, 0] * T[EDGE_I, 1] - T[EDGE_J, 1] * T[EDGE_I, 0]  # [j i]
    ang = (t[EDGE_J] - t[EDGE_I])                                   # <i j>
    ang = np.where(np.abs(ang) < 1e-12, 1e-12, ang)

    z = (br / ang) * (C[EDGE_I] * C[EDGE_J])
    return z


def loss_pos(ts, pair, free_idx, u_free, tx, ty):
    T = reconstruct_full_tildes(ts, pair, free_idx, u_free)
    z = compute_z(ts, T, tx, ty)

    pos_count = int(np.count_nonzero(z > POS_EPS))
    zmin = float(np.min(z))

    abs_z = np.abs(z)
    scale = float(np.median(abs_z) + SCALE_EPS)

    loss = float(np.sum(_softplus_vec(-z, scale)) + BETA * np.sum(_softplus_vec((MARGIN - z), scale)))
    return loss, pos_count, zmin


def init_tx_ty_interval_bruteforce(ts, pair, u_free, rng, tries=60):
    """
    Interval brute-force init for (tx,ty).
    Samples tx/ty from moment-curve intervals to find better initial C_i sign patterns.
    """
    t = ts
    intervals = [(-np.inf, t[0])] + [(t[k], t[k + 1]) for k in range(N - 1)] + [(t[-1], np.inf)]

    free_idx = [i for i in range(N) if i not in pair]
    best = None
    best_key = None

    for _ in range(tries):
        a = int(rng.integers(0, N + 1))
        b = int(rng.integers(0, N + 1))
        tx = sample_in_interval(rng, intervals[a][0], intervals[a][1])
        ty = sample_in_interval(rng, intervals[b][0], intervals[b][1])
        if tx >= ty - 1e-6:
            continue

        loss, pos_count, zmin = loss_pos(ts, pair, free_idx, u_free, tx, ty)
        key = (pos_count, -loss)  # prioritize positivity then loss
        if (best is None) or (key > best_key):
            best = (tx, ty, loss, pos_count, zmin)
            best_key = key

    if best is None:
        tx = float(ts[0] - rng.uniform(1.0, 3.0))
        ty = float(ts[-1] + rng.uniform(1.0, 3.0))
        return tx, ty

    return float(best[0]), float(best[1])


def simulated_annealing_restart(ts, pair, rng, stop_event=None):
    free_idx = [i for i in range(N) if i not in pair]

    u = rng.uniform(-5.0, 5.0, size=(len(free_idx),)).astype(np.float64)
    tx, ty = init_tx_ty_interval_bruteforce(ts, pair, u, rng, tries=60)

    loss, pos_count, zmin = loss_pos(ts, pair, free_idx, u, tx, ty)

    best = {
        "ts": ts.tolist(),
        "pair": [int(pair[0]), int(pair[1])],
        "free_idx": [int(i) for i in free_idx],
        "u_free": u.tolist(),
        "tx": float(tx),
        "ty": float(ty),
        "loss": float(loss),
        "pos_count": int(pos_count),
        "zmin": float(zmin),
    }

    if pos_count == 15:
        return best

    T0 = max(1.0, loss)
    sigma_u = float(SIGMA_U_INIT)
    sigma_xy = float(SIGMA_XY_INIT)

    accept_ct = 0
    trial_ct = 0

    for step in range(1, STEPS_PER_RESTART + 1):
        if stop_event is not None and stop_event.is_set():
            return best

        frac = step / float(STEPS_PER_RESTART)
        temp = T0 * (T_END / T0) ** frac

        # propose
        u2 = u + rng.normal(0.0, sigma_u, size=u.shape)
        tx2 = tx + float(rng.normal(0.0, sigma_xy))
        ty2 = ty + float(rng.normal(0.0, sigma_xy))

        # enforce ordering tx<ty
        if tx2 >= ty2 - 1e-6:
            m = 0.5 * (tx2 + ty2)
            tx2 = m - 1e-3
            ty2 = m + 1e-3

        loss2, pos2, zmin2 = loss_pos(ts, pair, free_idx, u2, tx2, ty2)
        dE = loss2 - loss

        trial_ct += 1
        if (dE <= 0.0) or (rng.random() < math.exp(-dE / max(temp, 1e-12))):
            u, tx, ty, loss, pos_count, zmin = u2, tx2, ty2, loss2, pos2, zmin2
            accept_ct += 1

            if (pos_count > best["pos_count"]) or (pos_count == best["pos_count"] and loss < best["loss"]):
                best.update({
                    "u_free": u.tolist(),
                    "tx": float(tx),
                    "ty": float(ty),
                    "loss": float(loss),
                    "pos_count": int(pos_count),
                    "zmin": float(zmin),
                })
                if pos_count == 15:
                    return best

        # mild step-size adaptation
        if step % 200 == 0:
            acc = accept_ct / max(trial_ct, 1)
            if acc < 0.20:
                sigma_u *= 0.85
                sigma_xy *= 0.85
            elif acc > 0.55:
                sigma_u *= 1.15
                sigma_xy *= 1.15
            accept_ct = 0
            trial_ct = 0

    return best


def worker_chunk(task_id, restarts, seed, stop_event=None):
    rng = np.random.default_rng(seed)

    best_overall = None
    for _ in range(restarts):
        if stop_event is not None and stop_event.is_set():
            break

        try:
            ts = sample_ts(rng)
            pair = pick_pair(rng)
            out = simulated_annealing_restart(ts, pair, rng, stop_event=stop_event)
        except Exception:
            # keep worker alive even if a rare numeric corner occurs
            continue

        if best_overall is None:
            best_overall = out
        else:
            if (out["pos_count"] > best_overall["pos_count"]) or (
                out["pos_count"] == best_overall["pos_count"] and out["loss"] < best_overall["loss"]
            ):
                best_overall = out

        if out["pos_count"] == 15 and out["zmin"] >= ZMIN_TRIGGER:
            return {"found": True, "candidate": out, "best": best_overall}

    return {"found": False, "candidate": None, "best": best_overall}


def check_active_charts_and_save(candidate):
    """
    Heavy check (Sage):
      - reconstruct tildes using solve_conservation_pair over RealField(200)
      - verify 15/15 positivity in hi precision
      - verify active charts under require_polytope=True
      - save RESULTS/golden_kinematics.json
    """
    from sage.all import RealField, vector, load

    # Lazy load (keeps hot loop light)
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
        ang = br(lambdas_sage[i], lambdas_sage[j])   # <i j>
        sq  = br(tildes_recon[j], tildes_recon[i])   # [j i]
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

    # Chart check via evaluator
    tf_final = {int(i): [1.0, float(u_free[k])] for k, i in enumerate(free_idx)}
    xs_final = [1.0, float(tx)]
    ys_final = [1.0, float(ty)]

    active_charts = 0
    hit_roots = []
    for roots in ALL_ROOTS:
        val = evaluator.evaluate_chart_fiber_gauge(
            roots, ts, tf_final, xs_final, ys_final, require_polytope=True
        )
        if abs(val) > 1e-20:
            active_charts += 1
            hit_roots.append(roots)

    if active_charts <= 0:
        print("  [reject] 15/15 but no active charts under polytope requirement.")
        return False

    os.makedirs("RESULTS", exist_ok=True)
    out = {
        "ts": ts,
        "pair": [int(pair[0]), int(pair[1])],
        "free_idx": [int(i) for i in free_idx],
        "u_free": [float(x) for x in u_free],
        "tildes_free": {str(i): tf_final[i] for i in tf_final},  # JSON keys become strings
        "x_s": xs_final,
        "y_s": ys_final,
        "z_values": z_vals,
        "zmin": float(zmin),
        "active_charts": int(active_charts),
        "active_roots": [list(r) for r in hit_roots],
        "meta": {
            "TOTAL_RESTARTS": int(TOTAL_RESTARTS),
            "STEPS_PER_RESTART": int(STEPS_PER_RESTART),
            "CHUNK_RESTARTS": int(CHUNK_RESTARTS),
            "POS_EPS": float(POS_EPS),
            "ZMIN_TRIGGER": float(ZMIN_TRIGGER),
            "timestamp": float(time.time()),
        },
    }

    with open("RESULTS/golden_kinematics.json", "w") as f:
        json.dump(out, f, indent=2)

    print(f"  ✨ GOLDEN FIND! Active Charts: {active_charts}")
    print("  Saved: RESULTS/golden_kinematics.json")
    return True


def serial_search():
    print("[fallback] Running SERIAL search (parallel pool unavailable).")
    rng = np.random.default_rng(SEED)
    best_global = {"pos_count": 0, "loss": float("inf")}
    started = time.time()

    for k in range(TOTAL_RESTARTS):
        ts = sample_ts(rng)
        pair = pick_pair(rng)
        cand = simulated_annealing_restart(ts, pair, rng, stop_event=None)

        if (cand["pos_count"] > best_global["pos_count"]) or (
            cand["pos_count"] == best_global["pos_count"] and cand["loss"] < best_global["loss"]
        ):
            best_global = {"pos_count": int(cand["pos_count"]), "loss": float(cand["loss"])}
            print(f"[progress] best_pos={best_global['pos_count']}/15 best_loss={best_global['loss']:.3e} t={time.time()-started:.1f}s")

        if cand["pos_count"] == 15 and cand["zmin"] >= ZMIN_TRIGGER:
            print(f"[candidate] Found 15/15 (loss={cand['loss']:.3e}, zmin={cand['zmin']:.3e})")
            if check_active_charts_and_save(cand):
                return

    print("Search exhausted with no golden kinematics found.")
    print(f"Best observed: {best_global['pos_count']}/15 (loss={best_global['loss']:.3e})")


def parallel_search():
    print("Starting FAST Parallel SA search for all-positive z_ij...")
    print(f"  TOTAL_RESTARTS={TOTAL_RESTARTS}, STEPS_PER_RESTART={STEPS_PER_RESTART}")
    print(f"  MAX_WORKERS={MAX_WORKERS}, CHUNK_RESTARTS={CHUNK_RESTARTS}")
    os.makedirs("RESULTS", exist_ok=True)

    import multiprocessing as mp

    num_chunks = int(math.ceil(TOTAL_RESTARTS / float(CHUNK_RESTARTS)))
    workers = min(MAX_WORKERS, os.cpu_count() or 1)

    best_global = {"pos_count": 0, "loss": float("inf")}
    started = time.time()

    mgr = mp.Manager()
    stop_event = mgr.Event()

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = []
        for k in range(num_chunks):
            restarts = CHUNK_RESTARTS if (k < num_chunks - 1) else (TOTAL_RESTARTS - CHUNK_RESTARTS * (num_chunks - 1))
            futures.append(ex.submit(worker_chunk, k, restarts, int(SEED + 1000 * k), stop_event))

        done_chunks = 0
        for fut in as_completed(futures):
            done_chunks += 1
            res = fut.result()

            b = res.get("best", None)
            if b is not None:
                if (b["pos_count"] > best_global["pos_count"]) or (
                    b["pos_count"] == best_global["pos_count"] and b["loss"] < best_global["loss"]
                ):
                    best_global = {"pos_count": int(b["pos_count"]), "loss": float(b["loss"])}
                    elapsed = time.time() - started
                    print(f"[progress] best_pos={best_global['pos_count']}/15 best_loss={best_global['loss']:.3e} t={elapsed:.1f}s")

            if res["found"]:
                cand = res["candidate"]
                print(f"[candidate] Found 15/15 in worker chunk. loss={cand['loss']:.3e} zmin={cand['zmin']:.3e}")

                stop_event.set()

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
                    print("[continue] Candidate rejected by active-charts/polytope check. Continuing...")

            if PRINT_EVERY_CHUNK and (done_chunks % PRINT_EVERY_CHUNK == 0):
                elapsed = time.time() - started
                print(f"[status] chunks_done={done_chunks}/{num_chunks}, elapsed={elapsed:.1f}s")

    print("Search exhausted with no golden kinematics found.")
    print(f"Best observed: {best_global['pos_count']}/15 (loss={best_global['loss']:.3e})")


def main():
    # On some CoCalc kernels / restricted environments, process pools can fail.
    # If that happens, fall back to serial search automatically.
    try:
        parallel_search()
    except Exception as e:
        print(f"[warn] Parallel search failed ({type(e).__name__}: {e}). Falling back to serial.")
        serial_search()


if __name__ == "__main__":
    main()

