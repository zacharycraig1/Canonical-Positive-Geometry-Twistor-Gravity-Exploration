# src/atlas/search_all_positive_z_FAST_PARALLEL_v4.sage
#
# FAST + PARALLEL Simulated Annealing search for n=6 all-positive z_ij (15/15),
# with EXACT momentum conservation solved in closed form at every step.
#
# v4 (tuned for ~16 cores / 16GB RAM):
# - Auto-tuned MAX_WORKERS default (safe for RAM) if left as None.
# - Replaces Manager().Event with shared mp.Value stop flag (lower overhead).
# - Reuses a preallocated t1 array inside the per-restart evaluator (reduces allocations in hot loop).
# - Multi-seed start: tries several u_free seeds and keeps the best before SA (higher hit-rate per restart).
# - Deterministic scan over all interval pairs for (tx,ty) (targets correct C-sign pattern quickly).
# - Energy shaping encourages increasing pos_count quickly (final correctness unchanged).
# - "Polish mode" automatically triggers extra refinement steps when reaching pos_count>=14.
# - Forces BLAS/OpenMP threads=1 per worker (prevents oversubscription).
#
# Output:
#   RESULTS/golden_kinematics.json
#
# Run:
#   sage src/atlas/search_all_positive_z_FAST_PARALLEL_v4.sage

import os

# ---- Prevent BLAS/OpenMP oversubscription in multi-process mode ----
for _k in [
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"
]:
    os.environ.setdefault(_k, "1")

import json, math, time
import itertools
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

# =========================
# CONFIG
# =========================
N = 6

# Overall budget
TOTAL_RESTARTS     = 5000      # for 16 cores, more restarts often beats deeper single runs
STEPS_PER_RESTART  = 1200      # main SA steps per restart
CHUNK_RESTARTS     = 50        # reduces scheduler overhead; still early-stops fast

SEED               = int(3)

# If None, auto-tune for ~16GB RAM: use at most (cores-2) and cap at 12
MAX_WORKERS        = 16 # Use all 16 cores as requested

# SA schedule
T_END              = 1e-3

# Softplus loss
MARGIN             = 1e-6
BETA               = 0.5
SCALE_EPS          = 1e-12

# Proposals
SIGMA_U_INIT       = 0.9
SIGMA_XY_INIT      = 0.6

# Strict positivity threshold
POS_EPS            = 0.0

# Only forward 15/15 candidates with zmin >= trigger to the heavy Sage evaluator
ZMIN_TRIGGER       = 1e-12

# Energy shaping: encourages higher pos_count earlier
POS_LAMBDA         = 6.0

# Multi-seed start (try K different u seeds per restart; pick best)
U_SEEDS_PER_RESTART = 8

# Polish phase (when we hit 14/15 or better)
POLISH_STEPS       = 1800
POLISH_SIGMA_SCALE = 0.35

PRINT_EVERY_CHUNK  = 1

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
    ts = np.sort(rng.uniform(1.0, 10.0, size=N).astype(np.float64))
    for k in range(1, N):
        if ts[k] <= ts[k - 1]:
            ts[k] = ts[k - 1] + 1e-6
    return ts


def interval_bounds(ts, idx):
    if idx == 0:
        return (-np.inf, float(ts[0]))
    if idx == N:
        return (float(ts[N-1]), np.inf)
    return (float(ts[idx - 1]), float(ts[idx]))


def interval_mid(ts, idx):
    lo, hi = interval_bounds(ts, idx)
    if not np.isfinite(lo) and np.isfinite(hi):
        return float(hi - 2.0)
    if np.isfinite(lo) and not np.isfinite(hi):
        return float(lo + 2.0)
    return float(0.5 * (lo + hi))


def sample_in_interval(rng, lo, hi, base_scale=2.0):
    if not np.isfinite(lo) and np.isfinite(hi):
        return float(hi - rng.uniform(0.5, base_scale))
    if np.isfinite(lo) and not np.isfinite(hi):
        return float(lo + rng.uniform(0.5, base_scale))
    if hi - lo < 1e-9:
        return float(lo)
    return float(rng.uniform(lo, hi))


def make_evaluator(ts, pair):
    """
    Fast evaluator for fixed (ts, pair).
    Reuses a preallocated t1 array to reduce allocations in the hot loop.
    """
    free_idx = [i for i in range(N) if i not in pair]
    a, b = int(pair[0]), int(pair[1])
    t = ts.astype(np.float64)

    ta, tb = float(t[a]), float(t[b])
    denom = tb - ta
    if abs(denom) < 1e-12:
        denom = 1e-12

    # Gauge: free tilde0 = 1
    S00 = float(len(free_idx))
    S10 = float(np.sum(t[free_idx]))

    # pair tilde0 constants
    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0

    # t0_const
    t0_const = np.zeros(N, dtype=np.float64)
    t0_const[free_idx] = 1.0
    t0_const[a] = a0
    t0_const[b] = b0

    free_t = t[free_idx].copy()

    ang = (t[EDGE_J] - t[EDGE_I])
    ang = np.where(np.abs(ang) < 1e-12, 1e-12, ang)
    ang_inv = 1.0 / ang

    # reused buffers (safe: single-threaded within process)
    t1 = np.zeros(N, dtype=np.float64)

    def eval_all(u_free, tx, ty):
        u = np.asarray(u_free, dtype=np.float64)

        # S01,S11 depend on u
        S01 = float(np.sum(u))
        S11 = float(np.sum(free_t * u))

        # pair tilde1
        b1 = (-S11 + ta * S01) / denom
        a1 = -S01 - b1

        # fill t1
        t1.fill(0.0)
        t1[free_idx] = u
        t1[a] = a1
        t1[b] = b1

        C = (tx - t) * (ty - t)
        Cprod = (C[EDGE_I] * C[EDGE_J])

        br = t0_const[EDGE_J] * t1[EDGE_I] - t1[EDGE_J] * t0_const[EDGE_I]
        z = (br * ang_inv) * Cprod

        pos_count = int(np.count_nonzero(z > POS_EPS))
        zmin = float(np.min(z))

        abs_z = np.abs(z)
        med = float(np.partition(abs_z, abs_z.size // 2)[abs_z.size // 2])
        scale = med + SCALE_EPS

        loss = float(np.sum(_softplus_vec(-z, scale)) + BETA * np.sum(_softplus_vec((MARGIN - z), scale)))
        energy = loss + POS_LAMBDA * (15 - pos_count) * scale
        return energy, loss, pos_count, zmin

    return free_idx, eval_all


def best_tx_ty_over_all_intervals(ts, rng, eval_all, u_free):
    """
    Scan all interval-pairs for tx<ty (fast; 21 combos for n=6) and choose best by (pos_count, loss).
    Also includes a few within-same-interval seeds.
    """
    best = None
    best_key = None

    mids = [interval_mid(ts, k) for k in range(N + 1)]

    # non-empty blocks ia<ib
    for ia in range(0, N + 1):
        tx0 = mids[ia]
        for ib in range(ia + 1, N + 1):
            ty0 = mids[ib]
            if tx0 >= ty0 - 1e-6:
                continue
            energy, loss, pos_count, zmin = eval_all(u_free, tx0, ty0)
            key = (pos_count, -loss)
            if (best is None) or (key > best_key):
                best = (tx0, ty0, energy, loss, pos_count, zmin)
                best_key = key

    # same-interval seeds
    for ia in range(0, N + 1):
        lo, hi = interval_bounds(ts, ia)
        tx = sample_in_interval(rng, lo, hi)
        ty = sample_in_interval(rng, lo, hi)
        if tx >= ty:
            tx, ty = min(tx, ty), max(tx, ty)
        if ty - tx < 1e-6:
            ty = tx + 1e-3
        energy, loss, pos_count, zmin = eval_all(u_free, tx, ty)
        key = (pos_count, -loss)
        if (best is None) or (key > best_key):
            best = (tx, ty, energy, loss, pos_count, zmin)
            best_key = key

    if best is None:
        tx = float(ts[0] - rng.uniform(1.0, 3.0))
        ty = float(ts[-1] + rng.uniform(1.0, 3.0))
        energy, loss, pos_count, zmin = eval_all(u_free, tx, ty)
        return tx, ty, energy, loss, pos_count, zmin

    return best


def choose_best_u_seed(ts, pair, rng):
    """
    Multi-seed initializer: sample several u seeds, for each choose best (tx,ty) via interval scan,
    then pick the overall best seed by (pos_count, loss).
    """
    free_idx, eval_all = make_evaluator(ts, pair)
    best = None
    best_key = None

    for _ in range(U_SEEDS_PER_RESTART):
        u = rng.uniform(-5.0, 5.0, size=(len(free_idx),)).astype(np.float64)
        tx, ty, energy, loss, pos_count, zmin = best_tx_ty_over_all_intervals(ts, rng, eval_all, u)
        key = (pos_count, -loss)
        if (best is None) or (key > best_key):
            best = (u, tx, ty, energy, loss, pos_count, zmin, free_idx, eval_all)
            best_key = key

        # early exit if already perfect
        if best is not None and best[5] == 15:
            break

    u, tx, ty, energy, loss, pos_count, zmin, free_idx, eval_all = best
    return free_idx, eval_all, u, tx, ty, energy, loss, pos_count, zmin


def sa_core(free_idx, eval_all, rng, u, tx, ty, energy, loss, pos_count, zmin, stop_flag=None, steps=STEPS_PER_RESTART, sigma_scale=1.0):
    """
    Core SA loop (shared by main and polish). Returns best dict for this phase.
    """
    best = {
        "u_free": u.tolist(),
        "tx": float(tx),
        "ty": float(ty),
        "energy": float(energy),
        "loss": float(loss),
        "pos_count": int(pos_count),
        "zmin": float(zmin),
    }

    T0 = max(1.0, energy)
    sigma_u = float(SIGMA_U_INIT) * sigma_scale
    sigma_xy = float(SIGMA_XY_INIT) * sigma_scale

    accept_ct = 0
    trial_ct = 0
    last_improve = 0

    for step in range(1, steps + 1):
        if stop_flag is not None and stop_flag.value == 1:
            return best

        frac = step / float(steps)
        temp = T0 * (T_END / T0) ** frac

        r = rng.random()
        u2 = u
        tx2, ty2 = tx, ty

        if r < 0.78:
            u2 = u + rng.normal(0.0, sigma_u, size=u.shape)
        elif r < 0.93:
            tx2 = tx + float(rng.normal(0.0, sigma_xy))
            ty2 = ty + float(rng.normal(0.0, sigma_xy))
        else:
            u2 = u + rng.normal(0.0, sigma_u, size=u.shape)
            tx2 = tx + float(rng.normal(0.0, sigma_xy))
            ty2 = ty + float(rng.normal(0.0, sigma_xy))

        if tx2 >= ty2 - 1e-6:
            m = 0.5 * (tx2 + ty2)
            tx2 = m - 1e-3
            ty2 = m + 1e-3

        energy2, loss2, pos2, zmin2 = eval_all(u2, tx2, ty2)
        dE = energy2 - energy

        trial_ct += 1
        if (dE <= 0.0) or (rng.random() < math.exp(-dE / max(temp, 1e-12))):
            u, tx, ty, energy, loss, pos_count, zmin = u2, tx2, ty2, energy2, loss2, pos2, zmin2
            accept_ct += 1

            if (pos_count > best["pos_count"]) or (pos_count == best["pos_count"] and loss < best["loss"]):
                best.update({
                    "u_free": u.tolist(),
                    "tx": float(tx),
                    "ty": float(ty),
                    "energy": float(energy),
                    "loss": float(loss),
                    "pos_count": int(pos_count),
                    "zmin": float(zmin),
                })
                last_improve = step
                if pos_count == 15:
                    return best

        if step % 250 == 0:
            acc = accept_ct / max(trial_ct, 1)
            if acc < 0.20:
                sigma_u *= 0.85
                sigma_xy *= 0.85
            elif acc > 0.55:
                sigma_u *= 1.15
                sigma_xy *= 1.15
            accept_ct = 0
            trial_ct = 0

        # macro jump on stagnation
        if (step - last_improve) >= 600 and best["pos_count"] <= 13:
            # keep u but reseed tx/ty via interval scan
            txJ, tyJ, eJ, lJ, pJ, zJ = best_tx_ty_over_all_intervals(ts_global, rng, eval_all, u)
            if pJ > pos_count or (pJ == pos_count and lJ < loss):
                tx, ty, energy, loss, pos_count, zmin = txJ, tyJ, eJ, lJ, pJ, zJ
            last_improve = step

        if best["pos_count"] >= 13 and step % 400 == 0:
            sigma_u *= 0.9
            sigma_xy *= 0.9

    return best


def simulated_annealing_restart(ts, pair, rng, stop_flag=None):
    global ts_global
    ts_global = ts  # used only for macro jump helper

    free_idx, eval_all, u, tx, ty, energy, loss, pos_count, zmin = choose_best_u_seed(ts, pair, rng)

    # main phase
    best_phase = sa_core(free_idx, eval_all, rng, u, tx, ty, energy, loss, pos_count, zmin, stop_flag=stop_flag, steps=STEPS_PER_RESTART, sigma_scale=1.0)

    out = {
        "ts": ts.tolist(),
        "pair": [int(pair[0]), int(pair[1])],
        "free_idx": [int(i) for i in free_idx],
        "u_free": best_phase["u_free"],
        "tx": best_phase["tx"],
        "ty": best_phase["ty"],
        "energy": best_phase["energy"],
        "loss": best_phase["loss"],
        "pos_count": best_phase["pos_count"],
        "zmin": best_phase["zmin"],
    }

    if out["pos_count"] == 15:
        return out

    # polish if close
    if out["pos_count"] >= 14:
        uP = np.asarray(out["u_free"], dtype=np.float64)
        txP, tyP = float(out["tx"]), float(out["ty"])
        eP, lP, pP, zP = eval_all(uP, txP, tyP)
        best_polish = sa_core(free_idx, eval_all, rng, uP, txP, tyP, eP, lP, pP, zP, stop_flag=stop_flag, steps=POLISH_STEPS, sigma_scale=POLISH_SIGMA_SCALE)

        if (best_polish["pos_count"] > out["pos_count"]) or (best_polish["pos_count"] == out["pos_count"] and best_polish["loss"] < out["loss"]):
            out.update({
                "u_free": best_polish["u_free"],
                "tx": best_polish["tx"],
                "ty": best_polish["ty"],
                "energy": best_polish["energy"],
                "loss": best_polish["loss"],
                "pos_count": best_polish["pos_count"],
                "zmin": best_polish["zmin"],
            })

    return out


def worker_chunk(task_id, restarts, seed, stop_flag):
    rng = np.random.default_rng(seed)
    best_overall = None

    for _ in range(restarts):
        if stop_flag.value == 1:
            break
        try:
            ts = sample_ts(rng)
            pair = pick_pair(rng)
            out = simulated_annealing_restart(ts, pair, rng, stop_flag=stop_flag)
        except Exception:
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
        "tildes_free": {str(i): tf_final[i] for i in tf_final},
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
            "POS_LAMBDA": float(POS_LAMBDA),
            "U_SEEDS_PER_RESTART": int(U_SEEDS_PER_RESTART),
            "POLISH_STEPS": int(POLISH_STEPS),
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

    for _ in range(TOTAL_RESTARTS):
        ts = sample_ts(rng)
        pair = pick_pair(rng)
        cand = simulated_annealing_restart(ts, pair, rng, stop_flag=None)

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
    print(f"  CHUNK_RESTARTS={CHUNK_RESTARTS}")

    os.makedirs("RESULTS", exist_ok=True)

    import multiprocessing as mp

    cores = os.cpu_count() or 1
    if MAX_WORKERS is None:
        # Safe default for ~16GB RAM: cap at 12 and leave 2 cores free
        workers = min(max(1, cores - 2), 12)
    else:
        workers = min(MAX_WORKERS, cores)

    # Shared stop flag using Manager (slower but safer for pickling)
    # mp.Value fails with fork/spawn issues in some envs unless passed carefully
    mgr = mp.Manager()
    stop_flag = mgr.Value('i', 0)

    num_chunks = int(math.ceil(TOTAL_RESTARTS / float(CHUNK_RESTARTS)))
    best_global = {"pos_count": 0, "loss": float("inf")}
    started = time.time()

    print(f"  cores={cores} -> workers={workers} (MAX_WORKERS={MAX_WORKERS})")

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = []
        for k in range(num_chunks):
            restarts = CHUNK_RESTARTS if (k < num_chunks - 1) else (TOTAL_RESTARTS - CHUNK_RESTARTS * (num_chunks - 1))
            futures.append(ex.submit(worker_chunk, k, restarts, int(SEED + 1000 * k), stop_flag))

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

                stop_flag.value = 1

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
    try:
        parallel_search()
    except Exception as e:
        print(f"[warn] Parallel search failed ({type(e).__name__}: {e}). Falling back to serial.")
        serial_search()


if __name__ == "__main__":
    main()

