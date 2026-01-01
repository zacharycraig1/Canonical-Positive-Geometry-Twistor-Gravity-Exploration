# src/atlas/search_all_positive_z_FAST_PARALLEL_v7_ALL_OPT.sage
#
# v7: "All the optimizations that should work" (CPU-first, optional GPU batch eval)
#
# What’s new vs v6:
#  A) Two-level scoring (MAJOR practical speedup):
#     - Always compute pos_count/zmin/neg_violation cheaply for all (u,combo).
#     - Only compute full Softplus loss (median scale + softplus) when pos_max >= LOSS_POS_THRESHOLD.
#     - For low pos regions, select using a fast proxy objective (sum of negative violations).
#
#  B) Refine (tx,ty) *within* the chosen interval pair:
#     - After choosing best combo (ia,ib), jitter tx,ty inside those intervals a few times
#       and keep the best by (pos, score). This often increases zmin and helps SA converge.
#
#  C) Near-miss intensification:
#     - If pos>=14, run additional short polish runs with perturbed starts.
#
#  D) Optional GPU acceleration (if CuPy is installed):
#     - Set USE_GPU=True to run **serial GPU mode** (1 worker) to avoid GPU contention.
#     - In GPU mode, the heavy CEM batch eval runs on GPU; correctness checks unchanged.
#
# Tuned defaults for ~16 cores / 16GB RAM.
#
# Output:
#   RESULTS/golden_kinematics.json
#
# Run:
#   sage src/atlas/search_all_positive_z_FAST_PARALLEL_v7_ALL_OPT.sage

import os
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
# USER TOGGLES
# =========================
USE_GPU = False         # if True and cupy installed -> run serial GPU mode
GPU_CEM_ONLY = True     # keep SA on CPU even in GPU mode (usually faster / less overhead)

# =========================
# CONFIG (16C / 16GB)
# =========================
N = 6

TOTAL_RESTARTS      = 900
SA_STEPS            = 700
POLISH_STEPS        = 1400

CHUNK_RESTARTS      = 20
SEED                = int(3)
MAX_WORKERS         = 16   # Use 16 cores

T_END               = 1e-3

# Softplus loss
MARGIN              = 1e-6
BETA                = 0.5
SCALE_EPS           = 1e-12

# Energy shaping / screening
POS_EPS             = 0.0
ZMIN_TRIGGER        = 1e-12
POS_LAMBDA          = 6.0

# Two-level scoring: only compute full softplus when we're close enough
LOSS_POS_THRESHOLD  = 13

# SA stepsizes
SIGMA_U             = 0.9
SIGMA_XY            = 0.6
POLISH_SIGMA_SCALE  = 0.35

# CEM seeding
CEM_POP             = 512
CEM_ITERS           = 10
CEM_ELITE_FRAC      = 0.15
CEM_STD_FLOOR       = 0.15
U_CLIP              = 12.0

# Interval refinement after picking (ia,ib)
REFINE_XY_TRIES     = 8
REFINE_XY_JITTER    = 0.12   # relative fraction of interval width for sampling (finite intervals)

# Near-miss intensification
EXTRA_POLISH_RUNS   = 3
EXTRA_POLISH_STEPS  = 600
EXTRA_SIGMA_SCALE   = 0.45

PRINT_EVERY_CHUNK   = 1

PAIR_BIAS = {
    (0, 1): 4.0,
    (4, 5): 4.0,
    (0, 5): 2.0,
    (1, 4): 2.0,
}

EDGE_PAIRS = [(i, j) for i in range(N) for j in range(i + 1, N)]
EDGE_I = np.array([i for (i, j) in EDGE_PAIRS], dtype=np.int64)
EDGE_J = np.array([j for (i, j) in EDGE_PAIRS], dtype=np.int64)
ALL_ROOTS = list(itertools.combinations(range(N), 3))  # 20 roots


# -------------------------
# Optional GPU backend
# -------------------------
try:
    import cupy as cp
    CUPY_OK = True
except Exception:
    cp = None
    CUPY_OK = False

def xp_backend(use_gpu):
    if use_gpu and CUPY_OK:
        return cp
    return np


def _softplus(x, scale, xp=np):
    xs = x / scale
    return scale * (xp.log1p(xp.exp(-xp.abs(xs))) + xp.maximum(xs, 0.0))


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


def interval_bounds(ts, idx):
    # intervals: 0=(-inf,t0), k=(t_{k-1},t_k), N=(t_{N-1},inf)
    if idx == 0:
        return (-np.inf, float(ts[0]))
    if idx == N:
        return (float(ts[N-1]), np.inf)
    return (float(ts[idx - 1]), float(ts[idx]))


def sample_tx_ty_in_intervals(ts, ia, ib, rng):
    lo_x, hi_x = interval_bounds(ts, ia)
    lo_y, hi_y = interval_bounds(ts, ib)

    def sample_in(lo, hi, side=2.0):
        if not np.isfinite(lo) and np.isfinite(hi):
            return float(hi - rng.uniform(0.5, side))
        if np.isfinite(lo) and not np.isfinite(hi):
            return float(lo + rng.uniform(0.5, side))
        if hi - lo < 1e-9:
            return float(lo)
        # jitter towards middle of interval
        width = hi - lo
        center = 0.5 * (lo + hi)
        rad = max(1e-6, REFINE_XY_JITTER * width)
        return float(rng.uniform(center - rad, center + rad))

    tx = sample_in(lo_x, hi_x)
    ty = sample_in(lo_y, hi_y)
    if tx >= ty - 1e-6:
        m = 0.5 * (tx + ty)
        tx = m - 1e-3
        ty = m + 1e-3
    return tx, ty


def interval_mid(ts, idx):
    if idx == 0:
        return float(ts[0] - 2.0)
    if idx == N:
        return float(ts[N-1] + 2.0)
    return float(0.5 * (ts[idx - 1] + ts[idx]))


def build_restart_struct(ts, pair, use_gpu=False):
    free_idx = [i for i in range(N) if i not in pair]
    a, b = int(pair[0]), int(pair[1])

    t = ts.astype(np.float64)
    ta, tb = float(t[a]), float(t[b])
    denom = tb - ta
    if abs(denom) < 1e-12:
        denom = 1e-12

    # gauge: free tilde0 = 1
    S00 = float(len(free_idx))
    S10 = float(np.sum(t[free_idx]))

    b0 = (-S10 + ta * S00) / denom
    a0 = -S00 - b0

    t0_const = np.zeros(N, dtype=np.float64)
    t0_const[free_idx] = 1.0
    t0_const[a] = a0
    t0_const[b] = b0

    free_t = t[free_idx].copy()

    ang = (t[EDGE_J] - t[EDGE_I])
    ang = np.where(np.abs(ang) < 1e-12, 1e-12, ang)
    ang_inv = 1.0 / ang

    mids = np.array([interval_mid(ts, k) for k in range(N + 1)], dtype=np.float64)

    tx_list = []
    ty_list = []
    ia_list = []
    ib_list = []
    Cprod_list = []

    for ia in range(N + 1):
        tx = float(mids[ia])
        for ib in range(ia + 1, N + 1):
            ty = float(mids[ib])
            if tx >= ty - 1e-6:
                continue
            C = (tx - t) * (ty - t)
            Cprod_edge = (C[EDGE_I] * C[EDGE_J]).astype(np.float64)
            ia_list.append(ia); ib_list.append(ib)
            tx_list.append(tx); ty_list.append(ty)
            Cprod_list.append(Cprod_edge)

    TX = np.array(tx_list, dtype=np.float64)
    TY = np.array(ty_list, dtype=np.float64)
    IA = np.array(ia_list, dtype=np.int16)
    IB = np.array(ib_list, dtype=np.int16)
    CPROD = np.stack(Cprod_list, axis=0)  # (M,15)
    M = CPROD.shape[0]

    if use_gpu:
        xp = cp
        # Move constant arrays to GPU once per restart.
        t0_const_x = xp.asarray(t0_const)
        ang_inv_x = xp.asarray(ang_inv)
        CPROD_x = xp.asarray(CPROD)
        free_t_x = xp.asarray(free_t)
    else:
        xp = np
        t0_const_x = t0_const
        ang_inv_x = ang_inv
        CPROD_x = CPROD
        free_t_x = free_t

    return {
        "ts": ts,
        "pair": pair,
        "free_idx": free_idx,
        "a": a,
        "b": b,
        "t": t,
        "ta": ta,
        "tb": tb,
        "denom": denom,
        "t0_const": t0_const,           # CPU copy (for SA eval)
        "ang_inv": ang_inv,             # CPU copy
        "free_t": free_t,               # CPU copy
        "TX": TX, "TY": TY, "IA": IA, "IB": IB,
        "CPROD": CPROD, "M": M,
        # GPU const arrays (or CPU refs)
        "xp": xp,
        "t0_const_x": t0_const_x,
        "ang_inv_x": ang_inv_x,
        "CPROD_x": CPROD_x,
        "free_t_x": free_t_x,
    }


def eval_batch_two_level(struct, U):
    """
    Evaluate U (K,4) across all combos (M) using:
      - Always compute pos and proxy_violation
      - Compute full softplus loss only when pos_max >= LOSS_POS_THRESHOLD
    Returns best-per-sample: best_pos, best_score, best_zmin, best_combo
    where score is full loss when close, else proxy violation.
    """
    xp = struct["xp"]
    free_idx = struct["free_idx"]
    a, b = struct["a"], struct["b"]
    ta, denom = struct["ta"], struct["denom"]

    free_t = struct["free_t_x"]
    t0_const = struct["t0_const_x"]
    ang_inv = struct["ang_inv_x"]
    CPROD = struct["CPROD_x"]
    M = struct["M"]

    K = U.shape[0]
    Ux = xp.asarray(U) if xp is cp else U

    S01 = xp.sum(Ux, axis=1)
    S11 = xp.sum(Ux * free_t[None, :], axis=1)

    b1 = (-S11 + ta * S01) / denom
    a1 = -S01 - b1

    # Build t1: (K,N)
    t1 = xp.zeros((K, N), dtype=xp.float64)
    t1[:, free_idx] = Ux
    t1[:, a] = a1
    t1[:, b] = b1

    br = t0_const[EDGE_J][None, :] * t1[:, EDGE_I] - t1[:, EDGE_J] * t0_const[EDGE_I][None, :]
    z = (br[:, None, :] * ang_inv[None, None, :]) * CPROD[None, :, :]   # (K,M,15)

    pos = xp.sum(z > POS_EPS, axis=2).astype(xp.int16)   # (K,M)
    zmin = xp.min(z, axis=2)                             # (K,M)

    # proxy violation: sum of negative parts (cheap)
    neg = xp.maximum(-z, 0.0)
    proxy = xp.sum(neg, axis=2)                          # (K,M)

    pos_max = xp.max(pos, axis=1)                        # (K,)
    mask = (pos == pos_max[:, None])

    # score init = proxy (but only on max-pos combos)
    INF = xp.asarray(np.inf) if xp is cp else np.inf
    score = xp.where(mask, proxy, INF)

    # Full loss only for samples with pos_max >= threshold
    if LOSS_POS_THRESHOLD is not None:
        good = (pos_max >= LOSS_POS_THRESHOLD)
    else:
        good = xp.ones_like(pos_max, dtype=bool)

    if xp is cp:
        good_cpu = cp.asnumpy(good)
    else:
        good_cpu = good

    # For each combo, compute full loss on the subset that needs it
    if np.any(good_cpu):
        # Move z to CPU only if using GPU_CEM_ONLY? No: keep on GPU if cupy, compute there.
        for m in range(M):
            if xp is cp:
                mk = mask[:, m] & good
                if not bool(cp.any(mk)):
                    continue
                z_m = z[mk, m, :]  # (Kk,15) on GPU
                abs_z = xp.abs(z_m)
                med = xp.partition(abs_z, abs_z.shape[1] // 2, axis=1)[:, abs_z.shape[1] // 2]
                scale = med + SCALE_EPS
                loss_m = xp.sum(_softplus(-z_m, scale[:, None], xp=xp), axis=1) + BETA * xp.sum(_softplus((MARGIN - z_m), scale[:, None], xp=xp), axis=1)
                score = score.copy()
                score[mk, m] = loss_m
            else:
                mk = mask[:, m] & good
                if not np.any(mk):
                    continue
                z_m = z[mk, m, :]  # (Kk,15) CPU
                abs_z = np.abs(z_m)
                med = np.partition(abs_z, abs_z.shape[1] // 2, axis=1)[:, abs_z.shape[1] // 2]
                scale = med + SCALE_EPS
                loss_m = np.sum(_softplus(-z_m, scale[:, None], xp=np), axis=1) + BETA * np.sum(_softplus((MARGIN - z_m), scale[:, None], xp=np), axis=1)
                score[mk, m] = loss_m

    best_combo = xp.argmin(score, axis=1).astype(xp.int32)
    idx = xp.arange(K)
    best_pos = pos[idx, best_combo].astype(xp.int32)
    best_score = score[idx, best_combo]
    best_zmin = zmin[idx, best_combo]

    if xp is cp:
        return (cp.asnumpy(best_pos), cp.asnumpy(best_score), cp.asnumpy(best_zmin), cp.asnumpy(best_combo))
    else:
        return (best_pos, best_score, best_zmin, best_combo)


def cem_seed(struct, rng):
    d = len(struct["free_idx"])
    assert d == 4

    mu = np.zeros((d,), dtype=np.float64)
    sig = np.ones((d,), dtype=np.float64) * 3.0

    best = None  # (pos, score, zmin, u, combo)

    elite_n = max(8, int(CEM_ELITE_FRAC * CEM_POP))

    for it in range(CEM_ITERS):
        U = mu[None, :] + sig[None, :] * rng.standard_normal(size=(CEM_POP, d))
        U = np.clip(U, -U_CLIP, U_CLIP)

        pos, score, zmin, combo = eval_batch_two_level(struct, U)

        # best sample in population: pos desc, score asc
        idx_best = np.lexsort((score, -pos))[0]
        cand = (int(pos[idx_best]), float(score[idx_best]), float(zmin[idx_best]), U[idx_best].copy(), int(combo[idx_best]))
        if best is None or (cand[0] > best[0]) or (cand[0] == best[0] and cand[1] < best[1]):
            best = cand
            if best[0] == 15:
                break

        # elites based on (pos, score)
        order = np.lexsort((score, -pos))
        elites = U[order[:elite_n], :]

        mu = np.mean(elites, axis=0)
        sig = np.std(elites, axis=0) + CEM_STD_FLOOR

    posB, scoreB, zminB, uB, comboB = best

    # Refine (tx,ty) inside selected intervals (ia,ib) for this uB
    ia = int(struct["IA"][comboB])
    ib = int(struct["IB"][comboB])
    tx0 = float(struct["TX"][comboB])
    ty0 = float(struct["TY"][comboB])

    # Small jitter refinement inside interval bounds
    best_tx, best_ty = tx0, ty0
    best_pos = posB
    best_score = scoreB
    best_zmin = zminB

    # Fast local scoring uses CPU eval_single (proxy if far) for stability in SA phase
    for _ in range(REFINE_XY_TRIES):
        tx, ty = sample_tx_ty_in_intervals(struct["ts"], ia, ib, rng)
        # quick proxy: compute pos + proxy (sum negative) (no full softplus unless close)
        p, sc, zm = quick_score_cpu(struct, uB, tx, ty)
        if (p > best_pos) or (p == best_pos and sc < best_score):
            best_pos, best_score, best_zmin = p, sc, zm
            best_tx, best_ty = tx, ty
        if best_pos == 15:
            break

    return uB, best_tx, best_ty, best_score, best_pos, best_zmin, comboB


def quick_score_cpu(struct, u, tx, ty):
    """
    Cheap score on CPU:
      - compute z
      - return (pos, score, zmin)
      - score is proxy if pos < LOSS_POS_THRESHOLD else full softplus loss
    """
    free_idx = struct["free_idx"]
    a, b = struct["a"], struct["b"]
    ta, denom = struct["ta"], struct["denom"]
    free_t = struct["free_t"]
    t0_const = struct["t0_const"]
    ang_inv = struct["ang_inv"]
    t = struct["t"]

    u = np.asarray(u, dtype=np.float64)
    S01 = float(np.sum(u))
    S11 = float(np.sum(u * free_t))

    b1 = (-S11 + ta * S01) / denom
    a1 = -S01 - b1

    t1 = np.zeros((N,), dtype=np.float64)
    t1[free_idx] = u
    t1[a] = a1
    t1[b] = b1

    C = (tx - t) * (ty - t)
    Cprod = (C[EDGE_I] * C[EDGE_J])

    br = t0_const[EDGE_J] * t1[EDGE_I] - t1[EDGE_J] * t0_const[EDGE_I]
    z = (br * ang_inv) * Cprod

    pos = int(np.count_nonzero(z > POS_EPS))
    zmin = float(np.min(z))

    if pos < LOSS_POS_THRESHOLD:
        score = float(np.sum(np.maximum(-z, 0.0)))
        return pos, score, zmin

    abs_z = np.abs(z)
    med = float(np.partition(abs_z, abs_z.size // 2)[abs_z.size // 2])
    scale = med + SCALE_EPS

    loss = float(np.sum(_softplus(-z, scale, xp=np)) + BETA * np.sum(_softplus((MARGIN - z), scale, xp=np)))
    return pos, loss, zmin


def eval_single(struct, u, tx, ty):
    """
    Full SA eval: energy uses full softplus loss (for stable polishing near feasibility),
    but we keep it cheap enough (15 numbers).
    """
    pos, loss_or_proxy, zmin = quick_score_cpu(struct, u, tx, ty)

    # If we're still far, convert proxy to a mild energy (still monotone in violations)
    if pos < LOSS_POS_THRESHOLD:
        scale = 1.0
        energy = float(loss_or_proxy + POS_LAMBDA * (15 - pos) * scale)
        return energy, float(loss_or_proxy), pos, zmin

    # near feasible: use loss
    abs_dummy_scale = 1.0
    energy = float(loss_or_proxy + POS_LAMBDA * (15 - pos) * abs_dummy_scale)
    return energy, float(loss_or_proxy), pos, zmin


def sa_polish(struct, rng, u0, tx0, ty0, stop_flag=None, steps=SA_STEPS, sigma_scale=1.0):
    u = np.asarray(u0, dtype=np.float64)
    tx, ty = float(tx0), float(ty0)

    energy, loss, pos, zmin = eval_single(struct, u, tx, ty)
    best = (pos, loss, zmin, u.copy(), tx, ty, energy)

    T0 = max(1.0, energy)
    sigma_u = SIGMA_U * sigma_scale
    sigma_xy = SIGMA_XY * sigma_scale

    accept_ct = 0
    trial_ct = 0

    for step in range(1, steps + 1):
        if stop_flag is not None and stop_flag.value == 1:
            break

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

        e2, l2, p2, z2 = eval_single(struct, u2, tx2, ty2)
        dE = e2 - energy

        trial_ct += 1
        if (dE <= 0.0) or (rng.random() < math.exp(-dE / max(temp, 1e-12))):
            u, tx, ty, energy, loss, pos, zmin = u2, tx2, ty2, e2, l2, p2, z2
            accept_ct += 1

            if (pos > best[0]) or (pos == best[0] and loss < best[1]):
                best = (pos, loss, zmin, u.copy(), tx, ty, energy)
                if pos == 15:
                    break

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

    posB, lossB, zminB, uB, txB, tyB, eB = best
    return {
        "pos_count": int(posB),
        "loss": float(lossB),
        "zmin": float(zminB),
        "u_free": uB.tolist(),
        "tx": float(txB),
        "ty": float(tyB),
        "energy": float(eB),
    }


def restart_once(rng, stop_flag=None, use_gpu=False):
    ts = sample_ts(rng)
    pair = pick_pair(rng)
    struct = build_restart_struct(ts, pair, use_gpu=use_gpu)

    u0, tx0, ty0, score0, pos0, zmin0, combo0 = cem_seed(struct, rng)

    out = sa_polish(struct, rng, u0, tx0, ty0, stop_flag=stop_flag, steps=SA_STEPS, sigma_scale=1.0)
    out.update({
        "ts": ts.tolist(),
        "pair": [int(pair[0]), int(pair[1])],
        "free_idx": [int(i) for i in struct["free_idx"]],
        "seed_combo_idx": int(combo0),
        "seed_combo_ia": int(struct["IA"][combo0]),
        "seed_combo_ib": int(struct["IB"][combo0]),
        "seed_score": float(score0),
    })

    if out["pos_count"] == 15:
        return out

    if out["pos_count"] >= 14:
        pol = sa_polish(struct, rng, np.asarray(out["u_free"]), out["tx"], out["ty"], stop_flag=stop_flag,
                        steps=POLISH_STEPS, sigma_scale=POLISH_SIGMA_SCALE)
        if (pol["pos_count"] > out["pos_count"]) or (pol["pos_count"] == out["pos_count"] and pol["loss"] < out["loss"]):
            out.update(pol)

    if out["pos_count"] >= 14 and out["pos_count"] < 15:
        u_base = np.asarray(out["u_free"], dtype=np.float64)
        tx_base = float(out["tx"]); ty_base = float(out["ty"])
        best = out
        for _ in range(EXTRA_POLISH_RUNS):
            u_try = u_base + rng.normal(0.0, 0.35, size=u_base.shape)
            tx_try = tx_base + float(rng.normal(0.0, 0.35))
            ty_try = ty_base + float(rng.normal(0.0, 0.35))
            if tx_try >= ty_try - 1e-6:
                m = 0.5*(tx_try+ty_try)
                tx_try = m - 1e-3
                ty_try = m + 1e-3
            pol2 = sa_polish(struct, rng, u_try, tx_try, ty_try, stop_flag=stop_flag,
                             steps=EXTRA_POLISH_STEPS, sigma_scale=EXTRA_SIGMA_SCALE)
            cand = dict(best)
            cand.update(pol2)
            if (cand["pos_count"] > best["pos_count"]) or (cand["pos_count"] == best["pos_count"] and cand["loss"] < best["loss"]):
                best = cand
            if best["pos_count"] == 15:
                return best
        return best

    return out


def worker_chunk(task_id, restarts, seed, stop_flag):
    rng = np.random.default_rng(seed)
    best_overall = None

    for _ in range(restarts):
        if stop_flag.value == 1:
            break
        try:
            out = restart_once(rng, stop_flag=stop_flag, use_gpu=False)
        except Exception:
            continue

        if best_overall is None:
            best_overall = out
        else:
            if (out["pos_count"] > best_overall["pos_count"]) or (out["pos_count"] == best_overall["pos_count"] and out["loss"] < best_overall["loss"]):
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
    out = dict(candidate)
    out["z_values"] = z_vals
    out["active_charts"] = int(active_charts)
    out["active_roots"] = [list(r) for r in hit_roots]
    out["meta"] = {
        "TOTAL_RESTARTS": int(TOTAL_RESTARTS),
        "SA_STEPS": int(SA_STEPS),
        "POLISH_STEPS": int(POLISH_STEPS),
        "CHUNK_RESTARTS": int(CHUNK_RESTARTS),
        "CEM_POP": int(CEM_POP),
        "CEM_ITERS": int(CEM_ITERS),
        "CEM_ELITE_FRAC": float(CEM_ELITE_FRAC),
        "LOSS_POS_THRESHOLD": int(LOSS_POS_THRESHOLD),
        "REFINE_XY_TRIES": int(REFINE_XY_TRIES),
        "EXTRA_POLISH_RUNS": int(EXTRA_POLISH_RUNS),
        "USE_GPU": bool(USE_GPU),
        "timestamp": float(time.time()),
    }

    with open("RESULTS/golden_kinematics.json", "w") as f:
        json.dump(out, f, indent=2)

    print(f"  ✨ GOLDEN FIND! Active Charts: {active_charts}")
    print("  Saved: RESULTS/golden_kinematics.json")
    return True


def cpu_parallel_search():
    print("Starting CPU parallel search (v7 ALL OPT) for all-positive z_ij...")
    print(f"  TOTAL_RESTARTS={TOTAL_RESTARTS}, SA_STEPS={SA_STEPS}, POLISH_STEPS={POLISH_STEPS}")
    print(f"  CEM_POP={CEM_POP}, CEM_ITERS={CEM_ITERS}, LOSS_POS_THRESHOLD={LOSS_POS_THRESHOLD}")
    print(f"  CHUNK_RESTARTS={CHUNK_RESTARTS}")

    os.makedirs("RESULTS", exist_ok=True)

    import multiprocessing as mp
    cores = os.cpu_count() or 1
    if MAX_WORKERS is None:
        workers = min(max(1, cores - 2), 12)
    else:
        workers = min(MAX_WORKERS, cores)

    # Use Manager for stop_flag to avoid pickle errors with mp.Value in some envs
    mgr = mp.Manager()
    stop_flag = mgr.Value('i', 0)
    
    num_chunks = int(math.ceil(TOTAL_RESTARTS / float(CHUNK_RESTARTS)))

    best_global = {"pos_count": 0, "loss": float("inf")}
    started = time.time()

    print(f"  cores={cores} -> workers={workers}")

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
                if (b["pos_count"] > best_global["pos_count"]) or (b["pos_count"] == best_global["pos_count"] and b["loss"] < best_global["loss"]):
                    best_global = {"pos_count": int(b["pos_count"]), "loss": float(b["loss"])}
                    elapsed = time.time() - started
                    print(f"[progress] best_pos={best_global['pos_count']}/15 best_loss={best_global['loss']:.3e} t={elapsed:.1f}s")

            if res["found"]:
                cand = res["candidate"]
                print(f"[candidate] Found 15/15 seed. loss={cand['loss']:.3e} zmin={cand['zmin']:.3e}")

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


def gpu_serial_search():
    if not CUPY_OK:
        print("[GPU] CuPy not available. Falling back to CPU parallel.")
        return cpu_parallel_search()

    print("Starting GPU-accelerated SERIAL search (v7 ALL OPT) ...")
    print("  NOTE: GPU mode uses 1 process to avoid GPU contention.")
    os.makedirs("RESULTS", exist_ok=True)

    rng = np.random.default_rng(SEED)
    best_global = {"pos_count": 0, "loss": float("inf")}
    started = time.time()

    for k in range(TOTAL_RESTARTS):
        out = restart_once(rng, stop_flag=None, use_gpu=True)

        if (out["pos_count"] > best_global["pos_count"]) or (out["pos_count"] == best_global["pos_count"] and out["loss"] < best_global["loss"]):
            best_global = {"pos_count": int(out["pos_count"]), "loss": float(out["loss"])}
            print(f"[progress] k={k} best_pos={best_global['pos_count']}/15 best_loss={best_global['loss']:.3e} t={time.time()-started:.1f}s")

        if out["pos_count"] == 15 and out["zmin"] >= ZMIN_TRIGGER:
            print(f"[candidate] Found 15/15 seed. loss={out['loss']:.3e} zmin={out['zmin']:.3e}")
            if check_active_charts_and_save(out):
                return

    print("GPU serial exhausted.")
    print(f"Best observed: {best_global['pos_count']}/15 (loss={best_global['loss']:.3e})")


def main():
    if USE_GPU:
        gpu_serial_search()
    else:
        cpu_parallel_search()


if __name__ == "__main__":
    main()


