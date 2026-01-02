import os
import json
import time
import numpy as np
import heapq
from sage.all import RealField, vector, load

# Load helpers
load("src/atlas/jacobian_fiber_gauge.sage")
load("src/atlas/solve_conservation_pair.sage")

# =========================
# CONFIG
# =========================
NUM_POINTS = 6
REFINE_SAMPLES = 200
U_CLIP_INIT = 50.0

# SciPy LP
try:
    from scipy.optimize import linprog
    SCIPY_OK = True
except ImportError:
    linprog = None
    SCIPY_OK = False

EDGE_PAIRS = [(i, j) for i in range(NUM_POINTS) for j in range(i + 1, NUM_POINTS)]
EDGE_I = np.array([i for (i, j) in EDGE_PAIRS], dtype=np.int64)
EDGE_J = np.array([j for (i, j) in EDGE_PAIRS], dtype=np.int64)
ALL_ROOTS = list(itertools.combinations(range(NUM_POINTS), 3))

def solve_lp_z_margin(C, s_edge, scale, u_clip_val):
    n_edges = C.shape[0]
    A_ub = np.zeros((n_edges, 5), dtype=np.float64)
    b_ub = np.zeros((n_edges,), dtype=np.float64)
    
    factor = -(scale * s_edge)[:, None]
    A_ub[:, :4] = factor * C
    A_ub[:, 4] = 1.0 
    
    c = np.array([0, 0, 0, 0, -1], dtype=np.float64)
    bounds = [(-u_clip_val, u_clip_val)] * 4 + [(0, None)]

    if SCIPY_OK:
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if not res.success: return False, None, -1.0
        x = res.x
        return (x[4] > 0), x[:4], float(x[4])
    return False, None, -1.0

def build_C_matrix(M1, t0):
    C = np.zeros((len(EDGE_PAIRS), 4), dtype=np.float64)
    for e, (i, j) in enumerate(EDGE_PAIRS):
        C[e, :] = t0[j] * M1[i, :] - t0[i] * M1[j, :]
    return C

def check_point_sage(candidate, u_opt, tx, ty, evaluator, RF):
    ts = [float(x) for x in candidate["ts"]]
    pair = (int(candidate["pair"][0]), int(candidate["pair"][1]))
    free_idx = [int(i) for i in candidate["free_idx"]]
    t0_free_vals = candidate["t0_free_vals"]
    
    lambdas = {i: vector(RF, [1, ts[i]]) for i in range(NUM_POINTS)}
    
    tildes_free = {}
    for k, idx in enumerate(free_idx):
        tildes_free[idx] = vector(RF, [float(t0_free_vals[k]), float(u_opt[k])])
        
    tildes = solve_conservation_pair(lambdas, tildes_free, pair, NUM_POINTS, RF)
    if tildes is None:
        return 0, -1.0, [], None, None

    xs = vector(RF, [1, tx])
    ys = vector(RF, [1, ty])
    
    def br(u, v): return u[0]*v[1] - u[1]*v[0]
    
    Cp = {i: br(lambdas[i], xs) * br(lambdas[i], ys) for i in range(NUM_POINTS)}
    
    pos_count = 0
    z_min = None
    z_vals = []
    
    for (i, j) in EDGE_PAIRS:
        ang = br(lambdas[i], lambdas[j])
        sq = br(tildes[j], tildes[i])
        val = (sq / ang) * Cp[i] * Cp[j]
        zf = float(val)
        z_vals.append(zf)
        if zf > 0: pos_count += 1
        z_min = zf if z_min is None else min(z_min, zf)
        
    return pos_count, z_min, z_vals, tildes, lambdas

def check_active_charts(evaluator, ts, tildes, tx, ty):
    xs = [1.0, float(tx)]
    ys = [1.0, float(ty)]
    # Convert tildes to dict of list for evaluator
    tf_final = {}
    for i in range(NUM_POINTS):
        tf_final[i] = [float(tildes[i][0]), float(tildes[i][1])]
        
    active_charts = 0
    hit_roots = []
    for roots in ALL_ROOTS:
        val = evaluator.evaluate_chart_fiber_gauge(roots, ts, tf_final, xs, ys, require_polytope=True)
        if abs(val) > 1e-20:
            active_charts += 1
            hit_roots.append(roots)
    return active_charts, hit_roots

def main():
    print("Starting V11 Refinement Phase")
    
    cand_file = "RESULTS/v11_top_candidates.json"
    if not os.path.exists(cand_file):
        print(f"No candidate file found at {cand_file}")
        return
        
    with open(cand_file, "r") as f:
        candidates = json.load(f)
        
    print(f"Loaded {len(candidates)} candidates.")
    
    RF = RealField(200)
    evaluator = FiberGaugeEvaluator(NUM_POINTS)
    
    rng = np.random.default_rng(1234)
    
    for idx, cand in enumerate(candidates):
        print(f"\nProcessing Candidate {idx+1}/{len(candidates)} (Init POS={cand.get('pos_float')}, mz={cand.get('mz_float'):.2e})")
        
        # Reconstruct Base
        ts = np.array(cand["ts"], dtype=np.float64)
        pair = cand["pair"]
        free_idx = cand["free_idx"]
        a, b = pair
        
        # Interval bounds
        ia, ib = cand["ia"], cand["ib"]
        Lx, Ly = cand.get("Lx", 2.0), cand.get("Ly", 2.0)
        
        # Construct intervals list again to get bounds
        intervals = []
        for k in range(NUM_POINTS + 1):
            if k == 0:
                intervals.append((-float('inf'), float(ts[0])))
            elif k == NUM_POINTS:
                intervals.append((float(ts[NUM_POINTS-1]), float('inf')))
            else:
                intervals.append((float(ts[k-1]), float(ts[k])))
                
        # Define ranges for ia, ib
        def get_range(idx, L_inf, U_inf, Lx, Ly):
            l, u = intervals[idx]
            if idx == 0:
                # (-inf, ts[0]). Target range: [ts[0]-Lx, ts[0]]
                # Actually, in v11 we sampled uniform dist from ts[0].
                # Let's just use [ts[0]-Lx, ts[0]-0.01] for safety or re-sample similarly
                return u - Lx, u
            elif idx == NUM_POINTS:
                return l, l + Ly
            else:
                return l, u
                
        La, Ua = get_range(ia, -np.inf, np.inf, Lx, Ly)
        Lb, Ub = get_range(ib, -np.inf, np.inf, Lx, Ly)
        
        # Local Search in (tx, ty) space
        best_local = None
        
        # Precompute C_matrix part
        # t0 is fixed for this candidate? Yes, t0_free_vals are fixed.
        # But u is optimized. 
        # Wait, the candidate has t0_free_vals. We need to reconstruct t0 vector.
        
        t = ts
        ta, tb = t[a], t[b]
        denom = tb - ta
        t_free = t[free_idx]
        t0_free_vals = np.array(cand["t0_free_vals"], dtype=np.float64)
        
        S00 = np.sum(t0_free_vals)
        S10 = np.sum(t_free * t0_free_vals)
        b0 = (-S10 + ta * S00) / denom
        a0 = -S00 - b0
        t0 = np.zeros(NUM_POINTS, dtype=np.float64)
        t0[free_idx] = t0_free_vals
        t0[a] = a0
        t0[b] = b0
        
        # M1
        B = (ta - t_free) / denom
        A = -np.ones_like(B) - B
        M1 = np.zeros((NUM_POINTS, 4), dtype=np.float64)
        for k, fi in enumerate(free_idx): M1[fi, k] = 1.0
        M1[a, :] = A
        M1[b, :] = B
        
        C_mat = build_C_matrix(M1, t0)
        
        print(f"  Local search with {REFINE_SAMPLES} samples...")
        
        for _ in range(REFINE_SAMPLES):
            # Sample tx, ty in bounds
            # For finite bounds, stay away from edges
            def sample_safe(L, U):
                w = U - L
                if w < 1e-4: return (L+U)/2
                return rng.uniform(L + 0.1*w, U - 0.1*w)
                
            tx = sample_safe(La, Ua)
            ty = sample_safe(Lb, Ub)
            
            if tx >= ty - 1e-6: continue
            
            # Compute signs & scale
            Ci = (tx - t) * (ty - t)
            s = np.sign(Ci[EDGE_I] * Ci[EDGE_J]).astype(np.int8)
            s[s == 0] = 1
            
            ang_diff = np.abs(t[EDGE_J] - t[EDGE_I])
            ang_diff = np.where(ang_diff < 1e-12, 1e-12, ang_diff)
            ang_inv = 1.0 / ang_diff
            Cprod = np.abs(Ci[EDGE_I] * Ci[EDGE_J])
            scale = Cprod * ang_inv
            
            # Solve LP
            ok, u, mz = solve_lp_z_margin(C_mat, s, scale, U_CLIP_INIT)
            
            if ok and mz > 0:
                if best_local is None or mz > best_local["mz"]:
                    best_local = {
                        "mz": mz,
                        "u": u,
                        "tx": tx,
                        "ty": ty
                    }
        
        if best_local is None:
            print("  Failed to find valid point in local search.")
            continue
            
        print(f"  Best local mz: {best_local['mz']:.2e}")
        
        # Hi-Prec Check
        pos, zmin, zvals, tildes, _ = check_point_sage(cand, best_local["u"], best_local["tx"], best_local["ty"], evaluator, RF)
        
        print(f"  Hi-Prec Check: POS={pos}/15, zmin={zmin:.3e}")
        
        if pos == 15 and zmin > 1e-16:
            # Check Charts
            print("  Checking charts...")
            n_charts, roots = check_active_charts(evaluator, [float(x) for x in ts], tildes, best_local["tx"], best_local["ty"])
            print(f"  Active Charts: {n_charts}")
            
            if n_charts > 0:
                print("  ✨ GOLDEN FIND!")
                # Save
                out = dict(cand)
                out["u_free"] = [float(x) for x in best_local["u"]]
                out["tx"] = float(best_local["tx"])
                out["ty"] = float(best_local["ty"])
                out["zmin"] = zmin
                out["pos"] = 15
                out["active_charts"] = n_charts
                out["active_roots"] = [list(r) for r in roots]
                
                with open("RESULTS/golden_kinematics_v11.json", "w") as f:
                    json.dump(out, f, indent=2)
                
                print("  Saved to RESULTS/golden_kinematics_v11.json. Exiting.")
                return 
        else:
            print("  Still failing robustness check.")

    print("Refinement complete. No golden point found in this batch.")

if __name__ == "__main__":
    main()


