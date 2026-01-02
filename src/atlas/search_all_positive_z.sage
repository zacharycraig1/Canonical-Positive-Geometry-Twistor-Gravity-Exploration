import sys
import os
import random as rnd
import json
import math
import itertools
import multiprocessing
import time
from sage.all import *

sys.path.append(os.getcwd())

# Move load() calls inside functions or handle imports carefully for multiprocessing
# Sage's load() might not work well with multiprocessing if called at module level for spawned processes
# But we can try. If it fails, we'll inline the logic or use sage.all properly.
# The previous files are pure sage/python so they should be fine if loaded in worker.
# To be safe, we will load them in the worker function.

def load_dependencies():
    # Only load if not already loaded (though load() usually re-executes)
    # We use a flag or check. But sage load is simple.
    # We will assume the worker process is fresh enough or just load.
    load("src/atlas/jacobian_fiber_gauge.sage")
    load("src/atlas/solve_conservation_pair.sage")

def softplus(x, scale=1.0):
    try:
        xs = x / scale
        if xs > 20: return x
        if xs < -20: return 0
        return scale * math.log(1 + math.exp(xs))
    except (OverflowError, ValueError):
        if x > 0: return x
        return 0

def compute_loss(z_vals, margin=1e-8, beta=0.1):
    abs_z = [abs(z) for z in z_vals]
    if not abs_z: return 0.0, 0
    med = sorted(abs_z)[len(abs_z)//2]
    scale = med + 1e-12
    loss = 0.0
    pos_count = 0
    for z in z_vals:
        if z > 0: pos_count += 1
        loss += softplus(-z, scale=scale)
        loss += beta * softplus(margin - z, scale=scale)
    return loss, pos_count

def get_z_values(lambdas, tildes, x_s, y_s, n):
    # This must match utilities
    def bracket_local(l1, l2):
         return l1[0]*l2[1] - l1[1]*l2[0]
         
    C = {}
    for i in range(n):
        C[i] = bracket_local(lambdas[i], x_s) * bracket_local(lambdas[i], y_s)
        
    z_vals = []
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            ang = bracket_local(lambdas[i], lambdas[j])
            sq = bracket_local(tildes[j], tildes[i]) 
            if abs(ang) < 1e-20: ang = 1e-20
            val = (sq / ang) * C[i] * C[j]
            z_vals.append(float(val))
            edges.append((i,j))
    return z_vals, edges

def get_interval_candidates(ts):
    n = len(ts)
    regions = []
    regions.append(ts[0] - 2.0)
    for i in range(n-1):
        regions.append( (ts[i] + ts[i+1]) / 2.0 )
    regions.append(ts[-1] + 2.0)
    candidates = []
    for i in range(len(regions)):
        for j in range(len(regions)): 
            candidates.append( (regions[i], regions[j]) )
    return candidates

def worker_search(args):
    """
    Worker function for parallel search.
    Args: (seed, num_restarts, steps_per_restart)
    """
    worker_id, num_restarts, steps_per_restart = args
    
    # Re-seed RNG
    rnd.seed(os.getpid() + time.time())
    set_random_seed(os.getpid())
    
    # Load Dependencies
    load("src/atlas/jacobian_fiber_gauge.sage")
    load("src/atlas/solve_conservation_pair.sage")
    
    n = 6
    RF = RealField(200)
    evaluator = FiberGaugeEvaluator(n)
    
    decay = math.exp(math.log(1e-4) / steps_per_restart)
    preferred_pairs = [(0,1), (4,5), (0,5), (2,3)]
    
    best_pos_local = 0
    
    for restart in range(num_restarts):
        # Init
        if rnd.random() < 0.7:
            pair = rnd.choice(preferred_pairs)
        else:
            all_pairs = list(itertools.combinations(range(n), 2))
            pair = rnd.choice(all_pairs)
            
        ts_raw = sorted([rnd.uniform(1, 10) for _ in range(n)])
        lambdas = {i: vector(RF, [1, ts_raw[i]]) for i in range(n)}
        
        free_indices = [i for i in range(n) if i not in pair]
        tildes_free = {}
        for idx in free_indices:
            u = rnd.uniform(-5, 5)
            tildes_free[idx] = vector(RF, [1, u])
            
        # Brute Force tx, ty
        candidates = get_interval_candidates(ts_raw)
        best_init_tx, best_init_ty = candidates[0]
        best_init_score = -float('inf')
        
        tildes_init = solve_conservation_pair(lambdas, tildes_free, pair, n, RF)
        
        if tildes_init is not None:
            for (cx, cy) in candidates:
                x_c = vector(RF, [1, cx])
                y_c = vector(RF, [1, cy])
                z, _ = get_z_values(lambdas, tildes_init, x_c, y_c, n)
                pc = sum(1 for v in z if v > 0)
                if pc > best_init_score:
                    best_init_score = pc
                    best_init_tx, best_init_ty = cx, cy
                    
        tx, ty = best_init_tx, best_init_ty
        
        current_tildes_free = tildes_free.copy()
        current_tx, current_ty = tx, ty
        
        tildes = solve_conservation_pair(lambdas, current_tildes_free, pair, n, RF)
        if tildes is None: continue
        
        x_s = vector(RF, [1, current_tx])
        y_s = vector(RF, [1, current_ty])
        z_vals, _ = get_z_values(lambdas, tildes, x_s, y_s, n)
        current_loss, current_pos = compute_loss(z_vals)
        
        temp = 1.0
        
        for step in range(steps_per_restart):
            # Update best local
            if current_pos > best_pos_local:
                best_pos_local = current_pos
                # Optional: Log occasionally? No, difficult from worker.
            
            # Check Golden
            if current_pos == 15:
                # Check Active Charts
                active_charts = 0
                all_roots = list(itertools.combinations(range(n), 3))
                hit_roots = []
                
                # Check ALL roots? Or lazy?
                # We need to find ONE to be "Golden", but having more is better.
                # Let's check all to report full stats.
                
                ts_final = ts_raw
                tf_final = {k: [float(v[0]), float(v[1])] for k,v in current_tildes_free.items()}
                xs_final = [1.0, float(current_tx)]
                ys_final = [1.0, float(current_ty)]
                
                valid_point = False
                
                for roots in all_roots:
                    val = evaluator.evaluate_chart_fiber_gauge(roots, ts_final, tf_final, xs_final, ys_final, require_polytope=True)
                    if abs(val) > 1e-20:
                        active_charts += 1
                        hit_roots.append(roots)
                        valid_point = True
                
                if valid_point:
                    result_data = {
                        "ts": ts_raw,
                        "tildes_free": tf_final,
                        "x_s": xs_final,
                        "y_s": ys_final,
                        "pair": (int(pair[0]), int(pair[1])),
                        "z_values": z_vals,
                        "active_charts": active_charts,
                        "active_roots": [list(r) for r in hit_roots]
                    }
                    return result_data # Success!
                    
            # Annealing Step
            next_tildes_free = {}
            for idx in free_indices:
                delta = vector(RF, [rnd.gauss(0, 0.1)*temp, rnd.gauss(0, 0.1)*temp])
                next_tildes_free[idx] = current_tildes_free[idx] + delta
            
            next_tx = current_tx + rnd.gauss(0, 0.5) * temp
            next_ty = current_ty + rnd.gauss(0, 0.5) * temp
            
            next_tildes = solve_conservation_pair(lambdas, next_tildes_free, pair, n, RF)
            if next_tildes is None: continue
            
            next_x_s = vector(RF, [1, next_tx])
            next_y_s = vector(RF, [1, next_ty])
            
            next_z, _ = get_z_values(lambdas, next_tildes, next_x_s, next_y_s, n)
            next_loss, next_pos = compute_loss(next_z)
            
            delta_loss = next_loss - current_loss
            if delta_loss < 0 or rnd.random() < math.exp(-delta_loss / temp):
                current_tildes_free = next_tildes_free
                current_tx, current_ty = next_tx, next_ty
                current_loss = next_loss
                current_pos = next_pos
                z_vals = next_z
            
            temp *= decay

    return best_pos_local # Failure (return best pos as int)

def search_parallel():
    print("Starting PARALLEL Optimized Kinematic Search...")
    
    num_cores = 16
    total_restarts = 1000 # Increase total search volume
    restarts_per_worker = total_restarts // num_cores
    steps_per_restart = 2000
    
    pool_args = [(i, restarts_per_worker, steps_per_restart) for i in range(num_cores)]
    
    print(f"Workers: {num_cores}, Restarts/Worker: {restarts_per_worker}, Total: {num_cores * restarts_per_worker}")
    
    # Use context manager for pool
    with multiprocessing.Pool(processes=num_cores) as pool:
        # We want to get results as they complete? 
        # Or just use map and check output.
        # imap_unordered is good.
        
        results = pool.imap_unordered(worker_search, pool_args)
        
        best_overall = 0
        
        for res in results:
            if isinstance(res, dict):
                # Found Golden!
                print("\n\n✨ GOLDEN CONFIGURATION FOUND! ✨")
                print(f"Active Charts: {res['active_charts']}")
                print(f"Roots: {res['active_roots'][:3]}...")
                
                os.makedirs("RESULTS", exist_ok=True)
                with open("RESULTS/golden_kinematics.json", "w") as f:
                    json.dump(res, f, indent=2)
                    
                pool.terminate()
                return
            else:
                # It returned an int (best pos)
                if res > best_overall:
                    best_overall = res
                print(f"Worker finished. Best local pos: {res}/15. Global Best: {best_overall}/15")
                
    print("\nSearch Exhausted. No golden configuration found.")

if __name__ == "__main__":
    search_parallel()
