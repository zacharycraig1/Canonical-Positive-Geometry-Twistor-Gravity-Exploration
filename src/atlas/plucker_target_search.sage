import sys
import os
import random as rnd
import json
import math
import numpy as np
import time
from scipy.optimize import minimize as scipy_minimize
from sage.all import *

sys.path.append(os.getcwd())

load("src/atlas/jacobian_fiber_gauge.sage")
load("src/atlas/solve_conservation_pair.sage")


def get_bracket(u, v):
    # u, v are 2-vectors (numpy or sage)
    return u[0]*v[1] - u[1]*v[0]

def softplus(x, scale=1.0):
    try:
        xs = x / scale
        if xs > 20: return x
        if xs < -20: return 0
        return scale * math.log(1 + math.exp(xs))
    except (OverflowError, ValueError):
        if x > 0: return x
        return 0

def search_plucker_softplus():
    print("Starting Plucker Softplus Search (Relaxed Conservation)...")
    
    n = 6
    RF = RealField(200)
    evaluator = FiberGaugeEvaluator(n)
    
    # ------------------------------------------------------------------
    # Parameters
    # ------------------------------------------------------------------
    num_attempts = 1000
    cons_weight = 1e5
    
    best_pos_count_global = 0
    
    for attempt in range(num_attempts):
        # 1. Initialize random Geometry
        ts_raw = sorted([rnd.uniform(1, 10) for _ in range(n)])
        L = np.array([[1.0, t] for t in ts_raw]) # n x 2
        
        tx = ts_raw[0] - rnd.uniform(1, 3)
        ty = ts_raw[-1] + rnd.uniform(1, 3)
        X = np.array([1.0, tx])
        Y = np.array([1.0, ty])
        
        # Precompute constants
        brackets_L = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                brackets_L[i,j] = get_bracket(L[i], L[j])
                
        C = np.zeros(n)
        for i in range(n):
            C[i] = get_bracket(L[i], X) * get_bracket(L[i], Y)
            
        # 2. Initial Guess for T
        # Use a "valid" point from solving conservation pair?
        # Or just random?
        # A valid point is a better start.
        # Pick random pair
        pair = (4, 5) # Default
        
        # Free tildes random
        tildes_free_sage = {}
        for i in range(n):
            if i in pair: continue
            u = rnd.uniform(-5, 5)
            tildes_free_sage[i] = vector(RF, [1, u])
            
        lambdas_sage = {i: vector(RF, [1, ts_raw[i]]) for i in range(n)}
        tildes_sage_init = solve_conservation_pair(lambdas_sage, tildes_free_sage, pair, n, RF)
        
        if tildes_sage_init is None:
            # Fallback to random initialization
            x0 = np.random.normal(0, 1, 2*n)
        else:
            # Flatten sage vectors to x0
            x0_list = []
            for i in range(n):
                x0_list.append(float(tildes_sage_init[i][0]))
                x0_list.append(float(tildes_sage_init[i][1]))
            x0 = np.array(x0_list)
            
        # 3. Objective Function
        def objective(x_params):
            T = x_params.reshape((n, 2))
            
            loss = 0.0
            
            # Conservation Penalty
            P_sum = np.dot(L.T, T)
            cons_loss = np.sum(P_sum**2)
            loss += cons_weight * cons_loss
            
            # Positivity Loss
            z_vals = []
            for i in range(n):
                for j in range(i+1, n):
                    sq = get_bracket(T[j], T[i]) # [j i]
                    ang = brackets_L[i,j]
                    # Avoid zero division in ang (though L is fixed and ordered so ang > 0 usually)
                    if abs(ang) < 1e-12: ang = 1e-12
                    
                    val = (sq / ang) * C[i] * C[j]
                    z_vals.append(val)
            
            # Compute softplus loss
            # Scale?
            abs_z = [abs(z) for z in z_vals]
            scale = 1.0
            if abs_z:
                 scale = sorted(abs_z)[len(abs_z)//2] + 1e-8
            
            margin = 1e-6
            beta = 0.5
            
            sp_loss = 0.0
            for z in z_vals:
                # Primary: Softplus(-z)
                sp_loss += softplus(-z, scale=scale)
                # Margin: Softplus(margin - z)
                sp_loss += beta * softplus(margin - z, scale=scale)
                
            loss += sp_loss
            
            return loss

        # Run Optimization
        res = scipy_minimize(objective, x0, method='L-BFGS-B', tol=1e-6, options={'maxiter': 2000})
        
        # 4. Check
        T_final = res.x.reshape((n, 2))
        P_sum = np.dot(L.T, T_final)
        cons_err = np.sum(P_sum**2)
        
        pos_count = 0
        z_vals = []
        for i in range(n):
            for j in range(i+1, n):
                sq = get_bracket(T_final[j], T_final[i])
                ang = brackets_L[i,j]
                val = (sq / ang) * C[i] * C[j]
                z_vals.append(val)
                if val > 0: pos_count += 1
                
        is_conserved = cons_err < 1e-8 # slightly loose for optimization result
        
        if attempt % 50 == 0:
            print(f"Attempt {attempt}: Pos={pos_count}/15, ConsErr={cons_err:.2e}, Loss={res.fun:.2e}")
            
        if pos_count > best_pos_count_global:
            best_pos_count_global = pos_count
            print(f"  New Best Positivity: {pos_count}/15 (ConsErr={cons_err:.2e})")
            
        if pos_count == 15 and is_conserved:
            print(f"  Candidate Found! 15/15, ConsErr={cons_err:.2e}")
            
            # Reconstruct exact conservation
            # We pick a pair and solve for it using the optimized free vars.
            # This ensures we are exactly on the manifold for the final check.
            
            # Pick pair that has least "disturbance" or just standard?
            pair = (4, 5) 
            free_indices = [0, 1, 2, 3]
            
            tildes_free_recon = {}
            for i in free_indices:
                tildes_free_recon[i] = vector(RF, [T_final[i,0], T_final[i,1]])
                
            lambdas_sage = {i: vector(RF, [1, ts_raw[i]]) for i in range(n)}
            tildes_recon = solve_conservation_pair(lambdas_sage, tildes_free_recon, pair, n, RF)
            
            if tildes_recon is None:
                print("  Reconstruction failed.")
                continue
                
            # Check recon positivity
            recon_pos_count = 0
            recon_z_vals = []
            
            x_sage = vector(RF, [1, tx])
            y_sage = vector(RF, [1, ty])
            C_sage = {}
            for i in range(n): C_sage[i] = bracket(lambdas_sage[i], x_sage) * bracket(lambdas_sage[i], y_sage)
            
            for i in range(n):
                for j in range(i+1, n):
                    ang = bracket(lambdas_sage[i], lambdas_sage[j])
                    sq = bracket(tildes_recon[j], tildes_recon[i])
                    val = (sq / ang) * C_sage[i] * C_sage[j]
                    recon_z_vals.append(float(val))
                    if val > 0: recon_pos_count += 1
            
            print(f"  Reconstructed Positivity: {recon_pos_count}/15")
            
            if recon_pos_count == 15:
                # Check Charts
                print("  Checking Active Charts on Reconstructed Point...")
                active_charts = 0
                all_roots = list(itertools.combinations(range(n), 3))
                hit_roots = []
                
                ts_final = ts_raw
                tf_final = {k: [float(v[0]), float(v[1])] for k,v in tildes_free_recon.items()}
                xs_final = [1.0, float(tx)]
                ys_final = [1.0, float(ty)]
                
                for roots in all_roots:
                    val = evaluator.evaluate_chart_fiber_gauge(roots, ts_final, tf_final, xs_final, ys_final, require_polytope=True)
                    if abs(val) > 1e-20:
                        active_charts += 1
                        hit_roots.append(roots)
                
                if active_charts > 0:
                    print(f"  ✨ GOLDEN FIND! Active Charts: {active_charts}")
                    
                    result_data = {
                        "ts": ts_final,
                        "tildes_free": tf_final,
                        "x_s": xs_final,
                        "y_s": ys_final,
                        "pair": (int(pair[0]), int(pair[1])),
                        "z_values": recon_z_vals,
                        "active_charts": active_charts,
                        "active_roots": [list(r) for r in hit_roots]
                    }
                    
                    with open("RESULTS/golden_kinematics.json", "w") as f:
                        json.dump(result_data, f, indent=2)
                    return
                else:
                     print("  Recon 15/15 -> No Active Charts.")
            else:
                print("  Lost positivity after reconstruction.")
                
    print("Search exhausted.")

if __name__ == "__main__":
    search_plucker_softplus()
