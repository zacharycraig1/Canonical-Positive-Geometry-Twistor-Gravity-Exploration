import sys
import os
import random as rnd
from sage.all import *
import itertools

sys.path.append(os.getcwd())

load("src/atlas/jacobian_fiber_gauge.sage")

def bracket(l1, l2):
    return l1[0]*l2[1] - l1[1]*l2[0]

def run_diagnostics(num_samples=200):
    print(f"Running Z-Sampling Diagnostics ({num_samples} samples)...")
    
    n = 6
    evaluator = FiberGaugeEvaluator(n)
    RF = RealField(200)
    
    # Track statistics
    stats = {
        "bitcounts": {}, # count -> frequency
        "active_charts_hist": {}, # count -> frequency
        "fully_positive": 0,
        "valid_samples": 0
    }
    
    # Pre-select some roots to check for active charts
    all_roots = list(itertools.combinations(range(n), 3))
    
    for k in range(num_samples):
        # Sampling Heuristic: Ordered Tildes (Strategy 3A)
        # 1. Lambdas
        ts = sorted([rnd.uniform(0, 10) for _ in range(n)])
        lambdas = {i: vector(RF, [1, ts[i]]) for i in range(n)}
        
        # 2. Tildes (Ordered on 0..3)
        # u_i strictly decreasing to get [ji] > 0
        us = sorted([rnd.uniform(0, 10) for _ in range(n-2)], reverse=True)
        ts_tilde_free = {i: vector(RF, [1, us[i]]) for i in range(n-2)}
        
        # Aux spinors
        x_s = vector(RF, [1, -2.0])
        y_s = vector(RF, [1, 12.0])
        
        # 3. Solve Conservation
        tildes = solve_conservation_generic(lambdas, ts_tilde_free, n, RF)
        
        if tildes is None:
            continue
            
        # Compute C
        C = {}
        for i in range(n): C[i] = bracket(lambdas[i], x_s) * bracket(lambdas[i], y_s)

        stats["valid_samples"] += 1
        
        z_vals = []
        pos_count = 0
        min_z = float('inf')
        max_z = float('-inf')
        
        # Check all 15 edges
        for i in range(n):
            for j in range(i+1, n):
                ang = bracket(lambdas[i], lambdas[j])
                # Note: [ji] convention fixed implies we use bracket(tildes[j], tildes[i])
                sq = bracket(tildes[j], tildes[i]) 
                if abs(ang) < RF(1e-25): ang = RF(1e-25)
                val = (sq / ang) * C[i] * C[j]
                
                z_vals.append(val)
                if val > 0: pos_count += 1
                else:
                    if k < 5: print(f"  Neg Z at {i},{j}: {val:.2e} (ang={ang:.2e}, sq={sq:.2e})")
                min_z = min(min_z, float(val))
                max_z = max(max_z, float(val))
                
        stats["bitcounts"][pos_count] = stats["bitcounts"].get(pos_count, 0) + 1
        if pos_count == 15:
            stats["fully_positive"] += 1
            
        # 2. Check Active Charts
        active_count = 0
        # Only check charts if we have reasonable positivity? 
        # No, check always to see correlation.
        
        for roots in all_roots:
            roots_t = tuple(sorted(list(roots)))
            # This uses the internal z calculation of the evaluator, which should now match the fix
            val = evaluator.evaluate_chart_fiber_gauge(roots_t, ts, ts_tilde_free, x_s, y_s, prec=200, require_polytope=True)
            if abs(val) > 0: # Non-zero means active
                active_count += 1
                
        stats["active_charts_hist"][active_count] = stats["active_charts_hist"].get(active_count, 0) + 1
        
        if k % 20 == 0:
            print(f"Sample {k}: Bitcount={pos_count}/15, Active Charts={active_count}, Range=[{min_z:.2e}, {max_z:.2e}]")

    print("\n--- Diagnostics Summary ---")
    print(f"Total Valid Samples: {stats['valid_samples']}/{num_samples}")
    print(f"Fully Positive Z Samples: {stats['fully_positive']}")
    
    print("\nBitcount Distribution:")
    for c in sorted(stats["bitcounts"].keys()):
        print(f"  {c}/15: {stats['bitcounts'][c]}")
        
    print("\nActive Charts Distribution:")
    for c in sorted(stats["active_charts_hist"].keys()):
        print(f"  {c} charts: {stats['active_charts_hist'][c]}")
        
    has_active = sum(k*v for k,v in stats['active_charts_hist'].items() if k > 0) > 0
    if stats['fully_positive'] > 0 or has_active:
        print("\n[SUCCESS] Found samples with positive Z and/or active charts.")
    else:
        print("\n[WARNING] Failed to find ideal samples. Check sign convention or sampling heuristic.")

if __name__ == "__main__":
    run_diagnostics()

