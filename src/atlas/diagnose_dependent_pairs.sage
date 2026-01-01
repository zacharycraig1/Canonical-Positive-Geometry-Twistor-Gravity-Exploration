import sys
import os
import random as rnd
import json
import itertools
from collections import defaultdict
from sage.all import *

sys.path.append(os.getcwd())

load("src/atlas/jacobian_fiber_gauge.sage")
load("src/atlas/solve_conservation_pair.sage")

def diagnose_dependent_pairs():
    print("Diagnosing Dependent Pairs for Positivity and Active Charts...")
    
    n = 6
    RF = RealField(200)
    evaluator = FiberGaugeEvaluator(n)
    
    # Pre-warm charts (optional, but good for timing)
    # We will just check a few charts or all? Checking all 20 is fine.
    all_roots = list(itertools.combinations(range(n), 3))
    
    pairs = list(itertools.combinations(range(n), 2))
    
    results = []
    
    # We stop if we find a "Golden Configuration"
    found_golden = False
    
    for pair_idx, (a, b) in enumerate(pairs):
        print(f"\nEvaluating Pair {pair_idx+1}/{len(pairs)}: ({a}, {b})")
        
        stats = {
            "pair": (int(a), int(b)),
            "samples": 0,
            "perfect_positivity": 0, # 15/15
            "high_positivity": 0,    # >= 14/15
            "active_charts_hits": 0, # > 0 active charts
            "max_positivity": 0,
            "avg_positivity": 0.0,
            "best_sample": None
        }
        
        positivity_sum = 0
        num_samples = 200 # Check enough points
        
        for s in range(num_samples):
            # 1. Sample lambdas (ordered t)
            # t_0 < t_1 < ... < t_5
            ts = sorted([rnd.uniform(1, 10) for _ in range(n)])
            lambdas = {i: vector(RF, [1, ts[i]]) for i in range(n)}
            
            # 2. Sample free tildes (ordered slopes for free indices?)
            # The heuristic "ordered slopes" applies to indices 0..3 usually.
            # If we change dependent pair, we should probably stick to 
            # "ordered slopes for the free indices" or just "ordered slopes for all i except a,b"?
            # The user note says: "u_0 > u_1 > u_2 > u_3" was the old heuristic.
            # Let's try to enforce "ordered slopes" for the *free* indices in their natural order.
            
            free_indices = [i for i in range(n) if i not in (a, b)]
            
            # Generate sorted slopes for free indices
            # u_vals = sorted([rnd.uniform(-5, 5) for _ in range(len(free_indices))], reverse=True)
            # Actually, let's just try random slopes first, or sorted.
            # The hypothesis is that SOME ordering works.
            # Let's try strictly decreasing slopes for the free ones as a start, 
            # mapping them to the free indices in order.
            
            u_vals = sorted([rnd.uniform(-2, 2) for _ in range(len(free_indices))], reverse=True)
            
            tildes_free = {}
            for k, idx in enumerate(free_indices):
                # tilde = [1, u] approx
                # Add some noise to first component too?
                # The user note says: "free tilde_i approx (1, u_i)"
                u = u_vals[k]
                # small noise
                noise = vector(RF, [rnd.gauss(0, 0.01), rnd.gauss(0, 0.01)])
                base = vector(RF, [1, u])
                tildes_free[idx] = base + noise
                
            # 3. Solve
            tildes = solve_conservation_pair(lambdas, tildes_free, (a, b), n, RF)
            
            if tildes is None: continue
            
            # 4. Check Positivity
            # We need C_i. For positivity of z_{ij}, we usually need C_i C_j > 0.
            # We can just compute z_{ij} directly.
            # To compute z, we need x, y spinors.
            # Let's pick x, y such that C_i are all same sign?
            # Or just pick x, y outside the range of t.
            
            tx = ts[0] - rnd.uniform(1, 3)
            ty = ts[-1] + rnd.uniform(1, 3)
            x_s = vector(RF, [1, tx])
            y_s = vector(RF, [1, ty])
            
            # Compute C
            C = {}
            for i in range(n):
                C[i] = bracket(lambdas[i], x_s) * bracket(lambdas[i], y_s)
                
            # Compute z for all edges
            # We need edge list.
            # Just loop i < j
            pos_count = 0
            all_z = []
            
            # 15 edges
            for i in range(n):
                for j in range(i+1, n):
                    ang = bracket(lambdas[i], lambdas[j])
                    sq = bracket(tildes[j], tildes[i]) # Note index order for [j i]
                    # z ~ [j i] / <i j> * C_i C_j
                    
                    val = (sq / ang) * C[i] * C[j]
                    if val > 0: pos_count += 1
                    all_z.append(val)
                    
            stats["max_positivity"] = max(stats["max_positivity"], pos_count)
            positivity_sum += pos_count
            
            if pos_count == 15:
                stats["perfect_positivity"] += 1
            if pos_count >= 14:
                stats["high_positivity"] += 1
                
            # 5. Check Active Charts (only if positivity is decent? or always?)
            # Checking charts is expensive-ish. Let's check if pos_count >= 14
            # or if we are desperate.
            # The goal is 15/15 AND active charts.
            
            active_charts = 0
            if pos_count >= 14:
                # Construct Z vector for MML (needs canonical order)
                # The evaluator helper _compute_z_vals does this given lambdas/tildes
                # But we have local code.
                # Let's use evaluator to check charts for this point.
                
                # Check a subset of charts or all?
                # If ANY chart is active, that's a win.
                for roots in all_roots:
                    # We need to pass data to evaluator.
                    # evaluate_chart_fiber_gauge computes everything from scratch.
                    # We can use that.
                    val = evaluator.evaluate_chart_fiber_gauge(roots, ts, tildes_free, x_s, y_s, prec=200, require_polytope=True)
                    if abs(val) > 1e-20: # Non-zero means active
                        active_charts += 1
                        
                if active_charts > 0:
                    stats["active_charts_hits"] += 1
                    
                    if pos_count == 15:
                        print(f"  ✨ GOLDEN SAMPLE FOUND! Pair ({a},{b}), 15/15 Positive, {active_charts} Active Charts.")
                        stats["best_sample"] = {
                            "ts": [float(x) for x in ts],
                            "tildes_free": {k: [float(v[0]), float(v[1])] for k,v in tildes_free.items()},
                            "x_s": [float(x) for x in x_s],
                            "y_s": [float(x) for x in y_s],
                            "pair": (int(a), int(b))
                        }
                        found_golden = True
                        break
        
        stats["samples"] = num_samples
        stats["avg_positivity"] = float(positivity_sum) / num_samples if num_samples > 0 else 0.0
        results.append(stats)
        
        print(f"  Avg Positivity: {stats['avg_positivity']:.2f}")
        print(f"  Perfect (15/15): {stats['perfect_positivity']}")
        print(f"  Active Charts Hits: {stats['active_charts_hits']}")
        
        if found_golden:
            break
            
    # Sort results
    results.sort(key=lambda x: (x["active_charts_hits"], x["perfect_positivity"], x["avg_positivity"]), reverse=True)
    
    # Convert Sage Integers to Python ints for JSON
    for r in results:
        r["active_charts_hits"] = int(r["active_charts_hits"])
        r["perfect_positivity"] = int(r["perfect_positivity"])
        r["high_positivity"] = int(r["high_positivity"])
        r["samples"] = int(r["samples"])
        r["max_positivity"] = int(r["max_positivity"])
        # Pair is already cast
    
    # Save Report
    os.makedirs("RESULTS", exist_ok=True)
    with open("RESULTS/dependent_pairs_report.json", "w") as f:
        json.dump(results, f, indent=2)
        
    print("\n" + "="*50)
    print("DIAGNOSIS COMPLETE")
    print("Top 3 Pairs:")
    for r in results[:3]:
        print(f"  Pair {r['pair']}: ActiveHits={r['active_charts_hits']}, Perfect={r['perfect_positivity']}, Avg={r['avg_positivity']:.2f}")
        
    if found_golden:
        print("\nSUCCESS: Found configuration with 15/15 positivity and active charts.")
        # Save the golden sample specifically
        with open("RESULTS/golden_kinematics.json", "w") as f:
            json.dump(results[0]["best_sample"], f, indent=2)
    else:
        print("\nWARNING: No golden configuration found yet.")

if __name__ == "__main__":
    diagnose_dependent_pairs()

