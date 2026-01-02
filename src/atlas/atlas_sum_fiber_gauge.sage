import sys
import os
import json
import random as rnd
from sage.all import *

sys.path.append(os.getcwd())

load("src/atlas/jacobian_fiber_gauge.sage")
load("src/atlas/check_residues_gate_b.sage")

def run_atlas_sum_verification():
    print("Running Full Atlas Sum Verification (Analytic Continuation / Chamber Restricted)...")
    
    # Load Coeffs
    coeff_path = "RESULTS/atlas_coeffs_exact.json"
    if not os.path.exists(coeff_path):
        print("Coefficients not found. Run solve_boundary_coefficients.sage first.")
        return
        
    with open(coeff_path, 'r') as f:
        data = json.load(f)
        
    coeffs = {}
    for item in data:
        roots = tuple(sorted(item['roots']))
        coeffs[roots] = item['coeff']
        
    # Load Best Chamber
    best_chamber = None
    if os.path.exists("RESULTS/best_chamber.txt"):
        with open("RESULTS/best_chamber.txt", "r") as f:
            best_chamber = int(f.read().strip())
        print(f"Restricting verification to chamber: {best_chamber}")
    else:
        print("No best chamber found. Verification might fail if mixing chambers.")
        
    n = 6
    evaluator = FiberGaugeEvaluator(n)
    all_roots = list(coeffs.keys())
    
    RF = RealField(200)
    
    print("\n--- Ratio Verification on New Points ---")
    
    num_points = 20
    valid_points = 0
    attempts = 0
    ratios = []
    
    while valid_points < num_points and attempts < num_points * 20:
        attempts += 1
        
        # Biased Sampling (Same as solver)
        ts = sorted([rnd.uniform(1, 9) for _ in range(n)])
        tx = ts[0] - rnd.uniform(1, 3)
        ty = ts[-1] + rnd.uniform(1, 3)
        x_s = vector(RF, [1, tx])
        y_s = vector(RF, [1, ty])
        
        ts_tilde_free = {}
        for i in range(n-2):
            noise = vector(RF, [rnd.gauss(0, 0.1), rnd.gauss(0, 0.1)])
            base = vector(RF, [1, -ts[i]])
            ts_tilde_free[i] = base + noise

        lambdas = {i: vector(RF, [1, ts[i]]) for i in range(n)}
        tildes = solve_conservation_generic(lambdas, ts_tilde_free, n, RF)
        if tildes is None: continue
        
        # Check Chamber
        if best_chamber is not None:
            key = evaluator.get_sign_pattern(ts, ts_tilde_free, x_s, y_s, RF)
            if key != best_chamber: continue
        
        # Atlas Sum
        atlas_sum = RF(0)
        active_charts = 0
        
        for roots in all_roots:
            c = coeffs.get(roots, 0)
            if c == 0: continue
            
            # Analytic continuation: require_polytope=False
            val = evaluator.evaluate_chart_fiber_gauge(roots, ts, ts_tilde_free, x_s, y_s, prec=200, require_polytope=False)
            if abs(val) > 1e-25:
                atlas_sum += c * val
                active_charts += 1
                
        # Target
        m_mhv = compute_M_MHV(n, lambdas, tildes, x_s, y_s, RF)
        
        if abs(m_mhv) > 1e-20:
            ratio = atlas_sum / m_mhv
            valid_points += 1
            ratios.append(float(ratio))
            print(f"Pt {valid_points}: Sum={float(atlas_sum):.2e}, MHV={float(m_mhv):.2e}, Ratio={float(ratio):.4f}, Active={active_charts}")
        
    # Stats
    if ratios:
        avg = sum(ratios) / len(ratios)
        dev = sum([(r - avg)**2 for r in ratios]) / len(ratios)
        print(f"\nAverage Ratio: {avg:.6f}")
        print(f"Variance: {dev:.2e}")
        
        stats = {
            "avg_ratio": avg,
            "variance": dev,
            "num_points": len(ratios),
            "chamber": best_chamber
        }
        with open("RESULTS/atlas_sum_ratio_stats.json", "w") as f:
            json.dump(stats, f, indent=2)
    else:
        print("No valid points found in chamber.")

if __name__ == "__main__":
    run_atlas_sum_verification()
