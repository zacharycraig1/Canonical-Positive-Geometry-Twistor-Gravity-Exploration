import sys
import os
import itertools
import json
import random as rnd
from collections import defaultdict
from sage.all import *

sys.path.append(os.getcwd())

load("src/atlas/jacobian_fiber_gauge.sage")
load("src/atlas/check_residues_gate_b.sage")

def solve_boundary_coefficients():
    print("Solving Boundary Coefficients (Analytic Continuation / Chamber Split)...")
    
    n = 6
    all_roots = list(itertools.combinations(range(n), 3))
    evaluator = FiberGaugeEvaluator(n)
    RF = RealField(200)
    
    # Storage for points by chamber
    # key -> list of { 'row': [...], 'target': val }
    chambers = defaultdict(list)
    
    num_points_target = 300
    points_generated = 0
    valid_points = 0
    
    print(f"Generating up to {num_points_target} valid points...")
    
    while valid_points < num_points_target and points_generated < num_points_target * 10:
        points_generated += 1
        
        # Use Biased Sampling to find consistent chambers
        # 1. Ordered ts
        ts = sorted([rnd.uniform(1, 9) for _ in range(n)])
        
        # 2. x, y outside range
        tx = ts[0] - rnd.uniform(1, 3)
        ty = ts[-1] + rnd.uniform(1, 3)
        
        x_s = vector(RF, [1, tx])
        y_s = vector(RF, [1, ty])
        
        # 3. Tildes: approx [1, -ti]
        ts_tilde_free = {}
        for i in range(n-2):
            # Add small noise
            noise = vector(RF, [rnd.gauss(0, 0.1), rnd.gauss(0, 0.1)])
            base = vector(RF, [1, -ts[i]])
            ts_tilde_free[i] = base + noise
            
        lambdas = {i: vector(RF, [1, ts[i]]) for i in range(n)}
        tildes = solve_conservation_generic(lambdas, ts_tilde_free, n, RF)
        
        if tildes is None: continue
        
        # Get chamber key
        key = evaluator.get_sign_pattern(ts, ts_tilde_free, x_s, y_s, RF)
        if key is None: continue
        
        # Evaluate all charts without polytope restriction
        row = []
        possible_singularity = False
        
        for roots in all_roots:
            # We use require_polytope=False
            val = evaluator.evaluate_chart_fiber_gauge(tuple(sorted(list(roots))), ts, ts_tilde_free, x_s, y_s, prec=200, require_polytope=False)
            
            # Check for singularity (huge values)
            if abs(val) > 1e20: 
                possible_singularity = True
                break
            row.append(val)
            
        if possible_singularity: continue
        
        # Compute Target MHV
        target = compute_M_MHV(n, lambdas, tildes, x_s, y_s, RF)
        if abs(target) > 1e20: continue # Skip singular targets
        
        chambers[key].append({
            'row': row,
            'target': target
        })
        valid_points += 1
        
        if valid_points % 50 == 0:
            print(f"  Generated {valid_points} points. Found {len(chambers)} chambers so far.")
            
    # Analyze Chambers
    print(f"\nFound {len(chambers)} chambers.")
    for k, pts in chambers.items():
        print(f"  Chamber {k}: {len(pts)} points")
        
    # Solve for each chamber with enough points
    # We have 20 charts, so we need >= 20 points. Ideally > 30 for overconstrained.
    
    coeffs_results = []
    fit_quality = []
    
    os.makedirs("RESULTS", exist_ok=True)
    
    best_chamber = None
    max_pts = 0
    
    for key, pts in chambers.items():
        if len(pts) < 25:
            continue
            
        if len(pts) > max_pts:
            max_pts = len(pts)
            best_chamber = key
            
        print(f"\nSolving for Chamber {key} ({len(pts)} points)...")
        
        A_rows = [p['row'] for p in pts]
        b_vals = [p['target'] for p in pts]
        
        A = matrix(RF, A_rows)
        b = vector(RF, b_vals)
        
        try:
            # Least squares
            AT = A.transpose()
            # (AT * A) x = AT * b
            lhs = AT * A
            rhs = AT * b
            
            x = lhs.solve_right(rhs)
            
            # Calculate Residuals
            residual_sq_sum = 0
            for i in range(len(pts)):
                pred = A[i] * x
                actual = b[i]
                diff = pred - actual
                residual_sq_sum += diff * diff
                
            rmse = sqrt(residual_sq_sum / len(pts))
            print(f"  RMSE: {rmse:.2e}")
            
            # Store Coefficients
            chamber_coeffs = []
            for i, roots in enumerate(all_roots):
                val = float(x[i])
                rounded = int(val) # Truncate or round?
                if abs(val - round(val)) < 0.1:
                    rounded = int(round(val))
                
                err = abs(val - rounded)
                chamber_coeffs.append({
                    "roots": [int(r) for r in roots],
                    "coeff": int(rounded),
                    "fit_val": float(val),
                    "is_integer": bool(err < 0.1)
                })
                
            coeffs_results.append({
                "chamber_key": int(key),
                "num_points": int(len(pts)),
                "rmse": float(rmse),
                "coeffs": chamber_coeffs
            })
            
            fit_quality.append({
                "chamber_key": int(key),
                "num_points": int(len(pts)),
                "rmse": float(rmse)
            })
            
        except Exception as e:
            print(f"  Solver failed for chamber {key}: {e}")
            
    # Save Results
    with open("RESULTS/coeffs_by_chamber.json", "w") as f:
        json.dump(coeffs_results, f, indent=2)
        
    with open("RESULTS/fit_quality_by_chamber.json", "w") as f:
        json.dump(fit_quality, f, indent=2)
        
    print(f"\nSaved results to RESULTS/coeffs_by_chamber.json")
    
    # If we found a best chamber, save it as the 'main' coefficients for step 2
    if best_chamber is not None:
        best_res = next(r for r in coeffs_results if r['chamber_key'] == best_chamber)
        
        with open("RESULTS/atlas_coeffs_exact.json", "w") as f:
            json.dump(best_res['coeffs'], f, indent=2)
        
        print(f"Saved best chamber ({best_chamber}) coefficients to RESULTS/atlas_coeffs_exact.json")
        with open("RESULTS/best_chamber.txt", "w") as f:
            f.write(str(best_chamber))

if __name__ == "__main__":
    solve_boundary_coefficients()
