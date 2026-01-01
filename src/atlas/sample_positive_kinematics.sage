import sys
import os
import json
import random as rnd
from sage.all import *

sys.path.append(os.getcwd())

load("src/atlas/jacobian_fiber_gauge.sage")
load("src/atlas/check_residues_gate_b.sage")

def sample_positive_kinematics():
    print("Sampling Positive Kinematics...")
    
    n = 6
    RF = RealField(200)
    evaluator = FiberGaugeEvaluator(n)
    
    found_any = False
    
    for attempt in range(500):
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
        
        # Check z signs
        key = evaluator.get_sign_pattern(ts, ts_tilde_free, x_s, y_s, RF)
        if key is None: continue
        
        # 15 bits should be 1
        target_mask = (1 << 15) - 1
        
        if key == target_mask:
            print(f"Success! Found all-positive z point at attempt {attempt}")
            found_any = True
            
            # Check polytope strictness
            # We check if ANY chart is active with require_polytope=True
            active_count = 0
            
            # Check a subset of charts to be fast
            check_roots = [(0,1,2), (0,1,3), (1,2,3), (3,4,5)] 
            
            for roots in check_roots:
                 # Need to re-instantiate evaluator or just call method?
                 # Method is stateless regarding roots (except cache).
                 val = evaluator.evaluate_chart_fiber_gauge(roots, ts, ts_tilde_free, x_s, y_s, prec=200, require_polytope=True)
                 if abs(val) > 1e-20:
                     active_count += 1
                     print(f"  Chart {roots} is ACTIVE (val={val:.2e})")
            
            if active_count > 0:
                print(f"  Point lies inside at least {active_count} checked charts.")
            else:
                print("  Point has positive z but might not be in checked charts (or evaluator gating failed).")
            
            # Save sample
            sample = {
                "ts": [float(t) for t in ts],
                "ts_tilde_free": {str(k): [float(x) for x in v] for k,v in ts_tilde_free.items()},
                "x": [float(x) for x in x_s],
                "y": [float(x) for x in y_s]
            }
            with open("RESULTS/positive_kinematics_sample.json", "w") as f:
                json.dump(sample, f, indent=2)
            break
            
    if not found_any:
        print("Failed to find all-positive sample after 500 attempts.")

if __name__ == "__main__":
    sample_positive_kinematics()


