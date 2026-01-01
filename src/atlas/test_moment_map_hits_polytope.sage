import sys
import os
import random as rnd
from sage.all import *


sys.path.append(os.getcwd())

# Import FiberGaugeEvaluator from existing codebase
# Assuming it is in src/atlas/jacobian_fiber_gauge.sage
# We use load because it's a sage script that might not be a python module
load("src/atlas/jacobian_fiber_gauge.sage")

def test_moment_map_gating():
    print("Testing Moment Map Gating with Direct Positive Z Inputs...")
    
    n = 6
    evaluator = FiberGaugeEvaluator(n)
    roots = (0, 1, 2)
    
    # Ensure chart data is loaded
    evaluator.get_chart_data(roots)
    
    # Access MML from the evaluator to get edge count
    data = evaluator.get_chart_data(roots)
    mml = data['mml']
    num_edges = mml.num_edges
    
    print(f"Chart Roots: {roots}")
    print(f"Number of edges: {num_edges}")
    
    num_samples = 100
    hits = 0
    
    RF = RealField(200) # Use high precision as in main code
    
    for i in range(num_samples):
        # 1. Generate random positive z
        # Log-uniform might be better to cover scales, or just uniform [0.1, 10]
        z_vals = [RF(rnd.uniform(0.1, 10.0)) for _ in range(num_edges)]
        
        # 2. Compute X directly
        # We access the internal method _compute_X_H_generic which we know exists
        try:
            X, _ = evaluator._compute_X_H_generic(mml, z_vals, RF)
            
            # 3. Check Polytope
            in_poly = evaluator.is_in_polytope(roots, X, RF)
            
            if in_poly:
                hits += 1
                
        except Exception as e:
            print(f"Error on sample {i}: {e}")
            continue
            
    print("-" * 50)
    print(f"Results for {num_samples} random positive Z samples:")
    print(f"Inside Polytope: {hits} ({float(hits)/num_samples*100:.1f}%)")
    
    if hits == 0:
        print("\nCRITICAL FAILURE: No positive Z landed in the polytope.")
        print("Possible causes:")
        print("1. Polytope inequalities are wrong (sign flip, wrong definition).")
        print("2. Moment map definition is inconsistent with polytope.")
        print("3. Vertices of polytope do not match moment map image.")
    else:
        print("\nSUCCESS: Positive Z values map to the interior of the polytope.")
        print("This confirms the geometry definition is consistent.")

if __name__ == "__main__":
    test_moment_map_gating()

