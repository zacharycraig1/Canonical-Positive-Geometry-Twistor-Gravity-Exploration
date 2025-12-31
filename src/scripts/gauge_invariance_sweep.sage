import sys
import os
import random as py_random
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.posgeom.forest_polytope import get_forest_polynomial
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import ang_bracket

def gauge_invariance_sweep():
    print("Checking Gauge Invariance (Reference Independence) of the Geometric Formula...")
    
    # Configuration
    n = 6
    roots = [0, 1, 2]
    num_gauge_trials = 10
    
    # 1. Fix Physical Kinematics (Momentum Twistors -> Spinors)
    print(f"Generating fixed kinematics for n={n}...")
    lambdas, tildes = sample_spinors_from_twistor(n=n, seed=999)
    
    # 2. Pre-compute Forest Polynomial
    print("Building Forest Polynomial...")
    F_poly = get_forest_polynomial(n, roots)
    
    # Store results
    values = []
    
    print(f"Sweeping over {num_gauge_trials} random reference spinors (x, y)...")
    
    for i in range(num_gauge_trials):
        # Generate random reference spinors x, y
        # We can simulate them as random complex vectors
        # Using simple integers/rationals for stability if possible, but floats are fine for invariance check
        # Let's use simple random integers to stay in QQ/ZZ if possible, or just floats
        
        # Random integers [-10, 10]
        x = [py_random.randint(-10, 10), py_random.randint(-10, 10)]
        y = [py_random.randint(-10, 10), py_random.randint(-10, 10)]
        
        # Avoid singular refs (collinear with each other or any lambda)
        # Simple check: x != k*y
        if x[0]*y[1] - x[1]*y[0] == 0: 
            continue
            
        # 3. Compute the Geometric Formula Value
        # M = (-1)^(n-1) * <01>^8 * F(z) / (Norm * prod(Ck^2))
        
        # a) Map to z_{ij}
        try:
            z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
        except ValueError:
            # Skip singular configs
            continue
            
        # b) Evaluate F(z)
        R_ring = F_poly.parent()
        eval_dict = {R_ring(k): v for k, v in z_map.items() if isinstance(k, str) and k.startswith('z_')}
        F_val = F_poly.subs(eval_dict).constant_coefficient()
        
        # c) Compute Denominator Factors
        # C_k for k not in roots
        def get_C(k, ref_x, ref_y):
            return (lambdas[k][0]*ref_x[1] - lambdas[k][1]*ref_x[0]) * \
                   (lambdas[k][0]*ref_y[1] - lambdas[k][1]*ref_y[0])
                   
        prod_C_sq = 1
        for k in range(n):
            if k not in roots:
                prod_C_sq *= get_C(k, x, y)**2
                
        # Norm factor (independent of x,y, but part of the formula)
        r1, r2, r3 = roots
        norm_factor = (ang_bracket(lambdas[r1], lambdas[r2]) * 
                       ang_bracket(lambdas[r2], lambdas[r3]) * 
                       ang_bracket(lambdas[r3], lambdas[r1]))**2
                       
        # Helicity factor <01>^8 (independent of x,y)
        h_factor = ang_bracket(lambdas[0], lambdas[1])**8
        
        sign = (-1)**(n-1)
        
        # Result
        if prod_C_sq == 0:
            continue
            
        M_geom = sign * h_factor * F_val / (norm_factor * prod_C_sq)
        values.append(M_geom)
        
        print(f"Gauge {i}: x={x}, y={y} -> M = {float(M_geom):.4e} + {float(M_geom.imag()):.4e}j")

    # 4. Check Invariance
    if not values:
        print("No valid trials.")
        return

    first_val = values[0]
    max_diff = 0
    
    for v in values[1:]:
        diff = abs(v - first_val)
        if diff > max_diff:
            max_diff = diff
            
    print(f"\nResults Summary:")
    print(f"Base Value: {float(first_val):.6e}")
    print(f"Max Relative Difference: {float(max_diff/abs(first_val)):.2e}")
    
    if max_diff / abs(first_val) < 1e-10:
        print("SUCCESS: Geometric formula is Gauge Invariant!")
    else:
        print("FAILURE: Formula depends on reference spinors.")

if __name__ == "__main__":
    gauge_invariance_sweep()

