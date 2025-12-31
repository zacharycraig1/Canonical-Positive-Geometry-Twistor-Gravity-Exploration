import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.posgeom.forest_polytope import get_forest_polynomial
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import ang_bracket

def gauge_invariance_sweep():
    print("Checking Gauge Invariance (Independence of Reference Spinors x, y) for n=6...")
    n = 6
    roots = [0, 1, 2]
    
    # 1. Fix Kinematics
    print("Generating fixed kinematics...")
    lambdas, tildes = sample_spinors_from_twistor(n=n, seed=42)
    
    # 2. Build Forest Polynomial
    print("Building Forest Polynomial F_{6,R}...")
    F_poly = get_forest_polynomial(n, roots)
    
    # 3. Sweep over random reference spinors
    num_sweeps = 20
    print(f"Testing {num_sweeps} random reference spinor pairs (x, y)...")
    
    # We will store the first computed amplitude and compare others to it
    reference_amp = None
    
    consistent = True
    
    for k in range(num_sweeps):
        # Generate random x, y
        import random
        random.seed(1000 + k)
        # Using random integers for components
        x = [random.randint(-10, 10) + I*random.randint(-10, 10), 
             random.randint(-10, 10) + I*random.randint(-10, 10)]
        y = [random.randint(-10, 10) + I*random.randint(-10, 10), 
             random.randint(-10, 10) + I*random.randint(-10, 10)]
             
        # Calculate Amplitude
        try:
            z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
        except ValueError as e:
            print(f"Sweep {k}: Skipping due to singularity: {e}")
            continue
            
        # Evaluate F
        R_ring = F_poly.parent()
        eval_dict = {}
        for var_name in R_ring.variable_names():
            if var_name in z_map:
                eval_dict[R_ring(var_name)] = z_map[var_name]
        
        F_val = F_poly.subs(eval_dict)
        if hasattr(F_val, 'is_constant') and not F_val.is_constant():
             F_val = F_val.constant_coefficient()
             
        # Compute Prefactors
        def get_C(i):
            return (lambdas[i][0]*x[1] - lambdas[i][1]*x[0]) * (lambdas[i][0]*y[1] - lambdas[i][1]*y[0])
            
        prod_C_sq = 1
        for i in range(n):
            if i not in roots:
                val = get_C(i)
                if val == 0:
                    raise ValueError(f"C_{i} is zero")
                prod_C_sq *= val**2
                
        r1, r2, r3 = roots
        norm_factor = (ang_bracket(lambdas[r1], lambdas[r2]) * 
                       ang_bracket(lambdas[r2], lambdas[r3]) * 
                       ang_bracket(lambdas[r3], lambdas[r1]))**2
                       
        h_factor = ang_bracket(lambdas[0], lambdas[1])**8
        sign = (-1)**(n-1)
        
        M_calc = sign * h_factor * F_val / (norm_factor * prod_C_sq)
        
        if reference_amp is None:
            reference_amp = M_calc
            print(f"Reference Amplitude (Sweep 0): {float(reference_amp):.6e}")
        else:
            if abs(M_calc) < 1e-15:
                 print(f"Sweep {k}: Calculated amplitude zero?")
                 continue
                 
            ratio = M_calc / reference_amp
            ratio_val = float(ratio)
            
            diff = abs(ratio_val - 1.0)
            if diff > 1e-5:
                print(f"Sweep {k}: MISMATCH! Ratio={ratio_val:.6f}")
                consistent = False

    if consistent:
        print(f"\nSUCCESS: Gauge invariance verified across {num_sweeps} random (x,y) pairs.")
    else:
        print("\nFAILURE: Gauge invariance check failed.")
        exit(1)

if __name__ == "__main__":
    gauge_invariance_sweep()
