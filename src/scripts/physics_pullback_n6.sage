import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.posgeom.forest_polytope_n6 import get_forest_polynomial_n6
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.chy_oracle.amplitude_spinor import ang_bracket

def physics_pullback_n6_check():
    print("Checking Exact Relation for n=6...")
    n = 6
    roots = [0, 1, 2]
    
    # Run 1 trial (slow?)
    
    # 1. Kinematics
    lambdas, tildes = sample_spinors_from_twistor(n=n, seed=42)
    x = [1, 0]
    y = [0, 1]
    
    # 2. Physics Variables
    z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
    
    # 3. Forest Polynomial F(z)
    print("Generating n=6 polynomial...")
    F_poly = get_forest_polynomial_n6(roots)
    print("Evaluating polynomial...")
    
    R = F_poly.parent()
    eval_dict = {R(k): v for k, v in z_map.items() if str(k).startswith('z_')}
    F_val = F_poly.subs(eval_dict)
    
    # 4. Amplitude
    print("Computing Amplitude...")
    M_amp, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=roots)
    if M_amp is None:
        print(f"Error: {status}")
        return
        
    # 5. Theoretical Prefactors
    def get_C(i):
        return (lambdas[i][0]*x[1] - lambdas[i][1]*x[0]) * (lambdas[i][0]*y[1] - lambdas[i][1]*y[0])
        
    prod_C_sq = 1
    for k in range(n):
        if k not in roots:
            prod_C_sq *= get_C(k)**2
            
    r1, r2, r3 = roots
    norm_factor = (ang_bracket(lambdas[r1], lambdas[r2]) * 
                   ang_bracket(lambdas[r2], lambdas[r3]) * 
                   ang_bracket(lambdas[r3], lambdas[r1]))**2
                   
    h_factor = ang_bracket(lambdas[0], lambdas[1])**8
    sign = (-1)**(n-1)
    
    M_expected = sign * h_factor * F_val / (norm_factor * prod_C_sq)
    
    # Comparison
    ratio = M_amp / M_expected
    print(f"Ratio = {float(ratio):.6f}  (M_amp={float(M_amp):.2e}, M_exp={float(M_expected):.2e})")
    
    if abs(float(ratio) - 1.0) < 1e-5:
        print("SUCCESS: Exact identity verified for n=6!")
    else:
        print("FAILURE.")

if __name__ == "__main__":
    physics_pullback_n6_check()



