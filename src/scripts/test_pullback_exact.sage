import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.scripts.physics_pullback_n6 import physics_pullback_n6_exact
from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.posgeom.forest_polytope import get_forest_polynomial
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.chy_oracle.amplitude_spinor import ang_bracket

def test_pullback_exact_harness():
    print("=== Exact Pullback Verification Harness ===")
    
    def check_n(n, roots=[0,1,2], trials=5):
        print(f"\nChecking n={n}...")
        F_poly = get_forest_polynomial(n, roots)
        consistent = True
        
        for t in range(trials):
            lambdas, tildes = sample_spinors_from_twistor(n=n, seed=200+t)
            x, y = [1, 0], [0, 1]
            
            try:
                z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
            except ValueError:
                continue
                
            # Eval F
            R_ring = F_poly.parent()
            eval_dict = {R_ring(k): v for k,v in z_map.items() if str(k) in R_ring.variable_names()}
            F_val = F_poly.subs(eval_dict)
            if hasattr(F_val, 'constant_coefficient'): F_val = F_val.constant_coefficient()
            
            # Eval Amp
            M_amp, _ = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=roots)
            if M_amp is None: continue
            
            # Predict
            def get_C(i):
                 return (lambdas[i][0]*x[1] - lambdas[i][1]*x[0]) * (lambdas[i][0]*y[1] - lambdas[i][1]*y[0])
            prod_C_sq = 1
            for k in range(n):
                if k not in roots: prod_C_sq *= get_C(k)**2
            
            r1, r2, r3 = roots
            norm = (ang_bracket(lambdas[r1], lambdas[r2]) * ang_bracket(lambdas[r2], lambdas[r3]) * ang_bracket(lambdas[r3], lambdas[r1]))**2
            hf = ang_bracket(lambdas[0], lambdas[1])**8
            sign = (-1)**(n-1)
            
            M_exp = sign * hf * F_val / (norm * prod_C_sq)
            
            if M_exp != 0:
                if M_amp / M_exp != 1:
                    print(f"  FAIL: Trial {t}, Ratio={M_amp/M_exp}")
                    consistent = False
            elif M_amp != 0:
                print(f"  FAIL: Trial {t}, Exp=0, Amp!=0")
                consistent = False
                
        if consistent: print(f"  PASS: n={n}")
        return consistent

    results = []
    results.append(check_n(4))
    results.append(check_n(5))
    results.append(check_n(6))
    
    if all(results):
        print("\nALL TESTS PASSED.")
        exit(0)
    else:
        print("\nSOME TESTS FAILED.")
        exit(1)

if __name__ == "__main__":
    test_pullback_exact_harness()
