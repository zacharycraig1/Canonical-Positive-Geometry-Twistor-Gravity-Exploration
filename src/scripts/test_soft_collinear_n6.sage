import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.posgeom.forest_polytope import get_forest_polynomial
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import ang_bracket

def test_soft_collinear_n6():
    print("Testing Soft and Collinear Limits for n=6...")
    n = 6
    roots = [0, 1, 2]
    
    # 1. Base Kinematics
    lambdas, tildes = sample_spinors_from_twistor(n=n, seed=55)
    x, y = [1, 0], [0, 1]
    
    # Pre-build Poly
    F_poly = get_forest_polynomial(n, roots)
    R_ring = F_poly.parent()
    
    def eval_M(L, T):
        try:
            z = eval_edge_vars_from_spinors(L, T, x, y)
        except ValueError:
            return None
        
        eval_dict = {R_ring(k): v for k,v in z.items() if str(k) in R_ring.variable_names()}
        val = F_poly.subs(eval_dict)
        if hasattr(val, 'constant_coefficient'): val = val.constant_coefficient()
        
        # Prefactors... (Simplified check: just check if M diverges or vanishes)
        return val

    # 2. Collinear Limit: 3 || 4
    print("\n--- Collinear Limit <34> -> 0 ---")
    eps = 1e-8
    L_coll = [list(l) for l in lambdas]
    L_coll[4] = [L_coll[3][0]+eps, L_coll[3][1]+eps]
    
    val_coll = eval_M(L_coll, tildes)
    print(f"M(eps) ~ {abs(val_coll):.4e}")
    
    if abs(val_coll) > 1e3:
        print("PASS: Diverges (Pole detected).")
    else:
        print("FAIL: Did not diverge significantly.")
        
    # 3. Soft Limit: k5 -> 0
    print("\n--- Soft Limit k5 -> 0 ---")
    # Scale lambda_5 and tilde_5 by sqrt(eps) -> k5 -> eps*k5
    eps_soft = 1e-8
    scale = sqrt(eps_soft)
    L_soft = [list(l) for l in lambdas]
    T_soft = [list(t) for t in tildes]
    L_soft[5] = [scale*v for v in L_soft[5]]
    T_soft[5] = [scale*v for v in T_soft[5]]
    
    val_soft = eval_M(L_soft, T_soft)
    print(f"M(soft) ~ {abs(val_soft):.4e}")
    
    # Gravity soft factor goes as 1/k^3? Or 1/k?
    # Weinberg soft factor S ~ 1.
    # Amplitude M_n -> S * M_{n-1}.
    # So M should be finite? No, S diverges as 1/k.
    # We expect divergence.
    
    if abs(val_soft) > 1e3:
        print("PASS: Diverges (Soft Pole detected).")
    else:
        print("FAIL: Did not diverge.")

if __name__ == "__main__":
    test_soft_collinear_n6()
