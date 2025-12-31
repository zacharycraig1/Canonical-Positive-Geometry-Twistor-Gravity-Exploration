import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.posgeom.forest_polytope import get_forest_polynomial
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import ang_bracket

def residue_factorization_n6():
    print("Checking Residue Factorization for n=6 (Collinear Limit z_34 -> 0)...")
    n = 6
    roots = [0, 1, 2]
    
    # 1. Generate base kinematics
    lambdas, tildes = sample_spinors_from_twistor(n=n, seed=42)
    x = [1, 0]
    y = [0, 1]
    
    # Build Polynomial once
    print("Building polynomial...")
    F_poly = get_forest_polynomial(n, roots)
    R_ring = F_poly.parent()
    print("Polynomial built.")
    
    # Define helper to evaluate M(eps)
    def get_M(eps):
        L = [list(l) for l in lambdas] 
        L[4] = [L[3][0] + eps, L[3][1] + eps] 
        
        # Eval edge vars
        try:
            z = eval_edge_vars_from_spinors(L, tildes, x, y)
        except ValueError:
            return 0
            
        # Eval F
        # Filter keys
        subs_dict = {R_ring(k): v for k,v in z.items() if str(k) in R_ring.variable_names()}
        val_F = F_poly.subs(subs_dict)
        
        # Handle full evaluation
        if hasattr(val_F, 'constant_coefficient'):
             val_F = val_F.constant_coefficient()
        
        # Prefactors
        pCsq = 1
        for i in range(n):
            if i not in roots:
                val = (L[i][0]*x[1] - L[i][1]*x[0]) * (L[i][0]*y[1] - L[i][1]*y[0])
                pCsq *= val**2
                
        nr = (ang_bracket(L[roots[0]], L[roots[1]]) * ang_bracket(L[roots[1]], L[roots[2]]) * ang_bracket(L[roots[2]], L[roots[0]]))**2
        hf = ang_bracket(L[0], L[1])**8
        sign = (-1)**(n-1)
        
        if nr * pCsq == 0: return float('inf')
        return sign * hf * val_F / (nr * pCsq)

    eps1 = 1e-5
    eps2 = 1e-8
    
    print(f"Evaluating at eps1={eps1}...")
    m1 = get_M(eps1)
    print(f"Evaluating at eps2={eps2}...")
    m2 = get_M(eps2)
    
    print(f"M(eps1) = {abs(m1):.4e}")
    print(f"M(eps2) = {abs(m2):.4e}")
    
    if m1 == 0 or m2 == 0:
        print("Amplitude is zero?")
        return
        
    ratio = abs(m2/m1)
    eps_ratio = eps1/eps2
    log_ratio = float(log(ratio)/log(eps_ratio))
    
    print(f"Estimated Pole Order: {log_ratio:.4f}")
    
    if abs(log_ratio - 3.0) < 0.2:
        print("SUCCESS: Cubic pole detected (Order ~ 3).")
    elif abs(log_ratio - 1.0) < 0.2:
        print("SUCCESS: Simple pole detected (Order ~ 1).")
    elif abs(log_ratio - 2.0) < 0.2:
        print("SUCCESS: Double pole detected (Order ~ 2).")
    else:
        print(f"Observed pole order {log_ratio:.4f}.")

if __name__ == "__main__":
    residue_factorization_n6()
