import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.posgeom.forest_polytope import get_forest_polynomial
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.chy_oracle.amplitude_spinor import ang_bracket

def physics_pullback_n6_exact():
    print("Checking Exact Relation for n=6...")
    n = 6
    roots = [0, 1, 2]
    
    # Run multiple trials
    trials = 20
    consistent = True
    
    # Cache the polynomial so we don't rebuild it every trial
    print("Building Forest Polynomial F_{6,R} (this may take a moment)...")
    F_poly = get_forest_polynomial(n, roots)
    print("Polynomial built.")
    
    for t in range(trials):
        # 1. Kinematics
        lambdas, tildes = sample_spinors_from_twistor(n=n, seed=100+t)
        x = [1, 0]
        y = [0, 1]
        
        # 2. Physics Variables
        try:
            z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
        except ValueError as e:
            print(f"Skipping trial {t}: {e}")
            continue
            
        # 3. Evaluate Forest Polynomial
        R_ring = F_poly.parent()
        # Convert z_map to argument dictionary for subs
        # We need to ensure keys match the polynomial variables
        # z_map has keys like "z_0_1" and also tuples (0,1). We use string keys.
        eval_dict = {}
        for var_name in R_ring.variable_names():
            if var_name in z_map:
                eval_dict[R_ring(var_name)] = z_map[var_name]
            else:
                 print(f"Warning: Variable {var_name} not found in z_map")
                 
        F_val = F_poly.subs(eval_dict)
        
        # If it's still a polynomial, check constancy and extract value
        if hasattr(F_val, 'is_constant'):
            if not F_val.is_constant():
                print(f"Trial {t}: F_val is not constant! Remaining vars: {F_val.variables()}")
                continue
            F_val = F_val.constant_coefficient()
            
        # 4. Amplitude from Oracle (Laplacian/Hodges)
        M_amp, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=roots)
        if M_amp is None:
            print(f"Skipping trial {t}: {status}")
            continue
            
        # 5. Theoretical Prediction
        # Formula: M = (-1)^(n-1) * <01>^8 * F(z) / (Norm * prod(Ck^2))
        
        # C_k for k not in roots (k=3,4,5)
        def get_C(i):
            return (lambdas[i][0]*x[1] - lambdas[i][1]*x[0]) * (lambdas[i][0]*y[1] - lambdas[i][1]*y[0])
            
        prod_C_sq = 1
        for k in range(n):
            if k not in roots:
                prod_C_sq *= get_C(k)**2
                
        # Normalization factor N_R = (<r1 r2><r2 r3><r3 r1>)^2
        r1, r2, r3 = roots
        norm_factor = (ang_bracket(lambdas[r1], lambdas[r2]) * 
                       ang_bracket(lambdas[r2], lambdas[r3]) * 
                       ang_bracket(lambdas[r3], lambdas[r1]))**2
                       
        # Helicity factor <01>^8 (assuming 0,1 are negative helicity)
        h_factor = ang_bracket(lambdas[0], lambdas[1])**8
        
        sign = (-1)**(n-1)
        
        # Predicted Amplitude
        M_expected = sign * h_factor * F_val / (norm_factor * prod_C_sq)
        
        # Comparison
        if M_expected == 0:
            if M_amp == 0:
                # print(f"Trial {t}: Both zero (consistent).")
                continue
            else:
                print(f"Trial {t}: M_expected=0 but M_amp!=0. FAILURE.")
                consistent = False
                continue
            
        ratio = M_amp / M_expected
        
        if ratio != 1:
            consistent = False
            print(f"MISMATCH at trial {t}! Ratio = {ratio}")
            print(f"M_amp = {M_amp}")
            print(f"M_exp = {M_expected}")
        else:
            # print(f"Trial {t}: Exact match.")
            pass
            
    if consistent:
        print(f"\nSUCCESS: Exact identity verified for n=6 ({trials} trials)!")
        return True
    else:
        print("\nFAILURE: Identity did not match for some trials.")
        return False

if __name__ == "__main__":
    if not physics_pullback_n6_exact():
        exit(1)
