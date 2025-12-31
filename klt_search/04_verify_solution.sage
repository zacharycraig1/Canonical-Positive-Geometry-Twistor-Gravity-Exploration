#!/usr/bin/env sage
from sage.all import *
import numpy as np
import os

def verify_solution():
    print("Step 4: Verifying Solution")
    
    if not os.path.exists("klt_search/solution_x.npy"):
        print("Error: Solution not found.")
        return
        
    x = np.load("klt_search/solution_x.npy")
    basis_data = load("klt_search/basis_metadata.sobj")
    monomials = basis_data['monomials']
    perms = basis_data['perms']
    
    # Generate new test points
    num_test = 20
    print(f"Testing on {num_test} new points...")
    
    load('src/hodges.sage')
    load('src/klt.sage')
    
    errors = []
    
    for i in range(num_test):
        twistor = MomentumTwistor(n=6, check_domain=False)
        
        try:
            # Oracle: Hodges Reduced * Normalization
            load('src/kinematics_map.sage')
            load('src/spinor_helicity.sage')
            
            lambdas, tilde_lambdas, _ = extract_spinors_from_twistor(twistor)
            if lambdas is None: continue
            adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
            
            amp_hodges, _ = hodges_6pt_mhv_reduced(adapter)
            if amp_hodges is None: continue
            
            ang01 = twistor.get_angle(0, 1)
            if ang01 == 0: continue
            norm_factor = - (ang01 ** 8)
            target = float(amp_hodges * norm_factor)
            
            # Prediction
            # Reconstruct Ansatz value
            
            # Compute vars
            s_vars = []
            pairs = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0)]
            for a,b in pairs: s_vars.append(mandelstam_invariant(twistor, a, b))
            triples = [(0,1,2), (1,2,3), (2,3,4)]
            for a,b,c in triples: 
                s1 = mandelstam_invariant(twistor, a, b)
                s2 = mandelstam_invariant(twistor, b, c)
                s3 = mandelstam_invariant(twistor, c, a)
                s_vars.append(s1 + s2 + s3)
                
            # Monomial values
            mon_vals = []
            for exps in monomials:
                val = 1
                for idx, e in enumerate(exps):
                    if e > 0: val *= (s_vars[idx] ** e)
                mon_vals.append(float(val))
            mon_vec = np.array(mon_vals)
                
            # PT values
            pt_vals = []
            for p in perms:
                order = [0] + list(p) + [4, 5]
                pt_vals.append(float(parke_taylor_6pt_mhv(twistor, order)))
            pt_vec = np.array(pt_vals)
            
            # Predict
            # Row = kron(pt_cross, mon_vec)
            # Pred = Row . x
            # x is shape (total_unknowns,)
            # Reshape x to (len_pt_cross, len_mon)
            # x_mat = x.reshape((36, 220))
            # Pred = pt_cross . (x_mat . mon_vec)
            
            pt_cross = np.outer(pt_vec, pt_vec).flatten()
            
            # Efficient computation
            # prediction = sum( pt_cross[k] * sum(x[k, m] * mon[m]) )
            
            # Reshape x once
            if i == 0:
                x_reshaped = x.reshape((36, len(monomials)))
                
            # Term 1: Weights vector w = x_reshaped @ mon_vec
            weights = x_reshaped @ mon_vec
            
            # Term 2: Dot with PT cross
            prediction = np.dot(pt_cross, weights)
            
            err = abs(prediction - target)
            rel_err = err / abs(target) if target != 0 else err
            
            errors.append(rel_err)
            print(f"Point {i}: Target={target:.4e}, Pred={prediction:.4e}, RelErr={rel_err:.4e}")
            
        except Exception as e:
            print(f"Point {i}: Error {e}")
            continue
            
    avg_rel_err = np.mean(errors)
    print(f"\nAverage Relative Error: {avg_rel_err:.4e}")
    if avg_rel_err < 1e-3:
        print("SUCCESS: Ansatz generalizes to new data!")
    else:
        print("FAILURE: Ansatz does not generalize.")

if __name__ == "__main__":
    verify_solution()

