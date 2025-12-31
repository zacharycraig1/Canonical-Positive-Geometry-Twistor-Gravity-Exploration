#!/usr/bin/env sage
from sage.all import *
load('src/hodges.sage') # For MomentumTwistor
load('src/klt.sage')
load('src/kinematics_map.sage')
load('src/spinor_helicity.sage')

def check_klt_symmetry():
    print("Checking KLT Amplitude Symmetry")
    
    twistor = MomentumTwistor(n=6, check_domain=False)
    lambdas, tilde_lambdas, _ = extract_spinors_from_twistor(twistor)
    if lambdas is None: return
    adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
    
    # 1. Base Calculation
    def m_func(tw, i, j):
        return adapter.get_angle(i,j) * adapter.get_square(i,j)
        
    base_val, reason = gravity_6pt_mhv_klt(adapter, m_func)
    print(f"Base Value: {base_val}")
    
    # 2. Permuted Calculation (Swap 0 and 1)
    # We need to permute the ADAPTER inputs
    # Swap indices 0 and 1 in lambdas and tildelambdas
    l_perm = list(lambdas)
    l_perm[0], l_perm[1] = l_perm[1], l_perm[0]
    
    tl_perm = list(tilde_lambdas)
    tl_perm[0], tl_perm[1] = tl_perm[1], tl_perm[0]
    
    adapter_perm = SpinorHelicityAdapter(l_perm, tl_perm)
    
    # Recalculate
    # Note: KLT implementation uses fixed legs 0, 4, 5 internally?
    # No, gravity_6pt_mhv_klt uses fixed_leg_1=0, fixed_leg_5=4, fixed_leg_6=5.
    # If we swap particles 0 and 1 in the input, the physics should change?
    # No, M(1,2,...) = M(2,1,...).
    # But we feed the permuted spinors.
    # So the function should return the same value (since M is symmetric).
    
    perm_val, reason = gravity_6pt_mhv_klt(adapter_perm, m_func)
    print(f"Permuted (0<->1): {perm_val}")
    
    ratio = perm_val / base_val if base_val != 0 else 0
    print(f"Ratio: {ratio}")
    
    if abs(ratio - 1) < 1e-5:
        print("SYMMETRIC")
    else:
        print("NOT SYMMETRIC")

if __name__ == "__main__":
    check_klt_symmetry()

