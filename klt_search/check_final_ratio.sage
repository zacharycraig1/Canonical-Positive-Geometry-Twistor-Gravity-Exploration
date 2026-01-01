#!/usr/bin/env sage
from sage.all import *
load('src/hodges.sage')
load('src/klt.sage')
load('src/kinematics_map.sage')
load('src/spinor_helicity.sage')

def check_ratio():
    print("Checking Ratio: KLT / (HodgesRed * <01>^8)")
    
    for i in range(5):
        twistor = MomentumTwistor(n=6, check_domain=False)
        
        # 1. KLT
        # Use mandelstam_invariant from hodges.sage (uses 4-bracket s)
        # Or should we use physical s from extracted spinors?
        # KLT kernel definition usually uses s_ij = (ki+kj)^2.
        # This should be physical s.
        # Let's use extraction-based s.
        
        lambdas, tilde_lambdas, _ = extract_spinors_from_twistor(twistor)
        if lambdas is None: continue
        adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
        
        # Define physical mandelstam for KLT
        def mandelstam_phys(tw, i, j):
            return adapter.get_angle(i,j) * adapter.get_square(i,j)
            
        klt_val, reason = gravity_6pt_mhv_klt(adapter, mandelstam_phys)
        if reason != "ok": continue
        
        # 2. Hodges Reduced
        # Uses adapter
        h_val, reason = hodges_6pt_mhv_reduced(adapter)
        if h_val is None: continue
        
        # 3. Normalization
        ang01 = adapter.get_angle(0, 1)
        norm = -(ang01 ** 8)
        
        target = h_val * norm
        
        if target == 0: continue
        
        ratio = klt_val / target
        print(f"Point {i}: Ratio = {float(ratio):.6f}")

if __name__ == "__main__":
    check_ratio()








