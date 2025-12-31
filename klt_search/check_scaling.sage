#!/usr/bin/env sage
from sage.all import *
load('src/hodges.sage')
load('src/klt.sage')

def check_scaling():
    print("Checking Scaling Properties (Z -> 2*Z)")
    
    twistor = MomentumTwistor(n=6, check_domain=False)
    
    # 1. Base Values
    h_base, _ = hodges_6pt_mhv(twistor)
    k_base, _ = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
    
    if h_base is None or k_base is None:
        print("Base computation failed")
        return
        
    # 2. Scaled Values
    # Create new twistor with 2*Z
    Z_scaled = [2 * z for z in twistor.Z]
    twistor_scaled = MomentumTwistor(n=6, Z=Z_scaled, check_domain=False)
    
    h_scaled, _ = hodges_6pt_mhv(twistor_scaled)
    k_scaled, _ = gravity_6pt_mhv_klt(twistor_scaled, mandelstam_invariant)
    
    # 3. Ratios
    ratio_h = h_scaled / h_base
    ratio_k = k_scaled / k_base
    
    print(f"Hodges Ratio (2Z/Z): {float(ratio_h):.4e}")
    print(f"KLT Ratio (2Z/Z):    {float(ratio_k):.4e}")
    
    # Analyze Powers of 2
    pow_h = log(abs(ratio_h)) / log(2)
    pow_k = log(abs(ratio_k)) / log(2)
    
    print(f"Hodges Scaling Power: {float(pow_h):.2f}")
    print(f"KLT Scaling Power:    {float(pow_k):.2f}")
    
    diff_pow = pow_k - pow_h
    print(f"Difference in Powers: {float(diff_pow):.2f}")

if __name__ == "__main__":
    check_scaling()





