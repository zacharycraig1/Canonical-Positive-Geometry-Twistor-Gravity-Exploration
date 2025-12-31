import sys
import os
import argparse
from sage.all import *

# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src'))

# Import from robust script
# This assumes phaseC2_factorization_robust.py is in the python path or same dir
try:
    from phaseC2_factorization_robust import (
        MomentumTwistor, hodges_6pt_mhv, get_s123, 
        bcfw_shift_spinors, spinors_to_twistors, solve_s123_bcfw
    )
except ImportError:
    # Try local import if running from scripts dir
    sys.path.append(os.path.join(os.getcwd(), 'src', 'scripts'))
    from phaseC2_factorization_robust import (
        MomentumTwistor, hodges_6pt_mhv, get_s123, 
        bcfw_shift_spinors, spinors_to_twistors, solve_s123_bcfw
    )

def analyze_scaling(seed, shift_a=0, shift_b=3):
    print(f"Analyzing scaling for Seed {seed}, Shift ({shift_a}, {shift_b})...")
    
    tw = MomentumTwistor(n=6, seed=seed)
    if not tw.domain_ok:
        print("Domain error in base kinematics.")
        return

    lambdas = [tw.get_lambda(i) for i in range(6)]
    tilde_lambdas = [tw.get_tilde_lambda(i) for i in range(6)]
    
    if any(x is None for x in tilde_lambdas):
        print("Invalid tilde lambdas.")
        return

    # Solve z*
    z_star = solve_s123_bcfw(lambdas, tilde_lambdas, shift_a, shift_b)
    if z_star is None:
        print("s123 is constant.")
        return
        
    print(f"z* = {z_star}")
    
    # Probe
    epsilons = [QQ(1)/10**k for k in range(3, 9)]
    
    print(f"{'eps':<10} | {'s123':<12} | {'|M|':<12} | {'M*s':<12}")
    print("-" * 50)
    
    scaling_data = []
    
    for eps in epsilons:
        z_probe = z_star + eps
        L_probe, Lt_probe = bcfw_shift_spinors(lambdas, tilde_lambdas, shift_a, shift_b, z_probe)
        tw_probe = spinors_to_twistors(L_probe, Lt_probe)
        
        if tw_probe is None:
            print(f"{float(eps):<10.1e} | Twistor Recon Failed")
            continue
            
        s_val = get_s123(tw_probe)
        M_val, reason = hodges_6pt_mhv(tw_probe)
        
        if M_val is None:
            print(f"{float(eps):<10.1e} | {str(s_val):<12} | Error: {reason}")
        else:
            abs_m = float(abs(M_val))
            s_float = float(s_val)
            prod = abs_m * s_float
            print(f"{float(eps):<10.1e} | {s_float:<12.4e} | {abs_m:<12.4e} | {prod:<12.4e}")
            
            scaling_data.append((s_float, abs_m))
            
    # Alpha
    if len(scaling_data) >= 2:
        s1, m1 = scaling_data[-1]
        s2, m2 = scaling_data[0]
        import math
        try:
            alpha = (math.log(m1) - math.log(m2)) / (math.log(abs(s1)) - math.log(abs(s2)))
            print(f"\nEstimated Alpha: {alpha:.4f}")
        except:
            print("\nCould not estimate alpha.")

if __name__ == "__main__":
    # Hardcoded default or use seeds from previous run
    # Example seed 5001 had alpha ~ 0
    analyze_scaling(seed=5001)



