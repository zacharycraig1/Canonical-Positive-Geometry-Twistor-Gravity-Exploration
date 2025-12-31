import sys
import os
# sys.path.append(os.getcwd())

from sage.all import *
load("src/kinematics_map.sage")
load("src/spinor_sampling.sage")

def test_channel_identities():
    print("Running Phase 1: Channel Identity Tests")
    print("---------------------------------------")
    
    num_samples = 50
    failures = 0
    
    for i in range(num_samples):
        res = sample_spinor_helicity_conserving(n=6, seed=i*100)
        if res is None: continue
        
        lambdas, tilde_lambdas = res
        
        # 1. Get s_ij
        sij = spinors_to_sij(lambdas, tilde_lambdas)
        
        # Check symmetry
        for a in range(6):
            for b in range(6):
                val_ab = sij.get((a,b), 0)
                val_ba = sij.get((b,a), 0)
                if val_ab != val_ba:
                    print(f"FAILURE: s_{a}{b} != s_{b}{a}")
                    failures += 1
                if a == b and val_ab != 0:
                     print(f"FAILURE: s_{a}{a} != 0")
                     failures += 1

        # 2. Get channels
        channels = spinors_to_channels(lambdas, tilde_lambdas)
        
        # Check sum consistency
        # s_ijk = s_ij + s_jk + s_ki
        for k_tuple, val in channels.items():
            if len(k_tuple) == 3:
                a, b, c = k_tuple
                sum_pairs = sij[(a,b)] + sij[(b,c)] + sij[(a,c)]
                if val != sum_pairs:
                    print(f"FAILURE: Channel sum mismatch for {k_tuple}")
                    failures += 1
                    
        if failures > 0: break
        
    if failures == 0:
        print(f"SUCCESS: {num_samples} samples verified for channel identities.")
    else:
        print("FAILED channel identities.")

if __name__ == "__main__":
    test_channel_identities()

