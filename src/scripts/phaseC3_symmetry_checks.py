import sys
import os
import random as rnd
import json
from sage.all import *

sys.path.append(os.path.join(os.getcwd(), 'src'))
from chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor

root_dir = os.getcwd()
hodges_path = os.path.join(root_dir, 'src', 'hodges.sage')
load(hodges_path)

def run_symmetry_check_v2():
    print("Phase C3 - Symmetry Check (Corrected)")
    
    # We only check permutations that preserve the {0, 1} set (negative helicities)
    # OR the {2,3,4,5} set (positive helicities).
    # i.e. S2 x S4 subgroup.
    
    # Actually, gravity amplitude is permutation symmetric if we track helicities.
    # But our formula implicitly assigns h=-2 to leg 0 (first arg) and leg 1 (second arg).
    # So if we want to check M is invariant under swapping 2 and 4, we swap inputs 2 and 4.
    # If we want to check swapping 0 and 1, we swap inputs 0 and 1.
    # This covers S2 x S4.
    
    log_data = []
    
    num_seeds = 10
    perms_per_seed = 50
    
    perms_passed = 0
    perms_total = 0
    
    seed_counter = 9100
    
    for _ in range(num_seeds):
        seed = seed_counter
        seed_counter += 1
        
        tw = MomentumTwistor(n=6, seed=seed)
        lambdas = [tw.get_lambda(i) for i in range(6)]
        tilde_lambdas = [tw.get_tilde_lambda(i) for i in range(6)]
        
        if any(x is None for x in tilde_lambdas): continue
        
        M_base, reason = hodges_6pt_mhv_spinor(lambdas, tilde_lambdas)
        if M_base is None: continue
        
        # S2 x S4 permutations
        # indices 0,1 map to 0,1
        # indices 2,3,4,5 map to 2,3,4,5
        
        sub_results = {'pass': 0, 'fail': 0}
        
        for __ in range(perms_per_seed):
            p01 = [0, 1]
            rnd.shuffle(p01)
            
            p2345 = [2, 3, 4, 5]
            rnd.shuffle(p2345)
            
            p = p01 + p2345
            
            L_p = [lambdas[i] for i in p]
            Lt_p = [tilde_lambdas[i] for i in p]
            
            M_p, reason_p = hodges_6pt_mhv_spinor(L_p, Lt_p)
            
            if M_p is None: continue
            
            ratio = M_p / M_base
            if abs(ratio - 1) < 1e-9:
                sub_results['pass'] += 1
            else:
                sub_results['fail'] += 1
                
        perms_passed += sub_results['pass']
        perms_total += sub_results['pass'] + sub_results['fail']
        
        print(f"Seed {seed}: Passed {sub_results['pass']}/{sub_results['pass']+sub_results['fail']}")
        log_data.append({'seed': seed, 'results': sub_results})
        
    print(f"\nTotal Passed: {perms_passed}/{perms_total}")
    
    with open('phaseC3_symmetry_v2.json', 'w') as f:
        json.dump(log_data, f, indent=2)

if __name__ == "__main__":
    run_symmetry_check_v2()
