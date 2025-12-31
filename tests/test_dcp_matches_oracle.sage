import sys
import os
# sys.path.append(os.getcwd())

from sage.all import *
load("src/dcp_eval.sage")
load("src/oracle_gravity_mhv6.sage")
load("src/spinor_sampling.sage")

def test_dcp_matches_oracle():
    print("Running Phase 4: DCP Form vs Oracle Test")
    print("----------------------------------------")
    
    # 1. Setup Form Infrastructure
    C = generate_channels_n6()
    triples = generate_triples(C)
    print(f"Generated {len(C)} channels and {len(triples)} basis triples.")
    
    # 2. Create a Dummy Candidate
    # Just pick the first triple (0,1,2) => coefficient 1.
    # Corresponds to 1 / (s_{ch0} * s_{ch1} * s_{ch2})
    # ch0 = C[0] = (1,2) => s_01
    # ch1 = C[1] = (1,3) => s_02
    # ch2 = C[2] = (1,4) => s_03
    # This is likely NOT the gravity amplitude, but it tests the pipeline.
    
    candidate = vector(QQ, len(triples))
    candidate[0] = 1 
    
    # 3. Evaluate on Samples
    num_samples = 5
    
    import random
    random.seed(int(101))
    
    for i in range(num_samples):
        res = sample_spinor_helicity_conserving(n=6, seed=random.randint(1,10000))
        if res is None: continue
        l, lt = res
        
        # Evaluate Form
        form_val = eval_form_on_spinors(candidate, l, lt, C, triples)
        
        # Evaluate Oracle (Hodges)
        oracle_res = oracle_M6(l, lt)
        if oracle_res["value"] is None: continue
        oracle_val = oracle_res["value"]
        
        # Compare
        ratio = form_val / oracle_val if oracle_val != 0 else 0
        
        print(f"Sample {i}: Form={float(abs(form_val)):.2e}, Oracle={float(abs(oracle_val)):.2e}, Ratio={ratio}")
        
    print("\nTest complete. (Mismatch expected as candidate is dummy)")

if __name__ == "__main__":
    test_dcp_matches_oracle()






