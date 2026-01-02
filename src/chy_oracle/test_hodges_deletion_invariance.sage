import sys
import os
from sage.all import *

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical
from itertools import combinations

def test_deletion_invariance():
    print("Testing Hodges Deletion Invariance...")
    
    n = 6
    seed = 123
    lambdas, tildes = sample_spinors_from_twistor(seed=seed, n=n)
    
    # Reference value using (0, 1, 2)
    ref_val, status = hodges_npt_mhv_canonical(lambdas, tildes, (0, 1), deletion_set=(0, 1, 2))
    if status != "ok":
        print(f"Failed to compute reference: {status}")
        return

    print(f"Reference value (0,1,2): {ref_val}")
    
    # Check all combinations
    triples = list(combinations(range(n), 3))
    
    match_count = 0
    fail_count = 0
    
    for triple in triples:
        val, status = hodges_npt_mhv_canonical(lambdas, tildes, (0, 1), deletion_set=triple)
        
        if status != "ok":
            print(f"  Triple {triple}: Error {status}")
            continue
            
        if val == ref_val:
            match_count += 1
        else:
            print(f"  Triple {triple}: Mismatch! Got {val}, Expected {ref_val}")
            if val == -ref_val:
                print("    (It is negative of reference)")
            fail_count += 1
            
    print(f"\nResults: {match_count} matches, {fail_count} failures out of {len(triples)} triples.")
    
    if fail_count == 0:
        print("PASS: Deletion invariance verified.")
    else:
        print("FAIL: Deletion invariance broken.")

if __name__ == "__main__":
    test_deletion_invariance()




