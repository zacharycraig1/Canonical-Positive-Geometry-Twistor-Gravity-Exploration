import sys
import os
# sys.path.append(os.getcwd())

from sage.all import *

# Load modules (Sage style)
# Assuming running from repo root
load("src/spinor_sampling.sage")
load("src/hodges.sage")
load("src/kinematics_map.sage")
load("tools/certify_sample.sage")

import random

def test_hodges_invariance():
    print("Running Phase 0: Hodges Invariance Test (Reference Spinors & Deletion Sets)")
    print("-------------------------------------------------------------------------")
    
    num_samples = 50
    failures = 0
    skips = 0
    
    # Random seed for overall test suite
    set_random_seed(12345)
    
    for i in range(num_samples):
        # 1. Sample kinematics
        seed_sample = random.randint(1, 1000000)
        res = sample_spinor_helicity_conserving(n=6, seed=seed_sample)
        if res is None:
            skips += 1
            continue
            
        lambdas, tilde_lambdas = res
        adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
        
        # 2. Compute baseline (default refs, default deletion set)
        baseline, reason = hodges_6pt_mhv_reduced(adapter)
        if baseline is None:
            # Skip degenerate samples
            # print(f"Sample {i}: Skipped ({reason})")
            skips += 1
            continue
            
        # 3. Test Reference Spinor Invariance
        # Try 5 random pairs of reference spinors
        for k in range(5):
            ref_x = vector(QQ, [QQ(random.randint(1, 10)), QQ(random.randint(1, 10))])
            ref_y = vector(QQ, [QQ(random.randint(1, 10)), QQ(random.randint(1, 10))])
            
            # Check non-degeneracy of refs locally
            # (In a real test we might just retry sampling refs if they fail, 
            # but hodges_6pt_mhv_reduced handles validation)
            
            val, reason_ref = hodges_6pt_mhv_reduced(adapter, ref_spinors=(ref_x, ref_y))
            
            if val is None:
                # If refs are bad, just skip this specific ref check
                continue
                
            if val != baseline:
                print(f"FAILURE: Reference spinor invariance violation at sample {i}")
                print(f"Baseline: {baseline}")
                print(f"Current:  {val}")
                print(f"Diff:     {baseline - val}")
                
                # Dump certificate
                cert = {
                    "type": "invariance_failure_ref_spinor",
                    "seed": seed_sample,
                    "lambdas": lambdas,
                    "tilde_lambdas": tilde_lambdas,
                    "baseline": baseline,
                    "comparison_value": val,
                    "ref_x": ref_x,
                    "ref_y": ref_y,
                    "diff": baseline - val
                }
                fpath = save_certificate(cert, "fail_ref_inv")
                print(f"Certificate saved to {fpath}")
                failures += 1
                break # Stop checking refs for this sample

        if failures > 0: break

        # 4. Test Deletion Set Invariance
        # Try 5 random deletion sets (disjoint rows/cols)
        # rows_delete (size 3), cols_delete (size 3)
        # For n=6, disjoint sets of size 3 partition the indices.
        # Actually Hodges formula requires |Φ|^{rst}_{ijk}.
        # i,j,k are rows removed. r,s,t are cols removed.
        # They don't strictly have to be disjoint from each other in the sense of {i,j,k} vs {r,s,t},
        # but the standard reduced form is defined for any deletion set of size 3.
        # The formula includes sigma(ijk, rst).
        
        all_indices = [0, 1, 2, 3, 4, 5]
        import itertools
        
        # We'll just pick random subsets of size 3 for rows and cols
        # There are (6 choose 3) = 20 choices for rows, 20 for cols. 400 combinations.
        # We test 5.
        
        for k in range(5):
            rows_del = tuple(sorted(random.sample(all_indices, 3)))
            cols_del = tuple(sorted(random.sample(all_indices, 3)))
            
            val, reason_del = hodges_6pt_mhv_reduced(adapter, deletion_set=(rows_del, cols_del))
            
            if val is None:
                continue
                
            if val != baseline:
                print(f"FAILURE: Deletion set invariance violation at sample {i}")
                print(f"Baseline: {baseline}")
                print(f"Current:  {val} (rows={rows_del}, cols={cols_del})")
                print(f"Diff:     {baseline - val}")
                
                cert = {
                    "type": "invariance_failure_deletion_set",
                    "seed": seed_sample,
                    "lambdas": lambdas,
                    "tilde_lambdas": tilde_lambdas,
                    "baseline": baseline,
                    "comparison_value": val,
                    "rows_delete": rows_del,
                    "cols_delete": cols_del,
                    "diff": baseline - val
                }
                fpath = save_certificate(cert, "fail_del_inv")
                print(f"Certificate saved to {fpath}")
                failures += 1
                break
        
        if failures > 0: break
        
        if (i+1) % 10 == 0:
            print(f"Verified {i+1} samples...")

    if failures == 0:
        print(f"SUCCESS: Passed {num_samples} samples (skipped {skips}). Invariance holds.")
    else:
        print("FAILED: See certificates.")

if __name__ == "__main__":
    test_hodges_invariance()

