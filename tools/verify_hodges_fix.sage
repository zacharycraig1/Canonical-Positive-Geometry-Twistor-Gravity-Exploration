#!/usr/bin/env sage
# =============================================================================
# VERIFY HODGES FIX: Test deletion independence after adding sigma
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('tools/diagnostic_probe_1_hodges_invariance.sage')

def verify_fix():
    """Verify Hodges is now deletion-independent."""
    print("="*70)
    print("VERIFYING HODGES FIX: Deletion Independence Test")
    print("="*70)
    
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    baseline_rows = [0, 1, 2]
    baseline_cols = [3, 4, 5]
    baseline_x = 0
    baseline_y = 5
    
    baseline_result = hodges_with_deletion(twistor, baseline_rows, baseline_cols, baseline_x, baseline_y)
    baseline_H = baseline_result[0]
    
    print(f"Baseline H = {baseline_H}")
    
    # Test multiple deletions
    test_cases = [
        ([0,1,2], [3,4,5], 0, 5),  # Same as baseline
        ([0,1,2], [0,1,2], 0, 5),  # Symmetric
        ([0,1,3], [2,4,5], 0, 5),  # Mixed
        ([0,2,4], [1,3,5], 0, 5),  # Mixed
        ([3,4,5], [3,4,5], 0, 5),  # Other symmetric
    ]
    
    print(f"\nTesting different deletions (all should give same H):")
    all_match = True
    for rows_del, cols_del, x_ref, y_ref in test_cases:
        result = hodges_with_deletion(twistor, list(rows_del), list(cols_del), x_ref, y_ref)
        H_val = result[0]
        
        if H_val is None:
            print(f"  Failed: rows={rows_del}, cols={cols_del}")
            continue
        
        ratio = H_val / baseline_H if baseline_H != 0 else None
        match = (ratio == 1) if ratio is not None else False
        
        print(f"  rows={rows_del}, cols={cols_del}: H={H_val}, ratio={ratio}, match={match}")
        
        if not match:
            all_match = False
    
    print(f"\n{'='*70}")
    if all_match:
        print("*** SUCCESS: All deletions give the same result!")
    else:
        print("*** FAILURE: Deletions still give different results")
        print("Need to check sigma computation or other factors")
    print(f"{'='*70}")

if __name__ == '__main__':
    verify_fix()










