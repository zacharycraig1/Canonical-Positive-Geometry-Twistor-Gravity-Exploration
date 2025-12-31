#!/usr/bin/env sage
# =============================================================================
# VERIFY COMPENSATION PATTERN ACROSS MULTIPLE DELETIONS
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('tools/diagnostic_probe_1_hodges_invariance.sage')
load('src/hodges_sigma.sage')

def verify_compensation():
    """Verify if ratio = sigma_ratio * c_ratio for multiple deletions."""
    print("="*70)
    print("VERIFYING COMPENSATION PATTERN")
    print("="*70)
    
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    x_ref = 0
    y_ref = 5
    
    # Baseline
    baseline_rows = [0, 1, 2]
    baseline_cols = [3, 4, 5]
    baseline_result = hodges_with_deletion(twistor, baseline_rows, baseline_cols, x_ref, y_ref)
    baseline_H = baseline_result[0]
    
    # Test cases
    test_cases = [
        ([0,1,2], [0,1,2], "symmetric"),
        ([0,1,3], [2,4,5], "mixed1"),
        ([0,2,4], [1,3,5], "mixed2"),
    ]
    
    print(f"\nBaseline: rows={baseline_rows}, cols={baseline_cols}, H={baseline_H}\n")
    
    all_match = True
    for test_rows, test_cols, name in test_cases:
        test_result = hodges_with_deletion(twistor, test_rows, test_cols, x_ref, y_ref)
        test_H = test_result[0]
        
        if test_H is None:
            continue
        
        ratio = test_H / baseline_H if baseline_H != 0 else None
        
        # Compute c ratios
        def get_c_factor(legs):
            if len(legs) != 3:
                return None
            i, j, k = legs
            aij = twistor.get_angle(i, j)
            ajk = twistor.get_angle(j, k)
            aki = twistor.get_angle(k, i)
            if aij == 0 or ajk == 0 or aki == 0:
                return None
            return QQ(1) / (aij * ajk * aki)
        
        c_baseline_rows = get_c_factor(baseline_rows)
        c_baseline_cols = get_c_factor(baseline_cols)
        c_test_rows = get_c_factor(test_rows)
        c_test_cols = get_c_factor(test_cols)
        
        if c_baseline_rows is None or c_baseline_cols is None or c_test_rows is None or c_test_cols is None:
            continue
        
        c_ratio = (c_test_rows * c_test_cols) / (c_baseline_rows * c_baseline_cols)
        
        # Compute sigma ratios
        sigma_baseline = hodges_sigma(baseline_rows, baseline_cols, 6)
        sigma_test = hodges_sigma(test_rows, test_cols, 6)
        sigma_ratio = sigma_test / sigma_baseline
        
        # Combined compensation
        compensation = sigma_ratio * c_ratio
        
        # Check if ratio matches compensation
        match = (ratio == compensation) if ratio is not None else False
        
        print(f"{name}: rows={test_rows}, cols={test_cols}")
        print(f"  H = {test_H}")
        print(f"  ratio = {ratio}")
        print(f"  sigma_ratio = {sigma_ratio}")
        print(f"  c_ratio = {c_ratio}")
        print(f"  compensation (sigma * c) = {compensation}")
        print(f"  Match: {match}")
        print()
        
        if not match:
            all_match = False
            # Check what's missing
            missing = ratio / compensation if compensation != 0 else None
            print(f"  Missing factor = {missing}")
            print()
    
    print(f"{'='*70}")
    if all_match:
        print("*** SUCCESS: All ratios match sigma * c compensation!")
        print("The compensation factors are correct!")
    else:
        print("*** PARTIAL: Some ratios don't match")
        print("Need to identify additional factor")
    print(f"{'='*70}")

if __name__ == '__main__':
    verify_compensation()






