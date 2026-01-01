#!/usr/bin/env sage
# =============================================================================
# DERIVE MISSING COMPENSATION FACTOR EMPIRICALLY
# =============================================================================
# Compare two deletions that should give the same result
# Solve for the missing factor
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('tools/diagnostic_probe_1_hodges_invariance.sage')

def derive_factor():
    """Derive the missing factor by comparing two deletions."""
    print("="*70)
    print("DERIVING MISSING COMPENSATION FACTOR")
    print("="*70)
    
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    x_ref = 0
    y_ref = 5
    
    # Baseline deletion
    baseline_rows = [0, 1, 2]
    baseline_cols = [3, 4, 5]
    baseline_result = hodges_with_deletion(twistor, baseline_rows, baseline_cols, x_ref, y_ref)
    baseline_H = baseline_result[0]
    
    # Test deletion
    test_rows = [0, 1, 3]
    test_cols = [2, 4, 5]
    test_result = hodges_with_deletion(twistor, test_rows, test_cols, x_ref, y_ref)
    test_H = test_result[0]
    
    print(f"\nBaseline: rows={baseline_rows}, cols={baseline_cols}")
    print(f"  H = {baseline_H}")
    print(f"\nTest: rows={test_rows}, cols={test_cols}")
    print(f"  H = {test_H}")
    
    ratio = test_H / baseline_H if baseline_H != 0 else None
    print(f"\nRatio = {ratio}")
    
    # Try to identify what factor would make them equal
    # If test_H = baseline_H * missing_factor, then missing_factor = ratio
    # But we want to see if ratio relates to angle brackets
    
    print(f"\nAnalyzing what factor relates these deletions...")
    
    # Compute c factors for both
    # Baseline: c_{012} and c_{345}
    c_baseline_rows = QQ(1) / (twistor.get_angle(0,1) * twistor.get_angle(1,2) * twistor.get_angle(2,0))
    c_baseline_cols = QQ(1) / (twistor.get_angle(3,4) * twistor.get_angle(4,5) * twistor.get_angle(5,3))
    
    # Test: c_{013} and c_{245}
    c_test_rows = QQ(1) / (twistor.get_angle(0,1) * twistor.get_angle(1,3) * twistor.get_angle(3,0))
    c_test_cols = QQ(1) / (twistor.get_angle(2,4) * twistor.get_angle(4,5) * twistor.get_angle(5,2))
    
    c_ratio = (c_test_rows * c_test_cols) / (c_baseline_rows * c_baseline_cols)
    
    print(f"  c_baseline = {c_baseline_rows} * {c_baseline_cols}")
    print(f"  c_test = {c_test_rows} * {c_test_cols}")
    print(f"  c_ratio = {c_ratio}")
    print(f"  ratio / c_ratio = {ratio / c_ratio if c_ratio != 0 else 'N/A'}")
    
    # Check sigma
    load('src/hodges_sigma.sage')
    sigma_baseline = hodges_sigma(baseline_rows, baseline_cols, 6)
    sigma_test = hodges_sigma(test_rows, test_cols, 6)
    sigma_ratio = sigma_test / sigma_baseline
    
    print(f"\n  sigma_baseline = {sigma_baseline}")
    print(f"  sigma_test = {sigma_test}")
    print(f"  sigma_ratio = {sigma_ratio}")
    print(f"  ratio / sigma_ratio = {ratio / sigma_ratio if sigma_ratio != 0 else 'N/A'}")
    
    # Combined
    combined_ratio = sigma_ratio * c_ratio
    print(f"\n  sigma_ratio * c_ratio = {combined_ratio}")
    print(f"  ratio / (sigma * c_ratio) = {ratio / combined_ratio if combined_ratio != 0 else 'N/A'}")
    
    # Check if there's a factor related to reference legs
    xy_bracket = twistor.get_angle(x_ref, y_ref)
    print(f"\n  <x y> = <{x_ref} {y_ref}> = {xy_bracket}")
    print(f"  Maybe need factor like <x y>^k?")
    print(f"  ratio * <xy>^2 = {ratio * (xy_bracket ** 2) if xy_bracket != 0 else 'N/A'}")
    print(f"  ratio * <xy>^4 = {ratio * (xy_bracket ** 4) if xy_bracket != 0 else 'N/A'}")

if __name__ == '__main__':
    derive_factor()









