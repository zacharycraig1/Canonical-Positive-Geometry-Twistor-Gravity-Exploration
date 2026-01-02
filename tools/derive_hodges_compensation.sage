#!/usr/bin/env sage
# =============================================================================
# DERIVE HODGES COMPENSATION FACTOR EMPIRICALLY
# =============================================================================
# Compare different deletions to find the missing compensation factor
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')

# Use the diagnostic probe function
load('tools/diagnostic_probe_1_hodges_invariance.sage')

def derive_compensation():
    """Derive the compensation factor empirically."""
    print("="*70)
    print("DERIVING HODGES COMPENSATION FACTOR")
    print("="*70)
    
    # Use one sample
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    baseline_rows = [0, 1, 2]
    baseline_cols = [3, 4, 5]
    baseline_x = 0
    baseline_y = 5
    
    baseline_result = hodges_with_deletion(twistor, baseline_rows, baseline_cols, baseline_x, baseline_y)
    baseline_H = baseline_result[0]
    
    print(f"Baseline H = {baseline_H}")
    
    # Test a few different deletions
    test_cases = [
        ([0,1,2], [0,1,2], 0, 5),  # Symmetric deletion
        ([0,1,3], [2,4,5], 0, 5),  # Mixed deletion
    ]
    
    print(f"\nAnalyzing compensation factors...")
    for rows_del, cols_del, x_ref, y_ref in test_cases:
        result = hodges_with_deletion(twistor, list(rows_del), list(cols_del), x_ref, y_ref)
        H_val = result[0]
        info = result[1] if len(result) > 1 else {}
        
        if H_val is None or H_val == 0:
            continue
        
        ratio = H_val / baseline_H if baseline_H != 0 else None
        
        print(f"\nDeletion: rows={rows_del}, cols={cols_del}, ref=({x_ref},{y_ref})")
        print(f"  H = {H_val}")
        print(f"  H / baseline = {ratio}")
        
        # Try to factor the ratio
        if ratio is not None:
            print(f"  Ratio factors: {ratio.factor()}")
            
            # Check if ratio relates to angle brackets
            # For rows_del=[0,1,2], cols_del=[0,1,2]:
            # Maybe need factor like <0,1><1,2><2,0> / <3,4><4,5><5,3>?
            if rows_del == [0,1,2] and cols_del == [0,1,2]:
                # Try: ratio should be related to c factors
                c_baseline_rows = QQ(1) / (twistor.get_angle(0,1) * twistor.get_angle(1,2) * twistor.get_angle(2,0))
                c_baseline_cols = QQ(1) / (twistor.get_angle(3,4) * twistor.get_angle(4,5) * twistor.get_angle(5,3))
                c_test_rows = QQ(1) / (twistor.get_angle(0,1) * twistor.get_angle(1,2) * twistor.get_angle(2,0))
                c_test_cols = QQ(1) / (twistor.get_angle(0,1) * twistor.get_angle(1,2) * twistor.get_angle(2,0))
                
                c_ratio = (c_test_rows * c_test_cols) / (c_baseline_rows * c_baseline_cols)
                print(f"  c_ratio = {c_ratio}")
                print(f"  ratio / c_ratio = {ratio / c_ratio if c_ratio != 0 else 'N/A'}")

if __name__ == '__main__':
    derive_compensation()










