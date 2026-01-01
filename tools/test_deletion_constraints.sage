#!/usr/bin/env sage
# =============================================================================
# TEST: Should reference legs be excluded from deleted sets?
# =============================================================================
# Hypothesis: Reference legs (x,y) should NOT be in deleted sets
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('tools/diagnostic_probe_1_hodges_invariance.sage')

def test_deletion_constraints():
    """Test if enforcing x,y not in deleted sets fixes the issue."""
    print("="*70)
    print("TESTING: Enforce reference legs NOT in deleted sets")
    print("="*70)
    
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    # Baseline: x=0, y=5, so deleted sets should NOT include 0 or 5
    x_ref = 0
    y_ref = 5
    
    print(f"\nReference legs: x={x_ref}, y={y_ref}")
    print(f"Testing deletions that EXCLUDE {x_ref} and {y_ref} from deleted sets...")
    
    # All ways to choose 3 rows and 3 cols that don't include x or y
    from itertools import combinations
    all_indices = list(range(6))
    valid_rows = [list(c) for c in combinations(all_indices, 3) if x_ref not in c and y_ref not in c]
    valid_cols = [list(c) for c in combinations(all_indices, 3) if x_ref not in c and y_ref not in c]
    
    print(f"Valid row deletions (excluding {x_ref},{y_ref}): {len(valid_rows)}")
    print(f"Valid col deletions (excluding {x_ref},{y_ref}): {len(valid_cols)}")
    
    # Test a few combinations
    test_cases = [
        ([1,2,3], [1,2,3]),  # Symmetric, excludes 0,5
        ([1,2,3], [2,3,4]),  # Mixed, excludes 0,5
        ([1,2,4], [2,3,4]),  # Mixed, excludes 0,5
    ]
    
    baseline_rows = [1, 2, 3]
    baseline_cols = [1, 2, 3]
    baseline_result = hodges_with_deletion(twistor, baseline_rows, baseline_cols, x_ref, y_ref)
    baseline_H = baseline_result[0]
    
    if baseline_H is None:
        print("ERROR: Baseline failed")
        return
    
    print(f"\nBaseline: rows={baseline_rows}, cols={baseline_cols}, H={baseline_H}")
    
    results = []
    for rows_del, cols_del in test_cases:
        # Verify x,y not in deleted sets
        if x_ref in rows_del or y_ref in rows_del:
            print(f"  SKIP: rows={rows_del} includes reference leg")
            continue
        if x_ref in cols_del or y_ref in cols_del:
            print(f"  SKIP: cols={cols_del} includes reference leg")
            continue
        
        result = hodges_with_deletion(twistor, rows_del, cols_del, x_ref, y_ref)
        H_val = result[0]
        
        if H_val is None:
            print(f"  FAILED: rows={rows_del}, cols={cols_del}")
            continue
        
        ratio = H_val / baseline_H if baseline_H != 0 else None
        match = (ratio == 1) if ratio is not None else False
        
        results.append({
            'rows': rows_del,
            'cols': cols_del,
            'H': H_val,
            'ratio': ratio,
            'match': match
        })
        
        print(f"  rows={rows_del}, cols={cols_del}: H={H_val}, ratio={ratio}, match={match}")
    
    print(f"\n{'='*70}")
    if results:
        all_match = all(r['match'] for r in results)
        if all_match:
            print("*** SUCCESS: All deletions (excluding refs) give same result!")
        else:
            print("*** PARTIAL: Some deletions still differ")
            print("This suggests the issue is NOT just about refs in deleted sets")
    print(f"{'='*70}")

if __name__ == '__main__':
    test_deletion_constraints()









