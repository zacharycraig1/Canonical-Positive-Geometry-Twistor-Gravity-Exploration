#!/usr/bin/env sage
# =============================================================================
# TEST: Hodges with fixed reference legs (not in deleted sets)
# =============================================================================
# Strategy: Always use x,y that are NOT in deleted sets
# =============================================================================

from sage.all import *
import sys
import os
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('tools/diagnostic_probe_1_hodges_invariance.sage')

def hodges_with_smart_refs(twistor, rows_delete, cols_delete):
    """
    Compute Hodges with reference legs chosen to NOT be in deleted sets.
    """
    n = 6
    # Choose x,y that are NOT in deleted sets
    all_indices = list(range(n))
    available = [i for i in all_indices if i not in rows_delete and i not in cols_delete]
    
    if len(available) < 2:
        # Fallback: use any two that aren't both in deleted sets
        row_set = set(rows_delete)
        col_set = set(cols_delete)
        # Try to find two indices where at least one is not in both sets
        candidates = [i for i in all_indices if i not in (row_set & col_set)]
        if len(candidates) >= 2:
            x_ref, y_ref = candidates[0], candidates[1]
        else:
            # Last resort: use first two available
            x_ref, y_ref = 0, 5
            if x_ref in rows_delete and x_ref in cols_delete:
                x_ref = 1
            if y_ref in rows_delete and y_ref in cols_delete:
                y_ref = 4
    else:
        x_ref, y_ref = available[0], available[1]
    
    return hodges_with_deletion(twistor, rows_delete, cols_delete, x_ref, y_ref)

def test_smart_refs():
    """Test if choosing refs outside deleted sets fixes deletion independence."""
    print("="*70)
    print("TESTING: Smart Reference Leg Selection")
    print("="*70)
    print("Strategy: Always choose x,y that are NOT in deleted sets")
    print("="*70)
    
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    # Baseline
    baseline_rows = [0, 1, 2]
    baseline_cols = [3, 4, 5]
    baseline_result = hodges_with_smart_refs(twistor, baseline_rows, baseline_cols)
    baseline_H = baseline_result[0]
    
    if baseline_H is None:
        print("ERROR: Baseline failed")
        return
    
    print(f"\nBaseline: rows={baseline_rows}, cols={baseline_cols}")
    print(f"  H = {baseline_H}")
    
    # Test different deletions
    test_deletions = [
        ([0,1,2], [3,4,5]),  # Baseline
        ([0,1,2], [0,1,2]),  # Symmetric
        ([0,1,3], [2,4,5]),  # Mixed
        ([0,2,4], [1,3,5]),  # Mixed
        ([1,2,3], [1,2,3]),  # Other symmetric
        ([1,2,4], [2,3,4]),  # Mixed
    ]
    
    print(f"\nTesting different deletions (with smart ref selection):")
    results = []
    for rows_del, cols_del in test_deletions:
        result = hodges_with_smart_refs(twistor, list(rows_del), list(cols_del))
        H_val = result[0]
        info = result[1] if len(result) > 1 else {}
        
        if H_val is None:
            print(f"  rows={rows_del}, cols={cols_del}: FAILED ({info.get('error', 'unknown')})")
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
        unique_H = list(set([r['H'] for r in results]))
        
        print(f"Total valid: {len(results)}")
        print(f"Unique H values: {len(unique_H)}")
        
        if all_match:
            print(f"\n*** SUCCESS: All deletions give the same result!")
            print(f"This fixes the deletion independence issue!")
            return True
        else:
            print(f"\n*** PARTIAL: {sum(1 for r in results if r['match'])}/{len(results)} match")
            print(f"Still need additional compensation")
            return False
    print(f"{'='*70}")

if __name__ == '__main__':
    test_smart_refs()









