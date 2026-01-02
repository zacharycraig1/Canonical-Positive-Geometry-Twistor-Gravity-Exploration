#!/usr/bin/env sage
# =============================================================================
# TEST: Use fixed reference legs (0,5) and only test deletions that exclude them
# =============================================================================

from sage.all import *
import sys
import os
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('tools/diagnostic_probe_1_hodges_invariance.sage')

def test_fixed_refs():
    """Test with fixed reference legs (0,5) and deletions that exclude them."""
    print("="*70)
    print("TESTING: Fixed Reference Legs (0,5) - Only Valid Deletions")
    print("="*70)
    print("Strategy: x=0, y=5 fixed. Only test deletions that exclude {0,5}")
    print("="*70)
    
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    x_ref = 0
    y_ref = 5
    
    # Only consider deletions that exclude 0 and 5
    all_indices = list(range(6))
    valid_triples = [list(c) for c in combinations(all_indices, 3) if 0 not in c and 5 not in c]
    
    print(f"\nValid deletion triples (excluding 0,5): {len(valid_triples)}")
    print(f"Examples: {valid_triples[:3]}")
    
    if len(valid_triples) < 2:
        print("ERROR: Not enough valid deletions")
        return
    
    # Baseline
    baseline_rows = valid_triples[0]
    baseline_cols = valid_triples[0]
    baseline_result = hodges_with_deletion(twistor, baseline_rows, baseline_cols, x_ref, y_ref)
    baseline_H = baseline_result[0]
    
    if baseline_H is None:
        print("ERROR: Baseline failed")
        return
    
    print(f"\nBaseline: rows={baseline_rows}, cols={baseline_cols}, H={baseline_H}")
    
    # Test other valid deletions
    print(f"\nTesting other valid deletions:")
    results = []
    for i, rows_del in enumerate(valid_triples[:5]):  # Test first 5
        for j, cols_del in enumerate(valid_triples[:5]):
            if i == 0 and j == 0:
                continue  # Skip baseline
            
            result = hodges_with_deletion(twistor, rows_del, cols_del, x_ref, y_ref)
            H_val = result[0]
            
            if H_val is None:
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
            
            status = "✓" if match else "✗"
            print(f"  {status} rows={rows_del}, cols={cols_del}: ratio={ratio}")
    
    print(f"\n{'='*70}")
    if results:
        all_match = all(r['match'] for r in results)
        unique_H = list(set([r['H'] for r in results]))
        
        print(f"Total tested: {len(results)}")
        print(f"Unique H values: {len(unique_H)}")
        print(f"Matches: {sum(1 for r in results if r['match'])}/{len(results)}")
        
        if all_match:
            print(f"\n*** SUCCESS: All valid deletions give the same result!")
            return True
        else:
            print(f"\n*** Still have differences")
            print(f"This means the issue is NOT just about refs in deleted sets")
            return False
    print(f"{'='*70}")

if __name__ == '__main__':
    test_fixed_refs()










