#!/usr/bin/env sage
# =============================================================================
# TEST REFERENCE INDEPENDENCE
# =============================================================================
# Fix one deletion, vary reference legs, verify results match
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('tools/diagnostic_probe_1_hodges_invariance.sage')

def test_reference_independence():
    """Test if Hodges is reference-independent for fixed deletion."""
    print("="*70)
    print("TESTING REFERENCE INDEPENDENCE")
    print("="*70)
    
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    # Test with baseline deletion
    rows_del = [0, 1, 2]
    cols_del = [3, 4, 5]
    
    print(f"\nFixed deletion: rows={rows_del}, cols={cols_del}")
    print(f"Testing different reference leg pairs...")
    
    # Test different reference pairs
    ref_pairs = [
        (0, 5),  # Baseline
        (0, 4),
        (1, 5),
        (2, 5),
        (0, 3),
        (1, 4),
    ]
    
    results = []
    for x_ref, y_ref in ref_pairs:
        result = hodges_with_deletion(twistor, rows_del, cols_del, x_ref, y_ref)
        H_val = result[0]
        
        if H_val is None:
            print(f"  ref=({x_ref},{y_ref}): FAILED")
            continue
        
        results.append({
            'ref': (x_ref, y_ref),
            'H': H_val
        })
        print(f"  ref=({x_ref},{y_ref}): H = {H_val}")
    
    if not results:
        print("ERROR: No valid results")
        return
    
    # Check if all H values are the same
    H_values = [r['H'] for r in results]
    unique_H = list(set(H_values))
    
    print(f"\n{'='*70}")
    print("ANALYSIS:")
    print(f"{'='*70}")
    print(f"Total valid computations: {len(results)}")
    print(f"Unique H values: {len(unique_H)}")
    
    if len(unique_H) == 1:
        print(f"\n*** SUCCESS: Reference independent!")
        print(f"All reference pairs give H = {unique_H[0]}")
        return True
    else:
        print(f"\n*** FAILURE: Reference dependent")
        print(f"Unique H values:")
        for i, h in enumerate(unique_H[:5]):
            print(f"  {i+1}: {h}")
        
        # Check ratios
        baseline_H = H_values[0]
        print(f"\nRatios relative to baseline (ref={results[0]['ref']}):")
        for r in results:
            ratio = r['H'] / baseline_H if baseline_H != 0 else None
            print(f"  ref={r['ref']}: ratio = {ratio}")
        
        return False

def test_deletion_with_ref_in_deleted():
    """Test if issue occurs when reference legs are in deleted sets."""
    print(f"\n{'='*70}")
    print("TESTING: Reference legs in deleted sets")
    print(f"{'='*70}")
    
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    # Test case where reference leg is in deleted set
    # rows_del = [0,1,2], so if x=0, then x is deleted
    rows_del = [0, 1, 2]
    cols_del = [0, 1, 2]  # Symmetric deletion
    
    print(f"\nDeletion: rows={rows_del}, cols={cols_del}")
    print(f"Note: This deletion includes leg 0, which we'll use as reference")
    
    ref_pairs = [
        (0, 5),  # x=0 is in deleted rows
        (1, 5),  # x=1 is in deleted rows
        (2, 5),  # x=2 is in deleted rows
        (3, 5),  # x=3 is NOT in deleted rows
    ]
    
    results = []
    for x_ref, y_ref in ref_pairs:
        result = hodges_with_deletion(twistor, rows_del, cols_del, x_ref, y_ref)
        H_val = result[0]
        
        if H_val is None:
            print(f"  ref=({x_ref},{y_ref}): FAILED")
            continue
        
        results.append({
            'ref': (x_ref, y_ref),
            'x_in_deleted': x_ref in rows_del,
            'H': H_val
        })
        print(f"  ref=({x_ref},{y_ref}), x in deleted: {x_ref in rows_del}, H = {H_val}")
    
    if results:
        unique_H = list(set([r['H'] for r in results]))
        print(f"\nUnique H values: {len(unique_H)}")
        if len(unique_H) > 1:
            print(f"This deletion shows reference dependence!")

if __name__ == '__main__':
    success1 = test_reference_independence()
    test_deletion_with_ref_in_deleted()






