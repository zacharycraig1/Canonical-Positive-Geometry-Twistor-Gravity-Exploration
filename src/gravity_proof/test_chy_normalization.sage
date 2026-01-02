#!/usr/bin/env sage
# =============================================================================
# Test CHY Normalization Fix
# =============================================================================
# Tests that adding the <12>^8 helicity factor fixes the CHY vs Hodges discrepancy.

from sage.all import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

load("src/gravity_proof/chamber_sum_analysis.sage")

def test_normalization_fix(seed=42):
    """Test that CHY with helicity factor matches Hodges."""
    print("="*60)
    print("TESTING CHY NORMALIZATION FIX")
    print("="*60)
    print(f"Seed: {seed}")
    
    analysis = ChamberSumAnalysis(seed=seed)
    analysis.solve_and_group_by_chamber()
    result = analysis.compute_per_chamber_contribution()
    hodges = analysis.compute_hodges_amplitude()
    
    print("")
    print("="*60)
    print("COMPARISON")
    print("="*60)
    
    if result and hodges.get('status') == 'success':
        chy_total = result['total_float']
        hodges_val = hodges['hodges_float']
        ratio = chy_total / hodges_val if hodges_val != 0 else float('inf')
        rel_diff = abs(chy_total - hodges_val) / abs(hodges_val) if hodges_val != 0 else float('inf')
        
        print(f"CHY (with helicity): {chy_total:.6e}")
        print(f"Hodges:              {hodges_val:.6e}")
        print(f"Ratio CHY/Hodges:    {ratio:.6f}")
        print(f"Relative diff:       {rel_diff:.6e}")
        
        if rel_diff < 1e-8:
            print("")
            print("SUCCESS: CHY matches Hodges!")
            return True
        else:
            print("")
            print("Still have discrepancy")
            return False
    else:
        print("ERROR: Computation failed")
        return False

if __name__ == "__main__":
    # Test multiple seeds
    seeds = [42, 123, 456]
    results = []
    
    for seed in seeds:
        print("")
        print("#"*60)
        result = test_normalization_fix(seed)
        results.append(result)
        print("#"*60)
    
    print("")
    print("="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Passed: {sum(results)}/{len(results)}")

