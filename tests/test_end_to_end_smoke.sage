#!/usr/bin/env sage
# =============================================================================
# TEST: End-to-End Smoke Test
# =============================================================================
# Quick verification that all components work together
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load modules
load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare.sage')


def test_end_to_end():
    """End-to-end test: sample, compute, compare."""
    print("="*70)
    print("TEST: End-to-End Smoke")
    print("="*70)
    
    n_tests = 10
    passed = 0
    failed = 0
    
    for seed in range(n_tests):
        # Sample
        Z = sample_positive_Z_moment_curve(n=6, seed=seed)
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
        
        if not twistor.domain_ok:
            print(f"Seed {seed}: Domain violation - {twistor.domain_reason}")
            continue
        
        # Compute Hodges
        H_result = hodges_6pt_mhv(twistor)
        H = H_result[0] if isinstance(H_result, tuple) else H_result
        H_reason = H_result[1] if isinstance(H_result, tuple) and len(H_result) > 1 else "ok"
        
        if H is None:
            print(f"Seed {seed}: Hodges None - {H_reason}")
            failed += 1
            continue
        
        # Compute KLT
        A_result = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
        A = A_result[0] if isinstance(A_result, tuple) else A_result
        A_reason = A_result[1] if isinstance(A_result, tuple) and len(A_result) > 1 else "ok"
        
        if A is None:
            print(f"Seed {seed}: KLT None - {A_reason}")
            failed += 1
            continue
        
        # Compare
        is_equal, ratio, diff = exact_equality_test(H, A)
        
        if is_equal:
            passed += 1
        else:
            print(f"Seed {seed}: Mismatch - ratio={ratio}")
            failed += 1
    
    print(f"\nResults: {passed} passed, {failed} failed out of {n_tests}")
    success = failed == 0 and passed == n_tests
    print(f"TEST: {'PASSED' if success else 'FAILED'}")
    
    return success


if __name__ == '__main__':
    test_end_to_end()






