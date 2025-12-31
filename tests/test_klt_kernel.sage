#!/usr/bin/env sage
# =============================================================================
# TEST: KLT Kernel Verification
# =============================================================================
# Verify the KLT momentum kernel implementation
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


def test_klt_kernel_symmetry():
    """Test properties of the KLT kernel."""
    print("="*70)
    print("TEST: KLT Kernel Properties")
    print("="*70)
    
    Z = sample_positive_Z_moment_curve(n=6, seed=42)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    if not twistor.domain_ok:
        print("Domain violation, cannot test")
        return False
    
    permuted_set = [1, 2, 3]
    from itertools import permutations
    all_perms = list(permutations(permuted_set))
    
    # Test: kernel exists for all permutation pairs
    kernel_values = {}
    for alpha in all_perms:
        for beta in all_perms:
            alpha_list = list(alpha)
            beta_list = list(beta)
            S = klt_momentum_kernel_6pt(alpha_list, beta_list, twistor, mandelstam_invariant)
            kernel_values[(alpha, beta)] = S
    
    none_count = sum(1 for v in kernel_values.values() if v is None)
    print(f"Kernel computations: {len(kernel_values)}")
    print(f"None values: {none_count}")
    
    # Test: identity permutation
    identity = tuple(permuted_set)
    S_id = kernel_values.get((identity, identity))
    print(f"S[identity|identity] = {S_id}")
    
    print(f"\nTEST: {'PASSED' if none_count == 0 else 'FAILED'}")
    return none_count == 0


def test_klt_sum_structure():
    """Test that KLT sum has correct structure (36 terms for n=6)."""
    print("="*70)
    print("TEST: KLT Sum Structure")
    print("="*70)
    
    Z = sample_positive_Z_moment_curve(n=6, seed=42)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    if not twistor.domain_ok:
        print("Domain violation, cannot test")
        return False
    
    # Count non-zero contributions
    from itertools import permutations
    permuted_set = [1, 2, 3]
    all_perms = list(permutations(permuted_set))
    
    contributions = []
    for alpha in all_perms:
        alpha = list(alpha)
        order_alpha = [4, 5] + alpha + [0]
        A_alpha = parke_taylor_6pt_mhv(twistor, order_alpha)
        
        for beta in all_perms:
            beta = list(beta)
            order_beta = [0] + beta + [4, 5]
            A_beta = parke_taylor_6pt_mhv(twistor, order_beta)
            
            S = klt_momentum_kernel_6pt(alpha, beta, twistor, mandelstam_invariant)
            
            if A_alpha is not None and A_beta is not None and S is not None:
                contributions.append(A_alpha * S * A_beta)
    
    print(f"Non-zero contributions: {len(contributions)}")
    print(f"Expected (6×6): 36")
    
    total = sum(contributions)
    print(f"Total KLT sum: {total}")
    
    passed = len(contributions) == 36
    print(f"\nTEST: {'PASSED' if passed else 'FAILED'}")
    return passed


if __name__ == '__main__':
    test_klt_kernel_symmetry()
    print()
    test_klt_sum_structure()






