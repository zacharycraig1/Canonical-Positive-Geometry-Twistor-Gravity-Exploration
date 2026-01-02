#!/usr/bin/env sage
# =============================================================================
# DIAGNOSTIC PROBE 2: KLT Basis Invariance
# =============================================================================
# Test if KLT amplitude is invariant under different fixed leg choices
# =============================================================================

from sage.all import *
import sys
import os
from itertools import permutations, combinations

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')

def klt_with_fixed_legs(twistor, mandelstam_func, fixed_legs):
    """
    Compute KLT with specific fixed legs.
    
    fixed_legs: tuple of 3 fixed leg indices (0-based)
    Returns: (amplitude, info_dict)
    """
    n = 6
    if len(fixed_legs) != 3:
        return (None, {'error': 'bad_fixed_legs'})
    
    fixed_set = set(fixed_legs)
    permuted_set = [i for i in range(n) if i not in fixed_set]
    
    if len(permuted_set) != 3:
        return (None, {'error': 'bad_permuted_set'})
    
    # Canonical ordering: lexicographic
    all_perms = sorted(list(permutations(permuted_set)))
    
    total = QQ(0)
    contributions = []
    
    for alpha in all_perms:
        alpha = list(alpha)
        # A(fixed[2], fixed[3], alpha, fixed[1])
        # Standard KLT: A(n-1, n, alpha, 1) * S[alpha|beta] * A(1, beta, n-1, n)
        order_alpha = [fixed_legs[2], fixed_legs[0]] + alpha + [fixed_legs[1]]
        A_alpha = parke_taylor_6pt_mhv(twistor, order_alpha)
        if A_alpha is None:
            continue
        
        for beta in all_perms:
            beta = list(beta)
            # A(fixed[1], beta, fixed[2], fixed[3])
            order_beta = [fixed_legs[1]] + beta + [fixed_legs[2], fixed_legs[0]]
            A_beta = parke_taylor_6pt_mhv(twistor, order_beta)
            if A_beta is None:
                continue
            
            # KLT kernel: S[alpha|beta] with pivot = fixed_legs[1]
            # Need to adapt kernel for different pivot
            S = klt_momentum_kernel_6pt(alpha, beta, twistor, mandelstam_func)
            if S is None:
                continue
            
            contrib = A_alpha * S * A_beta
            total += contrib
            contributions.append({
                'alpha': alpha,
                'beta': beta,
                'contrib': str(contrib)
            })
    
    return (total, {
        'fixed_legs': fixed_legs,
        'permuted_set': permuted_set,
        'perm_order': [list(p) for p in all_perms],
        'num_contribs': len(contributions)
    })

def test_klt_basis_invariance():
    """Test KLT invariance under different fixed leg choices."""
    print("="*70)
    print("DIAGNOSTIC PROBE 2: KLT Basis Invariance")
    print("="*70)
    
    # Use one generic sample
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    if not twistor.domain_ok:
        print("ERROR: Sample not in domain")
        return
    
    # Baseline: standard fixed legs (0, 4, 5) for (1, 5, 6)
    baseline_fixed = (0, 4, 5)
    
    baseline_result = klt_with_fixed_legs(twistor, mandelstam_invariant, baseline_fixed)
    baseline_KLT = baseline_result[0]
    
    if baseline_KLT is None:
        print(f"ERROR: Baseline computation failed: {baseline_result[1]}")
        return
    
    print(f"\nBaseline:")
    print(f"  Fixed legs: {baseline_fixed}")
    print(f"  KLT = {baseline_KLT}")
    
    # Test different fixed leg choices
    print(f"\nTesting different fixed leg choices...")
    
    # All ways to choose 3 fixed legs from 6
    all_fixed_choices = list(combinations(range(6), 3))
    
    # Test a few key choices
    test_fixed = [
        (0, 4, 5),  # Baseline: (1, 5, 6)
        (0, 3, 5),  # (1, 4, 6)
        (0, 4, 3),  # (1, 5, 4)
        (1, 4, 5),  # (2, 5, 6)
        (0, 1, 5),  # (1, 2, 6)
    ]
    
    results = []
    for fixed_legs in test_fixed:
        result = klt_with_fixed_legs(twistor, mandelstam_invariant, fixed_legs)
        KLT_val = result[0]
        info = result[1] if len(result) > 1 else {}
        
        if KLT_val is None:
            print(f"  Failed: fixed={fixed_legs}, error={info.get('error', 'unknown')}")
            continue
        
        ratio = KLT_val / baseline_KLT if baseline_KLT != 0 else None
        results.append({
            'fixed_legs': fixed_legs,
            'KLT': KLT_val,
            'ratio': ratio,
            'info': info
        })
        
        print(f"  Fixed legs: {fixed_legs}")
        print(f"    KLT = {KLT_val}")
        print(f"    KLT / baseline = {ratio}")
    
    # Analyze results
    print(f"\n{'='*70}")
    print("ANALYSIS:")
    print(f"{'='*70}")
    
    if not results:
        print("ERROR: No valid results")
        return
    
    ratios = [r['ratio'] for r in results]
    unique_ratios = list(set(ratios))
    
    print(f"Total valid computations: {len(results)}")
    print(f"Unique KLT/baseline ratios: {len(unique_ratios)}")
    
    if len(unique_ratios) == 1 and unique_ratios[0] == 1:
        print(f"\n*** SUCCESS: KLT is basis invariant!")
        print(f"All ratios = 1 (exact)")
    else:
        print(f"\n*** FAILURE: KLT is NOT basis invariant")
        print(f"Unique ratios: {unique_ratios[:5]}")
        print(f"\nThis explains the 6 classes - KLT basis dependence!")
        
        # Check if ratios match the 6 classes
        print(f"\nChecking if ratios match known 6 classes...")
        known_ratios = [
            QQ(3172174131043478298326403, 2000000000000000000000000),
            QQ(10453863930907666325843479461658434627489, 17047398293643451000000000000000000000000),
        ]
        for r in results:
            print(f"  fixed={r['fixed_legs']}: ratio={r['ratio']}")

if __name__ == '__main__':
    test_klt_basis_invariance()










