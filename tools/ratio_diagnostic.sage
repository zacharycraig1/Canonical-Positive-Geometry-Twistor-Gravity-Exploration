#!/usr/bin/env sage
# =============================================================================
# RATIO DIAGNOSTIC: Analyze and eliminate ratio variation
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare.sage')

def analyze_ratio_pattern(n_samples=50):
    """Analyze ratio patterns to identify convention mismatches."""
    print("="*70)
    print("RATIO PATTERN ANALYSIS")
    print("="*70)
    
    ratios = []
    samples = []
    
    for seed in range(n_samples):
        Z = sample_positive_Z_moment_curve(n=6, seed=seed)
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
        
        if not twistor.domain_ok:
            continue
        
        H_reduced = hodges_6pt_mhv_reduced(twistor)
        A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
        
        H = H_reduced[0] if isinstance(H_reduced, tuple) else H_reduced
        A = A_klt[0] if isinstance(A_klt, tuple) else A_klt
        
        if H is None or A is None or H == 0:
            continue
        
        ratio = A / H
        ratios.append(ratio)
        samples.append({
            'seed': seed,
            'ratio': ratio,
            'H': H,
            'A': A
        })
    
    print(f"Collected {len(ratios)} valid ratios")
    
    # Check if constant
    unique_ratios = list(set(ratios))
    print(f"Unique ratios: {len(unique_ratios)}")
    
    if len(unique_ratios) == 1:
        print(f"SUCCESS: Constant ratio = {unique_ratios[0]}")
        return True, unique_ratios[0]
    
    # Analyze pattern
    print("\nAnalyzing ratio variation...")
    
    # Check if ratios are related by simple factors
    # Try to see if ratio depends on seed pattern
    ratio_by_seed_mod = {}
    for i, s in enumerate(samples):
        mod = s['seed'] % 10
        if mod not in ratio_by_seed_mod:
            ratio_by_seed_mod[mod] = []
        ratio_by_seed_mod[mod].append(s['ratio'])
    
    print(f"Ratio groups by seed mod 10:")
    for mod in sorted(ratio_by_seed_mod.keys()):
        unique_in_group = len(set(ratio_by_seed_mod[mod]))
        print(f"  mod {mod}: {len(ratio_by_seed_mod[mod])} samples, {unique_in_group} unique ratios")
    
    # Try to see if ratio correlates with angle brackets
    if len(samples) >= 3:
        s0, s1, s2 = samples[0], samples[1], samples[2]
        print(f"\nSample ratios:")
        print(f"  Seed {s0['seed']}: {s0['ratio']}")
        print(f"  Seed {s1['seed']}: {s1['ratio']}")
        print(f"  Seed {s2['seed']}: {s2['ratio']}")
        
        # Check if ratios are proportional to some invariant
        Z0 = sample_positive_Z_moment_curve(n=6, seed=s0['seed'])
        tw0 = MomentumTwistor(n=6, Z=Z0, check_domain=True)
        
        # Try common invariants
        ang_01_0 = tw0.get_angle(0, 1)
        ang_12_0 = tw0.get_angle(1, 2)
        prod_cyclic_0 = QQ(1)
        for i in range(6):
            j = (i+1) % 6
            prod_cyclic_0 *= tw0.get_angle(i, j)
        
        print(f"\nFor seed {s0['seed']}:")
        print(f"  <01> = {ang_01_0}")
        print(f"  <12> = {ang_12_0}")
        print(f"  ∏<i,i+1> = {prod_cyclic_0}")
        print(f"  Ratio = {s0['ratio']}")
    
    return False, None


def test_alternative_conventions():
    """Test alternative convention choices to find constant ratio."""
    print("\n" + "="*70)
    print("TESTING ALTERNATIVE CONVENTIONS")
    print("="*70)
    
    # Test different Hodges minor choices
    print("\n1. Testing different Hodges minor selections...")
    
    seed = 42
    Z = sample_positive_Z_moment_curve(n=6, seed=seed)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    if not twistor.domain_ok:
        print("Domain violation")
        return
    
    # Standard: delete (0,1,2) rows and (3,4,5) cols
    H_std = hodges_6pt_mhv_reduced(twistor)
    H_std_val = H_std[0] if isinstance(H_std, tuple) else H_std
    
    # Alternative: delete (0,1,2) rows and (0,1,2) cols
    # (This would be symmetric deletion)
    print(f"Standard reduced: {H_std_val}")
    
    # Test KLT with different helicity assignments
    print("\n2. Testing KLT with different helicity...")
    A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
    A_klt_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
    print(f"KLT (neg helicity 0,1): {A_klt_val}")
    
    if H_std_val is not None and A_klt_val is not None and H_std_val != 0:
        ratio = A_klt_val / H_std_val
        print(f"Ratio: {ratio}")
        print(f"Ratio (float): {float(ratio):.10e}")


if __name__ == '__main__':
    success, constant = analyze_ratio_pattern(n_samples=50)
    if not success:
        test_alternative_conventions()










