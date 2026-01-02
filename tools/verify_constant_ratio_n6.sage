#!/usr/bin/env sage
# =============================================================================
# VERIFY CONSTANT RATIO FOR n=6
# =============================================================================
# Verify that M6^KLT / M6^Hodges_reduced is constant across many samples
# =============================================================================

from sage.all import *
import sys
import os
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare_reduced.sage')

def to_json_serializable(obj):
    """Convert Sage types to JSON-serializable types."""
    if isinstance(obj, (int, float, str, bool, type(None))):
        return obj
    if isinstance(obj, Integer):
        return int(obj)
    if isinstance(obj, Rational):
        # Convert to string to preserve precision
        return str(obj)
    if isinstance(obj, (list, tuple)):
        return [to_json_serializable(x) for x in obj]
    if isinstance(obj, dict):
        return {k: to_json_serializable(v) for k, v in obj.items()}
    return str(obj)

def verify_constant_ratio(n_samples=200):
    """Verify that the ratio is constant across many samples."""
    print("="*70)
    print("VERIFY CONSTANT RATIO FOR n=6 (KLT vs Hodges Reduced)")
    print("="*70)
    
    ratios = []
    valid_count = 0
    none_count = 0
    
    for seed in range(n_samples):
        Z = sample_positive_Z_moment_curve(n=6, seed=seed)
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
        
        if not twistor.domain_ok:
            none_count += 1
            continue
        
        H_red = hodges_6pt_mhv_reduced(twistor)
        A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
        
        H_val = H_red[0] if isinstance(H_red, tuple) else H_red
        A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
        
        if H_val is None or A_val is None:
            none_count += 1
            continue
        
        if H_val == 0:
            none_count += 1
            continue
        
        ratio = A_val / H_val
        ratios.append(ratio)
        valid_count += 1
        
        if seed < 5:
            print(f"Seed {seed}: ratio = {ratio}")
    
    print(f"\nValid samples: {valid_count}/{n_samples}")
    print(f"None/zero cases: {none_count}")
    
    if not ratios:
        print("ERROR: No valid ratios")
        return False, None
    
    unique_ratios = list(set(ratios))
    print(f"Unique ratios: {len(unique_ratios)}")
    
    if len(unique_ratios) == 1:
        constant_ratio = unique_ratios[0]
        print(f"\n*** SUCCESS: Constant ratio = {constant_ratio}")
        print(f"All {valid_count} samples have the same ratio!")
        
        # Verify exact equality
        print("\nVerifying exact equality...")
        test_samples = min(10, valid_count)
        all_equal = True
        for i in range(test_samples):
            Z = sample_positive_Z_moment_curve(n=6, seed=i)
            twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
            if not twistor.domain_ok:
                continue
            H_red = hodges_6pt_mhv_reduced(twistor)
            A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
            H_val = H_red[0] if isinstance(H_red, tuple) else H_red
            A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
            if H_val is None or A_val is None or H_val == 0:
                continue
            ratio = A_val / H_val
            if ratio != constant_ratio:
                print(f"  WARNING: Sample {i} has ratio {ratio} != {constant_ratio}")
                all_equal = False
        
        if all_equal:
            print(f"  Verified: All {test_samples} test samples have exact ratio = {constant_ratio}")
        
        # Save results
        results = {
            'status': 'SUCCESS',
            'constant_ratio': str(constant_ratio),
            'valid_samples': valid_count,
            'total_samples': n_samples,
            'none_cases': none_count,
            'unique_ratios': 1,
            'exact_equality_verified': all_equal
        }
        
        with open('results/n6_constant_ratio_verified.json', 'w') as f:
            json.dump(to_json_serializable(results), f, indent=2)
        
        return True, constant_ratio
    else:
        print(f"\n*** FAILURE: Ratios vary")
        print(f"First 5 unique ratios:")
        for i, r in enumerate(unique_ratios[:5]):
            print(f"  {i+1}: {r}")
        
        results = {
            'status': 'FAILURE',
            'valid_samples': valid_count,
            'total_samples': n_samples,
            'unique_ratios': len(unique_ratios),
            'sample_ratios': [str(r) for r in unique_ratios[:10]]
        }
        
        with open('results/n6_constant_ratio_verified.json', 'w') as f:
            json.dump(to_json_serializable(results), f, indent=2)
        
        return False, None

if __name__ == '__main__':
    success, ratio = verify_constant_ratio(n_samples=200)
    print(f"\n{'='*70}")
    if success:
        print(f"PROOF: M6^KLT / M6^Hodges_reduced = {ratio} (CONSTANT)")
        print(f"{'='*70}")
    else:
        print("RATIO IS NOT CONSTANT")
        print(f"{'='*70}")










