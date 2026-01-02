#!/usr/bin/env sage
# =============================================================================
# VERIFY CANONICAL CONVENTION: Fixed KLT ordering + Fixed Hodges deletion
# =============================================================================
# After enforcing canonical conventions, verify we get exactly 1 constant ratio
# =============================================================================

from sage.all import *
import sys
import os
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')

def to_json_serializable(obj):
    """Convert Sage types to JSON-serializable types."""
    if isinstance(obj, (int, float, str, bool, type(None))):
        return obj
    if isinstance(obj, Integer):
        return int(obj)
    if isinstance(obj, Rational):
        return str(obj)
    if isinstance(obj, (list, tuple)):
        return [to_json_serializable(x) for x in obj]
    if isinstance(obj, dict):
        return {k: to_json_serializable(v) for k, v in obj.items()}
    return str(obj)

def verify_canonical(n_samples=200):
    """Verify constant ratio with canonical conventions enforced."""
    print("="*70)
    print("VERIFYING CANONICAL CONVENTION")
    print("="*70)
    print("\nCanonical conventions:")
    print("  KLT: Lexicographically sorted permutations")
    print("  Hodges: Fixed deletion rows=[0,1,2], cols=[3,4,5]")
    print("  Reference legs: x=0, y=5")
    print("="*70)
    
    ratios = []
    seed_to_ratio = {}  # Fix: use dict to track seed->ratio mapping
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
        
        if H_val is None or A_val is None or H_val == 0:
            none_count += 1
            continue
        
        ratio = A_val / H_val
        ratios.append(ratio)
        seed_to_ratio[seed] = ratio  # Store mapping
        valid_count += 1
        
        if seed < 5:
            print(f"Seed {seed}: ratio = {ratio}")
    
    print(f"\nValid samples: {valid_count}/{n_samples}")
    print(f"None/zero cases: {none_count}")
    
    if not ratios:
        print("ERROR: No valid ratios")
        return False, None, {}
    
    unique_ratios = list(set(ratios))
    print(f"Unique ratios: {len(unique_ratios)}")
    
    if len(unique_ratios) == 1:
        constant_ratio = unique_ratios[0]
        print(f"\n*** SUCCESS: Constant ratio = {constant_ratio}")
        print(f"All {valid_count} samples have the same ratio!")
        
        results = {
            'status': 'SUCCESS',
            'constant_ratio': str(constant_ratio),
            'valid_samples': valid_count,
            'total_samples': n_samples,
            'unique_ratios': 1
        }
        
        with open('results/canonical_verified.json', 'w') as f:
            json.dump(to_json_serializable(results), f, indent=2)
        
        return True, constant_ratio, seed_to_ratio
    else:
        print(f"\n*** Still have {len(unique_ratios)} unique ratios")
        print(f"First 5 unique ratios:")
        for i, r in enumerate(unique_ratios[:5]):
            print(f"  {i+1}: {r}")
        
        # Group by ratio to analyze
        ratio_to_seeds = {}
        for seed, ratio in seed_to_ratio.items():
            ratio_str = str(ratio)
            if ratio_str not in ratio_to_seeds:
                ratio_to_seeds[ratio_str] = []
            ratio_to_seeds[ratio_str].append(seed)
        
        print(f"\nRatio class distribution:")
        for i, (ratio_str, seeds) in enumerate(ratio_to_seeds.items()):
            print(f"  Class {i+1}: {len(seeds)} samples, seeds: {seeds[:5]}...")
        
        results = {
            'status': 'PARTIAL',
            'valid_samples': valid_count,
            'total_samples': n_samples,
            'unique_ratios': len(unique_ratios),
            'ratio_classes': {str(k): len(v) for k, v in ratio_to_seeds.items()}
        }
        
        with open('results/canonical_verified.json', 'w') as f:
            json.dump(to_json_serializable(results), f, indent=2)
        
        return False, None, seed_to_ratio

if __name__ == '__main__':
    success, ratio, seed_map = verify_canonical(n_samples=200)
    print(f"\n{'='*70}")
    if success:
        print(f"PROOF: M6^KLT / M6^Hodges = {ratio} (CONSTANT)")
        print(f"Canonical convention verified!")
    else:
        print(f"Still have multiple ratio classes - need further investigation")
    print(f"{'='*70}")










