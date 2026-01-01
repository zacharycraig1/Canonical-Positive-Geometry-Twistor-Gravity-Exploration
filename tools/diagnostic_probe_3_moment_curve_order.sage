#!/usr/bin/env sage
# =============================================================================
# DIAGNOSTIC PROBE 3: Moment-Curve Ordering Hypothesis
# =============================================================================
# Test if 6 ratio classes correlate with ordering of permuted legs by t_i
# =============================================================================

from sage.all import *
import sys
import os
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')

def get_moment_curve_params(Z):
    """Extract moment-curve parameters t_i from Z."""
    # Z_i = (1, t_i, t_i^2, t_i^3), so t_i = Z_i[1]
    return [Z[i][1] for i in range(len(Z))]

def argsort(values):
    """Return indices that would sort the values."""
    return sorted(range(len(values)), key=lambda i: values[i])

def test_moment_curve_order():
    """Test if ratio classes correlate with moment-curve ordering."""
    print("="*70)
    print("DIAGNOSTIC PROBE 3: Moment-Curve Ordering Hypothesis")
    print("="*70)
    print("\nHypothesis: 6 classes = 6 orderings of permuted legs {1,2,3} by t_i")
    print("="*70)
    
    # Permuted set in KLT: {1,2,3} (0-based)
    permuted_indices = [1, 2, 3]
    
    data = []
    
    for seed in range(200):
        Z = sample_positive_Z_moment_curve(n=6, seed=seed)
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
        
        if not twistor.domain_ok:
            continue
        
        # Get moment-curve parameters
        t_vals = get_moment_curve_params(Z)
        
        # Get t values for permuted legs
        t_permuted = [t_vals[i] for i in permuted_indices]
        
        # Compute ordering permutation
        # argsort gives indices in sorted order, we want the permutation that sorts
        sorted_indices = argsort(t_permuted)
        # Convert to permutation of {0,1,2} representing order of {1,2,3}
        perm = tuple(sorted_indices)  # This is the permutation that sorts
        
        # Also get the actual sorted order as a tuple of leg labels
        sorted_legs = tuple([permuted_indices[i] for i in sorted_indices])
        
        # Compute ratio
        H_red = hodges_6pt_mhv_reduced(twistor)
        A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
        
        H_val = H_red[0] if isinstance(H_red, tuple) else H_red
        A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
        
        if H_val is None or A_val is None or H_val == 0:
            continue
        
        ratio = A_val / H_val
        ratio_str = str(ratio)
        
        data.append({
            'seed': seed,
            't_vals': t_vals,
            't_permuted': t_permuted,
            'perm': perm,
            'sorted_legs': sorted_legs,
            'ratio': ratio,
            'ratio_str': ratio_str
        })
    
    print(f"\nCollected {len(data)} valid samples")
    
    # Group by ratio
    ratio_groups = defaultdict(list)
    for d in data:
        ratio_groups[d['ratio_str']].append(d)
    
    print(f"Found {len(ratio_groups)} unique ratios")
    
    # Analyze correlation
    print(f"\n{'='*70}")
    print("ANALYSIS: Ratio Class vs Moment-Curve Ordering")
    print(f"{'='*70}")
    
    for i, (ratio_str, samples) in enumerate(ratio_groups.items()):
        print(f"\nRatio Class {i+1}: {len(samples)} samples")
        print(f"  Ratio: {ratio_str[:50]}...")
        
        # Check ordering patterns
        perms = [s['perm'] for s in samples]
        sorted_legs_list = [s['sorted_legs'] for s in samples]
        
        unique_perms = list(set(perms))
        unique_sorted_legs = list(set(sorted_legs_list))
        
        print(f"  Unique perms (argsort indices): {unique_perms}")
        print(f"  Unique sorted leg orders: {unique_sorted_legs}")
        
        if len(unique_perms) == 1:
            print(f"  *** All samples have same ordering!")
        else:
            print(f"  Ordering varies: {len(unique_perms)} different orderings")
    
    # Check one-to-one mapping
    print(f"\n{'='*70}")
    print("ONE-TO-ONE MAPPING TEST:")
    print(f"{'='*70}")
    
    # For each ratio class, check if it has a unique ordering
    ratio_to_ordering = {}
    for ratio_str, samples in ratio_groups.items():
        perms = [s['perm'] for s in samples]
        unique_perms = list(set(perms))
        if len(unique_perms) == 1:
            ratio_to_ordering[ratio_str] = unique_perms[0]
        else:
            ratio_to_ordering[ratio_str] = None
    
    if all(v is not None for v in ratio_to_ordering.values()):
        orderings = list(ratio_to_ordering.values())
        if len(orderings) == len(set(orderings)):
            print(f"\n*** SUCCESS: One-to-one mapping ratio_class ↔ ordering!")
            print(f"Each of the {len(ratio_groups)} ratio classes maps to a unique ordering")
            print(f"\nThis explains the 6 classes - moment-curve ordering!")
        else:
            print(f"\nMultiple ratio classes map to same ordering")
    else:
        print(f"\nSome ratio classes have multiple orderings")
    
    # Show the mapping
    print(f"\nMapping:")
    for ratio_str, ordering in ratio_to_ordering.items():
        if ordering is not None:
            print(f"  Ratio class {ratio_str[:30]}... → ordering {ordering}")

if __name__ == '__main__':
    test_moment_curve_order()









