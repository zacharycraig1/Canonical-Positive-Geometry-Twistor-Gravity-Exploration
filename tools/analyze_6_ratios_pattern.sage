#!/usr/bin/env sage
# =============================================================================
# ANALYZE THE 6 UNIQUE RATIOS PATTERN
# =============================================================================

from sage.all import *
import sys
import os
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')

def analyze_ratio_pattern(n_samples=200):
    """Analyze the pattern of the 6 unique ratios."""
    print("="*70)
    print("ANALYZING 6 UNIQUE RATIOS PATTERN")
    print("="*70)
    
    ratio_to_seeds = defaultdict(list)
    ratio_to_indices = defaultdict(list)
    
    for seed in range(n_samples):
        Z = sample_positive_Z_moment_curve(n=6, seed=seed)
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
        ratio_str = str(ratio)
        ratio_to_seeds[ratio_str].append(seed)
        ratio_to_indices[ratio_str].append(seed % 6)  # Check if related to index mod 6
    
    print(f"\nFound {len(ratio_to_seeds)} unique ratios")
    
    # Analyze clustering
    for i, (ratio_str, seeds) in enumerate(ratio_to_seeds.items()):
        print(f"\nRatio {i+1}: {len(seeds)} samples")
        print(f"  Ratio value: {ratio_str[:80]}...")
        print(f"  Seed range: {min(seeds)} to {max(seeds)}")
        print(f"  Seed mod 6 pattern: {sorted(set(ratio_to_indices[ratio_str]))}")
        
        # Check if seeds cluster
        seeds_sorted = sorted(seeds)
        clusters = []
        current_cluster = [seeds_sorted[0]]
        for s in seeds_sorted[1:]:
            if s - current_cluster[-1] <= 10:  # Within 10 of previous
                current_cluster.append(s)
            else:
                clusters.append(current_cluster)
                current_cluster = [s]
        clusters.append(current_cluster)
        
        print(f"  Clusters: {len(clusters)} (largest: {max(len(c) for c in clusters)} samples)")
    
    # Try to find a pattern in the ratios themselves
    ratios_list = [QQ(r) for r in ratio_to_seeds.keys()]
    print(f"\nAnalyzing ratio relationships...")
    
    # Check if ratios are related by simple factors
    for i, r1 in enumerate(ratios_list):
        for j, r2 in enumerate(ratios_list):
            if i < j and r1 != 0 and r2 != 0:
                factor = r2 / r1
                # Check if factor is a simple rational
                if factor.denominator() < 1000:  # Simple denominator
                    print(f"  Ratio {i+1} / Ratio {j+1} = {factor}")
    
    return ratio_to_seeds

if __name__ == '__main__':
    analyze_ratio_pattern(n_samples=200)






