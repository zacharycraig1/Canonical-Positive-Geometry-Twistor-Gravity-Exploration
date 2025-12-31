#!/usr/bin/env sage
# =============================================================================
# ANALYZE SOURCE OF 6 RATIO VALUES
# =============================================================================
# Check if it's related to KLT fixed legs, moment-curve pattern, or something else
# =============================================================================

from sage.all import *
import sys
import os
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')

def analyze_ratio_source(n_samples=200):
    """Analyze what causes the 6 ratio values."""
    print("="*70)
    print("ANALYZING SOURCE OF 6 RATIO VALUES")
    print("="*70)
    
    # Collect data with various properties
    data = []
    
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
        
        # Collect properties
        # Check moment-curve t values
        t_vals = []
        for i in range(6):
            # Z_i = (1, t, t^2, t^3), so t = Z_i[1]
            t_vals.append(Z[i][1])
        
        # Check if t values have a pattern
        t_diffs = [t_vals[i+1] - t_vals[i] for i in range(5)]
        
        data.append({
            'seed': seed,
            'ratio': ratio,
            'seed_mod_6': seed % 6,
            'seed_mod_7': seed % 7,
            't_vals': t_vals,
            't_diffs': t_diffs,
            'H': H_val,
            'A': A_val
        })
    
    print(f"Collected {len(data)} valid samples")
    
    # Group by ratio
    ratio_groups = defaultdict(list)
    for d in data:
        ratio_str = str(d['ratio'])
        ratio_groups[ratio_str].append(d)
    
    print(f"\nFound {len(ratio_groups)} unique ratios")
    
    # Analyze each ratio group
    for i, (ratio_str, samples) in enumerate(ratio_groups.items()):
        print(f"\nRatio {i+1}: {len(samples)} samples")
        print(f"  Ratio value: {ratio_str[:60]}...")
        
        # Check seed mod patterns
        seed_mods_6 = [s['seed_mod_6'] for s in samples]
        mod6_counts = defaultdict(int)
        for m in seed_mods_6:
            mod6_counts[m] += 1
        print(f"  Seed mod 6 distribution: {dict(mod6_counts)}")
        
        # Check if t_diffs have a pattern
        if samples:
            first_t_diffs = samples[0]['t_diffs']
            print(f"  First sample t_diffs: {[float(d) for d in first_t_diffs[:3]]}...")
        
        # Check seed ranges
        seeds = [s['seed'] for s in samples]
        print(f"  Seed range: {min(seeds)} to {max(seeds)}")
    
    # Check if ratios are related by a simple factor
    ratios_list = sorted([QQ(r) for r in ratio_groups.keys()], key=lambda x: abs(x))
    print(f"\nAnalyzing ratio relationships...")
    print(f"Ratio values (sorted by magnitude):")
    for i, r in enumerate(ratios_list):
        print(f"  {i+1}: {r}")
    
    # Check if they're related by simple factors
    if len(ratios_list) >= 2:
        base_ratio = ratios_list[0]
        print(f"\nRatios relative to first:")
        for i, r in enumerate(ratios_list[1:], 1):
            factor = r / base_ratio
            print(f"  Ratio {i+1} / Ratio 1 = {factor}")
            # Check if factor is a simple rational
            if factor.denominator() < 100:
                print(f"    -> Simple factor: {factor}")

if __name__ == '__main__':
    analyze_ratio_source(n_samples=200)






