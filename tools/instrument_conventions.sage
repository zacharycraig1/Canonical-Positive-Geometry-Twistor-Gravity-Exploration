#!/usr/bin/env sage
# =============================================================================
# INSTRUMENT CONVENTIONS: Log all discrete choices for each sample
# =============================================================================
# Goal: Identify what causes the 6 ratio classes
# =============================================================================

from sage.all import *
import sys
import os
import json
from collections import defaultdict
from itertools import permutations

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

def gravity_6pt_mhv_klt_instrumented(twistor, mandelstam_func):
    """
    KLT with instrumentation to log permutation ordering.
    """
    permuted_set = [1, 2, 3]
    fixed_leg_1 = 0
    fixed_leg_5 = 4
    fixed_leg_6 = 5
    
    # CRITICAL: Use lexicographically sorted permutations for canonical ordering
    all_perms = sorted(list(permutations(permuted_set)))
    
    # Log the permutation order
    perm_order = [list(p) for p in all_perms]
    
    total = QQ(0)
    contributions = []
    
    for alpha in all_perms:
        alpha = list(alpha)
        order_alpha = [fixed_leg_5, fixed_leg_6] + alpha + [fixed_leg_1]
        A_alpha = parke_taylor_6pt_mhv(twistor, order_alpha)
        if A_alpha is None:
            continue
        
        for beta in all_perms:
            beta = list(beta)
            order_beta = [fixed_leg_1] + beta + [fixed_leg_5, fixed_leg_6]
            A_beta = parke_taylor_6pt_mhv(twistor, order_beta)
            if A_beta is None:
                continue
            
            S = klt_momentum_kernel_6pt(alpha, beta, twistor, mandelstam_func)
            if S is None:
                continue
            
            contrib = A_alpha * S * A_beta
            total += contrib
            contributions.append({
                'alpha': alpha,
                'beta': beta,
                'A_alpha': str(A_alpha),
                'A_beta': str(A_beta),
                'S': str(S),
                'contrib': str(contrib)
            })
    
    return (total, "ok", {
        'perm_order': perm_order,
        'fixed_legs': [fixed_leg_1, fixed_leg_5, fixed_leg_6],
        'contributions': contributions
    })

def instrument_sample(seed, n_samples=200):
    """Instrument a single sample to log all convention choices."""
    Z = sample_positive_Z_moment_curve(n=6, seed=seed)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    
    if not twistor.domain_ok:
        return None
    
    # KLT with instrumentation
    klt_result = gravity_6pt_mhv_klt_instrumented(twistor, mandelstam_invariant)
    A_klt = klt_result[0]
    klt_info = klt_result[2] if len(klt_result) > 2 else {}
    
    # Hodges reduced
    H_red = hodges_6pt_mhv_reduced(twistor)
    H_val = H_red[0] if isinstance(H_red, tuple) else H_red
    
    if A_klt is None or H_val is None or H_val == 0:
        return None
    
    ratio = A_klt / H_val
    
    # Log Hodges choices
    hodges_info = {
        'rows_deleted': [0, 1, 2],
        'cols_deleted': [3, 4, 5],
        'rows_kept': [3, 4, 5],
        'cols_kept': [0, 1, 2],
        'reference_legs': [0, 5]
    }
    
    return {
        'seed': seed,
        'ratio': str(ratio),
        'A_klt': str(A_klt),
        'H_hodges': str(H_val),
        'klt_info': klt_info,
        'hodges_info': hodges_info
    }

def analyze_convention_classes(n_samples=200):
    """Instrument samples and group by ratio to find convention pattern."""
    print("="*70)
    print("INSTRUMENTING CONVENTIONS TO IDENTIFY 6 RATIO CLASSES")
    print("="*70)
    
    samples = []
    for seed in range(n_samples):
        result = instrument_sample(seed)
        if result:
            samples.append(result)
        if (seed + 1) % 50 == 0:
            print(f"  Processed {seed + 1}/{n_samples} samples...")
    
    print(f"\nCollected {len(samples)} valid samples")
    
    # Group by ratio
    ratio_groups = defaultdict(list)
    for s in samples:
        ratio_groups[s['ratio']].append(s)
    
    print(f"Found {len(ratio_groups)} unique ratios")
    
    # Analyze each group
    analysis = {}
    for i, (ratio_str, group_samples) in enumerate(ratio_groups.items()):
        print(f"\n{'='*70}")
        print(f"RATIO CLASS {i+1}: {len(group_samples)} samples")
        print(f"Ratio value: {ratio_str[:60]}...")
        
        # Check KLT permutation ordering (should be same for all)
        first_klt_info = group_samples[0]['klt_info']
        perm_order = first_klt_info.get('perm_order', [])
        print(f"KLT permutation order: {perm_order}")
        
        # Check if all samples have same KLT info
        all_same_klt = all(
            s['klt_info'].get('perm_order', []) == perm_order
            for s in group_samples
        )
        print(f"All samples have same KLT perm order: {all_same_klt}")
        
        # Check Hodges info (should be same for all)
        first_hodges = group_samples[0]['hodges_info']
        print(f"Hodges deletion: rows={first_hodges['rows_deleted']}, cols={first_hodges['cols_deleted']}")
        
        # Check seed pattern
        seeds = [s['seed'] for s in group_samples]
        print(f"Seed range: {min(seeds)} to {max(seeds)}")
        print(f"Seed mod 6: {sorted(set(s % 6 for s in seeds))}")
        
        analysis[ratio_str] = {
            'count': len(group_samples),
            'seeds': seeds[:10],  # First 10
            'klt_perm_order': perm_order,
            'hodges_deletion': first_hodges
        }
    
    # Save analysis
    with open('results/convention_analysis.json', 'w') as f:
        json.dump(to_json_serializable(analysis), f, indent=2)
    
    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE")
    print(f"Results saved to results/convention_analysis.json")
    print(f"{'='*70}")
    
    return analysis

if __name__ == '__main__':
    analyze_convention_classes(n_samples=200)










