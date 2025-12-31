#!/usr/bin/env sage
# =============================================================================
# RATIO RECOGNIZER FOR n=6
# =============================================================================
# Goal: Find kinematic factor F(Z) such that:
#   M6^KLT / (F(Z) * M6^Hodges) = constant
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare_reduced.sage')

def compute_invariant_features(twistor):
    """Compute candidate invariant features for ratio fitting."""
    features = {}
    
    # Cyclic angle brackets
    cyclic_prod = QQ(1)
    for i in range(6):
        j = (i + 1) % 6
        ang = twistor.get_angle(i, j)
        if ang == 0:
            return None
        cyclic_prod *= ang
    features['cyclic_prod'] = cyclic_prod
    
    # Products of specific angle brackets
    features['ang_01'] = twistor.get_angle(0, 1)
    features['ang_12'] = twistor.get_angle(1, 2)
    features['ang_23'] = twistor.get_angle(2, 3)
    features['ang_34'] = twistor.get_angle(3, 4)
    features['ang_45'] = twistor.get_angle(4, 5)
    features['ang_50'] = twistor.get_angle(5, 0)
    
    # Some four-brackets
    features['four_0123'] = twistor.get_four_bracket(0, 1, 2, 3)
    features['four_1234'] = twistor.get_four_bracket(1, 2, 3, 4)
    features['four_2345'] = twistor.get_four_bracket(2, 3, 4, 5)
    
    return features

def test_factor_candidate(twistor, factor_func):
    """Test if a candidate factor makes the ratio constant."""
    # Functions already loaded via load() at top
    
    H_red = hodges_6pt_mhv_reduced(twistor)
    A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
    
    H_val = H_red[0] if isinstance(H_red, tuple) else H_red
    A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
    
    if H_val is None or A_val is None or H_val == 0:
        return None
    
    try:
        F = factor_func(twistor)
        if F is None or F == 0:
            return None
        normalized_ratio = A_val / (F * H_val)
        return normalized_ratio
    except:
        return None

def find_normalization_factor(n_samples=50):
    """Try to find a normalization factor that makes ratio constant."""
    print("="*70)
    print("RATIO RECOGNIZER FOR n=6")
    print("="*70)
    
    # Collect data
    data = []
    for seed in range(n_samples):
        Z = sample_positive_Z_moment_curve(seed=seed)
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
        if not twistor.domain_ok:
            continue
        
        features = compute_invariant_features(twistor)
        if features is None:
            continue
        
        H_red = hodges_6pt_mhv_reduced(twistor)
        A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
        
        H_val = H_red[0] if isinstance(H_red, tuple) else H_red
        A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
        
        if H_val is None or A_val is None or H_val == 0:
            continue
        
        ratio = A_val / H_val
        data.append({
            'seed': seed,
            'twistor': twistor,
            'features': features,
            'ratio': ratio,
            'H': H_val,
            'A': A_val
        })
    
    print(f"Collected {len(data)} valid samples")
    
    if len(data) < 5:
        print("ERROR: Not enough valid samples")
        return None
    
    # Analyze ratio pattern
    ratios = [d['ratio'] for d in data]
    unique_ratios = list(set(ratios))
    print(f"Unique ratios (raw): {len(unique_ratios)}")
    
    # Try candidate factors
    candidates = [
        ("1", lambda t: QQ(1)),
        ("cyclic_prod", lambda t: compute_invariant_features(t)['cyclic_prod'] if compute_invariant_features(t) else None),
        ("cyclic_prod^2", lambda t: (compute_invariant_features(t)['cyclic_prod'])**2 if compute_invariant_features(t) else None),
        ("1/cyclic_prod", lambda t: QQ(1) / compute_invariant_features(t)['cyclic_prod'] if compute_invariant_features(t) and compute_invariant_features(t)['cyclic_prod'] != 0 else None),
        ("ang_01*ang_12", lambda t: compute_invariant_features(t)['ang_01'] * compute_invariant_features(t)['ang_12'] if compute_invariant_features(t) else None),
    ]
    
    best_factor = None
    best_unique_count = len(unique_ratios)
    
    for factor_name, factor_func in candidates:
        normalized_ratios = []
        for d in data:
            norm_ratio = test_factor_candidate(d['twistor'], factor_func)
            if norm_ratio is not None:
                normalized_ratios.append(norm_ratio)
        
        if normalized_ratios:
            unique_norm = list(set(normalized_ratios))
            print(f"\nFactor '{factor_name}': {len(unique_norm)} unique ratios from {len(normalized_ratios)} samples")
            if len(unique_norm) < best_unique_count:
                best_unique_count = len(unique_norm)
                best_factor = (factor_name, factor_func)
                if len(unique_norm) == 1:
                    print(f"  *** SUCCESS: Constant ratio = {unique_norm[0]}")
                    return factor_name, factor_func, unique_norm[0]
    
    if best_factor:
        print(f"\nBest factor: '{best_factor[0]}' reduces to {best_unique_count} unique ratios")
        return best_factor[0], best_factor[1], None
    else:
        print("\nNo factor found that reduces ratio variation")
        return None

if __name__ == '__main__':
    result = find_normalization_factor(n_samples=50)
    if result:
        factor_name, factor_func, constant = result
        print(f"\n{'='*70}")
        if constant:
            print(f"FOUND CONSTANT RATIO!")
            print(f"Factor: {factor_name}")
            print(f"Constant = {constant}")
        else:
            print(f"Best factor: {factor_name}")
        print(f"{'='*70}")

