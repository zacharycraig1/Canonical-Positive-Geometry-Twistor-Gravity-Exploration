# test_sign_in_positive_region.sage
"""
KEY HYPOTHESIS: In the "moment curve" region where all w_ij = [ij]/⟨ij⟩ > 0,
the gravity amplitude has a definite sign.

This would be a crucial step toward identifying the positive geometry.

For the forest polynomial expansion:
    det(L̃^R) = Σ_F ∏_{(i,j) ∈ E(F)} (-w_ij C_i C_j)

If all w_ij > 0 and C_i factors are positive, each term has sign (-1)^|E(F)|.

For n=6 with 3 roots, each forest has exactly n-3 = 3 edges.
So each term is (-1)^3 × (positive) = negative.
Therefore det(L̃^R) < 0 in the positive region.

This script tests:
1. Can we sample kinematics where all w_ij > 0?
2. Does the amplitude have a consistent sign there?
"""

from sage.all import *
from itertools import combinations
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical

print("="*70)
print("TESTING SIGN STRUCTURE IN POSITIVE REGION")
print("="*70)

def angle_bracket(lambdas, i, j):
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]

def square_bracket(lambdas_tilde, i, j):
    return lambdas_tilde[i][0] * lambdas_tilde[j][1] - lambdas_tilde[i][1] * lambdas_tilde[j][0]

def compute_weights(lambdas, lambdas_tilde):
    """Compute all w_ij = [ij]/⟨ij⟩"""
    n = len(lambdas)
    weights = {}
    for i in range(n):
        for j in range(i+1, n):
            ab = angle_bracket(lambdas, i, j)
            sb = square_bracket(lambdas_tilde, i, j)
            if ab != 0:
                weights[(i,j)] = sb / ab
            else:
                weights[(i,j)] = None  # undefined
    return weights

def all_weights_positive(weights):
    """Check if all w_ij > 0"""
    for (i,j), w in weights.items():
        if w is None or w <= 0:
            return False
    return True

def count_positive_weights(weights):
    """Count how many w_ij > 0"""
    count = 0
    for (i,j), w in weights.items():
        if w is not None and w > 0:
            count += 1
    return count

# Try to find kinematics in the positive region
print("\n[1] Searching for kinematics with all w_ij > 0...")

positive_samples = []
amplitude_signs = []

for seed in range(1, 1001):
    try:
        kin = SpinorKinematics.random_rational(n=6, seed=seed)
        lambdas = kin.lambdas
        lambdas_tilde = kin.tilde_lambdas
        
        weights = compute_weights(lambdas, lambdas_tilde)
        pos_count = count_positive_weights(weights)
        
        if pos_count >= 10:  # At least 10 out of 15 positive
            M_hodges, status = hodges_npt_mhv_canonical(lambdas, lambdas_tilde, (0, 1))
            
            if M_hodges is not None:
                ab12 = angle_bracket(lambdas, 0, 1)
                M_reduced = M_hodges / (ab12**8)
                
                sign = "+" if M_reduced > 0 else "-" if M_reduced < 0 else "0"
                
                if pos_count == 15:  # All positive - true positive region
                    positive_samples.append({
                        'seed': seed,
                        'M_reduced': float(M_reduced),
                        'sign': sign,
                        'weights': {f"({i},{j})": float(w) for (i,j), w in weights.items() if w is not None}
                    })
                    amplitude_signs.append(sign)
                    print(f"  Seed {seed}: ALL POSITIVE, M_reduced = {float(M_reduced):.6e} ({sign})")
                    
    except Exception as e:
        continue

print(f"\n[2] Found {len(positive_samples)} samples with all w_ij > 0")

if len(positive_samples) > 0:
    print("\n[3] Sign analysis:")
    pos_count = amplitude_signs.count("+")
    neg_count = amplitude_signs.count("-")
    zero_count = amplitude_signs.count("0")
    
    print(f"  Positive: {pos_count}")
    print(f"  Negative: {neg_count}")
    print(f"  Zero: {zero_count}")
    
    if pos_count > 0 and neg_count == 0:
        print(f"\n✓ DEFINITE SIGN: Amplitude is POSITIVE in the positive region")
    elif neg_count > 0 and pos_count == 0:
        print(f"\n✓ DEFINITE SIGN: Amplitude is NEGATIVE in the positive region")
    else:
        print(f"\n✗ MIXED SIGNS: No definite sign in the positive region")
else:
    print("\n[Alternative] Testing with partial positivity...")
    
    # Look for the best we can find
    best_pos_count = 0
    best_seed = None
    
    for seed in range(1, 501):
        try:
            kin = SpinorKinematics.random_rational(n=6, seed=seed)
            lambdas = kin.lambdas
            lambdas_tilde = kin.tilde_lambdas
            
            weights = compute_weights(lambdas, lambdas_tilde)
            pos_count = count_positive_weights(weights)
            
            if pos_count > best_pos_count:
                best_pos_count = pos_count
                best_seed = seed
                
        except:
            continue
    
    print(f"  Best found: seed={best_seed} with {best_pos_count}/15 positive weights")
    
    # Analyze this best case
    kin = SpinorKinematics.random_rational(n=6, seed=best_seed)
    lambdas = kin.lambdas
    lambdas_tilde = kin.tilde_lambdas
    weights = compute_weights(lambdas, lambdas_tilde)
    
    print("\n  Weight signs:")
    for (i,j), w in sorted(weights.items()):
        if w is not None:
            sign = "+" if w > 0 else "-"
            print(f"    w_{i}{j} = {float(w):.4f} ({sign})")

print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

if len(positive_samples) == 0:
    print("""
The "positive region" (all w_ij > 0) appears to be rare or empty for 
random rational kinematics. This suggests:

1. The positive region is measure-zero (lower-dimensional)
2. Special kinematics are needed to access it
3. The "moment curve" condition is very restrictive

This is consistent with the Phase V-X findings that the positive 
geometry approach faces fundamental obstacles.

For gravity, unlike Yang-Mills, the positive region may not be a 
simple convex polytope in kinematic space.
""")
else:
    if len(set(amplitude_signs)) == 1:
        print(f"""
✅ SIGNIFICANT FINDING:
In the positive region (all w_ij > 0), the amplitude has 
DEFINITE SIGN: {amplitude_signs[0]}

This supports the positive geometry interpretation:
- The amplitude could be the canonical form of a positive geometry
- The sign is consistent with the forest polynomial structure

Next step: Verify this holds for more samples and understand 
the boundary structure.
""")
    else:
        print("""
Mixed signs found in the positive region. This complicates the 
positive geometry interpretation.
""")

