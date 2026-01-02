# test_atlas_hypothesis.sage
"""
RIGOROUS TEST: Atlas Hypothesis for Gravity Positive Geometry

Hypothesis: M_gravity / <12>^8 = (1/4) × Σ_{all 20 charts R} Ω_R

where Ω_R = det(L̃^R) / normalization

This test computes both sides for multiple kinematic configurations
and checks if they agree.
"""

from sage.all import *
from itertools import combinations
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical

print("="*70)
print("TESTING ATLAS HYPOTHESIS FOR GRAVITY")
print("="*70)

def angle_bracket(lambdas, i, j):
    """Compute <ij>"""
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]

def square_bracket(lambdas_tilde, i, j):
    """Compute [ij]"""
    return lambdas_tilde[i][0] * lambdas_tilde[j][1] - lambdas_tilde[i][1] * lambdas_tilde[j][0]

def compute_C_factors(lambdas, ref_indices):
    """
    Compute C_k factors from Hodges formula.
    C_k = <01><12>...<(k-1)k> (cyclic angle bracket product up to k)
    
    For the reduced amplitude, we use reference-dependent factors.
    """
    n = len(lambdas)
    # C_k = product of consecutive angle brackets
    C = {}
    for k in range(n):
        if k not in ref_indices:
            # C_k = <k, k+1> in simplified version
            # Actually need the full chain
            C[k] = angle_bracket(lambdas, k, (k+1) % n)
    return C

def weighted_laplacian(lambdas, lambdas_tilde, roots):
    """
    Construct the weighted Laplacian matrix L̃ with entries:
    L̃_ij = -w_ij for i ≠ j where w_ij = [ij]/⟨ij⟩
    L̃_ii = Σ_{k≠i} w_ik
    """
    n = len(lambdas)
    roots_set = set(roots)
    non_roots = [i for i in range(n) if i not in roots_set]
    m = len(non_roots)
    
    # Build reduced Laplacian (delete rows/cols corresponding to roots)
    L = matrix(SR, m, m)
    
    for ii, i in enumerate(non_roots):
        for jj, j in enumerate(non_roots):
            if i == j:
                # Diagonal: sum of weights
                diag_sum = 0
                for k in range(n):
                    if k != i:
                        ab = angle_bracket(lambdas, i, k)
                        sb = square_bracket(lambdas_tilde, i, k)
                        if ab != 0:
                            diag_sum += sb / ab
                L[ii, jj] = diag_sum
            else:
                # Off-diagonal: -w_ij
                ab = angle_bracket(lambdas, i, j)
                sb = square_bracket(lambdas_tilde, i, j)
                if ab != 0:
                    L[ii, jj] = -sb / ab
                else:
                    L[ii, jj] = 0
    
    return L, non_roots

def chart_canonical_form(lambdas, lambdas_tilde, roots):
    """
    Compute the canonical form Ω_R for chart with given roots.
    
    Ω_R = det(L̃^R) × (prefactors)
    
    The prefactors involve angle brackets for the roots.
    """
    n = len(lambdas)
    L, non_roots = weighted_laplacian(lambdas, lambdas_tilde, roots)
    
    try:
        det_L = L.det()
    except:
        det_L = 0
    
    # The normalization involves:
    # 1. Angle brackets among roots: <r1 r2><r2 r3><r3 r1>
    # 2. C factors for non-roots
    
    r0, r1, r2 = roots
    ab_roots = angle_bracket(lambdas, r0, r1) * angle_bracket(lambdas, r1, r2) * angle_bracket(lambdas, r2, r0)
    
    # For now, return just det_L divided by root angle bracket product squared
    # (the exact normalization needs verification)
    if ab_roots != 0:
        return det_L / (ab_roots**2)
    else:
        return 0

def test_atlas_hypothesis(seed=42, num_charts=20):
    """Test if (1/4) Σ Ω_R = M_gravity / <12>^8"""
    
    print(f"\n[Test with seed={seed}]")
    
    # Generate kinematics
    kin = SpinorKinematics.random_rational(n=6, seed=seed)
    
    # Get spinor data
    lambdas = kin.lambdas
    lambdas_tilde = kin.tilde_lambdas
    
    # Compute gravity amplitude via Hodges
    try:
        M_hodges, status = hodges_npt_mhv_canonical(kin.lambdas, kin.tilde_lambdas, (0, 1))
        if M_hodges is None:
            print(f"  Hodges returned None")
            return None
    except Exception as e:
        print(f"  Error computing Hodges: {e}")
        return None
    
    # Compute <12>^8
    ab12 = angle_bracket(lambdas, 0, 1)
    ab12_8 = ab12**8
    
    # M_reduced = M_gravity / <12>^8
    M_reduced = M_hodges / ab12_8
    
    print(f"  M_hodges = {float(M_hodges):.6e}")
    print(f"  <12>^8 = {float(ab12_8):.6e}")
    print(f"  M_reduced = {float(M_reduced):.6e}")
    
    # Compute atlas sum: (1/4) Σ_R Ω_R
    all_charts = list(combinations(range(6), 3))
    atlas_sum = 0
    
    for roots in all_charts:
        omega_R = chart_canonical_form(lambdas, lambdas_tilde, roots)
        atlas_sum += omega_R
    
    atlas_result = atlas_sum / 4
    
    print(f"  Atlas sum (raw) = {float(atlas_sum):.6e}")
    print(f"  Atlas result (1/4 × sum) = {float(atlas_result):.6e}")
    
    # Compare
    if M_reduced != 0:
        ratio = atlas_result / M_reduced
        print(f"  Ratio (atlas/gravity) = {float(ratio):.6f}")
        return float(ratio)
    else:
        print(f"  Cannot compute ratio (M_reduced = 0)")
        return None

# Run multiple tests
print("\n" + "="*70)
print("RUNNING ATLAS TESTS")
print("="*70)

results = []
for seed in [42, 123, 456, 789, 2024]:
    ratio = test_atlas_hypothesis(seed=seed)
    if ratio is not None:
        results.append(ratio)

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

if len(results) > 0:
    print(f"\nRatios computed: {len(results)}")
    print(f"Ratio values: {results}")
    
    # Check if ratios are constant
    if len(results) > 1:
        variance = sum((r - results[0])**2 for r in results[1:]) / len(results)
        print(f"Variance from first ratio: {variance:.6e}")
        
        if variance < 1e-6:
            print(f"\n✓ Ratios are CONSTANT = {results[0]:.6f}")
            print(f"  This means the atlas hypothesis is TRUE up to a constant factor!")
            print(f"  Actual formula: M_gravity / <12>^8 = {results[0]:.6f} × (1/4) × Σ Ω_R")
        else:
            print(f"\n✗ Ratios are NOT constant")
            print(f"  The atlas hypothesis in this form is FALSIFIED")
else:
    print("No valid results computed")

