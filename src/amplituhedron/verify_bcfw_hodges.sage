#!/usr/bin/env sage
"""
Verify BCFW Cell Sum = Hodges Determinant
==========================================

THE KEY TEST: Does the BCFW cell sum equal the Hodges determinant?

If YES → We have a valid amplituhedron construction for gravity!
If NO → Analyze the discrepancy and iterate on the formula

This script performs systematic testing with:
- Exact rational arithmetic
- Multiple kinematic points
- Positivity region verification
"""

from sage.all import *
import sys
import os
import time

# Add project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.amplituhedron.momentum_twistor import MomentumTwistorData, generate_positive_twistors
from src.amplituhedron.bcfw_cells import amplitude_from_cells, gravity_mhv_bcfw_sum
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor


def reference_spinors():
    """Standard reference spinors."""
    return vector(QQ, [1, 0]), vector(QQ, [0, 1])


def hodges_from_spinors(n, seed, roots=(0, 1, 2)):
    """
    Compute Hodges determinant using the trusted laplacian_bridge.
    
    This is the ORACLE - the ground truth we're trying to match.
    """
    lambdas, tildes = sample_spinors_from_twistor(n=n, seed=seed)
    x, y = reference_spinors()
    
    M, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=roots)
    return M, lambdas, tildes


def hodges_from_momentum_twistors(tw, roots=(0, 1, 2)):
    """
    Compute Hodges determinant directly from momentum twistors.
    
    Uses the matrix formula:
    M_n = det(Φ) / ∏⟨i i+1⟩
    
    where Φ is the Hodges matrix.
    """
    n = tw.n
    
    # Reference legs (typically first and last)
    X, Y = roots[0], roots[-1]
    
    # Indices to keep in matrix (delete 3 for 6-point)
    kept_indices = [i for i in range(n) if i not in roots]
    d = len(kept_indices)
    
    if d == 0:
        return None
    
    # Build Hodges matrix
    Phi = matrix(QQ, d, d)
    
    for ii, i in enumerate(kept_indices):
        for jj, j in enumerate(kept_indices):
            if ii == jj:
                # Diagonal entry: -∑_{k≠i,X,Y} [ik]⟨kX⟩⟨kY⟩ / (⟨ik⟩⟨iX⟩⟨iY⟩)
                diag_sum = QQ(0)
                for k in range(n):
                    if k == i or k in roots:
                        continue
                    
                    ik_sq = tw.square_bracket(i, k)
                    if ik_sq is None:
                        continue
                    
                    ik_ang = tw.angle(i, k)
                    iX_ang = tw.angle(i, X)
                    iY_ang = tw.angle(i, Y)
                    kX_ang = tw.angle(k, X)
                    kY_ang = tw.angle(k, Y)
                    
                    if ik_ang == 0 or iX_ang == 0 or iY_ang == 0:
                        continue
                    
                    contrib = ik_sq * kX_ang * kY_ang / (ik_ang * iX_ang * iY_ang)
                    diag_sum -= contrib
                
                Phi[ii, jj] = diag_sum
            else:
                # Off-diagonal: [ij]/⟨ij⟩
                ij_ang = tw.angle(i, j)
                if ij_ang == 0:
                    return None
                
                ij_sq = tw.square_bracket(i, j)
                if ij_sq is None:
                    return None
                
                Phi[ii, jj] = ij_sq / ij_ang
    
    # Compute determinant
    try:
        det_Phi = Phi.det()
    except:
        return None
    
    # Denominator: ∏⟨i i+1⟩
    denom = QQ(1)
    for i in range(n):
        ip1 = (i + 1) % n
        ang = tw.angle(i, ip1)
        if ang == 0:
            return None
        denom *= ang
    
    return det_Phi / denom


def test_bcfw_equals_hodges_spinor(n_trials=50):
    """
    Test BCFW = Hodges using spinor-helicity variables.
    
    Uses the trusted laplacian_bridge as oracle.
    """
    print("="*70)
    print("TEST: BCFW Cell Sum = Hodges (Spinor-Helicity)")
    print("="*70)
    
    n = 6
    roots = (0, 1, 2)
    
    matches = 0
    constant_ratios = []
    failures = []
    
    for trial in range(n_trials):
        seed = 1000 + trial * 7  # Spread out seeds
        
        try:
            # Get Hodges from oracle
            hodges, lambdas, tildes = hodges_from_spinors(n, seed, roots)
            
            if hodges is None or hodges == 0:
                continue
            
            # Create momentum twistor from same seed for comparison
            tw = MomentumTwistorData(n=n, seed=seed)
            
            if tw.is_singular():
                continue
            
            # Compute BCFW sum
            bcfw_sum, cells, contribs = amplitude_from_cells(tw)
            
            if bcfw_sum is None:
                continue
            
            # Compare
            if bcfw_sum == hodges:
                matches += 1
            else:
                if hodges != 0:
                    ratio = bcfw_sum / hodges
                    constant_ratios.append(ratio)
                    failures.append({
                        'seed': seed,
                        'bcfw': bcfw_sum,
                        'hodges': hodges,
                        'ratio': ratio
                    })
        
        except Exception as e:
            continue
    
    # Report
    print(f"\nResults: {matches} exact matches out of {n_trials} trials")
    
    if constant_ratios:
        # Check if ratios are constant
        first = constant_ratios[0]
        all_same = all(r == first for r in constant_ratios)
        
        if all_same:
            print(f"[PROMISING] All ratios are identical: {first}")
            print("→ BCFW formula correct up to constant normalization!")
        else:
            # Check numeric similarity
            numeric_ratios = [float(r) for r in constant_ratios]
            mean_r = sum(numeric_ratios) / len(numeric_ratios)
            var_r = sum((r - mean_r)**2 for r in numeric_ratios) / len(numeric_ratios)
            std_r = var_r ** 0.5
            
            print(f"Ratio statistics:")
            print(f"  Mean: {mean_r:.6e}")
            print(f"  Std:  {std_r:.6e}")
            print(f"  CV:   {std_r/abs(mean_r):.6e}")
            
            if std_r / abs(mean_r) < 0.01:
                print("[PROMISING] Ratios are approximately constant!")
    
    return matches, constant_ratios, failures


def test_bcfw_equals_hodges_twistor(n_trials=50):
    """
    Test BCFW = Hodges using momentum twistor variables only.
    
    This tests consistency within the twistor framework.
    """
    print("\n" + "="*70)
    print("TEST: BCFW Cell Sum = Hodges (Momentum Twistors)")
    print("="*70)
    
    n = 6
    roots = (0, 1, 2)
    
    matches = 0
    failures = []
    
    for trial in range(n_trials):
        seed = 2000 + trial * 11
        
        try:
            tw = MomentumTwistorData(n=n, seed=seed)
            
            if tw.is_singular():
                continue
            
            # Hodges from twistors
            hodges = hodges_from_momentum_twistors(tw, roots)
            
            if hodges is None or hodges == 0:
                continue
            
            # BCFW sum
            bcfw_sum, cells, contribs = amplitude_from_cells(tw)
            
            if bcfw_sum is None:
                continue
            
            # Compare
            if bcfw_sum == hodges:
                matches += 1
                print(f"  [MATCH] Seed {seed}: bcfw = hodges = {hodges}")
            else:
                ratio = bcfw_sum / hodges if hodges != 0 else None
                failures.append({
                    'seed': seed,
                    'bcfw': bcfw_sum,
                    'hodges': hodges,
                    'ratio': ratio
                })
        
        except Exception as e:
            continue
    
    print(f"\nResults: {matches} exact matches out of {n_trials} trials")
    
    if failures and len(failures) <= 5:
        print("\nSample failures:")
        for f in failures[:5]:
            print(f"  Seed {f['seed']}: ratio = {f['ratio']}")
    
    return matches, failures


def analyze_discrepancy():
    """
    Detailed analysis of BCFW vs Hodges discrepancy.
    """
    print("\n" + "="*70)
    print("DISCREPANCY ANALYSIS")
    print("="*70)
    
    n = 6
    roots = (0, 1, 2)
    
    # Find a good test point
    for seed in range(100):
        tw = MomentumTwistorData(n=n, seed=seed)
        
        if tw.is_singular():
            continue
        
        hodges = hodges_from_momentum_twistors(tw, roots)
        
        if hodges is None or hodges == 0:
            continue
        
        bcfw_sum, cells, contribs = amplitude_from_cells(tw)
        
        if bcfw_sum is None:
            continue
        
        print(f"\nTest point: seed = {seed}")
        print(f"  Hodges:    {hodges}")
        print(f"  BCFW sum:  {bcfw_sum}")
        
        if hodges != 0:
            ratio = bcfw_sum / hodges
            print(f"  Ratio:     {ratio}")
            
            # Check if ratio is a simple expression
            print(f"\nAnalyzing ratio structure...")
            
            # Try to express as product of brackets
            # ratio * hodges = bcfw_sum
            # If ratio = ⟨ij⟩^k / ⟨kl⟩^m for some i,j,k,l,m
            # we can identify the missing factor
            
            print(f"  Ratio numerator:   {ratio.numerator() if hasattr(ratio, 'numerator') else 'N/A'}")
            print(f"  Ratio denominator: {ratio.denominator() if hasattr(ratio, 'denominator') else 'N/A'}")
        
        print(f"\nCell contributions:")
        for i, c in enumerate(contribs[:5]):
            print(f"  Cell {i}: {c['contribution']}")
        
        break
    
    return


def main():
    """Run all verification tests."""
    print("="*70)
    print("BCFW = HODGES VERIFICATION SUITE")
    print("="*70)
    print()
    
    t_start = time.time()
    
    # Test 1: Spinor-helicity based comparison
    matches1, ratios1, failures1 = test_bcfw_equals_hodges_spinor(n_trials=30)
    
    # Test 2: Pure momentum twistor comparison
    matches2, failures2 = test_bcfw_equals_hodges_twistor(n_trials=30)
    
    # Test 3: Detailed discrepancy analysis
    analyze_discrepancy()
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    total_matches = matches1 + matches2
    
    if total_matches > 0:
        print(f"[SUCCESS] Found {total_matches} exact matches!")
        print("→ BCFW cell structure is promising")
    else:
        print("[INVESTIGATING] No exact matches found")
        print("→ Formula needs refinement")
        
        if ratios1:
            # Check ratio pattern
            numeric = [float(r) for r in ratios1[:10]]
            print(f"\nSample ratios: {numeric}")
    
    print(f"\nTotal time: {time.time() - t_start:.1f}s")


if __name__ == '__main__':
    main()


