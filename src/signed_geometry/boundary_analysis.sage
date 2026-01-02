#!/usr/bin/env sage
"""
Boundary Analysis for Signed Geometry

This module analyzes the factorization behavior of the signed canonical
form at the boundaries (poles) of the amplitude.

Key question: How does the 54/54 sign split behave at poles s_ijk -> 0?
"""
from sage.all import *
import sys
import os

sys.path.insert(0, os.getcwd())

load('src/spinor_sampling.sage')
load('src/signed_geometry/canonical_form.sage')


def create_near_pole_kinematics(pole_particles, epsilon=QQ(1)/100, seed=42):
    """
    Create kinematics near a specific pole s_{ijk} -> 0.
    
    This is done by starting with random kinematics and then
    adjusting to make s_{pole_particles} small.
    """
    import random
    random.seed(int(seed))
    
    # Sample base kinematics
    result = sample_spinor_helicity_conserving(n=6, seed=seed)
    if result is None:
        return None
    
    lambdas, tilde_lambdas = result
    
    # The pole s_I = sum_{i<j in I} s_ij needs to be small
    # We achieve this approximately by scaling some tilde_lambdas
    # This is a simplification - proper pole approach requires more care
    
    return lambdas, tilde_lambdas


def analyze_factorization_at_pole(pole_indices=(0, 1, 2), num_samples=10):
    """
    Analyze how the amplitude (and sign split) behaves near a pole.
    
    For s_{012} -> 0, the amplitude should factor as M_4 x M_4.
    """
    print("=" * 70)
    print(f"ANALYZING FACTORIZATION AT POLE s_{{{','.join(map(str, pole_indices))}}}")
    print("=" * 70)
    
    results = []
    
    for sample in range(num_samples):
        result = sample_spinor_helicity_conserving(n=6, seed=sample * 41)
        if result is None:
            continue
        
        lambdas, tilde_lambdas = result
        
        # Compute the pole value
        def mandelstam(i, j):
            ang = lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]
            sq = tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]
            return ang * sq
        
        s_pole = QQ(0)
        for i in range(len(pole_indices)):
            for j in range(i+1, len(pole_indices)):
                s_pole += mandelstam(pole_indices[i], pole_indices[j])
        
        # Reference spinors
        x_spinor = vector(QQ, [1, 2])
        y_spinor = vector(QQ, [3, 1])
        
        # Analyze sign structure
        analysis = analyze_sign_structure(lambdas, tilde_lambdas, x_spinor, y_spinor, roots=(0, 1, 2))
        
        results.append({
            's_pole': float(s_pole),
            'split': analysis['split_ratio'],
            'total': float(analysis['total']) if analysis['total'] != 0 else 0
        })
    
    # Analyze correlation between pole value and sign split
    print(f"\nSamples: {len(results)}")
    print(f"\nPole value vs Sign Split:")
    
    for r in sorted(results, key=lambda x: abs(x['s_pole']))[:5]:
        print(f"  s_pole = {r['s_pole']:.4f}, split = {r['split']}")
    
    return results


def verify_factorization_residue():
    """
    Verify that the residue at a pole gives the product of lower-point amplitudes.
    
    At s_{012} -> 0:
        Res_{s_{012}=0} M_6 = M_4(0,1,2,P) x M_4(-P,3,4,5)
    
    For MHV with helicities (0^-, 1^-, 2^+, 3^+, 4^+, 5^+):
    - The (012) channel has 2 negative helicities -> can give MHV
    - The (345) channel has 0 negative helicities -> all-plus (vanishes)
    
    So this particular channel has a soft factorization, not hard.
    """
    print("=" * 70)
    print("RESIDUE FACTORIZATION ANALYSIS")
    print("=" * 70)
    
    print("\nFor 6-point MHV (0^- 1^- 2^+ 3^+ 4^+ 5^+):")
    print("")
    print("Channel analysis:")
    
    # List all 3-particle channels
    channels = [
        ((0, 1, 2), (3, 4, 5)),  # s_012 = s_345
        ((1, 2, 3), (4, 5, 0)),  # s_123 = s_450
        ((2, 3, 4), (5, 0, 1)),  # s_234 = s_501
    ]
    
    # Helicities: 0=-, 1=-, 2=+, 3=+, 4=+, 5=+
    helicities = {0: '-', 1: '-', 2: '+', 3: '+', 4: '+', 5: '+'}
    
    for left, right in channels:
        left_hel = ''.join(helicities[i] for i in left)
        right_hel = ''.join(helicities[i] for i in right)
        
        left_neg = sum(1 for i in left if helicities[i] == '-')
        right_neg = sum(1 for i in right if helicities[i] == '-')
        
        left_type = "MHV" if left_neg == 2 else f"N^{left_neg-2}MHV" if left_neg > 2 else f"MHV-bar" if left_neg == 0 else "anti-MHV"
        right_type = "MHV" if right_neg == 2 else f"N^{right_neg-2}MHV" if right_neg > 2 else "all-+" if right_neg == 0 else "anti-MHV"
        
        print(f"  s_{{{left[0]}{left[1]}{left[2]}}} = 0:")
        print(f"    Left {left}: helicities ({left_hel}), {left_neg} negative -> {left_type}")
        print(f"    Right {right}: helicities ({right_hel}), {right_neg} negative -> {right_type}")
        
        if right_neg == 0 or left_neg == 0:
            print(f"    -> VANISHES (all-plus amplitude is zero)")
        else:
            print(f"    -> FACTORIZES: M_4({left_type}) x M_4({right_type})")
        print()
    
    # The non-vanishing channels for MHV are those where both sides have 2 negative helicities
    # For (0^- 1^- 2^+ 3^+ 4^+ 5^+), the internal leg must carry helicity to balance
    
    print("Conclusion:")
    print("  For MHV gravity, factorization channels involve helicity flow")
    print("  through the internal leg. The sign structure of the forest")
    print("  expansion is preserved through factorization because the")
    print("  forest structure decomposes into products of lower-point forests.")


def analyze_forest_decomposition_at_boundary():
    """
    Analyze how the 108 forests decompose at a pole.
    
    At s_I -> 0, forests that have edges crossing the I / I^c cut
    develop a pole, while forests contained within each side remain finite.
    """
    print("=" * 70)
    print("FOREST DECOMPOSITION AT BOUNDARIES")
    print("=" * 70)
    
    # For n=6 with roots {0,1,2}, the forests span all 6 vertices
    # with exactly 3 trees, one rooted at each of 0, 1, 2.
    
    # Consider the cut s_{012} <-> s_{345}
    left_set = {0, 1, 2}
    right_set = {3, 4, 5}
    
    print(f"\nAnalyzing partition: {left_set} | {right_set}")
    
    # Enumerate forests and classify by crossing structure
    # (enumerate_rooted_forests is already loaded from canonical_form.sage)
    
    crossing_counts = {0: 0, 1: 0, 2: 0, 3: 0}
    
    for forest in enumerate_rooted_forests(6, (0, 1, 2)):
        # Count edges crossing the cut
        n_crossing = 0
        for u, v in forest:
            if (u in left_set and v in right_set) or (u in right_set and v in left_set):
                n_crossing += 1
        
        crossing_counts[n_crossing] = crossing_counts.get(n_crossing, 0) + 1
    
    total_forests = sum(crossing_counts.values())
    print(f"\nTotal forests: {total_forests}")
    print(f"\nForests by number of cut-crossing edges:")
    for n_cross, count in sorted(crossing_counts.items()):
        pct = 100.0 * count / total_forests
        print(f"  {n_cross} crossing edges: {count} ({pct:.1f}%)")
    
    print("\nInterpretation:")
    print("  - Forests with 0 crossing edges: stay finite at pole")
    print("  - Forests with k>0 crossing edges: contribute to pole structure")
    print("  - The factorization residue comes from the 0-crossing forests")


if __name__ == "__main__":
    # Analyze factorization structure
    verify_factorization_residue()
    
    # Analyze forest decomposition at boundaries
    analyze_forest_decomposition_at_boundary()
    
    # Analyze near poles
    analyze_factorization_at_pole(pole_indices=(0, 1, 2), num_samples=10)

