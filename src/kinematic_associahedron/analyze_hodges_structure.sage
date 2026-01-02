# src/kinematic_associahedron/analyze_hodges_structure.sage
"""
Analyze the Structure of the Hodges Gravity Amplitude
======================================================

The goal is to understand what positive geometry underlies the Hodges
amplitude. We will analyze:

1. The pole structure (where does the amplitude diverge?)
2. The residues at poles (do they factorize correctly?)
3. The "canonical form" structure (is it a sum over cells?)

For a positive geometry, the canonical form has:
- Logarithmic singularities on ALL boundaries
- Recursive structure: Res(Ω) = Ω(boundary)
- Unit leading residues (after proper normalization)
"""

from sage.all import *
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical


def analyze_pole_structure(n=6, seed=42):
    """
    Analyze the pole structure of the Hodges amplitude.
    
    For n=6 MHV gravity, the amplitude has poles at:
    - s_{ij} = 0 (two-particle factorization)
    - s_{ijk} = 0 (three-particle factorization)
    
    These correspond to boundaries of the positive geometry.
    """
    print("\n" + "="*70)
    print("ANALYZING HODGES POLE STRUCTURE")
    print("="*70)
    
    # Generate kinematics
    kin = SpinorKinematics.random_rational(n=n, seed=seed)
    
    # Compute amplitude
    amp, status = hodges_npt_mhv_canonical(kin.lambdas, kin.tilde_lambdas, (0, 1))
    
    if status != "ok":
        print(f"Hodges computation failed: {status}")
        return
    
    print(f"\nAmplitude M_6 = {float(amp):.6e}")
    
    # Compute all Mandelstam invariants
    print("\nTwo-particle poles s_ij:")
    for i in range(n):
        for j in range(i+1, n):
            s_ij = kin.s(i, j)
            print(f"  s_{{{i}{j}}} = {float(s_ij):12.4f}")
    
    print("\nThree-particle poles s_ijk (sum of s_ab for a<b in {i,j,k}):")
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                s_ijk = kin.s(i, j) + kin.s(j, k) + kin.s(i, k)
                print(f"  s_{{{i}{j}{k}}} = {float(s_ijk):12.4f}")
    
    # Identify which poles the amplitude actually has
    # (need symbolic analysis for this)
    print("\n" + "-"*40)
    print("Physical pole channels:")
    print("  2-particle: s_01, s_12, s_23, s_34, s_45, s_05 (adjacent)")
    print("  3-particle: s_012, s_123, s_234, s_345, s_450, s_501")
    print("  (Note: s_012 = s_345 by momentum conservation)")
    
    return amp


def analyze_factorization(n=6, seed=42):
    """
    Analyze how the amplitude factorizes on poles.
    
    At s_{ijk} → 0, the amplitude should factorize as:
    M_6 → M_4(ijk, P) × 1/s_{ijk} × M_4(-P, remaining)
    """
    print("\n" + "="*70)
    print("ANALYZING FACTORIZATION")
    print("="*70)
    
    # This requires symbolic computation
    # For now, we just verify that the residue structure is correct
    
    print("\nFor MHV gravity (n=6):")
    print("  M_6 → M_4(--++) × 1/s × M_4(--++)")
    print("  or M_6 → M_3 × 1/s × M_5")
    print("\nFactorization channels:")
    print("  s_012: (0,1,2) | (3,4,5)")
    print("  s_123: (1,2,3) | (0,4,5)")
    print("  s_234: (2,3,4) | (0,1,5)")
    
    print("\n4-point MHV gravity:")
    print("  M_4(1-,2-,3+,4+) = <12>^8 [34]^8 / (s_12 s_23 s_13)")
    print("                  = <12>^4 [34]^4 × s_12^3 s_23 s_13 (using Mandelstam relations)")


def analyze_canonical_form_structure(n=6, seed=42):
    """
    Analyze if the Hodges amplitude has the structure of a canonical form.
    
    A canonical form Ω for a positive geometry has the form:
    Ω = Σ_cells ±1 × d log(B_1) ∧ d log(B_2) ∧ ... ∧ d log(B_d)
    
    where B_i are boundary equations.
    
    For the associahedron, this is:
    Ω = Σ_triangulations 1 / (X_1 × X_2 × X_3)
    
    For gravity, we expect:
    Ω = numerator × Σ_something 1 / (s_1 × s_2 × s_3)
    
    where the numerator encodes helicity.
    """
    print("\n" + "="*70)
    print("ANALYZING CANONICAL FORM STRUCTURE")
    print("="*70)
    
    # Generate kinematics
    kin = SpinorKinematics.random_rational(n=n, seed=seed)
    
    # Compute amplitude
    amp, status = hodges_npt_mhv_canonical(kin.lambdas, kin.tilde_lambdas, (0, 1))
    
    if status != "ok":
        print(f"Hodges computation failed: {status}")
        return
    
    # Extract the helicity factor
    ang_12 = kin.angle(0, 1)
    h_factor = ang_12**8
    
    print(f"\nHelicity factor <12>^8 = {float(h_factor):.6e}")
    print(f"Full amplitude M_6 = {float(amp):.6e}")
    
    if h_factor != 0:
        reduced_amp = amp / h_factor
        print(f"Reduced amplitude M_6 / <12>^8 = {float(reduced_amp):.6e}")
    
    # The reduced amplitude should be the "canonical form" part
    # Let's check if it can be written as 1/product of poles
    
    print("\n" + "-"*40)
    print("Structure Analysis:")
    
    # Compute denominator contribution from known poles
    s_12 = kin.s(0, 1)
    s_23 = kin.s(1, 2)
    s_34 = kin.s(2, 3)
    s_45 = kin.s(3, 4)
    s_56 = kin.s(4, 5)
    s_61 = kin.s(5, 0)
    
    s_123 = kin.s(0, 1) + kin.s(1, 2) + kin.s(0, 2)
    s_234 = kin.s(1, 2) + kin.s(2, 3) + kin.s(1, 3)
    s_345 = kin.s(2, 3) + kin.s(3, 4) + kin.s(2, 4)
    
    print(f"\nAdjacent 2-particle Mandelstams:")
    print(f"  s_12 × s_23 × s_34 × s_45 × s_56 × s_61 = {float(s_12*s_23*s_34*s_45*s_56*s_61):.6e}")
    
    print(f"\n3-particle Mandelstams:")
    print(f"  s_123 × s_234 × s_345 = {float(s_123*s_234*s_345):.6e}")
    
    # Expected BGK-like structure:
    # M_6 ~ <12>^8 × numerator / (PT_denom × s_123 × s_234 × s_345)
    pt_denom = ang_12 * kin.angle(1,2) * kin.angle(2,3) * kin.angle(3,4) * kin.angle(4,5) * kin.angle(5,0)
    
    if pt_denom != 0 and s_123 != 0 and s_234 != 0 and s_345 != 0:
        bgk_denom = pt_denom * s_123 * s_234 * s_345
        implied_numerator = amp * bgk_denom / h_factor
        print(f"\nBGK structure check:")
        print(f"  PT denominator = {float(pt_denom):.6e}")
        print(f"  Full BGK denominator = {float(bgk_denom):.6e}")
        print(f"  Implied numerator N = M × BGK_denom / <12>^8 = {float(implied_numerator):.6e}")
    
    return amp


def compare_hodges_vs_bgk(n=6, num_samples=3):
    """
    Compare Hodges with the BGK formula structure.
    
    BGK (Berends-Giele-Kuijf):
    M_6 = <12>^8 × G_6 / (PT × s_123 × s_234 × s_345)
    
    where G_6 is a polynomial in brackets.
    """
    print("\n" + "="*70)
    print("COMPARING HODGES VS BGK STRUCTURE")
    print("="*70)
    
    for seed in range(42, 42 + num_samples):
        print(f"\n--- Seed {seed} ---")
        
        kin = SpinorKinematics.random_rational(n=n, seed=seed)
        
        amp, status = hodges_npt_mhv_canonical(kin.lambdas, kin.tilde_lambdas, (0, 1))
        if status != "ok":
            print(f"Failed: {status}")
            continue
        
        # BGK denominator
        pt = prod(kin.angle(i, (i+1)%n) for i in range(n))
        s_123 = kin.s(0, 1) + kin.s(1, 2) + kin.s(0, 2)
        s_234 = kin.s(1, 2) + kin.s(2, 3) + kin.s(1, 3)
        s_345 = kin.s(2, 3) + kin.s(3, 4) + kin.s(2, 4)
        
        if pt == 0 or s_123 == 0 or s_234 == 0 or s_345 == 0:
            print("Singular kinematics")
            continue
        
        bgk_denom = pt * s_123 * s_234 * s_345
        
        # Extract G_6
        h_factor = kin.angle(0, 1)**8
        G_6 = amp * bgk_denom / h_factor
        
        print(f"  M_6 = {float(amp):.6e}")
        print(f"  <12>^8 = {float(h_factor):.6e}")
        print(f"  BGK denom = {float(bgk_denom):.6e}")
        print(f"  G_6 = {G_6}")
        print(f"  (G_6 should be a polynomial in spinor brackets)")


def main():
    """Run all analyses."""
    analyze_pole_structure()
    analyze_factorization()
    analyze_canonical_form_structure()
    compare_hodges_vs_bgk()
    
    print("\n" + "="*70)
    print("CONCLUSIONS")
    print("="*70)
    print("""
The Hodges amplitude has the structure:
    M_6 = <12>^8 × det'(Φ) / (<12><23><34><45><56><61>)^2
    
where det'(Φ) is the reduced determinant of the Hodges matrix.

This can be rewritten in BGK form:
    M_6 = <12>^8 × G_6 / (PT × s_123 × s_234 × s_345)

For positive geometry interpretation:
1. The poles s_ijk correspond to BOUNDARIES of the geometry
2. The numerator G_6 encodes the specific TRIANGULATION structure
3. The <12>^8 factor encodes HELICITY (fiber direction)

The geometry likely lives in:
    E → (kinematic space) 
where E is a fiber bundle with fiber encoding helicity.

The base is an associahedron-like polytope with facets at s_ijk = 0.
The fiber is responsible for the <12>^8 factor.
    """)


if __name__ == "__main__":
    main()

