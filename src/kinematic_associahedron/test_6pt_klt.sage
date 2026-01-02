# src/kinematic_associahedron/test_6pt_klt.sage
"""
Test KLT at 6-point with the Sign Fix
======================================

From 4-point testing, we found: KLT = -Hodges (ratio = -1 exactly)

This is the expected (-1)^{n+1} factor for n=4: (-1)^5 = -1

For n=6: (-1)^7 = -1

So we expect the same relationship: KLT = -Hodges

Let's verify this!
"""

from sage.all import *
from itertools import permutations
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical


def angle_bracket(lambdas, i, j):
    """Compute <ij>."""
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def square_bracket(tilde_lambdas, i, j):
    """Compute [ij]."""
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def mandelstam(lambdas, tilde_lambdas, i, j):
    """Compute s_ij = <ij>[ji]."""
    return angle_bracket(lambdas, i, j) * square_bracket(tilde_lambdas, j, i)


def ym_amplitude_6pt(lambdas, order, neg_helicity=(0, 1)):
    """
    Compute 6-point MHV Yang-Mills amplitude.
    
    A_6 = <ab>^4 / (<order[0] order[1]> ... <order[5] order[0]>)
    """
    a, b = neg_helicity
    numerator = angle_bracket(lambdas, a, b) ** 4
    
    if numerator == 0:
        return None
    
    denom = QQ(1)
    n = len(order)
    for i in range(n):
        j = (i + 1) % n
        bracket = angle_bracket(lambdas, order[i], order[j])
        if bracket == 0:
            return None
        denom *= bracket
    
    return numerator / denom


def klt_kernel_6pt(alpha, beta, lambdas, tilde_lambdas):
    """
    Compute KLT momentum kernel S[α|β] for n=6.
    
    S[α|β] = Π_{i∈α'} (s_{1,α_i} + Σ_{j<i, θ_β(α_j,α_i)=1} s_{α_j,α_i})
    
    where α' = {α_0, α_1, α_2} (the permuted subset)
    """
    alpha = list(alpha)
    beta = list(beta)
    
    # Position map for beta
    pos_beta = {beta[i]: i for i in range(len(beta))}
    
    # θ_β(a,b) = 1 if a comes after b in β
    def theta(a, b):
        return 1 if pos_beta.get(a, -1) > pos_beta.get(b, -1) else 0
    
    kernel = QQ(1)
    
    for i in range(len(alpha)):
        # s_{1, α_i} where 1 is index 0
        s_1_ai = mandelstam(lambdas, tilde_lambdas, 0, alpha[i])
        term = s_1_ai
        
        for j in range(i):
            if theta(alpha[j], alpha[i]):
                s_aj_ai = mandelstam(lambdas, tilde_lambdas, alpha[j], alpha[i])
                term += s_aj_ai
        
        kernel *= term
    
    return kernel


def gravity_6pt_klt(lambdas, tilde_lambdas, neg_helicity=(0, 1)):
    """
    Compute 6-point MHV gravity amplitude via KLT double copy.
    
    M_6 = Σ_{α,β} A_YM(1,α,5,6) × S[α|β] × A_YM(1,β,6,5)
    """
    n = 6
    permuted_set = [1, 2, 3]  # Internal particles
    all_perms = list(permutations(permuted_set))
    
    total = QQ(0)
    
    for alpha in all_perms:
        # A_YM(1,α,5,6) → order [0] + list(alpha) + [4, 5]
        order_alpha = [0] + list(alpha) + [4, 5]
        A_alpha = ym_amplitude_6pt(lambdas, order_alpha, neg_helicity)
        if A_alpha is None:
            continue
        
        for beta in all_perms:
            # A_YM(1,β,6,5) → order [0] + list(beta) + [5, 4]
            order_beta = [0] + list(beta) + [5, 4]
            A_beta = ym_amplitude_6pt(lambdas, order_beta, neg_helicity)
            if A_beta is None:
                continue
            
            # KLT kernel
            S = klt_kernel_6pt(alpha, beta, lambdas, tilde_lambdas)
            
            # Add contribution
            total += A_alpha * S * A_beta
    
    return total


def test_6pt_klt_vs_hodges():
    """Test 6-point KLT against Hodges."""
    print("\n" + "="*70)
    print("TESTING 6-POINT KLT vs HODGES")
    print("="*70)
    print("\nExpected: KLT = -Hodges (from 4-point analysis)")
    print("Convention: (-1)^{n+1} = (-1)^7 = -1\n")
    
    ratios = []
    
    for seed in range(42, 52):
        kin = SpinorKinematics.random_rational(n=6, seed=seed)
        lambdas = kin.lambdas
        tilde_lambdas = kin.tilde_lambdas
        
        # KLT
        M_klt = gravity_6pt_klt(lambdas, tilde_lambdas, (0, 1))
        
        # Hodges
        M_hodges, status = hodges_npt_mhv_canonical(lambdas, tilde_lambdas, (0, 1))
        
        if M_hodges is None:
            print(f"Seed {seed}: Hodges failed [{status}]")
            continue
        
        print(f"Seed {seed}:")
        print(f"  M_6 (KLT):    {float(M_klt):>15.6e}")
        print(f"  M_6 (Hodges): {float(M_hodges):>15.6e}")
        
        if M_hodges != 0:
            ratio = M_klt / M_hodges
            print(f"  Ratio: {float(ratio):>15.10f}")
            ratios.append(ratio)
            
            if abs(float(ratio) + 1.0) < 1e-10:
                print(f"  ✓ KLT = -Hodges (as expected)")
            elif abs(float(ratio) - 1.0) < 1e-10:
                print(f"  ✓ KLT = +Hodges")
            else:
                print(f"  ✗ Unexpected ratio")
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    if ratios:
        # Check if all ratios are the same
        first_ratio = ratios[0]
        all_same = all(r == first_ratio for r in ratios)
        
        if all_same:
            print(f"✓ All ratios are identical: {float(first_ratio):.10f}")
            if float(first_ratio) == -1.0:
                print(f"✓ KLT = -Hodges confirmed!")
                print(f"\nThe sign comes from (-1)^{{n+1}} = (-1)^7 = -1")
            elif float(first_ratio) == 1.0:
                print(f"✓ KLT = Hodges confirmed!")
        else:
            print(f"✗ Ratios are NOT identical!")
            print(f"  Ratios: {[float(r) for r in ratios]}")
    else:
        print("No valid samples to compare")


if __name__ == "__main__":
    test_6pt_klt_vs_hodges()

