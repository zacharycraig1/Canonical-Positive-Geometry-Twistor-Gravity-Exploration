# src/kinematic_associahedron/debug_klt_6pt.sage
"""
Debug KLT at 6-point
====================

The 4-point case shows KLT = -Hodges exactly.
The 6-point case shows varying ratios.

This suggests either:
1. The kernel formula is wrong for n=6
2. The amplitude orderings are inconsistent
3. There's an additional factor we're missing

Let's debug systematically.
"""

from sage.all import *
from itertools import permutations
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical


def angle_bracket(lambdas, i, j):
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def square_bracket(tilde_lambdas, i, j):
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def mandelstam(lambdas, tilde_lambdas, i, j):
    """s_ij = <ij>[ji]"""
    return angle_bracket(lambdas, i, j) * square_bracket(tilde_lambdas, j, i)


def ym_amplitude(lambdas, order, neg_helicity=(0, 1)):
    """YM MHV amplitude: <ab>^4 / prod <order[i] order[i+1]>"""
    a, b = neg_helicity
    numer = angle_bracket(lambdas, a, b) ** 4
    if numer == 0:
        return None
    
    denom = QQ(1)
    n = len(order)
    for i in range(n):
        j = (i + 1) % n
        bracket = angle_bracket(lambdas, order[i], order[j])
        if bracket == 0:
            return None
        denom *= bracket
    
    return numer / denom


def klt_kernel_standard(alpha, beta, lambdas, tilde_lambdas):
    """
    Standard KLT kernel for n=6.
    
    S[α|β] = Π_{i∈α} (s_{1,α_i} + Σ_{j<i, θ_β(α_j,α_i)} s_{α_j,α_i})
    
    where θ_β(a,b) = 1 if a comes after b in β.
    """
    alpha = list(alpha)
    beta = list(beta)
    
    pos_beta = {beta[i]: i for i in range(len(beta))}
    
    def theta(a, b):
        return 1 if pos_beta.get(a, -1) > pos_beta.get(b, -1) else 0
    
    kernel = QQ(1)
    
    for i in range(len(alpha)):
        s_1_ai = mandelstam(lambdas, tilde_lambdas, 0, alpha[i])
        term = s_1_ai
        
        for j in range(i):
            if theta(alpha[j], alpha[i]):
                s_aj_ai = mandelstam(lambdas, tilde_lambdas, alpha[j], alpha[i])
                term += s_aj_ai
        
        kernel *= term
    
    return kernel


def klt_kernel_alternative(alpha, beta, lambdas, tilde_lambdas):
    """
    Alternative KLT kernel: S[β|α] instead of S[α|β].
    
    Some conventions have the indices swapped.
    """
    return klt_kernel_standard(beta, alpha, lambdas, tilde_lambdas)


def test_klt_variations(seed=42):
    """Test different KLT conventions."""
    print("\n" + "="*70)
    print(f"TESTING KLT VARIATIONS (seed={seed})")
    print("="*70)
    
    kin = SpinorKinematics.random_rational(n=6, seed=seed)
    lambdas = kin.lambdas
    tilde_lambdas = kin.tilde_lambdas
    
    # Hodges amplitude
    M_hodges, status = hodges_npt_mhv_canonical(lambdas, tilde_lambdas, (0, 1))
    if M_hodges is None:
        print(f"Hodges failed: {status}")
        return
    
    print(f"\nHodges amplitude: {float(M_hodges):.6e}")
    
    permuted_set = [1, 2, 3]
    all_perms = list(permutations(permuted_set))
    
    # Test different conventions
    conventions = []
    
    # Convention 1: Standard (what we have)
    # A(1,α,5,6) × S[α|β] × A(1,β,6,5)
    total_1 = QQ(0)
    for alpha in all_perms:
        order_alpha = [0] + list(alpha) + [4, 5]
        A_alpha = ym_amplitude(lambdas, order_alpha)
        if A_alpha is None:
            continue
        
        for beta in all_perms:
            order_beta = [0] + list(beta) + [5, 4]
            A_beta = ym_amplitude(lambdas, order_beta)
            if A_beta is None:
                continue
            
            S = klt_kernel_standard(alpha, beta, lambdas, tilde_lambdas)
            total_1 += A_alpha * S * A_beta
    
    conventions.append(("A(1,α,5,6) × S[α|β] × A(1,β,6,5)", total_1))
    
    # Convention 2: S[β|α] instead
    total_2 = QQ(0)
    for alpha in all_perms:
        order_alpha = [0] + list(alpha) + [4, 5]
        A_alpha = ym_amplitude(lambdas, order_alpha)
        if A_alpha is None:
            continue
        
        for beta in all_perms:
            order_beta = [0] + list(beta) + [5, 4]
            A_beta = ym_amplitude(lambdas, order_beta)
            if A_beta is None:
                continue
            
            S = klt_kernel_alternative(alpha, beta, lambdas, tilde_lambdas)
            total_2 += A_alpha * S * A_beta
    
    conventions.append(("A(1,α,5,6) × S[β|α] × A(1,β,6,5)", total_2))
    
    # Convention 3: Both amplitudes same ordering
    total_3 = QQ(0)
    for alpha in all_perms:
        order_alpha = [0] + list(alpha) + [4, 5]
        A_alpha = ym_amplitude(lambdas, order_alpha)
        if A_alpha is None:
            continue
        
        for beta in all_perms:
            order_beta = [0] + list(beta) + [4, 5]  # SAME as alpha
            A_beta = ym_amplitude(lambdas, order_beta)
            if A_beta is None:
                continue
            
            S = klt_kernel_standard(alpha, beta, lambdas, tilde_lambdas)
            total_3 += A_alpha * S * A_beta
    
    conventions.append(("A(1,α,5,6) × S[α|β] × A(1,β,5,6)", total_3))
    
    # Convention 4: Different fixed legs for kernel
    # Use n=6, so fixed legs 1, 6 (indices 0, 5)
    # Permuted: {2,3,4,5} but for n=6 with 3 permuted: still {2,3,4}
    # Actually, the kernel pivot is particle 1 (index 0)
    # Let's try a different approach: use s_{n,α_i} instead of s_{1,α_i}
    
    def klt_kernel_pivot_n(alpha, beta, lambdas, tilde_lambdas):
        """Kernel with pivot at particle n instead of 1."""
        alpha = list(alpha)
        beta = list(beta)
        
        pos_beta = {beta[i]: i for i in range(len(beta))}
        
        def theta(a, b):
            return 1 if pos_beta.get(a, -1) > pos_beta.get(b, -1) else 0
        
        kernel = QQ(1)
        
        for i in range(len(alpha)):
            s_n_ai = mandelstam(lambdas, tilde_lambdas, 5, alpha[i])  # s_{6,α_i}
            term = s_n_ai
            
            for j in range(i):
                if theta(alpha[j], alpha[i]):
                    s_aj_ai = mandelstam(lambdas, tilde_lambdas, alpha[j], alpha[i])
                    term += s_aj_ai
            
            kernel *= term
        
        return kernel
    
    total_4 = QQ(0)
    for alpha in all_perms:
        order_alpha = [0] + list(alpha) + [4, 5]
        A_alpha = ym_amplitude(lambdas, order_alpha)
        if A_alpha is None:
            continue
        
        for beta in all_perms:
            order_beta = [0] + list(beta) + [5, 4]
            A_beta = ym_amplitude(lambdas, order_beta)
            if A_beta is None:
                continue
            
            S = klt_kernel_pivot_n(alpha, beta, lambdas, tilde_lambdas)
            total_4 += A_alpha * S * A_beta
    
    conventions.append(("Kernel pivot at n=6", total_4))
    
    # Print results
    print(f"\nConvention tests:")
    for name, total in conventions:
        ratio = total / M_hodges if M_hodges != 0 else "N/A"
        print(f"\n  {name}")
        print(f"    KLT = {float(total):.6e}")
        print(f"    Ratio KLT/Hodges = {float(ratio):.6f}" if isinstance(ratio, Rational) else f"    Ratio = {ratio}")


def search_for_correct_convention():
    """Systematically search for the correct KLT convention."""
    print("\n" + "="*70)
    print("SEARCHING FOR CORRECT KLT CONVENTION")
    print("="*70)
    
    # Test across multiple seeds
    for seed in range(42, 47):
        test_klt_variations(seed)


def check_mandelstam_convention(seed=42):
    """Check if the Mandelstam convention is consistent."""
    print("\n" + "="*70)
    print("CHECKING MANDELSTAM CONVENTION")
    print("="*70)
    
    kin = SpinorKinematics.random_rational(n=6, seed=seed)
    lambdas = kin.lambdas
    tilde_lambdas = kin.tilde_lambdas
    
    print("\nMandelstam invariants s_ij = <ij>[ji]:")
    for i in range(6):
        for j in range(i+1, 6):
            ang = angle_bracket(lambdas, i, j)
            sq_ij = square_bracket(tilde_lambdas, i, j)
            sq_ji = square_bracket(tilde_lambdas, j, i)
            s_ij = ang * sq_ji
            s_ij_alt = ang * sq_ij
            
            print(f"  s_{i}{j}: <{i}{j}>={float(ang):.4f}, [{i}{j}]={float(sq_ij):.4f}, [{j}{i}]={float(sq_ji):.4f}")
            print(f"       s_ij = <ij>[ji] = {float(s_ij):.4f}")
            print(f"       s_ij_alt = <ij>[ij] = {float(s_ij_alt):.4f}")
    
    # Verify momentum conservation: Σ s_ij = 0 for each i
    print("\nMomentum conservation check (Σ_j s_ij should be 0):")
    for i in range(6):
        total = QQ(0)
        for j in range(6):
            if i != j:
                total += mandelstam(lambdas, tilde_lambdas, i, j)
        print(f"  Σ_j s_{i}j = {float(total):.10f}")
    
    # Check with kin.s method
    print("\nComparing with SpinorKinematics.s():")
    for i in range(6):
        for j in range(i+1, 6):
            my_s = mandelstam(lambdas, tilde_lambdas, i, j)
            kin_s = kin.s(i, j)
            print(f"  s_{i}{j}: my={float(my_s):.6f}, kin.s={float(kin_s):.6f}, match={my_s == kin_s}")


if __name__ == "__main__":
    check_mandelstam_convention()
    # search_for_correct_convention()

