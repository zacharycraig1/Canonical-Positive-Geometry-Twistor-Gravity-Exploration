# src/kinematic_associahedron/gravity_from_klt.sage
"""
Gravity Positive Geometry via KLT Double Copy
==============================================

This module computes the 6-point MHV gravity amplitude using the
KLT (Kawai-Lewellen-Tye) relations and compares with the Hodges formula.

The key equation is:

    M_n^gravity = Σ_{α,β} A_n^YM(α) × S[α|β] × A_n^YM(β)

For MHV:
    A_n^YM(α) = <12>^4 / (<α_1 α_2> <α_2 α_3> ... <α_n α_1>)

The positive geometry interpretation:
- Each A_YM corresponds to the canonical form of an associahedron
- The KLT kernel S[α|β] provides the "gluing" between two copies
- Gravity = "square" of Yang-Mills (geometrically)

Reference:
- Kawai, Lewellen, Tye (1986) - Original KLT relations
- Bern, Carrasco, Johansson - BCJ duality
- Arkani-Hamed et al. (1711.09102) - Positive geometry
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
    """Compute angle bracket <ij>."""
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def square_bracket(tilde_lambdas, i, j):
    """Compute square bracket [ij]."""
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def mandelstam_from_spinors(lambdas, tilde_lambdas, i, j):
    """Compute Mandelstam invariant s_ij = <ij>[ij]."""
    return angle_bracket(lambdas, i, j) * square_bracket(tilde_lambdas, i, j)


def parke_taylor_amplitude(lambdas, order, neg_helicity=(0, 1)):
    """
    Compute Yang-Mills MHV amplitude (Parke-Taylor formula).
    
    A_n^MHV = <ab>^4 / (<order[0] order[1]> <order[1] order[2]> ... <order[n-1] order[0]>)
    
    Args:
        lambdas: List of lambda spinors
        order: Color ordering as list of indices
        neg_helicity: Tuple of negative helicity particle indices
    
    Returns:
        Rational amplitude value, or None if singular
    """
    n = len(order)
    
    # Numerator: <ab>^4 where a, b are negative helicity particles
    a, b = neg_helicity
    numerator = angle_bracket(lambdas, a, b) ** 4
    if numerator == 0:
        return None
    
    # Denominator: cyclic product of angle brackets
    denominator = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        bracket = angle_bracket(lambdas, order[i], order[j])
        if bracket == 0:
            return None
        denominator *= bracket
    
    return numerator / denominator


def klt_kernel(alpha, beta, lambdas, tilde_lambdas):
    """
    Compute KLT momentum kernel S[α|β] for n=6.
    
    S[α|β] = Π_{i∈α'} (s_{1,α_i} + Σ_{j<i, θ_β(α_j,α_i)=1} s_{α_j,α_i})
    
    Args:
        alpha: Permutation of {1,2,3} (internal particles)
        beta: Permutation of {1,2,3}
        lambdas, tilde_lambdas: Spinor data
    
    Returns:
        Rational kernel value
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
        s_1_ai = mandelstam_from_spinors(lambdas, tilde_lambdas, 0, alpha[i])
        term = s_1_ai
        
        for j in range(i):
            if theta(alpha[j], alpha[i]):
                s_aj_ai = mandelstam_from_spinors(lambdas, tilde_lambdas, alpha[j], alpha[i])
                term += s_aj_ai
        
        kernel *= term
    
    return kernel


def gravity_via_klt(lambdas, tilde_lambdas, neg_helicity=(0, 1)):
    """
    Compute 6-point MHV gravity amplitude via KLT double copy.
    
    M_6 = Σ_{α,β} A_YM(1,α,5,6) × S[α|β] × A_YM(1,β,6,5)
    
    Args:
        lambdas: List of lambda spinors
        tilde_lambdas: List of tilde-lambda spinors
        neg_helicity: Tuple of negative helicity particle indices
    
    Returns:
        Tuple (amplitude, status)
    """
    n = 6
    permuted_set = [1, 2, 3]  # Internal particles (0-indexed: 2,3,4 → 1,2,3)
    all_perms = list(permutations(permuted_set))
    
    total = QQ(0)
    
    for alpha in all_perms:
        # A_YM(1,α,5,6) → order [0] + list(alpha) + [4, 5]
        order_alpha = [0] + list(alpha) + [4, 5]
        A_alpha = parke_taylor_amplitude(lambdas, order_alpha, neg_helicity)
        if A_alpha is None:
            continue
        
        for beta in all_perms:
            # A_YM(1,β,6,5) → order [0] + list(beta) + [5, 4]
            order_beta = [0] + list(beta) + [5, 4]
            A_beta = parke_taylor_amplitude(lambdas, order_beta, neg_helicity)
            if A_beta is None:
                continue
            
            # KLT kernel
            S = klt_kernel(alpha, beta, lambdas, tilde_lambdas)
            
            # Add contribution
            total += A_alpha * S * A_beta
    
    return (total, "ok")


def verify_klt_equals_hodges(num_samples=5, seed_start=42):
    """
    Verify that KLT double copy matches Hodges formula.
    """
    print("\n" + "="*70)
    print("VERIFYING KLT = HODGES FOR 6-POINT MHV GRAVITY")
    print("="*70)
    
    matches = 0
    tested = 0
    
    for k in range(num_samples):
        seed = seed_start + k
        print(f"\n--- Sample {k+1} (seed={seed}) ---")
        
        # Generate kinematics
        kin = SpinorKinematics.random_rational(n=6, seed=seed)
        lambdas = kin.lambdas
        tilde_lambdas = kin.tilde_lambdas
        
        # Compute via KLT
        klt_amp, klt_status = gravity_via_klt(lambdas, tilde_lambdas, neg_helicity=(0, 1))
        
        if klt_status != "ok":
            print(f"  KLT failed: {klt_status}")
            continue
        
        # Compute via Hodges
        hodges_amp, hodges_status = hodges_npt_mhv_canonical(lambdas, tilde_lambdas, (0, 1))
        
        if hodges_status != "ok":
            print(f"  Hodges failed: {hodges_status}")
            continue
        
        tested += 1
        
        print(f"  KLT amplitude:    {float(klt_amp):.10e}")
        print(f"  Hodges amplitude: {float(hodges_amp):.10e}")
        
        if hodges_amp == 0:
            print(f"  Hodges is zero, cannot compare ratio")
            continue
        
        ratio = klt_amp / hodges_amp
        print(f"  Ratio KLT/Hodges: {float(ratio):.10f}")
        
        # Check if they match (should be exactly 1 or a simple constant)
        if ratio == 1:
            print(f"  ✓ EXACT MATCH!")
            matches += 1
        elif abs(float(ratio) - 1.0) < 1e-10:
            print(f"  ✓ Match (numerical)")
            matches += 1
        else:
            print(f"  ✗ MISMATCH (ratio = {ratio})")
            
            # Check if it's a simple rational factor
            try:
                if ratio.denominator() == 1 or ratio.numerator() == 1:
                    print(f"    (ratio is simple: {ratio.numerator()}/{ratio.denominator()})")
            except:
                pass
    
    print(f"\n{'='*70}")
    print(f"SUMMARY: {matches}/{tested} samples matched")
    print(f"{'='*70}")
    
    return matches == tested


def analyze_klt_geometry(seed=42):
    """
    Analyze the geometric structure of the KLT double copy.
    """
    print("\n" + "="*70)
    print("KLT GEOMETRY ANALYSIS")
    print("="*70)
    
    # Generate kinematics
    kin = SpinorKinematics.random_rational(n=6, seed=seed)
    lambdas = kin.lambdas
    tilde_lambdas = kin.tilde_lambdas
    
    permuted_set = [1, 2, 3]
    all_perms = list(permutations(permuted_set))
    
    print(f"\nNumber of permutations: {len(all_perms)}")
    print(f"Permuted set: {permuted_set}")
    
    # Compute YM amplitude vector
    print(f"\nYang-Mills amplitudes A_YM(1,α,5,6):")
    A_left = []
    for alpha in all_perms:
        order = [0] + list(alpha) + [4, 5]
        A = parke_taylor_amplitude(lambdas, order, (0, 1))
        A_left.append(A if A else QQ(0))
        print(f"  α={alpha}: A = {float(A) if A else 'None':.6e}")
    
    print(f"\nYang-Mills amplitudes A_YM(1,β,6,5):")
    A_right = []
    for beta in all_perms:
        order = [0] + list(beta) + [5, 4]
        A = parke_taylor_amplitude(lambdas, order, (0, 1))
        A_right.append(A if A else QQ(0))
        print(f"  β={beta}: A = {float(A) if A else 'None':.6e}")
    
    # Compute KLT kernel matrix
    print(f"\nKLT Kernel Matrix S[α|β]:")
    S = matrix(QQ, len(all_perms), len(all_perms))
    for i, alpha in enumerate(all_perms):
        for j, beta in enumerate(all_perms):
            S[i, j] = klt_kernel(alpha, beta, lambdas, tilde_lambdas)
    
    for i, alpha in enumerate(all_perms):
        row = [f"{float(S[i,j]):10.2f}" for j in range(len(all_perms))]
        print(f"  {alpha}: [{', '.join(row)}]")
    
    # Compute gravity as matrix product
    A_left_vec = vector(QQ, A_left)
    A_right_vec = vector(QQ, A_right)
    
    gravity_klt = A_left_vec * S * A_right_vec
    
    print(f"\nGravity via KLT = A_left · S · A_right")
    print(f"  = {float(gravity_klt):.10e}")
    
    # Compare with Hodges
    hodges_amp, _ = hodges_npt_mhv_canonical(lambdas, tilde_lambdas, (0, 1))
    print(f"\nHodges amplitude: {float(hodges_amp):.10e}")
    
    if hodges_amp != 0:
        ratio = gravity_klt / hodges_amp
        print(f"Ratio: {float(ratio):.10f}")
    
    # Analyze kernel structure
    print(f"\n{'='*50}")
    print("KERNEL PROPERTIES")
    print(f"{'='*50}")
    
    print(f"det(S) = {S.det()}")
    print(f"rank(S) = {S.rank()}")
    
    # Check symmetry
    S_T = S.transpose()
    is_symmetric = (S == S_T)
    print(f"S symmetric: {is_symmetric}")
    
    if not is_symmetric:
        S_sym = (S + S_T) / 2
        S_antisym = (S - S_T) / 2
        print(f"||S - S^T|| (max entry) = {max(abs(x) for row in (S - S_T).rows() for x in row)}")
    
    return {
        'A_left': A_left_vec,
        'A_right': A_right_vec,
        'S_matrix': S,
        'gravity_klt': gravity_klt,
        'hodges': hodges_amp
    }


if __name__ == "__main__":
    # First verify KLT = Hodges
    verify_klt_equals_hodges(num_samples=5)
    
    # Then analyze geometry
    analyze_klt_geometry(seed=42)

