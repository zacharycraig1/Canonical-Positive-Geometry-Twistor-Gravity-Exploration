# src/kinematic_associahedron/test_4pt_klt.sage
"""
Test KLT at 4-point
====================

For n=4, there's only ONE scattering equation solution, making it the
simplest case to debug the KLT relation.

4-point MHV gravity amplitude:
    M_4 = <12>^8 [34]^8 / (s_12 s_23 s_13)

Or equivalently:
    M_4 = s_12^3 s_23 s_13 × |<12>[12]|^4 / (s_12 s_23 s_13)
        = s_12^2 × <12>^4 [12]^4 × ... (complex)

Simpler form using Mandelstam:
    M_4^MHV = <12>^4 [34]^4 × s_12 / s_23

Wait, let's use the correct form.

For 4-point MHV gravity (1-, 2-, 3+, 4+):
    M_4 = <12>^8 [34]^4 / (s s t)

where s = s_12 = s_34, t = s_23 = s_14, u = s_13 = s_24, and s + t + u = 0.

Actually the simplest is:
    M_4 = <12>^4 [34]^4 × s / t

Using momentum conservation and on-shell conditions.

Let's verify this against KLT.

For n=4 KLT:
    M_4 = A_4(1,2,3,4) × s_12 × A_4(1,2,4,3)

where A_4 is the YM amplitude.

A_4^MHV(1-,2-,3+,4+) = <12>^4 / (<12><23><34><41>)
"""

from sage.all import *
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.kinematics.spinors import SpinorKinematics


def angle_bracket(lambdas, i, j):
    """Compute <ij>."""
    return lambdas[i][0] * lambdas[j][1] - lambdas[i][1] * lambdas[j][0]


def square_bracket(tilde_lambdas, i, j):
    """Compute [ij]."""
    return tilde_lambdas[i][0] * tilde_lambdas[j][1] - tilde_lambdas[i][1] * tilde_lambdas[j][0]


def mandelstam(lambdas, tilde_lambdas, i, j):
    """Compute s_ij = <ij>[ji]."""
    # Note: s_ij = <ij>[ij] in some conventions, <ij>[ji] in others
    # Standard: s_ij = (p_i + p_j)^2 = 2 p_i.p_j = <ij>[ji]
    return angle_bracket(lambdas, i, j) * square_bracket(tilde_lambdas, j, i)


def ym_amplitude_4pt(lambdas, order, neg_helicity=(0, 1)):
    """
    Compute 4-point MHV Yang-Mills amplitude.
    
    A_4 = <ab>^4 / (<order[0] order[1]> <order[1] order[2]> <order[2] order[3]> <order[3] order[0]>)
    """
    a, b = neg_helicity
    numerator = angle_bracket(lambdas, a, b) ** 4
    
    denom = QQ(1)
    n = len(order)
    for i in range(n):
        j = (i + 1) % n
        bracket = angle_bracket(lambdas, order[i], order[j])
        if bracket == 0:
            return None
        denom *= bracket
    
    return numerator / denom


def gravity_4pt_direct(lambdas, tilde_lambdas, neg_helicity=(0, 1)):
    """
    Compute 4-point MHV gravity amplitude directly.
    
    From the directive (line 564):
    For MHV 4-point: M_4 = ⟨12⟩^4 [34]^4 / (s_12 s_23)
    
    where particles 1,2 have negative helicity and 3,4 have positive helicity.
    """
    a, b = neg_helicity
    c, d = 2, 3  # Positive helicity particles
    
    ang_ab = angle_bracket(lambdas, a, b)
    sq_cd = square_bracket(tilde_lambdas, c, d)
    
    s_12 = mandelstam(lambdas, tilde_lambdas, 0, 1)
    s_23 = mandelstam(lambdas, tilde_lambdas, 1, 2)
    
    if s_12 == 0 or s_23 == 0:
        return None
    
    # M_4^MHV = <12>^4 [34]^4 / (s_12 s_23)
    return (ang_ab**4) * (sq_cd**4) / (s_12 * s_23)


def gravity_4pt_klt(lambdas, tilde_lambdas, neg_helicity=(0, 1)):
    """
    Compute 4-point gravity via KLT.
    
    For n=4, the KLT relation is:
    M_4 = A_4(1,2,3,4) × s_12 × A_4(1,2,4,3)
    
    Note: Only ONE permutation of the middle particles (just {3} vs {4}).
    The KLT kernel for n=4 is just s_12.
    """
    # For n=4, permuted set is {2} (indices 1 in 0-indexed, or particle 3 in 1-indexed)
    # Wait, let me reconsider.
    
    # Standard KLT for n=4:
    # Fixed legs: 1 and 4 (indices 0 and 3)
    # Permuted set: {2, 3} (indices 1, 2)
    # But for n=4, there's only ONE permutation of the middle that matters
    
    # Actually, for n=4:
    # M_4 = A_4(1,2,3,4) × s_12 × A_4(1,2,4,3)
    # 
    # With ordering (0,1,2,3) for first, (0,1,3,2) for second
    
    order_left = [0, 1, 2, 3]  # (1,2,3,4)
    order_right = [0, 1, 3, 2]  # (1,2,4,3) - swapped last two
    
    A_left = ym_amplitude_4pt(lambdas, order_left, neg_helicity)
    A_right = ym_amplitude_4pt(lambdas, order_right, neg_helicity)
    
    if A_left is None or A_right is None:
        return None
    
    s_12 = mandelstam(lambdas, tilde_lambdas, 0, 1)
    
    # KLT: M_4 = A_left × s_12 × A_right
    return A_left * s_12 * A_right


def generate_4pt_kinematics(seed=42):
    """
    Generate 4-point kinematics satisfying momentum conservation.
    
    For 4 particles: p_1 + p_2 + p_3 + p_4 = 0
    """
    set_random_seed(seed)
    
    # Use momentum twistor approach for guaranteed conservation
    # For n=4, generate 4 random twistors
    Z = []
    for _ in range(4):
        z = vector(QQ, [QQ.random_element(num_bound=5, den_bound=3) for _ in range(4)])
        Z.append(z)
    
    # Extract lambdas
    lambdas = [vector(QQ, [z[0], z[1]]) for z in Z]
    
    # Compute tilde_lambdas from twistors using the formula
    def get_angle(i, j):
        return Z[i][0]*Z[j][1] - Z[i][1]*Z[j][0]
    
    n = 4
    tilde_lambdas = []
    
    for i in range(n):
        im1 = (i - 1) % n
        ip1 = (i + 1) % n
        
        mu_i = vector(QQ, [Z[i][2], Z[i][3]])
        mu_im1 = vector(QQ, [Z[im1][2], Z[im1][3]])
        mu_ip1 = vector(QQ, [Z[ip1][2], Z[ip1][3]])
        
        ang_im1_i = get_angle(im1, i)
        ang_i_ip1 = get_angle(i, ip1)
        ang_ip1_im1 = get_angle(ip1, im1)
        
        denom = ang_im1_i * ang_i_ip1
        if denom == 0:
            # Degenerate configuration
            return None, None
        
        numerator = ang_i_ip1 * mu_im1 + ang_ip1_im1 * mu_i + ang_im1_i * mu_ip1
        tilde_lambda = numerator / denom
        tilde_lambdas.append(tilde_lambda)
    
    return lambdas, tilde_lambdas


def test_4pt_klt_vs_direct():
    """Test KLT against direct formula for 4-point."""
    print("\n" + "="*70)
    print("TESTING 4-POINT KLT vs DIRECT FORMULA")
    print("="*70)
    
    num_tests = 10
    matches = 0
    
    for seed in range(42, 42 + num_tests):
        lambdas, tilde_lambdas = generate_4pt_kinematics(seed)
        
        if lambdas is None:
            print(f"Seed {seed}: Degenerate kinematics, skipping")
            continue
        
        # Compute via KLT
        M_klt = gravity_4pt_klt(lambdas, tilde_lambdas, (0, 1))
        
        if M_klt is None:
            print(f"Seed {seed}: KLT computation failed")
            continue
        
        # Compute direct
        M_direct = gravity_4pt_direct(lambdas, tilde_lambdas, (0, 1))
        
        if M_direct is None:
            print(f"Seed {seed}: Direct computation failed")
            continue
        
        print(f"\nSeed {seed}:")
        print(f"  M_4 (KLT):    {float(M_klt):.10e}")
        print(f"  M_4 (direct): {float(M_direct):.10e}")
        
        if M_direct != 0:
            ratio = M_klt / M_direct
            print(f"  Ratio: {float(ratio):.10f}")
            
            if ratio == 1:
                print(f"  ✓ EXACT MATCH")
                matches += 1
            else:
                print(f"  ✗ Mismatch")
        else:
            print(f"  Direct is zero, cannot compare")
    
    print(f"\n{'='*70}")
    print(f"SUMMARY: {matches}/{num_tests} exact matches")
    print(f"{'='*70}")
    
    return matches


def analyze_4pt_structure():
    """Analyze the structure of 4-point amplitudes."""
    print("\n" + "="*70)
    print("ANALYZING 4-POINT AMPLITUDE STRUCTURE")
    print("="*70)
    
    seed = 42
    lambdas, tilde_lambdas = generate_4pt_kinematics(seed)
    
    if lambdas is None:
        print("Degenerate kinematics")
        return
    
    print("\nSpinor brackets:")
    for i in range(4):
        for j in range(i+1, 4):
            ang = angle_bracket(lambdas, i, j)
            sq = square_bracket(tilde_lambdas, i, j)
            s = mandelstam(lambdas, tilde_lambdas, i, j)
            print(f"  <{i}{j}> = {float(ang):10.4f},  [{i}{j}] = {float(sq):10.4f},  s_{i}{j} = {float(s):10.4f}")
    
    # Verify momentum conservation: s + t + u = 0
    s = mandelstam(lambdas, tilde_lambdas, 0, 1)  # s_12
    t = mandelstam(lambdas, tilde_lambdas, 1, 2)  # s_23
    u = mandelstam(lambdas, tilde_lambdas, 0, 2)  # s_13
    
    print(f"\nMomentum conservation check:")
    print(f"  s_12 + s_23 + s_13 = {float(s + t + u):.10f} (should be 0)")
    
    # YM amplitudes
    print("\nYang-Mills amplitudes:")
    A_1234 = ym_amplitude_4pt(lambdas, [0,1,2,3], (0,1))
    A_1243 = ym_amplitude_4pt(lambdas, [0,1,3,2], (0,1))
    
    print(f"  A_4(1,2,3,4) = {float(A_1234):.10e}")
    print(f"  A_4(1,2,4,3) = {float(A_1243):.10e}")
    
    # KLT gravity
    M_klt = A_1234 * s * A_1243
    print(f"\nKLT gravity:")
    print(f"  M_4 = A(1234) × s_12 × A(1243)")
    print(f"      = {float(A_1234):.6e} × {float(s):.6f} × {float(A_1243):.6e}")
    print(f"      = {float(M_klt):.10e}")
    
    # Direct formula
    M_direct = gravity_4pt_direct(lambdas, tilde_lambdas, (0,1))
    print(f"\nDirect formula:")
    print(f"  M_4 = <12>^4 [34]^4 s_12 / s_23")
    ang_12 = angle_bracket(lambdas, 0, 1)
    sq_34 = square_bracket(tilde_lambdas, 2, 3)
    print(f"      = {float(ang_12**4):.6e} × {float(sq_34**4):.6e} × {float(s):.6f} / {float(t):.6f}")
    print(f"      = {float(M_direct):.10e}")
    
    if M_direct != 0:
        print(f"\nRatio KLT/Direct: {float(M_klt/M_direct):.10f}")


def gravity_4pt_hodges(lambdas, tilde_lambdas, neg_helicity=(0, 1)):
    """
    Compute 4-point gravity using Hodges formula.
    
    Import from the working Hodges module.
    """
    from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical
    return hodges_npt_mhv_canonical(lambdas, tilde_lambdas, neg_helicity)


def test_all_4pt_formulas():
    """Compare all 4-point gravity formulas."""
    print("\n" + "="*70)
    print("COMPARING ALL 4-POINT GRAVITY FORMULAS")
    print("="*70)
    
    for seed in range(43, 53):
        lambdas, tilde_lambdas = generate_4pt_kinematics(seed)
        
        if lambdas is None:
            continue
        
        print(f"\n--- Seed {seed} ---")
        
        # KLT
        M_klt = gravity_4pt_klt(lambdas, tilde_lambdas, (0, 1))
        
        # Direct
        M_direct = gravity_4pt_direct(lambdas, tilde_lambdas, (0, 1))
        
        # Hodges
        M_hodges, status = gravity_4pt_hodges(lambdas, tilde_lambdas, (0, 1))
        
        print(f"  M_4 (KLT):    {float(M_klt) if M_klt else 0:>15.6e}")
        print(f"  M_4 (direct): {float(M_direct) if M_direct else 0:>15.6e}")
        if M_hodges is not None:
            print(f"  M_4 (Hodges): {float(M_hodges):>15.6e} [{status}]")
        else:
            print(f"  M_4 (Hodges): None [{status}]")
        
        # Compare ratios
        if M_klt and M_hodges and M_hodges != 0:
            ratio_klt_hodges = M_klt / M_hodges
            print(f"  Ratio KLT/Hodges:    {float(ratio_klt_hodges):.10f}")
        
        if M_direct and M_hodges and M_hodges != 0:
            ratio_direct_hodges = M_direct / M_hodges
            print(f"  Ratio Direct/Hodges: {float(ratio_direct_hodges):.10f}")


if __name__ == "__main__":
    # analyze_4pt_structure()
    # test_4pt_klt_vs_direct()
    test_all_4pt_formulas()

