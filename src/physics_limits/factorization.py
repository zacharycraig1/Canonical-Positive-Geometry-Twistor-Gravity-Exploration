import sys
import os
from sage.all import *

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import hodges_npt_mhv_spinor
from src.physics_limits.bcfw import bcfw_shift_spinors, solve_bcfw_pole, get_channel_s, get_momentum_matrix, decompose_momentum_spinors

def get_amplitude_3pt(lambdas, tilde_lambdas, helicities):
    """
    Compute 3pt gravity amplitude explicitly.
    helicities: list of integers +/- 2.
    """
    h_sum = sum(helicities)
    
    def ang(i, j):
        return lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
        
    def sq(i, j):
        return tilde_lambdas[i][0]*tilde_lambdas[j][1] - tilde_lambdas[i][1]*tilde_lambdas[j][0]

    if h_sum == -2: # MHV (--+)
        minus = [i for i, h in enumerate(helicities) if h == -2]
        plus = [i for i, h in enumerate(helicities) if h == 2]
        if len(minus) != 2: return QQ(0)
        a, b = minus
        c = plus[0]
        num = ang(a, b)**6
        den = (ang(a, c) * ang(b, c))**2
        if den == 0: return QQ(0)
        return num / den
        
    elif h_sum == 2: # Anti-MHV (++-)
        plus = [i for i, h in enumerate(helicities) if h == 2]
        minus = [i for i, h in enumerate(helicities) if h == -2]
        if len(plus) != 2: return QQ(0)
        a, b = plus
        c = minus[0]
        num = sq(a, b)**6
        den = (sq(a, c) * sq(b, c))**2
        if den == 0: return QQ(0)
        return num / den
        
    return QQ(0)

def get_amplitude_mhv(lambdas, tilde_lambdas, negative_indices):
    n = len(lambdas)
    if n == 3:
        helicities = [2]*n
        for i in negative_indices:
            helicities[i] = -2
        return get_amplitude_3pt(lambdas, tilde_lambdas, helicities)

    # Hodges wrapper
    # Try different deletion sets if one fails
    for delete in [(0, 1, 2), (n-3, n-2, n-1), (0, 2, 4)]:
        val, status = hodges_npt_mhv_spinor(lambdas, tilde_lambdas, neg=negative_indices, delete=delete)
        if status == "ok":
            return val
            
    return None

def get_amplitude_anti_mhv(lambdas, tilde_lambdas, positive_indices):
    # Swap L and Lt
    return get_amplitude_mhv(tilde_lambdas, lambdas, positive_indices)

def factorization_check_n6_channel13(seed=42):
    print(f"Running Factorization Check for N=6, Channel s_{{13}} (Seed {seed})...")
    
    # 1. Setup Kinematics (MHV: 0-, 1- others +)
    n = 6
    lambdas, tildes = sample_spinors_from_twistor(seed=seed, n=n)
    
    channel = [1, 3]
    
    # Shift: Need to shift one in channel, one out.
    # 1 is in. 0 is out.
    # Shift 0, 1.
    shift_a, shift_b = 0, 1
    
    z_star = solve_bcfw_pole(lambdas, tildes, shift_a, shift_b, channel)
    print(f"z_star: {z_star}")
    
    if z_star is None:
        print("Error: Could not find pole.")
        return

    # 3. Residue
    epsilon = QQ(1)/QQ(10000000)
    z_probe = z_star + epsilon
    L_probe, Lt_probe = bcfw_shift_spinors(lambdas, tildes, shift_a, shift_b, z_probe)
    
    M6_probe = get_amplitude_mhv(L_probe, Lt_probe, negative_indices=(0, 1))
    s_probe = get_channel_s(L_probe, Lt_probe, channel)
    
    residue_numeric = M6_probe * s_probe
    print(f"Computed Residue (Numeric Limit): {residue_numeric.n()}")
    
    # 4. Exact Factorization
    L_star, Lt_star = bcfw_shift_spinors(lambdas, tildes, shift_a, shift_b, z_star)
    P_mat = get_momentum_matrix(L_star, Lt_star, channel)
    lam_P, lt_P = decompose_momentum_spinors(P_mat)
    
    # P = p1 + p3.
    # R has {1, 3, -P}. P flows out of R?
    # Usually residue is sum M_L(P) M_R(-P).
    # If channel is [1, 3], let's call this the "Right" side.
    # P_channel = p1 + p3.
    # R: {1, 3, -P_channel}. Sum p = 0.
    # L: {0, 2, 4, 5, P_channel}. Sum p = 0.
    
    lam_minusP = lam_P
    lt_minusP = -lt_P
    
    # Sum over h = +/- 2.
    total_residue_exact = QQ(0)
    
    # 1. h = -2 (P is minus)
    # R (1, 3, -P): 1 is -, 3 is +. -P is + (since P is -).
    # Wait. Helicities:
    # If P is (-), then -P is (+).
    # R has {1-, 3+, (-P)+}. 1 minus.
    # 3pt with 1 minus is Anti-MHV (++-). Correct.
    # L (0, 2, 4, 5, P): 0-, P-. 2 minus. MHV.
    
    print("\n  Term h=-2 (P is minus):")
    L_lambdas_m = [L_star[0], L_star[2], L_star[4], L_star[5], lam_P]
    L_tildes_m = [Lt_star[0], Lt_star[2], Lt_star[4], Lt_star[5], lt_P]
    M_L_m = get_amplitude_mhv(L_lambdas_m, L_tildes_m, negative_indices=(0, 4))
    print(f"    M_L (MHV): {M_L_m}")
    
    R_lambdas_m = [L_star[1], L_star[3], lam_minusP]
    R_tildes_m = [Lt_star[1], Lt_star[3], lt_minusP]
    M_R_m = get_amplitude_anti_mhv(R_lambdas_m, R_tildes_m, positive_indices=(1, 2))
    print(f"    M_R (Anti-MHV): {M_R_m}")
    
    if M_L_m is not None and M_R_m is not None:
        term_m = M_L_m * M_R_m
        total_residue_exact += term_m
        
    # 2. h = +2 (P is plus)
    # P is +. -P is -.
    # R has {1-, 3+, -P-}. 2 minus. MHV.
    # L has {0-, P+}. 1 minus. Anti-MHV (for 5pt? No. 5pt Anti-MHV has 3 minus. 1 minus is zero).
    # So this term should be zero.
    
    print("\n  Term h=+2 (P is plus):")
    L_lambdas_p = [L_star[0], L_star[2], L_star[4], L_star[5], lam_P]
    L_tildes_p = [Lt_star[0], Lt_star[2], Lt_star[4], Lt_star[5], lt_P]
    # L is Anti-MHV? No, 1 minus.
    # Let's try to compute MHV anyway? It will be 0.
    # Or Anti-MHV?
    # If L is zero, skip.
    print(f"    L (1 minus) -> Zero")
    
    # But R is MHV.
    R_lambdas_p = [L_star[1], L_star[3], lam_minusP]
    R_tildes_p = [Lt_star[1], Lt_star[3], lt_minusP]
    # Negatives: 1 (index 0), -P (index 2).
    M_R_p = get_amplitude_mhv(R_lambdas_p, R_tildes_p, negative_indices=(0, 2))
    print(f"    M_R (MHV): {M_R_p}")
    
    print(f"  Total Exact: {total_residue_exact.n()}")
    
    ratio = residue_numeric / total_residue_exact
    print(f"Ratio (Numeric/Exact): {ratio.n()}")
    
    if abs(ratio - 1) < 1e-3:
        print("SUCCESS: Factorization confirmed.")
    elif abs(ratio + 1) < 1e-3:
        print("SUCCESS: Factorization confirmed (with sign -1).")
    else:
        print("FAILURE: Mismatch.")

if __name__ == "__main__":
    for i in range(3):
        print(f"\n--- Seed {i} ---")
        try:
            factorization_check_n6_channel13(seed=i)
        except Exception as e:
            print(f"Error with seed {i}: {e}")
            import traceback
            traceback.print_exc()
