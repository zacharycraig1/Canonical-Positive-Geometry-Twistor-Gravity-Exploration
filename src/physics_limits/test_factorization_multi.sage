import sys
import os
from sage.all import *

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import hodges_npt_mhv_spinor, ang_bracket, sq_bracket
# Use the canonical Hodges
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical
from src.physics_limits.bcfw import bcfw_shift_spinors, solve_bcfw_pole, get_channel_s, get_momentum_matrix, decompose_momentum_spinors

def get_amplitude_3pt(lambdas, tilde_lambdas, helicities):
    # Same as before, but ensure signs are correct.
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

def get_amplitude_mhv_canonical(lambdas, tilde_lambdas, negative_indices):
    n = len(lambdas)
    if n == 3:
        helicities = [2]*n
        for i in negative_indices:
            helicities[i] = -2
        return get_amplitude_3pt(lambdas, tilde_lambdas, helicities)
        
    val, status = hodges_npt_mhv_canonical(lambdas, tilde_lambdas, negative_indices)
    if status != "ok":
        return None
    return val

def get_amplitude_anti_mhv_canonical(lambdas, tilde_lambdas, positive_indices):
    # Swap L and Lt
    return get_amplitude_mhv_canonical(tilde_lambdas, lambdas, positive_indices)

def check_factorization_multi(seed=42):
    print(f"M2: Checking Multi-Particle Factorization (Seed {seed})...")
    
    # 1. Kinematics N=6
    n = 6
    lambdas, tildes = sample_spinors_from_twistor(seed=seed, n=n)
    
    # Channel [0, 1, 2].
    channel = [0, 1, 2]
    # Shift 0, 3. (0 in, 3 out).
    shift_a, shift_b = 0, 3
    
    print(f"  Channel: {channel}")
    print(f"  Shift: {shift_a}, {shift_b}")
    
    # 2. Solve Pole
    z_star = solve_bcfw_pole(lambdas, tildes, shift_a, shift_b, channel)
    if z_star is None:
        print("  Error: No pole found.")
        return
    print(f"  z_star: {z_star}")
    
    # 3. Numeric Residue Extraction (Multiple Epsilon)
    residues = []
    epsilons = [QQ(1)/10000, QQ(1)/100000, QQ(1)/1000000]
    
    print("  Computing Numeric Residue limit...")
    for eps in epsilons:
        z_probe = z_star + eps
        L_probe, Lt_probe = bcfw_shift_spinors(lambdas, tildes, shift_a, shift_b, z_probe)
        
        M6 = get_amplitude_mhv_canonical(L_probe, Lt_probe, negative_indices=(0, 1))
        
        # We need s(z).
        # Note: We calibrated s normalization to match physical <ij>[ji].
        # But get_channel_s uses det(P).
        # If det(P) = - physical_s (for 3pt), we might have a sign.
        # But residues are defined on s(z).
        # Res = Lim s(z) M(z).
        # If we stick to det(P) for both pole condition and residue definition, it's self-consistent.
        # The normalization of s matters for comparison with M_L * M_R.
        # M_L M_R usually assumes "standard" propagator 1/P^2.
        # If P^2 (det) = +/- P^2 (phys), we might pick up a sign.
        
        s_val = get_channel_s(L_probe, Lt_probe, channel)
        res = M6 * s_val
        residues.append(res)
        # print(f"    eps={float(eps):.1e}: Res={res.n()}")
        
    # Check stability
    res_final = residues[-1]
    diff = abs(residues[-1] - residues[-2])
    print(f"  Convergence Diff: {diff.n()}")
    
    # If residue is small, it might be zero.
    if abs(res_final) < 1e-6:
        print("  PASS: Residue is effectively zero (consistent with selection rules).")
    else:
        print("  FAIL: Numeric residue is non-zero.")
    L_star, Lt_star = bcfw_shift_spinors(lambdas, tildes, shift_a, shift_b, z_star)
    P_mat = get_momentum_matrix(L_star, Lt_star, channel)
    
    # Decompose P.
    # Enforce BCFW consistency:
    # If shift is on 0 (in L), 3 (in R).
    # L side has {0, 1, 2, -P}. (Outgoing P convention usually P_L = -P_R).
    # Wait. P_channel = p0 + p1 + p2.
    # P flows out of L? No, p0+p1+p2 are outgoing. So P is incoming to L?
    # Standard: sum p_out = 0.
    # p0+p1+p2 + P_internal = 0. => P_internal = -(p0+p1+p2).
    # Then P_internal flows INTO R?
    # R side: p3+p4+p5 + (-P_internal) = 0.
    
    # We define P = p0+p1+p2.
    # L has {0, 1, 2, -P}.
    # R has {3, 4, 5, P}.
    
    # BCFW Internal Spinors:
    # Shifted leg 0 is in L.
    # lambda_P = lambda_0(z_star).
    # tilde_lambda_P is fixed by P = lambda_P tilde_lambda_P.
    
    lam_P, lt_P = decompose_momentum_spinors(P_mat, ref_lambda=L_star[0])
    
    # L side: {0, 1, 2, -P}.
    # -P has spinors {lam_P, -lt_P}.
    # R side: {3, 4, 5, P}.
    # P has spinors {lam_P, lt_P}.
    
    # Sum over h = +/- 2.
    # N=6 MHV (0-, 1-).
    # L side: 0-, 1-, 2+, -P?.
    # R side: 3+, 4+, 5+, P?.
    
    # h = +2 (Internal P is +):
    # R side: {3+, 4+, 5+, P+}. All plus -> 0 (Gravity requires 2 minus for MHV? No, Anti-MHV requires 2 plus?)
    # 4pt R: {3, 4, 5, P}.
    # If all +, M4(++++) = 0.
    # So h=+2 (P is +) is zero.
    
    # h = -2 (Internal P is -):
    # R side: {3+, 4+, 5+, P-}. One minus.
    # 4pt with 1 minus? M4(+++-) = 0?
    # Yes, MHV needs 2 minus. Anti-MHV needs 2 plus.
    # 4pt: MHV (--++), Anti-MHV (++--).
    # (+++-) is zero.
    
    # Wait. If both h terms are zero, residue is zero.
    # But numeric residue was non-zero?
    
    # Let's re-check helicity counting.
    # N=6 (0-, 1-).
    # L: {0-, 1-, 2+, -P}.
    # R: {3+, 4+, 5+, P}.
    
    # Term 1: P is h=-2.
    # R has P-. So {3+, 4+, 5+, -}. One minus.
    # Is M4(+++-) zero?
    # Yes, for Gravity.
    
    # Term 2: P is h=+2.
    # R has P+. So {3+, 4+, 5+, +}. All plus. Zero.
    
    # What went wrong?
    # Maybe my channel selection is bad for this helicity config.
    # For s_012 (0, 1 are minus).
    # 0, 1 both in L.
    # Then L has 2 minus from external, plus P.
    # If P is +, L has 2 minus. M4(--++). Non-zero.
    # But R needs to be non-zero. R has {3, 4, 5, P}. All + if P is +.
    # If P is -, L has 3 minus? M4(---+) is zero.
    # And R has 1 minus. M4(+++-) is zero.
    
    # CONCLUSION: s_012 residue MUST be zero for MHV(0,1).
    # If numeric residue is non-zero, we have a bug or I misunderstood helicity.
    # Let's check numeric residue.
    
    if abs(res_final) > 1e-10:
        print("  FAIL: Numeric residue is non-zero, but helicity arguments imply zero.")
        print("        Does s_012 have a pole? Yes, propagator exists.")
        print("        But numerator must vanish.")
    else:
        print("  PASS: Numeric residue is zero (consistent with selection rules).")
        
    # We need a NON-ZERO channel to test factorization.
    # Split the minuses.
    # Put 0 in L, 1 in R.
    # Channel {0, 2, 3}.
    # L: {0-, 2+, 3+, -P}.
    # R: {1-, 4+, 5+, P}.
    
    # If P is +:
    # L has {-, +, +, -}. (Using -P is -). Two minus. MHV!
    # R has {-, +, +, +}. One minus. Zero.
    
    # If P is -:
    # L has {-, +, +, +}. One minus. Zero.
    # R has {-, +, +, -}. Two minus. MHV!
    
    # So for {0, 2, 3}, only P=(h=-2) contributes?
    # Wait.
    # If P is -, then P particle has helicity -2.
    # In R: {..., P}. So P is -2. R has {1-, P-}. Two minus. MHV.
    # In L: {..., -P}. -P has helicity +2.
    # L has {0-, (-P)+}. One minus. Zero.
    
    # If P is + (helicity +2).
    # R: {1-, P+}. One minus. Zero.
    # L: {0-, (-P)-}. Two minus. MHV.
    
    # So we need one MHV and one Anti-MHV?
    # 4pt Anti-MHV is (++--).
    # But here we have 4pt sub-amplitudes.
    # 4pt is special. MHV = Anti-MHV.
    # M4(--++) is MHV. M4(++--) is Anti-MHV.
    # They are the same class.
    
    # So we need a split where both sides are MHV (or Anti-MHV).
    # Try Channel {0, 2, 3}.
    # L: {0-, 2+, 3+, -P}.
    # R: {1-, 4+, 5+, P}.
    
    # Case P is + (h=2):
    # R: {1-, 4+, 5+, P+}. 1 minus. Zero.
    # L: {0-, 2+, 3+, (-P)-}. 2 minus. MHV.
    # Product = 0 * MHV = 0.
    
    # Case P is - (h=-2):
    # R: {1-, 4+, 5+, P-}. 2 minus. MHV.
    # L: {0-, 2+, 3+, (-P)+}. 1 minus. Zero.
    # Product = 0.
    
    # Is there ANY channel that is non-zero?
    # Maybe 3-particle channel is bad?
    # For N=6, L=4pt, R=4pt.
    # Total minus = 2.
    # L needs 2 minus. R needs 2 minus. (Because P and -P carry opposite helicity).
    # If P is +, -P is -.
    # Total minus in L+R = (minuses in L) + (minuses in R).
    # L has (-P)-. R has P+.
    # Sum = (minus_L_ext + 1) + (minus_R_ext).
    # We want Sum = 4 (for two MHV sub-amps).
    # minus_L_ext + 1 + minus_R_ext = 4.
    # minus_tot_ext + 1 = 4 => minus_tot_ext = 3.
    # But we only have 2 minus external!
    
    # This implies N=6 MHV (2 minus) NEVER factorizes into two MHV 4-point amplitudes.
    # It must involve Anti-MHV?
    # 4pt Anti-MHV has 2 plus (so 2 minus).
    # Wait. M4 is N=4. MHV has 2 minus. Anti-MHV has 2 minus.
    # For N=4, MHV = Anti-MHV.
    # So 4pt always has 2 minus.
    
    # So for N=6 MHV to factorize into 4pt x 4pt:
    # We need L (2 minus) and R (2 minus).
    # Total minus count:
    # L: {ext_L} U {-P}.
    # R: {ext_R} U {P}.
    # #minus(L) = #minus(ext_L) + helicity(-P is minus)
    # #minus(R) = #minus(ext_R) + helicity(P is minus)
    
    # If P is +, -P is -.
    # #minus(L) = #minus(ext_L) + 1.
    # #minus(R) = #minus(ext_R) + 0.
    # Total required: 2 + 2 = 4.
    # #minus(ext_L) + 1 + #minus(ext_R) = 4.
    # #minus(ext_L) + #minus(ext_R) = 3.
    # But total external minus is 2.
    # 2 != 3.
    
    # If P is -, -P is +.
    # #minus(L) = #minus(ext_L) + 0.
    # #minus(R) = #minus(ext_R) + 1.
    # Total: #minus(ext_L) + #minus(ext_R) + 1 = 4 => 3.
    # Still 2 != 3.
    
    # Conclusion: N=6 MHV Gravity DOES NOT have poles in 3-particle channels (s_ijk) corresponding to 4pt x 4pt factorization!
    # Those residues are strictly zero.
    # N=6 MHV only has poles in 2-particle channels (s_ij)?
    # Where it splits into 3pt x 5pt.
    
    # Let's check 3pt x 5pt.
    # L=3pt, R=5pt.
    # L needs 2 minus (MHV) or 1 minus (Anti-MHV).
    # R needs 2 minus (MHV).
    
    # Case 1: L is Anti-MHV (1 minus). R is MHV (2 minus).
    # P is -. -P is +.
    # #minus(L) = #minus(ext_L) + 0 = 1 => #minus(ext_L) = 1.
    # #minus(R) = #minus(ext_R) + 1 = 2 => #minus(ext_R) = 1.
    # Total ext minus = 1 + 1 = 2.
    # This works!
    # So we need a 3-particle channel that is actually a 2-particle channel?
    # No, s_ij is a pole. s_ijk (3-particle) for n=6 means 3 particles on one side, 3 on other.
    # So 3pt x 3pt?
    # No, n=6 split into 3+3?
    # Yes, L=3+P (4pt), R=3-P (4pt).
    # Wait.
    # 6 particles. 3 left, 3 right.
    # L has 3 external + 1 internal = 4 particles.
    # R has 3 external + 1 internal = 4 particles.
    # This is the 4pt x 4pt case which we proved is zero.
    
    # So s_ijk poles are absent in MHV6?
    # If true, then my multi-particle factorization test should yield 0.
    # And we should verify that.
    
    print("  Analysis: N=6 MHV should have NO 3-particle poles (s_ijk).")
    print(f"  Numeric result: {res_final.n()}")
    
    if abs(res_final) < 1e-6:
        print("  PASS: Residue is zero as expected.")
    else:
        print("  FAIL: Residue is non-zero! (Unexpected for MHV6)")
        
if __name__ == "__main__":
    check_factorization_multi()

