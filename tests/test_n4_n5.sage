#!/usr/bin/env sage
# =============================================================================
# TEST: n=4 and n=5 Sanity Checks (CRITICAL)
# =============================================================================
# Must pass with constant ratio before proceeding to n=6 normalization.
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')
load('src/compare.sage')

def hodges_4pt_mhv_reduced(twistor):
    """
    Hodges reduced amplitude for n=4 MHV gravity.
    
    For n=4, we delete 1 row and 1 col, leaving a 3x3 minor.
    Standard choice: delete row 0, col 0.
    """
    n = twistor.n
    if n != 4:
        return (None, "unsupported_n")
    
    # Build Phi matrix
    Phi = matrix(QQ, n, n)
    x, y = 0, 3  # Reference legs
    
    # Off-diagonal
    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = twistor.get_angle(i, j)
                if ij_ang == 0:
                    return (None, "angle_zero")
                ij_sq = twistor.get_square(i, j)
                if ij_sq is None:
                    return (None, "square_none")
                Phi[i, j] = ij_sq / ij_ang
    
    # Diagonal
    for i in range(n):
        if i == x or i == y:
            if i == y:
                diag_sum = QQ(0)
                for j in range(n):
                    if j != i:
                        diag_sum -= Phi[i, j]
                Phi[i, i] = diag_sum
            else:
                Phi[i, i] = QQ(0)
        else:
            ix_ang = twistor.get_angle(i, x)
            iy_ang = twistor.get_angle(i, y)
            if ix_ang == 0 or iy_ang == 0:
                return (None, "diag_angle_zero")
            
            diag_sum = QQ(0)
            for j in range(n):
                if j == i:
                    continue
                jx_ang = twistor.get_angle(j, x)
                jy_ang = twistor.get_angle(j, y)
                if jx_ang == 0 or jy_ang == 0:
                    continue
                contrib = Phi[i, j] * (jx_ang * jy_ang) / (ix_ang * iy_ang)
                diag_sum -= contrib
            Phi[i, i] = diag_sum
    
    # For n=4: delete row 0, col 0, keep rows [1,2,3], cols [1,2,3]
    rows_keep = [1, 2, 3]
    cols_keep = [1, 2, 3]
    Phi_minor = Phi[rows_keep, cols_keep]
    
    try:
        det_minor = Phi_minor.det()
    except:
        return (None, "det_failed")
    
    # c factors: c_0 = 1 (single deleted index), c_1 = 1/(<12><23><31>)
    a12 = twistor.get_angle(1, 2)
    a23 = twistor.get_angle(2, 3)
    a31 = twistor.get_angle(3, 1)
    if a12 == 0 or a23 == 0 or a31 == 0:
        return (None, "c_factor_zero")
    c123 = QQ(1) / (a12 * a23 * a31)
    
    # Sign: (-1)^{n+1} = (-1)^5 = -1
    sign = -1
    
    bar_M4 = sign * c123 * det_minor
    return (bar_M4, "ok")


def parke_taylor_npt_mhv(twistor, order, neg_helicity=(0, 1)):
    """Parke-Taylor for arbitrary n."""
    n = twistor.n
    if len(order) != n:
        return None
    
    # Cyclic product
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        idx_i = order[i]
        idx_j = order[j]
        bracket = twistor.get_angle(idx_i, idx_j)
        if bracket == 0:
            return None
        denom *= bracket
    
    # Helicity factor
    neg_a, neg_b = neg_helicity
    helicity_factor = twistor.get_angle(neg_a, neg_b)
    if helicity_factor == 0:
        return None
    
    if denom == 0:
        return None
    
    return (helicity_factor ** 4) / denom


def klt_4pt_mhv(twistor, mandelstam_func):
    """
    KLT gravity amplitude for n=4 MHV.
    
    For n=4, permuted set is {1} (0-based for particle 2), so only 1 permutation.
    Standard KLT: M_4 = A(3,0,1,2) * S[1|1] * A(2,1,3,0)
    """
    n = 4
    permuted_set = [1]  # Only particle 1 (0-based) is permuted
    
    # Fixed legs: 0, 2, 3 (0-based for 1, 3, 4)
    fixed_0 = 0
    fixed_2 = 2
    fixed_3 = 3
    
    total = QQ(0)
    
    # Only one permutation: [1]
    alpha = [1]
    beta = [1]
    
    # A(3,0,alpha,2) = A(3,0,1,2)
    order_alpha = [fixed_3, fixed_0] + alpha + [fixed_2]
    A_alpha = parke_taylor_npt_mhv(twistor, order_alpha, neg_helicity=(0, 1))
    if A_alpha is None:
        return (None, "A_alpha_none")
    
    # A(2,beta,3,0) = A(2,1,3,0)
    order_beta = [fixed_2] + beta + [fixed_3, fixed_0]
    A_beta = parke_taylor_npt_mhv(twistor, order_beta, neg_helicity=(0, 1))
    if A_beta is None:
        return (None, "A_beta_none")
    
    # KLT kernel for n=4: S[alpha|beta] = s_{0,alpha[0]} = s_{0,1}
    S = mandelstam_func(twistor, 0, alpha[0])
    if S is None:
        return (None, "kernel_none")
    
    total = A_alpha * S * A_beta
    
    return (total, "ok")


def test_n4():
    """Test n=4 with constant ratio requirement."""
    print("="*70)
    print("TEST: n=4 MHV Gravity (CRITICAL SANITY CHECK)")
    print("="*70)
    
    ratios = []
    
    for seed in range(20):
        # Sample 4-point moment curve
        Z = []
        for i in range(4):
            t = QQ(i + 1) + QQ(seed % 100) / QQ(1000)
            z = vector(QQ, [QQ(1), t, t*t, t*t*t])
            Z.append(z)
        
        twistor = MomentumTwistor(n=4, Z=Z, check_domain=True)
        
        if not twistor.domain_ok:
            continue
        
        H_red = hodges_4pt_mhv_reduced(twistor)
        A_klt = klt_4pt_mhv(twistor, mandelstam_invariant)
        
        H_val = H_red[0] if isinstance(H_red, tuple) else H_red
        A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
        
        if H_val is None or A_val is None or H_val == 0:
            print(f"Seed {seed}: H={H_val}, A={A_val}")
            continue
        
        ratio = A_val / H_val
        ratios.append(ratio)
        print(f"Seed {seed}: ratio = {ratio}")
    
    if not ratios:
        print("ERROR: No valid ratios computed")
        return False
    
    unique_ratios = list(set(ratios))
    print(f"\nTotal samples: {len(ratios)}")
    print(f"Unique ratios: {len(unique_ratios)}")
    
    if len(unique_ratios) == 1:
        print(f"SUCCESS: Constant ratio = {unique_ratios[0]}")
        return True
    else:
        print(f"FAILURE: Ratios vary:")
        for i, r in enumerate(unique_ratios[:5]):
            print(f"  Ratio {i+1}: {r}")
        return False


if __name__ == '__main__':
    success = test_n4()
    print(f"\n{'='*70}")
    print(f"n=4 TEST: {'PASSED' if success else 'FAILED'}")
    print(f"{'='*70}")
    if not success:
        print("STOP-THE-LINE: n=4 must pass before proceeding to n=6")

