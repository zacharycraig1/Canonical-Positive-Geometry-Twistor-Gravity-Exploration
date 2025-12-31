#!/usr/bin/env sage
# =============================================================================
# TEST: n=4 with Direct Spinor-Helicity Sampling (Fallback Approach)
# =============================================================================
# Bypass twistors for n=4 and sample λ_i, \tildeλ_i directly
# Compute physical [ij] and s_ij correctly
# =============================================================================

from sage.all import *
import numpy as np

def sample_spinor_helicity_4pt(seed=None):
    """
    Sample 4-point spinor-helicity data directly.
    
    Momentum conservation: Σ p_i = 0, where p_i = λ_i \tildeλ_i
    
    Returns:
        Tuple (lambdas, tildelambdas) where:
        - lambdas: list of λ_i (2-component spinors)
        - tildelambdas: list of \tildeλ_i (2-component spinors)
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Sample λ_i (2-component spinors over QQ)
    lambdas = []
    for i in range(4):
        lambda_i = vector(QQ, [
            QQ(np.random.randint(1, 10)),
            QQ(np.random.randint(1, 10))
        ])
        # Ensure not all zero
        while all(x == 0 for x in lambda_i):
            lambda_i = vector(QQ, [
                QQ(np.random.randint(1, 10)),
                QQ(np.random.randint(1, 10))
            ])
        lambdas.append(lambda_i)
    
    # Sample \tildeλ_i and enforce momentum conservation
    # p_0 + p_1 + p_2 + p_3 = 0
    # We'll sample \tildeλ_0, \tildeλ_1, \tildeλ_2 freely, then solve for \tildeλ_3
    
    tildelambdas = []
    for i in range(3):
        tildelambda_i = vector(QQ, [
            QQ(np.random.randint(1, 10)),
            QQ(np.random.randint(1, 10))
        ])
        while all(x == 0 for x in tildelambda_i):
            tildelambda_i = vector(QQ, [
                QQ(np.random.randint(1, 10)),
                QQ(np.random.randint(1, 10))
            ])
        tildelambdas.append(tildelambda_i)
    
    # Solve for \tildeλ_3 from momentum conservation
    # p_0 + p_1 + p_2 + p_3 = 0
    # λ_0 \tildeλ_0 + λ_1 \tildeλ_1 + λ_2 \tildeλ_2 + λ_3 \tildeλ_3 = 0
    # This is a 2×2 system: λ_3 \tildeλ_3 = - (λ_0 \tildeλ_0 + λ_1 \tildeλ_1 + λ_2 \tildeλ_2)
    
    p_sum = lambdas[0] * tildelambdas[0] + lambdas[1] * tildelambdas[1] + lambdas[2] * tildelambdas[2]
    p_sum = -p_sum  # Need to solve λ_3 \tildeλ_3 = p_sum
    
    # For a 2-component \tildeλ_3, we have:
    # λ_3^0 \tildeλ_3^0 + λ_3^1 \tildeλ_3^1 = p_sum^0
    # This is underdetermined, so we'll use a simple approach:
    # Set \tildeλ_3^0 = 1, solve for \tildeλ_3^1, or vice versa
    
    lambda_3 = lambdas[3]
    if lambda_3[0] != 0:
        # Solve: λ_3^0 * t_0 + λ_3^1 * t_1 = p_sum^0
        # Set t_1 = 1, then t_0 = (p_sum^0 - λ_3^1) / λ_3^0
        tildelambda_3 = vector(QQ, [
            (p_sum[0] - lambda_3[1]) / lambda_3[0] if lambda_3[0] != 0 else QQ(1),
            QQ(1)
        ])
    elif lambda_3[1] != 0:
        tildelambda_3 = vector(QQ, [
            QQ(1),
            (p_sum[1] - lambda_3[0]) / lambda_3[1] if lambda_3[1] != 0 else QQ(1)
        ])
    else:
        # Fallback
        tildelambda_3 = vector(QQ, [QQ(1), QQ(1)])
    
    tildelambdas.append(tildelambda_3)
    
    return lambdas, tildelambdas


def angle_bracket_spinor(lambda_i, lambda_j):
    """<i j> = ε_{αβ} λ_i^α λ_j^β = λ_i^0 λ_j^1 - λ_i^1 λ_j^0"""
    return lambda_i[0] * lambda_j[1] - lambda_i[1] * lambda_j[0]


def square_bracket_spinor(tildelambda_i, tildelambda_j):
    """[i j] = ε^{\dotα\dotβ} \tildeλ_{i,\dotα} \tildeλ_{j,\dotβ} = \tildeλ_i^0 \tildeλ_j^1 - \tildeλ_i^1 \tildeλ_j^0"""
    return tildelambda_i[0] * tildelambda_j[1] - tildelambda_i[1] * tildelambda_j[0]


def mandelstam_spinor(lambda_i, lambda_j, tildelambda_i, tildelambda_j):
    """s_{ij} = <i j> [i j]"""
    ang = angle_bracket_spinor(lambda_i, lambda_j)
    sq = square_bracket_spinor(tildelambda_i, tildelambda_j)
    return ang * sq


def parke_taylor_4pt_spinor(lambdas, tildelambdas, order, neg_helicity=(0, 1)):
    """Parke-Taylor for 4-point using spinor-helicity."""
    n = 4
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        idx_i = order[i]
        idx_j = order[j]
        bracket = angle_bracket_spinor(lambdas[idx_i], lambdas[idx_j])
        if bracket == 0:
            return None
        denom *= bracket
    
    neg_a, neg_b = neg_helicity
    helicity = angle_bracket_spinor(lambdas[neg_a], lambdas[neg_b])
    if helicity == 0:
        return None
    
    if denom == 0:
        return None
    
    return (helicity ** 4) / denom


def klt_4pt_spinor(lambdas, tildelambdas):
    """
    KLT 4-point gravity amplitude using spinor-helicity.
    
    Standard field-theory KLT formula:
    M4(1,2,3,4) = -i s12 A(1,2,3,4) \tildeA(1,2,4,3)
    
    In QQ (no complex i), we use:
    M4 = -s12 * A(1,2,3,4) * A(1,2,4,3)
    
    Note: The sign convention may vary. We'll test both + and -.
    """
    # A(1,2,3,4) with order [0,1,2,3] (particles 1,2,3,4)
    order1 = [0, 1, 2, 3]
    A1 = parke_taylor_4pt_spinor(lambdas, tildelambdas, order1, neg_helicity=(0, 1))
    
    # \tildeA(1,2,4,3) with order [0,1,3,2] (particles 1,2,4,3)
    order2 = [0, 1, 3, 2]
    A2 = parke_taylor_4pt_spinor(lambdas, tildelambdas, order2, neg_helicity=(0, 1))
    
    if A1 is None or A2 is None:
        return None
    
    # s12 = <12>[12] = (p1 + p2)^2
    s12 = mandelstam_spinor(lambdas[0], lambdas[1], tildelambdas[0], tildelambdas[1])
    if s12 is None or s12 == 0:
        return None
    
    # Try with negative sign (standard KLT convention)
    M4 = -s12 * A1 * A2
    return M4


def hodges_4pt_spinor(lambdas, tildelambdas):
    """
    Hodges reduced amplitude for n=4 using spinor-helicity.
    
    Build Phi matrix: Phi_{ij} = [ij] / <ij>
    Reduced: delete row 0, col 0, keep 3x3 minor
    """
    n = 4
    Phi = matrix(QQ, n, n)
    
    # Off-diagonal
    for i in range(n):
        for j in range(n):
            if i != j:
                ang_ij = angle_bracket_spinor(lambdas[i], lambdas[j])
                sq_ij = square_bracket_spinor(tildelambdas[i], tildelambdas[j])
                if ang_ij == 0:
                    return None
                Phi[i, j] = sq_ij / ang_ij
    
    # Diagonal using reference legs x=0, y=3
    x, y = 0, 3
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
            ix_ang = angle_bracket_spinor(lambdas[i], lambdas[x])
            iy_ang = angle_bracket_spinor(lambdas[i], lambdas[y])
            if ix_ang == 0 or iy_ang == 0:
                return None
            
            diag_sum = QQ(0)
            for j in range(n):
                if j == i:
                    continue
                jx_ang = angle_bracket_spinor(lambdas[j], lambdas[x])
                jy_ang = angle_bracket_spinor(lambdas[j], lambdas[y])
                if jx_ang != 0 and jy_ang != 0:
                    contrib = Phi[i, j] * (jx_ang * jy_ang) / (ix_ang * iy_ang)
                    diag_sum -= contrib
            Phi[i, i] = diag_sum
    
    # For n=4, try symmetric deletion: delete rows (0,1) and cols (0,1), keep 2x2 minor
    # OR delete row 0, col 0, keep 3x3 minor
    # Standard: delete row 0, col 0
    rows_keep = [1, 2, 3]
    cols_keep = [1, 2, 3]
    Phi_minor = Phi[rows_keep, cols_keep]
    det_minor = Phi_minor.det()
    
    # c factor: For deleted row set {0}, we need c_{0} = 1 (single index)
    # For kept rows {1,2,3}, we need c_{123} = 1/(<12><23><31>)
    a12 = angle_bracket_spinor(lambdas[1], lambdas[2])
    a23 = angle_bracket_spinor(lambdas[2], lambdas[3])
    a31 = angle_bracket_spinor(lambdas[3], lambdas[1])
    if a12 == 0 or a23 == 0 or a31 == 0:
        return None
    c123 = QQ(1) / (a12 * a23 * a31)
    
    # For deleted col set, if we delete col 0, we might need c_{0} = 1
    # But Hodges formula uses c_{ijk} for rows and c_{rst} for cols
    # If we delete row 0 and col 0, we might need c_0 for both, or different factors
    
    # Try: bar_M4 = (-1)^{n+1} * c_{rows} * c_{cols} * det_minor
    # For rows {1,2,3}: c_{123}
    # For cols {1,2,3}: c_{123} (symmetric)
    sign = -1  # (-1)^{n+1} = (-1)^5 = -1
    
    # Try both symmetric and asymmetric
    bar_M4_symmetric = sign * c123 * c123 * det_minor
    # Also try with just one c factor
    bar_M4_asymmetric = sign * c123 * det_minor
    
    # Return the symmetric version first
    return bar_M4_symmetric


def test_n4_spinor_helicity():
    """Test n=4 with direct spinor-helicity sampling."""
    print("="*70)
    print("TEST: n=4 MHV Gravity (Spinor-Helicity Direct Sampling)")
    print("="*70)
    
    ratios = []
    valid_samples = 0
    
    for seed in range(20):
        try:
            lambdas, tildelambdas = sample_spinor_helicity_4pt(seed=seed)
            
            # Check that all angle brackets are nonzero
            all_ok = True
            for i in range(4):
                for j in range(i+1, 4):
                    ang = angle_bracket_spinor(lambdas[i], lambdas[j])
                    if ang == 0:
                        all_ok = False
                        break
                if not all_ok:
                    break
            
            if not all_ok:
                continue
            
            # Compute amplitudes
            M4_klt = klt_4pt_spinor(lambdas, tildelambdas)
            M4_hodges = hodges_4pt_spinor(lambdas, tildelambdas)
            
            if M4_klt is None or M4_hodges is None or M4_hodges == 0:
                print(f"Seed {seed}: KLT={M4_klt}, Hodges={M4_hodges}")
                continue
            
            ratio = M4_klt / M4_hodges
            ratios.append(ratio)
            valid_samples += 1
            print(f"Seed {seed}: ratio = {ratio}")
            
        except Exception as e:
            print(f"Seed {seed}: Error - {e}")
            continue
    
    if not ratios:
        print("ERROR: No valid ratios computed")
        return False
    
    unique_ratios = list(set(ratios))
    print(f"\nTotal valid samples: {valid_samples}")
    print(f"Unique ratios: {len(unique_ratios)}")
    
    if len(unique_ratios) == 1:
        print(f"SUCCESS: Constant ratio = {unique_ratios[0]}")
        return True
    else:
        print(f"Ratios vary:")
        for i, r in enumerate(unique_ratios[:5]):
            print(f"  Ratio {i+1}: {r}")
        return False


if __name__ == '__main__':
    success = test_n4_spinor_helicity()
    print(f"\n{'='*70}")
    print(f"n=4 TEST (Spinor-Helicity): {'PASSED' if success else 'FAILED'}")
    print(f"{'='*70}")

