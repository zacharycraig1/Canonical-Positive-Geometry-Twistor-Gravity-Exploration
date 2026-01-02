#!/usr/bin/env sage
# =============================================================================
# SYSTEMATIC n=4 CONVENTION TESTING
# =============================================================================
# Test all combinations of:
# - KLT sign conventions
# - KLT orderings
# - Hodges deletion patterns
# - Hodges normalization factors
# =============================================================================

from sage.all import *
load('tests/test_n4_spinor_helicity.sage')

def test_convention_combo(lambdas, tildelambdas, klt_sign, klt_order2, hodges_deletion, hodges_c_symmetric):
    """Test one convention combination."""
    # KLT
    order1 = [0, 1, 2, 3]
    A1 = parke_taylor_4pt_spinor(lambdas, tildelambdas, order1, neg_helicity=(0, 1))
    A2 = parke_taylor_4pt_spinor(lambdas, tildelambdas, klt_order2, neg_helicity=(0, 1))
    if A1 is None or A2 is None:
        return None
    s12 = mandelstam_spinor(lambdas[0], lambdas[1], tildelambdas[0], tildelambdas[1])
    if s12 is None or s12 == 0:
        return None
    M4_klt = klt_sign * s12 * A1 * A2
    
    # Hodges
    n = 4
    Phi = matrix(QQ, n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                ang_ij = angle_bracket_spinor(lambdas[i], lambdas[j])
                sq_ij = square_bracket_spinor(tildelambdas[i], tildelambdas[j])
                if ang_ij == 0:
                    return None
                Phi[i, j] = sq_ij / ang_ij
    
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
    
    # Deletion pattern
    if hodges_deletion == "row0_col0":
        rows_keep = [1, 2, 3]
        cols_keep = [1, 2, 3]
    elif hodges_deletion == "row0_col1":
        rows_keep = [1, 2, 3]
        cols_keep = [0, 2, 3]
    elif hodges_deletion == "row1_col1":
        rows_keep = [0, 2, 3]
        cols_keep = [0, 2, 3]
    else:
        rows_keep = [1, 2, 3]
        cols_keep = [1, 2, 3]
    
    Phi_minor = Phi[rows_keep, cols_keep]
    det_minor = Phi_minor.det()
    
    # c factors
    a12 = angle_bracket_spinor(lambdas[1], lambdas[2])
    a23 = angle_bracket_spinor(lambdas[2], lambdas[3])
    a31 = angle_bracket_spinor(lambdas[3], lambdas[1])
    if a12 == 0 or a23 == 0 or a31 == 0:
        return None
    c123 = QQ(1) / (a12 * a23 * a31)
    
    if hodges_c_symmetric:
        c_factor = c123 * c123
    else:
        c_factor = c123
    
    sign = -1
    bar_M4 = sign * c_factor * det_minor
    
    if bar_M4 == 0:
        return None
    
    ratio = M4_klt / bar_M4
    return ratio

def test_all_conventions():
    """Test all convention combinations."""
    print("="*70)
    print("SYSTEMATIC n=4 CONVENTION TESTING")
    print("="*70)
    
    # Test samples
    test_seeds = [1, 2, 3, 4, 5]
    
    # Convention options
    klt_signs = [-1, 1]
    klt_order2s = [[0, 1, 3, 2], [0, 3, 1, 2], [0, 2, 1, 3]]
    hodges_deletions = ["row0_col0", "row0_col1", "row1_col1"]
    hodges_c_symmetrics = [True, False]
    
    best_combo = None
    best_constant_count = 0
    
    for klt_sign in klt_signs:
        for klt_order2 in klt_order2s:
            for hodges_deletion in hodges_deletions:
                for hodges_c_symmetric in hodges_c_symmetrics:
                    ratios = []
                    combo_name = f"KLT_sign={klt_sign}, order2={klt_order2}, del={hodges_deletion}, c_sym={hodges_c_symmetric}"
                    
                    for seed in test_seeds:
                        try:
                            lambdas, tildelambdas = sample_spinor_helicity_4pt(seed=seed)
                            
                            # Check angle brackets
                            all_ok = True
                            for i in range(4):
                                for j in range(i+1, 4):
                                    if angle_bracket_spinor(lambdas[i], lambdas[j]) == 0:
                                        all_ok = False
                                        break
                                if not all_ok:
                                    break
                            if not all_ok:
                                continue
                            
                            ratio = test_convention_combo(lambdas, tildelambdas, klt_sign, klt_order2, hodges_deletion, hodges_c_symmetric)
                            if ratio is not None:
                                ratios.append(ratio)
                        except:
                            continue
                    
                    if ratios:
                        unique_ratios = list(set(ratios))
                        if len(unique_ratios) == 1:
                            print(f"\n*** SUCCESS: {combo_name}")
                            print(f"  Constant ratio = {unique_ratios[0]}")
                            return combo_name, unique_ratios[0]
                        elif len(unique_ratios) < len(ratios):
                            print(f"\nPartial: {combo_name}")
                            print(f"  {len(unique_ratios)} unique ratios from {len(ratios)} samples")
                            if len(unique_ratios) > best_constant_count:
                                best_constant_count = len(unique_ratios)
                                best_combo = combo_name
    
    print(f"\nNo constant ratio found. Best: {best_combo} with {best_constant_count} unique ratios")
    return None, None

if __name__ == '__main__':
    combo, ratio = test_all_conventions()
    if combo:
        print(f"\n{'='*70}")
        print(f"FOUND CONSTANT RATIO: {combo}")
        print(f"Ratio = {ratio}")
        print(f"{'='*70}")










