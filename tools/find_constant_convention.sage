#!/usr/bin/env sage
# =============================================================================
# FIND CONSTANT CONVENTION BY TESTING ALL DELETION PATTERNS
# =============================================================================

from sage.all import *
import sys
import os
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/klt.sage')

def test_hodges_deletion_pattern(twistor, rows_delete, cols_delete):
    """Test Hodges with specific deletion pattern."""
    n = 6
    Phi = matrix(QQ, n, n)
    x, y = 0, 5
    
    # Build Phi matrix
    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = twistor.get_angle(i, j)
                if ij_ang == 0:
                    return None
                ij_sq = twistor.get_square(i, j)
                if ij_sq is None:
                    return None
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
                return None
            diag_sum = QQ(0)
            for j in range(n):
                if j == i:
                    continue
                jx_ang = twistor.get_angle(j, x)
                jy_ang = twistor.get_angle(j, y)
                if jx_ang != 0 and jy_ang != 0:
                    contrib = Phi[i, j] * (jx_ang * jy_ang) / (ix_ang * iy_ang)
                    diag_sum -= contrib
            Phi[i, i] = diag_sum
    
    # Deletion pattern
    rows_keep = [i for i in range(n) if i not in rows_delete]
    cols_keep = [i for i in range(n) if i not in cols_delete]
    
    if len(rows_keep) != 3 or len(cols_keep) != 3:
        return None
    
    Phi_minor = Phi[rows_keep, cols_keep]
    try:
        det_minor = Phi_minor.det()
    except:
        return None
    
    # c factors
    if len(rows_delete) == 3:
        i, j, k = rows_delete
        aij = twistor.get_angle(i, j)
        ajk = twistor.get_angle(j, k)
        aki = twistor.get_angle(k, i)
        if aij == 0 or ajk == 0 or aki == 0:
            return None
        c_rows = QQ(1) / (aij * ajk * aki)
    else:
        return None
    
    if len(cols_delete) == 3:
        r, s, t = cols_delete
        ars = twistor.get_angle(r, s)
        ast = twistor.get_angle(s, t)
        atr = twistor.get_angle(t, r)
        if ars == 0 or ast == 0 or atr == 0:
            return None
        c_cols = QQ(1) / (ars * ast * atr)
    else:
        return None
    
    sign = -1
    bar_M6 = sign * c_rows * c_cols * det_minor
    return bar_M6

def find_best_convention(n_samples=200):
    """Test all deletion patterns to find one that gives constant ratio."""
    print("="*70)
    print("FINDING CONSTANT CONVENTION (Testing All Deletion Patterns)")
    print("="*70)
    
    # Generate test samples
    test_samples = []
    for seed in range(n_samples):
        Z = sample_positive_Z_moment_curve(n=6, seed=seed)
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
        if twistor.domain_ok:
            test_samples.append((seed, twistor))
    
    print(f"Generated {len(test_samples)} valid test samples")
    
    # Test standard first
    print("\nTesting STANDARD pattern: rows_del=[0,1,2], cols_del=[3,4,5]...")
    standard_ratios = []
    for seed, twistor in test_samples:
        H_red = hodges_6pt_mhv_reduced(twistor)
        A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
        H_val = H_red[0] if isinstance(H_red, tuple) else H_red
        A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
        if H_val is not None and A_val is not None and H_val != 0:
            standard_ratios.append(A_val / H_val)
    
    unique_standard = list(set(standard_ratios))
    print(f"  Standard pattern: {len(unique_standard)} unique ratios from {len(standard_ratios)} samples")
    
    # Test key symmetric patterns
    key_patterns = [
        ([0,1,2], [0,1,2], "symmetric_012"),
        ([3,4,5], [3,4,5], "symmetric_345"),
        ([0,1,2], [3,4,5], "standard"),
        ([0,1,3], [2,4,5], "mixed_1"),
        ([0,2,4], [1,3,5], "mixed_2"),
    ]
    
    best_pattern = None
    best_unique = len(unique_standard)
    best_ratio = None
    
    for rows_del, cols_del, name in key_patterns:
        print(f"\nTesting {name}: rows_del={rows_del}, cols_del={cols_del}...")
        pattern_ratios = []
        
        for seed, twistor in test_samples:
            H_alt = test_hodges_deletion_pattern(twistor, rows_del, cols_del)
            A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
            H_val = H_alt
            A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
            if H_val is not None and A_val is not None and H_val != 0:
                pattern_ratios.append(A_val / H_val)
        
        if pattern_ratios:
            unique_pattern = list(set(pattern_ratios))
            print(f"  {name}: {len(unique_pattern)} unique ratios from {len(pattern_ratios)} samples")
            
            if len(unique_pattern) < best_unique:
                best_unique = len(unique_pattern)
                best_pattern = (rows_del, cols_del, name)
                best_ratio = unique_pattern[0] if len(unique_pattern) == 1 else None
                
                if len(unique_pattern) == 1:
                    print(f"  *** SUCCESS: Constant ratio = {unique_pattern[0]}")
                    return best_pattern, unique_pattern[0]
    
    print(f"\nBest pattern: {best_pattern[2]} with {best_unique} unique ratios")
    return best_pattern, best_ratio

if __name__ == '__main__':
    pattern, ratio = find_best_convention(n_samples=200)
    print(f"\n{'='*70}")
    if ratio:
        print(f"FOUND CONSTANT RATIO!")
        print(f"Pattern: {pattern[2]}")
        print(f"Rows delete: {pattern[0]}, Cols delete: {pattern[1]}")
        print(f"Constant ratio = {ratio}")
    else:
        print(f"Best pattern: {pattern[2] if pattern else 'None'}")
        print(f"Unique ratios: {len(set([ratio])) if ratio else 'N/A'}")
    print(f"{'='*70}")









