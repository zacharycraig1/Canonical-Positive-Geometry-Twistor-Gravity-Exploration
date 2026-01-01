#!/usr/bin/env sage
# =============================================================================
# IDENTIFY CONVENTION CLASSES FOR 6 RATIO VALUES
# =============================================================================
# Goal: Determine what discrete choice causes each of the 6 ratio values
# =============================================================================

from sage.all import *
import sys
import os
from collections import defaultdict

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
    # c for rows_delete
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
    
    # c for cols_delete
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

def identify_convention_classes(n_samples=50):
    """Identify what convention choice causes each ratio value."""
    print("="*70)
    print("IDENTIFYING CONVENTION CLASSES FOR 6 RATIO VALUES")
    print("="*70)
    
    # Test with standard implementation first
    print("\n1. Testing with STANDARD Hodges (rows_delete=[0,1,2], cols_delete=[3,4,5])...")
    standard_ratios = defaultdict(list)
    
    for seed in range(n_samples):
        Z = sample_positive_Z_moment_curve(n=6, seed=seed)
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
        if not twistor.domain_ok:
            continue
        
        H_red = hodges_6pt_mhv_reduced(twistor)
        A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
        
        H_val = H_red[0] if isinstance(H_red, tuple) else H_red
        A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
        
        if H_val is None or A_val is None or H_val == 0:
            continue
        
        ratio = A_val / H_val
        ratio_str = str(ratio)
        standard_ratios[ratio_str].append(seed)
    
    print(f"   Found {len(standard_ratios)} unique ratios with standard implementation")
    for i, (ratio_str, seeds) in enumerate(standard_ratios.items()):
        print(f"   Ratio {i+1}: {len(seeds)} samples, seeds: {seeds[:5]}...")
    
    # Test alternative deletion patterns
    print("\n2. Testing alternative deletion patterns...")
    
    # All ways to choose 3 rows and 3 cols from 6
    from itertools import combinations
    all_rows = list(combinations(range(6), 3))
    all_cols = list(combinations(range(6), 3))
    
    # Test a few key patterns
    test_patterns = [
        ([0,1,2], [3,4,5]),  # Standard
        ([0,1,2], [0,1,2]),  # Symmetric
        ([3,4,5], [3,4,5]),  # Other symmetric
        ([0,1,3], [2,4,5]),  # Mixed
    ]
    
    for rows_del, cols_del in test_patterns:
        pattern_name = f"rows_del={rows_del}, cols_del={cols_del}"
        print(f"\n   Testing {pattern_name}...")
        
        pattern_ratios = defaultdict(list)
        for seed in range(min(20, n_samples)):
            Z = sample_positive_Z_moment_curve(n=6, seed=seed)
            twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
            if not twistor.domain_ok:
                continue
            
            H_alt = test_hodges_deletion_pattern(twistor, list(rows_del), list(cols_del))
            A_klt = gravity_6pt_mhv_klt(twistor, mandelstam_invariant)
            
            H_val = H_alt
            A_val = A_klt[0] if isinstance(A_klt, tuple) else A_klt
            
            if H_val is None or A_val is None or H_val == 0:
                continue
            
            ratio = A_val / H_val
            ratio_str = str(ratio)
            pattern_ratios[ratio_str].append(seed)
        
        print(f"      Found {len(pattern_ratios)} unique ratios")
        if len(pattern_ratios) == 1:
            print(f"      *** CONSTANT RATIO: {list(pattern_ratios.keys())[0][:60]}...")
            return pattern_name, list(pattern_ratios.keys())[0]
    
    # Check if ratio correlates with seed mod 6
    print("\n3. Checking correlation with seed mod 6...")
    for ratio_str, seeds in standard_ratios.items():
        mod6_counts = defaultdict(int)
        for s in seeds:
            mod6_counts[s % 6] += 1
        print(f"   Ratio {ratio_str[:40]}...: mod6 distribution = {dict(mod6_counts)}")
    
    return None, None

if __name__ == '__main__':
    pattern, ratio = identify_convention_classes(n_samples=50)
    if pattern:
        print(f"\n{'='*70}")
        print(f"FOUND CONSTANT RATIO WITH PATTERN: {pattern}")
        print(f"Ratio = {ratio}")
        print(f"{'='*70}")
    else:
        print(f"\n{'='*70}")
        print("No single pattern found that gives constant ratio")
        print("Need to investigate further (KLT fixed legs, reference legs, etc.)")
        print(f"{'='*70}")









