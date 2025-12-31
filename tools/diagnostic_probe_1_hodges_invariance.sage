#!/usr/bin/env sage
# =============================================================================
# DIAGNOSTIC PROBE: Hodges Invariance Tests (Reference & Deletion)
# =============================================================================
# Implements:
# Probe 1: Reference independence (fixed disjoint deletion, varying reference spinors)
# Probe 2: Deletion independence (fixed reference spinors, all 20 disjoint complement pairs)
# =============================================================================

from sage.all import *
import sys
import os
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/sampling.sage')
load('src/hodges.sage')
load('src/hodges_sigma.sage')

def sample_reference_spinors_random(twistor, seed_offset=0):
    """Sample random reference spinors λ_x, λ_y."""
    import random
    random.seed(int(42 + seed_offset))
    
    max_attempts = 100
    for attempt in range(max_attempts):
        # Random rational components
        lambda_x = vector(QQ, [QQ(random.randint(-10, 10)), QQ(random.randint(-10, 10))])
        lambda_y = vector(QQ, [QQ(random.randint(-10, 10)), QQ(random.randint(-10, 10))])
        
        # Ensure non-zero and independent
        if lambda_x == 0 or lambda_y == 0: continue
        
        all_ok = True
        for i in range(twistor.n):
            lambda_i = vector(QQ, [twistor.Z[i][0], twistor.Z[i][1]])
            ang_ix = lambda_i[0] * lambda_x[1] - lambda_i[1] * lambda_x[0]
            ang_iy = lambda_i[0] * lambda_y[1] - lambda_i[1] * lambda_y[0]
            if ang_ix == 0 or ang_iy == 0:
                all_ok = False
                break
        
        if all_ok:
            return (lambda_x, lambda_y)
    return (None, None)

def ang_vec(a, b):
    """Angle bracket for 2-vectors."""
    return a[0] * b[1] - a[1] * b[0]

def compute_hodges_explicit(twistor, rows_delete, cols_delete, lambda_x, lambda_y):
    """
    Compute Hodges reduced amplitude explicitly with given refs and deletion.
    """
    n = 6
    Phi = matrix(QQ, n, n)
    
    lambdas = []
    for i in range(n):
        lambdas.append(vector(QQ, [twistor.Z[i][0], twistor.Z[i][1]]))
    
    # Off-diagonal
    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = twistor.get_angle(i, j)
                if ij_ang == 0: return None
                ij_sq = twistor.get_square(i, j)
                if ij_sq is None: return None
                Phi[i, j] = ij_sq / ij_ang
    
    # Diagonal
    for i in range(n):
        lambda_i = lambdas[i]
        ang_ix = ang_vec(lambda_i, lambda_x)
        ang_iy = ang_vec(lambda_i, lambda_y)
        
        diag_sum = QQ(0)
        for j in range(n):
            if j == i: continue
            lambda_j = lambdas[j]
            ang_jx = ang_vec(lambda_j, lambda_x)
            ang_jy = ang_vec(lambda_j, lambda_y)
            
            contrib = Phi[i, j] * (ang_jx * ang_jy) / (ang_ix * ang_iy)
            diag_sum -= contrib
        Phi[i, i] = diag_sum
        
    rows_keep = [i for i in range(n) if i not in rows_delete]
    cols_keep = [i for i in range(n) if i not in cols_delete]
    
    Phi_minor = Phi[rows_keep, cols_keep]
    det_minor = Phi_minor.det()
    
    # c factors
    i, j, k = rows_delete
    c_rows = QQ(1) / (twistor.get_angle(i, j) * twistor.get_angle(j, k) * twistor.get_angle(k, i))
    
    r, s, t = cols_delete
    c_cols = QQ(1) / (twistor.get_angle(r, s) * twistor.get_angle(s, t) * twistor.get_angle(t, r))
    
    sigma = hodges_sigma(rows_delete, cols_delete, n)
    sign = -1 # (-1)^(n+1)
    
    return sign * sigma * c_rows * c_cols * det_minor

def probe_1_reference_invariance(twistor):
    """PROBE 1: Test reference independence for one fixed deletion."""
    print("\n[PROBE 1] Reference Invariance Test")
    print("Settings: Deletion rows=[0,1,2], cols=[3,4,5] (Fixed)")
    
    baseline_rows = [0, 1, 2]
    baseline_cols = [3, 4, 5]
    
    values = []
    
    for k in range(20):
        lx, ly = sample_reference_spinors_random(twistor, seed_offset=k*10)
        if lx is None: continue
        
        val = compute_hodges_explicit(twistor, baseline_rows, baseline_cols, lx, ly)
        if val is not None:
            values.append(val)
            
    unique_vals = list(set(values))
    print(f"Tested {len(values)} reference pairs.")
    print(f"Unique values: {len(unique_vals)}")
    if len(unique_vals) == 1:
        print("RESULT: PASS (All values identical)")
        return True, unique_vals[0]
    else:
        print("RESULT: FAIL")
        print(f"Values found: {unique_vals[:5]}")
        return False, None

def probe_2_deletion_invariance(twistor, valid_ref_pair):
    """PROBE 2: Test deletion independence across all 20 disjoint complement pairs."""
    print("\n[PROBE 2] Deletion Invariance Test")
    
    lx, ly = valid_ref_pair
    print(f"Using fixed reference spinors from Probe 1.")
    
    # Generate all 20 combinations of 3 rows from 6
    all_rows = list(combinations(range(6), 3))
    # For each, J is the complement
    deletions = []
    full_set = set(range(6))
    for r_del in all_rows:
        c_del = tuple(sorted(list(full_set - set(r_del))))
        deletions.append((list(r_del), list(c_del)))
        
    print(f"Testing {len(deletions)} disjoint complement pairs (I, J=I^c).")
    
    results = []
    for r_del, c_del in deletions:
        val = compute_hodges_explicit(twistor, r_del, c_del, lx, ly)
        if val is not None:
            results.append(val)
            
    unique_vals = list(set(results))
    print(f"Computed {len(results)}/{len(deletions)} deletions.")
    print(f"Unique values: {len(unique_vals)}")
    
    if len(unique_vals) == 1:
        print("RESULT: PASS (All deletions give identical result)")
        return True
    else:
        print("RESULT: FAIL")
        print(f"Unique values found (first 5): {unique_vals[:5]}")
        
        # Analyze ratios relative to first
        baseline = results[0]
        ratios = [val/baseline for val in unique_vals]
        print(f"Ratios relative to first: {ratios[:10]}")
        return False

def run_diagnostics():
    print("="*70)
    print("DIAGNOSTIC PROBE SUITE: Hodges Invariance")
    print("="*70)
    
    # 1. Kinematics
    Z = sample_positive_Z_moment_curve(n=6, seed=0)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)
    if not twistor.domain_ok:
        print("Fatal: Sample not in domain.")
        return

    # 2. Run Probe 1
    p1_pass, baseline_val = probe_1_reference_invariance(twistor)
    
    if not p1_pass:
        print("\nSTOPPING: Reference invariance failed. Fix diagonal/kinematics first.")
        return

    # 3. Run Probe 2 (using a valid ref pair from Probe 1)
    # Just sample one fresh pair that works
    lx, ly = sample_reference_spinors_random(twistor, seed_offset=999)
    if lx is None:
        print("Error sampling ref pair for Probe 2")
        return
        
    p2_pass = probe_2_deletion_invariance(twistor, (lx, ly))
    
    if not p2_pass:
        print("\nDIAGNOSIS: Reference invariance holds, but Deletion invariance fails.")
        print("LIKELY CAUSE: Kinematics (Z) are not valid momentum twistors corresponding")
        print("to a momentum-conserving 4D configuration, or [ij] definition is inconsistent.")
        print("ACTION: Implement spinor-helicity sampling with momentum conservation.")

if __name__ == '__main__':
    run_diagnostics()
