#!/usr/bin/env sage
# =============================================================================
# DIAGNOSTIC PROBE 2: Spinor Helicity Invariance
# =============================================================================
# Tests Hodges invariance using TRUE momentum-conserving spinor variables.
# Replaces Z-moment-curve sampling with direct lambda/tilde_lambda sampling.
# =============================================================================

from sage.all import *
import sys
import os
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/spinor_sampling.sage')
load('src/hodges_sigma.sage')

def sample_reference_spinors_random(lambdas, seed_offset=0):
    """Sample random reference spinors λ_x, λ_y."""
    import random
    random.seed(int(42 + seed_offset))
    
    max_attempts = 1000
    for attempt in range(max_attempts):
        lambda_x = vector(QQ, [QQ(random.randint(-10, 10)), QQ(random.randint(-10, 10))])
        lambda_y = vector(QQ, [QQ(random.randint(-10, 10)), QQ(random.randint(-10, 10))])
        
        if lambda_x == 0 or lambda_y == 0: continue
        
        all_ok = True
        for i in range(len(lambdas)):
            lambda_i = lambdas[i]
            # <i x>
            ang_ix = lambda_i[0] * lambda_x[1] - lambda_i[1] * lambda_x[0]
            # <i y>
            ang_iy = lambda_i[0] * lambda_y[1] - lambda_i[1] * lambda_y[0]
            if ang_ix == 0 or ang_iy == 0:
                all_ok = False
                break
        
        if all_ok:
            return (lambda_x, lambda_y)
            
    print("Warning: Failed to sample reference spinors after 1000 attempts.")
    return (None, None)

def ang_vec(a, b):
    return a[0] * b[1] - a[1] * b[0]

def compute_hodges_spinor(lambdas, tilde_lambdas, rows_delete, cols_delete, lambda_x, lambda_y):
    """
    Compute Hodges reduced amplitude using explicit spinor helicity variables.
    """
    n = len(lambdas)
    Phi = matrix(QQ, n, n)
    
    # Precompute angles and squares
    def get_angle(i, j):
        return ang_vec(lambdas[i], lambdas[j])
        
    def get_square(i, j):
        return ang_vec(tilde_lambdas[i], tilde_lambdas[j])
        
    # Off-diagonal
    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = get_angle(i, j)
                if ij_ang == 0: 
                    print(f"Angle <{i}{j}> is zero")
                    return None
                ij_sq = get_square(i, j)
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
    ang_ij = get_angle(i, j)
    ang_jk = get_angle(j, k)
    ang_ki = get_angle(k, i)
    if ang_ij == 0 or ang_jk == 0 or ang_ki == 0:
        print("c_rows denominator zero")
        return None
    c_rows = QQ(1) / (ang_ij * ang_jk * ang_ki)
    
    r, s, t = cols_delete
    ang_rs = get_angle(r, s)
    ang_st = get_angle(s, t)
    ang_tr = get_angle(t, r)
    if ang_rs == 0 or ang_st == 0 or ang_tr == 0:
        print("c_cols denominator zero")
        return None
    c_cols = QQ(1) / (ang_rs * ang_st * ang_tr)
    
    sigma = hodges_sigma(rows_delete, cols_delete, n)
    sign = -1 # (-1)^(n+1)
    
    return sign * sigma * c_rows * c_cols * det_minor

def probe_1_spinor_reference(lambdas, tilde_lambdas):
    """PROBE 1: Reference Invariance with Spinors."""
    print("\n[PROBE 1] Reference Invariance (Spinor Helicity)")
    baseline_rows = [0, 1, 2]
    baseline_cols = [3, 4, 5]
    
    values = []
    for k in range(10):
        lx, ly = sample_reference_spinors_random(lambdas, seed_offset=k*10)
        if lx is None: continue
        val = compute_hodges_spinor(lambdas, tilde_lambdas, baseline_rows, baseline_cols, lx, ly)
        if val is not None:
            values.append(val)
        else:
            print(f"Computation failed for seed offset {k*10}")
            
    unique_vals = list(set(values))
    print(f"Tested {len(values)} reference pairs. Unique values: {len(unique_vals)}")
    if len(unique_vals) == 1:
        print("RESULT: PASS")
        print(f"Value: {unique_vals[0]}")
        return True, unique_vals[0]
    else:
        print(f"RESULT: FAIL. Values: {unique_vals}")
        return False, None

def probe_2_spinor_deletion(lambdas, tilde_lambdas, valid_ref_pair):
    """PROBE 2: Deletion Invariance with Spinors (Disjoint Complements)."""
    print("\n[PROBE 2] Deletion Invariance (Spinor Helicity)")
    lx, ly = valid_ref_pair
    
    all_rows = list(combinations(range(6), 3))
    deletions = []
    full_set = set(range(6))
    for r_del in all_rows:
        c_del = tuple(sorted(list(full_set - set(r_del))))
        deletions.append((list(r_del), list(c_del)))
        
    print(f"Testing {len(deletions)} disjoint complement pairs.")
    
    results = []
    for r_del, c_del in deletions:
        val = compute_hodges_spinor(lambdas, tilde_lambdas, r_del, c_del, lx, ly)
        if val is not None:
            results.append(val)
        else:
            print(f"Failed for deletion {r_del}, {c_del}")
            
    unique_vals = list(set(results))
    print(f"Computed {len(results)} deletions. Unique values: {len(unique_vals)}")
    
    if len(unique_vals) == 1:
        print("RESULT: PASS (All identical)")
        print(f"Value: {unique_vals[0]}")
        return True
    else:
        print(f"RESULT: FAIL. Values: {unique_vals[:5]}")
        return False

def run_spinor_diagnostics():
    print("="*70)
    print("DIAGNOSTIC PROBE SUITE: SPINOR HELICITY")
    print("="*70)
    
    # 1. Sample Valid Kinematics
    res = sample_spinor_helicity_conserving(n=6, seed=0)
    if res is None:
        print("Fatal: Sampling failed.")
        return
    lambdas, tilde_lambdas = res
    print("Sampled valid 4D kinematics (momentum conserved).")
    
    # 2. Probe 1
    p1_pass, baseline_val = probe_1_spinor_reference(lambdas, tilde_lambdas)
    if not p1_pass: return
    
    # 3. Probe 2
    lx, ly = sample_reference_spinors_random(lambdas, seed_offset=999)
    if lx is None:
        print("Failed to sample ref spinors for Probe 2")
        return
    p2_pass = probe_2_spinor_deletion(lambdas, tilde_lambdas, (lx, ly))
    
    if p2_pass:
        print("\nCONCLUSION: Hodges formula works perfectly with valid 4D kinematics!")
    else:
        print("\nCONCLUSION: Still failing even with valid kinematics. Deep math/formula bug.")

if __name__ == '__main__':
    run_spinor_diagnostics()
