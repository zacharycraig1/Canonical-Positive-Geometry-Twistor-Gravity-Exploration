
from sage.all import *
import sys
import os
import numpy as np
from itertools import permutations, combinations

sys.path.append(os.getcwd())

# We define MomentumTwistor locally
load("src/hodges.sage")
load("klt_search/spinor_tools.sage")

def ts():
    import time
    return time.strftime("%H:%M:%S")

def print_log(msg):
    print(f"[{ts()}] {msg}")

# ROBUST HODGES
def robust_hodges_6pt_mhv(twistor, deleted_rows=None, ref_spinors=None):
    n = twistor.n
    if deleted_rows is None: deleted_rows = [0, 4, 5]
    deleted_rows = sorted(deleted_rows)
    rows_to_keep = [i for i in range(n) if i not in deleted_rows]
    cols_to_keep = rows_to_keep
    
    Phi = matrix(QQ, n, n)
    
    if ref_spinors is None:
        x_idx, y_idx = deleted_rows[0], deleted_rows[1]
        lambda_x = vector(QQ, [twistor.Z[x_idx][0], twistor.Z[x_idx][1]])
        lambda_y = vector(QQ, [twistor.Z[y_idx][0], twistor.Z[y_idx][1]])
    else:
        lambda_x, lambda_y = ref_spinors
        
    lambdas = []
    for i in range(n):
        lambdas.append(vector(QQ, [twistor.Z[i][0], twistor.Z[i][1]]))
        
    def ang_prod(l1, l2): return l1[0]*l2[1] - l1[1]*l2[0]

    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = twistor.get_angle(i, j)
                if ij_ang == 0: return None, "domain"
                ij_sq = twistor.get_square(i, j)
                if ij_sq is None: return None, "domain"
                Phi[i, j] = ij_sq / ij_ang
                
    for i in range(n):
        ix_ang = ang_prod(lambdas[i], lambda_x)
        iy_ang = ang_prod(lambdas[i], lambda_y)
        if ix_ang == 0 or iy_ang == 0:
            if i in rows_to_keep: return None, "domain_diag"
            Phi[i, i] = 0
            continue
        diag_sum = QQ(0)
        for j in range(n):
            if j == i: continue
            jx_ang = ang_prod(lambdas[j], lambda_x)
            jy_ang = ang_prod(lambdas[j], lambda_y)
            if jx_ang == 0 or jy_ang == 0: continue
            contrib = Phi[i, j] * (jx_ang * jy_ang) / (ix_ang * iy_ang)
            diag_sum -= contrib
        Phi[i, i] = diag_sum
        
    Phi_red = Phi[rows_to_keep, cols_to_keep]
    try: det_Phi_red = Phi_red.det()
    except: return None, "det_fail"
    
    r1, r2, r3 = deleted_rows
    norm_factor = QQ(1)
    norm_factor *= twistor.get_angle(r1, r2)
    norm_factor *= twistor.get_angle(r2, r3)
    norm_factor *= twistor.get_angle(r3, r1)
    if norm_factor == 0: return None, "norm_fail"
    norm_factor = norm_factor ** 2
    
    det_prime_Phi = det_Phi_red / norm_factor
    neg_factor = twistor.get_angle(0, 1) ** 8
    
    return det_prime_Phi * neg_factor, "ok"

# KLT IMPLEMENTATION
def mandelstam_s(twistor, i, j):
    ang = twistor.get_angle(i, j)
    sq = twistor.get_square(i, j)
    if ang == 0 or sq is None: return QQ(0)
    return ang * sq

def klt_kernel_S3(alpha, beta, twistor):
    pos_beta = {val: idx for idx, val in enumerate(beta)}
    def theta(a, b): return 1 if pos_beta[a] > pos_beta[b] else 0
    def s(i, j): return mandelstam_s(twistor, i, j)

    i1, i2, i3 = alpha
    term1 = s(0, i1) 
    if theta(i1, i2): term1 += s(i1, i2)
    if theta(i1, i3): term1 += s(i1, i3)
    term2 = s(0, i2)
    if theta(i2, i3): term2 += s(i2, i3)
    term3 = s(0, i3)
    return term1 * term2 * term3

def get_pt(twistor, order):
    denom = QQ(1)
    n = len(order)
    for k in range(n):
        i = order[k]
        j = order[(k+1)%n]
        ang = twistor.get_angle(i, j)
        if ang == 0: return None
        denom *= ang
    num = twistor.get_angle(0, 1)**4
    return num / denom

def calculate_klt_6pt(twistor):
    perm_indices = [1, 2, 3]
    perms = list(permutations(perm_indices))
    M_KLT = QQ(0)
    for alpha in perms:
        order_L = [0] + list(alpha) + [4, 5]
        AL = get_pt(twistor, order_L)
        if AL is None: continue
        for beta in perms:
            order_R = [0] + list(beta) + [5, 4]
            AR = get_pt(twistor, order_R)
            if AR is None: continue
            try: S = klt_kernel_S3(list(alpha), list(beta), twistor)
            except: continue
            M_KLT += AL * S * AR
    return M_KLT

# DIAGNOSTICS
def check_near_singularity(twistor):
    H, msg = robust_hodges_6pt_mhv(twistor)
    if H is None: return False, msg
    return True, {}

def check_permutation_symmetry(twistor):
    del_set = [0, 4, 5]
    H0 = robust_hodges_6pt_mhv(twistor, deleted_rows=del_set)
    if isinstance(H0, tuple): H0 = H0[0]
    if H0 is None: return None
    
    # Proper symmetry check via spinors
    lambdas, t_lambdas = spinors_from_twistors(twistor)
    if None in t_lambdas: return None
    
    # Permute 3 <-> 4 (indices 2, 3)
    perm_map = [0, 1, 3, 2, 4, 5]
    lambdas_new = [lambdas[p] for p in perm_map]
    t_lambdas_new = [t_lambdas[p] for p in perm_map]
    
    Z_new = twistors_from_spinors(lambdas_new, t_lambdas_new)
    if Z_new is None: return None
    
    twistor_perm = MomentumTwistor(n=6, Z=Z_new, check_domain=False)
    
    H_perm = robust_hodges_6pt_mhv(twistor_perm, deleted_rows=del_set)
    if isinstance(H_perm, tuple): H_perm = H_perm[0]
    
    if H_perm is None: return None
    
    ratio = H_perm / H0
    diff = abs(ratio - 1.0)
    
    if diff > 1e-4:
        print_log(f"DEBUG: Symmetry Fail. H0={float(H0):.4e}, H_perm={float(H_perm):.4e}, Ratio={float(ratio):.4f}")
        
    return diff

def check_soft_limit(seed):
    import random
    random.seed(seed)
    Z = []
    for i in range(5):
         Z.append(vector(QQ, [
             random.randint(-10,10), random.randint(-10,10),
             random.randint(-10,10), random.randint(-10,10)
         ]))
    epsilon = QQ(1)/QQ(1000)
    noise = vector(QQ, [1, 2, 3, 4])
    Z5 = Z[4] + Z[0] + epsilon * noise
    Z.append(Z5)
    twistor = MomentumTwistor(n=6, Z=Z, check_domain=False)
    del_set = [0, 4, 5]
    H6 = robust_hodges_6pt_mhv(twistor, deleted_rows=del_set)
    K6 = calculate_klt_6pt(twistor)
    if isinstance(H6, tuple): H6 = H6[0]
    if H6 is None or H6 == 0: return None
    if K6 == 0: return None
    ratio = K6 / H6
    return float(ratio)

def run_diagnostics():
    print_log("Starting Phase B Diagnostics (Robust Hodges + Proper Symmetry)...")
    
    # B1
    print_log("--- B1: Precision / Singularity Filter ---")
    ratios = []
    del_set = [0, 4, 5]
    
    for i in range(20):
        twistor = MomentumTwistor(n=6, seed=i+1000)
        H = robust_hodges_6pt_mhv(twistor, deleted_rows=del_set)
        if isinstance(H, tuple): H = H[0]
        K = calculate_klt_6pt(twistor)
        if H is None or H == 0: continue
        r = K/H
        ratios.append(r)
        
    import numpy as np
    vals = [float(x) for x in ratios]
    print_log(f"All points (n={len(vals)}): Mean={np.mean(vals):.6f}, Std/Mean={np.std(vals)/np.mean(vals):.6f}")
    
    # B2
    print_log("--- B2: Permutation Symmetry ---")
    sym_diffs = []
    for i in range(5):
        twistor = MomentumTwistor(n=6, seed=i+2000)
        diff = check_permutation_symmetry(twistor)
        if diff is not None: sym_diffs.append(diff)
    if sym_diffs:
        print_log(f"Avg Symmetry Violation: {np.mean(sym_diffs):.2e}")
        if np.mean(sym_diffs) < 1e-10:
            print_log("[PASS] Hodges is symmetric")
        else:
            print_log("[FAIL] Hodges is NOT symmetric")
            
    # B3
    print_log("--- B3: Soft Limit Trend ---")
    soft_ratios = []
    for i in range(5):
        r = check_soft_limit(seed=int(i+3000))
        if r is not None: soft_ratios.append(r)
    if soft_ratios:
        print_log(f"Ratio in Soft Limit: {np.mean(soft_ratios):.6f}")
    
if __name__ == "__main__":
    run_diagnostics()
