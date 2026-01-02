
from sage.all import *
import numpy as np
from itertools import permutations, combinations

load("src/hodges.sage")

# Use robust Hodges
def robust_hodges_6pt_mhv(twistor, deleted_rows=None):
    n = twistor.n
    if deleted_rows is None: deleted_rows = [0, 4, 5]
    deleted_rows = sorted(deleted_rows)
    rows_to_keep = [i for i in range(n) if i not in deleted_rows]
    cols_to_keep = rows_to_keep
    
    Phi = matrix(QQ, n, n)
    x_idx, y_idx = deleted_rows[0], deleted_rows[1]
    lambda_x = vector(QQ, [twistor.Z[x_idx][0], twistor.Z[x_idx][1]])
    lambda_y = vector(QQ, [twistor.Z[y_idx][0], twistor.Z[y_idx][1]])
    
    lambdas = []
    for i in range(n):
        lambdas.append(vector(QQ, [twistor.Z[i][0], twistor.Z[i][1]]))
        
    def ang_prod(l1, l2): return l1[0]*l2[1] - l1[1]*l2[0]

    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = twistor.get_angle(i, j)
                if ij_ang == 0: return None
                ij_sq = twistor.get_square(i, j)
                if ij_sq is None: return None
                Phi[i, j] = ij_sq / ij_ang
                
    for i in range(n):
        ix_ang = ang_prod(lambdas[i], lambda_x)
        iy_ang = ang_prod(lambdas[i], lambda_y)
        if ix_ang == 0 or iy_ang == 0:
            if i in rows_to_keep: return None
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
    except: return None
    
    r1, r2, r3 = deleted_rows
    norm_factor = QQ(1)
    norm_factor *= twistor.get_angle(r1, r2)
    norm_factor *= twistor.get_angle(r2, r3)
    norm_factor *= twistor.get_angle(r3, r1)
    if norm_factor == 0: return None
    norm_factor = norm_factor ** 2
    
    det_prime_Phi = det_Phi_red / norm_factor
    neg_factor = twistor.get_angle(0, 1) ** 8
    
    return det_prime_Phi * neg_factor

# KLT
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

# Cross Ratios
def get_cross_ratios(twistor):
    # Standard DCI cross ratios for n=6
    # u1 = <1234><4561> / <1245><3461>
    # u2 = <2345><5612> / <2356><4512>
    # u3 = <3456><6123> / <3461><5623>
    
    def fb(i, j, k, l):
        return twistor.get_four_bracket(i-1, j-1, k-1, l-1)
        
    u1_num = fb(1,2,3,4) * fb(4,5,6,1)
    u1_den = fb(1,2,4,5) * fb(3,4,6,1)
    u1 = u1_num / u1_den if u1_den != 0 else None
    
    u2_num = fb(2,3,4,5) * fb(5,6,1,2)
    u2_den = fb(2,3,5,6) * fb(4,5,1,2)
    u2 = u2_num / u2_den if u2_den != 0 else None
    
    u3_num = fb(3,4,5,6) * fb(6,1,2,3)
    u3_den = fb(3,4,6,1) * fb(5,6,2,3)
    u3 = u3_num / u3_den if u3_den != 0 else None
    
    return u1, u2, u3

def identify_prefactor():
    print("Identifying Prefactor (Ratio vs Cross Ratios)...")
    
    data = []
    
    for seed in range(2000, 2020):
        twistor = MomentumTwistor(n=6, seed=seed)
        H = robust_hodges_6pt_mhv(twistor)
        K = calculate_klt_6pt(twistor)
        
        if H is None or H == 0: continue
        
        ratio = K / H
        u1, u2, u3 = get_cross_ratios(twistor)
        
        if u1 is None or u2 is None or u3 is None: continue
        
        data.append({
            'ratio': ratio,
            'u1': u1, 'u2': u2, 'u3': u3
        })
        
    print(f"Collected {len(data)} points.")
    if not data: return
    
    # Try candidates
    candidates = [
        ("1", lambda u1, u2, u3: 1),
        ("u1", lambda u1, u2, u3: u1),
        ("u2", lambda u1, u2, u3: u2),
        ("u3", lambda u1, u2, u3: u3),
        ("1/u1", lambda u1, u2, u3: 1/u1),
        ("u1*u2*u3", lambda u1, u2, u3: u1*u2*u3),
        ("1/(u1*u2*u3)", lambda u1, u2, u3: 1/(u1*u2*u3)),
        ("u1+u2+u3", lambda u1, u2, u3: u1+u2+u3),
        ("u1*u2+u2*u3+u3*u1", lambda u1, u2, u3: u1*u2+u2*u3+u3*u1)
    ]
    
    for name, func in candidates:
        ratios = []
        for d in data:
            val = func(d['u1'], d['u2'], d['u3'])
            if val == 0: continue
            # Check Ratio / Candidate
            r = d['ratio'] / val
            ratios.append(float(r))
            
        mean = np.mean(ratios)
        std = np.std(ratios)
        cv = std / abs(mean) if mean != 0 else float('inf')
        
        print(f"{name:<20}: Mean={mean:.4e}, CV={cv:.4e}")

if __name__ == "__main__":
    identify_prefactor()








