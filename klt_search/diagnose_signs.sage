
from sage.all import *
import numpy as np
from itertools import permutations

load("src/hodges.sage")

# Use robust Hodges from diagnose_variation (copied here for self-containment)
def robust_hodges_6pt_mhv(twistor, deleted_rows=None):
    n = twistor.n
    if deleted_rows is None: deleted_rows = [0, 4, 5] # Default to fixed legs
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

# KLT Helpers
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

def get_perm_sign(p):
    # p is tuple/list of indices
    # Bubble sort count
    lst = list(p)
    swaps = 0
    for i in range(len(lst)):
        for j in range(0, len(lst)-i-1):
            if lst[j] > lst[j+1]:
                lst[j], lst[j+1] = lst[j+1], lst[j]
                swaps += 1
    return 1 if swaps % 2 == 0 else -1

def sample_moment_curve(n=6, seed=None):
    if seed is not None: np.random.seed(seed)
    t = sorted([np.random.uniform(0, 10) for _ in range(n)])
    Z = []
    for val in t:
        Z.append(vector(QQ, [1, val, val**2, val**3]))
    return Z

def run_sign_diagnosis():
    print("Checking KLT sign conventions (Moment Curve)...")
    
    perm_indices = [1, 2, 3]
    perms = list(permutations(perm_indices))
    
    # Precompute signs for perms
    p_signs = {p: get_perm_sign(p) for p in perms}
    
    # Ansatzes
    ansatzes = [
        ("Unity", lambda a, b: 1),
        ("Sign(Alpha)", lambda a, b: p_signs[a]),
        ("Sign(Beta)", lambda a, b: p_signs[b]),
        ("Sign(Alpha)*Sign(Beta)", lambda a, b: p_signs[a] * p_signs[b]),
        ("Minus Unity", lambda a, b: -1)
    ]
    
    results = {name: [] for name, _ in ansatzes}
    
    for seed in range(1000, 1010): # 10 points
        Z = sample_moment_curve(n=6, seed=seed)
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=False)
        H = robust_hodges_6pt_mhv(twistor)
        if H is None or H == 0: continue
        
        # Compute KLT terms
        terms = {} # (alpha, beta) -> val
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
                terms[(alpha, beta)] = AL * S * AR
        
        for name, sign_func in ansatzes:
            K = QQ(0)
            for (alpha, beta), val in terms.items():
                s = sign_func(alpha, beta)
                K += s * val
            
            ratio = K / H
            results[name].append(float(ratio))
            
    # Analyze
    for name, ratios in results.items():
        if not ratios: continue
        mean = np.mean(ratios)
        std = np.std(ratios)
        cv = std / abs(mean) if mean != 0 else float('inf')
        print(f"{name:<25}: Mean={mean:.4e}, CV={cv:.4e}")
        
if __name__ == "__main__":
    run_sign_diagnosis()

