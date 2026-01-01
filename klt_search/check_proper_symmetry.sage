
from sage.all import *
import numpy as np
from itertools import permutations

load("src/hodges.sage")
load("klt_search/spinor_tools.sage")

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

def check_proper_symmetry():
    print("Checking Proper Spinor Permutation Symmetry...")
    
    for seed in range(3000, 3005):
        twistor = MomentumTwistor(n=6, seed=seed)
        
        # 1. Base result
        H0 = robust_hodges_6pt_mhv(twistor)
        if H0 is None: continue
        
        # 2. Extract spinors
        lambdas, t_lambdas = spinors_from_twistors(twistor)
        if None in t_lambdas: continue
        
        # 3. Permute 2 and 3 (indices)
        # Permutation: 0, 1, 3, 2, 4, 5
        perm_map = [0, 1, 3, 2, 4, 5]
        
        lambdas_new = [lambdas[p] for p in perm_map]
        t_lambdas_new = [t_lambdas[p] for p in perm_map]
        
        # 4. Reconstruct Twistors
        Z_new = twistors_from_spinors(lambdas_new, t_lambdas_new)
        if Z_new is None:
            print("Failed to reconstruct twistors")
            continue
            
        twistor_new = MomentumTwistor(n=6, Z=Z_new, check_domain=False)
        
        # 5. Calculate Permuted Amplitude
        H_perm = robust_hodges_6pt_mhv(twistor_new)
        if H_perm is None: continue
        
        ratio = H_perm / H0
        diff = abs(ratio - 1.0)
        
        print(f"Seed {seed}: H0={float(H0):.4e}, H_perm={float(H_perm):.4e}, Ratio={float(ratio):.4f}")

if __name__ == "__main__":
    check_proper_symmetry()







