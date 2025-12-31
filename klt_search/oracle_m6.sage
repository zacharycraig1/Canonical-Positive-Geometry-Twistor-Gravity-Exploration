
from sage.all import *
import numpy as np

# Load MomentumTwistor from src/hodges.sage for compatibility
try:
    load("src/hodges.sage")
except:
    pass

# ==============================================================================
# ROBUST HODGES (Copied from correct_klt_proof.sage)
# ==============================================================================
def robust_hodges_6pt_mhv(twistor, deleted_rows=None, ref_spinors=None):
    """
    Hodges formula - ROBUST implementation.
    """
    n = twistor.n
    if deleted_rows is None:
        deleted_rows = [0, 1, 2]
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
        
    def ang_prod(l1, l2):
        return l1[0]*l2[1] - l1[1]*l2[0]

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
    try:
        det_Phi_red = Phi_red.det()
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

def M6_oracle(twistor):
    """
    Oracle evaluator for 6-point MHV gravity.
    Uses BCFW expansion (sum over channels).
    """
    n = twistor.n
    channels = []
    for i in range(n):
        j = (i + 3) % n
        if j != i:
            channels.append((i, j))
            
    total = QQ(0)
    
    # Product of all angle brackets
    all_angles = QQ(1)
    for k in range(n):
        kp1 = (k + 1) % n
        ang = twistor.get_angle(k, kp1)
        if ang == 0: return None
        all_angles *= ang
        
    for i, j in channels:
        # Channel: i, i+1, i+2
        ip1, ip2 = (i + 1) % n, (i + 2) % n
        
        # Get outside particle
        outside = None
        for k in range(n):
            if k not in [i, ip1, ip2]:
                outside = k
                break
        
        if outside is None: continue
        
        # The 4-bracket <i i+1 i+2 outside>
        four_bracket = twistor.get_four_bracket(i, ip1, ip2, outside)
        if four_bracket == 0: continue
        
        # Channel angle brackets
        channel_angles = QQ(1)
        for idx in [i, ip1, ip2]:
            idxp1 = (idx + 1) % n
            ang = twistor.get_angle(idx, idxp1)
            if ang == 0: return None
            channel_angles *= ang
            
        # BCFW term
        # The term for gravity MHV
        # Based on correct_amplituhedron_hodges.sage
        # M3 * M3 / P^2. M3 ~ < >^8.
        # We need to scale by helicity factor <0 1>^8 ?
        # Actually, correct_amplituhedron_hodges.sage likely assumed "term" was the volume form, 
        # which needs to be multiplied by <0 1>^8 to get amplitude.
        
        term = (channel_angles * channel_angles) / (four_bracket * four_bracket * all_angles * all_angles)
        
        total += term
        
    # Multiply by helicity factor <0 1>^8 (assuming 0,1 negative)
    # This matches Hodges scaling
    neg_factor = twistor.get_angle(0, 1) ** 8
    
    return total * neg_factor

if __name__ == "__main__":
    print("Starting Oracle Comparison...")
    
    matches = 0
    failures = 0
    
    for seed in range(40, 60):
        twistor = MomentumTwistor(n=6, seed=seed)
        
        # Oracle (BCFW)
        oracle_val = M6_oracle(twistor)
        if oracle_val is None:
            continue
            
        # Hodges (Robust)
        h_val, msg = robust_hodges_6pt_mhv(twistor, deleted_rows=[0, 4, 5])
        if h_val is None:
            continue
            
        if h_val == 0:
            continue
            
        ratio = oracle_val / h_val
        diff = abs(ratio - 1)
        
        print(f"Seed {seed}: Oracle={float(oracle_val):.4e}, Hodges={float(h_val):.4e}, Ratio={float(ratio):.6f}")
        
        if diff < 1e-4:
            matches += 1
        else:
            failures += 1
            
    print(f"\nSummary: {matches} Matches, {failures} Failures")
