#!/usr/bin/env sage
# =============================================================================
# FINAL VERIFICATION: KLT vs Hodges Equivalence (Spinor Helicity)
# =============================================================================
# Computes KLT gravity and Hodges reduced amplitude using valid 4D kinematics.
# Checks if the ratio is constant across random samples.
# =============================================================================

from sage.all import *
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

load('src/spinor_sampling.sage')
load('src/hodges_sigma.sage')
load('src/klt.sage') # for parke_taylor_6pt_mhv, klt_momentum_kernel_6pt

# -----------------------------------------------------------------------------
# Re-implement Hodges for Spinors (copy from diagnostic probe)
# -----------------------------------------------------------------------------
def ang_vec(a, b):
    return a[0] * b[1] - a[1] * b[0]

def compute_hodges_spinor(lambdas, tilde_lambdas):
    """Compute Hodges reduced amplitude (baseline deletion) with spinors."""
    n = 6
    # Reference spinors (fixed random choice for invariance)
    lambda_x = vector(QQ, [1, 0]) 
    lambda_y = vector(QQ, [0, 1])
    # Note: In real run, ensure <ix>, <iy> != 0. 
    # For generic lambdas, standard basis is usually fine, but safer to sample.
    # Let's use simple ones and return None if bad.
    
    Phi = matrix(QQ, n, n)
    
    def get_angle(i, j): return ang_vec(lambdas[i], lambdas[j])
    def get_square(i, j): return ang_vec(tilde_lambdas[i], tilde_lambdas[j])
        
    # Off-diagonal
    for i in range(n):
        for j in range(n):
            if i != j:
                ij_ang = get_angle(i, j)
                if ij_ang == 0: return None
                ij_sq = get_square(i, j)
                Phi[i, j] = ij_sq / ij_ang
    
    # Diagonal
    for i in range(n):
        lambda_i = lambdas[i]
        ang_ix = ang_vec(lambda_i, lambda_x)
        ang_iy = ang_vec(lambda_i, lambda_y)
        if ang_ix == 0 or ang_iy == 0: return None
        
        diag_sum = QQ(0)
        for j in range(n):
            if j == i: continue
            lambda_j = lambdas[j]
            ang_jx = ang_vec(lambda_j, lambda_x)
            ang_jy = ang_vec(lambda_j, lambda_y)
            if ang_jx == 0 or ang_jy == 0: return None
            
            contrib = Phi[i, j] * (ang_jx * ang_jy) / (ang_ix * ang_iy)
            diag_sum -= contrib
        Phi[i, i] = diag_sum
        
    rows_delete = [0, 1, 2]
    cols_delete = [3, 4, 5]
    rows_keep = [3, 4, 5]
    cols_keep = [0, 1, 2]
    
    Phi_minor = Phi[rows_keep, cols_keep]
    det_minor = Phi_minor.det()
    
    # c factors
    i, j, k = rows_delete
    c_rows = QQ(1) / (get_angle(i, j) * get_angle(j, k) * get_angle(k, i))
    
    r, s, t = cols_delete
    c_cols = QQ(1) / (get_angle(r, s) * get_angle(s, t) * get_angle(t, r))
    
    sigma = hodges_sigma(rows_delete, cols_delete, n)
    sign = -1 
    
    return sign * sigma * c_rows * c_cols * det_minor

# -----------------------------------------------------------------------------
# KLT Gravity Implementation
# -----------------------------------------------------------------------------
def compute_klt_gravity(lambdas, tilde_lambdas):
    """Compute KLT gravity amplitude for n=6 MHV."""
    n = 6
    
    # 1. Compute Mandelstam invariants s_ij
    s = matrix(QQ, n, n)
    for i in range(n):
        for j in range(n):
            # s_ij = <ij>[ij]
            ang = ang_vec(lambdas[i], lambdas[j])
            sq = ang_vec(tilde_lambdas[i], tilde_lambdas[j])
            s[i, j] = ang * sq
            
    # 2. Compute Parke-Taylor amplitudes (A_YM)
    # Need to pass angle bracket function or object to parke_taylor_6pt_mhv
    # klt.sage's parke_taylor_6pt_mhv takes 'twistor' object with get_angle
    
    class MockTwistor:
        def __init__(self, ls): 
            self.lambdas = ls
            self.n = len(ls)
        def get_angle(self, i, j): return ang_vec(self.lambdas[i], self.lambdas[j])
        
    mock_twistor = MockTwistor(lambdas)
    
    # 3. Use klt.sage functions
    # gravity_6pt_mhv_klt expects a momentum twistor object to compute s_ij internally
    # but our s_ij comes from spinors.
    # We should reconstruct the KLT sum manually using s matrix and PT function.
    
    from itertools import permutations
    
    # Fixed legs
    # Left A: fixed 0, 4, 5 (indices 1,5,6 in 1-based) -> permute {1,2,3} (indices 2,3,4)
    # Right A: fixed 0, 4, 5 -> permute {1,2,3}
    # (Matches klt.sage convention)
    
    perm_set = [1, 2, 3] # Particles 2,3,4 in 0-based indexing
    
    total_amp = 0
    
    # S_KLT kernel depends on s matrix
    # We can reuse klt_momentum_kernel_6pt if we can adapt it to take 's'
    # klt.sage: klt_momentum_kernel_6pt(alpha, beta, twistor)
    # It calls mandelstam_invariant(i, j, twistor)
    
    # Let's monkey-patch mandelstam_invariant or rewrite kernel logic
    
    def get_s(i, j): return s[i, j]
    
    for alpha in permutations(perm_set):
        # A_YM(5, 6, alpha, 1) -> (4, 5, alpha, 0) 0-based
        # Order: [4, 5] + list(alpha) + [0]
        order_L = [4, 5] + list(alpha) + [0]
        # parke_taylor_6pt_mhv(twistor, order)
        amp_L = parke_taylor_6pt_mhv(mock_twistor, order_L)
        
        for beta in permutations(perm_set):
            # Kernel S[alpha|beta]
            # S = prod_{i=2}^{4} (s_{1, alpha_i} + sum_{j>i} theta(alpha_i, alpha_j in beta) * s_{alpha_i, alpha_j})
            # indices in formula are 1-based. 1 -> 0.
            # alpha is list of indices.
            
            # Re-implement kernel logic locally for clarity
            kernel = 1
            # Loop over elements in alpha (which are {1,2,3})
            # alpha = (a2, a3, a4) corresponding to i=2,3,4
            # We iterate k from 0 to 2
            
            for k in range(len(alpha)):
                # current particle index alpha[k]
                node_i = alpha[k]
                
                # Term 1: s_{1, alpha_i} -> s_{0, node_i}
                term = get_s(0, node_i)
                
                # Term 2: sum over j > k (future elements in alpha)
                for m in range(k + 1, len(alpha)):
                    node_j = alpha[m]
                    
                    # Check ordering in beta
                    # theta = 1 if node_i appears AFTER node_j in beta
                    pos_i_beta = beta.index(node_i)
                    pos_j_beta = beta.index(node_j)
                    
                    theta = 1 if pos_i_beta > pos_j_beta else 0
                    
                    term += theta * get_s(node_i, node_j)
                
                kernel *= term
            
            # A_YM(1, beta, 5, 6) -> (0, beta, 4, 5)
            order_R = [0] + list(beta) + [4, 5]
            amp_R = parke_taylor_6pt_mhv(mock_twistor, order_R)
            
            term = amp_L * kernel * amp_R
            total_amp += term
            
    return total_amp

# -----------------------------------------------------------------------------
# Main Verification Loop
# -----------------------------------------------------------------------------
def verify_equivalence():
    print("="*70)
    print("FINAL VERIFICATION: KLT vs Hodges (Spinor Helicity)")
    print("="*70)
    
    ratios = []
    
    for k in range(10):
        # Sample kinematics
        res = sample_spinor_helicity_conserving(n=6, seed=k*100)
        if res is None: continue
        lambdas, tilde_lambdas = res
        
        # Compute Hodges REDUCED
        H_red = compute_hodges_spinor(lambdas, tilde_lambdas)
        if H_red is None or H_red == 0: 
            print(f"Sample {k}: Hodges failed or zero")
            continue
            
        # Compute KLT (FULL AMPLITUDE)
        # Assumes negative helicity on legs 0 and 1
        K_full = compute_klt_gravity(lambdas, tilde_lambdas)
        
        # Hypothesis: KLT = C * <0 1>^8 * H_red
        # Check ratio = K_full / (H_red * <0 1>^8)
        
        ang01 = ang_vec(lambdas[0], lambdas[1])
        if ang01 == 0:
            print(f"Sample {k}: <0 1> is zero")
            continue
            
        ratio = K_full / (H_red * (ang01**8))
        ratios.append(ratio)
        
        print(f"Sample {k}: Ratio K / (H_red * <0 1>^8) = {ratio}")
        
    unique_ratios = list(set(ratios))
    print("\nSummary:")
    print(f"Tested {len(ratios)} samples.")
    print(f"Unique ratios: {unique_ratios}")
    
    if len(unique_ratios) == 1:
        print("\n*** DISCOVERY CONFIRMED ***")
        print(f"Constant Ratio: {unique_ratios[0]}")
        print("KLT Gravity is EQUIVALENT to Hodges Reduced Amplitude (up to this constant).")
    else:
        print("\n*** DISCREPANCY ***")
        print("Ratio is not constant.")

if __name__ == '__main__':
    verify_equivalence()

