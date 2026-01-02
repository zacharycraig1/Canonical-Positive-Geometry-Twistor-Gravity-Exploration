import sys
import os
from sage.all import *

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_twistor
from src.chy_oracle.amplitude_spinor import mandelstam_s

# Need a generic n-point parke_taylor
def parke_taylor_n_pt(twistor, order, neg_helicity=(0, 1)):
    n = len(order)
    if n != twistor.n:
        return None
        
    denom = QQ(1)
    for i in range(n):
        j = (i + 1) % n
        idx_i = order[i]
        idx_j = order[j]
        bracket = twistor.get_angle(idx_i, idx_j)
        if bracket == 0: return None
        denom *= bracket
        
    neg_a, neg_b = neg_helicity
    h_factor = twistor.get_angle(neg_a, neg_b)
    if h_factor == 0: return None
    
    return (h_factor ** 4) / denom

def verify_klt_n5():
    print("O1.c (n=5): Verifying KLT against permutation sum...")
    
    n = 5
    twistor = sample_twistor(seed=100, n=n)
    ls = [twistor.get_lambda(i) for i in range(n)]
    lts = []
    for i in range(n):
        lt = twistor.get_tilde_lambda(i)
        if lt is None:
            print("Singular")
            return
        lts.append(lt)
        
    def m_s(i, j): return mandelstam_s(ls, lts, i, j)
    
    # Formula for n=5
    # M_5 = (-1)^6 * sum_{alpha, beta in S2} A(1, alpha, 4, 5) S[alpha|beta] A(1, beta, 5, 4)
    # Fixed legs: 1(0), 4(3), 5(4).
    # Permutations of {2, 3} (indices 1, 2).
    # M5 = A(0, alpha, 3, 4) S A(0, beta, 4, 3)
    
    perms = [[1, 2], [2, 1]]
    
    # Kernel S[alpha|beta]
    # alpha = (a1, a2)
    # S = (s_{0, a1}) * (s_{0, a2} + theta(a1, a2; beta) s_{a1, a2})
    
    total_M5 = QQ(0)
    
    for alpha in perms:
        # A(0, alpha, 3, 4)
        order_a = [0] + alpha + [3, 4]
        A_a = parke_taylor_n_pt(twistor, order_a, neg_helicity=(0, 1))
        
        for beta in perms:
            # A(0, beta, 4, 3)
            order_b = [0] + beta + [4, 3]
            A_b = parke_taylor_n_pt(twistor, order_b, neg_helicity=(0, 1))
            
            # Kernel
            # 1. Term 1: s_{0, alpha[0]}
            val1 = m_s(0, alpha[0])
            
            # 2. Term 2: s_{0, alpha[1]} + theta * s_{alpha[0], alpha[1]}
            theta = 0
            # theta(alpha[0], alpha[1]; beta)
            # 1 if alpha[0] appears AFTER alpha[1] in beta? 
            # Check definition: theta(aj, ai) = 1 if aj before ai in beta (Wait, plan says 'aj before ai'?)
            # Plan O1.a: "theta(aj, ai) = 1 if aj appears BEFORE ai in the list beta? No, plan says 'alpha[j] appears BEFORE alpha[i]'? 
            # Let's check plan text: "theta(aj, ai) = 1 if aj appears *before* ai in beta".
            # Wait, usually it's "after".
            # Bern et al (1998): S[i1...ik | j1...jk] = prod_t (s_1it + sum_{q<t} theta(iq, it)_J s_iq it) where theta is regarding J ordering.
            # Let's try "after" first, as per some conventions (e.g. reverse ordering).
            # Plan says: "theta(alpha[j], alpha[i]) = 1 if alpha[j] appears *before* alpha[i] in beta" (Wait, re-read plan).
            # Plan says: "theta(alpha[j], alpha[i]) = 1 if alpha[j] appears **before** alpha[i] in beta".
            # Okay, let's use "before".
            
            # Check alpha[0] vs alpha[1] in beta.
            pos0 = beta.index(alpha[0])
            pos1 = beta.index(alpha[1])
            if pos0 < pos1:
                theta = 1
            else:
                theta = 0
            
            val2 = m_s(0, alpha[1]) + theta * m_s(alpha[0], alpha[1])
            
            S = val1 * val2
            
            total_M5 += A_a * S * A_b
            
    print(f"M5 (Calculated): {total_M5}")
    
    # Check permutation symmetry (swap 1 and 2, indices 1 and 2)
    # Oh wait, we summed over 1 and 2.
    # Check swapping 1 and 3? (indices 1 and 3)
    # Recompute with swapped indices.
    
    # Or just use Hodges N=5 as oracle (since we trust Hodges deletion invariance).
    from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical
    
    H5, status = hodges_npt_mhv_canonical(ls, lts, negative_indices=(0, 1))
    print(f"M5 (Hodges):     {H5}")
    
    ratio = H5 / total_M5
    print(f"Ratio (Hodges/KLT): {ratio}")
    
    if abs(ratio - 1) < 1e-9:
        print("PASS: KLT matches Hodges for N=5.")
    elif abs(ratio + 1) < 1e-9:
        print("PASS: KLT matches Hodges for N=5 (Sign flip).")
    else:
        print("FAIL: Mismatch.")

if __name__ == "__main__":
    verify_klt_n5()

