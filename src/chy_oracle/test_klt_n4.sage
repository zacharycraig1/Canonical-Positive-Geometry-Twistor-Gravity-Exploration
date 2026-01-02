import sys
import os
from sage.all import *

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_twistor
from src.chy_oracle.amplitude_spinor import mandelstam_s
from src.chy_oracle.klt import parke_taylor_6pt_mhv # We'll need a generic n-pt version

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

def verify_klt_n4():
    print("O1.c (n=4): Verifying KLT against ground truth formula...")
    
    # Ground truth: M4 = - s_12 A(1,2,3,4) A(1,2,4,3)
    # Using indices 0,1,2,3.
    # Fixed legs: 1 and n-1, n?
    # Usually: M4 = - s_12 A(1,2,3,4) A(1,2,4,3)
    # Let's check indices carefully.
    # 0, 1, 2, 3.
    # s_01 * A(0,1,2,3) * A(0,1,3,2) ?
    # Note: A(0,1,3,2) has swapped 2,3 (last two).
    
    n = 4
    twistor = sample_twistor(seed=100, n=n) # Seed 42 was singular for n=6, maybe 4 too
    ls = [twistor.get_lambda(i) for i in range(n)]
    lts = []
    for i in range(n):
        lt = twistor.get_tilde_lambda(i)
        if lt is None:
            print("Singular")
            return
        lts.append(lt)
        
    def m_s(i, j): return mandelstam_s(ls, lts, i, j)
    
    # Compute A(0,1,2,3)
    A_1234 = parke_taylor_n_pt(twistor, [0, 1, 2, 3], neg_helicity=(0, 1))
    
    # Compute A(0,1,3,2)
    A_1243 = parke_taylor_n_pt(twistor, [0, 1, 3, 2], neg_helicity=(0, 1))
    
    s_01 = m_s(0, 1)
    
    # Formula M4 = - s_01 * A(0,1,2,3) * A(0,1,3,2)
    M4_ref = - s_01 * A_1234 * A_1243
    
    print(f"M4 (Reference): {M4_ref}")
    
    # Now check if this matches general KLT kernel logic
    # KLT form: sum_{alpha, beta} A(1, alpha, n-1, n) S[alpha|beta] A(1, beta, n, n-1)
    # n=4. Fixed 1 (0), n-1 (2), n (3).
    # Permutations of {2..n-2} -> {1}.
    # alpha, beta in S_1 (just [1]).
    # M = A(0, 1, 2, 3) * S[1|1] * A(0, 1, 3, 2)
    # S[1|1] = - (s_01 + ...) ? No.
    # S[alpha|beta] formula:
    # Kernel for n=4 is just -s_01? Or s_01?
    # KLT kernel S[i|j] = k . p_i?
    
    # Using the formula in plan:
    # S[alpha|beta] = prod (s_{1, alpha} + ...)
    # Here pivot is 1 (index 0).
    # alpha = [1].
    # S = s_{0, 1}.
    
    # Wait, the KLT formula usually has M = (-1)^(n+1) ...
    # n=4 => (-1)^5 = -1.
    # So M = -1 * A * S * A'
    # M = - A(0,1,2,3) * s_01 * A(0,1,3,2).
    # Matches!
    
    # Now let's try to verify M4 is permutation symmetric?
    # Gravity should be symmetric in 0,1,2,3.
    # Let's swap 1 and 2 (indices 1 and 2).
    
    # A(0,2,1,3)
    A_1324 = parke_taylor_n_pt(twistor, [0, 2, 1, 3], neg_helicity=(0, 1))
    # A(0,2,3,1)
    A_1342 = parke_taylor_n_pt(twistor, [0, 2, 3, 1], neg_helicity=(0, 1))
    s_02 = m_s(0, 2)
    
    M4_swap = - s_02 * A_1324 * A_1342
    
    print(f"M4 (Swapped):   {M4_swap}")
    
    ratio = M4_ref / M4_swap
    print(f"Ratio: {ratio}")
    
    if abs(ratio - 1) < 1e-9:
        print("PASS: M4 KLT formula is permutation invariant.")
    else:
        print("FAIL: M4 KLT formula inconsistent.")

if __name__ == "__main__":
    verify_klt_n4()

