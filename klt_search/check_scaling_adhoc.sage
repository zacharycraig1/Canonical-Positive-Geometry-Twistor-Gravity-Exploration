
from sage.all import *
import sys
import os
sys.path.append(os.getcwd())

load("src/hodges.sage")
from itertools import permutations

# Copied from klt_search/diagnose_ranks.sage
def klt_momentum_kernel_6pt(alpha, beta, twistor):
    if len(alpha) != 3 or len(beta) != 3: return None
    pos_in_beta = {val: idx for idx, val in enumerate(beta)}
    def theta_beta(a, b):
        if a not in pos_in_beta or b not in pos_in_beta: return 0
        return 1 if pos_in_beta[a] > pos_in_beta[b] else 0
    
    s_0_a0 = mandelstam_invariant(twistor, 0, alpha[0])
    sum1 = s_0_a0
    if theta_beta(alpha[0], alpha[1]): sum1 += mandelstam_invariant(twistor, alpha[0], alpha[1])
    if theta_beta(alpha[0], alpha[2]): sum1 += mandelstam_invariant(twistor, alpha[0], alpha[2])
    
    s_0_a1 = mandelstam_invariant(twistor, 0, alpha[1])
    sum2 = s_0_a1
    if theta_beta(alpha[1], alpha[2]): sum2 += mandelstam_invariant(twistor, alpha[1], alpha[2])
    
    s_0_a2 = mandelstam_invariant(twistor, 0, alpha[2])
    return sum1 * sum2 * s_0_a2

def get_pt(twistor, order):
    # MHV PT: <0 1>^4 / prod <i i+1>
    denom = QQ(1)
    for k in range(6):
        denom *= twistor.get_angle(order[k], order[(k+1)%6])
    if denom == 0: return None
    num = twistor.get_angle(0, 1)**4
    return num / denom

def compute_klt(twistor):
    permuted = [1, 2, 3]
    perms = list(permutations(permuted))
    M_KLT = QQ(0)
    for alpha in perms:
        order_L = [0] + list(alpha) + [4, 5]
        AL = get_pt(twistor, order_L)
        if AL is None: continue
        for beta in perms:
            order_R = [0] + list(beta) + [5, 4]
            AR = get_pt(twistor, order_R)
            if AR is None: continue
            S = klt_momentum_kernel_6pt(list(alpha), list(beta), twistor)
            if S is None: continue
            M_KLT += AL * S * AR
    return M_KLT

def check_scaling():
    print("Checking Scaling Dimensions (Z -> 2*Z)...")
    
    # 1. Base point
    twistor = None
    for seed in range(100):
        t = MomentumTwistor(n=6, seed=seed)
        if t.domain_ok:
            twistor = t
            break
            
    if twistor is None:
        print("Could not find valid twistor")
        return

    H1 = hodges_6pt_mhv(twistor)
    if isinstance(H1, tuple): H1 = H1[0]
    
    K1 = compute_klt(twistor)
    
    if H1 is None or K1 is None:
        print("H1 or K1 is None")
        return
    
    print(f"Base: H={float(H1):.4e}, KLT={float(K1):.4e}")
    
    # 2. Scaled point
    Z_scaled = [z * 2 for z in twistor.Z]
    twistor2 = MomentumTwistor(n=6, Z=Z_scaled)
    
    H2 = hodges_6pt_mhv(twistor2)
    if isinstance(H2, tuple): H2 = H2[0]
    
    K2 = compute_klt(twistor2)
    
    print(f"Scaled: H={float(H2):.4e}, KLT={float(K2):.4e}")
    
    # 3. Ratios
    ratio_H = H2 / H1
    ratio_K = K2 / K1
    
    log2_H = log(ratio_H, 2)
    log2_K = log(ratio_K, 2)
    
    print(f"Scaling H: 2^{float(log2_H):.2f}")
    print(f"Scaling KLT: 2^{float(log2_K):.2f}")
    
    diff = log2_K - log2_H
    print(f"Difference in scaling power: {float(diff):.2f}")

if __name__ == "__main__":
    check_scaling()

