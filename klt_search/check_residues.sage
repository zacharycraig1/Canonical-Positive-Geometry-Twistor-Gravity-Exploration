
from sage.all import *
import sys
import os
sys.path.append(os.getcwd())

load("src/hodges.sage")
from itertools import permutations

# KLT Kernel (Standard Pivot 1)
def klt_kernel_S3(alpha, beta, twistor):
    pos_beta = {val: idx for idx, val in enumerate(beta)}
    def theta(a, b): return 1 if pos_beta[a] > pos_beta[b] else 0
    def s(i, j): return mandelstam_invariant(twistor, i, j)
    i1, i2, i3 = alpha
    t1 = s(0, i1) 
    if theta(i1, i2): t1 += s(i1, i2)
    if theta(i1, i3): t1 += s(i1, i3)
    t2 = s(0, i2)
    if theta(i2, i3): t2 += s(i2, i3)
    t3 = s(0, i3)
    return t1 * t2 * t3

def get_pt(twistor, order):
    denom = QQ(1)
    for k in range(6):
        denom *= twistor.get_angle(order[k], order[(k+1)%6])
    if denom == 0: return None
    num = twistor.get_angle(0, 1)**4
    return num / denom

def compute_klt(twistor):
    perms = list(permutations([1, 2, 3]))
    M_KLT = QQ(0)
    for alpha in perms:
        oL = [0] + list(alpha) + [4, 5]
        AL = get_pt(twistor, oL)
        if AL is None: continue
        for beta in perms:
            oR = [0] + list(beta) + [5, 4]
            AR = get_pt(twistor, oR)
            if AR is None: continue
            try:
                S = klt_kernel_S3(list(alpha), list(beta), twistor)
                M_KLT += AL * S * AR
            except: continue
    return M_KLT

def check_residues():
    print("Checking KLT/Hodges Ratio on Factorization Channels...")
    
    # We want to approach a pole.
    # For n=6, poles are s_ij, s_ijk.
    # s_123 = (p1+p2+p3)^2 -> <12>[12] + <13>[13] + <23>[23]
    # Hard to constrain random twistors to s_ijk=0 exactly without advanced geometry.
    # Instead, let's look at a soft collinear limit? s_12 -> 0.
    # s_12 -> 0 implies <12> -> 0 or [12] -> 0.
    
    # Construct a limit where <0 1> -> epsilon.
    # Z0 and Z1 close.
    
    for eps_pow in range(1, 10):
        epsilon = QQ(10)**(-eps_pow)
        
        # Base Z
        t = [QQ(k+1) for k in range(6)]
        import random
        random.seed(int(42))
        t = [x + QQ(random.randint(1,100))/1000 for x in t]
        Z = [vector(QQ, [1, v, v**2, v**3]) for v in t]
        
        # Modify Z1 to be Z0 + eps * (Z1_orig - Z0)
        # Actually just Z1 = Z0 + eps * V
        Z[1] = Z[0] + epsilon * (Z[1] - Z[0])
        
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=False)
        
        # Hodges
        H_res = hodges_6pt_mhv(twistor)
        H = H_res[0] if isinstance(H_res, tuple) else H_res
        
        # KLT
        K = compute_klt(twistor)
        
        if H == 0: continue
        
        ratio = K / H
        print(f"Eps = 1e-{eps_pow}: Ratio = {float(ratio):.6f}")

if __name__ == "__main__":
    check_residues()

