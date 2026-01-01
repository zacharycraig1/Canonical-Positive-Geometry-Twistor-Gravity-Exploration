
from sage.all import *
import sys
import os
sys.path.append(os.getcwd())

load("src/hodges.sage")
from itertools import permutations

def klt_kernel_S2(alpha, beta, twistor):
    # n=5. Perms of {1,2} (indices)
    # Pivot 0.
    # S[i1, i2 | j1, j2]
    
    # Formula:
    # Term 1: s(0, i1) + theta(i1, i2)s(i1, i2)
    # Term 2: s(0, i2)
    
    pos_beta = {val: idx for idx, val in enumerate(beta)}
    def theta(a, b): return 1 if pos_beta[a] > pos_beta[b] else 0
    def s(i, j): return mandelstam_invariant(twistor, i, j)
    
    i1, i2 = alpha
    
    t1 = s(0, i1)
    if theta(i1, i2): t1 += s(i1, i2)
    
    t2 = s(0, i2)
    
    return t1 * t2

def get_pt(twistor, order):
    n = len(order)
    denom = QQ(1)
    for k in range(n):
        denom *= twistor.get_angle(order[k], order[(k+1)%n])
    if denom == 0: return None
    # MHV n=5: <0 1>^4 (assuming 0,1 negative)
    num = twistor.get_angle(0, 1)**4
    return num / denom

def check_n5():
    print("Checking n=5 KLT vs Hodges...")
    
    # Hodges n=5
    # Remove 0,1,2. Keep 3,4.
    
    for i in range(5):
        # Moment curve n=5
        t = [QQ(k+1) for k in range(5)]
        import random
        random.seed(i)
        t = [x + QQ(random.randint(1,50))/500 for x in t]
        Z = [vector(QQ, [1, v, v**2, v**3]) for v in t]
        twistor = MomentumTwistor(n=5, Z=Z, check_domain=False)
        
        # Hodges
        # We need a custom Hodges call for n=5 or use existing if it supports n=5
        # The existing hodges_6pt_mhv is hardcoded for n=6 (deletion [0,1,2] -> keeps 3,4,5).
        # We need to adapt it.
        
        # Manually compute Hodges n=5
        ref = sample_reference_spinors(twistor)
        Phi = build_hodges_phi(twistor, ref)
        if isinstance(Phi, tuple): continue
        
        # Reduced det: delete 0,1,2. Keep 3,4.
        Phi_red = Phi[[3,4], [3,4]]
        det_red = Phi_red.det()
        
        norm = (twistor.get_angle(0,1)*twistor.get_angle(1,2)*twistor.get_angle(2,0))**2
        H = det_red / norm * twistor.get_angle(0,1)**8 # Helicity factor
        
        # KLT n=5
        # Sum over perms of {1, 2} (particles 2,3)
        # Fixed: 0, 3, 4 (particles 1, 4, 5)
        
        perms = list(permutations([1, 2]))
        M_KLT = QQ(0)
        
        for alpha in perms:
            # L: 0, alpha, 3, 4
            oL = [0] + list(alpha) + [3, 4]
            AL = get_pt(twistor, oL)
            
            for beta in perms:
                # R: 0, beta, 4, 3
                oR = [0] + list(beta) + [4, 3]
                AR = get_pt(twistor, oR)
                
                S = klt_kernel_S2(alpha, beta, twistor)
                
                M_KLT += AL * S * AR
                
        ratio = M_KLT / H
        print(f"Sample {i}: Ratio = {float(ratio):.6f}")

if __name__ == "__main__":
    check_n5()







