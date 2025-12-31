
from sage.all import *
import sys
import os
sys.path.append(os.getcwd())

load("src/hodges.sage")
from itertools import permutations

# Standard KLT kernel pivot 1
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

def find_factor():
    print("Gathering data for regression...")
    
    data = []
    
    # Generate 15 points
    for i in range(15):
        # Moment curve
        t = [QQ(k+1) for k in range(6)]
        import random
        random.seed(int(i+1000))
        t = [x + QQ(random.randint(1,100))/1000 for x in t]
        Z = [vector(QQ, [1, v, v**2, v**3]) for v in t]
        twistor = MomentumTwistor(n=6, Z=Z, check_domain=False)
        
        # Hodges (Del 0,1,2)
        H_res = hodges_6pt_mhv(twistor)
        H = H_res[0] if isinstance(H_res, tuple) else H_res
        
        # KLT
        perms = list(permutations([1,2,3]))
        M_KLT = QQ(0)
        for alpha in perms:
            oL = [0] + list(alpha) + [4, 5]
            AL = get_pt(twistor, oL)
            for beta in perms:
                oR = [0] + list(beta) + [5, 4] # Swap 4,5
                AR = get_pt(twistor, oR)
                S = klt_kernel_S3(list(alpha), list(beta), twistor)
                M_KLT += AL * S * AR
                
        if H == 0: continue
        ratio = M_KLT / H
        
        # Try to correlate with other cross ratios
        # Maybe 1 + u1 + u2 + u3?
        
        # Standard cross ratios:
        # u1 = s12 s45 / s123 s345 ? No, for 6pt.
        # u1 = s12 s45 / s123 s234 ...
        
        s = lambda i,j: mandelstam_invariant(twistor, i, j)
        s3 = lambda i,j,k: s(i,j) + s(j,k) + s(k,i)
        
        # OPE-like variables?
        
        data.append({
            'ratio': ratio,
            'idx': i
        })
        
        print(f"Sample {i}: R={float(ratio):.4f}")
        
    print("\n[FAIL] Simple regression failed. Trying constant search on n=5 again.")

if __name__ == "__main__":
    find_factor()

