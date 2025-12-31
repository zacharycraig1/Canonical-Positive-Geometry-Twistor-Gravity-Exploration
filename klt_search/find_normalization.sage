
from sage.all import *
import sys
import os
sys.path.append(os.getcwd())

# Load dependencies
load("src/hodges.sage")
from itertools import permutations

def ts():
    import time
    return time.strftime("%H:%M:%S")

def log_message(msg):
    print(f"[{ts()}] {msg}")

# KLT Kernel (Standard Pivot 1)
# S[alpha|beta] p1
def klt_kernel_S3(alpha, beta, twistor):
    # alpha, beta are permutations of {2,3,4} (indices 1,2,3)
    # Pivot is 1 (index 0)
    
    # Map values to indices for theta function
    pos_beta = {val: idx for idx, val in enumerate(beta)}
    
    def theta(a, b):
        return 1 if pos_beta[a] > pos_beta[b] else 0
        
    def s(i, j):
        val = mandelstam_invariant(twistor, i, j)
        if val is None: raise ValueError(f"s_{i}_{j} is None")
        return val

    # Formula from BDPR 9811140 eq 2.11 for n=6 (N=3 in sum)
    # S[i1,i2,i3 | j1,j2,j3]
    # alpha = [i1, i2, i3]
    # beta = [j1, j2, j3]
    
    i1, i2, i3 = alpha
    
    # Term 1: s(1, i1) + theta(i1, i2)*s(i1, i2) + theta(i1, i3)*s(i1, i3)
    term1 = s(0, i1) 
    if theta(i1, i2): term1 += s(i1, i2)
    if theta(i1, i3): term1 += s(i1, i3)
    
    # Term 2: s(1, i2) + theta(i2, i3)*s(i2, i3)
    term2 = s(0, i2)
    if theta(i2, i3): term2 += s(i2, i3)
    
    # Term 3: s(1, i3)
    term3 = s(0, i3)
    
    return term1 * term2 * term3

def get_pt(twistor, order):
    # A_MHV = <0 1>^4 / prod <i i+1>
    denom = QQ(1)
    for k in range(6):
        denom *= twistor.get_angle(order[k], order[(k+1)%6])
    if denom == 0: return None
    num = twistor.get_angle(0, 1)**4
    return num / denom

def find_exact_match():
    log_message("Searching for Exact KLT/Hodges Match (Variations)...")
    
    # Use moment curve for stability
    def get_mc(seed):
        seed = int(seed) # Ensure integer
        t = [QQ(i+1) for i in range(6)]
        import random
        random.seed(seed)
        t = [x + QQ(random.randint(1,100))/1000 for x in t]
        t.sort()
        Z = []
        for val in t:
             Z.append(vector(QQ, [1, val, val**2, val**3]))
        return MomentumTwistor(n=6, Z=Z)

    # Prepare data
    samples = []
    for i in range(5):
        samples.append(get_mc(i+500))
        
    # Variations
    variations = [
        ("Standard (Swap 4,5)", True, lambda a,b: 1, [0, 1, 2]),
        ("Standard (Del 0,4,5)", True, lambda a,b: 1, [0, 4, 5]),
        ("No Swap (Del 0,4,5)", False, lambda a,b: 1, [0, 4, 5]),
    ]
    
    for label, do_swap, sign_func, del_set in variations:
        log_message(f"Testing: {label} [Del {del_set}]")
        ratios = []
        
        for twistor in samples:
            # Hodges
            H_res = hodges_6pt_mhv(twistor, deletion_set=del_set)
            H = H_res[0] if isinstance(H_res, tuple) else H_res
            if H is None or H == 0: continue
            
            # KLT
            permuted = [1, 2, 3] # indices
            perms = list(permutations(permuted))
            
            M_KLT = QQ(0)
            
            for alpha in perms:
                # Left: 1, alpha, n-1, n -> 0, alpha, 4, 5
                order_L = [0] + list(alpha) + [4, 5]
                AL = get_pt(twistor, order_L)
                if AL is None: continue
                
                for beta in perms:
                    # Right: 1, beta, n, n-1 -> 0, beta, 5, 4 (if swap)
                    if do_swap:
                        order_R = [0] + list(beta) + [5, 4]
                    else:
                        order_R = [0] + list(beta) + [4, 5]
                        
                    AR = get_pt(twistor, order_R)
                    if AR is None: continue
                    
                    try:
                        S = klt_kernel_S3(list(alpha), list(beta), twistor)
                    except ValueError: continue
                    
                    sign = sign_func(alpha, beta)
                    M_KLT += sign * AL * S * AR
            
            if M_KLT != 0:
                ratios.append(M_KLT / H)
        
        # Analyze ratios
        if not ratios: continue
        first = ratios[0]
        # Check variance
        rel_diffs = [abs((r - first)/first) for r in ratios]
        max_diff = max(rel_diffs)
        
        log_message(f"  First: {float(first):.6f}")
        log_message(f"  Max Rel Diff: {float(max_diff):.2e}")
        
        if max_diff < 1e-5:
            log_message("  [SUCCESS] Constant Ratio!")
            
    log_message("Note: Variation of ~12% suggests the ratio depends on cross-ratios (weight 0 functions).")
    log_message("Next step: Regress Ratio vs {u1, u2, u3} to find the exact factor.")

if __name__ == "__main__":
    find_exact_match()
        
if __name__ == "__main__":
    find_exact_match()
