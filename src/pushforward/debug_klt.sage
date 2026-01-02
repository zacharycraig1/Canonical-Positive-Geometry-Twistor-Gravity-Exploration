import sys
import os
from sage.all import *

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_twistor
from src.chy_oracle.klt import parke_taylor_6pt_mhv, klt_momentum_kernel_6pt, gravity_6pt_mhv_klt
from src.chy_oracle.amplitude_spinor import mandelstam_s, hodges_npt_mhv_canonical

def debug():
    tw = sample_twistor(seed=100, n=6)
    ls = [tw.get_lambda(i) for i in range(6)]
    lts = []
    for i in range(6):
        lt = tw.get_tilde_lambda(i)
        if lt is None:
             print("Singular")
             return
        lts.append(lt)
    
    def m_s(t, i, j): return mandelstam_s(ls, lts, i, j)
    
    # Hodges
    H, status = hodges_npt_mhv_canonical(ls, lts, negative_indices=(0, 1))
    print(f"Hodges: {H}")

    # KLT
    K, status = gravity_6pt_mhv_klt(tw, m_s)
    print(f"KLT: {K}")
    
    if K != 0:
        print(f"Ratio: {H/K}")

    # alpha = (1, 2, 3)
    # beta = (1, 2, 3)
    alpha = (1, 2, 3)
    beta = (1, 2, 3)
    
    order_a = [0, 1, 2, 3, 4, 5]
    Aa = parke_taylor_6pt_mhv(tw, order_a)
    
    order_b = [0, 1, 2, 3, 5, 4] # swapped 5,4
    Ab = parke_taylor_6pt_mhv(tw, order_b)
    
    S = klt_momentum_kernel_6pt(alpha, beta, tw, m_s)
    
    print(f"Aa: {Aa}")
    print(f"Ab: {Ab}")
    print(f"S: {S}")
    print(f"Term (123|123): {Aa*S*Ab}")

if __name__ == "__main__":
    debug()

