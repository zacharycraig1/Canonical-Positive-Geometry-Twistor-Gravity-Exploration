import sys
from sage.all import *
import os

sys.path.append(os.path.join(os.getcwd(), 'src'))

from kinematics.spinors import SpinorKinematics
from chy.mhv_solution import sigma_mhv
from chy.scattering_eqs import detprime_phi
from chy.integrands import mhv_gravity_amplitude

root_dir = os.getcwd()
hodges_path = os.path.join(root_dir, 'src', 'hodges.sage')
load(hodges_path)

def check_symmetry():
    print("Checking Permutation Symmetry of M_CHY...")
    
    n = 6
    mt = MomentumTwistor(n, seed=100)
    lambdas = [mt.get_lambda(i) for i in range(n)]
    tilde_lambdas = [mt.get_tilde_lambda(i) for i in range(n)]
    
    theta = vector(QQ, [1, 0])
    eta = vector(QQ, [0, 1])
    chi = vector(QQ, [13, 7])
    
    # Base case: 0,1 negative.
    # We will swap particles 2 and 4.
    # 2 is in deletion set (0,1,2). 4 is kept (3,4,5).
    # If M is invariant, this must pass.
    
    # Calc M_original
    k1 = SpinorKinematics(n, lambdas, tilde_lambdas)
    sig1 = sigma_mhv(k1, theta, eta, chi)
    det1 = detprime_phi(sig1, k1)
    m1 = mhv_gravity_amplitude(sig1, k1, det1)
    
    print(f"M_original: {m1}")
    
    # Calc M_swapped (swap 2 and 4)
    L2 = list(lambdas)
    Lt2 = list(tilde_lambdas)
    L2[2], L2[4] = L2[4], L2[2]
    Lt2[2], Lt2[4] = Lt2[4], Lt2[2]
    
    k2 = SpinorKinematics(n, L2, Lt2)
    sig2 = sigma_mhv(k2, theta, eta, chi)
    det2 = detprime_phi(sig2, k2)
    m2 = mhv_gravity_amplitude(sig2, k2, det2)
    
    print(f"M_swapped:  {m2}")
    
    if m1 == m2:
        print("Symmetry 2 <-> 4 CONFIRMED.")
    else:
        print("Symmetry 2 <-> 4 BROKEN.")
        ratio = m1/m2
        print(f"Ratio: {ratio}")

if __name__ == "__main__":
    check_symmetry()

