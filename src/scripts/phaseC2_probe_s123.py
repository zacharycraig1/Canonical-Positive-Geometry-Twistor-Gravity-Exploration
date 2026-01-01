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

def probe_s123_pole():
    print("Phase C2: Probing s_123 pole behavior...")
    
    R = PolynomialRing(QQ, 'eps')
    eps = R.gen()
    
    # Construct kinematics where s_123 = eps
    # We can do this by taking generic kinematics and shifting one vector to satisfy P^2=eps.
    # Hard to do exactly in symbolic ring without introducing square roots.
    # Alternatively, use numerical probe with very small eps.
    # Or: P^2 = s12 + s13 + s23.
    # Set s12 = eps, s13=0, s23=0? No, s_ij must be physical.
    
    # Easier: Use random line Z(t) and find root of s_123(t).
    # Then expand around it.
    
    print("Finding root of s_123(t) on a line...")
    seed = 123
    mt_A = MomentumTwistor(6, seed=seed)
    mt_B = MomentumTwistor(6, seed=seed+1)
    
    # Define s_123(t)
    # s_123 = <12>[21] + <13>[31] + <23>[32]
    # This is a polynomial in t?
    # <ij> is deg 2. [ji] is deg 2?
    # No, <ij> is deg 1 in Z? No, quadratic in Z components.
    # Z(t) is linear. <ij>(t) is quadratic. [ji](t) is rational?
    # [ji] depends on tilde_lambda. tilde_lambda is rational in Z.
    # So s_123(t) is a rational function.
    
    # Let's just sample values and interpolate s_123(t) numerator?
    # Or just use Newton's method to find a root t* numerically.
    
    # We'll use high precision float.
    CIF = ComplexField(200)
    
    def get_s123(t_val):
        t_qq = QQ(real(t_val)) # approx
        # Reconstruct exactly if t is rational
        # But for root finding we need field.
        # Let's stick to QQ and find root in algebraic closure?
        pass

    # Let's just evaluate s_123 at t=0,1,2... interpolate, solve.
    points = []
    for k in range(10):
        t = QQ(k)
        Z_t = [mt_A.Z[i] + t * mt_B.Z[i] for i in range(6)]
        mt_t = MomentumTwistor(6, Z=Z_t, check_domain=False)
        L = [mt_t.get_lambda(i) for i in range(6)]
        Lt = [mt_t.get_tilde_lambda(i) for i in range(6)]
        if any(x is None for x in Lt): continue
        k_t = SpinorKinematics(6, L, Lt)
        
        s123 = k_t.s(0,1) + k_t.s(0,2) + k_t.s(1,2)
        points.append((t, s123))
        
    # Reconstruct rational function for s_123(t)
    # It might be complex.
    # Let's check if it crosses 0.
    
    # Actually, simpler check:
    # Does M blow up when s_123 is small?
    # Run a loop decreasing s_123.
    
    print("Numerical probe of s_123 -> 0")
    
    # We can force s_123 = 0 by solving constraint.
    # But checking M values near 0 is easier.
    
    # Scan t
    # Find t where s123 is small
    best_t = 0
    min_s = 999
    
    # Coarse scan
    for k in range(-100, 100):
        t = QQ(k)/10
        Z_t = [mt_A.Z[i] + t * mt_B.Z[i] for i in range(6)]
        mt_t = MomentumTwistor(6, Z=Z_t, check_domain=False)
        L = [mt_t.get_lambda(i) for i in range(6)]
        Lt = [mt_t.get_tilde_lambda(i) for i in range(6)]
        if any(x is None for x in Lt): continue
        k_t = SpinorKinematics(6, L, Lt)
        s123 = abs(k_t.s(0,1) + k_t.s(0,2) + k_t.s(1,2))
        if s123 < min_s:
            min_s = s123
            best_t = t
            
    print(f"Best t ~ {float(best_t):.4f} with s123 ~ {float(min_s):.4e}")
    
    # Refine
    t_curr = best_t
    for step in range(5):
        # Gradient descentish or just binary search
        pass
        
    # If s123 doesn't vanish, we can't test.
    # But s123 usually has roots.
    
    # Instead, let's use the known result:
    # Compute M at this t.
    # If M ~ 1/s123, then M * s123 ~ constant.
    # If M is regular, M ~ constant.
    
    # Let's verify M is finite.
    
    theta = vector(QQ, [1, 0])
    eta = vector(QQ, [0, 1])
    chi = vector(QQ, [1, 1])
    
    try:
        # Check M at the "small s123" point
        Z_t = [mt_A.Z[i] + best_t * mt_B.Z[i] for i in range(6)]
        mt_t = MomentumTwistor(6, Z=Z_t, check_domain=False)
        L = [mt_t.get_lambda(i) for i in range(6)]
        Lt = [mt_t.get_tilde_lambda(i) for i in range(6)]
        k_t = SpinorKinematics(6, L, Lt)
        
        sigmas = sigma_mhv(k_t, theta, eta, chi)
        det_p = detprime_phi(sigmas, k_t)
        m_val = mhv_gravity_amplitude(sigmas, k_t, det_p)
        
        print(f"M at s123={float(min_s):.4e} is {float(abs(m_val)):.4e}")
        
        # Now try to make s123 smaller?
        # Or just conclude based on magnitude.
        # If M ~ 10^30, likely pole.
        # If M ~ 10^0, likely finite.
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    probe_s123_pole()






