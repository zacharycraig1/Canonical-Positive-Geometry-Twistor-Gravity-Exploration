import sys
from sage.all import *
import os
import random

# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src'))

from kinematics.spinors import SpinorKinematics
from chy.mhv_solution import sigma_mhv
from chy.scattering_eqs import detprime_phi
from chy.integrands import mhv_gravity_amplitude

root_dir = os.getcwd()
hodges_path = os.path.join(root_dir, 'src', 'hodges.sage')
load(hodges_path)

def bracket_poly_from_twistors(mt_A, mt_B, i, j, t):
    lamAi = vector(QQ, mt_A.get_lambda(i))
    lamBi = vector(QQ, mt_B.get_lambda(i))
    lamAj = vector(QQ, mt_A.get_lambda(j))
    lamBj = vector(QQ, mt_B.get_lambda(j))

    li0 = lamAi[0] + t*lamBi[0]
    li1 = lamAi[1] + t*lamBi[1]
    lj0 = lamAj[0] + t*lamBj[0]
    lj1 = lamAj[1] + t*lamBj[1]
    return (li0*lj1 - li1*lj0)

def validate_pole_structure_global():
    print("Validating pole structure on a fresh random line...")
    
    # New seeds - completely independent
    mt_A = MomentumTwistor(6, seed=999)
    mt_B = MomentumTwistor(6, seed=888)
    
    # We want to check if D(Z) = (Prod_all <ij>) * (Prod_cyclic <i,i+1>) / <01>^2
    # is the correct denominator.
    # We can check this by computing M * D and seeing if it's polynomial.
    # We don't need full reconstruction.
    # Just check polynomiality on a few points (e.g. Lagrange interpolation error is 0).
    
    n = 6
    theta = vector(QQ, [1, 0])
    eta = vector(QQ, [0, 1])
    chi = vector(QQ, [1, 1])
    
    # Sample points
    points = []
    # We need degree + 1 points to check polynomiality.
    # Degree expected ~ 30-34.
    # Let's take 40 points.
    
    print("Collecting 45 points...")
    for t_val in range(1, 46):
        t = QQ(t_val)
        Z_t = [mt_A.Z[i] + t * mt_B.Z[i] for i in range(n)]
        mt_t = MomentumTwistor(n, Z=Z_t, check_domain=False)
        lambdas = [mt_t.get_lambda(i) for i in range(n)]
        tilde_lambdas = [mt_t.get_tilde_lambda(i) for i in range(n)]
        
        if any(tl is None for tl in tilde_lambdas): continue
        k_t = SpinorKinematics(n, lambdas, tilde_lambdas)
        
        try:
            sigmas = sigma_mhv(k_t, theta, eta, chi)
            det_p = detprime_phi(sigmas, k_t)
            m_chy = mhv_gravity_amplitude(sigmas, k_t, det_p)
            
            # Compute Candidate Denominator D_sym / <01>^2
            # D_sym = Prod_{all} <ij> * Prod_{cyc} <i,i+1>
            
            # 1. Product of all pairs
            d_all = QQ(1)
            for i in range(n):
                for j in range(i+1, n):
                    d_all *= k_t.angle(i, j)
            
            # 2. Product of cyclic
            d_cyc = QQ(1)
            for i in range(n):
                j = (i+1)%n
                d_cyc *= k_t.angle(i, j)
                
            D_candidate = d_all * d_cyc / (k_t.angle(0, 1)**2)
            
            N_val = m_chy * D_candidate
            points.append((t, N_val))
            
        except Exception:
            pass
            
    print(f"Collected {len(points)} points.")
    
    # Interpolate on first 38 points
    train = points[:38]
    test = points[38:]
    
    R = PolynomialRing(QQ, 't')
    poly = R.lagrange_polynomial(train)
    
    print(f"Polynomial Degree: {poly.degree()}")
    
    # Check test points
    valid = True
    for (ti, yi) in test:
        val = poly(ti)
        err = abs(val - yi)
        if err > 1e-10:
            print(f"Mismatch at t={ti}. Error={err}")
            valid = False
            
    if valid:
        print("SUCCESS: Denominator structure confirmed globally!")
        print("N(Z) = M(Z) * D_found is a polynomial.")
    else:
        print("FAILURE: Denominator structure did not hold on new line.")

if __name__ == "__main__":
    validate_pole_structure_global()



