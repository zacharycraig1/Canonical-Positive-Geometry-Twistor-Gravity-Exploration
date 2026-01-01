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

def debug_poles():
    print("Debugging Pole Orders for <23> vs <24>...")
    
    R = PolynomialRing(QQ, 'eps')
    eps = R.gen()
    
    def check_pair(a, b):
        # Deterministic setup
        set_random_seed(12345)
        
        # Base Z
        Z_base = []
        for k in range(6):
            Z_base.append([QQ(ZZ.random_element(1, 100)) for _ in range(4)])
            
        lambdas = []
        nu = vector(R, [1, 2])
        
        for k in range(6):
            if k == b:
                La = vector(R, [Z_base[a][0], Z_base[a][1]])
                Lb = La + eps * nu
                lambdas.append(Lb)
            else:
                lambdas.append(vector(R, [Z_base[k][0], Z_base[k][1]]))
                
        # Tildes
        def get_angle(i, j):
            return lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
            
        tilde_lambdas = []
        for i in range(6):
            im1 = (i - 1) % 6
            ip1 = (i + 1) % 6
            
            mu_i = vector(R, [Z_base[i][2], Z_base[i][3]])
            mu_im1 = vector(R, [Z_base[im1][2], Z_base[im1][3]])
            mu_ip1 = vector(R, [Z_base[ip1][2], Z_base[ip1][3]])
            
            ang_i_ip1 = get_angle(i, ip1)
            ang_ip1_im1 = get_angle(ip1, im1)
            ang_im1_i = get_angle(im1, i)
            
            denom = ang_im1_i * ang_i_ip1
            num = mu_im1 * ang_i_ip1 + mu_i * ang_ip1_im1 + mu_ip1 * ang_im1_i
            tilde_lambdas.append(num / denom)
            
        class GenericKinematics:
            def __init__(self, n, L, L_tilde):
                self.n = n
                self.lambdas = L
                self.tilde_lambdas = L_tilde
            def angle(self, i, j): return self.lambdas[i][0]*self.lambdas[j][1] - self.lambdas[i][1]*self.lambdas[j][0]
            def square(self, i, j): return self.tilde_lambdas[i][0]*self.tilde_lambdas[j][1] - self.tilde_lambdas[i][1]*self.tilde_lambdas[j][0]
            def s(self, i, j): return self.angle(i, j) * self.square(j, i)

        k_eps = GenericKinematics(6, lambdas, tilde_lambdas)
        
        # Refs
        def r_val(): return QQ(ZZ.random_element(1, 100))
        theta = vector(R, [r_val(), r_val()])
        eta = vector(R, [r_val(), r_val()])
        chi = vector(R, [1, r_val()])
        
        sigmas = sigma_mhv(k_eps, theta, eta, chi)
        det_p = detprime_phi(sigmas, k_eps)
        m_chy = mhv_gravity_amplitude(sigmas, k_eps, det_p)
        
        def valuation(expr):
            if expr == 0: return 999
            num = expr.numerator()
            den = expr.denominator()
            vn = num.valuation() if hasattr(num, 'valuation') else 0
            vd = den.valuation() if hasattr(den, 'valuation') else 0
            return vn - vd
            
        ord_val = valuation(m_chy)
        print(f"Pair <{a}{b}>: Order {ord_val}")
        
    check_pair(2, 3)
    check_pair(2, 4)

if __name__ == "__main__":
    debug_poles()






