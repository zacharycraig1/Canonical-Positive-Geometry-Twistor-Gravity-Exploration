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

def analyze_divisor_table():
    print("Phase C1: Generating Divisor Behavior Table...")
    
    R = PolynomialRing(QQ, 'eps')
    eps = R.gen()
    
    pairs = [(i, j) for i in range(6) for j in range(i+1, 6)]
    
    print(f"{'Divisor':<10} | {'Type':<10} | {'ord(M)':<8} | {'ord(N)':<8} | {'Notes'}")
    print("-" * 60)
    
    for (a, b) in pairs:
        # Construct generic deformation <ab> -> eps
        # lambda_b = lambda_a + eps * random_vec
        
        # Random base
        # Ensure generic enough to avoid other brackets vanishing
        set_random_seed(a*10 + b)
        
        Z_list_poly = []
        for k in range(6):
            Z_list_poly.append([QQ(ZZ.random_element(1, 100)) for _ in range(4)])
            
        # Fix lambda_a and lambda_b
        # Let lambda_a be random (already set)
        # Set lambda_b = lambda_a + eps * (1, 1) (or random)
        
        # We need to construct vectors in R
        lambdas = []
        for k in range(6):
            if k == b:
                # lambda_b = lambda_a + eps * nu
                # Pick nu not collinear to lambda_a
                nu = vector(R, [1, 2]) # Simple choice
                La = vector(R, [Z_list_poly[a][0], Z_list_poly[a][1]])
                Lb = La + eps * nu
                lambdas.append(Lb)
                # Update Z list for consistency (though we only use lambdas/tildes)
                Z_list_poly[b][0] = Lb[0]
                Z_list_poly[b][1] = Lb[1]
            else:
                lambdas.append(vector(R, [Z_list_poly[k][0], Z_list_poly[k][1]]))
                
        # Tilde lambdas
        def get_angle(i, j):
            return lambdas[i][0]*lambdas[j][1] - lambdas[i][1]*lambdas[j][0]
            
        tilde_lambdas = []
        for i in range(6):
            im1 = (i - 1) % 6
            ip1 = (i + 1) % 6
            
            mu_i = vector(R, [Z_list_poly[i][2], Z_list_poly[i][3]])
            mu_im1 = vector(R, [Z_list_poly[im1][2], Z_list_poly[im1][3]])
            mu_ip1 = vector(R, [Z_list_poly[ip1][2], Z_list_poly[ip1][3]])
            
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
                
            def angle(self, i, j):
                return self.lambdas[i][0]*self.lambdas[j][1] - self.lambdas[i][1]*self.lambdas[j][0]
                
            def square(self, i, j):
                return self.tilde_lambdas[i][0]*self.tilde_lambdas[j][1] - self.tilde_lambdas[i][1]*self.tilde_lambdas[j][0]
                
            def s(self, i, j):
                return self.angle(i, j) * self.square(j, i)

        k_eps = GenericKinematics(6, lambdas, tilde_lambdas)
        
        # References
        def r_val(): return QQ(ZZ.random_element(1, 100))
        theta = vector(R, [r_val(), r_val()])
        eta = vector(R, [r_val(), r_val()])
        chi = vector(R, [1, r_val()])

        try:
            sigmas = sigma_mhv(k_eps, theta, eta, chi)
            det_p = detprime_phi(sigmas, k_eps)
            m_chy = mhv_gravity_amplitude(sigmas, k_eps, det_p)
            
            # Denominator D
            d_all = R(1)
            for i in range(6):
                for j in range(i+1, 6):
                    d_all *= k_eps.angle(i, j)
            d_cyc = R(1)
            for i in range(6):
                d_cyc *= k_eps.angle(i, (i+1)%6)
            ang01 = k_eps.angle(0, 1)
            D_val = d_all * d_cyc / (ang01**2)
            
            N_val = m_chy * D_val
            
            def valuation(expr):
                if expr == 0: return 999
                num = expr.numerator()
                den = expr.denominator()
                vn = num.valuation() if hasattr(num, 'valuation') else 0
                vd = den.valuation() if hasattr(den, 'valuation') else 0
                return vn - vd
                
            ord_M = valuation(m_chy)
            ord_N = valuation(N_val)
            
            is_cyclic = (b == (a+1)%6) or (a == (b+1)%6)
            pair_type = "Cyclic" if is_cyclic else "Non-Cyc"
            if (a,b) == (0,1): pair_type = "Special"
            
            note = ""
            if ord_N > 0: note = f"Vanishes (deg {ord_N})"
            
            print(f"<{a}{b}>      | {pair_type:<10} | {ord_M:<8} | {ord_N:<8} | {note}")
            
        except Exception as e:
            print(f"<{a}{b}>      | Error: {e}")

if __name__ == "__main__":
    analyze_divisor_table()



