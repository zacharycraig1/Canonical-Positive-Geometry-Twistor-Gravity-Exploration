import sys
from sage.all import *
import os

# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src'))

from kinematics.spinors import SpinorKinematics
from chy.mhv_solution import sigma_mhv
from chy.scattering_eqs import detprime_phi
from chy.integrands import mhv_gravity_amplitude

# Load Hodges for MomentumTwistor structure
root_dir = os.getcwd()
hodges_path = os.path.join(root_dir, 'src', 'hodges.sage')
load(hodges_path)

def analyze_01_vanishing_order():
    """
    Phase C0: Determine exact vanishing order v_{01} of N(epsilon)
    as <01> -> epsilon.
    """
    print("Phase C0: Analyzing <01> vanishing order...")
    
    R = PolynomialRing(QQ, 'eps')
    eps = R.gen()
    
    seed = 42
    set_random_seed(seed)
    
    def rand_Z():
        # Avoid 0 to prevent accidental collinearity
        def nz():
            r = QQ.random_element(num_bound=10, den_bound=1)
            while r == 0: r = QQ.random_element(num_bound=10, den_bound=1)
            return r
        return [nz() for _ in range(4)]
        
    Z_list_poly = []
    
    # Z0: Fixed lambda=(1,0)
    z0 = rand_Z()
    z0[0], z0[1] = 1, 0
    
    # Z1: Fixed lambda=(1, eps) - Generic approach to <01>=0 (collinear)
    z1 = rand_Z()
    z1[0], z1[1] = 1, eps
    
    Z_list_poly.append(vector(R, z0))
    Z_list_poly.append(vector(R, z1))
    
    for i in range(2, 6):
        Z_list_poly.append(vector(R, rand_Z()))
        
    lambdas = []
    for i in range(6):
        lambdas.append(vector(R, [Z_list_poly[i][0], Z_list_poly[i][1]]))
    
    # Helper for angle brackets
    def get_angle_vec(L, i, j):
        return L[i][0]*L[j][1] - L[i][1]*L[j][0]
        
    tilde_lambdas = []
    for i in range(6):
        im1 = (i - 1) % 6
        ip1 = (i + 1) % 6
        
        mu_i = vector(R, [Z_list_poly[i][2], Z_list_poly[i][3]])
        mu_im1 = vector(R, [Z_list_poly[im1][2], Z_list_poly[im1][3]])
        mu_ip1 = vector(R, [Z_list_poly[ip1][2], Z_list_poly[ip1][3]])
        
        ang_i_ip1 = get_angle_vec(lambdas, i, ip1)
        ang_ip1_im1 = get_angle_vec(lambdas, ip1, im1)
        ang_im1_i = get_angle_vec(lambdas, im1, i)
        
        denom = ang_im1_i * ang_i_ip1
        if denom == 0:
            raise ValueError(f"Degenerate kinematics at i={i}")
            
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
    
    # 2. Evaluate M and N
    # Use generic references to avoid collisions
    # Randomize them to avoid accidental integer collisions
    # Ensure chi[0] != 0 to avoid collinearity with lambda_1=(0,eps)
    def r_val(): return QQ(ZZ.random_element(1, 100))
    
    theta = vector(R, [r_val(), r_val()])
    eta = vector(R, [r_val(), r_val()])
    chi = vector(R, [1, r_val()]) # Force chi[0]=1 != 0
    
    print(f"Chi: {chi}")
    print(f"Lambda_1: {k_eps.lambdas[1]}")
    
    try:
        sigmas = sigma_mhv(k_eps, theta, eta, chi)
        print("Sigmas computed.")
        # Debug collisions
        for i in range(6):
            for j in range(i+1, 6):
                diff = sigmas[i] - sigmas[j]
                if diff == 0:
                     print(f"Collision: sigma_{i} == sigma_{j}")
        
        det_p = detprime_phi(sigmas, k_eps)
        m_chy = mhv_gravity_amplitude(sigmas, k_eps, det_p)
        
        # Calculate D(epsilon)
        d_all = R(1)
        for i in range(6):
            for j in range(i+1, 6):
                d_all *= k_eps.angle(i, j)
                
        d_cyc = R(1)
        for i in range(6):
            d_cyc *= k_eps.angle(i, (i+1)%6)
            
        ang01 = k_eps.angle(0, 1)
        print(f"<01> = {ang01}")
        
        D_val = d_all * d_cyc / (ang01**2)
        
        N_val = m_chy * D_val
        
        def valuation(expr):
            if expr == 0: return float('inf')
            num = expr.numerator()
            den = expr.denominator()
            # Handle fraction field elements
            if hasattr(num, 'valuation'): 
                v_num = num.valuation()
            else:
                v_num = 0
                
            if hasattr(den, 'valuation'):
                v_den = den.valuation()
            else:
                v_den = 0
                
            return v_num - v_den
            
        v_01 = valuation(N_val)
        
        print(f"Vanishing order v_01 = {v_01}")
        
        # Interpret
        if v_01 == 0:
            print("Interpretation: <01>^2 is a real pole of M (N doesn't vanish).")
        elif v_01 > 0:
            print(f"Interpretation: N vanishes to order {v_01}. The <01>^2 pole is cancelled.")
            if v_01 >= 2:
                print("Specifically, <01>^2 in denominator is fully cancelled.")
        else:
            print(f"Interpretation: N has a pole at <01>=0 (Order {v_01}).")
            
        # Also check M order
        v_M = valuation(m_chy)
        print(f"Amplitude M order at <01>: {v_M}")
        
    except Exception as e:
        print(f"Failed Phase C0: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    analyze_01_vanishing_order()
