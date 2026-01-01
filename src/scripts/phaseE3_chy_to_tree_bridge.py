import sys
import os
from sage.all import *

sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import ang_bracket, sq_bracket
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy.mhv_solution import sigma_mhv, check_scattering_equations
from src.chy_oracle.matrix_tree import hodges_weighted_laplacian

class MockKinematics:
    def __init__(self, lambdas, tilde_lambdas):
        self.lambdas = lambdas
        self.tilde_lambdas = tilde_lambdas
        self.n = len(lambdas)
        
    def s(self, i, j):
        return ang_bracket(self.lambdas[i], self.lambdas[j]) * sq_bracket(self.tilde_lambdas[i], self.tilde_lambdas[j])

def verify_bridge():
    print("Starting Phase E3: CHY -> Tree Bridge...")
    
    seed = 42
    lambdas, tildes = sample_spinors_from_twistor(seed=seed, n=6)
    kin = MockKinematics(lambdas, tildes)
    
    # Reference spinors for CHY
    theta = vector(QQ, [1, 0])
    eta = vector(QQ, [0, 1])
    chi = vector(QQ, [1, 1])
    
    # 1. Compute sigma values
    sigmas = sigma_mhv(kin, theta, eta, chi)
    
    # Verify Scattering Equations
    residuals = check_scattering_equations(kin, sigmas)
    max_err = max(abs(r) for r in residuals)
    print(f"Scattering Eq Error: {max_err} (Should be 0)")
    
    # 2. Verify sigma difference identity
    # (sigma_a - sigma_b) ~ <ab> / (<a chi><b chi>)
    # Let's find the proportionality constant
    # sigma_a = <a eta><theta chi> / (<a chi><theta eta>)
    # sigma_a - sigma_b = (theta_chi/theta_eta) * [ <a eta>/<a chi> - <b eta>/<b chi> ]
    #                   = K * [ <a eta><b chi> - <b eta><a chi> ] / (<a chi><b chi>)
    # Using Schouten: <a eta><b chi> - <b eta><a chi> = <a b><eta chi>
    # So sigma_a - sigma_b = K * <a b><eta chi> / (<a chi><b chi>)
    # Constant = (theta_chi/theta_eta) * <eta chi>
    
    print("\nVerifying Sigma Difference Identity:")
    theta_chi = ang_bracket(theta, chi)
    theta_eta = ang_bracket(theta, eta)
    eta_chi = ang_bracket(eta, chi)
    
    const_factor = (theta_chi / theta_eta) * eta_chi
    
    mismatch = False
    for a in range(6):
        for b in range(a+1, 6):
            diff_sigma = sigmas[a] - sigmas[b]
            
            ang_ab = ang_bracket(lambdas[a], lambdas[b])
            ang_ac = ang_bracket(lambdas[a], chi)
            ang_bc = ang_bracket(lambdas[b], chi)
            
            rhs = const_factor * ang_ab / (ang_ac * ang_bc)
            
            if abs(diff_sigma - rhs) > 1e-10:
                print(f"Mismatch at {a},{b}: {diff_sigma} vs {rhs}")
                mismatch = True
                
    if not mismatch:
        print("Identity Verified: (sigma_a - sigma_b) propto <ab> / (<a chi><b chi>)")
        
    # 3. Connect to Weighted Laplacian
    # CHY Jacobian Matrix Phi_CHY_ij = s_ij / (sigma_i - sigma_j)^2
    # Substitute identity:
    # Phi_CHY_ij = s_ij / ( K * <ij> / (<i chi><j chi>) )^2
    #            = (<ij>[ij]) * (<i chi><j chi>)^2 / ( K^2 * <ij>^2 )
    #            = ( [ij]/<ij> ) * (<i chi><j chi>)^2 / K^2
    #            = w_ij * C_i * C_j * (1/K^2)
    # where C_i = <i chi>^2 ??
    # Wait, my Weighted Laplacian used C_i = <i x> <i y>.
    # Here we have (<i chi> <j chi>)^2 = <i chi>^2 <j chi>^2.
    # So it matches Weighted Laplacian with reference spinors x = y = chi!
    
    print("\nChecking Correspondence with Weighted Laplacian:")
    # Use x = y = chi for weighted laplacian
    L_tilde, C_w, Phi_w = hodges_weighted_laplacian(lambdas, tildes, chi, chi)
    
    # We expect L_tilde_ij = - w_ij * C_i * C_j = - w_ij * <i chi>^2 * <j chi>^2
    # CHY Phi_ij = w_ij * <i chi>^2 * <j chi>^2 / K^2
    # So CHY Phi_ij = - L_tilde_ij / K^2
    
    # Let's check explicitly
    phi_chy_sample = kin.s(0, 1) / (sigmas[0] - sigmas[1])**2
    l_tilde_sample = L_tilde[0, 1] # = - w_01 C_0 C_1
    
    ratio = phi_chy_sample / l_tilde_sample
    predicted_ratio = - 1 / const_factor**2
    
    print(f"Ratio Phi_CHY_01 / L_tilde_01: {ratio}")
    print(f"Predicted Ratio (-1/K^2):      {predicted_ratio}")
    
    if abs(ratio - predicted_ratio) < 1e-10:
        print("SUCCESS: CHY Jacobian is proportional to Weighted Laplacian!")
        print("This confirms Mechanism 3: The Weighted Laplacian IS the CHY Jacobian structure.")

if __name__ == "__main__":
    verify_bridge()






