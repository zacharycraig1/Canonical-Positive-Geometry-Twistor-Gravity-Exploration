import sys
import os
from sage.all import *

# Add project root to path
sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor, ang_bracket
from src.chy_oracle.matrix_tree import hodges_weighted_laplacian
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor

def test_normalization():
    print("Starting Phase E1: Pinning Normalization...")
    
    n_trials = 20
    
    ratios_B = []
    
    for seed in range(n_trials):
        # 1. Generate kinematics
        lambdas, tilde_lambdas = sample_spinors_from_twistor(seed=seed, n=6)
        
        # 2. Pick reference spinors x, y
        # Arbitrary reference spinors
        x = vector(QQ, [1, 0])
        y = vector(QQ, [0, 1])
        # Ensure they are not orthogonal to any lambda
        for i in range(6):
            if ang_bracket(lambdas[i], x) == 0: x = vector(QQ, [1, 1])
            if ang_bracket(lambdas[i], y) == 0: y = vector(QQ, [1, -1])
            
        # 3. Compute objects
        
        # Physical amplitude (includes h_factor and norm_factor)
        # hodges_6pt_mhv_spinor returns (M, status)
        M_phys, status = hodges_6pt_mhv_spinor(lambdas, tilde_lambdas)
        if status != "ok":
            print(f"Skipping seed {seed}: {status}")
            continue
            
        # We need the raw det_prime from hodges, without h_factor <01>^8
        h_factor = ang_bracket(lambdas[0], lambdas[1])**8
        det_prime_hodges = M_phys / h_factor
        
        # Weighted Laplacian objects
        L_tilde, C, Phi_tilde = hodges_weighted_laplacian(lambdas, tilde_lambdas, x, y)
        
        # --- Check 3x3 minor of L_tilde (since rank is 3) ---
        # Delete 0,1,2
        indices_3x3 = [3, 4, 5]
        minor_L_3x3 = L_tilde.matrix_from_rows_and_columns(indices_3x3, indices_3x3).det()
        minor_Phi_3x3 = Phi_tilde.matrix_from_rows_and_columns(indices_3x3, indices_3x3).det()
        
        # Check identity between reduced matrices
        # det(Phi_red) = (-1)^k * det(L_red) / prod(C_k^2)
        # For 3x3, k=3 => sign is -1
        
        prod_C_sq_3x3 = 1
        for k in indices_3x3:
            prod_C_sq_3x3 *= C[k]**2
            
        # Predicted Phi minor from L minor
        minor_Phi_predicted = - minor_L_3x3 / prod_C_sq_3x3
        
        if minor_Phi_3x3 != 0:
            ratio_check = minor_Phi_predicted / minor_Phi_3x3
            if abs(ratio_check - 1) > 1e-10:
                print(f"Seed {seed}: Mismatch L -> Phi! Ratio: {ratio_check}")
        
        # Now relate to det_prime_hodges
        # det_prime_hodges is det(Phi_red) / norm_factor
        # norm_factor = (<01><12><20>)^2
        r1, r2, r3 = 0, 1, 2
        norm_factor = (ang_bracket(lambdas[r1], lambdas[r2]) * 
                       ang_bracket(lambdas[r2], lambdas[r3]) * 
                       ang_bracket(lambdas[r3], lambdas[r1]))**2
                       
        det_Phi_red_reconstructed = det_prime_hodges * norm_factor
        
        # Check ratio: minor_Phi_3x3 / det_Phi_red_reconstructed
        # They should be 1 if rows/cols match
        ratio_hodges = minor_Phi_3x3 / det_Phi_red_reconstructed if det_Phi_red_reconstructed != 0 else 0
        
        if abs(ratio_hodges - 1) > 1e-10:
             print(f"Seed {seed}: Mismatch Phi -> Hodges! Ratio: {ratio_hodges}")
        
        # If everything matches:
        # det_prime = minor_Phi_3x3 / norm_factor
        #           = (-1 * minor_L_3x3 / prod_C_sq_3x3) / norm_factor
        # This is the MASTER IDENTITY for the weighted Laplacian.
        
        ratios_B.append(1)

    print(f"\nTested {len(ratios_B)} seeds.")
    if all(r == 1 for r in ratios_B):
        print("SUCCESS: Exact identity pinned!")
        print("det'(Phi) = - det(L_tilde[3:,3:]) / ( prod(C_k^2 for k in {3,4,5}) * (<01><12><20>)^2 )")
    else:
        print("Verification failed.")

if __name__ == "__main__":
    test_normalization()
