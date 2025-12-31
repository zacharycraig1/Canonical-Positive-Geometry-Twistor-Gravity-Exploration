from sage.all import *

load("src/boundary_sampler.sage")
load("src/spinor_sampling.sage")
load("src/kinematics_map.sage")

# from src.boundary_sampler import solve_boundary_shift, apply_shift
# from src.spinor_sampling import sample_spinor_helicity_conserving
# from src.kinematics_map import spinors_to_sij

import random

def debug_boundary():
    random.seed(int(42))
    res = sample_spinor_helicity_conserving(n=6, seed=123)
    lambdas, tilde = res
    
    S = {0, 1, 2}
    a, b = 0, 3
    
    # 1. Calculate s_S(0)
    sij = spinors_to_sij(lambdas, tilde)
    s_S_0 = sij[(0,1)] + sij[(0,2)] + sij[(1,2)]
    print(f"s_S(0) = {s_S_0}")
    
    # 2. Calculate matrix element manually
    # <0 | P_S | 3] = <0 1>[1 3] + <0 2>[2 3]
    term1 = (lambdas[0][0]*lambdas[1][1] - lambdas[0][1]*lambdas[1][0]) * \
            (tilde[1][0]*tilde[3][1] - tilde[1][1]*tilde[3][0])
            
    term2 = (lambdas[0][0]*lambdas[2][1] - lambdas[0][1]*lambdas[2][0]) * \
            (tilde[2][0]*tilde[3][1] - tilde[2][1]*tilde[3][0])
            
    mat_elem = term1 + term2
    print(f"Matrix elem <0|P|3] = {mat_elem}")
    
    # 3. Predict z_pole
    # s(z) = s(0) - z * mat_elem = 0 => z = s(0)/mat_elem
    z_pred = s_S_0 / mat_elem
    print(f"Predicted z_pole = {z_pred}")
    
    # 4. Check solve_boundary_shift result
    z_solved = solve_boundary_shift(lambdas, tilde, S, (a,b))
    print(f"Solved z_pole = {z_solved}")
    
    # 5. Apply shift and check s_S
    l_new, t_new = apply_shift(lambdas, tilde, a, b, z_solved)
    
    sij_new = spinors_to_sij(l_new, t_new)
    s_S_new = sij_new[(0,1)] + sij_new[(0,2)] + sij_new[(1,2)]
    print(f"s_S(z_pole) = {s_S_new}")
    
    if s_S_new == 0:
        print("SUCCESS: Boundary reached exactly.")
    else:
        print(f"FAILURE: Residual s_S = {s_S_new}")
        
    # Check linearity
    l_half, t_half = apply_shift(lambdas, tilde, a, b, z_solved/2)
    sij_half = spinors_to_sij(l_half, t_half)
    s_S_half = sij_half[(0,1)] + sij_half[(0,2)] + sij_half[(1,2)]
    print(f"s_S(z_pole/2) = {s_S_half}")
    print(f"Expected linear half: {s_S_0 - (z_solved/2)*mat_elem}")
    
    # Check if quadratic term exists?
    # P(z)^2 = P(0)^2 - z <a|P|b] + z^2 k^2.
    # k = |0>[3|. k^2 = <0 0>[3 3] = 0.
    # But wait, k^2 = s_k = 0?
    # k is a vector. k^2 = k_mu k^mu.
    # |a>[b| vector is v_mu = 1/2 <a|sigma_mu|b].
    # v^2 = 0? Yes, standard identity.
    # So it should be strictly linear.

if __name__ == "__main__":
    debug_boundary()

