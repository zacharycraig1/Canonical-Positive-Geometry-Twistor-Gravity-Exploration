import sys
import os
# sys.path.append(os.getcwd())

from sage.all import *
load("src/boundary_sampler.sage")
load("src/oracle_gravity_mhv6.sage")
load("src/kinematics_map.sage")

def test_factorization_oracle():
    print("Running Phase 3: Factorization Oracle Test")
    print("------------------------------------------")
    
    # Test 3-particle channel factorization
    # s_012 -> 0.
    # M_6 ~ M_L(0,1,2,-P) * (1/s_012) * M_R(P,3,4,5)
    # M_6 * s_012 should approach a constant (the residue).
    
    channel = (0, 1, 2)
    epsilons = [QQ(1)/100, QQ(1)/1000, QQ(1)/10000, QQ(1)/100000]
    
    # We use a fixed seed to track the SAME kinematic point approaching the boundary
    # Actually sample_near_boundary generates a NEW random point each time if we don't fix the base.
    # To test scaling, we must use the SAME base point and shift along the SAME line.
    
    # Let's manually do the sampling loop here to control the base point
    import random
    random.seed(int(42))
    base_seed = random.randint(1, 100000)
    
    # Generate base once
    res = sample_spinor_helicity_conserving(n=6, seed=base_seed)
    lambdas_0, tilde_0 = res
    
    # Find shift
    S_set = set(channel)
    a = 0 # in S
    b = 3 # not in S
    z_pole = solve_boundary_shift(lambdas_0, tilde_0, S_set, (a, b))
    
    print(f"Base point generated. z_pole = {z_pole}")
    
    prev_val = None
    
    print(f"\nChecking channel {channel} scaling:")
    print(f"{'Epsilon':<15} {'s_S':<15} {'M_hodges':<25} {'Product (Residue)':<25}")
    
    for eps in epsilons:
        z = z_pole + eps
        l_eps, lt_eps = apply_shift(lambdas_0, tilde_0, a, b, z)
        
        # Compute amplitude
        # We assume oracle handles potential large numbers
        res = oracle_M6(l_eps, lt_eps)
        if res["value"] is None:
            print(f"Error at eps={eps}: {res.get('error')}")
            continue
            
        m_val = res["value"]
        
        # Compute s_S
        sij = spinors_to_sij(l_eps, lt_eps)
        s_S = sij[(0,1)] + sij[(1,2)] + sij[(0,2)]
        
        product = m_val * s_S
        
        print(f"{float(eps):<15.2e} {float(s_S):<15.2e} {float(abs(m_val)):<25.2e} {product}")
        
        if prev_val is not None:
            # Check convergence
            # Using rational exact arithmetic, so we can see digits
            # product should be stable
            diff = abs(product - prev_val)
            # print(f"  Diff from prev: {float(diff):.2e}")
            pass
            
        prev_val = product

    # Check 2-particle channel (collinear limit)
    # s_01 -> 0
    # Gravity amplitude scaling: M ~ 1/s_01^3? Or 1/s_01?
    # MHV Gravity collinear limit:
    # If 0,1 collinear: 1/s <01> ...
    # Standard splitting function behavior.
    # Let's just check if M * s_01 is finite or if we need higher powers.
    # Hint: Gravity usually has 1/s poles?
    # Wait, splitting function for gravity is more singular?
    # No, tree amplitude poles are simple poles 1/P^2.
    # Collinear limit is a pole in s_ij.
    # So M * s_ij should be finite.
    
    channel2 = (0, 1)
    print(f"\nChecking channel {channel2} scaling:")
    
    # New base for this channel
    z_pole2 = solve_boundary_shift(lambdas_0, tilde_0, set(channel2), (0, 3))
    
    for eps in epsilons:
        z = z_pole2 + eps
        l_eps, lt_eps = apply_shift(lambdas_0, tilde_0, 0, 3, z)
        
        res = oracle_M6(l_eps, lt_eps)
        if res["value"] is None: continue
        m_val = res["value"]
        
        # s_01
        s_01 = ang_bracket(l_eps, 0, 1) * sq_bracket(lt_eps, 0, 1)
        
        product = m_val * s_01
        print(f"{float(eps):<15.2e} {float(s_01):<15.2e} {float(abs(m_val)):<25.2e} {product}")

if __name__ == "__main__":
    test_factorization_oracle()

