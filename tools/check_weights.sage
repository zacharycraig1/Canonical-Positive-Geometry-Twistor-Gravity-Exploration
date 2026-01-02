import sys
import os
# sys.path.append(os.getcwd())

from sage.all import *
load("src/oracle_gravity_mhv6.sage")
load("src/spinor_sampling.sage")

def check_weights():
    print("Checking Little Group Weights...")
    
    # Sample base point
    res = sample_spinor_helicity_conserving(n=6, seed=42)
    if res is None:
        print("Sampling failed")
        return
    lambdas, tildes = res
    
    # Compute base values
    res_base = oracle_M6(lambdas, tildes)
    val_hodges = res_base["M_hodges_reduced"]
    val_klt = res_base["M_klt"]
    
    if val_hodges is None or val_klt is None:
        print("Base computation failed")
        return
        
    print(f"Base Hodges: {val_hodges}")
    print(f"Base KLT:    {val_klt}")
    
    # Scale ONE particle's spinors to check individual weight
    # Particle 0
    # lambda_0 -> t * lambda_0
    # tilde_0 -> t^-1 * tilde_0
    # This preserves momentum p_0 = lambda*tilde
    # Gravity amplitudes scale as h_i = -2 for negative, +2 for positive helicity?
    # M(t*lam, t^-1*til) = t^(-2h) * M ?
    
    # Actually, let's just scale ALL lambdas by T and ALL tildes by 1/T.
    # Momentum is invariant.
    # M should scale by total weight.
    
    T = QQ(2)
    lambdas_scaled = [l * T for l in lambdas]
    tildes_scaled = [lt * (1/T) for lt in tildes] # actually we don't need to scale tildes to keep P invariant?
    # Wait. p = lam * til. If lam -> T lam, til -> 1/T til, then p -> p.
    # So P is invariant.
    
    res_scaled = oracle_M6(lambdas_scaled, tildes_scaled)
    h_scaled = res_scaled["M_hodges_reduced"]
    k_scaled = res_scaled["M_klt"]
    
    ratio_h = h_scaled / val_hodges
    ratio_k = k_scaled / val_klt
    
    print(f"\nScaling ALL: lambda->2*lambda, tilde->1/2*tilde (P invariant)")
    print(f"Hodges ratio: {ratio_h} = 2^{log(ratio_h, 2)}")
    print(f"KLT ratio:    {ratio_k} = 2^{log(ratio_k, 2)}")
    
    # Now scale uniformly: lambda -> T*lambda, tilde -> T*tilde
    # P -> T^2 * P.
    # Amplitude dimension? M ~ P^2 (dimensionless in some conventions, or mass dimension?)
    # KLT involves s_ij. Parke-Taylor involves 1/ang^2 ~ 1/T^2.
    
    lambdas_uni = [l * T for l in lambdas]
    tildes_uni = [lt * T for lt in tildes]
    
    res_uni = oracle_M6(lambdas_uni, tildes_uni)
    h_uni = res_uni["M_hodges_reduced"]
    k_uni = res_uni["M_klt"]
    
    r_uni_h = h_uni / val_hodges
    r_uni_k = k_uni / val_klt
    
    print(f"\nScaling Uniform: lambda->2*lambda, tilde->2*tilde (P -> 4*P)")
    print(f"Hodges ratio: {r_uni_h} = 2^{log(r_uni_h, 2)}")
    print(f"KLT ratio:    {r_uni_k} = 2^{log(r_uni_k, 2)}")
    
    # Check little group weight for Particle 0 (assumed negative helicity in KLT?)
    # lambda_0 -> T * lambda_0, tilde_0 -> T^-1 * tilde_0. (P0 fixed)
    # Check weight.
    l_lg = [v for v in lambdas]
    t_lg = [v for v in tildes]
    l_lg[0] = l_lg[0] * T
    t_lg[0] = t_lg[0] * (1/T)
    
    res_lg = oracle_M6(l_lg, t_lg)
    h_lg = res_lg["M_hodges_reduced"]
    k_lg = res_lg["M_klt"]
    
    if h_lg is not None:
        print(f"\nLittle Group Particle 0 (lambda->2*lambda, tilde->1/2*tilde):")
        print(f"Hodges ratio: {h_lg/val_hodges} = 2^{log(h_lg/val_hodges, 2)}")
        print(f"KLT ratio:    {k_lg/val_klt}    = 2^{log(k_lg/val_klt, 2)}")

if __name__ == "__main__":
    check_weights()










