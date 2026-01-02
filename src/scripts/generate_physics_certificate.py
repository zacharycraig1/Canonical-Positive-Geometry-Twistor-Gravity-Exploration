import sys
import os
import json
from sage.all import *

# Ensure we can import from src
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import hodges_npt_mhv_spinor
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical
from src.physics_limits.bcfw import bcfw_shift_spinors, solve_bcfw_pole, get_channel_s, get_momentum_matrix, decompose_momentum_spinors
from src.physics_limits.soft import soft_factor_gravity_plus

def check_factorization_n6_s012(seed):
    """Check factorization for s_012 (should be zero for MHV)."""
    lambdas, tildes = sample_spinors_from_twistor(seed=seed, n=6)
    channel = [0, 1, 2]
    # Shift 0, 3
    z_star = solve_bcfw_pole(lambdas, tildes, 0, 3, channel)
    if z_star is None: return {"error": "no_pole"}
    
    eps = QQ(1)/1000000
    z_probe = z_star + eps
    L, Lt = bcfw_shift_spinors(lambdas, tildes, 0, 3, z_probe)
    
    M6_res = hodges_npt_mhv_canonical(L, Lt, (0, 1))
    if M6_res[1] != "ok":
        return {"error": "M6_failed", "status": M6_res[1]}
    M6 = M6_res[0]
    
    s_val = get_channel_s(L, Lt, channel)
    res_numeric = M6 * s_val
    
    # We expect zero
    passed = abs(res_numeric) < 1e-6
    return {
        "channel": "s_012",
        "z_star": str(z_star),
        "residue_numeric": str(float(res_numeric)),
        "expected": "0",
        "pass": passed
    }

def check_soft_plus(seed):
    """Check soft limit h=+2."""
    n_hard = 5
    lambdas_hard, tildes_hard = sample_spinors_from_twistor(seed=seed, n=n_hard)
    M5_res = hodges_npt_mhv_canonical(lambdas_hard, tildes_hard, (0, 1))
    if M5_res[1] != "ok":
        return {"error": "M5_hard_failed"}
    M5_hard = M5_res[0]
    
    s_lambda = vector(QQ, [1, 2])
    s_tilde = vector(QQ, [3, 4])
    epsilon = QQ(1)/100000
    
    # Recoil
    ts = s_tilde
    t0 = tildes_hard[0]
    t1 = tildes_hard[1]
    det_01 = t0[0]*t1[1] - t0[1]*t1[0]
    rhs = -epsilon * ts
    det_rhs_t1 = rhs[0]*t1[1] - rhs[1]*t1[0]
    det_t0_rhs = t0[0]*rhs[1] - t0[1]*rhs[0]
    a = det_rhs_t1 / det_01
    b = det_t0_rhs / det_01
    
    lambdas_6 = []
    tildes_6 = []
    lambdas_6.append(lambdas_hard[0] + a * s_lambda)
    tildes_6.append(tildes_hard[0])
    lambdas_6.append(lambdas_hard[1] + b * s_lambda)
    tildes_6.append(tildes_hard[1])
    for i in range(2, 5):
        lambdas_6.append(lambdas_hard[i])
        tildes_6.append(tildes_hard[i])
    lambdas_6.append(epsilon * s_lambda)
    tildes_6.append(s_tilde)
    
    M6_res = hodges_npt_mhv_canonical(lambdas_6, tildes_6, (0, 1))
    if M6_res[1] != "ok":
        return {"error": "M6_soft_failed"}
    M6 = M6_res[0]
    
    x_ref = vector(QQ, [1, 0])
    y_ref = vector(QQ, [0, 1])
    S_val = soft_factor_gravity_plus(lambdas_6, tildes_6, 5, x_ref, y_ref)
    
    predicted = S_val * M5_hard
    ratio = M6 / predicted if predicted != 0 else 0
    
    passed = abs(ratio - 1) < 1e-3
    
    return {
        "type": "soft_plus",
        "ratio": str(float(ratio)),
        "pass": passed
    }

def generate_certificate():
    print("M4: Generating Physics Certificate...")
    results = {}
    
    # 1. Factorization Checks
    fact_results = []
    for seed in range(3):
        res = check_factorization_n6_s012(seed)
        fact_results.append(res)
    results["factorization_s012"] = fact_results
    
    # 2. Soft Checks
    soft_results = []
    for seed in range(3):
        res = check_soft_plus(seed)
        soft_results.append(res)
    results["soft_plus"] = soft_results
    
    # Write JSON
    os.makedirs("RESULTS", exist_ok=True)
    with open("RESULTS/physics_limits_n6.json", "w") as f:
        json.dump(results, f, indent=2)
        
    print("Certificate written to RESULTS/physics_limits_n6.json")
    
    # Summary
    all_pass = True
    for r in fact_results: 
        if "error" in r or not r.get("pass", False): all_pass = False
    for r in soft_results:
        if "error" in r or not r.get("pass", False): all_pass = False
        
    if all_pass:
        print("ALL TESTS PASSED.")
    else:
        print("SOME TESTS FAILED.")

if __name__ == "__main__":
    generate_certificate()

