import sys
import os
from sage.all import *

if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())

from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import ang_bracket, sq_bracket
from src.chy_oracle.hodges_reduced import hodges_npt_mhv_canonical

def soft_factor_gravity_plus(lambdas, tilde_lambdas, s_idx, x_ref, y_ref):
    """
    Compute Soft Factor for soft particle s with helicity +2.
    Using 2-reference spinor formula to avoid gauge artifacts.
    Usually formula involves one reference eta?
    S^+ = - Sum_{a != s} ([a s] / <a s>) * (<a eta>^2 / <s eta>^2) ?
    
    If we use two references, maybe we mean checking independence?
    Or is there a specific 2-ref formula?
    Standard soft factor S^+ = Sum_a [a s]/<a s> * <a x><a y> / <s x><s y> ?
    
    The Weinberg soft factor is:
    S = Sum_a (epsilon_s . p_a)^2 / (p_s . p_a)
    For h=+2, epsilon_s = |ref] <s| / [s ref] ? No.
    epsilon^+_mu = <q gamma_mu s] / <q s> / sqrt(2).
    epsilon^+ . p_a = <q a s] / <q s> / sqrt(2) = <q a> [a s] / <q s>.
    (eps . p_a)^2 = <q a>^2 [a s]^2 / <q s>^2.
    p_s . p_a = <s a> [a s] / 2.
    S = Sum_a 2 * (<q a>^2 [a s]^2 / <q s>^2) / (<s a> [a s])
      = Sum_a 2 * <q a>^2 [a s] / (<q s>^2 <s a>)
      = - Sum_a [a s]/<a s> * (<a q>^2 / <s q>^2).
      
    This matches the standard formula with ONE reference q.
    We can accept x_ref and y_ref and check they give same result.
    """
    def compute_with_ref(eta):
        total = QQ(0)
        lts = tilde_lambdas[s_idx]
        ls = lambdas[s_idx]
        sq_s_eta = lts[0]*eta[1] - lts[1]*eta[0] # This is [s eta]
        # Wait, formula has <s q> in denominator.
        # Reference q is a spinor (lambda-like).
        # epsilon^+ = |q] <s| / <s q>. Ref is q (lambda).
        # But we pass x_ref as vector. Is it lambda or tilde?
        # If we use auxiliary lambda q:
        # S^+ = - Sum [a s]/<a s> * <a q>^2 / <s q>^2.
        
        # Let's assume ref is lambda-type.
        ang_s_q = ls[0]*eta[1] - ls[1]*eta[0]
        if ang_s_q == 0: return None
        
        for a in range(len(lambdas)):
            if a == s_idx: continue
            
            # <a s>
            ang_as = lambdas[a][0]*ls[1] - lambdas[a][1]*ls[0]
            # [a s]
            sq_as = tilde_lambdas[a][0]*lts[1] - tilde_lambdas[a][1]*lts[0]
            # <a q>
            ang_a_q = lambdas[a][0]*eta[1] - lambdas[a][1]*eta[0]
            
            if ang_as == 0: continue
            
            term = (sq_as / ang_as) * (ang_a_q**2 / ang_s_q**2)
            total -= term
        return total

    S1 = compute_with_ref(x_ref)
    S2 = compute_with_ref(y_ref)
    
    if S1 is None or S2 is None: return None
    
    if abs(S1 - S2) > 1e-9:
        print(f"  Warning: Soft factor gauge dependent! Diff={abs(S1-S2)}")
        
    return S1

def check_soft_limit_n6_mhv_plus(seed=42):
    print(f"M3: Checking Soft Limit (h=+2) for N=6 MHV (Seed {seed})...")
    
    # 1. Hard Kinematics (N=5)
    n_hard = 5
    lambdas_hard, tildes_hard = sample_spinors_from_twistor(seed=seed, n=n_hard)
    
    # Check M5
    M5_hard_val, status = hodges_npt_mhv_canonical(lambdas_hard, tildes_hard, negative_indices=(0, 1))
    if status != "ok":
        print("  M5 computation failed")
        return
    print(f"  M5: {M5_hard_val}")

    # 2. Add Soft Particle (index 5)
    s_lambda = vector(QQ, [1, 2])
    s_tilde = vector(QQ, [3, 4])
    
    # Holomorphic Soft Limit: lambda_s -> epsilon * lambda_s
    epsilon = QQ(1)/100000 # 1e-5
    
    # Recoil logic
    # p_s = |eps lam> [tilde|.
    # p_s = eps * |lam> [tilde|.
    # We absorb recoil into 0 and 1.
    ts = s_tilde
    t0 = tildes_hard[0]
    t1 = tildes_hard[1]
    
    det_01 = t0[0]*t1[1] - t0[1]*t1[0] # [0 1]
    rhs = -epsilon * ts
    det_rhs_t1 = rhs[0]*t1[1] - rhs[1]*t1[0]
    det_t0_rhs = t0[0]*rhs[1] - t0[1]*rhs[0]
    
    a = det_rhs_t1 / det_01
    b = det_t0_rhs / det_01
    
    lambdas_6 = []
    tildes_6 = []
    
    # 0
    lambdas_6.append(lambdas_hard[0] + a * s_lambda)
    tildes_6.append(tildes_hard[0])
    
    # 1
    lambdas_6.append(lambdas_hard[1] + b * s_lambda)
    tildes_6.append(tildes_hard[1])
    
    # 2, 3, 4
    for i in range(2, 5):
        lambdas_6.append(lambdas_hard[i])
        tildes_6.append(tildes_hard[i])
        
    # 5 (Soft)
    lambdas_6.append(epsilon * s_lambda)
    tildes_6.append(s_tilde)
    
    # Evaluate M6
    M6_val, status = hodges_npt_mhv_canonical(lambdas_6, tildes_6, negative_indices=(0, 1))
    
    # Soft Factor
    # S^+ uses HARD lambdas (shifted?)?
    # Usually Leading Soft Factor uses unshifted hard particles, or shifted?
    # The difference is O(epsilon).
    # But we should be precise.
    # The amplitude M6(eps) ~ S^+ M5 + O(eps^0).
    # Since S^+ ~ 1/eps, we check eps * M6 -> eps * S * M5.
    
    # Using shifted lambdas 0,1 inside S makes it more accurate?
    # Let's use lambdas_6 (which includes shifts).
    # Note: lambda_6[5] = epsilon * s_lambda.
    # S factor formula uses <a s>.
    # <a s> is proportional to epsilon.
    # [a s] is O(1).
    # S^+ ~ [as]/<as> ~ 1/eps.
    
    x_ref = vector(QQ, [1, 0])
    y_ref = vector(QQ, [0, 1])
    
    S_val = soft_factor_gravity_plus(lambdas_6, tildes_6, 5, x_ref, y_ref)
    
    # Prediction
    predicted = S_val * M5_hard_val
    
    print(f"  M6: {M6_val.n()}")
    print(f"  Pred: {predicted.n()}")
    
    if predicted == 0:
        print("  Predicted is 0.")
        return
        
    ratio = M6_val / predicted
    print(f"  Ratio: {ratio.n()}")
    
    if abs(ratio - 1) < 1e-3:
        print("  PASS: Soft limit confirmed.")
    else:
        print("  FAIL: Mismatch.")

def compute_soft_ratio_n6(epsilon, seed=42):
    """
    Compute ratio M6 / (S * M5) for a given epsilon.
    Returns (ratio, status).
    """
    n_hard = 5
    lambdas_hard, tildes_hard = sample_spinors_from_twistor(seed=seed, n=n_hard)
    M5_hard_val, status = hodges_npt_mhv_canonical(lambdas_hard, tildes_hard, negative_indices=(0, 1))
    if status != "ok": return None, f"M5_{status}"

    s_lambda = vector(QQ, [1, 2])
    s_tilde = vector(QQ, [3, 4])
    
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
    
    M6_val, status = hodges_npt_mhv_canonical(lambdas_6, tildes_6, negative_indices=(0, 1))
    if status != "ok": return None, f"M6_{status}"
    
    x_ref = vector(QQ, [1, 0])
    y_ref = vector(QQ, [0, 1])
    S_val = soft_factor_gravity_plus(lambdas_6, tildes_6, 5, x_ref, y_ref)
    
    if S_val is None: return None, "S_val_none"
    
    predicted = S_val * M5_hard_val
    if predicted == 0: return None, "zero_prediction"
    
    return M6_val / predicted, "ok"

if __name__ == "__main__":
    check_soft_limit_n6_mhv_plus()
