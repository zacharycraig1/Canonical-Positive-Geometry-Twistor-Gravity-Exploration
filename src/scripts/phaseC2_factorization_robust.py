import sys
import os
import json
import time
from sage.all import *

# Ensure src is in path
sys.path.append(os.path.join(os.getcwd(), 'src'))

# Load Hodges module for reference if needed
root_dir = os.getcwd()
hodges_path = os.path.join(root_dir, 'src', 'hodges.sage')
load(hodges_path)

# --- Helpers ---

def ang_bracket(la, lb):
    return la[0]*lb[1] - la[1]*lb[0]

def sq_bracket(lta, ltb):
    # [i j] = u0 v1 - u1 v0
    return lta[0]*ltb[1] - lta[1]*ltb[0]

def get_s_spinor(lambdas, tilde_lambdas, i, j):
    """Computes s_ij = <ij>[ji]"""
    ang = ang_bracket(lambdas[i], lambdas[j])
    sq = sq_bracket(tilde_lambdas[i], tilde_lambdas[j])
    return ang * sq

def get_s_multi_spinor(lambdas, tilde_lambdas, indices):
    """
    Computes s_S = (sum p_i)^2 = sum_{i<j in S} 2 p_i.p_j
    2 p_i.p_j = <ij>[ji] (or -<ij>[ij])
    We use get_s_spinor(i,j) which returns <ij>[ji] implicitly.
    """
    s_total = QQ(0)
    indices = sorted(list(indices))
    for i in range(len(indices)):
        for j in range(i+1, len(indices)):
            idx1 = indices[i]
            idx2 = indices[j]
            s_total += get_s_spinor(lambdas, tilde_lambdas, idx1, idx2)
    return s_total

def bcfw_shift_spinors(lambdas, tilde_lambdas, a, b, z):
    """
    Apply BCFW shift:
    la -> la + z lb
    lb_t -> lb_t - z la_t
    """
    new_L = [list(l) for l in lambdas]
    new_Lt = [list(lt) for lt in tilde_lambdas]
    
    # la
    new_L[a][0] += z * lambdas[b][0]
    new_L[a][1] += z * lambdas[b][1]
    
    # lb_t
    new_Lt[b][0] -= z * tilde_lambdas[a][0]
    new_Lt[b][1] -= z * tilde_lambdas[a][1]
    
    # Convert back to vectors
    return [vector(QQ, l) for l in new_L], [vector(QQ, lt) for lt in new_Lt]

def solve_s_channel_bcfw(lambdas, tilde_lambdas, channel_indices, shift_a, shift_b):
    """
    Find z* such that s_channel(z*) = 0.
    """
    # Evaluate at z=0 and z=1
    def get_s(z_val):
        L, Lt = bcfw_shift_spinors(lambdas, tilde_lambdas, shift_a, shift_b, z_val)
        return get_s_multi_spinor(L, Lt, channel_indices)

    s0 = get_s(QQ(0))
    s1 = get_s(QQ(1))
    
    slope = s1 - s0
    if slope == 0:
        return None # Constant
        
    z_star = -s0 / slope
    return z_star

def hodges_6pt_mhv_spinor(lambdas, tilde_lambdas, deletion_set=None):
    """
    Hodges formula for 6-point MHV gravity implemented directly with spinors.
    """
    n = 6
    if len(lambdas) != n:
        return None, "n_must_be_6"
        
    # Build Phi matrix
    Phi = matrix(QQ, n, n)
    
    # Off-diagonal: Phi_{ij} = [i j] / <i j>
    for i in range(n):
        for j in range(n):
            if i != j:
                ang = ang_bracket(lambdas[i], lambdas[j])
                sq = sq_bracket(tilde_lambdas[i], tilde_lambdas[j])
                
                if ang == 0:
                    return None, "domain_violation_angle_bracket_offdiag"
                    
                Phi[i, j] = sq / ang
                
    # Diagonal: Reference spinor formula
    lx = vector(QQ, [1, 0])
    ly = vector(QQ, [0, 1])
    
    # Check if valid references
    for i in range(n):
        if ang_bracket(lambdas[i], lx) == 0:
            lx = vector(QQ, [1, 1]) # fallback
        if ang_bracket(lambdas[i], ly) == 0:
            ly = vector(QQ, [1, -1]) # fallback

    for i in range(n):
        ang_ix = ang_bracket(lambdas[i], lx)
        ang_iy = ang_bracket(lambdas[i], ly)
        
        if ang_ix == 0 or ang_iy == 0:
             return None, "domain_violation_reference_spinor"
             
        diag_sum = QQ(0)
        for j in range(n):
            if j == i: continue
            
            ang_jx = ang_bracket(lambdas[j], lx)
            ang_jy = ang_bracket(lambdas[j], ly)
            
            term = Phi[i, j] * (ang_jx * ang_jy) / (ang_ix * ang_iy)
            diag_sum -= term
            
        Phi[i, i] = diag_sum
        
    # Reduced determinant
    if deletion_set is None:
        rows_to_delete = [0, 1, 2]
    else:
        rows_to_delete = deletion_set
        
    cols_to_delete = rows_to_delete
    
    rows_keep = [r for r in range(n) if r not in rows_to_delete]
    cols_keep = [c for c in range(n) if c not in cols_to_delete]
    
    Phi_red = Phi[rows_keep, cols_keep]
    det_Phi_red = Phi_red.det()
    
    # Normalization factor
    r1, r2, r3 = rows_to_delete
    ang_12 = ang_bracket(lambdas[r1], lambdas[r2])
    ang_23 = ang_bracket(lambdas[r2], lambdas[r3])
    ang_31 = ang_bracket(lambdas[r3], lambdas[r1])
    
    norm_factor = (ang_12 * ang_23 * ang_31)**2
    if norm_factor == 0:
        return None, "domain_violation_norm_factor"
        
    det_prime = det_Phi_red / norm_factor
    
    # Helicity factor <0 1>^8
    h_factor = ang_bracket(lambdas[0], lambdas[1])**8
    
    return det_prime * h_factor, "ok"

def run_all_channels_classification():
    print("Phase C2 (Redux) - Channel Classification [Optimized]")
    
    # Channels: All 3-particle subsets of {0..5}
    # But s_S = s_{complement}. So {0,1,2} is same channel as {3,4,5}.
    # We only need to check half of them.
    # Total subsets of size 3 from 6 is 6C3 = 20.
    # We have 10 unique channels (pairs of complement sets).
    
    from itertools import combinations
    indices = set(range(6))
    
    # Generate unique channels
    seen_channels = set()
    channels_to_test = []
    
    for comb in combinations(indices, 3):
        c = tuple(sorted(comb))
        complement = tuple(sorted(list(indices - set(c))))
        
        if c in seen_channels or complement in seen_channels:
            continue
            
        seen_channels.add(c)
        seen_channels.add(complement)
        channels_to_test.append(c)
        
    print(f"Testing {len(channels_to_test)} unique 3-particle channels: {channels_to_test}")
    
    log_data = []
    
    epsilons = [QQ(1)/10**k for k in range(3, 9)]
    num_trials_per_channel = 3
    
    seed_counter = 6000
    
    for channel in channels_to_test:
        print(f"\n=== Channel {channel} ===")
        
        # Pick shift indices (a, b) such that ONE is in channel, ONE is out.
        # This ensures momentum flows through the channel cut.
        # If both in or both out, s_channel is invariant (constant).
        
        in_set = set(channel)
        out_set = indices - in_set
        
        shift_a = list(in_set)[0]
        shift_b = list(out_set)[0]
        
        print(f"Shift indices: ({shift_a}, {shift_b})")
        
        valid_trials = 0
        while valid_trials < num_trials_per_channel:
            seed = seed_counter
            seed_counter += 1
            
            tw = MomentumTwistor(n=6, seed=seed)
            if not tw.domain_ok: continue
            
            lambdas = [tw.get_lambda(i) for i in range(6)]
            tilde_lambdas = [tw.get_tilde_lambda(i) for i in range(6)]
            if any(x is None for x in tilde_lambdas): continue
            
            # Solve z*
            z_star = solve_s_channel_bcfw(lambdas, tilde_lambdas, channel, shift_a, shift_b)
            if z_star is None:
                continue # constant channel?
                
            # Probe
            scaling_results = []
            for eps in epsilons:
                z_probe = z_star + eps
                L_probe, Lt_probe = bcfw_shift_spinors(lambdas, tilde_lambdas, shift_a, shift_b, z_probe)
                
                s_val = get_s_multi_spinor(L_probe, Lt_probe, channel)
                M_val, reason = hodges_6pt_mhv_spinor(L_probe, Lt_probe)
                
                if M_val is None:
                    scaling_results.append({'eps': float(eps), 's': float(s_val) if s_val else None, 'M': None, 'error': reason})
                else:
                    scaling_results.append({'eps': float(eps), 's': float(s_val), 'M': float(abs(M_val))})
            
            # Estimate alpha
            valid_probes = [d for d in scaling_results if d['M'] is not None and d['s'] != 0 and d['M'] != 0]
            alpha_est = "N/A"
            classification = "Unknown"
            
            if len(valid_probes) >= 2:
                p1 = valid_probes[-1]
                p2 = valid_probes[0]
                import math
                try:
                    logs1 = math.log(abs(p1['s']))
                    logs2 = math.log(abs(p2['s']))
                    logm1 = math.log(p1['M'])
                    logm2 = math.log(p2['M'])
                    if abs(logs1 - logs2) > 1e-9:
                        alpha = (logm1 - logm2) / (logs1 - logs2)
                        alpha_est = f"{alpha:.4f}"
                        
                        if abs(alpha + 1) < 0.1:
                            classification = "Pole (1/s)"
                        elif abs(alpha) < 0.1:
                            classification = "Finite (s^0)"
                        elif alpha > 0.5:
                            classification = "Zero (s^1 or higher)"
                        else:
                            classification = f"Other ({alpha:.2f})"
                except:
                    pass
            
            print(f"Trial {valid_trials+1}: z*={z_star}, alpha={alpha_est} -> {classification}")
            
            log_data.append({
                'channel': str(channel),
                'seed': seed,
                'z_star': str(z_star),
                'alpha_est': alpha_est,
                'classification': classification,
                'scaling_data': scaling_results
            })
            
            valid_trials += 1

    with open('channel_classification_results.json', 'w') as f:
        json.dump(log_data, f, indent=2)
    print("\nResults saved to channel_classification_results.json")

if __name__ == "__main__":
    run_all_channels_classification()
