import sys
import os
import json
import math
import random as rnd
from sage.all import *

sys.path.append(os.getcwd())

# Import MomentumTwistor from kinematics_samples
from src.chy_oracle.kinematics_samples import MomentumTwistor, sample_spinors_from_twistor

from src.chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor, ang_bracket
from src.chy_oracle.matrix_tree import hodges_weighted_laplacian

def measure_slopes(boundary_pair, n_seeds=3):
    i_idx, j_idx = boundary_pair
    print(f"\nAnalyzing Boundary <{i_idx},{j_idx}>...")
    
    epsilons = [QQ(1)/10**k for k in range(3, 8)] # 1e-3 .. 1e-7
    
    slopes = []
    
    for seed in range(n_seeds):
        # 1. Generate singular configuration
        # Use simple twistor construction: Z_j = Z_i => <ij>=0
        mt = MomentumTwistor(n=6, seed=seed+1000)
        
        # Force Z_j to have same spinor part as Z_i
        Z = [list(z) for z in mt.Z]
        Z[j_idx][0] = Z[i_idx][0]
        Z[j_idx][1] = Z[i_idx][1]
        
        # Perturbation direction V
        V = [rnd.randint(-10, 10) for _ in range(4)]
        # Ensure V spinor part is not parallel to Z_i spinor part
        # i.e. <i, V> != 0
        ang_iV = Z[i_idx][0]*V[1] - Z[i_idx][1]*V[0]
        if ang_iV == 0: V[1] += 1
        V = vector(QQ, V)
        
        data_points = []
        
        for eps in epsilons:
            # Perturb Z_j
            Z_eps = [vector(QQ, z) for z in Z]
            Z_eps[j_idx] = Z_eps[j_idx] + eps * V
            
            mt_eps = MomentumTwistor(n=6, Z=Z_eps)
            
            # Check angle
            ang_val = mt_eps.get_angle(i_idx, j_idx)
            if ang_val == 0: continue
            
            # Get spinors
            lambdas = [mt_eps.get_lambda(k) for k in range(6)]
            tildes = [mt_eps.get_tilde_lambda(k) for k in range(6)]
            
            if any(t is None for t in tildes): continue
            
            # 1. M_phys
            M_phys, status = hodges_6pt_mhv_spinor(lambdas, tildes)
            if status != "ok" or M_phys == 0: continue
            
            # 2. Geometric Object (Minor)
            # Use ref x, y
            x = vector(QQ, [1, 0])
            y = vector(QQ, [0, 1])
            try:
                L_tilde, C, _ = hodges_weighted_laplacian(lambdas, tildes, x, y)
            except ValueError: continue
            
            # Minor deleting {0,1,2}
            indices = [3, 4, 5]
            minor_det = L_tilde.matrix_from_rows_and_columns(indices, indices).det()
            
            # 3. Prefactors
            # Global factor G = <01>^8
            # Norm factor N = (<01><12><20>)^2
            # Vertex factor P = prod_{k in 3,4,5} C_k^2
            
            G = ang_bracket(lambdas[0], lambdas[1])**8
            
            N_val = (ang_bracket(lambdas[0], lambdas[1]) * 
                     ang_bracket(lambdas[1], lambdas[2]) * 
                     ang_bracket(lambdas[2], lambdas[0]))**2
                     
            P_val = 1
            for k in indices:
                P_val *= C[k]**2
                
            # Full Reconstructed M = - G * minor_det / (P * N)
            
            data_points.append({
                'log_ang': float(log(abs(ang_val))),
                'log_M': float(log(abs(M_phys))),
                'log_det': float(log(abs(minor_det))),
                'log_G': float(log(abs(G))),
                'log_N': float(log(abs(N_val))),
                'log_P': float(log(abs(P_val)))
            })
            
        # Debug print
        if data_points:
            print(f"  Seed {seed} debug:")
            p_start = data_points[0]
            p_end = data_points[-1]
            print(f"    Eps large (1e-3): log_ang={p_start['log_ang']:.2f}, log_M={p_start['log_M']:.2f}")
            print(f"    Eps small (1e-7): log_ang={p_end['log_ang']:.2f}, log_M={p_end['log_M']:.2f}")
            dx = p_start['log_ang'] - p_end['log_ang']
            dy = p_start['log_M'] - p_end['log_M']
            print(f"    dx={dx:.2f}, dy={dy:.2f}, slope={dy/dx:.2f}")

        # Fit slopes
        if len(data_points) >= 2:
            p1 = data_points[0] # largest eps (approx) - actually list is sorted large eps to small eps (1e-3 to 1e-7)
            p2 = data_points[-1] # smallest eps
            
            dx = p1['log_ang'] - p2['log_ang'] # Note: 1e-3 has larger log than 1e-7. 
            # Wait, epsilons are 1e-3, 1e-4... so p1 is 1e-3, p2 is 1e-7.
            # ang ~ eps. log(ang1) > log(ang2). dx > 0.
            
            seed_slopes = {}
            if abs(dx) > 1e-6:
                for key in ['M', 'det', 'G', 'N', 'P']:
                    dy = p1[f'log_{key}'] - p2[f'log_{key}']
                    seed_slopes[key] = dy / dx
                slopes.append(seed_slopes)
                
    # Average slopes
    if not slopes:
        return None
        
    avg_slopes = {}
    for key in slopes[0].keys():
        vals = [s[key] for s in slopes]
        avg_slopes[key] = sum(vals) / len(vals)
        
    return avg_slopes

def run_pole_audit():
    results = []
    
    # Check Adjacent <0,1>
    print("Checking Adjacent Pole <0,1>...")
    slopes_adj = measure_slopes((0, 1))
    if slopes_adj:
        print("  Slopes:", slopes_adj)
        results.append({"type": "adjacent", "pair": [0,1], "slopes": slopes_adj})
    
    # Check Non-Adjacent <0,2>
    print("Checking Non-Adjacent Pole <0,2>...")
    slopes_non = measure_slopes((0, 2))
    if slopes_non:
        print("  Slopes:", slopes_non)
        results.append({"type": "non_adjacent", "pair": [0,2], "slopes": slopes_non})
        
    # Check Kept Index Pole <3,4>
    print("Checking Kept Index Pole <3,4>...")
    slopes_kept = measure_slopes((3, 4))
    if slopes_kept:
        print("  Slopes:", slopes_kept)
        results.append({"type": "kept_index", "pair": [3,4], "slopes": slopes_kept})

    with open("results/phaseF_pole_orders.json", "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    run_pole_audit()
