import sys
import os
from sage.all import *
import json
import random

# Add project root to path
sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor, ang_bracket
from src.chy_oracle.matrix_tree import hodges_weighted_laplacian
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor

def check_reference_independence():
    print("Starting Phase F1.2: Reference Independence Check...")
    
    n_seeds = 3
    n_refs_per_seed = 20
    
    results_log = []
    
    for seed in range(n_seeds):
        print(f"\nProcessing seed {seed}...")
        
        # 1. Sample Kinematics
        try:
            lambdas, tildes = sample_spinors_from_twistor(seed=seed, n=6)
        except Exception as e:
            print(f"  Failed to sample kinematics: {e}")
            continue
            
        # 2. Compute Physical Amplitude (Ground Truth)
        M_phys, status = hodges_6pt_mhv_spinor(lambdas, tildes)
        if status != "ok":
            print(f"  Skipping seed {seed}: {status}")
            continue
            
        # Precompute static factors for the identity
        # M = - <01>^8 * det(L_3x3) / ( prod(C^2) * norm )
        h_factor = ang_bracket(lambdas[0], lambdas[1])**8
        norm_factor = (ang_bracket(lambdas[0], lambdas[1]) * 
                       ang_bracket(lambdas[1], lambdas[2]) * 
                       ang_bracket(lambdas[2], lambdas[0]))**2
        
        # We need M_rec / M_phys == 1
        # M_rec = - h_factor * det(L) / (prod(C^2) * norm)
        
        seed_failures = 0
        max_dev = 0
        
        for ref_idx in range(n_refs_per_seed):
            # 3. Sample random reference spinors x, y
            # Use random rationals
            x = vector(QQ, [QQ(random.randint(-10, 10)), QQ(random.randint(-10, 10))])
            y = vector(QQ, [QQ(random.randint(-10, 10)), QQ(random.randint(-10, 10))])
            
            # Avoid 0 vector
            if x == 0: x = vector(QQ, [1, 0])
            if y == 0: y = vector(QQ, [0, 1])
            
            try:
                L_tilde, C, _ = hodges_weighted_laplacian(lambdas, tildes, x, y)
            except ValueError:
                # Orthogonal to some particle
                continue
                
            # Compute minor (delete 0,1,2 -> keep 3,4,5)
            indices = [3, 4, 5]
            det_L = L_tilde.matrix_from_rows_and_columns(indices, indices).det()
            
            prod_C_sq = 1
            for k in indices:
                prod_C_sq *= C[k]**2
                
            if prod_C_sq == 0: continue
            
            M_rec = - det_L * h_factor / (prod_C_sq * norm_factor)
            
            if M_phys == 0:
                print("  M_phys is 0 (unlikely for generic kinematics)")
                continue
                
            ratio = M_rec / M_phys
            dev = abs(ratio - 1)
            
            if dev > max_dev:
                max_dev = dev
                
            if dev > 1e-10:
                seed_failures += 1
                print(f"  FAIL at ref {ref_idx}: ratio={float(ratio)}")
                print(f"  x={x}, y={y}")
                
        results_log.append({
            "seed": seed,
            "max_deviation": float(max_dev),
            "failures": seed_failures,
            "refs_tested": n_refs_per_seed
        })
        
        if seed_failures == 0:
            print(f"  Seed {seed} PASSED (Max dev: {max_dev})")
        else:
            print(f"  Seed {seed} FAILED ({seed_failures} mismatches)")
            
    # Save log
    with open("reference_independence_results.json", "w") as f:
        json.dump(results_log, f, indent=2)
        
    print("\nVerification Complete. Log saved.")

if __name__ == "__main__":
    check_reference_independence()







